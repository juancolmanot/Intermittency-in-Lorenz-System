#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <omp.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "../include/lorenz.h"
#include "../../../common/include/file_handle.h"
#include "../../../common/include/ini.h"
#include "../../../common/include/parameters.h"
#include "../../../common/include/stats.h"
#include "../../../common/include/linear_algebra.h"
#include "../../../common/include/gsl_utilities.h"

int main(int argc, char *argv[]) {

    // ===============================================================================
    // Check for right passed arguments.
    // ===============================================================================
    if (argc != 4){
        perror("Wrong amount of arguments.");
        printf("Usage: ./run-script.sh script_to_run write_to_file.dat params_file.ini\n");
        exit(EXIT_FAILURE);
    }

    // ===============================================================================
    // Load arguments into variables.
    // ===============================================================================
    const char *write_filename = argv[2];
    const char *params_file = argv[3];

    // ===============================================================================
    // We load parameters in params1 variable
    // ===============================================================================
    Parameters5 params;
    load_parameters_from_file(params_file, &params, handler5);

    // ===============================================================================
    // File to write to
    // ===============================================================================
    FILE *f = open_file(write_filename);

    // ===============================================================================
    // Integration parameters.
    // ===============================================================================
    double h = params.h;
    double atol, rtol;
    atol = params.atol;
    rtol = params.rtol;
    double t, t_transient, t_stationary;
    t = 0.0;
    t_transient = params.t_transient;
    t_stationary = t_transient + params.t_stationary;

    // ===============================================================================
    // Initial state.
    // ===============================================================================
    size_t n = 3;
    double x0[n];
    x0[0] = 2.0;
    x0[1] = -1.0;
    x0[2] = 150.0;
    double perturbation_range[n];
    for (unsigned int i = 0; i < n; i++) {
        perturbation_range[i] = x0[i] * 0.01;
    }

    // ===============================================================================
    // Threads.
    // ===============================================================================
    unsigned int num_threads = 4;  // Number of threads (adjust if necessary)

    // ===============================================================================
    // Reinjection parameters.
    // ===============================================================================
    double yf = 41.2861;
    //double clam = 1.85;

    // ===============================================================================
    // Mapping count parameters.
    // ===============================================================================
    unsigned int map_target_per_thread = 2000;
    unsigned int total_map_target = map_target_per_thread * num_threads;
    double ymap[total_map_target][2];//, zmap[total_map_target][2];

    // ===============================================================================
    // Integrate stationary state.
    // ===============================================================================
    omp_set_num_threads((int)num_threads);
    #pragma omp parallel
    {
        // Get the thread ID and number of threads
        unsigned int thread_id = (unsigned int)omp_get_thread_num();
        unsigned int thread_num = (unsigned int)omp_get_num_threads();

        // ===============================================================================
        // Instanciate parameters for system.
        // ===============================================================================
        lorenz_par p = {params.sigma, params.rho, params.beta};

        // ===============================================================================
        // Initialize random number generator
        // ===============================================================================
        gsl_rng *rng = gsl_rng_alloc(gsl_rng_default);
        gsl_rng_set(rng, (long unsigned int)time(NULL) + (long unsigned int)thread_num);

        // ===============================================================================
        // Apply perturbation to initial conditions
        // ===============================================================================
        double x[n];
        for (size_t i = 0; i < n; ++i) {
            x[i] = x0[i] + gsl_ran_gaussian(rng, perturbation_range[i]);
        }

        // ===============================================================================
        // Preparing integrator.
        // ===============================================================================
        const gsl_odeiv2_step_type *integrator = gsl_odeiv2_step_rkf45;
        gsl_odeiv2_step *s = gsl_odeiv2_step_alloc(integrator, n);
        gsl_odeiv2_control *c = gsl_odeiv2_control_y_new(atol, rtol);
        gsl_odeiv2_evolve *e = gsl_odeiv2_evolve_alloc(n);
        gsl_odeiv2_system sys = {lorenz_gsl, jaclorenz_gsl, n, &p};

        // ===============================================================================
        // Integrate transient.
        // ===============================================================================    
        while (t < t_transient) {
            int status = gsl_odeiv2_evolve_apply(e, c, s, &sys, &t, t_stationary, &h, x);
            
            if (status != GSL_SUCCESS) {
                printf("error: %d\n", status);
            }
        }

        // ===============================================================================
        // Arrays for interpolation and map.
        // ===============================================================================
        long double xfit[n], yfit[n];//, zfit[n];
        double yreg[2];//, zreg[2];
        yreg[0] = yreg[1] = 0.0;
        //zreg[0] = zreg[1] = 0.0;
        double xp = 0;
        long double ciy[n];//, ciz[n];

        // ===============================================================================
        // Fill first three points.
        // ===============================================================================
        for (unsigned int i = 1; i < 3; i++) {
            int status = gsl_odeiv2_evolve_apply(e, c, s, &sys, &t, t_stationary, &h, x);
            xfit[i] = x[0];
            yfit[i] = x[1];
            //zfit[i] = x[2];
            if (status != GSL_SUCCESS) {
                printf("error initializing, return value=%d\n", status);
                break;
            }
        }

        time_t tprev, tnow;
        time(&tprev);
        unsigned int local_map_count = 0;

        while (local_map_count < map_target_per_thread) {
            int status = gsl_odeiv2_evolve_apply(e, c, s, &sys, &t, t_stationary, &h, x);
            
            xfit[0] = xfit[1];
            xfit[1] = xfit[2];
            xfit[2] = x[0];

            yfit[0] = yfit[1];
            yfit[1] = yfit[2];
            yfit[2] = x[1];

            //zfit[0] = zfit[1];
            //zfit[1] = zfit[2];
            //zfit[2] = x[2];

            if (xfit[0] < xp) {
                if (xfit[1] > xp) {
                    quadratic_regression(xfit, yfit, 3, ciy);
                    //quadratic_regression(xfit, zfit, 3, ciz);
                    yreg[1] = (double)(ciy[0] * xp * xp + ciy[1] * xp + ciy[2]);
                    //zreg[1] = (double)(ciz[0] * xp * xp + ciz[1] * xp + ciz[2]);
                    //if (yreg[0] - yf < 20.0 && yreg[0] - yf > -20.0){
                    //    if (yreg[1] - yf < 20.0 && yreg[1] - yf > -20.0){
                    unsigned int index = thread_id * map_target_per_thread + local_map_count;
                    ymap[index][0] = yreg[0];
                    ymap[index][1] = yreg[1];
                            //zmap[index][0] = zreg[0];
                            //zmap[index][1] = zreg[1];
                    local_map_count++;
//                        }
//                    }
                    yreg[0] = yreg[1];
                    //zreg[0] = zreg[1];
                }
            }

            time(&tnow);
            if (difftime(tnow, tprev) > 2) {
                #pragma omp critical
                {
                    printf("Thread: %d, ", thread_id);
                    printf("map count: %d, ", local_map_count);
                    printf("%% completed: %3.2f %%\n", (double)local_map_count * 100.0 / (double)map_target_per_thread);
                }
                tprev = tnow;
            }

            if (status != GSL_SUCCESS) {
                printf("error: %d\n", status);
                break;
            }
        }

        // ===============================================================================
        // Free integrator objects.
        // ===============================================================================
        gsl_odeiv2_evolve_free(e);
        gsl_odeiv2_control_free(c); 
        gsl_odeiv2_step_free(s);
    }

    // ===============================================================================
    // Write results to file.
    // ===============================================================================
    for (unsigned int i = 10; i < total_map_target; i++) {
        fprintf(f, "%5.5f %5.5f\n",
            ymap[i][0],
            ymap[i][1]
        );
    }

    // ===============================================================================
    // Close file.
    // ===============================================================================    
    fclose(f);

    // ===============================================================================
    // Finish statement.
    // ===============================================================================
    printf("Process finished. Results stored in %s\n", write_filename);
    return 0;
}   
