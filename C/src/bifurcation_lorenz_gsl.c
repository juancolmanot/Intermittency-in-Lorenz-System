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
    Parameters5a params;
    load_parameters_from_file(params_file, &params, handler5a);

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
    // Rho parameter.
    // ===============================================================================
    double rho_init, rho_final, delta_rho, rho_i;
    unsigned int n_rho = params.n_rho;
    rho_init = params.rho_init;
    rho_final = params.rho_final;
    delta_rho = (rho_final - rho_init) / (double)n_rho;

    // ===============================================================================
    // System syze and perturbation.
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
    // Poincaré map parameters and fixed point.
    // ===============================================================================
    double yf = 41.2861;
    double dymax = 5.0;
    unsigned int states_per_thread = 1000;
    unsigned int total_states = states_per_thread * num_threads * n_rho;
    double yreinj[total_states][3];
    // ===============================================================================
    // Integrate stationary state.
    // ===============================================================================
    omp_set_num_threads((int)num_threads);
    
    #pragma omp parallel
    {
        // Loop through rho values
        for (unsigned int j = 0; j < n_rho; j++){
            // Get the thread ID and number of threads
            unsigned int thread_id = (unsigned int)omp_get_thread_num();
            unsigned int thread_num = (unsigned int)omp_get_num_threads();

            // Print status
            #pragma omp critical
            {
                printf("Thread: %d, n_rho: %d of %d\n", thread_id, j + 1, n_rho);
                
            }

            // ===============================================================================
            // Initialize random number generator
            // ===============================================================================
            gsl_rng *rng = gsl_rng_alloc(gsl_rng_default);
            gsl_rng_set(rng, (long unsigned int)time(NULL) + (long unsigned int)thread_num);

            // ===============================================================================
            // Restore initial state.
            // ===============================================================================
            x0[0] = 2.0;
            x0[1] = -1.0;
            x0[2] = 150.0;

            // ===============================================================================
            // Instanciate parameters for system.
            // ===============================================================================
            rho_i = rho_init + delta_rho * (double)j;
            lorenz_par p = {params.sigma, rho_i, params.beta};

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
            long double xfit[n], yfit[n];
            double yreg[2];
            yreg[0] = yreg[1] = 0.0;
            double xp = 0;
            long double ci[n];

            // ===============================================================================
            // Fill first three points.
            // ===============================================================================
            for (unsigned int i = 1; i < 3; i++) {
                int status = gsl_odeiv2_evolve_apply(e, c, s, &sys, &t, t_stationary, &h, x);
                xfit[i] = x[0];
                yfit[i] = x[1];
                if (status != GSL_SUCCESS) {
                    printf("error initializing, return value=%d\n", status);
                    break;
                }
            }

            unsigned int local_state_count = 0;
            unsigned long int iteration_count = 0;
            while (local_state_count < states_per_thread) {
                iteration_count++;
                int status = gsl_odeiv2_evolve_apply(e, c, s, &sys, &t, t_stationary, &h, x);
                
                xfit[0] = xfit[1];
                xfit[1] = xfit[2];
                xfit[2] = x[0];

                yfit[0] = yfit[1];
                yfit[1] = yfit[2];
                yfit[2] = x[1];

                if (xfit[0] < xp && xfit[1] > xp) {
                    quadratic_regression(xfit, yfit, 3, ci);
                    yreg[1] = (double)(ci[0] * xp * xp + ci[1] * xp + ci[2]);
                    if (yreg[1] >= yf - dymax && yreg[1] <= yf + dymax) {
                        if (yreg[0] >= yf - dymax && yreg[0] <= yf + dymax) {
                            unsigned int index = thread_id * (states_per_thread * n_rho) + j * states_per_thread + local_state_count;
                            yreinj[index][0] = rho_i;
                            yreinj[index][1] = yreg[0];
                            yreinj[index][2] = yreg[1];
                            local_state_count++;
                        }
                    }
                    yreg[0] = yreg[1];
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
    }

    // ===============================================================================
    // Write results to file.
    // ===============================================================================
    for (unsigned int i = 0; i < total_states; i++) {
        if (fabs(yreinj[i][1] - yreinj[i][2]) < 1e-4){
            fprintf(f, "%12.5E %12.5E %12.5E\n",
            yreinj[i][0],
            yreinj[i][1],
            yreinj[i][2]
            );
        }
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
