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
    Parameters6 params;
    load_parameters_from_file(params_file, &params, handler6);

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
    double perturbation_range = 1.0;

    // ===============================================================================
    // Threads.
    // ===============================================================================
    unsigned int num_threads = 4;  // Number of threads (adjust if necessary)

    // ===============================================================================
    // Reinjection parameters.
    // ===============================================================================
    double yf = 41.2861;
    double clam = 1.85;
    double ymin, ymax;
    ymin = params.xmin;
    ymax = params.xmax;
    unsigned int rtarget_per_thread = params.rtarget;
    unsigned int total_rtarget = rtarget_per_thread * num_threads;
    double yreinj[total_rtarget][2];

    // ===============================================================================
    // Integrate stationary state.
    // ===============================================================================
    omp_set_num_threads((int)num_threads);
    int stop_signal = 0;
    #pragma omp parallel shared(stop_signal)
    {
        // Get the thread ID and number of threads
        unsigned int thread_id = (unsigned int)omp_get_thread_num();
        unsigned int thread_num = (unsigned int)omp_get_num_threads();

        // ===============================================================================
        // Error and overflow contention parameters.
        // ===============================================================================
        double ymax_allow = 300.0;
        long unsigned int steps = 0, max_steps = params.max_steps;
        double max_runtime = 40.0 * 60.0; //(45 minutes)
        time_t start_time;
        time(&start_time);

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
            x[i] = x0[i] + gsl_ran_gaussian(rng, perturbation_range);
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
                printf("GSL error: %s\n", gsl_strerror(status));
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
                printf("GSL error: %s\n", gsl_strerror(status));
                break;
            }
        }

        time_t tprev, tnow;
        time(&tprev);
        unsigned int local_rcount = 0;
        time_t current_time;
        time(&current_time);
        double elapsed_time = difftime(current_time, start_time);
        while (local_rcount < rtarget_per_thread) {
            int status = gsl_odeiv2_evolve_apply(e, c, s, &sys, &t, t_stationary, &h, x);
            
            time_t current_time;
            time(&current_time);
            elapsed_time = difftime(current_time, start_time);
            #pragma omp atomic
            steps++;

            xfit[0] = xfit[1];
            xfit[1] = xfit[2];
            xfit[2] = x[0];

            yfit[0] = yfit[1];
            yfit[1] = yfit[2];
            yfit[2] = x[1];

            if (xfit[0] < xp) {
                if (xfit[1] > xp) {
                    quadratic_regression(xfit, yfit, 3, ci);
                    yreg[1] = (double)(ci[0] * xp * xp + ci[1] * xp + ci[2]);
                    if (yreg[0] > ymin && yreg[0] < ymax){
                        if (yreg[1] >= yf - clam && yreg[1] <= yf + clam) {
                            if (yreg[0] < yf - clam || yreg[0] > yf + clam) {
                                if (fabs(yreg[0]) < ymax_allow && fabs(yreg[1]) < ymax_allow) {
                                    unsigned int index = thread_id * rtarget_per_thread + local_rcount;
                                    yreinj[index][0] = yreg[0];
                                    yreinj[index][1] = yreg[1];
                                    local_rcount++;
                                }
                                else {
                                    printf("Unrealistic value encountered: yreg[0]=%f, yreg[1]=%f\n", yreg[0], yreg[1]);
                                    stop_signal = 1;
                                }
                            }
                        }
                    }
                    yreg[0] = yreg[1];
                }
            }

            time(&tnow);
            if (difftime(tnow, tprev) > 2) {
                #pragma omp critical
                {
                    printf("Thread: %d, ", thread_id);
                    printf("reinject count: %d, ", local_rcount);
                    printf("%% completed: %3.2f ,%% ", (double)local_rcount * 100.0 / (double)rtarget_per_thread);
                    printf("%% steps realized: %3.2f%%, ", (double)steps * 100.0 / (double)max_steps);
                    printf("%% time elapsed: %3.2f%%\n", elapsed_time  * 100.0 / max_runtime);
                }
                tprev = tnow;
            }

            if (status != GSL_SUCCESS) {
                printf("error: %d\n", status);
                printf("GSL error: %s\n", gsl_strerror(status));
                stop_signal = 1;
            }

            if (elapsed_time > max_runtime) {
                #pragma omp critical
                {
                    printf("Maximum runtime exceeded. Exiting.\n");
                    stop_signal = 1;
                }
            }

            if (steps >= max_steps) {
                #pragma omp critical
                {
                    printf("Maximum steps exceeded. Exiting.\n");
                    stop_signal = 1;
                }
            }

            if (stop_signal) {
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
    for (unsigned int i = 0; i < total_rtarget; i++) {
        fprintf(f, "%12.5E %12.5E %12.5E %12.5E\n",
            yreinj[i][0],
            yreinj[i][1],
            yreinj[i][0] - yf,
            yreinj[i][1] - yf
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