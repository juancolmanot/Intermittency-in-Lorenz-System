#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <omp.h>
#include <float.h>
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
    // Common min and max array.
    // ===============================================================================
    double common_min_values[n][4], common_max_values[n][4];
    double min_minimorum[n], max_maximorum[n];
    // Fill max and min arrays.
    for (size_t i = 0; i < n; i++) {
        min_minimorum[i] = DBL_MAX;
        max_maximorum[i] = -DBL_MAX;
        for (size_t j = 0; j < num_threads; j++) {
            common_min_values[i][j] = DBL_MAX;
            common_max_values[i][j] = -DBL_MAX;    
        }
    }

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
        gsl_rng_set(rng, (long unsigned int)omp_get_wtime() + (long unsigned int)thread_num);

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

        double tprev = omp_get_wtime();

        while (t < t_stationary) {
            int status = gsl_odeiv2_evolve_apply(e, c, s, &sys, &t, t_stationary, &h, x);

            /*
            double tnow = omp_get_wtime();
            if (tnow - tprev > 2.0) {
                #pragma omp critical
                {
                    printf("Thread: [%d] - ", thread_id);
                    printf(" progress: [%4.2f%%]\n", t * 100 / t_stationary);
                }
                tprev = tnow;
            }
            */
            double tnow = omp_get_wtime();
            if (tnow - tprev > 2.0 && thread_id == 0) { // Only the master thread prints progress
                double progress = t * 100 / t_stationary;
                printf("Progress: [%4.2f%%]\n", progress);
                fflush(stdout); // Ensure the output is flushed immediately
                tprev = tnow;
            }

            // Set max and min values.
            for (size_t i = 0; i < n; i++) {
                if (x[i] < common_min_values[i][thread_id]) {
                    common_min_values[i][thread_id] = x[i];
                }
                else if (x[i] > common_max_values[i][thread_id]) {
                    common_max_values[i][thread_id] = x[i];
                }
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
    // Determin overall max and min.
    // ===============================================================================
    for (size_t i = 0; i < num_threads; i++) {
        for (size_t j = 0; j < n; j++) {
            if (common_min_values[j][i] < min_minimorum[j]) {
                min_minimorum[j] = common_min_values[j][i];
            }
            if (common_max_values[j][i] > max_maximorum[j]) {
                max_maximorum[j] = common_max_values[j][i];
            }
        }
    }

    // ===============================================================================
    // Write results to file.
    // ===============================================================================
    for (size_t i = 0; i < n; i++) {
        fprintf(f, "%12.5E %12.5E\n",
            min_minimorum[i],
            max_maximorum[i]
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
