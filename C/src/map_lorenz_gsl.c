#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
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
    // Instanciate parameters for system.
    // ===============================================================================
    size_t n = 3;
    lorenz_par p = {params.sigma, params.rho, params.beta};

    // ===============================================================================
    // Initial state.
    // ===============================================================================
    double x[n];
    x[0] = 2.0;
    x[1] = -1.0;
    x[2] = 150.0;

    // ===============================================================================
    // Integration parameters.
    // ===============================================================================
    unsigned int n_steps = params.n_steps;
    double h = params.h;
    double atol, rtol;
    atol = params.atol;
    rtol = params.rtol;
    double t, t_transient, t_stationary;
    t = 0.0;
    t_transient = params.t_transient;
    t_stationary = params.t_stationary;

    // ===============================================================================
    // Preparing integrator.
    // ===============================================================================
    gsl_odeiv2_system sys = {lorenz_gsl, jaclorenz_gsl, n, &p};
    gsl_odeiv2_driver *d = gsl_odeiv2_driver_alloc_y_new(&sys, gsl_odeiv2_step_rkf45, h, atol, rtol);

    // ===============================================================================
    // Integrate transient.
    // ===============================================================================    
    int status = gsl_odeiv2_driver_apply(d, &t, t_transient, x);

    if (status != GSL_SUCCESS) {
        printf("error: %d\n", status);
    }

    // ===============================================================================
    // Arrays for interpolation and map.
    // ===============================================================================
    long double xfit[n], yfit[n];
    double yreg;
    unsigned int mtarget = 100000, mcount = 0;
    unsigned int icount = 0;
    double xp = 0;
    long double ci[n];

    // ===============================================================================
    // Fill first three points.
    // ===============================================================================
    t = 0.0;
    h = params.h;
    d = gsl_odeiv2_driver_alloc_y_new(&sys, gsl_odeiv2_step_rkf45, h, atol, rtol);
    for (unsigned int i = 0; i < 3; i++) {
        double ti = i * h;
        int status = gsl_odeiv2_driver_apply(d, &t, ti, x);
        xfit[i] = x[0];
        yfit[i] = x[1];

        if (status != GSL_SUCCESS) {
            printf("error: %d\n", status);
            break;
        }
    }

    // ===============================================================================
    // Reset t and h and integrate stationary state.
    // ===============================================================================
    t = 0.0;
    h = params.h;
    time_t tprev, tnow;
    time(&tprev);
    d = gsl_odeiv2_driver_alloc_y_new(&sys, gsl_odeiv2_step_rkf45, h, atol, rtol);
    while (mcount < mtarget) {
        icount++;
        double ti = icount * t_stationary / (double)n_steps;
        int status = gsl_odeiv2_driver_apply(d, &t, ti, x);
        
        xfit[0] = xfit[1];
        xfit[1] = xfit[2];
        xfit[2] = x[0];

        yfit[0] = yfit[1];
        yfit[1] = yfit[2];
        yfit[2] = x[1];

        if (xfit[0] < xp) {
            if (xfit[1] > xp) {
                quadratic_regression(xfit, yfit, 3, ci);
                yreg = (double)(ci[0] * xp * xp + ci[1] * xp + ci[2]);
                fprintf(f, "%d %5.5f\n", mcount, yreg);
                mcount++;
            }
        }

        time(&tnow);
        if (difftime(tnow, tprev) > 2) {
            printf("map count: %d\n", mcount);
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
    fclose(f);
    gsl_odeiv2_driver_free(d);
    // ===============================================================================
    // Finish statement.
    // ===============================================================================
    printf("Process finished. Results stored in %s\n", write_filename);
    return 0;
}   
