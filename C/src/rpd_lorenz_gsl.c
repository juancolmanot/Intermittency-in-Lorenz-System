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
    double h = params.h;
    double atol, rtol;
    atol = params.atol;
    rtol = params.rtol;
    double t, t_transient, t_stationary;
    t = 0.0;
    t_transient = params.t_transient;
    t_stationary = t_transient + params.t_stationary;

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

    // ===============================================================================
    // Reinjection parameters.
    // ===============================================================================
    double yf = 41.2861;
    double clam = 1.85;
    unsigned int rtarget = 40000, rcount = 0;
    double yreinj[rtarget];

    // ===============================================================================
    // Integrate stationary state.
    // ===============================================================================
    time_t tprev, tnow;
    time(&tprev);
    while (rcount < rtarget) {
        int status = gsl_odeiv2_evolve_apply(e, c, s, &sys, &t, t_stationary, &h, x);
        
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
                if (yreg[1] >= yf - clam && yreg[1] <= yf + clam) {
                    if (yreg[0] < yf - clam || yreg[0] > yf + clam) {
                        yreinj[rcount] = yreg[1];
                        rcount++;
                    }
                }
                yreg[0] = yreg[1];
            }
        }

        time(&tnow);
        if (difftime(tnow, tprev) > 2) {
            printf("reinject count: %d\n", rcount);
            tprev = tnow;
        }

        if (status != GSL_SUCCESS) {
            printf("error: %d\n", status);
            break;
        }
    }

    // ===============================================================================
    // RPD computing.
    // ===============================================================================
    unsigned int nbins = 500;
    double dy = (2 * clam) / (double)(nbins - 1);
    double *bins = calloc(nbins, sizeof(double));
    for (unsigned int i = 0; i < nbins; i++) {
        bins[i] = yf - clam + dy * (double)i;
    }
    double *rpd = calloc(nbins, sizeof(double));
    stats_histogram_double(rpd, bins, yreinj, rtarget, nbins);
    for (unsigned int i = 0; i < nbins; i++) {
        fprintf(f, "%5.5f %5.5f\n", bins[i], rpd[i] / (double)(rtarget));
    }

    // ===============================================================================
    // Free integrator objects.
    // ===============================================================================
    free(bins);
    free(rpd);
    gsl_odeiv2_evolve_free(e);
    gsl_odeiv2_control_free(c);
    gsl_odeiv2_step_free(s);
    fclose(f);
    // ===============================================================================
    // Finish statement.
    // ===============================================================================
    printf("Process finished. Results stored in %s\n", write_filename);
    return 0;
}   