#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <stdbool.h>
#include "../include/lorenz.h"
#include "../include/intermittency_lorenz.h"
#include "../../../common/include/file_handle.h"
#include "../../../common/include/parameters.h"
#include "../../../common/include/odeintegration.h"
#include "../../../common/include/ini.h"
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
    int n = 3;
    long double p[n];
    p[0] = params.sigma;
    p[1] = params.rho;
    p[2] = params.beta;

    // ===============================================================================
    // Initial state.
    // ===============================================================================
    long double x0[n];
    x0[0] = 2.0;
    x0[1] = -1.0;
    x0[2] = 150.0;

    // ===============================================================================
    // Integration parameters.
    // ===============================================================================
    long double h0 = params.h;
    long double atol, rtol;
    atol = params.atol;
    rtol = params.rtol;
    long double t0, t_transient;
    t0 = 0.0;
    t_transient = params.t_transient;

    RKF45State state;
    state.t = t0;
    state.h = h0;
    state.x = x0;
    state.n = n;

    // ===============================================================================
    // Integration of transient.
    // ===============================================================================
    rkf45_integrate(&state, lorenz, p, t_transient, atol, rtol);

    // ===============================================================================
    // Arrays for interpolation and mapping.
    // ===============================================================================
    long double xfit[n], yfit[n], zfit[n];
    long double xp = 0;
    long double ciy[n], ciz[n];
    long unsigned int mcount = 0;
    long double yreg[2];
    yreg[0] = yreg[1] = 0.0;

    // ===============================================================================
    // Arrays and values for reinjection.
    // ===============================================================================
    unsigned int rtarget = 10000, rcount = 0;
    long double yf = 41.2861;
    long double clam = 1.85;
    long double *xreinj = calloc(rtarget, sizeof(long double));
    time_t tprev, tnow;
    time(&tprev);

    // ===============================================================================
    // Integration.
    // ===============================================================================
    while (rcount < rtarget) {
        rkf45_step(&state, lorenz, p);
        xfit[0] = xfit[1];
        xfit[1] = xfit[2];
        xfit[2] = state.x[0];

        yfit[0] = yfit[1];
        yfit[1] = yfit[2];
        yfit[2] = state.x[1];

        zfit[0] = zfit[1];
        zfit[1] = zfit[2];
        zfit[2] = state.x[2];

        if (xfit[0] < xp) {
        	if (xfit[1] > xp) {
                mcount++;
                time(&tnow);
                if (difftime(tnow, tprev) > 2) {
                    printf("map count: %ld, reinj count: %d\n", mcount, rcount);
                    tprev = tnow;
                }
        		quadratic_regression(xfit, yfit, 3, ciy);
        		quadratic_regression(xfit, zfit, 3, ciz);
        		yreg[1] = ciy[0] * xp * xp + ciy[1] * xp + ciy[2];
                if (reinjection_1d(yreg[0], yreg[1], yf, clam) == 1) {
                    xreinj[rcount] = yreg[1];
                    rcount++;
                }
                yreg[0] = yreg[1];
        	}
        }
    }

    // ===============================================================================
    // RPD computing.
    // ===============================================================================
    unsigned int nbins = 500;
    long double dx = (2 * clam) / (long double)(nbins - 1);
    long double *bins = calloc(nbins, sizeof(long double));
    for (unsigned int i = 0; i < nbins; i++) {
        bins[i] = yf - clam + dx * (long double)i;
    }
    long double *rpd = calloc(nbins, sizeof(long double));
    stats_histogram(rpd, bins, xreinj, rtarget, nbins);
    for (unsigned int i = 0; i < nbins; i++) {
        fprintf(f, "%5.5Lf %5.5Lf\n", bins[i], rpd[i]);
    }

    free(bins);
    free(rpd);
    free(xreinj);
    fclose(f);
    // ===============================================================================
    // Finish statement.
    // ===============================================================================
    printf("Process finished. Results stored in %s\n", write_filename);
    return 0;
}