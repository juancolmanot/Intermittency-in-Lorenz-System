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

#define BIN_COUNT 100 // Define number of bins of domain.

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
    unsigned int n_steps = (unsigned int)params.max_steps;
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
    // State array.
    // ===============================================================================
    double *x_state, *y_state, *z_state;
    x_state = calloc(n_steps, sizeof(double));
    y_state = calloc(n_steps, sizeof(double));
    z_state = calloc(n_steps, sizeof(double));

    // ===============================================================================
    // Reset t and h and integrate stationary state.
    // ===============================================================================
    t = 0.0;
    h = params.h;
    time_t tprev, tnow;
    time(&tprev);
    for (unsigned int i = 1; i <= n_steps; i++) {
        double ti = (double)i * t_stationary / (double)n_steps;
        int status = gsl_odeiv2_driver_apply(d, &t, ti, x);
        x_state[i] = x[0];
        y_state[i] = x[1];
        z_state[i] = x[2];

        time(&tnow);
        if (difftime(tnow, tprev) > 2) {
            printf("progress: %d/%d [%5.5f%%]\n",
                i, n_steps, ((double)i / (double)n_steps));
            tprev = tnow;
        }

        if (status != GSL_SUCCESS) {
            printf("error: %d\n", status);
            break;
        }
    }

    // ===============================================================================
    // Find minimum and maximum values.
    // ===============================================================================
    double xmin, xmax, ymin, ymax, zmin, zmax;
    xmin = la_min_d(x_state, n_steps);
    ymin = la_min_d(y_state, n_steps);
    zmin = la_min_d(z_state, n_steps);
    xmax = la_max_d(x_state, n_steps);
    ymax = la_max_d(y_state, n_steps);
    zmax = la_max_d(z_state, n_steps);

    // ===============================================================================
    // Bins array.
    // ===============================================================================
    int bins[BIN_COUNT][BIN_COUNT][BIN_COUNT] = {0};
    double x_bin_size = (xmax - xmin) / BIN_COUNT;
    double y_bin_size = (ymax - ymin) / BIN_COUNT;
    double z_bin_size = (zmax - zmin) / BIN_COUNT;    

    // ===============================================================================
    // Fill bins array.
    // ===============================================================================
    for (unsigned int i = 0; i < n_steps;) {
        int x_bin = (int)((x_state[i] - xmin) / x_bin_size);
        int y_bin = (int)((y_state[i] - ymin) / y_bin_size);
        int z_bin = (int)((z_state[i] - zmin) / z_bin_size);
        
        if (x_bin >= 0 && x_bin < BIN_COUNT &&
            y_bin >= 0 && y_bin < BIN_COUNT &&
            z_bin >= 0 && z_bin < BIN_COUNT){
            bins[x_bin][y_bin][z_bin]++;
        }
    }

    // ===============================================================================
    // Write to file.
    // ===============================================================================
    for (int i = 0; i < BIN_COUNT; i++) {
        for (int j = 0; j < BIN_COUNT; j++) {
            for (int k = 0; k < BIN_COUNT; k++) {
                if (bins[i][j][k] > 0) {
                    double x_center = xmin + (i + 0.5) * x_bin_size;
                    double y_center = ymin + (j + 0.5) * y_bin_size;
                    double z_center = zmin + (k + 0.5) * z_bin_size;
                    fprintf(f, "%12.5E %12.5E %12.5E %d\n",
                        x_center,
                        y_center,
                        z_center,
                        bins[i][j][k]
                    );
                }
            }
        }
    }

    // ===============================================================================
    // Free integrator objects.
    // ===============================================================================
    gsl_odeiv2_driver_free(d);
    fclose(f);
    // ===============================================================================
    // Finish statement.
    // ===============================================================================
    printf("Process finished. Results stored in %s\n", write_filename);
    return 0;
}   
