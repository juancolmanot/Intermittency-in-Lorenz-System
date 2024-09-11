#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <omp.h>

#define BIN_COUNT 100 // Define number of bins of domain.


// Lorenz system struct.
typedef struct {
    double sigma;
    double rho;
    double beta;
} lorenz_par;


int main(void) {

    // ===============================================================================
    // Thread config
    // ===============================================================================
    int num_threads = omp_get_max_threads();

    // GSL random number generator
    gsl_rng_env_setup();
    gsl_rng *rng = gsl_rng_alloc(gsl_rng_default);
    gsl_rng_set(rng, time(NULL));
    
    // ===============================================================================
    // System parameters
    // ===============================================================================
    double sigma = 10.0;
    double rho = 166.07;
    double beta = 8.0 / 3.0;
    double delta = 1e-4;

    // ===============================================================================
    // Instanciate parameters for system.
    // ===============================================================================
    size_t n = 3;
    lorenz_par p = {sigma, rho, beta};

    // ===============================================================================
    // Maximum and minimum values.
    // ===============================================================================
    double xmin = -5.83170e1, xmax = 5.95057e1;
    double ymin = -1.37485e2, ymax = 1.42208e2;
    double zmin = 2.86333e2, zmax = 2.93883e2;

    // ===============================================================================
    // Binning parameters
    // ===============================================================================
    double x_bin_size = (xmax - xmin) / BIN_COUNT;
    double y_bin_size = (ymax - ymin) / BIN_COUNT;
    double z_bin_size = (zmax - zmin) / BIN_COUNT;

    // ===============================================================================
    // Integration parameters.
    // ===============================================================================
    unsigned int n_steps = 1000000;
    double t_stationary = 100000.0;
    double t_transient = 1000.0;
    double h0 = 1e.6;
    double atol = 1e-12;
    double rtol = 1e-14;

    // ===============================================================================
    // Bin array.
    // ===============================================================================
    int bins[BIN_COUNT][BIN_COUNT][BIN_COUNT] = {0};

    #pragma omp parallel num_threads(num_threads)
    {

        // Private random number generator per thread
        gsl_rng *thread_rng = gsl_rng_alloc(gsl_rng_default);
        gsl_rng_set(thread_rng, time(NULL) + omp_get_thread_num());

        // ===============================================================================
        // Initial state.
        // ===============================================================================
        double x[n];
        x[0] = 2.0 + gsl_rng_uniform(thread_rng) * delta;
        x[1] = -1.0 + gsl_rng_uniform(thread_rng) * delta;
        x[2] = 150.0 + gsl_rng_uniform(thread_rng) * delta;
        double h = h0;

        // ===============================================================================
        // Preparing integrator.
        // ===============================================================================
        gsl_odeiv2_system sys = {lorenz_gsl, jaclorenz_gsl, n, &p};
        gsl_odeiv2_driver *d = gsl_odeiv2_driver_alloc_y_new(&sys, gsl_odeiv2_step_rkf45, h, atol, rtol);

        // ===============================================================================
        // Integrate transient.
        // ===============================================================================    
        double t = 0.0
        int status = gsl_odeiv2_driver_apply(d, &t, t_transient, x);

        if (status != GSL_SUCCESS) {
            printf("error: %d\n", status);
        }


        // ===============================================================================
        // Reset t and h and integrate stationary state.
        // ===============================================================================
        t = 0.0;
        h = h0;

        for (unsigned int i = 1; i <= n_steps; i++) {
            double ti = (double)i * t_stationary / (double)n_steps;
            int status = gsl_odeiv2_driver_apply(d, &t, ti, x);

            int x_bin = (int)((x[0] - xmin) / x_bin_size);
            int y_bin = (int)((x[1] - ymin) / y_bin_size);
            int z_bin = (int)((x[2] - zmin) / z_bin_size);

            if (x_bin >= 0 && x_bin < BIN_COUNT &&
                y_bin >= 0 && y_bin < BIN_COUNT &&
                z_bin >= 0 && z_bin < BIN_COUNT) {
                #pragma omp atomic
                bins[x_bin][y_bin][z_bin]++;
            }

            if (status != GSL_SUCCESS) {
                printf("error: %d\n", status);
                break;
            }
        }

        gsl_odeiv2_driver_free(d);
        gsl_rng_free(thread_rng);
        free(x_state);
        free(y_state);
        free(z_state);
    }

    // ===============================================================================
    // Write to file.
    // ===============================================================================
    FILE *f = fopen("../datafiles/lorenz_invariant_measure_1.dat", "w");
    if (f == NULL) {
        fprintf(stderr, "Error: Could not open file for writing.\n");
        return EXIT_FAILURE;
    }

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
    // Free objects.
    // ===============================================================================
    fclose(f);

    return 0;
}
