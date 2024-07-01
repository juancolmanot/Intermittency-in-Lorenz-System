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
    // Reset t and h and integrate stationary state.
    // ===============================================================================
    t = 0.0;
    h = params.h;
    for (unsigned int i = 1; i <= n_steps; i++) {
        double ti = i * t_stationary / (double)n_steps;
        int status = gsl_odeiv2_driver_apply(d, &t, ti, x);
        fprintf(f, "%5.6f %5.6f %5.6f %5.6f\n", t, x[0], x[1], x[2]);
        if (status != GSL_SUCCESS) {
            printf("error: %d\n", status);
            break;
        }
    }

    // ===============================================================================
    // Free integrator objects.
    // ===============================================================================
    fclose(f);
    // ===============================================================================
    // Finish statement.
    // ===============================================================================
    printf("Process finished. Results stored in %s\n", write_filename);
    return 0;
}   
