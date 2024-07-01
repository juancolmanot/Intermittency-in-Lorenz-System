#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include "../include/lorenz.h"
#include "../../../common/include/file_handle.h"
#include "../../../common/include/parameters.h"
#include "../../../common/include/odeintegration.h"
#include "../../../common/include/ini.h"

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
    long double t0, t1;
    t0 = 0.0;
    t1 = 5.0;

    RKF45State state;
    state.t = t0;
    state.h = h0;
    state.x = x0;
    state.n = n;

    // ===============================================================================
    // Integration.
    // ===============================================================================
    while (state.t < t1) {
        if (state.t + state.h > t1) {
            state.h = t1 - state.t;
        }
        rkf45_step_adap(&state, lorenz, p, atol, rtol);
        fprintf(f, "%5.8Lf %5.8Lf %5.8Lf %5.8Lf\n",
            state.t,
            state.x[0],
            state.x[1],
            state.x[2]
        );
    }

    fclose(f);
    // ===============================================================================
    // Finish statement.
    // ===============================================================================
    printf("Process finished. Results stored in %s\n", write_filename);
    return 0;
}