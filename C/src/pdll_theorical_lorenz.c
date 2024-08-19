#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <math.h>
#include <time.h>
#include "../include/lorenz.h"
#include "../../../common/include/file_handle.h"
#include "../../../common/include/progress_handle.h"
#include "../../../common/include/parameters.h"
#include "../../../common/include/ini.h"
#include "../../../common/include/stats.h"
#include "../../../common/include/linear_algebra.h"
#include "../../../common/include/gsl_utilities.h"
#include "../../../common/include/dynamical_systems.h"
#include "../../../common/include/intermittency.h"

int main(int argc, char *argv[]) {

    // ===============================================================================
    // Check for right passed arguments.
    // ===============================================================================
    if (argc != 5){
        printf("Wrong amount of arguments, you passed %d, but 5 are required\n.", argc);
        printf("Usage: ./run-script.sh script_to_run pdll_write_file rpd_read_file");
        printf(" coeffs_read file.\n");
        exit(EXIT_FAILURE);
    }

    // ===============================================================================
    // Load arguments into variables.
    // ===============================================================================
    const char *pdll_write_file = argv[2];
    const char *rpd_read_file = argv[3];
    const char *coeffs_read_file = argv[4];

    // ===============================================================================
    // File to write to
    // ===============================================================================
    FILE *f = open_file(pdll_write_file);

    // ===============================================================================
    // Store data from rpd file.
    // ===============================================================================
    unsigned int rpd_rows = 0, rpd_cols = 0;
    double **rpd_data;
    read_data_file_unsigned_double(rpd_read_file, &rpd_data, &rpd_rows, &rpd_cols);
    double x[rpd_rows], rpd[rpd_rows];
    for (unsigned int i = 0; i < rpd_rows; i++) {
        x[i] = rpd_data[i][2];
        rpd[i] = rpd_data[i][3];
    }

    // ===============================================================================
    // Store data from coeffs fit file.
    // ===============================================================================
    unsigned int coeffs_rows = 0, coeffs_cols = 0;
    double **coeffs_data;
    read_data_file_unsigned_double(coeffs_read_file, &coeffs_data, &coeffs_rows, &coeffs_cols);
    double a, b, c;
    a = coeffs_data[0][0];
    b = coeffs_data[0][1];
    c = coeffs_data[0][2];

    // ===============================================================================
    // Size of laminar region.
    // ===============================================================================
    double yf = 41.2861;
    double clam = 1.85;

    // ===============================================================================
    // Compute pdll.
    // ===============================================================================
    unsigned int nlams = 20;
    double pdll[nlams], ll[nlams], local_map = 0.0, x_i = 0.0;
    for (unsigned int i = 0; i < nlams; i++) {
        x_i = type_I_xl_2((double)i, a, b, c, clam);
        local_map = a * x[i] * x[i] + (b - 1) * x[i] + c;
        pdll[i] = rpd[i] * fabs(local_map);
        fprintf(f, "%5.5f %5.5f %5.5f %5.5f %5.5f %5.5f %5.5f\n",
            x[i],
            ll[i],
            pdll[i],
            rpd[i],
            local_map,
            x_i,
            atan(x_i)
        );
    }

    // ===============================================================================
    // Free memory.
    // ===============================================================================
    for (unsigned int i = 0; i < rpd_rows; i++) {
        free(rpd_data[i]);
    }
    for (unsigned int i = 0; i < coeffs_rows; i++) {
        free(coeffs_data[i]);
    }
    free(rpd_data);
    free(coeffs_data);
    fclose(f);
    
    // ===============================================================================
    // Finish statement.
    // ===============================================================================
    printf("Process finished. Results stored in %s.\n", pdll_write_file);
    return 0;
}