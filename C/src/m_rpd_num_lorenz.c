#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <stdbool.h>
#include <complex.h>
#include "../include/lorenz.h"
#include "../../../common/include/file_handle.h"
#include "../../../common/include/progress_handle.h"
#include "../../../common/include/parameters.h"
#include "../../../common/include/ini.h"
#include "../../../common/include/stats.h"
#include "../../../common/include/linear_algebra.h"
#include "../../../common/include/dynamical_systems.h"
#include "../../../common/include/gsl_utilities.h"
#include "../../../common/include/multipurpose_routines.h"
#include "../../../common/include/intermittency.h"

int main(int argc, char *argv[]) {

    // ===============================================================================
    // Check for right passed arguments.
    // ===============================================================================
    if (argc != 6){
        printf("Wrong amount of arguments, you passed %d, but 5 are required\n.", argc);
        printf("Usage: ./run-script.sh script_to_run m_write_file");
        printf(" rpd_write_files reinjection_data.dat params_file.ini\n");
        exit(EXIT_FAILURE);
    }

    // ===============================================================================
    // Load arguments into variables.
    // ===============================================================================
    const char *m_write_file = argv[2];
    const char *rpd_write_file = argv[3];
    const char *reinjection_data = argv[4];
    const char *params_file = argv[5];

    // ===============================================================================
    // Parameters.
    // ===============================================================================
    Parameters8 params;
    load_parameters_from_file(params_file, &params, handler8);

    // ===============================================================================
    // Store data from reinjection points file.
    // ===============================================================================
    unsigned int rows = 0, cols = 0;
    long double **data;
    read_data_file_unsigned(reinjection_data, &data, &rows, &cols);

    // ===============================================================================
    // Load reinjection points.
    // ===============================================================================
    long double *xreinj = calloc(rows, sizeof(long double));
    for (unsigned int i = 0; i < rows; i++) {
        xreinj[i] = data[i][3];
    }
    // Sort x array.
    quicksort_long_unsigned(xreinj, rows);

    // ===============================================================================
    // Compute numerical M(x).
    // ===============================================================================
    long double *M_tot = calloc(rows, sizeof(long double));
    compute_m_function(xreinj, M_tot, rows);

    // ===============================================================================
    // Compute numerical RPD(x).
    // ===============================================================================
    unsigned int nbins_total = 0;
    for (int i = 0; i < params.n_regions; i++) {
        nbins_total += params.region_bins[i];
    }
    long double *rpd_total = calloc(nbins_total, sizeof(long double));
    long double *bins_total = calloc(nbins_total, sizeof(long double));
    compute_rpd_numerical(xreinj, rows, bins_total, rpd_total, nbins_total, 1.0);

    // ===============================================================================
    // Open files to complete M(x) and RPD(x) and erase their content.
    // ===============================================================================
    FILE *fm = fopen(m_write_file, "w"), *frpd = fopen(rpd_write_file, "w");
    for (unsigned int i = 1; i < rows; i++) {
        fprintf(fm, "%12.5LE %12.5LE\n", xreinj[i], M_tot[i]);
    }
    for (unsigned int i = 1; i < nbins_total; i++) {
        fprintf(frpd, "%12.5LE %12.5LE\n", bins_total[i], rpd_total[i]);
    }
    fclose(fm);
    fclose(frpd);

    // ===============================================================================
    // Free alocated memory.
    // ===============================================================================
    for (unsigned int i = 0; i < rows; i++) {
        free(data[i]);
    }
    free(data);
    free(xreinj);
    free(M_tot);
    free(rpd_total);
    free(bins_total);

    // ===============================================================================
    // Finsh statement.
    // ===============================================================================
    printf("Results stored in %s & %s\n", m_write_file, rpd_write_file);
    printf("Process finished.\n");
    return 0;
}