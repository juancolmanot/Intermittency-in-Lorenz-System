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
    if (argc != 5){
        printf("Wrong amount of arguments, you passed %d, but 5 are required\n.", argc);
        printf("Usage: ./run.sh script_to_run rpd_write_file");
        printf(" m_read_file.dat params_file.ini\n");
        exit(EXIT_FAILURE);
    }

    // ===============================================================================
    // Load arguments into variables.
    // ============================================================ ===================
    const char *rpd_write_file = argv[2];
    const char *m_read_file = argv[3];
    const char *params_file = argv[4];

    // ===============================================================================
    // Parameters.
    // ===============================================================================
    Parameters7a params;
    load_parameters_from_file(params_file, &params, handler7a);

    // ===============================================================================
    // Store data from m numeric files.
    // ===============================================================================
    unsigned int rows_m = 0, cols_m = 0;
    long double **data_m;
    read_data_file_unsigned(m_read_file, &data_m, &rows_m, &cols_m);


    // ===============================================================================
    // Get arrays of x and M(x).
    // ===============================================================================
    long double *x_tot = calloc(rows_m, sizeof(long double));
    long double *M_tot = calloc(rows_m, sizeof(long double));
    for (size_t i = 0; i < rows_m; i++) {
        x_tot[i] = data_m[i][0];
        M_tot[i] = data_m[i][1];
    }
    
    // ===============================================================================
    // Open files to complete M(x) and RPD(x) and erase their content.
    // ===============================================================================
    char file_rpd_all[256];
    sprintf(file_rpd_all, "%s_rpd_all.dat", rpd_write_file);
    FILE *frpd_all = fopen(file_rpd_all, "w");
    fclose(frpd_all);

    for (int i = 0; i < params.n_regions; i++){
        
        // ===============================================================================
        // Get region's arrays for x and M(x).
        // ===============================================================================
        long double *x_region = calloc(100, sizeof(long double));
        unsigned int n_region = 0;
        long double xmin = params.xmins[i], xmax = params.xmaxs[i];
        get_region_array(x_tot, xmin, xmax, rows_m, &x_region, &n_region);
        
        // ===============================================================================
        // Get indeces for xmin and xmax in region.
        // ===============================================================================
        unsigned int xminloc = 0, xmaxloc = 0;
        xminloc = get_value_loc(x_region, x_region[0], n_region);
        xmaxloc = get_value_loc(x_region, x_region[n_region - 1], n_region);
        
        // Load portion of M(x) inside this M_region:
        long double *M_region = calloc(n_region, sizeof(long double));
        for (unsigned int j = xminloc; j < xmaxloc + 1; j++) {
            M_region[j - xminloc] = M_tot[j];
        }

        // ===============================================================================
        // Compute slope m and constant b from M(x) function.
        // ===============================================================================
        long double mi = 0.0, bi = 0.0;
        linear_regression(x_region, M_region, n_region, &mi, &bi);
        
        // ===============================================================================
        // Compute Theoretical RPD(x).
        // ===============================================================================
        long double alpha = 0.0;
        alpha = (2.0 * mi - 1.0) / (1.0 - mi);
        printf("subregion: %d - m: %Lf - alpha: %Lf\n", i, mi, alpha);
        long double *rpd_theo = calloc(n_region, sizeof(long double));
        long double xci = params.xc[i];
        // Normalize rpd.
        long double int_rpd = 0.0;
        long double lbound = x_region[0], ubound = x_region[n_region - 1];
        int_rpd = rpd_theoretical_integral(alpha, xci, lbound, ubound);
        int_rpd = fabsl(int_rpd);
        long double bnorm = (long double)(params.wi * n_region) / ((long double)(rows_m) * int_rpd);
        // Compute rpd theoretical
        for (unsigned int j = 0; j < n_region; j++) {
            rpd_theo[j] = bnorm * powl(fabsl(x_region[j] - xci), alpha);
        }

        // ===============================================================================
        // Write results to files.
        // ===============================================================================
        char file_rpd[256];
        sprintf(file_rpd, "%s_rpd_%d.dat", rpd_write_file, i);
        FILE *frpd = fopen(file_rpd, "w");
        for (unsigned int j = 1; j < n_region; j++) {
            fprintf(frpd, "%12.7LE %12.7LE\n",
                x_region[j],
                rpd_theo[j]
            );
        }

        printf("Results stored in %s.\n", file_rpd);
        printf("===================================================================================\n");
        fclose(frpd);
        free(rpd_theo);
        free(x_region);
        free(M_region);
    }

    // ===============================================================================
    // Free memory.
    // ===============================================================================
    for (unsigned int i = 0; i < rows_m; i++) {
        free(data_m[i]);
    }
    free(data_m);
    free(M_tot);
    free(x_tot);

    // ===============================================================================
    // Finsh statement.
    // ===============================================================================
    printf("Process finished.\n");
    return 0;
}