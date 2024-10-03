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
        printf("Usage: ./run-script.sh script_to_run m_rpd_write_file");
        printf(" reinjection_data.dat params_file.ini\n");
        exit(EXIT_FAILURE);
    }

    // ===============================================================================
    // Load arguments into variables.
    // ===============================================================================
    const char *m_rpd_write_file = argv[2];
    const char *reinjection_data = argv[3];
    const char *params_file = argv[4];

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
    char file_m_all[256], file_rpd_all[256];
    sprintf(file_m_all, "%s_m_all.dat", m_rpd_write_file);
    sprintf(file_rpd_all, "%s_rpd_all.dat", m_rpd_write_file);
    FILE *fm_all = fopen(file_m_all, "w"), *frpd_all = fopen(file_rpd_all, "w");
    fclose(fm_all);
    fclose(frpd_all);

    for (int i = 0; i < params.n_regions; i++){
        
        // ===============================================================================
        // Get region's array of data.
        // ===============================================================================
        long double *x_region = calloc(100, sizeof(long double));
        unsigned int n_region = 0;
        long double xmin = params.xmins_fit[i], xmax = params.xmaxs_fit[i];
        get_region_array(xreinj, xmin, xmax, rows, &x_region, &n_region);

        // ===============================================================================
        // Get region M values, from M_tot array.
        // ===============================================================================
        unsigned int idxmin = get_value_loc(xreinj, x_region[0], rows);
        unsigned int idxmax = get_value_loc(xreinj, x_region[n_region - 1], rows);

        // Load portion of M(x) inside this M_region:
        long double *M_region = calloc(n_region, sizeof(long double));
        for (unsigned int j = idxmin; j < idxmax + 1; j++) {
            M_region[j - idxmin] = M_tot[j];
        }
        
        // ===============================================================================
        // Compute M(x) slope and theoric M function.
        // ===============================================================================
        long double m, b;
        linear_regression(x_region, M_region, n_region, &m, &b);

        // Now use the domain values for the theoretical function
        long double *x_region_d = calloc(100, sizeof(long double));
        unsigned int n_region_d = 0;
        long double xmin_d = params.xmins_domain[i], xmax_d = params.xmaxs_domain[i];
        get_region_array(xreinj, xmin_d, xmax_d, rows, &x_region_d, &n_region_d);

        long double *M_theoretical = calloc(n_region_d, sizeof(long double));
        long double *M_numerical = calloc(n_region_d, sizeof(long double));
        unsigned int idxmin_d = get_value_loc(xreinj, x_region_d[0], rows);
        unsigned int idxmax_d = get_value_loc(xreinj, x_region_d[n_region_d - 1], rows);
        for (unsigned int j = 0; j < n_region_d; j++) {
            M_theoretical[j] = b + m * x_region_d[j];
        }
        for (unsigned int j = idxmin_d; j < idxmax_d + 1; j++) {
            M_numerical[j - idxmin_d] = M_tot[j];
        }

        // Compute alpha.
        long double alpha = 0.0;
        alpha = (2 * m - 1) / (1 - m);
        /*
        if (m > 0.51 && m < 1.0) {
            alpha = (-2 * m - 1) / (1 + m);
        }
        else if (m < 0.51 && m > 0.0) {
            alpha = (2 * m - 1) / (1 - m);
        }
        else if (fabsl(m) > 1) {
            alpha = (2 * m - 1) / (1 - m);
        }
        */

        // ===============================================================================
        // Compute region theorical RPD(x).
        // ===============================================================================
        unsigned int nbins_i = params.region_bins[i];
        long double *bins_i = calloc(nbins_i, sizeof(long double));
        long double *rpd_i = calloc(nbins_i, sizeof(long double));
        long double xci = params.xc[i];
        // Normalize rpd.
        long double int_rpd = 0.0;
        long double lbound = x_region_d[0], ubound = x_region_d[n_region_d - 1];
        int_rpd = rpd_theoretical_integral(alpha, xci, lbound, ubound);
        int_rpd = fabsl(int_rpd);
        long double bnorm = (long double)(n_region_d)/ ((long double)(rows) * int_rpd);
        // Compute rpd theoretical
        for (unsigned int j = 0; j < nbins_i; j++) {
            bins_i[j] = x_region_d[0] + (x_region_d[n_region_d - 1] - x_region_d[0]) * (long double)j / (long double)(nbins_i - 1);
            rpd_i[j] = bnorm * powl(fabsl(bins_i[j] - xci), alpha);
        }

        // Compute region RPD(x) numerical;
        long double *rpd_region = calloc(nbins_i, sizeof(long double));
        long double wi = (long double)n_region_d / (long double)rows;
        compute_rpd_numerical(x_region_d, n_region_d, bins_i, rpd_region, nbins_i, wi);
        
        printf("m: %Lf, alpha: %Lf\n", m, alpha);
        printf("%% points in region: %Lf\n", (long double)n_region_d / (long double)rows);
        printf("===================================================================================\n");

        // ===============================================================================
        // Write results to files.
        // ===============================================================================
        char file_m[256], file_rpd[256];
        sprintf(file_m, "%s_m_%d.dat", m_rpd_write_file, i);
        sprintf(file_rpd, "%s_rpd_%d.dat", m_rpd_write_file, i);
        FILE *fm = fopen(file_m, "w");
        FILE *frpd = fopen(file_rpd, "w");
        // Open general files to append.
        FILE *fm_all = fopen(file_m_all, "a");
        FILE *frpd_all = fopen(file_rpd_all, "a");
        for (unsigned int j = 0; j < idxmax_d - idxmin_d; j++) {
            fprintf(fm, "%d %12.5LE %12.5LE %12.5LE\n",
                j,
                x_region_d[j],
                M_numerical[j],
                M_theoretical[j]
            );
            fprintf(fm_all, "%12.5LE %12.5LE %12.5LE\n",
                x_region_d[j],
                M_numerical[j],
                M_theoretical[j]
            );
        }
        for (unsigned int j = 1; j < nbins_i - 1; j++) {
            fprintf(frpd, "%d %12.5LE %12.5LE %12.5LE\n",
                j,
                bins_i[j],
                rpd_region[j],
                rpd_i[j]
            );
            fprintf(frpd_all, "%12.5LE %12.5LE %12.5LE\n",
                bins_i[j],
                rpd_region[j],
                rpd_i[j]
            );
        }
        printf("Results stored in %s & %s.\n", file_m, file_rpd);
        fclose(fm);
        fclose(frpd);
        fclose(fm_all);
        fclose(frpd_all);
        free(x_region);
        free(x_region_d);
        free(M_region);
        free(rpd_region);
        free(M_theoretical);
        free(M_numerical);
        free(rpd_i);
        free(bins_i);
    }

    // ===============================================================================
    // Finsh statement.
    // ===============================================================================
    printf("Process finished.\n");
    return 0;
}