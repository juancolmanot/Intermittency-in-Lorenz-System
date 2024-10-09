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
    Parameters9 params;
    load_parameters_from_file(params_file, &params, handler9);

    // ===============================================================================
    // Store data from reinjection points file.
    // ===============================================================================
    unsigned int rows = 0, cols = 0;
    long double **data;
    read_data_file_unsigned(reinjection_data, &data, &rows, &cols);

    // ===============================================================================
    // Load reinjection points.
    // ===============================================================================
    long double *xprev = calloc(rows, sizeof(long double));
    long double *xreinj = calloc(rows, sizeof(long double));
    long double *xreinj_sort = calloc(rows, sizeof(long double));
    for (unsigned int i = 0; i < rows; i++) {
        xprev[i] = data[i][0];
        xreinj[i] = data[i][3];
        xreinj_sort[i] = data[i][3];
    }

    // Sort x array.
    quicksort_long_unsigned(xreinj_sort, rows);

    // ===============================================================================
    // Compute numerical M(x).
    // ===============================================================================
    long double *M_tot = calloc(rows, sizeof(long double));
    compute_m_function(xreinj_sort, M_tot, rows);

    // ===============================================================================
    // Compute numerical RPD(x).
    // ===============================================================================
    unsigned int nbins_total = 0;
    for (int i = 0; i < params.n_regions; i++) {
        nbins_total += (unsigned int)params.bins[i];
    }
    long double *rpd_total = calloc(nbins_total, sizeof(long double));
    long double *bins_total = calloc(nbins_total, sizeof(long double));
    compute_rpd_numerical(xreinj_sort, rows, bins_total, rpd_total, nbins_total, 1.0);

    
    // ===============================================================================
    // Open files to complete M(x) and RPD(x) and erase their content.
    // ===============================================================================
    char file_m_all[256], file_rpd_all[256];
    sprintf(file_m_all, "%s_prev_region_m_all.dat", m_rpd_write_file);
    sprintf(file_rpd_all, "%sprev_region_rpd_all.dat", m_rpd_write_file);
    FILE *fm_all = fopen(file_m_all, "w"), *frpd_all = fopen(file_rpd_all, "w");
    for (unsigned int i = 0; i < rows; i++) {
        fprintf(fm_all, "%12.5LE %12.5LE\n", xreinj_sort[i], M_tot[i]);
    }
    for (unsigned int i = 0; i < nbins_total; i++) {
        fprintf(fm_all, "%12.5LE %12.5LE\n", bins_total[i], rpd_total[i]);
    }
    fclose(fm_all);
    fclose(frpd_all);

    for (int i = 0; i < params.n_regions; i++){
        
        // ===============================================================================
        // Get region's array of data for previous to reinjected points.
        // ===============================================================================
        long double *x_prev_reg = calloc(100, sizeof(long double));
        unsigned int n_region = 0;
        long double xmin = params.xmins_regions[i], xmax = params.xmaxs_regions[i];
        get_region_array_scrambled(xprev, xmin, xmax, rows, &x_prev_reg, &n_region);
        
        // ===============================================================================
        // Compute numerical M function associated with given region.
        // ===============================================================================
        unsigned int *indeces = calloc(100, sizeof(unsigned int));
        unsigned int nidx = 0;
        
        get_region_indeces(xprev, xmin, xmax, rows, &indeces, &nidx);
        
        // Load portion of M(x) inside this M_region:
        long double *x_region = calloc(n_region, sizeof(long double));
        long double *M_region = calloc(n_region, sizeof(long double));
        long double Mi = 0.0;
        
        for (unsigned int j = 0; j < nidx; j++) {
            x_region[j] = xreinj[indeces[j]];
        }
                
        quicksort_long_unsigned(x_region, n_region);
        for (unsigned int j = 0; j < nidx; j++) {
            Mi += x_region[j];
            M_region[j] = Mi / (long double)(j + 1);
        }
        
        
        // ===============================================================================
        // Compute region numerical RPD(x).
        // ===============================================================================
        unsigned int nbins_i = (unsigned int) ((double)nidx / 100.0) ;
        long double *bins_i = calloc(nbins_i, sizeof(long double));
        long double *rpd_region = calloc(nbins_i, sizeof(long double));
        long double wi = (long double)n_region / (long double)rows;
        compute_rpd_numerical(x_region, n_region, bins_i, rpd_region, nbins_i, wi);
        
        printf("%% points in region: %5.2Lf%% - ", (long double)n_region * 100.0 / (long double)rows);
        printf("Number of points in region: %d\n", nidx);
        printf("===================================================================================\n");

        // ===============================================================================
        // Write results to files.
        // ===============================================================================
        char file_m[256], file_rpd[256];
        sprintf(file_m, "%s_m_%d_prev_region.dat", m_rpd_write_file, i);
        sprintf(file_rpd, "%s_rpd_%d_prev_region.dat", m_rpd_write_file, i);
        FILE *fm = fopen(file_m, "w");
        FILE *frpd = fopen(file_rpd, "w");
        for (unsigned int j = 1; j < n_region; j++) {
            fprintf(fm, "%12.7LE %12.7LE\n",
                x_region[j],
                M_region[j]
            );
        }
        for (unsigned int j = 1; j < nbins_i - 1; j++) {
            fprintf(frpd, "%12.7LE %12.7LE\n",
                bins_i[j],
                rpd_region[j]
            );
        }
        printf("Results stored in %s & %s.\n", file_m, file_rpd);
        printf("===================================================================================\n");
        fclose(fm);
        fclose(frpd);
        free(bins_i);
        free(rpd_region);
        free(indeces);
        free(x_region);
        free(M_region);
        free(x_prev_reg);
    }

    // ===============================================================================
    // Free memory.
    // ===============================================================================
    for (unsigned int i = 0; i < rows; i++) {
        free(data[i]);
    }
    free(data);
    free(xprev);
    free(xreinj);
    free(xreinj_sort);
    free(M_tot);
    free(rpd_total);
    free(bins_total);

    // ===============================================================================
    // Finsh statement.
    // ===============================================================================
    printf("Process finished.\n");
    return 0;
}