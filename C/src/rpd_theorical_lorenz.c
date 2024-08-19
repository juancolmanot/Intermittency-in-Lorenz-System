#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include "../include/lorenz.h"
#include "../../../common/include/file_handle.h"
#include "../../../common/include/parameters.h"
#include "../../../common/include/ini.h"
#include "../../../common/include/stats.h"
#include "../../../common/include/linear_algebra.h"
#include "../../../common/include/gsl_utilities.h"

int main(int argc, char *argv[]) {

    // ===============================================================================
    // Check for right passed arguments.
    // ===============================================================================
    if (argc != 7){
        printf("Wrong amount of arguments, you passed %d, but 7 are required\n.", argc);
        printf("Usage: ./run-script.sh script_to_run rpd_write_file");
        printf(" rpd_theoric_write_file reinject_data.dat");
        printf(" reinject_count_data.dat params_file.ini\n");
        exit(EXIT_FAILURE);
    }

    // ===============================================================================
    // Load arguments into variables.
    // ===============================================================================
    const char *rpd_write_file = argv[2];
    const char *rpd_theoric_write_file = argv[3];
    const char *reinject_data = argv[4];
    const char *reinject_count_data = argv[5];
    const char *params_file = argv[6];

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
    read_data_file_unsigned(reinject_data, &data, &rows, &cols);

    // ===============================================================================
    // Store data from reinjection count file.
    // ===============================================================================
    unsigned int rows_c = 0, cols_c = 0;
    long double **count_data;
    read_data_file_unsigned(reinject_count_data, &count_data, &rows_c, &cols_c);
    unsigned int count_tot = 0, count_zone = 0;
    count_zone = (unsigned int)count_data[params.current_zone][2];
    for (unsigned int i = 0; i < rows_c; i++) {
        count_tot += (unsigned int)count_data[i][2];
    }

    // ===============================================================================
    // Array to store numerical and theoric regions:
    //  - The structure of it is: A[2 * nregions][ni] (ni could be either the number of
    //  points for the M function or the number of bins for the RPD).
    // ===============================================================================
    long double **RPD, *b_norm, *alphas, *x_fp, **xdomain;
    RPD = (long double **)calloc((size_t)(3 * params.n_regions), sizeof(long double));
    xdomain = (long double **)calloc((size_t)(params.n_regions), sizeof(long double));
    for (int i = 0; i < params.n_regions; i++) {
        xdomain[i] = (long double *)calloc(2, sizeof(long double));
    }
    alphas = calloc((size_t)params.n_regions, sizeof(long double));
    b_norm = calloc((size_t)params.n_regions, sizeof(long double));
    x_fp = calloc((size_t)params.n_regions, sizeof(long double));

    // ===============================================================================
    // Load reinjected points.
    // ===============================================================================
    long double *xreinj = calloc(rows, sizeof(long double));
    for (unsigned int j = 0; j < rows; j++) {
        xreinj[j] = data[j][3];
    }

    // ===============================================================================
    // Loop through regions of M function.
    // ===============================================================================
    for (int i = 0; i < params.n_regions; i++) {

        // ===============================================================================
        // Sort reinjected points.
        // ===============================================================================
        quicksort_long_unsigned(xreinj, rows);
        
        // ===============================================================================
        // Compute M function and get slope
        // ===============================================================================
        long double xstart_fit, xend_fit;
        xstart_fit = (long double)params.xmins_fit[i];
        xend_fit = (long double)params.xmaxs_fit[i];
        xdomain[i][0] = (long double)params.xmins_domain[i];
        xdomain[i][1] = (long double)params.xmaxs_domain[i];
        unsigned int npoints_fit = 0, npoints_domain = 0;
        for (unsigned int j = 0; j < rows; j++) {
            if (xreinj[j] >= xdomain[i][0] && xreinj[j] < xdomain[i][1]) {
                npoints_fit++;
            }
            if (xreinj[j] >= xdomain[i][0] && xreinj[j] < xdomain[i][1]) {         
                npoints_domain++;
            }
        }
        long double Mx[npoints_fit], x_dom[npoints_domain];
        long double Mi = 0.0;
        long double *xfit, *Mfit;
        xfit = malloc(npoints_fit * sizeof(long double));
        Mfit = malloc(npoints_fit * sizeof(long double));
        long double mi = 0.0, bi = 0.0;
        unsigned int k = 0, l = 0;
        for (unsigned int j = 0; j < rows; j++) {
            Mi += xreinj[j];
            if (xreinj[j] >= xstart_fit && xreinj[j] < xend_fit) {
                xfit[k] = xreinj[j];
                Mx[k] = Mi / (long double)(j + 1);
                Mfit[k] = Mx[k];
                k++;
            }
            if (xreinj[j] >= xdomain[i][0] && xreinj[j] < xdomain[i][1]) {         
                x_dom[l] = xreinj[j];
                l++;
            }
        }


        linear_regression(xfit, Mfit, npoints_fit, &mi, &bi);
        
        // ===============================================================================
        // Compute numerical RPD.
        // ===============================================================================
        unsigned int nbins = params.region_bins[i];
        long double bins[nbins - 2], rpd[nbins - 2];
        long double xmin, xmax, dx;
        long double wi = (long double)npoints_domain / (long double)count_tot;
        xmin = la_min(x_dom, npoints_domain);
        xmax = la_max(x_dom, npoints_domain);
        dx = (xmax - xmin) / (long double)(nbins - 1);
        for (unsigned int j = 1; j < nbins - 1; j++) {
            bins[j - 1] = xmin + dx * (long double)j;
        }
        stats_histogram(rpd, bins, x_dom, npoints_domain, nbins - 1);
        
        // ===============================================================================
        // Normalize RPD.
        // ===============================================================================
        unsigned int mtecarlopoints = 10000;
        long double integral_rpd = montecarlo_integration_long(
            bins,
            rpd,
            nbins - 2,
            mtecarlopoints
        );
        long double const_norm = wi / integral_rpd;
        for (unsigned int j = 1; j < nbins - 3; j++) {
            rpd[j] = rpd[j] * const_norm;
        }

        // ===============================================================================
        // Compute theoretical RPD.
        // ===============================================================================
        long double alpha = (2 * mi - 1) / (1 - mi);
        long double xci = (long double)params.xc[i];
        long double rpd_theoric[nbins - 2];
        for (unsigned int j = 1; j < nbins - 3; j++) {
            rpd_theoric[j] = powl(fabsl(bins[j] - xci), alpha);
        }
        long double integral_rpd_theoric = montecarlo_integration_long(
            bins,
            rpd_theoric,
            nbins - 2,
            mtecarlopoints
        );
        long double b_rpd = wi / integral_rpd_theoric;
        alphas[i] = alpha;
        b_norm[i] = b_rpd;
        x_fp[i] = xci;

        // ===============================================================================
        // Allocate corresponding pointers for storing results.
        // ===============================================================================
        for (int j = i * 3; j < (i + 1) * 3; j++) {
            RPD[j] = (long double *)malloc((nbins - 2) * sizeof(long double));
        }
        
        // ===============================================================================
        // Write into given arrays
        // ===============================================================================
        for (unsigned int k = 0; k < nbins - 3; k++) {
            RPD[(size_t)(i * 3)][k] = bins[k];
            RPD[(size_t)(i * 3 + 1)][k] = rpd[k];
            RPD[(size_t)(i * 3 + 2)][k] = b_rpd * rpd_theoric[k];
        }
        // ===============================================================================
        // Free memory.
        // ===============================================================================
        free(xreinj);
        free(xfit);
        free(Mfit);
    }

    // ===============================================================================
    // Compute Numerical whole RPD.
    // ===============================================================================
    unsigned int nbins_tot = 0;
    for (int j = 0; j < params.n_regions; j++){
        nbins_tot += params.region_bins[j];
    }
    long double bins_tot[nbins_tot - 2], rpd_tot[nbins_tot - 2];
    long double xmin_tot, xmax_tot, dx_tot;
    xmin_tot = la_min(xreinj, rows);
    xmax_tot = la_max(xreinj, rows);
    dx_tot = (xmax_tot - xmin_tot) / (long double)(nbins_tot - 1);
    for (unsigned int j = 1; j < nbins_tot - 1; j++) {
        bins_tot[j - 1] = xmin_tot + dx_tot * (long double)j;
    }
    stats_histogram(rpd_tot, bins_tot, xreinj, rows, nbins_tot - 1);

    // ===============================================================================
    // Normalize whole RPD.
    // ===============================================================================
    unsigned int mtecarlopoints = 10000;
    long double integral_rpd_tot = montecarlo_integration_long(
        bins_tot,
        rpd_tot,
        nbins_tot - 2,
        mtecarlopoints
    );
    long double const_norm_tot = 1.0 / integral_rpd_tot;
    for (unsigned int j = 1; j < nbins_tot - 3; j++) {
        rpd_tot[j] = rpd_tot[j] * const_norm_tot;
    }
    // ===============================================================================
    // Write into given files
    // ===============================================================================
    FILE *rpd_file, *rpd_theoric_file;
    rpd_file = open_file(rpd_write_file);
    rpd_theoric_file = open_file(rpd_theoric_write_file);
    fprintf(rpd_theoric_file, "%12s %12s %12s %12s %12s %12s\n", 
        "zone",
        "b",
        "alpha",
        "xc",
        "x_start",
        "x_end"
    );
    unsigned int rpd_idx = 0;
    for (int i = 0; i < params.n_regions; i++) {
        for (unsigned int j = 1; j < params.region_bins[i] - 3; j++) {
            rpd_idx++;
            fprintf(rpd_file, "%5.10Lf %5.10Lf %5.10Lf %5.10Lf\n",
                RPD[i * 3][j],
                RPD[(i * 3) + 1][j],
                RPD[(i * 3) + 2][j],
                rpd_tot[rpd_idx]
            );
        }
        fprintf(rpd_theoric_file, "%12d %12.3LE %12.3LE %12.3LE %12.3LE %12.3LE\n",
            i,
            b_norm[i],
            alphas[i],
            x_fp[i],
            xdomain[i][0],
            xdomain[i][1]
        );
    }

    // ===============================================================================
    // Free memory.
    // ===============================================================================
    for (unsigned int i = 0; i < rows; i++) {
        free(data[i]);
    }
    for (int i = 0; i < (3 * params.n_regions); i++) {
        free(RPD[i]);
    }
    for (int i = 0; i < params.n_regions; i++){
        free(xdomain[i]);
    }
    free(b_norm);
    free(x_fp);
    free(alphas);
    free(data);
    free(RPD);
    free(xdomain);
    free(params.xmins_fit);
    free(params.xmaxs_fit);
    free(params.xmins_domain);
    free(params.xmaxs_domain);
    free(params.region_bins);
    free(params.xc);
    // ===============================================================================
    // Finish statement.
    // ===============================================================================
    printf("Process finished. Results stored in %s and %s.\n", rpd_write_file, rpd_theoric_write_file);
    return 0;
}