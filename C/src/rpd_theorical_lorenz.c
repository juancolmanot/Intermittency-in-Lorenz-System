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
    long double *RPD, *b_norm, *alphas, *x_fp, **xdomain;
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
    // Count total bins.
    // ===============================================================================
    unsigned int nbins_tot = 0;
    for (int j = 0; j < params.n_regions; j++){
        nbins_tot += params.region_bins[j];
    }
    RPD = calloc(nbins_tot, sizeof(long double));

    // ===============================================================================
    // Compute Numerical whole RPD.
    // ===============================================================================
    long double bins_tot[nbins_tot - 2], rpd_tot[nbins_tot - 2];
    bins_tot[0] = rpd_tot[0] = 0.0;
    long double xmin_tot, xmax_tot, dx_tot;
    xmin_tot = la_min(xreinj, rows);
    xmax_tot = la_max(xreinj, rows);
    dx_tot = (xmax_tot - xmin_tot) / (long double)(nbins_tot - 1);
    for (unsigned int j = 1; j < nbins_tot - 1; j++) {
        bins_tot[j - 1] = xmin_tot + dx_tot * (long double)j;
        rpd_tot[j - 1] = 0.0;
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

    for (unsigned int j = 1; j < nbins_tot - 1; j++) {
        rpd_tot[j] = rpd_tot[j] * const_norm_tot;
    }

    // ===============================================================================
    // Sort reinjected points.
    // ===============================================================================
    quicksort_long_unsigned(xreinj, rows);

    // ===============================================================================
    // Loop through regions of M function.
    // ===============================================================================
    for (int i = 0; i < params.n_regions; i++) {
        
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
            if (xreinj[j] >= xstart_fit && xreinj[j] < xend_fit) {
                npoints_fit++;
            }
            if (xreinj[j] >= xdomain[i][0] && xreinj[j] < xdomain[i][1]) {
                npoints_domain++;
            }
        }
        long double Mx[npoints_fit], xdom[npoints_domain];
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
                xdom[l] = xreinj[j];
                l++;
            }
        }

        linear_regression(xfit, Mfit, npoints_fit, &mi, &bi);

        // ===============================================================================
        // Compute theoretical RPD.
        // ===============================================================================
        unsigned int nbins = params.region_bins[i];
        long double bins[nbins], rpd_theoric[nbins];
        long double xmin, xmax, dx;
        long double wi = (long double)npoints_domain / (long double)count_tot;
        xmin = la_min(xdom, npoints_domain);
        xmax = la_max(xdom, npoints_domain);
        dx = (xmax - xmin) / (long double)(nbins - 1);
        long double alpha = (2 * mi - 1) / (1 - mi);
        long double xci = params.xc[i];

        for (unsigned int j = 0; j < nbins; j++) {
            bins[j] = 0.0;
            bins[j] = xmin + dx * (long double)j;
            rpd_theoric[j] = powl(fabsl(bins[j] - xci), alpha);
        }
        long double integral_rpd_theoric = 0;
        integral_rpd_theoric = montecarlo_integration_long(
            bins,
            rpd_theoric,
            nbins - 2,
            mtecarlopoints
        );
        long double b_rpd = wi / integral_rpd_theoric;
        alphas[i] = alpha;
        b_norm[i] = b_rpd;
        x_fp[i] = xci;
        printf("%d %d %Lf\n", npoints_domain, count_tot, wi);
        printf("%d %Lf %Lf %Lf %Lf\n", i, mi, alpha, b_rpd, integral_rpd_theoric);

        // ===============================================================================
        // Write into given arrays
        // ===============================================================================
        for (unsigned int k = 0; k < nbins_tot; k++) {
            if (bins_tot[k] > xmin && bins_tot[k] < xmax){
                RPD[k] = b_rpd * powl(fabsl(bins_tot[k] - xci), alpha);
            }
        }
        // ===============================================================================
        // Free memory.
        // ===============================================================================
        free(xfit);
        free(Mfit);
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
    for (int i = 0; i < params.n_regions; i++) {
        fprintf(rpd_theoric_file, "%12d %12.3LE %12.3LE %12.3LE %12.3LE %12.3LE\n",
            i,
            b_norm[i],
            alphas[i],
            x_fp[i],
            xdomain[i][0],
            xdomain[i][1]
        );
    }

    for (unsigned int i = 1; i < nbins_tot; i++) {
        if ((bins_tot[i] == 0.0 && rpd_tot[i] == 0.0) || (bins_tot[i] == 0.0 && RPD[i] == 0.0)) {
            printf("zeroes\n");
        }
        else {
            fprintf(rpd_file, "%12.5LE %12.5LE %12.5LE\n",
                bins_tot[i],
                rpd_tot[i],
                RPD[i]
            );
        }
    }


    // ===============================================================================
    // Free memory.
    // ===============================================================================
    for (unsigned int i = 0; i < rows; i++) {
        free(data[i]);
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
    free(xreinj);
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