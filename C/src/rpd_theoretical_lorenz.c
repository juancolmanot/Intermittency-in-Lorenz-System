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
    if (argc != 5){
        printf("Wrong amount of arguments, you passed %d, but 5 are required\n.", argc);
        printf("Usage: ./run-script.sh script_to_run rpd_write_file");
        printf(" reinject_data.dat params_file.ini\n");
        exit(EXIT_FAILURE);
    }

    // ===============================================================================
    // Load arguments into variables.
    // ===============================================================================
    const char *rpd_write_file = argv[2];
    const char *reinject_data = argv[3];
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
    read_data_file_unsigned(reinject_data, &data, &rows, &cols);

    // ===============================================================================
    // Array to store numerical and theoric regions:
    //  - The structure of it is: A[2 * nregions][ni] (ni could be either the number of
    //  points for the M function or the number of bins for the RPD).
    // ===============================================================================
    long double *alphas, *x_fp, **xdomain;
    alphas = calloc((size_t)params.n_regions, sizeof(long double));
    x_fp = calloc((size_t)params.n_regions, sizeof(long double));
    xdomain = (long double **)calloc((size_t)(params.n_regions), sizeof(long double));
    for (int i = 0; i < params.n_regions; i++) {
        xdomain[i] = (long double *)calloc(2, sizeof(long double));
    }

    // ===============================================================================
    // Load reinjected points.
    // ===============================================================================
    long double *xreinj = calloc(rows, sizeof(long double));
    for (unsigned int j = 0; j < rows; j++) {
        xreinj[j] = data[j][3];
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
        long double x_dom[npoints_domain];
        long double Mi = 0.0;
        long double *xfit, *Mfit;
        xfit = calloc(npoints_fit, sizeof(long double));
        Mfit = calloc(npoints_fit, sizeof(long double));
        long double mi = 0.0, bi = 0.0;
        unsigned int k = 0, l = 0;
        for (unsigned int j = 0; j < rows; j++) {
            Mi += xreinj[j];
            if (xreinj[j] >= xstart_fit && xreinj[j] < xend_fit) {
                xfit[k] = xreinj[j];
                Mfit[k] = Mi / (long double)(j + 1);
                k++;
            }
            if (xreinj[j] >= xdomain[i][0] && xreinj[j] < xdomain[i][1]) {         
                x_dom[l] = xreinj[j];
                l++;
            }
        }

        linear_regression(xfit, Mfit, npoints_fit, &mi, &bi);
        
        // ===============================================================================
        // Compute theoretical RPD.
        // ===============================================================================
        long double alpha = (2 * mi - 1) / (1 - mi);
        long double xci = 0.0;
        if (alpha < -1.0) {
            xci = xdomain[i][1];
        }
        else if (alpha >= -1.0) {
            xci = xdomain[i][0];
        }

        alphas[i] = alpha;
        x_fp[i] = xci;
        printf("%d %Lf %Lf\n", i, xci, alpha);
        
        // ===============================================================================
        // Free memory.
        // ===============================================================================
        free(xfit);
        free(Mfit);
    }
 
    // ===============================================================================
    // Compute Numerical RPD.
    // ===============================================================================
    unsigned int nbins_tot = 0;
    for (int j = 0; j < params.n_regions; j++){
        nbins_tot += params.region_bins[j];
    }
    long double *bins_tot, *rpd_tot;
    bins_tot = calloc(nbins_tot, sizeof(long double));
    rpd_tot = calloc(nbins_tot, sizeof(long double));
    long double xmin_tot, xmax_tot, dx_tot;
    xmin_tot = la_min(xreinj, rows);
    xmax_tot = la_max(xreinj, rows);
    dx_tot = (xmax_tot - xmin_tot) / (long double)(nbins_tot - 1);
    for (unsigned int j = 1; j < nbins_tot; j++) {
        bins_tot[j - 1] = xmin_tot + dx_tot * (long double)j;
    }
    stats_histogram(rpd_tot, bins_tot, xreinj, rows, nbins_tot);

    // ===============================================================================
    // Normalize Numerical RPD.
    // ===============================================================================
    unsigned int mtecarlopoints = 100000;
    long double integral_rpd_tot = montecarlo_integration_long(
        bins_tot,
        rpd_tot,
        nbins_tot,
        mtecarlopoints
    );
    long double const_norm_tot = 1.0 / integral_rpd_tot;
    for (unsigned int j = 0; j < nbins_tot; j++) {
        rpd_tot[j] *= const_norm_tot;
    }

    // ===============================================================================
    // Compute Theoretical RPD.
    // ===============================================================================
    unsigned int nbins = nbins_tot;
    long double *bins_theoric, *rpd_theoric;
    bins_theoric = calloc(nbins, sizeof(long double));
    rpd_theoric = calloc(nbins, sizeof(long double));
    long double xmin, xmax, dx;
    xmin = la_min(xreinj, rows);
    xmax = la_max(xreinj, rows);
    dx = (xmax - xmin) / (long double)(nbins - 1);
    for (unsigned int j = 1; j < nbins; j++) {
        bins_theoric[j - 1] = xmin + dx * (long double)j;
    }

    for (unsigned int i = 0; i < nbins; i++) {
        
        for (int j = 0; j < params.n_regions; j++){
            rpd_theoric[i] += powl(fabsl(bins_theoric[i] - x_fp[j]), alphas[j]);
        }
        
        //rpd_theoric[i] = powl(fabsl(bins_theoric[i] - x_fp[4]), alphas[4]);
    }

    // ===============================================================================
    // Normalize theoretical RPD.
    // ===============================================================================
    /*
    long double integral_rpd_theoric = montecarlo_integration_long(
        bins_theoric,
        rpd_theoric,
        nbins_tot,
        mtecarlopoints
    );
    long double const_norm_theoric = 1.0L / integral_rpd_theoric;
    for (unsigned int j = 0; j < nbins_tot; j++) {
        rpd_theoric[j] *= const_norm_theoric;
    }
    */

    // ===============================================================================
    // Write into given files
    // ===============================================================================
    FILE *rpd_file = open_file(rpd_write_file);
        
    for (unsigned int i = 0; i < nbins_tot; i++){
        fprintf(rpd_file, "%12.3LE %12.3LE %12.3LE %12.3LE\n",
            bins_tot[i],
            bins_theoric[i],
            rpd_tot[i],
            rpd_theoric[i]
        );
    }

    // ===============================================================================
    // Free memory.
    // ===============================================================================
    for (unsigned int i = 0; i < rows; i++) {
        free(data[i]);
    }
    for (unsigned int i = 0; i < params.n_regions; i++){
        free(xdomain[i]);
    }
    free(xdomain);
    free(bins_tot);
    free(bins_theoric);
    free(rpd_tot);
    free(rpd_theoric);
    free(xreinj);
    free(x_fp);
    free(alphas);
    free(data);
    free(params.xmins_fit);
    free(params.xmaxs_fit);
    free(params.xmins_domain);
    free(params.xmaxs_domain);
    free(params.region_bins);
    free(params.xc);

    fclose(rpd_file);

    // ===============================================================================
    // Finish statement.
    // ===============================================================================
    printf("Process finished. Results stored in %s.\n", rpd_write_file);
    return 0;
}