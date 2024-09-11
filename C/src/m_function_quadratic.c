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
    if (argc != 6){
        printf("Wrong amount of arguments, you passed %d, but 6 are required\n.", argc);
        printf("Usage: ./run-script.sh script_to_run m_write_file");
        printf(" m_theoric_write_file reinject_data.dat");
        printf(" params_file.ini\n");
        exit(EXIT_FAILURE);
    }

    // ===============================================================================
    // Load arguments into variables.
    // ===============================================================================
    const char *m_write_file = argv[2];
    const char *m_theoric_write_file = argv[3];
    const char *reinject_data = argv[4];
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
    read_data_file_unsigned(reinject_data, &data, &rows, &cols);

    // ===============================================================================
    // Array to store numerical and theoric regions:
    //  - The structure of it is: A[2 * nregions][ni] (ni could be either the number of
    //  points for the M function or the number of bins for the RPD).
    // ===============================================================================
    long double **M, *bm, *am, **xdomain;
    M = (long double **)calloc((size_t)(3 * params.n_regions), sizeof(long double));
    bm = (long double*)calloc((size_t)(params.n_regions), sizeof(long double));
    am = (long double*)calloc((size_t)(params.n_regions), sizeof(long double));
    xdomain = (long double **)calloc((size_t)(params.n_regions), sizeof(long double));
    for (int i = 0; i < params.n_regions; i++) {
        xdomain[i] = (long double *)calloc(2, sizeof(long double));
    }

    // ===============================================================================
    // Array to store counts of points per region and bins per región.
    // ===============================================================================
    unsigned int counts_region[params.n_regions][2];

    // ===============================================================================
    // Loop through regions of M function.
    // ===============================================================================
    for (int i = 0; i < params.n_regions; i++) {
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
        
        long double Mx[npoints_fit], x_dom[npoints_domain], M_dom[npoints_domain];
        long double Mi = 0.0;
        counts_region[i][0] = npoints_domain;
        counts_region[i][1] = npoints_fit;
        long double *xfit, *Mfit;
        xfit = calloc(npoints_fit, sizeof(long double));
        Mfit = calloc(npoints_fit,   sizeof(long double));
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
                M_dom[l] = Mi / (long double)(j + 1);
                l++;
            }
        }
        
        linear_regression(xfit, Mfit, npoints_fit, &mi, &bi);
        
        bm[i] = bi;
        am[i] = mi;

        // ===============================================================================
        // Allocate corresponding pointers for storing results.
        // ===============================================================================
        for (int j = i * 3; j < (i + 1) * 3; j++) {
            M[j] = (long double *)malloc(npoints_domain * sizeof(long double));
        }

        // ===============================================================================
        // Write into given arrays
        // ===============================================================================
        for (unsigned int k = 0; k < npoints_fit; k++) {
            M[(size_t)(i * 3)][k] = xfit[k];
            M[(size_t)(i * 3 + 1)][k] = Mfit[k];
            M[(size_t)(i * 3 + 2)][k] = bi + mi * xfit[k];
        }
        // ===============================================================================
        // Free memory.
        // ===============================================================================
        free(xreinj);
        free(xfit);
        free(Mfit);
    }
    // ===============================================================================
    // Write into given files
    // ===============================================================================
    FILE *m_file, *m_theoric_file;
    m_file = open_file(m_write_file);
    m_theoric_file = open_file(m_theoric_write_file);
    fprintf(m_theoric_file, "%12s %12s %12s %12s %12s\n", 
        "zone",
        "b",
        "m",
        "x_start",
        "x_end"
    );
    for (int i = 0; i < params.n_regions; i++) {
        for (unsigned int j = 0; j < counts_region[i][1]; j++) {
            fprintf(m_file, "%5.10Lf %5.10Lf %5.10Lf\n",
                M[i * 3][j],
                M[(i * 3) + 1][j],
                M[(i * 3) + 2][j]
            );  
        }
        fprintf(m_theoric_file, "%12d %12.3LE %12.3LE %12.3LE %12.3LE\n",
            i,
            bm[i],
            am[i],
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
        free(M[i]);
    }
    for (int i = 0; i < params.n_regions; i++){
        free(xdomain[i]);
    }
    free(bm);
    free(am);
    free(data);
    free(M);
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
    printf("Process finished. Results stored in %s and %s.\n", m_write_file, m_theoric_write_file);
    return 0;
}