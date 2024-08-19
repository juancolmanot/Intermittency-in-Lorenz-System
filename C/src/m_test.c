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
        printf("Usage: ./run-script.sh script_to_run m_write_base rpd_write_base");
        printf(" reinject_data.dat reinject_count_data.dat params_file.ini\n");
        exit(EXIT_FAILURE);
    }

    // ===============================================================================
    // Load arguments into variables.
    // ===============================================================================
    const char *m_write_base = argv[2];
    const char *rpd_write_base = argv[3];
    const char *reinject_data = argv[4];
    const char *reinject_count_data = argv[5];
    const char *params_file = argv[6];

    // ===============================================================================
    // Parameters.
    // ===============================================================================
    Parameters7 params;
    load_parameters_from_file(params_file, &params, handler7);

    // ===============================================================================
    // Store data from reinjection points file.
    // ===============================================================================
    unsigned int rows = 0, cols = 0;
    double **data;
    read_data_file_unsigned_double(reinject_data, &data, &rows, &cols);

    // ===============================================================================
    // Store data from reinjection count file.
    // ===============================================================================
    unsigned int rows_c = 0, cols_c = 0;
    double **count_data;
    read_data_file_unsigned_double(reinject_count_data, &count_data, &rows_c, &cols_c);
    unsigned int count_tot = 0, count_zone = 0;
    for (unsigned int i = 0; i < rows_c; i++) {
        count_tot += (unsigned int)count_data[i][2];
    }
    count_zone = (unsigned int)count_data[params.current_zone][params.current_zone];

    // ===============================================================================
    // Loop through regions of M function.
    // ===============================================================================
    for (int i = 0; i < params.n_regions; i++) {
        // ===============================================================================
        // Load reinjected points.
        // ===============================================================================
        double *xreinj = calloc(rows, sizeof(double));
        for (unsigned int j = 0; j < rows; j++) {
            xreinj[j] = data[j][1];
        }
        // ===============================================================================
        // Sort reinjected points.
        // ===============================================================================
        quicksort_double_unsigned(xreinj, rows);
        // ===============================================================================
        // Compute M function
        // ===============================================================================
        double Mx[rows];
        double Mi = 0.0;
        for (unsigned int j = 0; j < rows; j++) {
            Mi += xreinj[j];
            Mx[j] = Mi / (double)(j + 1);
        }
        // ===============================================================================
        // Get M slope
        // ===============================================================================
        double xstart, xend;
        xstart = params.xmins[i];
        xend = params.xmaxs[i];
        printf("%f %f\n", xstart, xend);
        unsigned int npoints = 0;
        for (unsigned int j = 0; j < rows; j++) {
            if (xreinj[j] >= xstart && xreinj[j] < xend) {
                npoints++;
            }
        }

        FILE *test_m = fopen("../datafiles/m_test_fit.dat", "w");
        double mi, bi;
        double xfit[npoints], Mfit[npoints];
        unsigned int k = 0;
        for (unsigned int j = 0; j < rows; j++) {
            if (xreinj[j] >= xstart && xreinj[j] < xend) {
                xfit[k] = xreinj[j];
                Mfit[k] = Mx[j];
                k++;
            }
        }
        gsl_linear_regression(xfit, Mfit, npoints, &mi, &bi);
        printf("mi %f\n", mi);s
        for (unsigned int j = 0; j < npoints; j++) {
            fprintf(test_m, "%f %f %f\n", xfit[j], Mfit[j], bi + mi * xfit[j]);
        }
        fclose(test_m);
        // ===============================================================================
        // Compute numerical RPD.
        // ===============================================================================
        /*unsigned int nbins = 100;
        double bins[nbins], rpd[nbins];
        double xmin, xmax, dx;
        double wi = (double)count_zone / (double)count_tot;
        xmin = la_min_d(xreinj, rows);
        xmax = la_max_d(xreinj, rows);
        dx = (xmax - xmin) / (double)(nbins - 1);
        for (unsigned int j = 0; j < nbins; j++) {
            bins[j] = xmin + dx * (double)j;
        }
        stats_histogram_double(rpd, bins, xreinj, rows, nbins);
        // ===============================================================================
        // Normalize RPD.
        // ===============================================================================
        unsigned int mtecarlopoints = 1000000;
        double integral_rpd = montecarlo_integration(bins, rpd, nbins, mtecarlopoints);
        double const_norm = wi / integral_rpd;
        for (unsigned int j = 0; j < nbins; j++) {
            rpd[j] = rpd[j] * const_norm;
        }
        // ===============================================================================
        // Compute theoretical RPD.
        // ===============================================================================
        double alpha = (2 * mi - 1) / (1 - mi);
        double xci = params.xc[i];
        double rpd_theoric[nbins];
        printf("%f %f %f\n", xci, mi, alpha);
        for (unsigned int j = 0; j < nbins; j++) {
            rpd_theoric[j] = pow(fabs(xci - bins[j]), alpha);
            //printf("%f\n", pow(xci - bins[j], alpha));
        }
        double integral_rpd_theoric = montecarlo_integration(bins, rpd_theoric, nbins, mtecarlopoints);
        double b_rpd = wi / integral_rpd_theoric;
        printf("%f  %f\n", integral_rpd, integral_rpd_theoric);
        // ===============================================================================
        // Write into given files
        // ===============================================================================
        char m_filename[256], rpd_filename[256];
        sprintf(m_filename, "%s_%d.dat", m_write_base, i);
        sprintf(rpd_filename, "%s_%d.dat", rpd_write_base, i);
        FILE *m_file, *rpd_file;
        m_file = open_file(m_filename);
        rpd_file = open_file(rpd_filename);
        if (m_file && rpd_file) {
            for (unsigned int j = 0; j < rows; j++) {
                fprintf(m_file, "%5.5f %5.5f %5.5f\n", xreinj[j], Mx[j], bi + mi * xreinj[j]);
            }

            for (unsigned int j = 0; j < nbins; j++) {
                fprintf(rpd_file, "%5.5f %5.5f %5.5f\n", bins[j], rpd[j], b_rpd * rpd_theoric[j]);
            }
        }
        else {
            perror("Failed to open files for writing m and/or rpd.");
        }
        // ===============================================================================
        // Free memory.
        // ===============================================================================
        free(xreinj);
        fclose(m_file);
        fclose(rpd_file);*/
    }
    // ===============================================================================
    // Free memory.
    // ===============================================================================
    for (unsigned int i = 0; i < rows; i++) {
        free(data[i]);
    }
    free(data);
    free(params.xmins);
    free(params.xmaxs);
    free(params.xc);

    // ===============================================================================
    // Finish statement.
    // ===============================================================================
    /*for (int i = 0; i < params.n_regions; i++) {
        char m_filename[256], rpd_filename[256];
        sprintf(m_filename, "%s_%d.dat", m_write_base, i);
        sprintf(rpd_filename, "%s_%d.dat", rpd_write_base, i);
        printf("Process finished. Results stored in %s and %s.\n", m_filename, rpd_filename);
    }*/
    return 0;
}