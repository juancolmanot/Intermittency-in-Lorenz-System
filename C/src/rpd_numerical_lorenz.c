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

int main(int argc, char *argv[]) {

    // ===============================================================================
    // Check for right passed arguments.
    // ===============================================================================
    if (argc != 4){
        printf("Wrong amount of arguments, you passed %d, but 4 are required\n.", argc);
        printf("Usage: ./run-script.sh script_to_run write_to_file.dat read_from_file.dat\n");
        exit(EXIT_FAILURE);
    }

    // ===============================================================================
    // Load arguments into variables.
    // ===============================================================================
    const char *write_filename = argv[2];
    const char *read_filename = argv[3];
    
    // ===============================================================================
    // Store data from read file.
    // ===============================================================================
    unsigned int rows = 0, cols = 0;
    double **data;
    read_data_file_unsigned_double(read_filename, &data, &rows, &cols);

    // ===============================================================================
    // File to write to
    // ===============================================================================
    FILE *f = open_file(write_filename);

    // ===============================================================================
    // Load reinjected points.
    // ===============================================================================
    double yf = 41.2861;
    double *xreinj1 = calloc(rows, sizeof(double));
    double *xreinj2 = calloc(rows, sizeof(double));
    double *xreinjabs = calloc(rows, sizeof(double));
    for (unsigned int i = 0; i < rows; i++) {
        xreinj1[i] = data[i][1];
        xreinj2[i] = data[i][1] - yf;
        xreinjabs[i] = fabs(yf - data[i][1]);
    }

    // ===============================================================================
    // Compute numerical RPD
    // ===============================================================================
    unsigned int nbins = 800;
    double bins1[nbins], rpd1[nbins], bins2[nbins], rpd2[nbins], binsabs[nbins], rpdabs[nbins];
    double xmin1, xmax1, dx1;
    double xmin2, xmax2, dx2;
    double xminabs, xmaxabs, dxabs;
    xmin1 = la_min_d(xreinj1, rows);
    xmax1 = la_max_d(xreinj1, rows);
    dx1 = (xmax1 - xmin1) / (double)(nbins - 1);
    xmin2 = la_min_d(xreinj2, rows);
    xmax2 = la_max_d(xreinj2, rows);
    dx2 = (xmax2 - xmin2) / (double)(nbins - 1);
    xminabs = la_min_d(xreinjabs, rows);
    xmaxabs = la_max_d(xreinjabs, rows);
    dxabs = (xmaxabs - xminabs) / (double)(nbins - 1);
    for (unsigned int i = 0; i < nbins; i++) {
        bins1[i] = xmin1 + dx1 * (double)i;
        bins2[i] = xmin2 + dx2 * (double)i;
        binsabs[i] = xminabs + dxabs * (double)i;
        rpd1[i] = 0;
        rpd2[i] = 0;
        rpdabs[i] = 0;
    }
    stats_histogram_double(rpd1, bins1, xreinj1, rows, nbins);
    stats_histogram_double(rpd2, bins2, xreinj2, rows, nbins);
    stats_histogram_double(rpdabs, binsabs, xreinjabs, rows, nbins);

    // Normalize RPDs.
    double b1 = 0.0, b2 = 0.0, babs = 0.0;
    b1 = normalize_histogram_double(bins1, rpd1, nbins, 1.0);
    b2 = normalize_histogram_double(bins2, rpd2, nbins, 1.0);
    babs = normalize_histogram_double(binsabs, rpdabs, nbins, 1.0);

    for (unsigned int i = 0; i < nbins; i++) {
        if (rpd1[i] < rows && rpd2[i] < rows && rpdabs[i] < rows) {
            fprintf(f, "%12.5E %12.5E %12.5E %12.5E %12.5E %12.5E\n",
                bins1[i], b1 * rpd1[i],
                bins2[i], b2 * rpd2[i],
                binsabs[i], babs * rpdabs[i]
            );
        }
    }

    // ===============================================================================
    // Free memory.
    // ===============================================================================
    for (unsigned int i = 0; i < rows; i++) {
        free(data[i]);
    }
    free(data);
    free(xreinj1);
    free(xreinj2);
    free(xreinjabs);
    fclose(f);

    // ===============================================================================
    // Finish statement.
    // ===============================================================================
    printf("Process finished. Results stored in %s\n", write_filename);
    return 0;
}