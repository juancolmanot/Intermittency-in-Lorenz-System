    #include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <omp.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "../include/lorenz.h"
#include "../../../common/include/file_handle.h"
#include "../../../common/include/ini.h"
#include "../../../common/include/parameters.h"
#include "../../../common/include/stats.h"
#include "../../../common/include/linear_algebra.h"
#include "../../../common/include/gsl_utilities.h"

int main(int argc, char *argv[]) {

    // ===============================================================================
    // Check for right passed arguments.
    // ===============================================================================
    if (argc != 4){
        perror("Wrong amount of arguments.");
        printf("Usage: ./run-script.sh script_to_run write_to_file.dat read_from_file.dat\n");
        exit(EXIT_FAILURE);
    }

    // ===============================================================================
    // Load arguments into variables.
    // ===============================================================================
    const char *write_filename = argv[2];
    const char *read_file = argv[3];
    // ===============================================================================
    // File to write to
    // ===============================================================================
    FILE *f = open_file(write_filename);

    // ===============================================================================
    // Store data from laminar lengths file.
    // ===============================================================================
    unsigned int laminar_rows = 0, laminar_cols = 0;
    double **laminar_data;
    read_data_file_unsigned_double(read_file, &laminar_data, &laminar_rows, &laminar_cols);
    double ylaminar[laminar_rows];
    for (unsigned int i = 0; i < laminar_rows; i++) {
        ylaminar[i] = laminar_data[i][0];
    }

    // ===============================================================================
    // Compute numerical pdll.
    // ===============================================================================
    double ylaminar_sorted[laminar_rows];
    for (unsigned int i = 0; i < laminar_rows; i++) {
        ylaminar_sorted[i] = ylaminar[i];
    }
    quicksort_double_unsigned(ylaminar_sorted, laminar_rows);
    double lmax, lmin;
    lmin = ylaminar_sorted[0];
    lmax = ylaminar_sorted[laminar_rows - 1];
    unsigned int nlam = (unsigned int)(lmax - lmin);
    double pdll[nlam];
    for (unsigned int i = 0; i < nlam - 1; i++) {
        for (unsigned int j = 0; j < laminar_rows; j++) {
            if (ylaminar_sorted[j] >= lmin + (double)i && ylaminar_sorted[j] < lmin + (double)(i + 1)) {
                pdll[i]++;
            }    
        }
    }

    // ===============================================================================
    // Write results to files.
    // ===============================================================================
    for (unsigned int i = 0; i < nlam - 1; i++){
        fprintf(f, "%5.5f %5.5f\n", lmin + (double)i, pdll[i]);
    }

    // ===============================================================================
    // Close file.
    // ===============================================================================    
    fclose(f);

    // ===============================================================================
    // Finish statement.
    // ===============================================================================
    printf("Process finished. Results stored in %s.\n", write_filename);
    return 0;
}
