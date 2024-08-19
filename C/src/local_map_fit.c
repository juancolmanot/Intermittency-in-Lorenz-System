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
    if (argc != 5){
        perror("Wrong amount of arguments.");
        printf("Usage: ./run-script.sh script_to_run write_to_file_map.dat");
        printf(" write_to_file_coeffs.dat read_from_file.dat\n");
        exit(EXIT_FAILURE);
    }

    // ===============================================================================
    // Load arguments into variables.
    // ===============================================================================
    const char *write_filename_map = argv[2];
    const char *write_filename_coeffs = argv[3];
    const char *read_filename = argv[4];

    // ===============================================================================
    // Store data from read file.
    // ===============================================================================
    unsigned int rows = 0, cols = 0;
    double **data;
    read_data_file_unsigned_double(read_filename, &data, &rows, &cols);

    // ===============================================================================
    // File to write to
    // ===============================================================================
    FILE *f1, *f2;
    f1 = open_file(write_filename_map);
    f2 = open_file(write_filename_coeffs);

    // ===============================================================================
    // Fit map
    // ===============================================================================
    double xfit[rows], yfit[rows], ci[3], xfit_sort[rows];
    for (unsigned int i = 0; i < rows; i++) {
        xfit[i] = data[i][2];
        yfit[i] = data[i][3];
        xfit_sort[i] = data[i][2];
    }

    gsl_regression_quadratic(xfit, yfit, rows, ci);
    printf("%5.5f %5.5f %5.5f\n", ci[0], ci[1], ci[2]);
    // ===============================================================================
    // Write results to file.
    // ===============================================================================
    quicksort_double_unsigned(xfit_sort, rows);
    for (unsigned int i = 0; i < rows; i++) {
        fprintf(f1, "%5.5f %5.5f %5.5f %5.5f\n",
            xfit[i],
            yfit[i],
            xfit_sort[i],
            ci[2] * xfit_sort[i] * xfit_sort[i] + ci[1] * xfit_sort[i] + ci[0]
        );
    }

    fprintf(f2, "%5.5f %5.5f %5.5f\n", ci[2], ci[1], ci[0]);

    // ===============================================================================
    // Free memory.
    // ===============================================================================
    for (unsigned int i = 0; i < rows; i++) {
        free(data[i]);
    }
    free(data);

    // ===============================================================================
    // Close file.
    // ===============================================================================    
    fclose(f1);
    fclose(f2);

    // ===============================================================================
    // Finish statement.
    // ===============================================================================
    printf("Process finished. Results stored in %s\n & %s\n", write_filename_map, write_filename_coeffs);
    return 0;
}   
