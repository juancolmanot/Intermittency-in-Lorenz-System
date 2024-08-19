#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include "../include/lorenz.h"
#include "../../../common/include/file_handle.h"
#include "../../../common/include/stats.h"
#include "../../../common/include/linear_algebra.h"
#include "../../../common/include/calculus.h"

int main(int argc, char *argv[]) {

    // ===============================================================================
    // Check for right passed arguments.
    // ===============================================================================
    if (argc != 4){
        printf("Wrong amount of arguments, you passed %d, but 5 are required\n.", argc);
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
    double *xreinj = calloc(rows, sizeof(double));
    for (unsigned int i = 0; i < rows; i++) {
        xreinj[i] = data[i][1];
    }

    // ===============================================================================
    // Sort reinjected points.
    // ===============================================================================
    quicksort_double_unsigned(xreinj, rows);

    // ===============================================================================
    // Compute M function.
    // ===============================================================================
    double Mx[rows], Mxi = 0.0;
    for (unsigned int i = 0; i < rows; i++) {
        Mxi += xreinj[i];
        Mx[i] = Mxi / (double)(i + 1);
    }

    // ===============================================================================
    // Compute local slope.
    // ===============================================================================
    unsigned int sample_size = 100;
    double *slopes = local_slope(Mx, xreinj, rows, sample_size);
    for (unsigned int i = 1; i < rows - 1; i++) {
        fprintf(f, "%5.5f %5.5f %5.5f\n", xreinj[i], Mx[i], slopes[i]);
    }


    // ===============================================================================
    // Free memory.
    // ===============================================================================
    for (unsigned int i = 0; i < rows; i++) {
        free(data[i]);
    }
    free(data);
    free(xreinj);
    free(slopes);
    fclose(f);

    // ===============================================================================
    // Finish statement.
    // ===============================================================================
    printf("Process finished. Results stored in %s\n", write_filename);
    return 0;
}