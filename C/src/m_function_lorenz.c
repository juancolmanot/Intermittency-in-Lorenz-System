#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include "../include/lorenz.h"
#include "../../../common/include/file_handle.h"
#include "../../../common/include/stats.h"
#include "../../../common/include/linear_algebra.h"

int main(int argc, char *argv[]) {

    // ===============================================================================
    // Check for right passed arguments.
    // ===============================================================================
    if (argc != 4){
        printf("Wrong amount of arguments, you passed %d, but 4 are required\n.", argc);
        printf("Usage: ./run.sh script_to_run write_to_file.dat read_from_file.dat\n");
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
    // Fixed point.
    // ===============================================================================
    double yf = 41.2861;

    // ===============================================================================
    // Load reinjected points.
    // ===============================================================================
    double *xreinj = calloc(rows, sizeof(double));
    double *xreinj_fp = calloc(rows, sizeof(double));
    for (unsigned int i = 0; i < rows; i++) {
        xreinj[i] = data[i][1];
        xreinj_fp[i] = data[i][1] -  yf;
    }

    // ===============================================================================
    // Sort reinjected points.
    // ===============================================================================
    bubble_sort_double_unsigned(xreinj, rows);
    bubble_sort_double_unsigned(xreinj_fp, rows);

    // ===============================================================================
    // Compute M function
    // ===============================================================================
    double Mi = 0.0, Mi_fp = 0.0;
    for (unsigned int i = 0; i < rows; i++) {
        Mi += xreinj[i];
        Mi_fp += xreinj_fp[i];
        fprintf(f, "%12.5E %12.5E %12.5E %12.5E %12.5E %12.5E\n",
            xreinj[i],
            Mi / (double)(i + 1),
            xreinj_fp[i],
            (Mi_fp / (double)(i + 1)),
            xreinj[i] - xreinj[0],
            (Mi / (double)(i + 1)) - (xreinj[0] / 1.0)
        );
    }

    // ===============================================================================
    // Free memory.
    // ===============================================================================
    for (unsigned int i = 0; i < rows; i++) {
        free(data[i]);
    }
    free(data);
    free(xreinj);
    free(xreinj_fp);
    fclose(f);

    // ===============================================================================
    // Finish statement.
    // ===============================================================================
    printf("Process finished. Results stored in %s\n", write_filename);
    return 0;
}