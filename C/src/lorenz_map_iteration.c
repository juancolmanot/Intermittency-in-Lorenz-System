#include <stdlib.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../../../common/include/file_handle.h"

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
    // Load reinjected points.
    // ===============================================================================
    size_t n_map = 4;
    for (unsigned int i = 0; i < rows - n_map; i++) {
        if (data[i][0] != 0.0 && data[i + n_map][0]) {
            fprintf(f, "%12.5E %12.5E %12.5E\n",
                (double)i,
                data[i][0],
                data[i + n_map][0]
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
    fclose(f);

    // ===============================================================================
    // Finish statement.
    // ===============================================================================
    printf("Process finished. Results stored in %s\n", write_filename);
    return 0;
}