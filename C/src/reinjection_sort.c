#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include "../include/lorenz.h"
#include "../../../common/include/file_handle.h"
#include "../../../common/include/parameters.h"
#include "../../../common/include/ini.h"

int main(int argc, char *argv[]) {

    // ===============================================================================
    // Check for right passed arguments.
    // ===============================================================================
    if (argc != 5){
        printf("Wrong amount of arguments. You passed %d, but 5 are required.\n", argc);
        printf("Usage: ./run-script.sh script_to_run write_to_file.dat ");
        printf("reinjected_file.dat regions_file.dat\n");
        exit(EXIT_FAILURE);
    }

    // ===============================================================================
    // Load arguments into variables.
    // ===============================================================================
    const char *write_filename = argv[2];
    const char *reinjections_file = argv[3];
    const char *regions_file = argv[4];

    // ===============================================================================
    // Store data from read files.
    // ===============================================================================
    unsigned int rows1 = 0, cols1 = 0;
    double **reinjections_data;
    read_data_file_unsigned_double(reinjections_file, &reinjections_data, &rows1, &cols1);

    unsigned int rows2 = 0, cols2 = 0;
    double **regions_data;
    read_data_file_unsigned_double(regions_file, &regions_data, &rows2, &cols2);

    // ===============================================================================
    // Sort reinjected points by region.
    // ===============================================================================
    for (unsigned int i = 0; i < rows1; i++) {
        double xn = reinjections_data[i][0];
        double xn1 = reinjections_data[i][1];
        
        for (unsigned int j = 0; j < rows2; j++) {
            double xmin = regions_data[j][0];
            double xmax = regions_data[j][1];

            if (xn > xmin && xn < xmax) {
                char region_filename[256];
                sprintf(region_filename, "%s_%d.dat", write_filename, j + 1);
                FILE *region_file = fopen(region_filename, "a");
                if (region_file) {
                    fprintf(region_file, "%12.5E %12.5E\n", xn, xn1);
                    fclose(region_file);
                } else {
                    perror("Failed to open region file");
                }
                break;
            }
        }
    }

    // ===============================================================================
    // Freeing allocated memory.
    // ===============================================================================
    for (unsigned int i = 0; i < rows1; i++) {
        free(reinjections_data[i]);
    }
    for (unsigned int i = 0; i < rows2; i++) {
        free(regions_data[i]);
    }
    free(reinjections_data);
    free(regions_data);
    
    // ===============================================================================
    // Finish statement.
    // ===============================================================================
    printf("Process finished. Results stored in %s\n", write_filename);
    return 0;
}