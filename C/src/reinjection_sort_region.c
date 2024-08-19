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
    if (argc != 6){
        printf("Wrong amount of arguments. You passed %d, but 6 are required.\n", argc);
        printf("Usage: ./run-script.sh script_to_run write_to_file_base write_counts_filename.dat ");
        printf("reinjected_file.dat regions_file.dat\n");
        exit(EXIT_FAILURE);
    }

    // ===============================================================================
    // Load arguments into variables.
    // ===============================================================================
    const char *write_filename = argv[2];
    const char *write_counts_filename = argv[3];
    const char *reinjections_file = argv[4];
    const char *regions_file = argv[5];

    // ===============================================================================
    // Files to write to
    // ===============================================================================
    FILE *f = open_file(write_counts_filename);

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
    // Array to count reinjections by zone.
    // ===============================================================================
    unsigned int regions_count[rows2];
    for (unsigned int i = 0; i < rows2; i++) {
        regions_count[i] = 0;
    }

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
                regions_count[j]++;
                char region_filename[256];
                sprintf(region_filename, "%s_%d.dat", write_filename, j);
                FILE *region_file = fopen(region_filename, "a");
                if (region_file) {
                    fprintf(region_file, "%5.5f %5.5f\n", xn, xn1);
                    fclose(region_file);
                } else {
                    perror("Failed to open region file");
                }
                break;
            }
        }
    }

    // ===============================================================================
    // Write to file amount of reinjections per region.
    // ===============================================================================
    for (unsigned int i = 0; i < rows2; i++) {
        fprintf(f, "%5.5f %5.5f %d\n", regions_data[i][0], regions_data[i][1], regions_count[i]);
        printf("reinjection count region %d: %d\n", i, regions_count[i]);
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
    fclose(f);
    
    // ===============================================================================
    // Finish statement.
    // ===============================================================================
    printf("Process finished. Results stored in %s\n", write_counts_filename);
    return 0;
}