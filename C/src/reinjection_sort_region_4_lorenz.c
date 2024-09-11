#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include "../include/lorenz.h"
#include "../../../common/include/file_handle.h"
#include "../../../common/include/parameters.h"
#include "../../../common/include/ini.h"

int main(int argc, char *argv[]) {

    /*
    ==================================================================================
      This script sorts the points of the reinjection map in region 4 of the domain
     by how they are located in relationship with two lines, put there to isolate the
     effect of each structure in the RPD.
    ==================================================================================
    */

    // ===============================================================================
    // Check for right passed arguments.
    // ===============================================================================
    if (argc != 7){
        printf("Wrong amount of arguments. You passed %d, but 6 are required.\n", argc);
        printf("Usage: ./run-script.sh script_to_run write_to_file_1 write_to_file_2");
        printf("write_to_file_3 functions_file.dat reinjections_file.dat\n");
        exit(EXIT_FAILURE);
    }

    // ===============================================================================
    // Load arguments into variables.
    // ===============================================================================
    const char *write_filename_1 = argv[2];
    const char *write_filename_2 = argv[3];
    const char *write_filename_3 = argv[4];
    const char *functions_file = argv[5];
    const char *reinjections_file = argv[6];

    // ===============================================================================
    // Store data from read files.
    // ===============================================================================
    unsigned int rows1 = 0, cols1 = 0;
    double **reinjections_data;
    read_data_file_unsigned_double(reinjections_file, &reinjections_data, &rows1, &cols1);

    unsigned int rows2 = 0, cols2 = 0;
    double **functions_data;
    read_data_file_unsigned_double(functions_file, &functions_data, &rows2, &cols2);

    // ===============================================================================
    // Load lines parameters.
    // ===============================================================================
    double a1, a2, b1, b2;
    a1 = functions_data[0][0];
    b1 = functions_data[0][1];
    a2 = functions_data[1][0];
    b2 = functions_data[1][1];

    // ===============================================================================
    // Sort reinjected points by region.
    // ===============================================================================
    for (unsigned int i = 0; i < rows1; i++) {
        double xn = reinjections_data[i][0];
        double xn1 = reinjections_data[i][3];
        double z = reinjections_data[i][4];
        double fx1, fx2;
        fx1 = a1 * xn + b1;
        fx2 = a2 * xn + b2;
        
        if (xn1 > fx1) {
            FILE *region_file_1 = fopen(write_filename_1, "a");
            if (region_file_1) {
                fprintf(region_file_1, "%12.5E %12.5E %12.5E\n", xn, xn1, z);
                fclose(region_file_1);
            } else {
                perror("Failed to open region file");
            }
        }
        else if (xn1 < fx1 && xn1 > fx2) {
            FILE *region_file_2 = fopen(write_filename_2, "a");
            if (region_file_2) {
                fprintf(region_file_2, "%12.5E %12.5E %12.5E\n", xn, xn1, z);
                fclose(region_file_2);
            } else {
                perror("Failed to open region file");
            }
        }
        else if (xn1 < fx2) {
            FILE *region_file_3 = fopen(write_filename_3, "a");
            if (region_file_3) {
                fprintf(region_file_3, "%12.5E %12.5E %12.5E\n", xn, xn1, z);
                fclose(region_file_3);
            } else {
                perror("Failed to open region file");
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
        free(functions_data[i]);
    }
    free(reinjections_data);
    free(functions_data);
    
    // ===============================================================================
    // Finish statement.
    // ===============================================================================
    printf("Process finished. Results stored in %s, %s and %s.\n",
        write_filename_1,
        write_filename_2,
        write_filename_3
    );
    return 0;
}