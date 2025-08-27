#define _CRT_SECURE_NO_DEPRECATE
#include <stdio.h>
#include <stdlib.h>
#include "surface_solns.h"

void LoadCustomParamsFromFile(double p_custom[], int file_num, char params_dir[], int max_elements) {

    if (file_num > 9999 || file_num < -999) {
        printf("Error: File number cannot exceed 4 digits\n");
        return;
    }

    // Get base name of file
    char basename[35];
    snprintf(basename, sizeof(basename), "custom_mirror_array_params_%d.txt", file_num);

    // Concatenate to the absolute path of the file
    char filename[547];  // 512 + 35
    snprintf(filename, sizeof(filename), "%s%s", params_dir, basename);

    FILE *file = fopen(filename, "r");
    if (!file) {
        printf("Error: Could not open file %s\n", filename);
        return;
    }

    int n_slices_per_col, n_cols, surface_type;
    double dx, dy, gx_width, gx_depth, gy_width, gy_depth;

    // Read first line (integers)
    if (fscanf(file, "%d %d %d\n", &n_slices_per_col, &n_cols, &surface_type) != 3) {
        printf("Error: Invalid format on first line\n");
        fclose(file);
        return;
    }

    // Read second line (doubles)
    if (fscanf(file, "%lf %lf %lf %lf %lf %lf\n", &dx, &dy, &gx_width, &gx_depth, &gy_width, &gy_depth) != 6) {
        printf("Error: Invalid format on second line\n");
        fclose(file);
        return;
    }

    // Calculate length of array based on the total number of slices
    int num_slices = n_slices_per_col * n_cols;
    int array_size = 9 + num_slices * 5; // First 9 entries + slice parameters

    if (array_size > max_elements) {
        printf("Error: Number of entries exceeds maximum limit of %d\n", max_elements);
        fclose(file);
        return;
    }

    // Set all values to zero initially
    for (int i = 0; i < array_size; i++) {
        p_custom[i] = 0;
    }

    // Store first 9 parameters
    p_custom[0] = (double) n_slices_per_col;
    p_custom[1] = (double) n_cols;
    p_custom[2] = (double) surface_type;
    p_custom[3] = dx;
    p_custom[4] = dy;
    p_custom[5] = gx_width;
    p_custom[6] = gx_depth;
    p_custom[7] = gy_width;
    p_custom[8] = gy_depth;

    // Read slice parameters
    double alpha, beta, gamma, cv, k;
    int base_idx;
    for (int i = 0; i < num_slices; i++) {
        if (fscanf(file, "%lf %lf %lf %lf %lf\n", &alpha, &beta, &gamma, &cv, &k) == 5) {
            base_idx = 9 + 5 * i;
            p_custom[base_idx] = alpha;
            p_custom[base_idx + 1] = beta;
            p_custom[base_idx + 2] = gamma;
            p_custom[base_idx + 3] = cv;
            p_custom[base_idx + 4] = k;
        } else {
            break; // Stop reading if fewer than expected values are found
        }
    }

    fclose(file);
    // Make sure to free the array when you're done!
}