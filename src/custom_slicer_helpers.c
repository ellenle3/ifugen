#include <stdio.h>
#include <stdlib.h>

double* load_slice_params_file(int file_num, int *array_size) {

    if (file_num > 9999 || file_num < -999) {
        printf("Error: File number cannot exceed 4 digits\n");
        return NULL;
    }
    char filename[36]; // 11 characters for file name + 4 for number
    snprintf(filename, sizeof(filename), "custom_mirror_array_params_%d.txt", file_num);
    
    FILE *file = fopen(filename, "r");
    if (!file) {
        printf("Error: Could not open file %s\n", filename);
        return NULL;
    }

    int n_slices_per_col, n_cols, surface_type;
    double dx, dy, gx_width, gx_depth, gy_width, gy_depth;

    // Read first line (integers)
    if (fscanf(file, "%d %d %d", &n_slices_per_col, &n_cols, &surface_type) != 3) {
        printf("Error: Invalid format on first line\n");
        fclose(file);
        return NULL;
    }

    // Read second line (doubles)
    if (fscanf(file, "%lf %lf %lf %lf %lf %lf", &dx, &dy, &gx_width, &gx_depth, &gy_width, &gy_depth) != 6) {
        printf("Error: Invalid format on second line\n");
        fclose(file);
        return NULL;
    }

    // Calculate length of array based on the total number of slices
    int num_slices = n_slices_per_col * n_cols;
    *array_size = 9 + num_slices * 5; // First 9 entries + slice parameters

    // Allocate memory for parameter array
    double *slice_param_arr = (double *)malloc(*array_size * sizeof(double));
    if (!slice_param_arr) {
        printf("Error: Memory allocation failed\n");
        fclose(file);
        return NULL;
    }

    // Store first 9 parameters
    slice_param_arr[0] = (double) n_slices_per_col;
    slice_param_arr[1] = (double) n_cols;
    slice_param_arr[2] = (double) surface_type;
    slice_param_arr[3] = dx;
    slice_param_arr[4] = dy;
    slice_param_arr[5] = gx_width;
    slice_param_arr[6] = gx_depth;
    slice_param_arr[7] = gy_width;
    slice_param_arr[8] = gy_depth;

    // Read slice parameters

    for (int i = 0; i < *array_size; i++) {
        slice_param_arr[i] = 0.0;
    }

    int i = 0;
    while (i < num_slices * 5) {
        double alpha, beta, gamma, c, k;
        if (fscanf(file, "%lf %lf %lf %lf %lf", &alpha, &beta, &gamma, &c, &k) == 5) {
            slice_param_arr[9 + i] = alpha;
            slice_param_arr[10 + i] = beta;
            slice_param_arr[11 + i] = gamma;
            slice_param_arr[12 + i] = c;
            slice_param_arr[13 + i] = k;
            i += 5;
        } else {
            break; // Stop reading if fewer than expected values are found
        }
    }

    fclose(file);
    return slice_param_arr;
}