#include <stdio.h>
#include <stdlib.h>

void LoadCustomParamsFromFile(double *custom_slice_params, int file_num, char params_dir[], int max_elements, int max_path_length) {

    if (file_num > 9999 || file_num < -999) {
        printf("Error: File number cannot exceed 4 digits\n");
        return;
    }

    // Get base name of file
    int max_basename_length = 35; // 31 characters for file name + 4 for number
    char basename[max_basename_length];
    snprintf(basename, sizeof(basename), "custom_mirror_array_params_%d.txt", file_num);

    // Concatenate to the absolute path of the file
    int max_filename_length = max_path_length + max_basename_length;
    char filename[max_filename_length];
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
        custom_slice_params[i] = 0;
    }

    // Store first 9 parameters
    custom_slice_params[0] = (double) n_slices_per_col;
    custom_slice_params[1] = (double) n_cols;
    custom_slice_params[2] = (double) surface_type;
    custom_slice_params[3] = dx;
    custom_slice_params[4] = dy;
    custom_slice_params[5] = gx_width;
    custom_slice_params[6] = gx_depth;
    custom_slice_params[7] = gy_width;
    custom_slice_params[8] = gy_depth;

    // Read slice parameters
    double alpha, beta, gamma, cv, k;
    int base_idx;
    for (int i = 0; i < num_slices; i++) {
        if (fscanf(file, "%lf %lf %lf %lf %lf\n", &alpha, &beta, &gamma, &cv, &k) == 5) {
            base_idx = 9 + 5 * i;
            custom_slice_params[base_idx] = alpha;
            custom_slice_params[base_idx + 1] = beta;
            custom_slice_params[base_idx + 2] = gamma;
            custom_slice_params[base_idx + 3] = cv;
            custom_slice_params[base_idx + 4] = k;
        } else {
            break; // Stop reading if fewer than expected values are found
        }
    }

    fclose(file);
    // Make sure to free the array when you're done!
}

void GetSliceParamsCustom(double* alpha, double* beta, double* gamma, double* cv, double* k,
int slice_num, int col_num, double custom_slice_params[]) {

    int n_slices_per_col = (int) custom_slice_params[0];
    int start_idx = 9 + 5 * (col_num * n_slices_per_col + slice_num);

    *alpha = custom_slice_params[start_idx];
    *beta = custom_slice_params[start_idx + 1];
    *gamma = custom_slice_params[start_idx + 2];
    *cv = custom_slice_params[start_idx + 3];
    *k = custom_slice_params[start_idx + 4];
}