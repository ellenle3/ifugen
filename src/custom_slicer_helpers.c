#define _CRT_SECURE_NO_DEPRECATE
#include <stdio.h>
#include <stdlib.h>
#include "custom_slicer_helpers.h"

#define NUM_PARAMS_PER_SLICE 10

void LoadCustomParamsFromFile(double p_custom[], int file_num, char params_dir[], int max_elements) {
    if (file_num > 9999 || file_num < -999) {
        printf("Error: File number cannot exceed 4 digits\n");
        return;
    }
    char basename[35];
    snprintf(basename, sizeof(basename), "p_custom_%d.txt", file_num);

    char filename[547];  // 512 + 35
    snprintf(filename, sizeof(filename), "%s%s", params_dir, basename);

    FILE *file = fopen(filename, "r");
    if (!file) {
        printf("Error: Could not open file %s\n", filename);
        return;
    }

    int n_rows, n_cols;
    int surface_type;
    double dx, dy, gx_width, gx_depth, gy_width, gy_depth;

    // Read first line
    if (fscanf(file, "%d %d %d\n", &n_rows, &n_cols, &surface_type) != 3) {
        printf("Error: Invalid format on first line\n");
        fclose(file);
        return;
    }

    // Read second line
    if (fscanf(file, "%lf %lf %lf %lf %lf %lf\n", &dx, &dy, &gx_width, &gx_depth, &gy_width, &gy_depth) != 6) {
        printf("Error: Invalid format on second line\n");
        fclose(file);
        return;
    }

    int num_slices = n_rows * n_cols;
    int array_size = 9 + n_rows + NUM_PARAMS_PER_SLICE * num_slices;
    if (array_size > max_elements) {
        printf("Error: Number of entries exceeds maximum limit of %d\n", max_elements);
        fclose(file);
        return;
    }

    // Zero-initialize array
    for (int i = 0; i < array_size; i++) p_custom[i] = 0.0;

    // Store first 9 entries
    p_custom[0] = (double)n_rows;
    p_custom[1] = (double)n_cols;
    p_custom[2] = (double)surface_type;
    p_custom[3] = dx;
    p_custom[4] = dy;
    p_custom[5] = gx_width;
    p_custom[6] = gx_depth;
    p_custom[7] = gy_width;
    p_custom[8] = gy_depth;

    // Read third line: u values
    for (int i = 0; i < n_rows; i++) {
        double u_val;
        if (fscanf(file, "%lf", &u_val) == 1) {
            p_custom[9 + i] = u_val;
        } else {
            p_custom[9 + i] = 0.0; // missing values default to zero
        }
    }

    // Read slice parameters directly into p_custom
    double alpha, beta, gamma, cv, k, zp, syx, syz, sxy, sxz;
    int base_idx;
    for (int slice_idx = 0; slice_idx < num_slices; slice_idx++) {
        if (fscanf(file, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",
                &alpha, &beta, &gamma, &cv, &k, &zp, &syx, &syz, &sxy, &sxz) == 10) {
            base_idx = 9 + n_rows + NUM_PARAMS_PER_SLICE * slice_idx;
            p_custom[base_idx]     = alpha;
            p_custom[base_idx + 1] = beta;
            p_custom[base_idx + 2] = gamma;
            p_custom[base_idx + 3] = cv;
            p_custom[base_idx + 4] = k;
            p_custom[base_idx + 5] = zp;
            p_custom[base_idx + 6] = syx;
            p_custom[base_idx + 7] = syz;
            p_custom[base_idx + 8] = sxy;
            p_custom[base_idx + 9] = sxz;
        } else {
            break; // stop if fewer than 10 values found
        }
    }

    fclose(file);
}

SLICE_PARAMS GetSliceParamsCustom(int row_num, int col_num, double p_custom[]) {
    SLICE_PARAMS slice = {0};  // Initialize all fields to zero
    int n_rows = (int)p_custom[0];
    int n_cols = (int)p_custom[1];

    // Safety check
    if (row_num < 0 || row_num >= n_rows || col_num < 0 || col_num >= n_cols) {
        printf("Error: Requested slice (%d, %d) out of bounds (%d, %d)\n",
            row_num, col_num, n_rows, n_cols);
        return slice;  // Return zeroed slice on error
    }

    int u_start_idx = 9;
    int slice_params_start_idx = 9 + n_rows;

    // Compute the 1D index of the slice in the flattened array (row-major)
    int row_idx = row_num + col_num * n_rows;

    // Get the u value for this row
    slice.u = p_custom[u_start_idx + row_num];

    // Get the 10 slice parameters
    int base_idx = slice_params_start_idx + NUM_PARAMS_PER_SLICE * row_idx;
    slice.alpha = p_custom[base_idx];
    slice.beta  = p_custom[base_idx + 1];
    slice.gamma = p_custom[base_idx + 2];
    slice.cv    = p_custom[base_idx + 3];
    slice.k     = p_custom[base_idx + 4];
    slice.zp    = p_custom[base_idx + 5];
    slice.syx   = p_custom[base_idx + 6];
    slice.syz   = p_custom[base_idx + 7];
    slice.sxy   = p_custom[base_idx + 8];
    slice.sxz   = p_custom[base_idx + 9];

    return slice;
}
