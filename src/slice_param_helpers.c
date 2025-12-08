#define _USE_MATH_DEFINES
#define _CRT_SECURE_NO_DEPRECATE
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "slice_param_helpers.h"

// assumes that p_custom[] is malloced for MAX_ELEMENTS
void LoadCustomParamsFromFile(double p_custom[], int file_num, char params_dir[]) {
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

    int n_rows, n_cols, is_linear;
    int surface_type;
    double dx, dy, gx_width, gx_depth, gy_width, gy_depth, f;

    // Read first line
    if (fscanf(file, "%d %d %d %d\n", &n_rows, &n_cols, &surface_type, &is_linear) != 4) {
        printf("Error: Invalid format on first line\n");
        fclose(file);
        return;
    }

    // Read second line
    if (fscanf(file, "%lf %lf %lf %lf %lf %lf %lf\n", &dx, &dy, &gx_width, &gx_depth, &gy_width, &gy_depth, &f) != 7) {
        printf("Error: Invalid format on second line\n");
        fclose(file);
        return;
    }

    int num_slices = n_rows * n_cols;
    int array_size = NUM_BASE_PARAMS + n_rows + NUM_PARAMS_PER_SLICE * num_slices;
    if (array_size > MAX_ELEMENTS) {
        printf("Error: Number of entries exceeds maximum limit of %d\n", MAX_ELEMENTS);
        fclose(file);
        return;
    }

    // Zero-initialize array
    for (int i = 0; i < array_size; i++) p_custom[i] = 0.0;

    if (f <= 0) {
        printf("Error: f must be positive. Setting to 100.\n");
        f = 100.0;
    }

    // Store first 10 entries
    p_custom[0] = (double)n_rows;
    p_custom[1] = (double)n_cols;
    p_custom[2] = 1;
    p_custom[3] = (double)surface_type;
    p_custom[4] = dx;
    p_custom[5] = dy;
    p_custom[6] = gx_width;
    p_custom[7] = gx_depth;
    p_custom[8] = gy_width;
    p_custom[9] = gy_depth;
    //number of params above must match NUM_BASE_PARAMS

    // Read third line: u values
    for (int i = 0; i < n_rows; i++) {
        double u_val;
        if (fscanf(file, "%lf", &u_val) == 1) {
            p_custom[NUM_BASE_PARAMS + i] = u_val;
        } else {
            p_custom[NUM_BASE_PARAMS + i] = 0.0; // missing values default to zero
        }
    }

    // Read slice parameters directly into p_custom
    double alpha, beta, gamma, cv, k, zp, syx, syz, sxy, sxz, theta, szx, szy;
    int base_idx;
    for (int slice_idx = 0; slice_idx < num_slices; slice_idx++) {
        if (fscanf(file, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",
                &alpha, &beta, &gamma, &cv, &k, &zp, &syx, &syz, &sxy, &sxz, &theta, &szx, &szy) == NUM_PARAMS_PER_SLICE) {
            base_idx = NUM_BASE_PARAMS + n_rows + NUM_PARAMS_PER_SLICE * slice_idx;

            if (is_linear) {
                // If params are defined in linear space, convert to angles first
                if (cv != 0) f = 1 / (2 * cv);
                Conic2DOffAxisAngle(&alpha, &beta, cv, k, beta, alpha);
                alpha *= 180.0 / M_PI;
                beta  *= 180.0 / M_PI;
                gamma = atan( gamma / f ) * 180.0 / M_PI; 
            }

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
            p_custom[base_idx + 10] = theta;
            p_custom[base_idx + 11] = szx;
            p_custom[base_idx + 12] = szy;
        } else {
            break; // stop if fewer than 10 values found
        }
    }

    fclose(file);
}

SLICE_PARAMS GetSliceParams(int slice_num, int col_num, double p_custom[]) {
    SLICE_PARAMS slice = {0};  // Initialize all fields to zero
    int n_rows = (int)p_custom[0];
    int n_cols = (int)p_custom[1];
    int n_each = (int)p_custom[2];

    // Safety check
    if (slice_num < 0 || slice_num >= n_rows * n_each || col_num < 0 || col_num >= n_cols) {
        printf("Requested slice (%d, %d) out of bounds (%d, %d)\n",
            slice_num, col_num, n_rows * n_each, n_cols);
        return slice;  // Return zeroed slice on error
    }

    int u_start_idx = NUM_BASE_PARAMS;
    int slice_params_start_idx = NUM_BASE_PARAMS + n_rows;

    // Compute the 1D index of the slice in the flattened array. slice_num is the
    // local index within the column while slice_idx is the global index.
    int slice_idx = slice_num + col_num * n_rows * n_each;
    int row_num = floor((double)slice_num / n_each);

    // Get the u value for this row
    slice.u = p_custom[u_start_idx + row_num];

    // Get the slice parameters
    int base_idx = slice_params_start_idx + NUM_PARAMS_PER_SLICE * slice_idx;
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
    slice.theta  = p_custom[base_idx + 10];
    slice.szx    = p_custom[base_idx + 11];
    slice.szy    = p_custom[base_idx + 12];

    return slice;
}

double GetUForRow(int row_num, double p_custom[]) {
    return p_custom[NUM_BASE_PARAMS + row_num];
}

// Check whether all rows have the same u value for triangle generation
int IsAllRowsAligned(double p_custom[]) {
    int n_rows = (int)p_custom[0];
    double first_u = p_custom[NUM_BASE_PARAMS];
    for (int i = 1; i < n_rows; i++) {
        if (p_custom[NUM_BASE_PARAMS + i] != first_u) {
            return 0; // Found a row with different u value
        }
    }
    return 1; // All rows have the same u value
}

GRID_PARAMS_BASIC MakeBasicParamsFromCustom(double p_custom[]) {
    GRID_PARAMS_BASIC p_basic;
    p_basic.n_rows = (int)p_custom[0];
    p_basic.n_cols = (int)p_custom[1];
    p_basic.n_each = (int)p_custom[2];
    p_basic.surface_type = (int)p_custom[3];
    p_basic.dx = p_custom[4];
    p_basic.dy = p_custom[5];
    p_basic.gx_width = p_custom[6];
    p_basic.gx_depth = p_custom[7];
    p_basic.gy_width = p_custom[8];
    p_basic.gy_depth = p_custom[9];
    return p_basic;
}

// Validate image slicer parameters, modifying illegal parameters as needed.
int ValidateBasicParams(GRID_PARAMS_BASIC* p) {
    // Keep track of whether we had to change any parameters
    int is_valid = 1;

    if (!(p->surface_type==0 || p->surface_type==1)) { p->surface_type=0; is_valid = 0; }
    if (p->n_cols < 1) { p->n_cols = 1; is_valid = 0; }
    if (p->n_rows < 1) { p->n_rows = 1; is_valid = 0; }
    if (p->n_each < 1) { p->n_each = 1; is_valid = 0; }

    // No need to check angles because we will convert them to be between -180
    // and 180 degrees...

    if (p->dx <= 0) { p->dx = 1; is_valid = 0; }
    if (p->dy <= 0) { p->dy = 1; is_valid = 0; }
    if (p->gx_width < 0) { p->gx_width = 0; is_valid = 0; }
    if (p->gy_width < 0) { p->gy_width = 0; is_valid = 0; }

    // Gap depths can be whatever
    // There are also no limitations on cv and k
    return is_valid;
}

int ValidateSlicerParamsAngular(IMAGE_SLICER_PARAMS_ANGULAR* p) {

    double* p_custom = (double*)malloc(MAX_ELEMENTS * sizeof(double));
    MakeSliceParamsArrayAngular(p_custom, *p);
    GRID_PARAMS_BASIC p_basic = MakeBasicParamsFromCustom(p_custom);
    int is_valid = ValidateBasicParams(&p_basic);

    p->n_cols = p_basic.n_cols;
    p->n_rows = p_basic.n_rows;
    p->n_each = p_basic.n_each;
    p->surface_type = p_basic.surface_type;
    p->dx = p_basic.dx;
    p->dy = p_basic.dy;
    p->gx_width = p_basic.gx_width;
    p->gx_depth = p_basic.gx_depth;
    p->gy_width = p_basic.gy_width;
    p->gy_depth = p_basic.gy_depth;
    
    if (!(p->angle_mode==0 || p->angle_mode==1 || p->angle_mode==2 || p->angle_mode==3)){
        p->angle_mode = 0; is_valid = 0;
        }

    free(p_custom);

    return is_valid;
}

int ValidateSlicerParamsLinear(IMAGE_SLICER_PARAMS_LINEAR* p) {

    double* p_custom = (double*)malloc(MAX_ELEMENTS * sizeof(double));
    MakeSliceParamsArrayLinear(p_custom, *p);
    GRID_PARAMS_BASIC p_basic = MakeBasicParamsFromCustom(p_custom);
    int is_valid = ValidateBasicParams(&p_basic);

    p->n_cols = p_basic.n_cols;
    p->n_rows = p_basic.n_rows;
    p->n_each = p_basic.n_each;
    p->surface_type = p_basic.surface_type;
    p->dx = p_basic.dx;
    p->dy = p_basic.dy;
    p->gx_width = p_basic.gx_width;
    p->gx_depth = p_basic.gx_depth;
    p->gy_width = p_basic.gy_width;
    p->gy_depth = p_basic.gy_depth;
    
    if (!(p->angle_mode==0 || p->angle_mode==1 || p->angle_mode==2 || p->angle_mode==3)){
        p->angle_mode = 0; is_valid = 0;
        }
    if (p->f <= 0) {
        p->f = 100.0; is_valid = 0;
    }

    free(p_custom);

    return is_valid;
}

// Checks whether every member in each of the IMAGER_SLICER_PARAMS structs are
// equivalent
int IsParametersEqualAngular(IMAGE_SLICER_PARAMS_ANGULAR p1, IMAGE_SLICER_PARAMS_ANGULAR p2) {
   if (
        p1.surface_type == p2.surface_type &&
        p1.n_each == p2.n_each &&
        p1.n_rows == p2.n_rows &&
        p1.n_cols == p2.n_cols &&
        p1.angle_mode == p2.angle_mode &&

        p1.dalpha == p2.dalpha &&
        p1.dbeta == p2.dbeta &&
        p1.dgamma == p2.dgamma &&
        p1.gamma_offset == p2.gamma_offset &&

        p1.azps == p2.azps &&
        p1.dsyx == p2.dsyx &&
        p1.dsyz == p2.dsyz &&
        p1.dsxy == p2.dsxy &&
        p1.dsxz == p2.dsxz &&
        p1.du == p2.du &&

        p1.alpha_cen == p2.alpha_cen &&
        p1.beta_cen == p2.beta_cen &&
        p1.gamma_cen == p2.gamma_cen &&
        p1.syx_cen == p2.syx_cen &&
        p1.syz_cen == p2.syz_cen &&
        p1.sxy_cen == p2.sxy_cen &&
        p1.sxz_cen == p2.sxz_cen &&
        p1.u_cen == p2.u_cen &&

        p1.dx == p2.dx &&
        p1.dy == p2.dy &&
        p1.cv == p2.cv &&
        p1.k == p2.k &&

        p1.gx_width == p2.gx_width &&
        p1.gx_depth == p2.gx_depth &&
        p1.gy_width == p2.gy_width &&
        p1.gy_depth == p2.gy_depth
   ) return 1;
   return 0;
}

int IsParametersEqualLinear(IMAGE_SLICER_PARAMS_LINEAR p1, IMAGE_SLICER_PARAMS_LINEAR p2) {
   if (
        p1.surface_type == p2.surface_type &&
        p1.n_each == p2.n_each &&
        p1.n_rows == p2.n_rows &&
        p1.n_cols == p2.n_cols &&
        p1.angle_mode == p2.angle_mode &&

        p1.dy0 == p2.dy0 &&
        p1.dx0 == p2.dx0 &&
        p1.dd == p2.dd &&
        p1.d_offset == p2.d_offset &&

        p1.azps == p2.azps &&
        p1.dsyx == p2.dsyx &&
        p1.dsyz == p2.dsyz &&
        p1.dsxy == p2.dsxy &&
        p1.dsxz == p2.dsxz &&
        p1.du == p2.du &&

        p1.y0_cen == p2.y0_cen &&
        p1.x0_cen == p2.x0_cen &&
        p1.d_cen == p2.d_cen &&
        p1.syx_cen == p2.syx_cen &&
        p1.syz_cen == p2.syz_cen &&
        p1.sxy_cen == p2.sxy_cen &&
        p1.sxz_cen == p2.sxz_cen &&
        p1.u_cen == p2.u_cen &&

        p1.dx == p2.dx &&
        p1.dy == p2.dy &&
        p1.cv == p2.cv &&
        p1.k == p2.k &&
        p1.f == p2.f &&

        p1.gx_width == p2.gx_width &&
        p1.gx_depth == p2.gx_depth &&
        p1.gy_width == p2.gy_width &&
        p1.gy_depth == p2.gy_depth
   ) return 1;
   return 0;
}

// Calculates u for the given row index. This is computed in a separate function
// as it is needed to calculate the column index of a slice.
void MakeSliceParamsArrayAngular(double p_custom[], IMAGE_SLICER_PARAMS_ANGULAR p) {
    int n_each = p.n_each;
    int n_rows = p.n_rows;
    int n_cols = p.n_cols;
    int num_slices = n_rows * n_cols;
    int array_size = NUM_BASE_PARAMS + n_rows + NUM_PARAMS_PER_SLICE * num_slices;
    if (array_size > MAX_ELEMENTS) {
        printf("Error: Number of entries exceeds maximum limit of %d\n", MAX_ELEMENTS);
        return;
    }

    p_custom[0] = (double)n_rows;
    p_custom[1] = (double)n_cols;
    p_custom[2] = (double)n_each;
    p_custom[3] = (double)p.surface_type;
    p_custom[4] = p.dx;
    p_custom[5] = p.dy;
    p_custom[6] = p.gx_width;
    p_custom[7] = p.gx_depth;
    p_custom[8] = p.gy_width;
    p_custom[9] = p.gy_depth;

    int slice_num, slice_idx, base_idx;

    // Calculate parameters for each slice
    for (int col = 0; col < n_cols; col++) {

        for (int row = 0; row < n_rows; row++) {

            for (int subidx = 0; subidx < n_each; subidx++) {
                slice_num = subidx + row * n_each;
                slice_idx = slice_num + col * n_rows * n_each;
                SLICE_PARAMS pslice = GetSliceParamsAngular(slice_num, col, p);
                if (subidx == 0 && col == 0) {
                    // Store u values only once per row
                    p_custom[NUM_BASE_PARAMS + row] = pslice.u;
                }

                base_idx = NUM_BASE_PARAMS + n_rows + NUM_PARAMS_PER_SLICE * slice_idx;
                p_custom[base_idx]     = pslice.alpha;
                p_custom[base_idx + 1] = pslice.beta;
                p_custom[base_idx + 2] = pslice.gamma;
                p_custom[base_idx + 3] = pslice.cv;
                p_custom[base_idx + 4] = pslice.k;
                p_custom[base_idx + 5] = pslice.zp;
                p_custom[base_idx + 6] = pslice.syx;
                p_custom[base_idx + 7] = pslice.syz;
                p_custom[base_idx + 8] = pslice.sxy;
                p_custom[base_idx + 9] = pslice.sxz;

                // no rotations about z-axis, this parameters should never be accessed
                p_custom[base_idx + 10] = 0.0;
                p_custom[base_idx + 11] = 0.0;
                p_custom[base_idx + 12] = 0.0;
            }
        }
    }
}

void MakeSliceParamsArrayLinear(double p_custom[], IMAGE_SLICER_PARAMS_LINEAR p) {
    int n_each = p.n_each;
    int n_rows = p.n_rows;
    int n_cols = p.n_cols;
    int num_slices = n_rows * n_cols;
    int array_size = NUM_BASE_PARAMS + n_rows + NUM_PARAMS_PER_SLICE * num_slices;
    if (array_size > MAX_ELEMENTS) {
        printf("Error: Number of entries exceeds maximum limit of %d\n", MAX_ELEMENTS);
        return;
    }

    p_custom[0] = (double)n_rows;
    p_custom[1] = (double)n_cols;
    p_custom[2] = (double)n_each;
    p_custom[3] = (double)p.surface_type;
    p_custom[4] = p.dx;
    p_custom[5] = p.dy;
    p_custom[6] = p.gx_width;
    p_custom[7] = p.gx_depth;
    p_custom[8] = p.gy_width;
    p_custom[9] = p.gy_depth;

    int slice_num, slice_idx, base_idx;

    // Calculate parameters for each slice
    for (int col = 0; col < n_cols; col++) {

        for (int row = 0; row < n_rows; row++) {

            for (int subidx = 0; subidx < n_each; subidx++) {
                slice_num = subidx + row * n_each;
                slice_idx = slice_num + col * n_rows * n_each;
                SLICE_PARAMS pslice = GetSliceParamsLinear(slice_num, col, p);
                if (subidx == 0 && col == 0) {
                    // Store u values only once per row
                    p_custom[NUM_BASE_PARAMS + row] = pslice.u;
                }

                base_idx = NUM_BASE_PARAMS + n_rows + NUM_PARAMS_PER_SLICE * slice_idx;
                p_custom[base_idx]     = pslice.alpha;
                p_custom[base_idx + 1] = pslice.beta;
                p_custom[base_idx + 2] = pslice.gamma;
                p_custom[base_idx + 3] = pslice.cv;
                p_custom[base_idx + 4] = pslice.k;
                p_custom[base_idx + 5] = pslice.zp;
                p_custom[base_idx + 6] = pslice.syx;
                p_custom[base_idx + 7] = pslice.syz;
                p_custom[base_idx + 8] = pslice.sxy;
                p_custom[base_idx + 9] = pslice.sxz;

                // no rotations about z-axis, this parameters should never be accessed
                p_custom[base_idx + 10] = 0.0;
                p_custom[base_idx + 11] = 0.0;
                p_custom[base_idx + 12] = 0.0;
            }
        }
    }
}
double CalcZpFromGamma(double gamma, double cv, double k) {
    if (cv == 0) {
        return 0;
    }
    double gamma_rad = gamma * M_PI / 180.0;
    double sing = sin(gamma_rad);
    if (cv == 0) {
        return 0;
    }
    return 1/cv * sing*sing / 2;
}

static SLICE_PARAMS GetSliceParamsAngular(int slice_num, int col_num, IMAGE_SLICER_PARAMS_ANGULAR p) {
    // Get row number and the subindex of the slice within that row
    int row_num = floor((double)slice_num / p.n_each);
    int slice_num_row = slice_num % p.n_each;
    double gamma_extra;
    SLICE_PARAMS pslice = {0};

    double row_mid = (p.n_rows % 2 == 0) ? (p.n_rows - 1) / 2.0 : p.n_rows / 2;
    double col_mid = (p.n_cols % 2 == 0) ? (p.n_cols - 1) / 2.0 : p.n_cols / 2;
    double slice_mid = (p.n_each % 2 == 0) ? (p.n_each - 1) / 2.0 : p.n_each / 2;
    double offset_row = row_num - row_mid;
    double offset_col = col_num - col_mid;

    pslice.alpha = p.alpha_cen + p.dalpha * offset_row;
    pslice.syx = p.syx_cen + p.dsyx * offset_row;
    pslice.syz = p.syz_cen + p.dsyz * offset_row;
    gamma_extra = p.gamma_offset * offset_row;
    
    pslice.beta = p.beta_cen + p.dbeta * offset_col;
    pslice.sxy = p.sxy_cen + p.dsxy * offset_col;
    pslice.sxz = p.sxz_cen + p.dsxz * offset_col;

    // Get the angles of the bottom- and top-most slices of the central row,
    // If n_rows is even, there are 2 rows straddling the x=0 center line,
    // Set the "central row" to the one above the x-axis (+y direction),
    double gamma_bot, gamma_top;
    gamma_bot = p.gamma_cen - p.dgamma * slice_mid;
    gamma_top = p.gamma_cen + p.dgamma * slice_mid;

    // Get u value a similar way
    double u_extra;
    pslice.u = p.u_cen;

    // Determine offsets in gamma. First, check whether the mode allows the extra
    // offsets to stack or not. If no, then set gamma_extra to repeat every 2 rows.
    if (p.angle_mode == 2 || p.angle_mode == 3) {
        // gamma_offset does not stack
        gamma_extra = (row_num % 2 == 0) ? -p.gamma_offset / 2.0 : p.gamma_offset / 2.0;
        u_extra = (row_num % 2 == 0) ? -p.du / 2.0 : p.du / 2.0;
    }
    else {
        u_extra = p.du * offset_row;
    }
    pslice.u += u_extra;

    // Staircase mode, need to alternate which direction gamma is incremented in
    if (p.angle_mode == 0 || p.angle_mode == 2) {
        if (row_num % 2 == 0) {
            // If the row is even, the angles are the same as the central row
            pslice.gamma = gamma_bot + slice_num_row * p.dgamma + gamma_extra;
        } else {
            // If odd, the top/bottom angles are flipped
            pslice.gamma = gamma_top - slice_num_row * p.dgamma + gamma_extra;
        }
    // Not staircase - copy the same angle pattern as the central row
    } else if (p.angle_mode == 1 || p.angle_mode == 3) {
        pslice.gamma = gamma_bot + slice_num_row * p.dgamma + gamma_extra;
    }

    pslice.zp = p.azps * CalcZpFromGamma(pslice.gamma, p.cv, p.k);

    // Constant values for all slices
    pslice.cv = p.cv;
    pslice.k  = p.k;

    return pslice;
}

static SLICE_PARAMS GetSliceParamsLinear(int slice_num, int col_num, IMAGE_SLICER_PARAMS_LINEAR p) {
    // Get row number and the subindex of the slice within that row
    int row_num = floor((double)slice_num / p.n_each);
    int slice_num_row = slice_num % p.n_each;
    SLICE_PARAMS pslice = {0};

    double row_mid = (p.n_rows % 2 == 0) ? (p.n_rows - 1) / 2.0 : p.n_rows / 2;
    double col_mid = (p.n_cols % 2 == 0) ? (p.n_cols - 1) / 2.0 : p.n_cols / 2;
    double slice_mid = (p.n_each % 2 == 0) ? (p.n_each - 1) / 2.0 : p.n_each / 2;
    double offset_row = row_num - row_mid;
    double offset_col = col_num - col_mid;

    double y0, x0, d_extra;
    double d = 0;

    y0 = p.y0_cen + p.dy0 * offset_row;
    pslice.syx = p.syx_cen + p.dsyx * offset_row;
    pslice.syz = p.syz_cen + p.dsyz * offset_row;
    d_extra = p.d_offset * offset_row;
    
    x0 = p.x0_cen + p.dx0 * offset_col;
    pslice.sxy = p.sxy_cen + p.dsxy * offset_col;
    pslice.sxz = p.sxz_cen + p.dsxz * offset_col;

    // Get the angles of the bottom- and top-most slices of the central row,
    // If n_rows is even, there are 2 rows straddling the x=0 center line,
    // Set the "central row" to the one above the x-axis (+y direction),
    double d_bot, d_top;
    d_bot = p.d_cen - p.dd * slice_mid;
    d_top = p.d_cen + p.dd * slice_mid;

    // Get u value a similar way
    double u_extra;
    pslice.u = p.u_cen;

    // Determine offsets in gamma. First, check whether the mode allows the extra
    // offsets to stack or not. If no, then set gamma_extra to repeat every 2 rows.
    if (p.angle_mode == 2 || p.angle_mode == 3) {
        // gamma_offset does not stack
        d_extra = (row_num % 2 == 0) ? -p.d_offset / 2.0 : p.d_offset / 2.0;
        u_extra = (row_num % 2 == 0) ? -p.du / 2.0 : p.du / 2.0;
    }
    else {
        u_extra = p.du * offset_row;
    }
    pslice.u += u_extra;

    // Staircase mode, need to alternate which direction gamma is incremented in
    if (p.angle_mode == 0 || p.angle_mode == 2) {
        if (row_num % 2 == 0) {
            // If the row is even, the angles are the same as the central row
            d = d_bot + slice_num_row * p.dd + d_extra;
        } else {
            // If odd, the top/bottom angles are flipped
            d = d_top - slice_num_row * p.dd + d_extra;
        }
    // Not staircase - copy the same angle pattern as the central row
    } else if (p.angle_mode == 1 || p.angle_mode == 3) {
        d = d_bot + slice_num_row * p.dd + d_extra;
    }
    // Would usually be atan(d/f) where f is the distance to the plane where it's
    // linearly spaced, but just set f=1 here. It will be linearly spaced no
    // matter the distance.
    pslice.gamma = atan( d ) * 180.0 / M_PI;

    double alpha, beta;
    Conic2DOffAxisAngle(&alpha, &beta, p.cv, p.k, x0, y0);
    pslice.alpha = alpha * 180.0 / M_PI;
    pslice.beta  = beta * 180.0 / M_PI;
    
    pslice.zp = p.azps * CalcZpFromGamma(pslice.gamma, p.cv, p.k);

    // Constant values for all slices
    pslice.cv = p.cv;
    pslice.k  = p.k;

    return pslice;
}