#define _USE_MATH_DEFINES
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "slicer_generation.h"
#include "custom_slicer_helpers.h"

/*
Slicer generation and ray tracing algorithm.

Ellen Lee
*/

void linspace(double *array, double start, double end, int n) {
    if (n <= 1) {
        array[0] = start;
        return;
    }
    double step = (end - start) / (n - 1);
    for (int i = 0; i < n; i++) {
        array[i] = start + i * step;
    }
}

// Validate image slicer parameters, modifying illegal parameters as needed.
int ValidateSlicerParams(IMAGE_SLICER_PARAMS *p) {
    // Keep track of whether we had to change any parameters
    int is_valid = 1;

    // Do not touch the custom flag!

    if (!(p->surface_type==0 || p->surface_type==1)) { p->surface_type=0; is_valid = 0; }
    if (p->n_cols < 1) { p->n_cols = 1; is_valid = 0; }
    if (p->n_rows < 1) { p->n_rows = 1; is_valid = 0; }
    if (p->n_each < 1) { p->n_each = 1; is_valid = 0; }
    if (!(p->angle_mode==0 || p->angle_mode==1 || p->angle_mode==2 || p->angle_mode==3)){
        p->angle_mode = 0; is_valid = 0;
        }

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

// Checks whether every member in each of the IMAGER_SLICER_PARAMS structs are
// equivalent
int IsParametersEqual(IMAGE_SLICER_PARAMS p1, IMAGE_SLICER_PARAMS p2) {
   if (
        p1.custom == p2.custom &&
        p1.surface_type == p2.surface_type &&
        p1.n_each == p2.n_each &&
        p1.n_rows == p2.n_rows &&
        p1.n_cols == p2.n_cols &&
        p1.angle_mode == p2.angle_mode &&

        p1.dalpha == p2.dalpha &&
        p1.dbeta == p2.dbeta &&
        p1.dgamma == p2.dgamma &&
        p1.gamma_offset == p2.gamma_offset &&

        p1.dzps == p2.dzps &&
        p1.dzp_col == p2.dzp_col &&
        p1.dzp_row == p2.dzp_row &&
        p1.dsyx == p2.dsyx &&
        p1.dsyz == p2.dsyz &&
        p1.dsxy == p2.dsxy &&
        p1.dsxz == p2.dsxz &&
        p1.du == p2.du &&

        p1.alpha_cen == p2.alpha_cen &&
        p1.beta_cen == p2.beta_cen &&
        p1.gamma_cen == p2.gamma_cen &&
        p1.zps_cen == p2.zps_cen &&
        p1.zp_cen == p2.zp_cen &&
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


// The struct p is needed to determine the boundaries betwen slices and overall
// size of the image slicer. Parameters relating to the angles, curvature, or conic
// of any individual slice will be ignored in lieu of the custom slice parameters
// if the custom flag is enabled.
//
// The important parameters (row and column numbers, slice dimensions, etc.) are
// read in from the text file. Store those parameters in a struct with this function.
IMAGE_SLICER_PARAMS MakeSlicerParamsFromCustom(double p_custom[]) {
    IMAGE_SLICER_PARAMS p;
    p.custom = 1;
    p.n_each = 1;
    p.n_rows = (int) p_custom[0];
    p.n_cols = (int) p_custom[1];
    p.surface_type = (int) p_custom[2];
    p.dx = p_custom[3];
    p.dy = p_custom[4];
    p.gx_width = p_custom[5];
    p.gx_depth = p_custom[6];
    p.gy_width = p_custom[7];
    p.gy_depth = p_custom[8];

    // The rest of these will not be touched, set to zero
    p.cv = 0;
    p.k = 0;
    p.angle_mode = 0;
    p.dalpha = 0;
    p.dbeta = 0;
    p.dgamma = 0;
    p.gamma_offset = 0;
    p.dzps = 0;
    p.dzp_col = 0;
    p.dzp_row = 0;
    p.dsyx = 0;
    p.dsyz = 0;
    p.dsxy = 0;
    p.dsxz = 0;
    p.du = 0;
    p.alpha_cen = 0;
    p.beta_cen = 0;
    p.gamma_cen = 0;
    p.zps_cen = 0;
    p.zp_cen = 0;
    p.syx_cen = 0;
    p.syz_cen = 0;
    p.sxy_cen = 0;
    p.sxz_cen = 0;
    p.u_cen = 0;

    return p;
}

void GetSurfaceFuncs(SAG_FUNC *sag_func, TRANSFER_DIST_FUNC *transfer_dist_func,
SURF_NORMAL_FUNC *surf_normal_func, CRITICAL_XY_FUNC *critical_xy_func, SLICE_PARAMS pslice, IMAGE_SLICER_PARAMS p) {
    if (pslice.cv == 0) {
        *sag_func = &TiltedPlaneSag;
        *transfer_dist_func = &TiltedPlaneTransfer;
        *surf_normal_func = &TiltedPlaneSurfaceNormal;
        *critical_xy_func = &TiltedPlaneCriticalXY;
    }
    else {
        *sag_func = &Conic2DSag;
        *transfer_dist_func = &Conic2DTransfer;
        *surf_normal_func = &Conic2DSurfaceNormal;
        *critical_xy_func = &Conic2DCriticalXY;
    }

}

// Sag generation and helpers for determining the parameters of individual slices

// The size of the image slicer along x is determined by the number of columns and
// the width of the gaps between them. Ditto for y except with the total number
// of slices within a column.
void GetSlicerSize(double *xsize, double *ysize, IMAGE_SLICER_PARAMS p) {
    int n_slices = p.n_each * p.n_rows;
    *ysize = n_slices * p.dy + (n_slices - 1) * p.gy_width;
    *xsize = p.n_cols * p.dx + (p.n_cols - 1) * p.gx_width;
}

// Column, slices indices of the slice determine its angle gamma. The indices are
// determined by the given x, y values for the image slicer parameters.
//
// Indexing for columns goes in the +x direction and slices goes in the +y direction.
// Looking at the image slicer face-on (facing the +z direction), indices (0, 0)
// correspond to the bottom left of the image slicer.
// 
// Gaps in the +x, +y directions are included as part of the index, e.g., slice 0
// includes the gap above the slice and column 0 includes the gap to the left of
// the column. The computation for the sag is cut off before the topmost and leftmost
// gaps, so in practice gaps can only appear between slices.
void GetSlicerIndex(int *col_num, int *slice_num, double x, double y, IMAGE_SLICER_PARAMS p, double p_custom[]) {
    double xsize, ysize;
    GetSlicerSize(&xsize, &ysize, p);
    *slice_num = (int) floor((y + ysize / 2) / (p.dy + p.gy_width));
    int row_num = *slice_num / p.n_each;
    // Column index depends on how much the row is shifted by (u)
    double u = GetUForRow(row_num, p, p_custom);
    *col_num = (int) floor((x - u + xsize / 2) / (p.dx + p.gx_width));
}


// Checks whether a point (x, y) is inside a gap. This is useful for sag and ray
// tracing computations.
void IsInsideSlicerGap(int *in_xgap, int *in_ygap, double x, double y, IMAGE_SLICER_PARAMS p, double p_custom[]) {
    double xsize, ysize;
    int col_num, slice_num;
    GetSlicerSize(&xsize, &ysize, p);
    GetSlicerIndex(&col_num, &slice_num, x, y, p, p_custom);

    // Check vertical (y) gap
    double ygap_bot = (slice_num + 1) * p.dy + slice_num * p.gy_width - ysize / 2;
    double ygap_top = ygap_bot + p.gy_width;
    *in_ygap = (y > ygap_bot && y <= ygap_top);

    // Check horizontal (x) gap with u shift
    int row_num = slice_num / p.n_each;
    double u = GetUForRow(row_num, p, p_custom);
    double xgap_left = (col_num + 1) * p.dx + col_num * p.gx_width - xsize / 2 + u;
    double xgap_right = xgap_left + p.gx_width;
    *in_xgap = (x > xgap_left && x <= xgap_right);
}

// If there are an even number of slices in a column, the "central" slice used
// for the paraxial ray trace is ambiguous. In this case we should allow the user
// to specify whether they want to slice above or below the center line. Same
// goes for the columns. The parameters p.active_x and p.active_y exist for this
// purpose.
void GetParaxialSliceIndex(int *col_num, int *slice_num, int active_x, int active_y, IMAGE_SLICER_PARAMS p) {
    *col_num = p.n_cols / 2;
    if (p.n_cols % 2 == 0 && active_x) (*col_num)++;

    int n_slices = p.n_each * p.n_rows;
    *slice_num = n_slices / 2;
    if (n_slices % 2 == 0 && active_y) (*slice_num)++;
}

// Computes the maximum and minimum row offsets for the image slicer. umin and umax
// are needed along with the global extrema of the sag to compute the ray trace.
void GetMinMaxU(double *umin, double *umax, IMAGE_SLICER_PARAMS p, double p_custom[]) {
    int i;
    
    // If using custom slice parameters, brute force search
    if (p.custom) {
        // double u_all[p.n_rows];  //u_all needs to be a constant length!
        // for (i = 0; i < p.n_rows; ++i) {
        //     u_all[i] = p_custom[9 + i];
        // }
        // *umin = u_all[0];
        // *umax = u_all[0];
        // for (i = 1; i < p.n_rows; ++i) {
        //     if (u_all[i] < *umin) *umin = u_all[i];
        //     if (u_all[i] > *umax) *umax = u_all[i];
        // }
        *umin = p_custom[9];
        *umax = p_custom[9];
        for (i = 1; i < p.n_rows; ++i) {
            if (p_custom[9 + i] < *umin) *umin = p_custom[9 + i];
            if (p_custom[9 + i] > *umax) *umax = p_custom[9 + i];
        }
        return;
    }

    int num_slices = p.n_each * p.n_rows;
    if (num_slices < 2) {
        // Only one slice, so no need to compare
        *umin = p.u_cen;
        *umax = p.u_cen;
        return;
    }

    // u does not alternate. Minimum and maximum must be at the first and last rows.
    if (p.angle_mode == 0 || p.angle_mode == 1) {
        *umin = get_u_for_row(0, p, p_custom);
        *umax = get_u_for_row(p.n_rows - 1, p, p_custom);
    // u alternates odd/even rows. 
    } else if (p.angle_mode == 2 || p.angle_mode == 3) {
        *umin = get_u_for_row(0, p, p_custom);
        *umax = get_u_for_row(1, p, p_custom);
    }

    // umin and umax need to be swapped if the change in u is negative.
    if (p.du < 0) {
        double tmp = *umin;
        *umin = *umax;
        *umax = tmp;
    }
}

// If p.custom, then use the slice parameters provided in the array p_custom.
// Otherwise, this array is ignored.
SLICE_PARAMS GetSliceParams(int slice_num, int col_num, IMAGE_SLICER_PARAMS p,
    double p_custom[]) {
    if (p.custom) {
        return GetSliceParamsCustom(slice_num, col_num, p_custom);
    }
    return GetSliceParamsStandard(slice_num, col_num, p);
}

// Calculates u for the given row index. This is computed in a separate function
// as it is needed to calculate the column index of a slice.
double GetUForRow(int row_num, IMAGE_SLICER_PARAMS p, double p_custom[]) {
    if (p.custom) {
        return p_custom[9 + row_num];
    }

    double u_extra;
    // Similar to GetSliceParamsStandard
    if (p.n_rows % 2 == 0) {
        u_extra = p.du * (row_num - (p.n_rows - 1) / 2.0);
    } else {
        u_extra = p.du * (row_num - p.n_rows / 2);
    }
    // angle_mode also determines whether u should stack because the way gamma is
    // constructed on the image slicer should determine how u is constructed on the
    // subpupil mirrors. This somewhat restricts which surfaces can be defined in
    // standard mode, but if this is insufficient the user can use custom mode instead.
    if (p.angle_mode == 2 || p.angle_mode == 3) {
        if (row_num % 2 == 0) {
            u_extra = -p.du / 2.0;
        } else {
            u_extra = p.du / 2.0;
        }
    }

    return p.u_cen + u_extra;
}

// For a given column, row indices determine the angles alpha, beta, and gamma.
// alpha and beta determine the behavior of one "section" of the image slicer.
// gamma defines the rotation of each slice within that section. Essentially,
// alpha results in rows of pupil images (e.g., IRTF/SPECTRE) and beta gives
// different columns (e.g., VLT/MUSE).
// 
// The p.angle_mode parameter allows the user to specify how the angles gamma switch
// between rows. If the mode is 1, then gamma is incremented the same way in every
// section. If 0, then sign of dgamma is flipped so that pattern is alternated
// like a staircase.
SLICE_PARAMS GetSliceParamsStandard(int slice_num, int col_num, IMAGE_SLICER_PARAMS p) {
    // Get row number and the subindex of the slice within that row
    int row_num = slice_num / p.n_each;
    int slice_num_row = slice_num % p.n_each;
    double gamma_extra;
    SLICE_PARAMS pslice;

    // Get u value
    pslice.u = get_u_for_row(row_num, p, NULL);  // assuming p_custom not needed

    double row_mid = (p.n_rows % 2 == 0) ? (p.n_rows - 1) / 2.0 : p.n_rows / 2;
    double col_mid = (p.n_cols % 2 == 0) ? (p.n_cols - 1) / 2.0 : p.n_cols / 2;
    double slice_mid = (p.n_each % 2 == 0) ? (p.n_each - 1) / 2.0 : p.n_each / 2;
    double offset_row = row_num - row_mid;
    double offset_col = col_num - col_mid;
    double offset_slice = slice_num_row - slice_mid;

    pslice.alpha = p.alpha_cen + p.dalpha * offset_row;
    pslice.syx = p.syx_cen + p.dsyx * offset_row;
    pslice.syz = p.syz_cen + p.dsyz * offset_row;
    pslice.zp  = p.zp_cen + p.dzp_row * offset_row;
    gamma_extra = p.gamma_offset * offset_row;
    
    pslice.beta = p.beta_cen + p.dbeta * offset_col;
    pslice.sxy = p.sxy_cen + p.dsxy * offset_col;
    pslice.sxz = p.sxz_cen + p.dsxz * offset_col;
    pslice.zp += p.dzp_col * offset_col;

    pslice.zp += p.zps_cen + p.dzps * offset_slice;

    // Get the angles of the bottom- and top-most slices of the central row,
    // If n_rows is even, there are 2 rows straddling the x=0 center line,
    // Set the "central row" to the one above the x-axis (+y direction),
    double gamma_bot, gamma_top;
    if (p.n_each % 2 == 0) {
        gamma_bot = p.gamma_cen - p.dgamma * (p.n_each - 1) / 2.0;
        gamma_top = p.gamma_cen + p.dgamma * (p.n_each - 1) / 2.0;
    } else {
        gamma_bot = p.gamma_cen - p.dgamma * (p.n_each / 2);
        gamma_top = p.gamma_cen + p.dgamma * (p.n_each / 2);
    }

    // Determine offsets in gamma. First, check whether the mode allows the extra
    // offsets to stack or not. If no, then set gamma_extra to repeat every 2 rows.
    if (p.angle_mode == 2 || p.angle_mode == 3) {
        // gamma_offset does not stack
        gamma_extra = (row_num % 2 == 0) ? -p.gamma_offset / 2.0 : p.gamma_offset / 2.0;
    }

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

    // Constant values for all slices
    pslice.cv = p.cv;
    pslice.k  = p.k;
}

SLICE_PARAMS GetSliceParamsCustom(int slice_num, int col_num, double p_custom[]) {

    int n_slices_per_col = (int) p_custom[0];
    int start_idx = 9 + 5 * (col_num * n_slices_per_col + slice_num);

    SLICE_PARAMS pslice;
    pslice.alpha = p_custom[start_idx];
    pslice.beta = p_custom[start_idx + 1];
    pslice.gamma = p_custom[start_idx + 2];
    pslice.cv = p_custom[start_idx + 3];
    pslice.k = p_custom[start_idx + 4];

    return pslice;
}


// The sag of the image slicer can be found by
double ImageSlicerSag(double x, double y, IMAGE_SLICER_PARAMS p, double p_custom[]) {
    // Check if out of bounds
    int col_num, slice_num;
    GetSlicerIndex(&col_num, &slice_num, x, y, p, p_custom);
    int row_num = slice_num / p.n_each;
    if (row_num < 0 || row_num >= p.n_rows || col_num < 0 || col_num >= p.n_cols) {
        return NAN;
    }

    // If inside a gap, return the gap depth unless at the outer edge (in which
    // case it should be considered out of bounds)
    int in_xgap, in_ygap;
    IsInsideSlicerGap(&in_xgap, &in_ygap, x, y, p, p_custom);
    if (in_xgap) {
        if (col_num == p.n_cols - 1) {
            return NAN;
        }
        return p.gx_depth;
    }
    if (in_ygap) {
        if (slice_num == p.n_each * p.n_rows - 1) {
            return NAN;
        }
        return p.gy_depth;
    }

    // Inside a slice. From the slice number, determine the angles
    SLICE_PARAMS pslice = GetSliceParams(slice_num, col_num, p, p_custom);
    SAG_FUNC sag_func;
    TRANSFER_DIST_FUNC transfer_dist_func;
    SURF_NORMAL_FUNC surf_normal_func;
    CRITICAL_XY_FUNC critical_xy_func;
    GetSurfaceFuncs(&sag_func, &transfer_dist_func, &surf_normal_func, &critical_xy_func, pslice, p);
    return sag_func(x, y, pslice);
}


// Ray tracing

double FindBoundedSliceExtremum(double x0, double y0, int mode, IMAGE_SLICER_PARAMS p, double p_custom[]) {
    // If x0, y0 are inside gaps then we don't need to go further
    int in_xgap, in_ygap;
    IsInsideSlicerGap(&in_xgap, &in_ygap, x0, y0, p, p_custom);
    if (in_xgap) return p.gx_depth;
    if (in_ygap) return p.gy_depth;

    // Not inside a gap. We need to compute the sag at the bounds and critical
    // point(s) of the slice.
    // First, compute the bounds of the current slice
    int col_num, slice_num;
    GetSlicerIndex(&col_num, &slice_num, x0, y0, p, p_custom);
    double xsize, ysize;
    GetSlicerSize(&xsize, &ysize, p);
    SLICE_PARAMS pslice = GetSliceParams(slice_num, col_num, p, p_custom);
    SAG_FUNC sag_func;
    TRANSFER_DIST_FUNC transfer_dist_func;
    SURF_NORMAL_FUNC surf_normal_func;
    CRITICAL_XY_FUNC critical_xy_func;
    GetSurfaceFuncs(&sag_func, &transfer_dist_func, &surf_normal_func, &critical_xy_func, pslice, p);

    double u = GetUForRow(slice_num / p.n_each, p, p_custom);
    double xlo = col_num * (p.dx + p.gx_width) - xsize / 2 + u;
    double xhi = xlo + p.dx;
    double ylo = slice_num * (p.dy + p.gy_width) - ysize / 2;
    double yhi = ylo + p.dy;

    // Compute critical point
    double xc, yc;
    critical_xy_func(&xc, &yc, pslice);

    // There are up to 5 points to compare depending on whether the critical point
    // is within bounds.
    double zsolns[5];
    zsolns[0] = sag_func(xlo, ylo, pslice);
    zsolns[1] = sag_func(xlo, yhi, pslice);
    zsolns[2] = sag_func(xhi, ylo, pslice);
    zsolns[3] = sag_func(xhi, yhi, pslice);

    int n_compare = 4; // Number of elements to compare so far

    // Check whether critical points are in bounds. If yes, compute the sag an
    // add to the array of points to compare.
    if (yc >= ylo && yc <= yhi && xc >= xlo && xc <= xhi) {
        zsolns[4] = sag_func(xc, yc, pslice);
        n_compare++;
    }
    
    // Compare potential solutions to get the maximum or minimum
    double result = zsolns[0];
    for (int i = 1; i < n_compare; i++) {
        if (mode) {
            if (zsolns[i] > result) result = zsolns[i]; // max
        } else {
            if (zsolns[i] < result) result = zsolns[i]; // min
        }
    }
    return result;
}


void FindSlicerGlobalExtrema(double *zmin, double *zmax, IMAGE_SLICER_PARAMS p, double p_custom[]) {
    // Determine which slice the global extrema are on by roughly sampling the
    // entire image slicer
    double xsize, ysize;
    GetSlicerSize(&xsize, &ysize, p);
    int nx = p.n_cols * 6;
    int ny = p.n_rows * p.n_each * 8;

    // Safeguard in case the user attempts to initialize an obscenely large number
    // of slices
    if (nx > 30000 || ny > 40000) {
        nx = 30000; // Implies 5,000 columns
        ny = 40000; // 5,000 slices per column...
    }

    // Generally the number of points shouldn't be a problem but dynamically allocate
    // memory for xpts and ypts just in case
    double *xpts = (double *)malloc(nx * sizeof(double));
    double *ypts = (double *)malloc(ny * sizeof(double));
    if (xpts == NULL || ypts == NULL) {
        // If memory allocation fails, print an error and give up
        fprintf(stderr, "Memory allocation failed.\n");
        return;
    }
    linspace(xpts, -xsize/2, xsize/2, nx);
    linspace(ypts, -ysize/2, ysize/2, ny);

    // Evaluate the grid, keeping track of the maximum and minimum
    double x0_max = 0, y0_max = 0, z0_max = -INFINITY;
    double x0_min = 0, y0_min = 0, z0_min = INFINITY;
    double z;
    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            z = ImageSlicerSag(xpts[i], ypts[j], p, p_custom);
            if (z > z0_max) {
                x0_max = xpts[i];
                y0_max = ypts[j];
                z0_max = z;
            } else if (z < z0_min) {
                x0_min = xpts[i];
                y0_min = ypts[j];
                z0_min = z;
            }
        }
    }
    
    free(xpts); free(ypts);

    // Use the estimated maximum and minimum to find the exact values
    *zmin = FindBoundedSliceExtremum(x0_min, y0_min, 0, p, p_custom);
    *zmax = FindBoundedSliceExtremum(x0_max, y0_max, 1, p, p_custom);
}

double TransferFunction(double t, RAY_IN ray_in, IMAGE_SLICER_PARAMS p, double p_custom[]) {
    double xs = ray_in.xt + t*ray_in.l;
    double ys = ray_in.yt + t*ray_in.m;
    double zs = t*ray_in.n;
    double sag = ImageSlicerSag(xs, ys, p, p_custom);
    return sag - zs;
}

RAY_BOUNDS GetRayBounds(RAY_IN ray_in, double umin, double umax, double zmin, double zmax,
    IMAGE_SLICER_PARAMS p, void *p_custom) {

    RAY_BOUNDS bounds;

    double xsize, ysize;
    double tmin, tmax;
    GetSlicerSize(&xsize, &ysize, p);

    if (fabs(ray_in.n) < 1e-13) {
        // Ray is nearly parallel to the z-axis

        if (fabs(ray_in.l) < 1e-13) {
            // Ray moves only along y
            double ymin = -ysize / 2.0;
            double ymax =  ysize / 2.0;
            if (ray_in.m < 0) {
                double tmp = ymin;
                ymin = ymax;
                ymax = tmp;
            }
            tmin = (ymin - ray_in.yt) / ray_in.m;
            bounds.xmin = ray_in.xt + tmin * ray_in.l;
            bounds.ymin = ymin;

            tmax = (ymax - ray_in.yt) / ray_in.m;
            bounds.xmax = ray_in.xt + tmax * ray_in.l;
            bounds.ymax = ymax;

        } else {
            // Ray moves in x (or diagonally)
            double xmin = -xsize / 2.0 + umin;
            double xmax =  xsize / 2.0 + umax;
            if (ray_in.l < 0) {
                double tmp = xmin;
                xmin = xmax;
                xmax = tmp;
            }
            tmin = (xmin - ray_in.xt) / ray_in.l;
            bounds.xmin = xmin;
            bounds.ymin = ray_in.yt + tmin * ray_in.m;

            tmax = (xmax - ray_in.xt) / ray_in.l;
            bounds.xmax = xmax;
            bounds.ymax = ray_in.yt + tmax * ray_in.m;
        }

    } else {
        // Ray has z component
        tmin = zmin / ray_in.n;
        bounds.xmin = ray_in.xt + tmin * ray_in.l;
        bounds.ymin = ray_in.yt + tmin * ray_in.m;

        tmax = zmax / ray_in.n;
        bounds.xmax = ray_in.xt + tmax * ray_in.l;
        bounds.ymax = ray_in.yt + tmax * ray_in.m;
    }

    // Get slicer indices for start and end
    GetSlicerIndex(&bounds.nc_min, &bounds.ns_min, bounds.xmin, bounds.ymin, p,
        p_custom);

    GetSlicerIndex(&bounds.nc_max, &bounds.ns_max, bounds.xmax, bounds.ymax, p,
        p_custom);

    bounds.sgnc = (bounds.xmin <= bounds.xmax) ? 1 : -1;
    bounds.sgns = (bounds.ymin <= bounds.ymax) ? 1 : -1;

    return bounds;
}

// Checks if the ray is bounds at least some of the time.
int IsRayInBounds(int nc_min, int ns_min, int nc_max, int ns_max, double umax, double umin, IMAGE_SLICER_PARAMS p) {
    int dcol = ceil( (umax - umin) / p.dx);
    int col_min = -dcol;
    int col_max = p.n_cols + dcol;

    if ( (nc_min < col_min && nc_max < col_min) || (nc_min >= col_max && nc_max >= col_max)) {
        return 0;  // x-value of the ray is too high or low
    }

    int n_sperc = p.n_each * p.n_rows;  // number of slices per column
    if ( (ns_min < 0 && ns_max < 0) || (ns_min >= n_sperc && ns_max >= n_sperc) ) {
        return 0;  // y-value of the ray is too high or low
    }
    
    return 1;
}

int IsSectionInBounds(int col_num, int row_num, double umin, double umax, IMAGE_SLICER_PARAMS p) {
    if (row_num < 0 || row_num >= p.n_rows) {
        return 0;
    }

    // Column indices require us to consider how each row is shifted
    int dcol = (int)ceil((umax - umin) / p.dx);
    if (col_num < -dcol || col_num >= p.n_cols + dcol) {
        return 0;
    }

    return 1;
}

void CheckSliceSolution(RAY_OUT *ray_out, double tol, RAY_IN ray_in, int ns_test, int nc_test,
    IMAGE_SLICER_PARAMS p, void *p_custom) {

    *ray_out = (RAY_OUT){NAN, NAN, NAN, NAN, NAN, NAN, NAN};

    double xt = ray_in.xt; double yt = ray_in.yt;
    double l = ray_in.l; double m = ray_in.m; double n = ray_in.n;

    // Check if this slice is a solution. Get params for this slice and
    // compute the transfer distance.
    SLICE_PARAMS pslice = GetSliceParams(ns_test, nc_test, p, p_custom);
    SAG_FUNC sag_func;
    TRANSFER_DIST_FUNC transfer_dist_func;
    SURF_NORMAL_FUNC surf_normal_func;
    CRITICAL_XY_FUNC critical_xy_func;
    GetSurfaceFuncs(&sag_func, &transfer_dist_func, &surf_normal_func, &critical_xy_func, pslice, p);

    double t = transfer_dist_func(xt, yt, l, m, n, pslice);
    double result = TransferFunction(t, ray_in, p, p_custom);

    // Is this a valid zero to the transfer function?
    if (fabs(result) < tol) {
        // Yes - found a solution!
        double xs = xt + t * l;
        double ys = yt + t * m;
        double zs = t * n;

        double ln, mn, nn;
        surf_normal_func(&ln, &mn, &nn, xs, ys, pslice, 1);

        // WAIT - Is the solution inside of a gap? If we're unlucky and zs is
        // equal to the gap depth then this may be the case.
        int in_xgap, in_ygap;
        IsInsideSlicerGap(&in_xgap, &in_ygap, xs, ys, p, p_custom);
        if (in_xgap || in_ygap) {
            return;
        }

        // Not in a gap. This slice is the solution.
        ray_out->xs = xs; ray_out->ys = ys; ray_out->zs = zs;
        ray_out->t  = t; ray_out->ln = ln; ray_out->mn = mn; ray_out->nn = nn;
    }

    return;
}
void CheckYWallCollision(RAY_OUT *ray_out, RAY_IN ray_in, int ns_test, int nc_test, int sgns,
    IMAGE_SLICER_PARAMS p, void *p_custom){

    *ray_out = (RAY_OUT){NAN, NAN, NAN, NAN, NAN, NAN, NAN};
    double xt = ray_in.xt; double yt = ray_in.yt;
    double l = ray_in.l; double m = ray_in.m; double n = ray_in.n;

    // Get y-dimension size of the slicer
    double xsize, ysize;
    GetSlicerSize(&xsize, &ysize, p);

    // Get sag function for the current column
    SLICE_PARAMS pslice = GetSliceParams(ns_test, nc_test, p, p_custom);
    SAG_FUNC sag_func;
    TRANSFER_DIST_FUNC transfer_dist_func;
    SURF_NORMAL_FUNC surf_normal_func;
    CRITICAL_XY_FUNC critical_xy_func;
    GetSurfaceFuncs(&sag_func, &transfer_dist_func, &surf_normal_func, &critical_xy_func, pslice, p);

    if (fabs(m) > 1e-13) {
        // Near wall
        double ynear = ((1 + sgns) / 2.0) * p.dy + ns_test * (p.dy + p.gy_width) - ysize / 2.0;
        double tnear = (ynear - yt) / m;
        double xnear = xt + tnear * l;
        double znear = tnear * n;

        // Far wall
        double yfar = ynear + sgns * p.gy_width;
        double tfar = (yfar - yt) / m;
        double xfar = xt + tfar * l;
        double zfar = tfar * n;

        // Sag of near slice
        SLICE_PARAMS pslice_near = GetSliceParams(ns_test, nc_test, p, p_custom);
        double znear_slice = sag_func(xnear, ynear, pslice_near);

        // Sag of far slice
        SLICE_PARAMS pslice_far = GetSliceParams(ns_test + sgns, nc_test, p, p_custom);
        double zfar_slice = sag_func(xfar, yfar, pslice_far);

        // Check near wall collision
        if ((znear_slice > p.gy_depth && p.gy_width > 0) &&
            (znear <= znear_slice && znear >= p.gy_depth)) {
            ray_out->xs = xnear; ray_out->ys = ynear; ray_out->zs = znear;
            ray_out->t = tnear; ray_out->ln = 0.0; ray_out->mn = -1.0 * sgns; ray_out->nn = 0.0;
            return;
        }

        // Check far wall collision
        double zcompare = (p.gy_width == 0) ? znear_slice : p.gy_depth;
        if ((zfar >= zfar_slice && zfar <= zcompare) ||
            (zfar <= zfar_slice && zfar >= zcompare)) {
            ray_out->xs = xfar; ray_out->ys = yfar; ray_out->zs = zfar;
            ray_out->t = tfar; ray_out->ln = 0.0; ray_out->mn = -1.0 * sgns; ray_out->nn = 0.0;
            return;
        }

        // gap
        if ((zfar >= p.gy_depth) && fabs(n) > 1e-13) {
            double t = p.gx_depth / n;
            ray_out->xs = xt + t * l; ray_out->ys = yt + t * m; ray_out->zs = p.gx_depth;
            ray_out->t  = t; ray_out->ln = 0.0; ray_out->mn = 0.0; ray_out->nn = -1.0;
            return;
        }
    }

    // Ray missed
    return;
}

void CheckXWallCollision(RAY_OUT *ray_out, RAY_IN ray_in, int ns_test, int nc_test, int sgnc,
    IMAGE_SLICER_PARAMS p, void *p_custom) {

    *ray_out = (RAY_OUT){NAN, NAN, NAN, NAN, NAN, NAN, NAN};

    double xt = ray_in.xt; double yt = ray_in.yt;
    double l = ray_in.l; double m = ray_in.m; double n = ray_in.n;

    double xsize, ysize;
    GetSlicerSize(&xsize, &ysize, p);

    SLICE_PARAMS pslice = GetSliceParams(ns_test, nc_test, p, p_custom);
    SAG_FUNC sag_func;
    TRANSFER_DIST_FUNC transfer_dist_func;
    SURF_NORMAL_FUNC surf_normal_func;
    CRITICAL_XY_FUNC critical_xy_func;
    GetSurfaceFuncs(&sag_func, &transfer_dist_func, &surf_normal_func, &critical_xy_func, pslice, p);

    if (fabs(l) > 1e-13) {
        double u = GetUForRow(ns_test / p.n_each, p, p_custom);

        // Ray coordinates on near wall
        double xnear = ((1 + sgnc) / 2.0) * p.dx + nc_test * (p.dx + p.gx_width) - xsize / 2.0 + u;
        double tnear = (xnear - xt) / l;
        double ynear = yt + tnear * m;
        double znear = tnear * n;

        // Ray coordinates on far wall
        double xfar = xnear + sgnc * p.gx_width;
        double tfar = (xfar - xt) / l;
        double yfar = yt + tfar * m;
        double zfar = tfar * n;

        // Sag of near and far slices
        SLICE_PARAMS pslice_near = GetSliceParams(ns_test, nc_test, p, p_custom);
        double znear_slice = sag_func(xnear, ynear, pslice_near);

        SLICE_PARAMS pslice_far = GetSliceParams(ns_test, nc_test + 1 * sgnc, p, p_custom);
        double zfar_slice = sag_func(xfar, yfar, pslice_far);

        // Check near wall collision
        if ((znear_slice > p.gx_depth && p.gx_width > 0) &&
            (znear <= znear_slice && znear >= p.gx_depth)) {
            ray_out->xs = xnear; ray_out->ys = ynear; ray_out->zs = znear;
            ray_out->t  = tnear; ray_out->ln = -1.0 * sgnc; ray_out->mn = 0.0; ray_out->nn = 0.0;
            return;
        }

        // Check far wall collision
        double zcompare = (p.gx_width == 0) ? znear_slice : p.gx_depth;
        if ((zfar >= zfar_slice && zfar <= zcompare) ||
            (zfar <= zfar_slice && zfar >= zcompare)) {
            ray_out->xs = xfar; ray_out->ys = yfar; ray_out->zs = zfar;
            ray_out->t  = tfar; ray_out->ln = -1.0 * sgnc; ray_out->mn = 0.0; ray_out->nn = 0.0;
            return;
        }

        // gap
        if ((zfar >= p.gx_depth) && fabs(n) > 1e-13) {
            double t = p.gy_depth / n;
            ray_out->xs = xt + t * l; ray_out->ys = yt + t * m; ray_out->zs = p.gy_depth;
            ray_out->t  = t; ray_out->ln = 0.0; ray_out->mn = 0.0; ray_out->nn = -1.0;
            return;
        }
    }

    return;
}

void CalcNextCoords(double *x_next, double *y_next, int *code, RAY_IN ray_in, int sgnc, int sgns, int nc_test, int nr_test,
    double x_test, double y_test, double xmax, double ymax, IMAGE_SLICER_PARAMS p,
    void *p_custom) {

    double l = ray_in.l;
    double m = ray_in.m;

    // Incoming ray is straight-on in z-direction
    if (fabs(l) < 1e-13 && fabs(m) < 1e-13) {
        *x_next = xmax;
        *y_next = ymax;
        *code = 0;
        return;
    }

    // Avoid division by zero
    if (fabs(l) < 1e-13) l = 1e-13;
    if (fabs(m) < 1e-13) m = 1e-13;

    double xsize, ysize;
    GetSlicerSize(&xsize, &ysize, p);

    double u = GetUForRow(nr_test, p, p_custom);

    // Compute next column intersection
    double x_nextcol = (nc_test + (1 + sgnc) / 2.0) * (p.dx + p.gx_width) - xsize / 2.0 + u;
    double y_nextcol = y_test + (x_nextcol - x_test) * m / l;

    // Compute next row intersection
    double y_nextrow = (nr_test + (1 + sgns) / 2.0) * p.n_each * (p.dy + p.gy_width) - ysize / 2.0;
    double x_nextrow = x_test + (y_nextrow - y_test) * l / m;

    // Compute distances to next options
    double dtocol = hypot(x_nextcol - x_test, y_nextcol - y_test);
    double dtorow = hypot(x_nextrow - x_test, y_nextrow - y_test);
    double dtomax = hypot(xmax - x_test, ymax - y_test);

    if (dtomax <= dtocol && dtomax <= dtorow) {
        *x_next = xmax;
        *y_next = ymax;
        *code = 0;
    } else if (dtocol <= dtorow && dtocol <= dtomax) {
        *x_next = x_nextcol;
        *y_next = y_nextcol;
        *code = 1;
    } else {
        *x_next = x_nextrow;
        *y_next = y_nextrow;
        *code = 2;
    }
}

void CalcNumSlicesToCheck(int sgnc, int sgns, int nc_test, int nr_test,
                          double x_test, double y_test,
                          double x_next, double y_next, int code,
                          IMAGE_SLICER_PARAMS p, void *p_custom,
                          int *n_stocheck, int *nc_new, int *nr_new) {
    int ns1 = 0, ns2 = 0;

    // Get slice indices (assuming GetSlicerIndex returns column and slice indices)
    int dummy_nc;

    GetSlicerIndex(&dummy_nc, &ns1, x_test, y_test, p, p_custom);
    GetSlicerIndex(&dummy_nc, &ns2,
                   x_next - 1e-13 * sgnc,
                   y_next - 1e-13 * sgns,
                   p, p_custom);

    int diff = abs(ns2 - ns1) + 1;
    *n_stocheck = (diff < p.n_each) ? diff : p.n_each;

    switch(code) {
        case 0: // stay in this section
            *nc_new = nc_test;
            *nr_new = nr_test;
            break;
        case 1: // next column
            *nc_new = nc_test + sgnc;
            *nr_new = nr_test;
            break;
        case 2: // next row
            *nc_new = nc_test;
            *nr_new = nr_test + sgns;
            break;
        default:
            // fallback: stay in place
            *nc_new = nc_test;
            *nr_new = nr_test;
            break;
    }
}

void RayTraceSlicer(RAY_OUT *ray_out, RAY_IN ray_in, double zmin, double zmax, double umin, double umax,
                    int trace_walls, IMAGE_SLICER_PARAMS p, double p_custom[]) {
    // Tolerance for accepting the transfer distance as valid
    double tol = 1e-12;

    // Initialize ray_out with NaNs
    *ray_out = (RAY_OUT){NAN, NAN, NAN, NAN, NAN, NAN, NAN};

    // Get starting and ending rows and columns
    RAY_BOUNDS bounds = GetRayBounds(ray_in, umin, umax, zmin, zmax, p, p_custom);
    int nc_min = bounds.nc_min; int ns_min = bounds.ns_min;
    int nc_max = bounds.nc_max; int ns_max = bounds.ns_max;
    int sgnc = bounds.sgnc; int sgns = bounds.sgns;
    double xmin = bounds.xmin; double ymin = bounds.ymin;
    double xmax = bounds.xmax; double ymax = bounds.ymax;

    if (!IsRayInBounds(nc_min, ns_min, nc_max, ns_max, umin, umax, p)) {
        return;  // ray_out is already filled with NaNs
    }

    // Ray is in bounds at least some of the time. Start from the min col and slice
    // indices and check solutions until we hit the max
    int nc_test = nc_min;
    int ns_test = ns_min;
    int nr_test = ns_min / p.n_each;
    double x_test = xmin; double y_test = ymin;

    // If starting index isn't in bounds, set it to where it will be in bounds
    while (!IsSectionInBounds(nc_test, nr_test, umin, umax, p)) {
        double x_next, y_next;
        int code;
        CalcNextCoords(&x_next, &y_next, &code, ray_in, sgnc, sgns, nc_test, nr_test, x_test, y_test,
                       xmax, ymax, p, p_custom);
        int n_stocheck;
        CalcNumSlicesToCheck(sgnc, sgns, nc_test, nr_test, x_test, y_test, x_next, y_next, code,
                             p, p_custom, &n_stocheck, &nc_test, &nr_test);
        x_test = x_next;
        y_test = y_next;
    }

    // Get slice index from test point
    int junk;
    GetSlicerIndex(&junk, &ns_test, x_test, y_test, p, p_custom);

    int nc_new = -1;
    int nr_new = -1;
    int is_same_section = 0;

    while (IsSectionInBounds(nc_test, nr_test, umin, umax, p) && !is_same_section) {
        // Number of slices to check in this section
        double x_next, y_next;
        int code;
        CalcNextCoords(&x_next, &y_next, &code, ray_in, sgnc, sgns, nc_test, nr_test, x_test, y_test,
                       xmax, ymax, p, p_custom);

        int n_stocheck;
        CalcNumSlicesToCheck(sgnc, sgns, nc_test, nr_test, x_test, y_test, x_next, y_next, code,
                             p, p_custom, &n_stocheck, &nc_new, &nr_new);
        x_test = x_next;
        y_test = y_next;

        // Iterate slices within this section
        for (int n = 0; n < n_stocheck; n++) {
            // Check whether each slice is a solution to the transfer equation
            CheckSliceSolution(ray_out, tol, ray_in, ns_test, nc_test, p, p_custom);

            if (!isnan(ray_out->t)) {
                return;  // Found a solution within a slice
            }

            // Before going to the next slice, check if there is a collision with
            // a wall. Skip if this is the last slice in the section.
            if (trace_walls && (ns_test % p.n_each != 0 && ns_test % p.n_each != p.n_each - 1)) {
                CheckYWallCollision(ray_out, ray_in, ns_test, nc_test, sgns, p, p_custom);
                if (!isnan(ray_out->t)) {
                    return;  // Found a wall collision
                }
            }

            // Increment the slice index
            ns_test += sgns;
        }

        // Before continuing, check walls between sections.
        is_same_section = (nc_new == nc_test && nr_new == nr_test);

        if (trace_walls && IsSectionInBounds(nc_new, nr_new, umin, umax, p)) {
            if (nc_new != nc_test) {
                // Section changed columns, check collision with x-wall
                CheckXWallCollision(ray_out, ray_in, ns_test, nc_test, sgnc, p, p_custom);
                if (!isnan(ray_out->t)) {
                    return;
                }
            } else if (nr_new != nr_test) {
                // Section changed rows, check collision with y-wall
                CheckYWallCollision(ray_out, ray_in, ns_test, nc_test, sgns, p, p_custom);
                if (!isnan(ray_out->t)) {
                    return;
                }
            } else {
                // Staying on the same section. Check both x- and y-walls...
                int in_xgap = 0, in_ygap = 0;
                IsInsideSlicerGap(&in_xgap, &in_ygap,
                                  xmax - 1e-13 * sgnc, y_test - 1e-13 * sgns,
                                  p, p_custom);

                if (in_xgap) {
                    CheckXWallCollision(ray_out, ray_in, ns_test, nc_test, sgnc, p, p_custom);
                    if (!isnan(ray_out->t)) {
                        return;
                    }
                } else if (in_ygap) {
                    CheckYWallCollision(ray_out, ray_in, ns_test, nc_test, sgns, p, p_custom);
                    if (!isnan(ray_out->t)) {
                        return;
                    }
                }
                // If neither in_xgap or in_ygap are true, then a slice intersection should
                // have already been found so we don't need to worry about this case.
            }
        }

        // If crossing columns, need to double check the last slice index in the
        // next section
        if (nc_new != nc_test) {
            ns_test -= sgns;
        }

        nc_test = nc_new;
        nr_test = nr_new;
    }

    // If none of that worked then this is a bizarre edge case, e.g., the direction
    // cosine is less than 1E-13 but it grazed off of a wall somehow. Sweep this ray
    // under the rug and say that it missed... ray_out will be filled with NaN values
    // at this point.
    return;  // 5 tumbleweeds roll by...
}

void ParaxialRayTraceSlicer(RAY_OUT *ray_out, double *l_out, double *m_out, double *n_out,
    RAY_IN *ray_in, double n1, double n2, int active_x, int active_y,
    IMAGE_SLICER_PARAMS p, double p_custom[]) 
{
    int nc, ns;
    GetParaxialSliceIndex(&nc, &ns, active_x, active_y, p);
    SLICE_PARAMS pslice = GetSliceParams(ns, nc, p, p_custom);

    // Transform the ray into local coordinates
    double coords[3]   = { ray_in->xt, ray_in->yt, 0.0 };
    double cosines[3]  = { ray_in->l,  ray_in->m,  ray_in->n };
    double coords_out[3], cosines_out[3];

    Conic2DTransformation(coords_out, coords, pslice, 1, 1);   // forward, translate
    Conic2DTransformation(cosines_out, cosines, pslice, 1, 0); // forward, no translate

    // Power calculation
    double power = (n2 - n1) * pslice.cv;
    double powerx, powery;
    if (p.surface_type == 0) {
        powerx = power;
        powery = power;
    } else if (p.surface_type == 1) {    // cylindrical
        powerx = power;
        powery = 0.0;
    } else {
        powerx = power;
        powery = power;
    }

    double x = coords_out[0];
    double y = coords_out[1];
    double z = coords_out[2];
    double l = cosines_out[0];
    double m = cosines_out[1];
    double n = cosines_out[2];

    if (fabs(n) < 1e-13) {
        *ray_out = (RAY_OUT){NAN,NAN,NAN,NAN,NAN,NAN,NAN};
        if (l_out) *l_out = NAN;
        if (m_out) *m_out = NAN;
        if (n_out) *n_out = NAN;
        return;
    }

    // Convert to slopes
    l /= n;
    m /= n;

    l = (n1 * l - x * powerx) / n2;
    m = (n1 * m - y * powery) / n2;

    // Convert back to direction cosines and normalize
    n = sqrt(1.0 / (1.0 + l*l + m*m));
    l *= n;
    m *= n;

    double coords_par[3]   = { x, y, z };
    double cosines_par[3]  = { l, m, n };

    // Transform back to global coordinates
    double coords_back[3], cosines_back[3];
    Conic2DTransformation(coords_back, coords_par, pslice, -1, 1);  // inverse, translate
    Conic2DTransformation(cosines_back, cosines_par, pslice, -1, 0); // inverse, no translate

    SAG_FUNC sag_func;
    TRANSFER_DIST_FUNC transfer_dist_func;
    SURF_NORMAL_FUNC surf_normal_func;
    CRITICAL_XY_FUNC critical_xy_func;
    GetSurfaceFuncs(&sag_func, &transfer_dist_func, &surf_normal_func,
                    &critical_xy_func, pslice, p);

    double ln, mn, nn;
    surf_normal_func(&ln, &mn, &nn, x, y, pslice, 1);

    // Fill ray_out with intersection coords and surface normals
    ray_out->xs = coords_back[0];
    ray_out->ys = coords_back[1];
    ray_out->zs = coords_back[2];
    ray_out->t  = 0.0;
    ray_out->ln = ln;
    ray_out->mn = mn;
    ray_out->nn = nn;

    // Return the modified (l,m,n) separately
    if (l_out) *l_out = l;
    if (m_out) *m_out = m;
    if (n_out) *n_out = n;
}
