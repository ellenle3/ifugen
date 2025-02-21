#define _USE_MATH_DEFINES
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "slicer_generation.h"
#include "surface_solns.h"
#include "custom_slicer_helpers.h"

/*
Slicer generation and ray tracing algorithm.

Ellen Lee
*/

// TODO: This function does not check that array is large enough to store n
// values! This could lead to unexpected behavior...
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

    if (!(p->cylinder==0 || p->cylinder==1)) { p->cylinder=0; is_valid = 0; }
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
        p1.cylinder == p2.cylinder &&
        p1.n_each == p2.n_each &&
        p1.n_rows == p2.n_rows &&
        p1.n_cols == p2.n_cols &&
        p1.angle_mode == p2.angle_mode &&
        p1.dalpha == p2.dalpha &&
        p1.dbeta == p2.dbeta &&
        p1.dgamma == p2.dgamma &&
        p1.gamma_offset == p2.gamma_offset &&
        p1.alpha_cen == p2.alpha_cen &&
        p1.beta_cen == p2.beta_cen &&
        p1.gamma_cen == p2.gamma_cen &&
        p1.dx == p2.dx &&
        p1.dy == p2.dy &&
        p1.gx_width == p2.gx_width &&
        p1.gx_depth == p2.gx_depth &&
        p1.gy_width == p2.gy_width &&
        p1.gy_depth == p2.gy_depth &&
        p1.cv == p2.cv &&
        p1.k == p2.k
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
IMAGE_SLICER_PARAMS MakeSlicerParamsFromCustom(double custom_slice_params[]) {
    IMAGE_SLICER_PARAMS p;
    p.custom = 1;
    p.n_each = 1;
    p.n_rows = (int) custom_slice_params[0];
    p.n_cols = (int) custom_slice_params[1];
    p.cylinder = (int) custom_slice_params[2];
    p.dx = custom_slice_params[3];
    p.dy = custom_slice_params[4];
    p.gx_width = custom_slice_params[5];
    p.gx_depth = custom_slice_params[6];
    p.gy_width = custom_slice_params[7];
    p.gy_depth = custom_slice_params[8];

    // The rest of these will not be touched, set to zero
    p.cv = 0;
    p.k = 0;
    p.angle_mode = 0;
    p.dalpha = 0;
    p.dbeta = 0;
    p.dgamma = 0;
    p.gamma_offset = 0;
    p.alpha_cen = 0;
    p.beta_cen = 0;
    p.gamma_cen = 0;
    return p;
}

void GetSurfaceFuncs(SAG_FUNC *sag_func, TRANSFER_DIST_FUNC *transfer_dist_func,
SURF_NORMAL_FUNC *surf_normal_func, CRITICAL_XY_FUNC *critical_xy_func, double cv, int is_cylinder) {
    // If the curvature is 0, then the image slicer is a plane
    if (cv == 0) {
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
void GetSlicerIndex(int *col_num, int *slice_num, double x, double y, IMAGE_SLICER_PARAMS p) {
    double xsize, ysize;
    GetSlicerSize(&xsize, &ysize, p);
    // Calculate the column and slice number based on x, y position
    *col_num = floor( (x + xsize / 2) / (p.dx + p.gx_width) );
    *slice_num = floor( (y + ysize / 2) / (p.dy + p.gy_width) );
}

// Checks whether a point (x, y) is inside a gap. This is useful for sag and ray
// tracing computations.
void IsInsideSlicerGap(int *in_xgap, int *in_ygap, double x, double y, IMAGE_SLICER_PARAMS p) {
    double xsize, ysize;
    int col_num, slice_num;
    GetSlicerSize(&xsize, &ysize, p);
    GetSlicerIndex(&col_num, &slice_num, x, y, p);
    
    // Check vertical (y) gap
    double ygap_bot = (slice_num + 1) * p.dy + slice_num * p.gy_width - ysize/2;
    double ygap_top = ygap_bot + p.gy_width;
    *in_ygap = (y > ygap_bot && y <= ygap_top);
    // Check horizontal (x) gap
    double xgap_left = (col_num + 1) * p.dx + col_num * p.gx_width - xsize/2;
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


// If p.custom, then use the slice parameters provided in the array custom_slice_params.
// Otherwise, this array is ignored.
void GetSliceParams(double* alpha, double* beta, double* gamma, double* cv, double* k,
int slice_num, int col_num, IMAGE_SLICER_PARAMS p, double custom_slice_params[]) {
    if (p.custom) {
        GetSliceParamsCustom(alpha, beta, gamma, cv, k, slice_num, col_num, custom_slice_params);
        return;
    }
    GetSliceParamsStandard(alpha, beta, gamma, cv, k, slice_num, col_num, p);
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
void GetSliceParamsStandard(double* alpha, double* beta, double* gamma, double* cv, double* k, int slice_num, int col_num, IMAGE_SLICER_PARAMS p) {
    // cv and k are determined by p and are the same for all slices
    *cv = p.cv; *k = p.k;
    
    // Get row number as well as the subindex of the slice on that row
    int row_num = slice_num / p.n_each;
    int slice_num_row = slice_num - row_num * p.n_each;
    double gamma_extra;

    // Set the angles alpha and beta depending on the row and column
    *alpha = 0; *beta = 0; *gamma = 0;
    if (p.n_rows % 2 == 0) {
        *alpha = p.alpha_cen + p.dalpha * (row_num - (p.n_rows - 1.) / 2);
        gamma_extra = p.gamma_offset * (row_num - (p.n_rows - 1.) / 2);
    }
    else {
        *alpha = p.alpha_cen + p.dalpha * (row_num - floor(p.n_rows / 2) );
        gamma_extra = p.gamma_offset * (row_num - floor(p.n_rows / 2) );
    }
    if (p.n_cols % 2 == 0) {
        *beta = p.beta_cen + p.dbeta * (col_num - (p.n_cols - 1.)/2);
    }
    else {
        *beta = p.beta_cen + p.dbeta * (col_num - floor(p.n_cols / 2) );
    }

    // Get the angles of the bottom- and top-most slices of the central row
    // If n_rows is even, there are 2 rows straddling the x=0 center line
    // Set the "central row" to the one above the x-axis (+y direction)
    double gamma_bot = 0;
    double gamma_top = 0;
    if (p.n_each % 2 == 0) {
        gamma_bot = p.gamma_cen - p.dgamma * (p.n_each - 1.) / 2;
        gamma_top = p.gamma_cen + p.dgamma * (p.n_each - 1.) / 2;
    }
    else {
        gamma_bot = p.gamma_cen - p.dgamma * floor(p.n_each / 2);
        gamma_top = p.gamma_cen + p.dgamma * floor(p.n_each / 2);
    }

    // Determine offsets in gamma. First, check whether the mode allows the extra
    // offsets to stack or not. If no, then set gamma_extra to repeat every 2 rows.
    if (p.angle_mode == 2 || p.angle_mode == 3) {
        // Offsets should not stack
        if (row_num % 2 == 0) gamma_extra = -p.gamma_offset / 2;
        else gamma_extra = p.gamma_offset / 2;
    }

    switch (p.angle_mode) {
        case 0:
        case 2:
            // If the row is even, the angles are the same as the central row
            // If odd, the top/bottom angles are flipped
            if (row_num % 2 == 0) {
                *gamma = gamma_bot + slice_num_row * p.dgamma + gamma_extra;
            }  
            else {
                *gamma = gamma_top - slice_num_row * p.dgamma + gamma_extra;
            }
            break;

        case 1:
        case 3:
            // Copy the same angle pattern as the central row
            *gamma = gamma_bot + slice_num_row * p.dgamma + gamma_extra;
            break;
    }
}

// The sag of the image slicer can be found by
double ImageSlicerSag(double x, double y, IMAGE_SLICER_PARAMS p, double custom_slice_params[]) {
    // Get dimensions of the image slicer
    double xsize, ysize;
    GetSlicerSize(&xsize, &ysize, p);
    // Check if (x, y) is out of bounds
    if (fabs(x) >= xsize/2 || fabs(y) >= ysize/2) {
        return NAN;
    }

    // Figure out which slice and column number we are on to determine gaps
    int col_num, slice_num;
    GetSlicerIndex(&col_num, &slice_num, x, y, p);
    // If inside a gap, return the gap depth instead of a curved surface
    int in_xgap, in_ygap;
    IsInsideSlicerGap(&in_xgap, &in_ygap, x, y, p);
    if (in_xgap) {
        return p.gx_depth;
    }
    if (in_ygap) {
        return p.gy_depth;
    }

    // Inside a slice. From the slice number, determine the angles
    double alpha, beta, gamma, cv, k;
    GetSliceParams(&alpha, &beta, &gamma, &cv, &k, slice_num, col_num, p, custom_slice_params);

    SAG_FUNC sag_func;
    TRANSFER_DIST_FUNC transfer_dist_func;
    SURF_NORMAL_FUNC surf_normal_func;
    CRITICAL_XY_FUNC critical_xy_func;
    GetSurfaceFuncs(&sag_func, &transfer_dist_func, &surf_normal_func, &critical_xy_func, cv, p.cylinder);
    return sag_func(x, y, cv, k, alpha, beta, gamma);
}


// Ray tracing

double FindBoundedSliceExtremum(double x0, double y0, int mode, IMAGE_SLICER_PARAMS p, double custom_slice_params[]) {
    // If x0, y0 are inside gaps then we don't need to go further
    int in_xgap, in_ygap;
    IsInsideSlicerGap(&in_xgap, &in_ygap, x0, y0, p);
    if (in_xgap) return p.gx_depth;
    if (in_ygap) return p.gy_depth;

    // Not inside a gap. We need to compute the sag at the bounds and critical
    // point(s) of the slice.
    // First, compute the bounds of the current slice
    int col_num, slice_num;
    GetSlicerIndex(&col_num, &slice_num, x0, y0, p);
    double xsize, ysize;
    GetSlicerSize(&xsize, &ysize, p);
    double alpha, beta, gamma, cv, k;
    GetSliceParams(&alpha, &beta, &gamma, &cv, &k, slice_num, col_num, p, custom_slice_params);
    SAG_FUNC sag_func;
    TRANSFER_DIST_FUNC transfer_dist_func;
    SURF_NORMAL_FUNC surf_normal_func;
    CRITICAL_XY_FUNC critical_xy_func;
    GetSurfaceFuncs(&sag_func, &transfer_dist_func, &surf_normal_func, &critical_xy_func, cv, p.cylinder);

    double xlo = col_num * (p.dx + p.gx_width) - xsize / 2;
    double xhi = xlo + p.dx;
    double ylo = slice_num * (p.dy + p.gy_width) - ysize / 2;
    double yhi = ylo + p.dy;

    // Compute critical point
    double xc, yc;
    critical_xy_func(&xc, &yc, cv, k, alpha, beta, gamma);

    // There are up to 5 points to compare depending on whether the critical point
    // is within bounds.
    double zsolns[5];
    zsolns[0] = sag_func(xlo, ylo, cv, k, alpha, beta, gamma);
    zsolns[1] = sag_func(xlo, yhi, cv, k, alpha, beta, gamma);
    zsolns[2] = sag_func(xhi, ylo, cv, k, alpha, beta, gamma);
    zsolns[3] = sag_func(xhi, yhi, cv, k, alpha, beta, gamma);

    int n_compare = 4; // Number of elements to compare so far

    // Check whether critical points are in bounds. If yes, compute the sag an
    // add to the array of points to compare.
    if (yc >= ylo && yc <= yhi && xc >= xlo && xc <= xhi) {
        zsolns[4] = sag_func(xc, yc, cv, k, alpha, beta, gamma);
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


void FindSlicerGlobalExtrema(double *zmin, double *zmax, IMAGE_SLICER_PARAMS p, double custom_slice_params[]) {
    // Determine which slice the global extrema are on by roughly sampling the
    // entire image slicer
    double xsize, ysize;
    GetSlicerSize(&xsize, &ysize, p);
    int nx = p.n_cols * 6;
    int ny = p.n_rows * p.n_each * 8;

    // Safeguard in case the user attempts to initialize an obscenely large number
    // of slices
    if (nx > 30000 || ny > 40000) {
        nx = 30000; // 5,000 columns???
        ny = 40000; // Implies 5,000 slices per column which seems excessive...
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
            z = ImageSlicerSag(xpts[i], ypts[j], p, custom_slice_params);
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
    *zmin = FindBoundedSliceExtremum(x0_min, y0_min, 0, p, custom_slice_params);
    *zmax = FindBoundedSliceExtremum(x0_max, y0_max, 1, p, custom_slice_params);
}

double TransferEquation(double t, double xt, double yt, double l, double m, double n,
IMAGE_SLICER_PARAMS p, double custom_slice_params[]) {
    double xs = xt + t*l;
    double ys = yt + t*m;
    double zs = t*n;
    double sag = ImageSlicerSag(xs, ys, p, custom_slice_params);
    return sag - zs;
}

void GetRayBounds(int *nc_min, int *ns_min, int *nc_max, int *ns_max, RAY_IN ray_in, double zmin, double zmax, IMAGE_SLICER_PARAMS p) {
    double xmin, xmax, ymin, ymax;
    double xsize, ysize;
    GetSlicerSize(&xsize, &ysize, p);
    
    if (fabs(ray_in.n) < 1e-13) {
        // Ray is moving perpendicular to the z-axis! This is an unusual case.
        // Set bounds on potential xs and ys to edges of the image slicer
        xmin = -xsize / 2; xmax = xsize / 2;
        ymin = -ysize / 2; ymax = ysize / 2;
        
        if (ray_in.l < 0) {xmin *= -1; xmax *= -1;}  // propagate right to left
        if (ray_in.m < 0) {ymin *= -1; ymax *= -1;}  // propagate bottom to top
    }
        
    else {
        // Get maximum and minimum possible values of xs and ys
        double tmin = zmin / ray_in.n;
        xmin = ray_in.xt + tmin*ray_in.l; ymin = ray_in.yt + tmin*ray_in.m;
        double tmax = zmax / ray_in.n;
        xmax = ray_in.xt + tmax*ray_in.l; ymax = ray_in.yt + tmax*ray_in.m;
    }
        

    // Get starting and ending rows and columns
    GetSlicerIndex(nc_min, ns_min, xmin, ymin, p);
    GetSlicerIndex(nc_max, ns_max, xmax, ymax, p);
}

int IsRayInBounds(int nc_min, int ns_min, int nc_max, int ns_max, IMAGE_SLICER_PARAMS p) {
    if ( (nc_min < 0 && nc_max < 0) || (nc_min >= p.n_cols && nc_max >= p.n_cols)) {
        return 0;  // x-value of the ray is too high or low
    }

    int n_sperc = p.n_each * p.n_rows;  // number of slices per column
    if ( (ns_min < 0 && ns_max < 0) || (ns_min >= n_sperc && ns_max >= n_sperc) ) {
        return 0;  // y-value of the ray is too high or low
    }
    
    return 1;
}

void RayTraceSlicer(RAY_OUT *ray_out, RAY_IN ray_in, double zmin, double zmax, int trace_walls, IMAGE_SLICER_PARAMS p, double custom_slice_params[]) {

    // Tolerance for accepting the transfer distance as valid
    double tol = 1e-11;
    ray_out->xs = ray_out->ys = ray_out->zs = NAN;
    ray_out->t = ray_out->ln = ray_out->mn = ray_out->nn = NAN;

    int nc_min, ns_min, nc_max, ns_max;
    GetRayBounds(&nc_min, &ns_min, &nc_max, &ns_max, ray_in, zmin, zmax, p);

    // Ray missed
    if (!IsRayInBounds(nc_min, ns_min, nc_max, ns_max, p)) return;

    // Ray is in bounds at least some of the time. Start from the min col and slice
    // indices and check solutions until we hit the max. Initialize the loop parameters...
    int nc_test = nc_min; int ns_test = ns_min;
    double x_test = ray_in.xt; double y_test = ray_in.yt;

    int dcol = abs(nc_max - nc_test);        // Number of col indices left to iterate
    int dslice = abs(ns_max - ns_test);      // Number of slice (row) indices left to iterate

    // Keep track of signs - which way to iterate. Set to zero for now.
    // If dcol is 0 then sgnc doesn't matter since we aren't iterating it anyway.
    // Same goes for dslice and sgns.
    int sgnc = (dcol > 0) ? (nc_max - nc_min) / dcol : 0;
    int sgns = (dslice > 0) ? (ns_max - ns_min) / dslice : 0;
        
    int slice_iter = 0;       // Keep track of whether we're on the first slice in a column to check walls
    
    double xsize, ysize;
    GetSlicerSize(&xsize, &ysize, p);
    double xt = ray_in.xt, yt = ray_in.yt;
    double l = ray_in.l, m = ray_in.m, n = ray_in.n;

    int n_sforc, in_xgap, in_ygap;
    double x_cross, y_cross, alpha, beta, gamma, cv, k, t_test, result;
    double xs, ys, zs, t, ln, mn, nn;
    double xnear, ynear, znear, tnear, xfar, yfar, zfar, tfar, znear_slice, zfar_slice, zcompare;
    SAG_FUNC sag_func;
    TRANSFER_DIST_FUNC transfer_dist_func;
    SURF_NORMAL_FUNC surf_normal_func;
    CRITICAL_XY_FUNC critical_xy_func;

    // Begin iterating through slices

    while (dcol >= 0) {
        
        // How many slices do we need to check on this column?

        // No column switching needed, check all remaining slices
        if (fabs(l) < 1e-13 || dcol == 0) {n_sforc = dslice;}

        // Ray can potentially switch columns
        else {
            // Crossover point between current column and next one
            x_cross = nc_test * (p.dx + p.gx_width) + (1 + sgnc) * p.dx / 2;
            // y-intercept between ray and the next column to check
            y_cross = yt + (x_cross - x_test) * m / l;
            // Number of slices to check
            n_sforc = ceil( fabs(y_cross - y_test) / (p.dy + p.gy_width) );
            x_test = x_cross; y_test = y_cross;
        }
            
        while (dslice >= 0 && n_sforc >= 0) {
            
            if (nc_test >= 0 && ns_test >= 0) {

                // Check if out of bounds
                GetSliceParams(&alpha, &beta, &gamma, &cv, &k, ns_test, nc_test, p, custom_slice_params);
                GetSurfaceFuncs(&sag_func, &transfer_dist_func, &surf_normal_func, &critical_xy_func, cv, p.cylinder);

                t_test = transfer_dist_func(xt, yt, l, m, n, cv, k, alpha, beta, gamma);
                result = TransferEquation(t_test, xt, yt, l, m, n, p, custom_slice_params);
                
                // Check whether the transfer distance of the current slice is a
                // valid zero of the transfer equation.
                if (fabs(result) < tol) {
                    // Yes - found a solution!
                    t = t_test;
                    xs = xt + t_test*l;
                    ys = yt + t_test*m;
                    zs = t_test*n;
                    surf_normal_func(&ln, &mn, &nn, xs, ys, cv, k, alpha, beta, gamma, 1);

                    // WAIT - Is the solution inside of a gap? If we're unlucky
                    // and zs is equal to the gap depth then this may be the case.
                    IsInsideSlicerGap(&in_xgap, &in_ygap, xs, ys, p);
                    if (in_xgap || in_ygap) {
                        dslice = -1; dcol = -1;
                        break;
                    }

                    ray_out->xs = xs; ray_out->ys = ys; ray_out->zs = zs; ray_out->t = t;
                    ray_out->ln = ln; ray_out->mn = mn; ray_out->nn = nn;
                    return;
                }

                // Check if the ray is hitting a wall after this slice
                if (trace_walls) {

                    // Going between columns since slice_iter==0. Skip if there
                    // are no columns left to iterate
                    if (slice_iter == 0 && fabs(l) > 1E-13 && dcol > 0) {
                        
                        xnear = (nc_test + 1) * p.dx - xsize / 2;
                        tnear = (xnear - xt) / l;
                        ynear = yt + tnear*m;
                        znear = tnear*n;
                        
                        xfar = (nc_test + 1) * p.dx + sgnc * p.gx_width - xsize / 2;
                        tfar = (xfar - xt) / l;
                        yfar = yt + tfar*m;
                        zfar = tfar*n;

                        znear_slice = sag_func(xnear, ynear, cv, k, alpha, beta, gamma);
                        GetSliceParams(&alpha, &beta, &gamma, &cv, &k, ns_test, nc_test + 1*sgnc, p, custom_slice_params);
                        zfar_slice = sag_func(xfar, yfar, cv, k, alpha, beta, gamma);

                        // Did it hit a near wall? This is only possible if the gap depth
                        // protrudes further in -z than the near edge
                        // THIS IS ALSO POSSIBLE IF THE RAY APPROACHES FROM THE WRONG SIDE
                        // OF THE IMAGE SLICER! But this should never happen under normal circumstances.
                        if ( (znear_slice > p.gx_depth && p.gx_width > 0) && (znear < znear_slice && znear > p.gx_depth) ) {
                            ray_out->xs = xnear; ray_out->ys = ynear; ray_out->zs = znear; ray_out->t = tnear;
                            ray_out->ln = -1*sgnc; ray_out->mn = 0; ray_out->nn = 0;
                            return;
                        }

                        // Did it hit a far wall?
                        zcompare = (p.gx_width == 0) ? znear_slice : p.gx_depth;
                        if ( (zfar > zfar_slice && zfar < zcompare) || (zfar < zfar_slice && zfar > zcompare) ) {
                            ray_out->xs = xfar; ray_out->ys = yfar; ray_out->zs = zfar; ray_out->t = tfar;
                            ray_out->ln = -1*sgnc; ray_out->mn = 0; ray_out->nn = 0;
                            return;
                        }
                    }
                
                    // Not going between columns, check walls along y-axis instead.
                    // Same as above...
                    else {
                        if (fabs(m) > 1E-13 && dslice > 0) {
                            ynear = (ns_test + 1) * p.dy - ysize / 2;
                            tnear = (ynear - yt) / m;
                            xnear = xt + tnear*l;
                            znear = tnear*n;
                            
                            yfar = (ns_test + 1) * p.dy + sgns * p.gy_width - ysize / 2;
                            tfar = (yfar - yt) / m;
                            xfar = xt + tfar*l;
                            zfar = tfar*n;
                            
                            znear_slice = sag_func(xnear, ynear, cv, k, alpha, beta, gamma);
                            GetSliceParams(&alpha, &beta, &gamma, &cv, &k, ns_test + 1*sgns, nc_test, p, custom_slice_params);
                            zfar_slice = sag_func(xfar, yfar, cv, k, alpha, beta, gamma);

                            if ( (znear_slice > p.gy_depth && p.gy_width > 0) && (znear < znear_slice && znear > p.gy_depth) ) {
                                ray_out->xs = xnear; ray_out->ys = ynear; ray_out->zs = znear; ray_out->t = tnear;
                                ray_out->ln = 0; ray_out->mn = -1*sgns; ray_out->nn = 0;
                                return;
                            }

                            // Did it hit a far wall?
                            zcompare = (p.gy_width == 0) ? znear_slice : p.gy_depth;
                            if ( (zfar > zfar_slice && zfar < zcompare) || (zfar < zfar_slice && zfar > zcompare) ) {
                                ray_out->xs = xfar; ray_out->ys = yfar; ray_out->zs = zfar; ray_out->t = tfar;
                                ray_out->ln = 0; ray_out->mn = -1*sgns; ray_out->nn = 0;
                                return;
                            }
                        }
                    }
                }
            }
            
            // Not a solution - increment slice and try again.
            slice_iter++;
            ns_test += sgns;
            dslice -= 1;
        }

        // The last slice index needs to be checked again when crossing over to
        // the next column
        slice_iter = 0;
        ns_test -= sgns;
        dslice++;
        // Go to the next column
        nc_test += sgnc;
        dcol -= 1;
    }

    // If none of the above worked then we probably hit a gap
    if (!trace_walls || fabs(n) < 1e-13) return;

    // Gaps between columns take precedence over gaps between slices in rows
    t_test = p.gx_depth / n;
    if ( fabs(TransferEquation(t_test, xt, yt, l, m, n, p, custom_slice_params)) < tol ) {
        t = t_test;
        ray_out->xs = xt + t*l; ray_out->ys = yt + t*m; ray_out->zs = p.gx_depth; ray_out->t = t;
        ray_out->ln = 0; ray_out->mn = 0; ray_out->nn = -1;
        return;
    }
    t_test = p.gy_depth / n;
    if ( fabs(TransferEquation(t_test, xt, yt, l, m, n, p, custom_slice_params)) < tol ) {
        t = t_test;
        ray_out->xs = xt + t*l; ray_out->ys = yt + t*m; ray_out->zs = p.gy_depth; ray_out->t = t;
        ray_out->ln = 0; ray_out->mn = 0; ray_out->nn = -1;
        return;
    }

    // If none of that worked then this is a bizarre edge case, e.g., the direction
    // cosine is less than 1E-13 but it grazed off of a wall somehow. Sweep this
    // ray under the rug and say that it missed...
    
    return;  // *5 tumbleweeds roll along in front of you*
}