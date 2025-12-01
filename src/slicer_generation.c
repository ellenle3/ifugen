#define _USE_MATH_DEFINES
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "slicer_generation.h"

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

void GetSurfaceFuncs(TRANSFER_DIST_FUNC *transfer_dist_func, SURF_NORMAL_FUNC *surf_normal_func,
    CRITICAL_XY_FUNC *critical_xy_func, TRANSFORMATION_FUNC *transform_func, SLICE_PARAMS pslice, IMAGE_SLICER_PARAMS_BASIC p) {
    if (pslice.cv == 0) {
        *transfer_dist_func = &PlaneTransfer;
        *surf_normal_func = &PlaneSurfaceNormal;
        *critical_xy_func = &PlaneCriticalXY;
        *transform_func = &PlaneTransformation;
        return;
    }
    if (p.surface_type == 1) {
        *transfer_dist_func = &CylinderTransfer;
        *surf_normal_func = &CylinderSurfaceNormal;
        *critical_xy_func = &CylinderCriticalXY;
        *transform_func = &CylinderTransformation;
        return;
    }
    *transfer_dist_func = &Conic2DTransfer;
    *surf_normal_func = &Conic2DSurfaceNormal;
    *critical_xy_func = &Conic2DCriticalXY;
    *transform_func = &Conic2DTransformation;
}

// Sag generation and helpers for determining the parameters of individual slices

// The size of the image slicer along x is determined by the number of columns and
// the width of the gaps between them. Ditto for y except with the total number
// of slices within a column.
void GetSlicerSize(double *xsize, double *ysize, IMAGE_SLICER_PARAMS_BASIC p) {
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
void GetSlicerIndex(int *col_num, int *slice_num, double x, double y, IMAGE_SLICER_PARAMS_BASIC p, double p_custom[]) {
    double xsize, ysize;
    GetSlicerSize(&xsize, &ysize, p);
    *slice_num = (int) floor((y + ysize / 2) / (p.dy + p.gy_width));
    int row_num = floor((double)*slice_num / p.n_each);
    // Column index depends on how much the row is shifted by (u)
    double u = GetUForRow(row_num, p_custom);
    *col_num = (int) floor((x - u + xsize / 2) / (p.dx + p.gx_width));
}


// Checks whether a point (x, y) is inside a gap. This is useful for sag and ray
// tracing computations.
void IsInsideSlicerGap(int *in_xgap, int *in_ygap, double x, double y, IMAGE_SLICER_PARAMS_BASIC p, double p_custom[]) {
    double xsize, ysize;
    int col_num, slice_num;
    GetSlicerSize(&xsize, &ysize, p);
    GetSlicerIndex(&col_num, &slice_num, x, y, p, p_custom);

    // Check vertical (y) gap
    double ygap_bot = (slice_num + 1) * p.dy + slice_num * p.gy_width - ysize / 2;
    double ygap_top = ygap_bot + p.gy_width;
    *in_ygap = (y > ygap_bot && y <= ygap_top);

    // Check horizontal (x) gap with u shift
    int row_num = floor((double)slice_num / p.n_each);
    double u = GetUForRow(row_num, p_custom);
    double xgap_left = (col_num + 1) * p.dx + col_num * p.gx_width - xsize / 2 + u;
    double xgap_right = xgap_left + p.gx_width;
    *in_xgap = (x > xgap_left && x <= xgap_right);
}

// If there are an even number of slices in a column, the "central" slice used
// for the paraxial ray trace is ambiguous. In this case we should allow the user
// to specify whether they want to slice above or below the center line. Same
// goes for the columns. The parameters p.active_x and p.active_y exist for this
// purpose.
void GetParaxialSliceIndex(int *col_num, int *slice_num, int active_x, int active_y, IMAGE_SLICER_PARAMS_BASIC p) {
    *col_num = p.n_cols / 2;
    if (p.n_cols % 2 == 0 && active_x) (*col_num)++;

    int n_slices = p.n_each * p.n_rows;
    *slice_num = n_slices / 2;
    if (n_slices % 2 == 0 && active_y) (*slice_num)++;
}

// Computes the maximum and minimum row offsets for the image slicer. umin and umax
// are needed along with the global extrema of the sag to compute the ray trace.
void GetMinMaxU(double *umin, double *umax, IMAGE_SLICER_PARAMS_BASIC p, double p_custom[]) {
    int i;
    *umin = p_custom[9];
    *umax = p_custom[9];
    for (i = 1; i < p.n_rows; ++i) {
        if (p_custom[9 + i] < *umin) *umin = p_custom[9 + i];
        if (p_custom[9 + i] > *umax) *umax = p_custom[9 + i];
    }
    return;
}

double ImageSlicerSag(double x, double y, IMAGE_SLICER_PARAMS_BASIC p, double p_custom[]) {
    // Check if out of bounds
    int col_num, slice_num;
    GetSlicerIndex(&col_num, &slice_num, x, y, p, p_custom);
    int row_num = floor((double)slice_num / p.n_each);
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
    SLICE_PARAMS pslice = GetSliceParams(slice_num, col_num, p_custom);
    TRANSFER_DIST_FUNC transfer_dist_func;
    SURF_NORMAL_FUNC surf_normal_func;
    CRITICAL_XY_FUNC critical_xy_func;
    TRANSFORMATION_FUNC transform_func;
    GetSurfaceFuncs(&transfer_dist_func, &surf_normal_func, &critical_xy_func, &transform_func, pslice, p);
    return SliceSag(x, y, pslice, transfer_dist_func, transform_func);
}


// Ray tracing

double FindBoundedSliceExtremum(double x0, double y0, int mode, IMAGE_SLICER_PARAMS_BASIC p, double p_custom[]) {
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
    SLICE_PARAMS pslice = GetSliceParams(slice_num, col_num, p_custom);
    TRANSFER_DIST_FUNC transfer_dist_func;
    SURF_NORMAL_FUNC surf_normal_func;
    CRITICAL_XY_FUNC critical_xy_func;
    TRANSFORMATION_FUNC transform_func;
    GetSurfaceFuncs(&transfer_dist_func, &surf_normal_func, &critical_xy_func, &transform_func, pslice, p);

    int row_num = floor((double)slice_num / p.n_each);
    double u = GetUForRow(row_num, p_custom);
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
    zsolns[0] = SliceSag(xlo, ylo, pslice, transfer_dist_func, transform_func);
    zsolns[1] = SliceSag(xlo, yhi, pslice, transfer_dist_func, transform_func);
    zsolns[2] = SliceSag(xhi, ylo, pslice, transfer_dist_func, transform_func);
    zsolns[3] = SliceSag(xhi, yhi, pslice, transfer_dist_func, transform_func);

    int n_compare = 4; // Number of elements to compare so far

    // Check whether critical points are in bounds. If yes, compute the sag an
    // add to the array of points to compare.
    if (yc >= ylo && yc <= yhi && xc >= xlo && xc <= xhi) {
        zsolns[4] = SliceSag(xc, yc, pslice, transfer_dist_func, transform_func);
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


void FindSlicerGlobalExtrema(double *zmin, double *zmax, IMAGE_SLICER_PARAMS_BASIC p, double p_custom[]) {
    // Determine which slice the global extrema are on by roughly sampling the
    // entire image slicer
    double xsize, ysize;
    GetSlicerSize(&xsize, &ysize, p);
    int nx = p.n_cols * 6;
    int ny = p.n_rows * p.n_each * 8;

    // Safeguard in case the user attempts to initialize an obscenely large number
    // of slices
    if (nx > 300000 || ny > 400000) {
        nx = 300000; // Implies 50,000 columns
        ny = 400000; // 50,000 slices per column...
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

double TransferFunction(double t, RAY_IN ray_in, IMAGE_SLICER_PARAMS_BASIC p, double p_custom[]) {
    double xs = ray_in.xt + t*ray_in.l;
    double ys = ray_in.yt + t*ray_in.m;
    double zs = t*ray_in.n;
    double sag = ImageSlicerSag(xs, ys, p, p_custom);
    return sag - zs;
}

RAY_BOUNDS GetRayBounds(RAY_IN ray_in, double umin, double umax, double zmin, double zmax,
    IMAGE_SLICER_PARAMS_BASIC p, void *p_custom) {

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
int IsRayInBounds(int nc_min, int ns_min, int nc_max, int ns_max, double umax, double umin, IMAGE_SLICER_PARAMS_BASIC p) {
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

int IsSectionInBounds(int col_num, int row_num, double umin, double umax, IMAGE_SLICER_PARAMS_BASIC p) {
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

int IsSectionValid(int col_num, int row_num, IMAGE_SLICER_PARAMS_BASIC p) {
    if (row_num < 0 || row_num >= p.n_rows) {
        return 0;
    }
    if (col_num < 0 || col_num >= p.n_cols) {
        return 0;
    }
    return 1;
}

void CheckSliceSolution(RAY_OUT *ray_out, double tol, RAY_IN ray_in, int ns_test, int nc_test,
    IMAGE_SLICER_PARAMS_BASIC p, void *p_custom) {

    *ray_out = (RAY_OUT){NAN, NAN, NAN, NAN, NAN, NAN, NAN};

    double xt = ray_in.xt; double yt = ray_in.yt;
    double l = ray_in.l; double m = ray_in.m; double n = ray_in.n;

    // Check if this slice is a solution. Get params for this slice and
    // compute the transfer distance.
    SLICE_PARAMS pslice = GetSliceParams(ns_test, nc_test, p_custom);
    TRANSFER_DIST_FUNC transfer_dist_func;
    SURF_NORMAL_FUNC surf_normal_func;
    CRITICAL_XY_FUNC critical_xy_func;
    TRANSFORMATION_FUNC transform_func;
    GetSurfaceFuncs(&transfer_dist_func, &surf_normal_func, &critical_xy_func, &transform_func, pslice, p);

    RAY_OUT ray_out_test = SliceRayTrace(ray_in, pslice, transfer_dist_func, surf_normal_func, transform_func, 1);
    double result = TransferFunction(ray_out_test.t, ray_in, p, p_custom);

    // Is this a valid zero to the transfer function?
    if (fabs(result) < tol) {
        // Yes - found a solution!

        // WAIT - Is the solution inside of a gap? If we're unlucky and zs is
        // equal to the gap depth then this may be the case.
        int in_xgap, in_ygap;
        IsInsideSlicerGap(&in_xgap, &in_ygap, ray_out_test.xs, ray_out_test.ys, p, p_custom);
        if (in_xgap || in_ygap) {
            return;
        }

        // Not in a gap. This slice is the solution.
        memcpy(ray_out, &ray_out_test, sizeof(RAY_OUT));
    }
}

void CheckYWallCollision(RAY_OUT *ray_out, RAY_IN ray_in, int ns_test, int nc_test, int sgns,
    IMAGE_SLICER_PARAMS_BASIC p, void *p_custom){

    *ray_out = (RAY_OUT){NAN, NAN, NAN, NAN, NAN, NAN, NAN};
    double xt = ray_in.xt; double yt = ray_in.yt;
    double l = ray_in.l; double m = ray_in.m; double n = ray_in.n;

    // Get y-dimension size of the slicer
    double xsize, ysize;
    GetSlicerSize(&xsize, &ysize, p);

    // Get sag function for the current column
    SLICE_PARAMS pslice = GetSliceParams(ns_test, nc_test, p_custom);
    TRANSFER_DIST_FUNC transfer_dist_func;
    SURF_NORMAL_FUNC surf_normal_func;
    CRITICAL_XY_FUNC critical_xy_func;
    TRANSFORMATION_FUNC transform_func;

    double t;

    if (fabs(l) <= 1e-13 && fabs(m) <= 1e-13 && p.gy_width > 0) {
        // Ray is going straight down z-axis. If you're checking for a wall collision
        // here then it must have gone into a gap.
        t = p.gy_depth / n;
        ray_out->xs = xt + t * l; ray_out->ys = yt + t * m; ray_out->zs = p.gy_depth;
        ray_out->t  = t; ray_out->ln = 0.0; ray_out->mn = 0.0; ray_out->nn = -1.0;
        return;

    }

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
        SLICE_PARAMS pslice_near = GetSliceParams(ns_test, nc_test, p_custom);
        GetSurfaceFuncs(&transfer_dist_func, &surf_normal_func, &critical_xy_func, &transform_func, pslice_near, p);
        double znear_slice = SliceSag(xnear, ynear, pslice_near, transfer_dist_func, transform_func);

        // Sag of far slice
        SLICE_PARAMS pslice_far = GetSliceParams(ns_test + sgns, nc_test, p_custom);
        GetSurfaceFuncs(&transfer_dist_func, &surf_normal_func, &critical_xy_func, &transform_func, pslice_far, p);
        double zfar_slice = SliceSag(xfar, yfar, pslice_far, transfer_dist_func, transform_func);

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
        if (zfar >= p.gy_depth && fabs(n) > 1e-13 && p.gy_width > 0) {
            t = p.gy_depth / n;
            ray_out->xs = xt + t * l; ray_out->ys = yt + t * m; ray_out->zs = p.gy_depth;
            ray_out->t  = t; ray_out->ln = 0.0; ray_out->mn = 0.0; ray_out->nn = -1.0;
            return;
        }
    }

    // Ray missed
    return;
}

void CheckXWallCollision(RAY_OUT *ray_out, RAY_IN ray_in, int ns_test, int nc_test, int sgnc,
    IMAGE_SLICER_PARAMS_BASIC p, void *p_custom) {

    *ray_out = (RAY_OUT){NAN, NAN, NAN, NAN, NAN, NAN, NAN};

    double xt = ray_in.xt; double yt = ray_in.yt;
    double l = ray_in.l; double m = ray_in.m; double n = ray_in.n;

    double xsize, ysize;
    GetSlicerSize(&xsize, &ysize, p);

    SLICE_PARAMS pslice = GetSliceParams(ns_test, nc_test, p_custom);
    TRANSFER_DIST_FUNC transfer_dist_func;
    SURF_NORMAL_FUNC surf_normal_func;
    CRITICAL_XY_FUNC critical_xy_func;
    TRANSFORMATION_FUNC transform_func;

    double t;
    
    if (fabs(l) <= 1e-13 && fabs(m) <= 1e-13 && p.gx_width > 0) {
        // Ray is going straight down z-axis. If you're checking for a wall collision
        // here then it must have gone into a gap.
        t = p.gx_depth / n;
        ray_out->xs = xt + t * l; ray_out->ys = yt + t * m; ray_out->zs = p.gx_depth;
        ray_out->t  = t; ray_out->ln = 0.0; ray_out->mn = 0.0; ray_out->nn = -1.0;
        return;

    }

    if (fabs(l) > 1e-13) {
        double u = GetUForRow(ns_test / p.n_each, p_custom);

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
        SLICE_PARAMS pslice_near = GetSliceParams(ns_test, nc_test, p_custom);
        GetSurfaceFuncs(&transfer_dist_func, &surf_normal_func, &critical_xy_func, &transform_func, pslice_near, p);
        double znear_slice = SliceSag(xnear, ynear, pslice_near, transfer_dist_func, transform_func);

        SLICE_PARAMS pslice_far = GetSliceParams(ns_test, nc_test + 1 * sgnc, p_custom);
        GetSurfaceFuncs(&transfer_dist_func, &surf_normal_func, &critical_xy_func, &transform_func, pslice_far, p);
        double zfar_slice = SliceSag(xfar, yfar, pslice_far, transfer_dist_func, transform_func);

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
        if (zfar >= p.gx_depth && fabs(n) > 1e-13 && p.gx_width > 0) {
            t = p.gx_depth / n;
            ray_out->xs = xt + t * l; ray_out->ys = yt + t * m; ray_out->zs = p.gx_depth;
            ray_out->t  = t; ray_out->ln = 0.0; ray_out->mn = 0.0; ray_out->nn = -1.0;
            return;
        }
    }

    return;
}

void CalcNextCoords(double *x_next, double *y_next, int *code, RAY_IN ray_in, int sgnc, int sgns, int nc_test, int nr_test,
    double x_test, double y_test, double xmax, double ymax, IMAGE_SLICER_PARAMS_BASIC p,
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

    double u = GetUForRow(nr_test, p_custom);

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
                          IMAGE_SLICER_PARAMS_BASIC p, void *p_custom,
                          int *n_stocheck, int *nc_new, int *nr_new) {
    int ns1 = 0, ns2 = 0;

    // Get slice indices (assuming GetSlicerIndex returns column and slice indices)
    int dummy_nc;

    GetSlicerIndex(&dummy_nc, &ns1, x_test, y_test, p, p_custom);
    GetSlicerIndex(&dummy_nc, &ns2,
                   x_next - 1e-15 * sgnc,
                   y_next - 1e-15 * sgns,
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

int IsLastSliceInSection(int ns_test, int sgns, IMAGE_SLICER_PARAMS_BASIC p) {
    if (sgns > 0) {
        return (ns_test % p.n_each == p.n_each - 1);
    } else if (sgns < 0) {
        return (ns_test % p.n_each == 0);
    }
    return 0;
}


void RayTraceSlicer(RAY_OUT *ray_out, RAY_IN ray_in, double zmin, double zmax, double umin, double umax,
                    int trace_walls, IMAGE_SLICER_PARAMS_BASIC p, double p_custom[]) {
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
        
        // Only check slices if the section is valid
        if (IsSectionValid(nc_test, nr_test, p)) {

            // Iterate slices within this section
            for (int n = 0; n < n_stocheck; n++) {
                // Check whether each slice is a solution to the transfer equation
                CheckSliceSolution(ray_out, tol, ray_in, ns_test, nc_test, p, p_custom);

                if (!isnan(ray_out->t)) {
                    return;  // Found a solution within a slice
                }

                // Before going to the next slice, check if there is a collision with
                // a wall. Skip if this is the last slice in the section.
                if (trace_walls && !IsLastSliceInSection(ns_test, sgns, p)) {
                    CheckYWallCollision(ray_out, ray_in, ns_test, nc_test, sgns, p, p_custom);
                    if (!isnan(ray_out->t)) {
                        return;  // Found a wall collision
                    }
                }

                // Increment the slice index
                ns_test += sgns;
            }

            // If crossing columns, need to double check the last slice index in the
            // next section
            if (nc_new != nc_test) {
                ns_test -= sgns;
            }
        }

        // Before continuing, check walls between sections.
        is_same_section = (nc_new == nc_test && nr_new == nr_test);

        // If both the current and next sections are invalid, skip wall checks
        if (IsSectionValid(nc_test, nr_test, p) || IsSectionValid(nc_new, nr_new, p)) {

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
                                    xmax - 1e-15 * sgnc, y_test - 1e-15 * sgns,
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
    IMAGE_SLICER_PARAMS_BASIC p, double p_custom[]) 
{
    int nc, ns;
    GetParaxialSliceIndex(&nc, &ns, active_x, active_y, p);
    SLICE_PARAMS pslice = GetSliceParams(ns, nc, p_custom);

    TRANSFER_DIST_FUNC transfer_dist_func;
    SURF_NORMAL_FUNC surf_normal_func;
    CRITICAL_XY_FUNC critical_xy_func;
    TRANSFORMATION_FUNC transform_func;
    GetSurfaceFuncs(&transfer_dist_func, &surf_normal_func, &critical_xy_func, &transform_func, pslice, p);

    // Transform the ray into local coordinates
    double coords[3]   = { ray_in->xt, ray_in->yt, 0.0 };
    double cosines[3]  = { ray_in->l,  ray_in->m,  ray_in->n };
    double coords_local[3], cosines_local[3];

    transform_func(coords_local, coords, pslice, -1, 1);   // forward, translate
    transform_func(cosines_local, cosines, pslice, -1, 0); // forward, no translate

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

    double x = coords_local[0];
    double y = coords_local[1];
    double z = coords_local[2];
    double l = cosines_local[0];
    double m = cosines_local[1];
    double n = cosines_local[2];

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

    cosines_local[0] = l;
    cosines_local[1] = m;
    cosines_local[2] = n;

    double surf_normal_local[3] = {0.0, 0.0, -1.0}; // Paraxial approx: normal is along -z

    // Transform back to global coordinates
    double coords_back[3], cosines_back[3], surf_normal_back[3];
    transform_func(coords_back, coords_local, pslice, 1, 1);  // inverse, translate
    transform_func(cosines_back, cosines_local, pslice, 1, 0); // inverse, no translate
    transform_func(surf_normal_back, surf_normal_local, pslice, 1, 0); // inverse, no translate

    // Fill ray_out with intersection coords and surface normals
    ray_out->xs = coords_back[0];
    ray_out->ys = coords_back[1];
    ray_out->zs = coords_back[2];
    ray_out->t  = 0.0;
    ray_out->ln = surf_normal_back[0];
    ray_out->mn = surf_normal_back[1];
    ray_out->nn = surf_normal_back[2];

    *l_out = cosines_back[0];
    *m_out = cosines_back[1];
    *n_out = cosines_back[2];
}
