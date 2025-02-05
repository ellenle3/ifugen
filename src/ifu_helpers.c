#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "ifu_helpers.h"

/*
Slicer generation and ray tracing algorithm.

Ellen Lee
*/


int CheckSlicerParams() {
    // Check that the parameters are valid.
    return 0;
}


double ConvertAngle(double t) {
    t = fmod(t, 360);
    if (t > 180) {
        return t - 360;
    }
    return t;
}


// Sag generation and helpers for determining the parameters of individual slices


void GetSliceAngles(double* alpha, double* beta, double* gamma, int slice_num, int col_num, IMAGE_SLICER_PARAMS p) {
    // Get row number as well as the subindex of the slice on that row
    int row_num = slice_num / p.n_each;
    int slice_num_row = slice_num - row_num * p.n_each;
    // Set the angles alpha and beta depending on the row and column
    *alpha = 0; *beta = 0; *gamma = 0;
    if (p.n_rows % 2 == 0) {
        *alpha = p.alpha_cen + p.dalpha * (row_num - (p.n_rows - 1.) / 2);
    }
    else {
        *alpha = p.alpha_cen + p.dalpha * (row_num - p.n_rows / 2);
    }
    if (p.n_cols % 2 == 0) {
        *beta = p.beta_cen + p.dbeta * (col_num - (p.n_cols - 1.)/2);
    }
    else {
        *beta = p.beta_cen + p.dbeta * (col_num - p.n_cols / 2);
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
        gamma_bot = p.gamma_cen - p.dgamma * (p.n_each / 2);
        gamma_top = p.gamma_cen + p.dgamma * (p.n_each / 2);
    }
    switch (p.mode) {
        case 0:
            // If the row is even, the angles are the same as the central row
            // If odd, the top/bottom angles are flipped
            if (row_num % 2 == 0) {
                *gamma = gamma_bot + slice_num_row * p.dgamma;
            }  
            else {
                *gamma = gamma_top - slice_num_row * p.dgamma;
            }
            break;
        case 1:
            // Copy the same angle pattern as the central row
            *gamma = gamma_bot + slice_num_row * p.dgamma;
            break;
    }
}


void GetSlicerSize(double *xsize, double *ysize, IMAGE_SLICER_PARAMS p) {
    int n_slices = p.n_each * p.n_rows;
    *ysize = n_slices * p.dy + (n_slices - 1) * p.gy_width;
    *xsize = p.n_cols * p.dx + (p.n_cols - 1) * p.gx_width;
}


void GetSlicerIndex(int *col_num, int *slice_num, double x, double y, IMAGE_SLICER_PARAMS p) {
    double xsize, ysize;
    GetSlicerSize(&xsize, &ysize, p);

    // Calculate the column and slice number based on x, y position
    *col_num = (x + xsize / 2) / (p.dx + p.gx_width);
    *slice_num = (y + ysize / 2) / (p.dy + p.gy_width);
}


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


void ImageSlicerSag(double *z, double x, double y, IMAGE_SLICER_PARAMS p) {
    // Get dimensions of the image slicer
    double xsize, ysize;
    GetSlicerSize(&xsize, &ysize, p);
    // Check if (x, y) is out of bounds
    if (fabs(x) >= xsize/2 || fabs(y) >= ysize/2) {
        *z = NAN;
        return;
    }

    // Figure out which slice and column number we are on to determine gaps
    int col_num, slice_num;
    GetSlicerIndex(&col_num, &slice_num, x, y, p);
    // If inside a gap, return the gap depth instead of a curved surface
    int in_xgap, in_ygap;
    IsInsideSlicerGap(&in_xgap, &in_ygap, x, y, p);
    if (in_xgap) {
        *z = p.gx_depth;
        return;
    }
    if (in_ygap) {
        *z = p.gy_depth;
        return;
    }

    // Inside a slice. From the slice number, determine the angles
    double alpha, beta, gamma;
    GetSliceAngles(&alpha, &beta, &gamma, slice_num, col_num, p);
    *z = Conic3DSag(x, y, p.cv, p.k, alpha, beta, gamma);
}


// Ray tracing


double FindBoundedSliceExtremum(double x0, double y0, int mode, IMAGE_SLICER_PARAMS p) {
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
    double alpha, beta, gamma;
    GetSliceAngles(&alpha, &beta, &gamma, slice_num, col_num, p);

    double xlo = col_num * (p.dx + p.gx_width) - xsize / 2;
    double xhi = xlo + p.dx;
    double ylo = slice_num * (p.dy + p.gy_width) - ysize / 2;
    double yhi = ylo + p.dy;

    // Compute critical points
    double xc1, xc2, yc;
    Conic3DCriticalXY(&xc1, &xc2, &yc, p.cv, p.k, alpha, beta, gamma);

    // There are up to 6 points to compare depending on whether the critical points
    // are within bounds. Get the edges first
    double zsolns[6];
    zsolns[0] = Conic3DSag(xlo, ylo, p.cv, p.k, alpha, beta, gamma);
    zsolns[1] = Conic3DSag(xlo, yhi, p.cv, p.k, alpha, beta, gamma);
    zsolns[2] = Conic3DSag(xhi, ylo, p.cv, p.k, alpha, beta, gamma);
    zsolns[3] = Conic3DSag(xhi, yhi, p.cv, p.k, alpha, beta, gamma);

    int n_compare = 4; // Number of elements to compare so far

    // Check whether critical points are in bounds. If yes, compute the sag an
    // add to the array of points to compare.
    if (!isnan(yc)) {
        if (yc >= ylo && yc <= yhi) {
            if (xc1 >= xlo && xc1 <= xhi) {
                zsolns[4] = Conic3DSag(xc1, yc, p.cv, p.k, alpha, beta, gamma);
                n_compare++;
            }
            if (!isnan(xc2)) {
                if (xc2 >= xlo && xc2 <= xhi) {
                    zsolns[5] = Conic3DSag(xc2, yc, p.cv, p.k, alpha, beta, gamma);
                    n_compare++;
                }
            }
        }
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


void FindSlicerGlobalExtrema(double *zmin, double *zmax, IMAGE_SLICER_PARAMS p) {
    // Determine which slice the global extrema are on by roughly sampling the
    // entire image slicer
    double xsize, ysize;
    GetSlicerSize(&xsize, &ysize, p);
    int nx = p.n_cols * 8;
    int ny = p.n_rows * p.n_each * 10;

    // Safeguard in case the user attempts to initialize an obscenely large number
    // of slices
    if (nx > 100000 || ny > 100000) {
        nx = 100000; 
        ny = 100000; // Implies >10,000 slices per column which seems excessive...
    }

    // Generally the number of points shouldn't be a problem but dynamically allocate
    // memory for xpts and ypts just in case
    double *xpts = (double *)malloc(nx * sizeof(double));
    double *ypts = (double *)malloc(ny * sizeof(double));
    if (xpts == NULL || ypts == NULL) {
        // If memory allocation fails, print an error and give up
        fprintf(stderr, "Memory allocation failed!\n");
        return;
    }

    for (int i = 0; i < nx; i++) {
        xpts[i] = -xsize / 2 + (xsize * i) / nx;
    }
    for (int i = 0; i < ny; i++) {
        ypts[i] = -ysize / 2 + (ysize * i) / ny;
    }

    // Evaluate the grid, keeping track of the maximum and minimum
    double x0_max = 0, y0_max = 0, z0_max = -INFINITY;
    double x0_min = 0, y0_min = 0, z0_min = INFINITY;
    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            double z = make_image_slicer(xpts[i], ypts[j], p);
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

    // Use the estimated maximum and minimum to find the exact values
    *zmin = FindBoundedSliceExtremum(x0_min, y0_min, 0, p);
    *zmax = FindBoundedSliceExtremum(x0_max, y0_max, 1, p);
}