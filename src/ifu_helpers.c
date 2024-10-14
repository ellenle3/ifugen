#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "ifu_helpers.h"

/*
Helper functions for generating an image slicer IFU.
Ellen Lee
*/


double ConvertAngle(double t) {
    t = fmod(t, 360);
    if (t > 180) {
        return t - 360;
    }
    return t;
}

double Conic3DSag(double x, double y, int derv, double cv, double k, double alpha, double beta, double gamma) {
    // Also keep track of the angle, which determines which solution of
    // the quadratic to use
    alpha = ConvertAngle(alpha) * M_PI/180;
    beta = ConvertAngle(beta) * M_PI/180;
    gamma = ConvertAngle(gamma) * M_PI/180;
    int sgn;
    if (fabs(gamma) <= M_PI/2) {sgn = 1;}
    else {sgn = -1;}

    // Determine the off-axis distance
    if (fabs(cv) > 1E-10) {
        // If curvature is very small, it's basically a plane so no need to shift y-coordinate
        double y0 = sin(alpha) / ( cv * (1+cos(alpha)) );
        double x0 = sin(beta) / ( cv * (1+cos(beta)) );
        y = y + y0;
        x = x + x0;
    }

    // Rotate about the x-axis
    double cosg = cos(gamma);
    double cosg2 = cosg*cosg;
    double sing = sin(gamma);
    double sing2 = sing*sing;
    double asol = cv*(sing2 + (k+1)*cosg2);
    double bsol = -2*cosg*(x*k*cv*sing + 1);
    double csol = 2*x*sing + cv*(y*y + x*x*cosg2 + (k+1)*x*x*sing2);

    if (bsol*bsol - 4*asol*csol < 0 || fabs(2*asol) < 1E-10) {
        // Invalid solution
        return 0;
    }
    return (2*csol) / (-bsol + sgn*sqrt(bsol*bsol - 4*asol*csol));
}


void SliceAngles(double* alpha, double* beta, double* gamma, int slice_num, int col_num, IMAGE_SLICER_PARAMS p) {
    // Get row number as well as the subindex of the slice on that row
    int row_num = slice_num / p.n_each;
    int slice_num_row = slice_num - row_num * p.n_each;
    // Set the angles alpha and beta depending on the row and column
    *alpha = 0; *beta = 0; *gamma = 0;
    if (p.n_rows % 2 == 0) {
        *alpha = p.alpha_cen + p.dalpha*(row_num - (p.n_rows-1.)/2);
    }
    else {
        *alpha = p.alpha_cen + p.dalpha*(row_num - p.n_rows/2);
    }
    if (p.n_cols % 2 == 0) {
        *beta = p.beta_cen + p.dbeta*(col_num - (p.n_cols-1.)/2);
    }
    else {
        *beta = p.beta_cen + p.dbeta*(col_num - p.n_cols/2);
    }
    // Get the angles of the bottom- and top-most slices of the central row
    // If n_rows is even, there are 2 rows straddling the x=0 center line
    // Set the "central row" to the one above the x-axis (+y direction)
    double gamma_bot = 0;
    double gamma_top = 0;
    if (p.n_each % 2 == 0) {
        gamma_bot = p.gamma_cen - p.dgamma*(p.n_each-1.)/2;
        gamma_top = p.gamma_cen + p.dgamma*(p.n_each-1.)/2;
    }
    else {
        gamma_bot = p.gamma_cen - p.dgamma*(p.n_each/2);
        gamma_top = p.gamma_cen + p.dgamma*(p.n_each/2);
    }
    switch (p.mode) {
        case 0:
            // If the row is even, the angles are the same as the central row
            // If odd, the top/bottom angles are flipped
            if (row_num % 2 == 0) {
            *gamma = gamma_bot + slice_num_row*p.dgamma;
            }  
            else {
            *gamma = gamma_top - slice_num_row*p.dgamma;
            }
            break;
        case 1:
            // Copy the same angle pattern as the central row
            *gamma = gamma_bot + slice_num_row*p.dgamma;
            break;
    }
}

int ImageSlicerSag(double *z, double x, double y, IMAGE_SLICER_PARAMS p) {
    // Get dimensions of the image slicer
    int n_slices = p.n_each * p.n_rows;
    double ysize = (n_slices*p.dy + (n_slices-1)*p.gy_width);
    double xsize = (p.n_cols*p.dx + (p.n_cols-1)*p.gx_width);
    if (fabs(x) >= xsize/2 || fabs(y) >= ysize/2) {
        // x, y coordinate is out of bounds (x=0, y=0 is in the middle of the slicer)
        *z = 0;
        return -1;
    }
    // Figure out which slice and column number we are on to determine gaps
    int slice_num = floor((y + ysize/2) / (p.dy + p.gy_width));
    int col_num = floor((x + xsize/2) / (p.dx + p.gx_width));
    // If inside a gap, return the gap depth instead of a curved surface
    double ygap_bot = (slice_num + 1) * p.dy + slice_num*p.gy_width - ysize/2;
    double ygap_top = ygap_bot + p.gy_width;
    if (y > ygap_bot && y <= ygap_top) {
        *z = p.gy_depth;
        return 0;
    }
    double xgap_left = (col_num + 1) * p.dx + col_num*p.gx_width - xsize/2;
    double xgap_right = xgap_left + p.gx_width;
    if (x > xgap_left && x <= xgap_right) {
        *z = p.gx_depth;
        return 0;
    }
    // Inside a slice. From the slice number, determine the angles
    double* alphaptr;
    double* betaptr;
    double* gammaptr;
    SliceAngles(alphaptr, betaptr, gammaptr, slice_num, col_num, p);
    double alpha = *alphaptr;
    double beta = *betaptr;
    double gamma = *gammaptr;
    *z = Conic3DSag(x, y, p.cv, p.k, 0, alpha, beta, gamma);
    return 0;
}