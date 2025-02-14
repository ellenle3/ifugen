#define _USE_MATH_DEFINES
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "surface_solns.h"

/*
Solutions for the sag, critical points, ray transfer distance, and surface normals
for different surface types.

Ellen Lee
*/

// Bounds the angle to be between -180 and 180 degrees.
double ConvertAngle(double t) {
    t = fmod(t, 360);
    if (t > 180) {
        return t - 360;
    }
    return t;
}

// 3D conic solutions - planes (cv -> 0) are handled differently, see next section

// Determine the off-axis distance. If the curvature is near-zero, this is
// basically a plane. Set a limit to how small cv can be. In practice the
// user should not set a non-zero off-axis distance if the surface is a plane.
void Conic3DOffAxisDistance(double *x0, double *y0, double cv, double alpha, double beta) {
    if (fabs(cv) < 1E-13) cv = 1E-13;
    *y0 = sin(alpha) / ( cv * (1 + cos(alpha)) );  // Angles must be in radians
    *x0 = sin(beta) / ( cv * (1 + cos(beta)) );
}

// This is the sag for a conic surface translated by x0 and y0, and then rotated
// about the global y-axis. The sag for this surface is solved by doing the
// transformations backward, i.e., undoing the rotation and then the translation.
// The reason for this is that the sag is easier to solve for this way, but it
// is equivalent. Refer to the derivation for more informaiton.
double Conic3DSag(double x, double y, double cv, double k, double alpha, double beta, double gamma) {
    alpha = ConvertAngle(alpha) * M_PI/180;
    beta = ConvertAngle(beta) * M_PI/180;
    gamma = ConvertAngle(gamma) * M_PI/180;

    // The value of gamma determines which solution of the quadratic to use
    int sgn;
    if (fabs(gamma) <= M_PI/2) sgn = 1;
    else sgn = -1;

    double x0, y0;
    Conic3DOffAxisDistance(&x0, &y0, cv, alpha, beta);

    if (k == -1 && gamma == 0) {
        return cv * (x*x + y*y) / 2;
    }

    // Rotate about the x-axis
    double cosg = cos(gamma), cosg2 = cosg*cosg;
    double sing = sin(gamma), sing2 = sing*sing;

    double asol = cv*(sing2 + (k+1)*cosg2);
    double bsol = 2*cv*sing*(x*k*cosg - x0) - 2*cosg;
    double csol = cv*k*x*x*sing2 - 2*x*sing + 2*cv*x0*x*cosg + cv*(x*x + x0*x0 + (y-y0)*(y-y0));

    // Invalid solution, sag is undefined here
    if (bsol*bsol - 4*asol*csol < 0 || fabs(2*asol) < 1E-13) return NAN;

    return (2*csol) / (-bsol + sgn*sqrt(bsol*bsol - 4*asol*csol));
}

// The transfer distance is obtained by solving the system of equations
//     x_s = x_t + t*l
//     y_s = y_t + t*m
//     z_s = t*n
// where the ray starts at (x_t, y_t, 0) with direction cosines (l, m, n) and is 
// propagated to the surface at (x_s, y_s, z_s). This function computes an analytic
// expression for t.
double Conic3DTransfer(double xt, double yt, double l, double m, double n, double cv, double k, double alpha, double beta, double gamma) {

    alpha = ConvertAngle(alpha) * M_PI/180;
    beta = ConvertAngle(beta) * M_PI/180;
    gamma = ConvertAngle(gamma) * M_PI/180;

    // The value of gamma determines which solution of the quadratic to use, as
    // with the sag
    int sgn;
    if (fabs(gamma) <= M_PI/2) sgn = 1;
    else sgn = -1;
        
    double x0, y0;
    Conic3DOffAxisDistance(&x0, &y0, cv, alpha, beta);

    double cosg = cos(gamma), cosg2 = cosg*cosg;
    double sing = sin(gamma), sing2 = sing*sing;
    
    double dsol = cv + cv*k*(n*n*cosg2 + l*l*sing2 + 2*l*n*sing*cosg);
    double fsol = cv*l*(xt*(1+k*sing2)+x0*cosg) - l*sing + cv*m*(yt-y0) + cv*n*sing*(k*xt*cosg-x0) - n*cosg;
    double gsol = cv*(xt*xt + x0*x0 + (yt-y0)*(yt-y0)) + 2*cv*x0*xt*cosg + xt*sing*(cv*k*xt*sing-2);

    // Ray missed this surface
    if (fsol*fsol - dsol*gsol < 0) return NAN;
        
    return gsol/(-fsol + sgn*sqrt(fsol*fsol - dsol*gsol));
}


// The normal vectors are obtained by taking the gradient of the sag function.
// We need these along with the transfer distance to perform the ray refraction.
void Conic3DSurfaceNormal(double *ln, double *mn, double *nn, double x, double y, double cv, double k, double alpha, double beta, double gamma, int normalize) {
    alpha = ConvertAngle(alpha) * M_PI/180;
    beta = ConvertAngle(beta) * M_PI/180;
    gamma = ConvertAngle(gamma) * M_PI/180;

    // The sag depends on the value of gamma, so preseve the sign here
    int sgn;
    if (fabs(gamma) <= M_PI/2) sgn = 1;
    else sgn = -1;

    double x0, y0;
    Conic3DOffAxisDistance(&x0, &y0, cv, alpha, beta);
    
    double cosg = cos(gamma), cosg2 = cosg*cosg;
    double sing = sin(gamma), sing2 = sing*sing;
    double asol = cv*(sing2 + (k+1)*cosg2);
    double bsol = 2*cv*sing*(x*k*cosg - x0) - 2*cosg;
    double csol = cv*k*x*x*sing2 - 2*x*sing + 2*cv*x0*x*cosg + cv*(x*x + x0*x0 + (y-y0)*(y-y0));

    double arg0 = bsol*bsol - 4*asol*csol;
    if (arg0 <= 0) {
        *ln = *mn = *nn = NAN;
        return;
    }
    double eta = sqrt(arg0);
    double denom = -bsol + sgn*eta;

    // Partial derivative with respect to x
    double arg1 = 4 * (cv*x*(1+k*sing2) + cv*x0*cosg - sing) / denom;
    double arg2 = 4 * csol / (denom*denom);
    double arg3 = -cv*k*sing*cosg;
    double arg4 = cv*k*bsol*sing*cosg - 2*asol*(cv*x*(1+k*sing2) + cv*x0*cosg - sing);
    double dervx = arg1 - arg2 * (arg3 + sgn * arg4 / eta);

    // ...with respect to y
    arg1 = 4*cv*(y-y0) / denom;
    arg2 = 2*asol*csol / (eta*denom);
    double dervy = arg1 * (1 + sgn*arg2);
    
    // Surface normals should be normalized. But if we do this, we lose the magnitude
    // of the derivatives. So we can choose to normalize or not.
    double norm = 1;
    if (normalize) norm = sqrt(dervx*dervx + dervy*dervy + 1);

    // d/dz is always equal to -1, no need to calculate it
    // The sign in Zemax's example seems to be opposite of Shannon (1997)...
    *ln = dervx / norm;
    *mn = dervy / norm;
    *nn = -1 / norm;
}

// Wrapper function to make accessing the partial derivative about x easier. This
// is only used for the critical point calculation. 
double Conic3DDervX(double x, double y, double cv, double k, double alpha, double beta, double gamma) {
    double dervx, dervy, dervz;
    Conic3DSurfaceNormal(&dervx, &dervy, &dervz, x, y, cv, k, alpha, beta, gamma, 0);
    return dervx;
}

// The critical points are used to find the bounded extrema within a slice. By
// symmetry, the y-coordinate of the critical point *must* be equal to the amount
// we shifted the conic (y0).
//
// The x-coordinate is less obvious because there is a rotation by gamma on top of the
// shift x0. Solving for dz/dx = 0 analytically is not feasible, so I use a secant
// method to find the root. This should be fast because the derivative is monotonic
// and well-behaved.
// 
// Also, we know that there can only be one critical point because an unrotated
// conic has only one critical point. For there to be two, the sag would not be
// a function anymore (multiple values of the sag for the same (x, y)).
void Conic3DCriticalXY(double *xc, double *yc, double cv, double k, double alpha, double beta, double gamma) {
    double tol = 1E-13;    // Tolerance for accepting root
    int niter_max = 50;    // Max number of iterations - under normal circumstances
                           // should converge quickly (10 iterations or less)

    // If the curvature is very small this is basically a plane. In this case
    // do not attempt to find critical points because the derivative is constant.
    if (fabs(cv) < 1E-13) {
        *yc = NAN; *xc = NAN;
        return;
    }

    alpha = ConvertAngle(alpha) * M_PI/180;
    beta = ConvertAngle(beta) * M_PI/180;
    gamma = ConvertAngle(gamma) * M_PI/180;

    double x0, y0;
    Conic3DOffAxisDistance(&x0, &y0, cv, alpha, beta);

    // Perform secant method to find xc; it should be around x0 if gamma is not
    // huge. Use 10% the radius of curvature as initial guesses to compute the
    // first secant line.
    double xc0 = x0 + 0.1 * 1/cv;
    double xc1 = x0 - 0.1 * 1/cv;
    double xc2, dervx0, dervx1;
    double err = 1;
    int i = 1;
    while (err < tol && i < niter_max) {
        dervx0 = Conic3DDervX(xc0, -y0, cv, k, alpha, beta, gamma);
        dervx1 = Conic3DDervX(xc1, -y0, cv, k, alpha, beta, gamma);
        xc2 = (xc0 * dervx1 - xc1 * dervx0) / (dervx1 - dervx0);
        err = fabs(xc2 - xc1);
        xc0 = xc1; xc1 = xc2;
    }
    // If it was unable to converge for some reason set xc to NAN...
    *xc = (err < tol) ? xc2 : NAN;
    *yc = -y0;
}

// Tilted plane solutions (cv = 0). The functions below match the call signatures
// of the functions for the conic surface. This allows us to pass these functions
// as arguments to other functions that are used to define the behavior of
// the entire image slicer (see slicer_generation.c).

// alpha, beta, and gamma are implemented a bit differently here because we cannot
// use the previous definition of the off-axis distance. Instead, we set alpha to
// be an rotation about the global y-axis, and then beta + gamma to be about the
// global x-axis. These are extrinsic rotations.
//
// The reason why beta and gamma are kept separate is that alpha and beta are
// used to define the behavior of a given "section" of the image slicer. gamma
// defines the tilt of each individual slice within that section.
double TiltedPlaneSag(double x, double y, double cv, double k, double alpha, double beta, double gamma) {
    alpha = ConvertAngle(alpha) * M_PI/180;
    beta = ConvertAngle(beta) * M_PI/180;
    gamma = ConvertAngle(gamma) * M_PI/180;
    
    double sina = sin(alpha), cosa = cos(alpha);
    double sinbg = sin(beta + gamma), cosbg = cos(beta + gamma);
    // Cap to prevent these from exploding
    if (fabs(cosa) < 1E-13) cosa = 1E-13;
    if (fabs(cosbg) < 1E-13) cosbg = 1E-13;

    return x * sinbg / (cosa*cosbg) - y * sina / cosa;
}

// See explanation for the transfer distance for the conic surface.
double TiltedPlaneTransfer(double xt, double yt, double l, double m, double n, double cv, double k, double alpha, double beta, double gamma) {
    alpha = ConvertAngle(alpha) * M_PI/180;
    beta = ConvertAngle(beta) * M_PI/180;
    gamma = ConvertAngle(gamma) * M_PI/180;
    
    double sina = sin(alpha), cosa = cos(alpha);
    double sinbg = sin(beta + gamma), cosbg = cos(beta + gamma);
    // Cap to prevent these from exploding
    if (fabs(cosa) < 1E-13) cosa = 1E-13;
    if (fabs(cosbg) < 1E-13) cosbg = 1E-13;

    double arg1 = xt * sinbg / (cosa * cosbg) - yt * sina / cosa;
    double arg2 = n - l * sinbg / (cosa * cosbg) + m * sina / cosa;
    
    if (fabs(arg2) < 1E-13) return NAN;
    return arg1 / arg2;
}

// See explanation for the surface normal for the conic surface.
void TiltedPlaneSurfaceNormal(double *ln, double *mn, double *nn, double x, double y, double cv, double k, double alpha, double beta, double gamma, int normalize) {
    alpha = ConvertAngle(alpha) * M_PI/180;
    beta = ConvertAngle(beta) * M_PI/180;
    gamma = ConvertAngle(gamma) * M_PI/180;
    
    double sina = sin(alpha), cosa = cos(alpha);
    double sinbg = sin(beta + gamma), cosbg = cos(beta + gamma);
    // Cap to prevent these from exploding
    if (fabs(cosa) < 1E-13) cosa = 1E-13;
    if (fabs(cosbg) < 1E-13) cosbg = 1E-13;

    double dervx = sinbg / (cosa * cosbg);
    double dervy = -sina / cosa;

    double norm = 1;
    if (normalize) norm = sqrt(dervx*dervx + dervy*dervy + 1);

    *ln = dervx / norm;
    *mn = dervy / norm;
    *nn = -1 / norm;
}

// Planes have no critical points, so this does nothing.
void TiltedPlaneCriticalXY(double *xc, double *yc, double cv, double k, double alpha, double beta, double gamma) {
     *xc = NAN; *yc = NAN;
}