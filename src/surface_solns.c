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

double CalcQuadraticDerv(int sgn, double A, double B, double C,
                         double dA, double dB, double dC) {

    double discrim = B * B - A * C;

    // Guard against invalid sqrt
    if (discrim < 0.0) return NAN;

    double sqrt_disc = sqrt(discrim);
    double Eta = -B + sgn * sqrt_disc;

    // Derivative of Eta
    double dEta = -dB + sgn * (2.0 * B * dB - C * dA - A * dC) / (2.0 * sqrt_disc);

    // Final derivative expression
    return (Eta * dC - C * dEta) / (Eta * Eta);
}


// 2D conic solutions - planes (cv -> 0) are handled differently, see next section

// Determine the off-axis distance. If the curvature is near-zero, this is
// basically a plane. Set a limit to how small cv can be. In practice the
// user should not set a non-zero off-axis distance if the surface is a plane.
void Conic2DOffAxisDistance(double *x0, double *y0, double c, double k, double alpha, double beta) {
    if (fabs(c) < 1E-13) c = 1E-13;

    // y0: off-axis distance in y direction
    if (alpha == 0.0) {
        *y0 = 0.0;
    } else {
        double tana = tan(alpha);
        double num = (k - 1.0) + sqrt(4.0 + tana * tana * (3.0 - k));
        double den = tana * tana + (1.0 + k);
        *y0 = (tana / (2.0 * c)) * (num / den);
    }

    // x0: off-axis distance in x direction
    if (beta == 0.0) {
        *x0 = 0.0;
    } else {
        double tanb = tan(beta);
        double num = (k - 1.0) + sqrt(4.0 + tanb * tanb * (3.0 - k));
        double den = tanb * tanb + (1.0 + k);
        *x0 = (tanb / (2.0 * c)) * (num / den);
    }
}

// This is the sag for a conic surface translated by x0 and y0, and then rotated
// about the global y-axis. The sag for this surface is solved by doing the
// transformations backward, i.e., undoing the rotation and then the translation.
// The reason for this is that the sag is easier to solve for this way, but it
// is equivalent. Refer to the derivation for more informaiton.
double Conic2DSag(double x, double y, SLICE_PARAMS pslice) {
    x -= pslice.u;

    double alpha = ConvertAngle(pslice.alpha) * M_PI / 180.0;
    double beta  = ConvertAngle(pslice.beta)  * M_PI / 180.0;
    double gamma = ConvertAngle(pslice.gamma) * M_PI / 180.0;

    int sgn = (fabs(gamma) <= M_PI / 2.0) ? 1 : -1;
    if (pslice.k < -1.0) sgn *= -1;

    double x0, y0;
    Conic2DOffAxisDistance(&x0, &y0, pslice.cv, pslice.k, alpha, beta);

    double v1 = pslice.syx + x0;
    double v2 = pslice.u - pslice.syx;
    double v3 = pslice.syz + pslice.zp;

    double cosg = cos(gamma);
    double cos2g = cos(2.0 * gamma);
    double sing = sin(gamma);
    double sin2g = sin(2.0 * gamma);

    double A, B, C;
    double cv = pslice.cv;

    if (pslice.k == -1.0) {
        A = cv * sing * sing;

        B = -(
            cosg
            + cv * (v2 + x) * cosg * sing
            + cv * sing * (v1 + v3 * sing)
        );

        C = (
            -2.0 * pslice.syz
            + cv * (v1 * v1 + (y + y0) * (y + y0))
            + 2.0 * (v3 + cv * v1 * (v2 + x)) * cosg
            + cv * (v2 + x) * (v2 + x) * cosg * cosg
            - 2.0 * (v2 - cv * v1 * v3 + x) * sing
            + cv * v3 * v3 * sing * sing
            + cv * v3 * (v2 + x) * sin2g
        );
    } else {
        A = cv * (1.0 + pslice.k) * (2.0 + pslice.k + pslice.k * cos2g);

        B = -(1.0 + pslice.k) * (
            -2.0 * (-1.0 + cv * (1.0 + pslice.k) * pslice.syz) * cosg
            + cv * (
                (2.0 + pslice.k) * v3
                + pslice.k * v3 * cos2g
                + 2.0 * v1 * sing
                - pslice.k * (v2 + x) * sin2g
            )
        );

        C = (1.0 + pslice.k) * (
            -4.0 * pslice.syz
            + 2.0 * cv * ((1.0 + pslice.k) * pslice.syz * pslice.syz + v1 * v1 + (y + y0) * (y + y0))
            + cv * (2.0 + pslice.k) * (v2 * v2 + v3 * v3 + 2.0 * v2 * x + x * x)
            + 4.0 * (v3 - cv * (1.0 + pslice.k) * pslice.syz * v3 + cv * v1 * (v2 + x)) * cosg
            - cv * pslice.k * (v2 - v3 + x) * (v2 + v3 + x) * cos2g
            + 4.0 * ((-1.0 + cv * (1.0 + pslice.k) * pslice.syz) * v2 + cv * v1 * v3 - x + cv * (1.0 + pslice.k) * pslice.syz * x) * sing
            - 2.0 * cv * pslice.k * v3 * (v2 + x) * sin2g
        );
    }

    double discrim = B * B - A * C;

    if (discrim < 0.0) return 0.0;

    return C / (-B + sgn * sqrt(discrim));
}


// The transfer distance is obtained by solving the system of equations
//     x_s = x_t + t*l
//     y_s = y_t + t*m
//     z_s = t*n
// where the ray starts at (x_t, y_t, 0) with direction cosines (l, m, n) and is 
// propagated to the surface at (x_s, y_s, z_s). This function computes an analytic
// expression for t.
double Conic2DTransfer(double xt, double yt, double l, double m, double n, SLICE_PARAMS pslice) {
    double alpha = ConvertAngle(pslice.alpha) * M_PI / 180.0;
    double beta  = ConvertAngle(pslice.beta)  * M_PI / 180.0;
    double gamma = ConvertAngle(pslice.gamma) * M_PI / 180.0;

    int sgn = (fabs(gamma) <= M_PI / 2.0) ? 1 : -1;
    if (pslice.k < -1.0) sgn *= -1;

    double x0, y0;
    Conic2DOffAxisDistance(&x0, &y0, pslice.cv, pslice.k, alpha, beta);

    double v1 = pslice.syx + x0;
    double v2 = pslice.u - pslice.syx;
    double v3 = pslice.syz + pslice.zp;

    double cosg = cos(gamma);
    double sing = sin(gamma);
    double cos2g = cos(2.0 * gamma);
    double sin2g = sin(2.0 * gamma);

    double D, F, G;
    double cv = pslice.cv;

    if (pslice.k == -1.0) {
        D = cv * (m * m + (l * cosg - n * sing) * (l * cosg - n * sing));

        F = (
            cv * m * (y0 + yt)
            + cv * l * (v2 + xt) * cosg * cosg
            - sing * (l + cv * n * v1 + cv * n * v3 * sing)
            - cosg * (n - cv * l * v1 + cv * (-l * v3 + n * (v2 + xt)) * sing)
        );

        G = (
            -2.0 * pslice.syz
            + cv * (v1 * v1 + (y0 + yt) * (y0 + yt))
            + 2.0 * (v3 + cv * v1 * (v2 + xt)) * cosg
            + cv * (v2 + xt) * (v2 + xt) * cosg * cosg
            - 2.0 * (v2 - cv * v1 * v3 + xt) * sing
            + cv * v3 * v3 * sing * sing
            + cv * v3 * (v2 + xt) * sin2g
        );
    } else {
        D = 2.0 * cv * pslice.k * (n * cosg + l * sing) * (n * cosg + l * sing)
            + 2.0 * cv * (l * l + m * m + n * n);

        F = (
            cv * (-(2.0 + pslice.k) * n * v3 + (2.0 + pslice.k) * l * (v2 + xt) + 2.0 * m * (y0 + yt))
            + 2.0 * (n * (-1.0 + cv * (1.0 + pslice.k) * pslice.syz) + cv * l * v1) * cosg
            - cv * pslice.k * (n * v3 + l * (v2 + xt)) * cos2g
            + 2.0 * (l * (-1.0 + cv * (1.0 + pslice.k) * pslice.syz) - cv * n * v1) * sing
            + cv * pslice.k * (-l * v3 + n * (v2 + xt)) * sin2g
        );

        G = (
            -4.0 * pslice.syz
            + 2.0 * cv * ((1.0 + pslice.k) * pslice.syz * pslice.syz + v1 * v1 + (y0 + yt) * (y0 + yt))
            + cv * (2.0 + pslice.k) * (v2 * v2 + v3 * v3 + 2.0 * v2 * xt + xt * xt)
            + 4.0 * (v3 - cv * (1.0 + pslice.k) * pslice.syz * v3 + cv * v1 * (v2 + xt)) * cosg
            - cv * pslice.k * (v2 - v3 + xt) * (v2 + v3 + xt) * cos2g
            + 4.0 * ((-1.0 + cv * (1.0 + pslice.k) * pslice.syz) * v2 + cv * v1 * v3 - xt + cv * (1.0 + pslice.k) * pslice.syz * xt) * sing
            - 2.0 * cv * pslice.k * v3 * (v2 + xt) * sin2g
        );
    }

    double discrim = F * F - D * G;

    if (discrim < 0.0)
        return NAN;

    return G / (-F + sgn * sqrt(discrim));
}


// The normal vectors are obtained by taking the gradient of the sag function.
// We need these along with the transfer distance to perform the ray refraction.
void Conic2DSurfaceNormal(double *ln, double *mn, double *nn, double x, double y, SLICE_PARAMS pslice, int normalize) {
    double alpha = ConvertAngle(pslice.alpha) * M_PI / 180.0;
    double beta  = ConvertAngle(pslice.beta)  * M_PI / 180.0;
    double gamma = ConvertAngle(pslice.gamma) * M_PI / 180.0;

    int sgn = (fabs(gamma) <= M_PI / 2.0) ? 1 : -1;
    if (pslice.k < -1.0) sgn *= -1;

    double x0, y0;
    Conic2DOffAxisDistance(&x0, &y0, pslice.cv, pslice.k, alpha, beta);

    double v1 = pslice.syx + x0;
    double v2 = pslice.u - pslice.syx;
    double v3 = pslice.syz + pslice.zp;

    double cosg = cos(gamma);
    double sing = sin(gamma);
    double cos2g = cos(2.0 * gamma);
    double sin2g = sin(2.0 * gamma);

    double A, B, C;
    double dAx = 0, dAy = 0, dBx = 0, dBy = 0, dCx = 0, dCy = 0;
    double cv = pslice.cv;
    double k = pslice.k;

    if (k == -1.0) {
        // Parabolic case
        A = cv * sing * sing;

        B = -(
            cosg
            + cv * (v2 + x) * cosg * sing
            + cv * sing * (v1 + v3 * sing)
        );

        C = (
            -2.0 * pslice.syz
            + cv * (v1 * v1 + (y + y0) * (y + y0))
            + 2.0 * (v3 + cv * v1 * (v2 + x)) * cosg
            + cv * (v2 + x) * (v2 + x) * cosg * cosg
            - 2.0 * (v2 - cv * v1 * v3 + x) * sing
            + cv * v3 * v3 * sing * sing
            + cv * v3 * (v2 + x) * sin2g
        );

        dBx = -cv * cosg * sing;
        dCx = -2.0 * sing + 2.0 * cv * cosg * (v1 + (v2 + x) * cosg + v3 * sing);
        dCy = 2.0 * cv * (y + y0);
    } else {
        // Non-parabolic conic
        A = cv * (1.0 + k) * (2.0 + k + k * cos2g);

        B = - (1.0 + k) * (
            -2.0 * (-1.0 + cv * (1.0 + k) * pslice.syz) * cosg
            + cv * (
                (2.0 + k) * v3
                + k * v3 * cos2g
                + 2.0 * v1 * sing
                - k * (v2 + x) * sin2g
            )
        );

        C = (1.0 + k) * (
            -4.0 * pslice.syz
            + 2.0 * cv * ((1.0 + k) * pslice.syz * pslice.syz + v1 * v1 + (y + y0) * (y + y0))
            + cv * (2.0 + k) * (v2 * v2 + v3 * v3 + 2.0 * v2 * x + x * x)
            + 4.0 * (v3 - cv * (1.0 + k) * pslice.syz * v3 + cv * v1 * (v2 + x)) * cosg
            - cv * k * (v2 - v3 + x) * (v2 + v3 + x) * cos2g
            + 4.0 * ((-1.0 + cv * (1.0 + k) * pslice.syz) * v2 + cv * v1 * v3 - x + cv * (1.0 + k) * pslice.syz * x) * sing
            - 2.0 * cv * k * v3 * (v2 + x) * sin2g
        );

        dBx = -cv * k * (1.0 + k) * sin2g;
        dCx = (1.0 + k) * (
            2.0 * cv * (2.0 + k) * (v2 + x)
            + 4.0 * cv * v1 * cosg
            - 2.0 * cv * k * (v2 + x) * cos2g
            - 2.0 * cv * k * (v2 - v3 + x) * cos2g
            - 4.0 * (-1.0 + cv * (1.0 + k) * pslice.syz) * sing
            + 2.0 * cv * k * v3 * sin2g
        );
        dCy = 4.0 * cv * (1.0 + k) * (y + y0);
    }

    double dervx = calc_quadratic_derv(sgn, A, B, C, dAx, dBx, dCx);
    double dervy = calc_quadratic_derv(sgn, A, B, C, dAy, dBy, dCy);

    double norm = 1.0;
    if (normalize) {
        norm = sqrt(dervx * dervx + dervy * dervy + 1.0);
    }

    *ln = dervx / norm;
    *mn = dervy / norm;
    *nn = -1.0 / norm;
}


// Wrapper function to make accessing the partial derivative about x easier. This
// is only used for the critical point calculation. 
double Conic2DDervX(double x, double y, SLICE_PARAMS pslice) {
    double dervx, dervy, dervz;
    Conic2DSurfaceNormal(&dervx, &dervy, &dervz, x, y, pslice, 0);
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
void Conic2DCriticalXY(double *xc, double *yc, SLICE_PARAMS pslice) {
    double tol = 1E-13;    // Tolerance for accepting root
    int niter_max = 50;    // Max number of iterations - under normal circumstances
                           // should converge quickly (10 iterations or less)

    double cv = pslice.cv;
    double k = pslice.k;
    double alpha = ConvertAngle(pslice.alpha) * M_PI / 180.0;
    double beta  = ConvertAngle(pslice.beta)  * M_PI / 180.0;
    double gamma = ConvertAngle(pslice.gamma) * M_PI / 180.0;

    // If the curvature is very small this is basically a plane. In this case
    // do not attempt to find critical points because the derivative is constant.
    if (fabs(cv) < 1E-13) {
        *yc = NAN; *xc = NAN;
        return;
    }
    
    double x0, y0;
    Conic2DOffAxisDistance(&x0, &y0, cv, k, alpha, beta);

    // Perform secant method to find xc; it should be around x0 if gamma is not
    // huge. Use 5% the radius of curvature as initial guesses to compute the
    // first secant line.
    double xc0 = x0 + 0.05 * 1/pslice.cv;
    double xc1 = x0 - 0.05 * 1/cv;
    double xc2, dervx0, dervx1;
    double err = 1;
    int i = 1;
    while (err < tol && i < niter_max) {
        dervx0 = Conic2DDervX(xc0, -y0, pslice);
        dervx1 = Conic2DDervX(xc1, -y0, pslice);
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
double TiltedPlaneSag(double x, double y, SLICE_PARAMS p) {
    // Convert angles to radians
    double cv = p.cv;
    double k = p.k;
    double alpha = ConvertAngle(p.alpha) * M_PI / 180.0;
    double beta  = ConvertAngle(p.beta)  * M_PI / 180.0;
    double gamma = ConvertAngle(p.gamma) * M_PI / 180.0;

    double cosa  = cos(alpha);
    double tana  = tan(alpha);
    double cosbg = cos(beta + gamma);
    double tanbg = tan(beta + gamma);

    // Cap to prevent explosion near grazing incidence
    if (fabs(cosa) < 1e-13)  cosa = 1e-13;
    if (fabs(cosbg) < 1e-13) cosbg = 1e-13;

    double seca   = 1.0 / cosa;
    double secbg  = 1.0 / cosbg;

    double z = (
        secbg * (p.sxz - p.syz - p.sxz * seca + (y - p.sxy) * tana)
        - (x - p.syx + p.u) * tanbg
        + p.syz
        + p.zp
    );

    return z;
}


// See explanation for the transfer distance for the conic surface.
double TiltedPlaneTransfer(double xt, double yt, double l, double m, double n, SLICE_PARAMS pslice) {
    double alpha = ConvertAngle(pslice.alpha) * M_PI / 180.0;
    double beta  = ConvertAngle(pslice.beta)  * M_PI / 180.0;
    double gamma = ConvertAngle(pslice.gamma) * M_PI / 180.0;

    double cosa  = cos(alpha);
    double tana  = tan(alpha);
    double cosbg = cos(beta + gamma);
    double tanbg = tan(beta + gamma);

    // Cap to prevent these from exploding
    if (fabs(cosa) < 1e-13)  cosa = 1e-13;
    if (fabs(cosbg) < 1e-13) cosbg = 1e-13;

    double seca  = 1.0 / cosa;
    double secbg = 1.0 / cosbg;

    double arg1 = secbg * (pslice.sxz - pslice.syz - pslice.sxz * seca + (yt - pslice.sxy) * tana);
    double arg2 = (xt - pslice.syx + pslice.u) * tanbg;
    double arg3 = pslice.syz + pslice.zp;

    double den = n - m * secbg * tana + l * tanbg;

    if (fabs(den) < 1e-13) {
        return NAN;  // handle division by zero gracefully
    }

    return (arg1 - arg2 + arg3) / den;
}


// See explanation for the surface normal for the conic surface.
void TiltedPlaneSurfaceNormal(double *ln, double *mn, double *nn, double x, double y, SLICE_PARAMS pslice, int normalize) {
    double alpha = ConvertAngle(pslice.alpha) * M_PI / 180.0;
    double beta  = ConvertAngle(pslice.beta)  * M_PI / 180.0;
    double gamma = ConvertAngle(pslice.gamma) * M_PI / 180.0;

    double cosa  = cos(alpha);
    double tana  = tan(alpha);
    double cosbg = cos(beta + gamma);
    double tanbg = tan(beta + gamma);

    // Cap to prevent these from exploding
    if (fabs(cosa) < 1e-13)  cosa  = 1e-13;
    if (fabs(cosbg) < 1e-13) cosbg = 1e-13;

    double secbg = 1.0 / cosbg;

    double dervx = -tanbg;
    double dervy = secbg * tana;

    double norm = 1.0;
    if (normalize) {
        norm = sqrt(dervx * dervx + dervy * dervy + 1.0);
    }

    *ln = dervx / norm;
    *mn = dervy / norm;
    *nn = -1.0;
}

// Planes have no critical points, so this does nothing.
void TiltedPlaneCriticalXY(double *xc, double *yc, SLICE_PARAMS pslice) {
     *xc = NAN; *yc = NAN;
}