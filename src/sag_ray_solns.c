#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "ifu_helpers.h"

/*
Solutions for the sag, critical points, ray transfer distance, and surface normals
for different surface types.

Ellen Lee
*/


double ConvertAngle(double t) {
    t = fmod(t, 360);
    if (t > 180) {
        return t - 360;
    }
    return t;
}

// 3D conic sag calculations - planes (cv -> 0) have different definitions of
// alpha and beta and are handled separately.


double Conic3DSag(double x, double y, double cv, double k, double alpha, double beta, double gamma) {
    // The value of gamma determines which solution of the quadratic to use
    alpha = ConvertAngle(alpha) * M_PI/180;
    beta = ConvertAngle(beta) * M_PI/180;
    gamma = ConvertAngle(gamma) * M_PI/180;
    int sgn;
    if (fabs(gamma) <= M_PI/2) {sgn = 1;}
    else {sgn = -1;}

    // Determine the off-axis distance. If the curvature is near-zero, this is
    // basically a plane. Set a limit to how small cv can be. In practice the
    // user should not set a non-zero off-axis distance if the surface is a plane.
    if (fabs(cv) < 1E-13) cv = 1E-13;

    double y0 = sin(alpha) / ( cv * (1 + cos(alpha)) );
    double x0 = sin(beta) / ( cv * (1 + cos(beta)) );
    y = y + y0;
    x = x + x0;

    // Rotate about the x-axis
    double cosg = cos(gamma);
    double cosg2 = cosg*cosg;
    double sing = sin(gamma);
    double sing2 = sing*sing;

    double asol = cv*(sing2 + (k+1)*cosg2);
    double bsol = -2*cosg*(x*k*cv*sing + 1);
    double csol = 2*x*sing + cv*(y*y + x*x*cosg2 + (k+1)*x*x*sing2);

    if (bsol*bsol - 4*asol*csol < 0 || fabs(2*asol) < 1E-10) {
        // Invalid solution, set sag to 0
        return 0;
    }
    return (2*csol) / (-bsol + sgn*sqrt(bsol*bsol - 4*asol*csol));
}


void Conic3DCriticalXY(double *xc1, double *xc2, double *yc, double cv, double k, double alpha, double beta, double gamma) {

    alpha = ConvertAngle(alpha) * M_PI/180;
    beta = ConvertAngle(beta) * M_PI/180;
    gamma = ConvertAngle(gamma) * M_PI/180;

    // If the curvature is very small this is basically a plane. In this case
    //  there are no critical points because the derivative is a constant.
    if (fabs(cv) < 1E-13) {
        *yc = NAN; *xc1 = NAN; *xc2 = NAN;
        return;
    }

    // Determine off axis distances.
    double y0 = sin(alpha) / ( cv * (1 + cos(alpha)) );
    double x0 = sin(beta) / ( cv * (1 + cos(beta)) );

    *yc = -y0;

    double sing = sin(gamma);
    double sin3g = sin(3 * gamma);
    double cosg = cos(gamma);
    double cos2g = cos(2 * gamma);

    if (k == -1) {
        *xc1 = -x0 + ( sing*(5 + cv*cv*y0*y0) + sin3g*(1 + cv*cv*y0*y0) ) / (cosg*cosg*-8*cv);
        *xc2 = NAN; // No second solution
    }
    else {
        double U = -sing / (cv*(1+k)*(-2-k+k*cos2g));
        double V = -2 - k + k*cos2g;
        double W = cosg*k* sqrt( (-2+2*cv*cv*y0*y0*(1+k)) * V );

        *xc1 = U * (V - W) - x0;
        *xc2 = U * (V + W) - x0;
    }
}


//