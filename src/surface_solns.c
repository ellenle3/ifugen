#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "surface_solns.h"

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

// 3D conic solutions - planes (cv -> 0) are handled differently, see next section

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
    double cosg = cos(gamma), cosg2 = cosg*cosg;
    double sing = sin(gamma), sing2 = sing*sing;

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

    double sing = sin(gamma), sin3g = sin(3 * gamma);
    double cosg = cos(gamma), cos2g = cos(2 * gamma);

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

double Conic3DTransfer(double xt, double yt, double l, double m, double n, double cv, double k, double alpha, double beta, double gamma) {
    alpha = ConvertAngle(alpha) * M_PI/180;
    beta = ConvertAngle(beta) * M_PI/180;
    gamma = ConvertAngle(gamma) * M_PI/180;
    int sgn;
    if (fabs(gamma) <= M_PI/2) {sgn = 1;}
    else {sgn = -1;}
        
    // Determine the off-axis distance.
    if (fabs(cv) < 1E-13) {cv = 1E-13;}
        
    double y0 = sin(alpha) / ( cv * (1 + cos(alpha)) );
    double x0 = sin(beta) / ( cv * (1 + cos(beta)) );
    yt = yt + y0;
    xt = xt + x0;

    double cosg = cos(gamma), cosg2 = cosg*cosg;
    double sing = sin(gamma), sing2 = sing*sing;
    
    double xi = cosg2 + (1+k)*sing2;
    double asol = cv*(sing2 + (1+k)*cosg2);
    double dsol = asol*n*n - 2*l*n*k*cv*sing*cosg + cv*l*l*xi + m*m*cv;
    double fsol = -xt*n*k*cv*sing*cosg - n*cosg + l*sing + cv*xt*l*xi + yt*m*cv;
    double gsol = 2*xt*sing + cv*xt*xt*xi + cv*yt*yt;

    // Ray missed this surface
    if (fsol*fsol - dsol*gsol < 0) {return NAN;}
        
    return gsol/(-fsol + sgn*sqrt(fsol*fsol - dsol*gsol));
}

void Conic3DSurfaceNormal(double *ln, double *mn, double *nn, double x, double y, double cv, double k, double alpha, double beta, double gamma) {
    // The value of gamma determines which solution of the quadratic to use
    alpha = ConvertAngle(alpha) * M_PI/180;
    beta = ConvertAngle(beta) * M_PI/180;
    gamma = ConvertAngle(gamma) * M_PI/180;
    int sgn;
    if (fabs(gamma) <= M_PI/2) {sgn = 1;}
    else {sgn = -1;}

    if (fabs(cv) < 1E-13) cv = 1E-13;

    double y0 = sin(alpha) / ( cv * (1 + cos(alpha)) );
    double x0 = sin(beta) / ( cv * (1 + cos(beta)) );
    y = y + y0;
    x = x + x0;
    
    double sq2 = sqrt(2);
    double cosg = cos(gamma), cos2g = cos(2*gamma);
    double sing = sin(gamma), sin2g = sin(2*gamma);
    double arg0 = 1 - cv*cv*( 2*(1+k)*x*x + (2+k)*y*y ) + (1-cv*cv*k*y*y)*cos2g - 4*cv*x*sing;

    // Surface is not defined here
    if (arg0 <= 0) {
        *ln = NAN; *mn = NAN; *nn = NAN;
        return;
    }

    double eta = sqrt(arg0);
    double num = ( eta*sq2 + sgn*2*cosg*(1+cv*k*x*sing) );
    double psi = 2 * cv / (num*num);
    
    // Partial derivative with respect to x
    double arg1 = cv*y*y + cv*x*x*(2+k-k*cos2g)/2 + 2*x*sing;
    double arg2 = k*sin2g - sgn*2*sq2*( cv*(1+k)*x + sing )/eta;
    double arg3 = 2*cv*(2+k)*x - 2*cv*k*x*cos2g + 4*sing;
    double arg4 = 2*cosg + cv*k*x*sin2g + sgn*eta*sq2;
    double dervx = -psi * arg1 * arg2 + arg3 / arg4;

    // ...with respect to y
    double dervy;
    if (sgn == -1) {dervy = -sq2 * cv * y / eta;}
    else {
        arg1 = 2*eta*sq2 + 4*cosg*(1+cv*k*x*sing);
        arg2 = cv*sq2*(2+k+k*cos2g) / eta;
        arg3 = cv*y*y + cv*x*x*cosg*cosg + x*sing*(2+cv*(1+k)*x*sing);
        dervy = y * psi * ( arg1 + arg2 * arg3 );
    }
        
    double norm = sqrt(dervx*dervx + dervy*dervy + 1);

    // d/dz is always equal to -1, no need to calculate it
    //  The sign in Zemax's example seems to be opposite of Shannon (1997)...
    *ln = dervx / norm;
    *mn = dervy / norm;
    *nn = -1 / norm;
}

// Tilted plane solutions (cv = 0)

double TiltedPlaneSag(double x, double y, double cv, double k, double alpha, double beta, double gamma) {
    alpha = ConvertAngle(alpha) * M_PI/180;
    beta = ConvertAngle(beta) * M_PI/180;
    gamma = ConvertAngle(gamma) * M_PI/180;
    
    double sina = sin(alpha), cosa = cos(alpha);
    double sinbg = sin(beta + gamma), cosbg = cos(beta + gamma);
    // Cap to prevent these from exploding
    if (fabs(cosa) < 1E-13) {cosa = 1E-13;}
    if (fabs(cosbg) < 1E-13) {cosbg = 1E-13;}

    return x * sinbg / (cosa*cosbg) - y * sina / cosa;
}

void TiltedPlaneCriticalXY(double *xc1, double *xc2, double *yc, double cv, double k, double alpha, double beta, double gamma) {
    // Planes have no critical points
    *yc = NAN; *xc1 = NAN; *xc2 = NAN;
}

double TiltedPlaneTransfer(double xt, double yt, double l, double m, double n, double cv, double k, double alpha, double beta, double gamma) {
    alpha = ConvertAngle(alpha) * M_PI/180;
    beta = ConvertAngle(beta) * M_PI/180;
    gamma = ConvertAngle(gamma) * M_PI/180;
    
    double sina = sin(alpha), cosa = cos(alpha);
    double sinbg = sin(beta + gamma), cosbg = cos(beta + gamma);
    // Cap to prevent these from exploding
    if (fabs(cosa) < 1E-13) {cosa = 1E-13;}
    if (fabs(cosbg) < 1E-13) {cosbg = 1E-13;}

    double arg1 = xt * sinbg / (cosa * cosbg) - yt * sina / cosa;
    double arg2 = n - l * sinbg / (cosa * cosbg) + m * sina / cosa;
    
    if (fabs(arg2) < 1E-13) return NAN;
    return arg1 / arg2;
}

void TiltedPlaneSurfaceNormal(double *ln, double *mn, double *nn, double x, double y, double cv, double k, double alpha, double beta, double gamma) {
    alpha = ConvertAngle(alpha) * M_PI/180;
    beta = ConvertAngle(beta) * M_PI/180;
    gamma = ConvertAngle(gamma) * M_PI/180;
    
    double sina = sin(alpha), cosa = cos(alpha);
    double sinbg = sin(beta + gamma), cosbg = cos(beta + gamma);
    // Cap to prevent these from exploding
    if (fabs(cosa) < 1E-13) {cosa = 1E-13;}
    if (fabs(cosbg) < 1E-13) {cosbg = 1E-13;}

    double dervx = sinbg / (cosa * cosbg);
    double dervy = -sina / cosa;
    double norm = sqrt(dervx*dervx + dervy*dervy + 1);

    *ln = dervx / norm;
    *mn = dervy / norm;
    *nn = -1 / norm;
}