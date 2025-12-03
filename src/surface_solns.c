#define _USE_MATH_DEFINES
#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include "surface_solns.h"

/*
Solutions for the sag, critical points, ray transfer distance, and surface normals
for different surface types.

Ellen Lee
*/

// Bounds the angle to be between -180 and 180 degrees.
static double ConvertAngle(double t) {
    t = fmod(t, 360);
    if (t > 180) {
        return t - 360;
    }
    return t;
}

// 4x4 matrix multiplication
static void Mat4Mul(double out[4][4], const double A[4][4], const double B[4][4]) {
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            out[i][j] = 0.0;
            for (int k = 0; k < 4; k++)
                out[i][j] += A[i][k] * B[k][j];
        }
    }
}

// Multiplies 4x4 matrix with a 4x1 vector
static void Mat4VecMul(double out[3], const double M[4][4], const double v[3], const double w) {
    double hv[4] = { v[0], v[1], v[2], w };

    double r[4];
    for (int i = 0; i < 4; i++) {
        r[i] = M[i][0] * hv[0] + M[i][1] * hv[1] +
               M[i][2] * hv[2] + M[i][3] * hv[3];
    }

    out[0] = r[0];
    out[1] = r[1];
    out[2] = r[2];
}

static void Mat4AffineInverse(double Minv[4][4], const double M[4][4]) {
    // Upper-left 3×3 is orthonormal (rotation only)
    // So inverse rotation = transpose
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++)
            Minv[i][j] = M[j][i];
    }

    // Translation is last column of M
    double tx = M[0][3];
    double ty = M[1][3];
    double tz = M[2][3];

    // Inverse translation = -R^T * t
    Minv[0][3] = -(Minv[0][0] * tx + Minv[0][1] * ty + Minv[0][2] * tz);
    Minv[1][3] = -(Minv[1][0] * tx + Minv[1][1] * ty + Minv[1][2] * tz);
    Minv[2][3] = -(Minv[2][0] * tx + Minv[2][1] * ty + Minv[2][2] * tz);

    // Last row
    Minv[3][0] = 0;
    Minv[3][1] = 0;
    Minv[3][2] = 0;
    Minv[3][3] = 1;
}


static void NoSurfaceNormal(double *ln, double *mn, double *nn, double x, double y, SLICE_PARAMS pslice, int normalize) {
    *ln = NAN;
    *mn = NAN;
    *nn = NAN;
}

RAY_IN ConvertRayInToLocal(RAY_IN ray_in, SLICE_PARAMS pslice, TRANSFORMATION_FUNC transform_func) {
    double coords_global[3] = { ray_in.xt, ray_in.yt, ray_in.zt };  // starting coordinates
    double cosines_global[3] = { ray_in.l, ray_in.m, ray_in.n };    // direction cosines
    double coords_local[3];
    double cosines_local[3];

    // The ray is transformed backwards in the reference frame of the surface,
    // hence why converting the ray into local coordinates requires the inverse
    // transformation.
    transform_func(coords_local, coords_global, pslice, -1, 1);   // inverse, translate
    transform_func(cosines_local, cosines_global, pslice, -1, 0); // inverse, no translate

    RAY_IN ray_in_local = {
        coords_local[0],
        coords_local[1],
        coords_local[2],
        cosines_local[0],
        cosines_local[1],
        cosines_local[2]
    };
    return ray_in_local;
}

RAY_OUT ConvertRayOutToGlobal(RAY_OUT ray_out, SLICE_PARAMS pslice, TRANSFORMATION_FUNC transform_func) {
    double coords_local[3] = { ray_out.xs, ray_out.ys, ray_out.zs };   // surface coordinates
    double normals_local[3] = { ray_out.ln, ray_out.mn, ray_out.nn };  // surface normals
    double coords_global[3];
    double normals_global[3];

    transform_func(coords_global, coords_local, pslice, 1, 1);   // forward, translate
    transform_func(normals_global, normals_local, pslice, 1, 0); // forward, no translate
    
    RAY_OUT ray_out_global = {
        coords_global[0],
        coords_global[1],
        coords_global[2],
        ray_out.t,
        normals_global[0],
        normals_global[1],
        normals_global[2]
    };
    return ray_out_global;
}

RAY_OUT SliceRayTrace(RAY_IN ray_in, SLICE_PARAMS pslice, TRANSFER_DIST_FUNC transfer_dist_func,
    SURF_NORMAL_FUNC surface_normal_func, TRANSFORMATION_FUNC transform_func, int normalize) {

    // Convert the ray into local coordinates
    RAY_IN ray_in_local = ConvertRayInToLocal(ray_in, pslice, transform_func);

    double xt = ray_in_local.xt;
    double yt = ray_in_local.yt;
    double zt = ray_in_local.zt;
    double l  = ray_in_local.l;
    double m  = ray_in_local.m;
    double n  = ray_in_local.n;

    // Ray transfer in local coordinates
    double t = transfer_dist_func(xt, yt, zt, l, m, n, pslice);
    double xs = xt + t * l;
    double ys = yt + t * m;
    double zs = zt + t * n;
    double ln, mn, nn;
    surface_normal_func(&ln, &mn, &nn, xs, ys, pslice, normalize);

    RAY_OUT ray_out_local = {
        xs,
        ys,
        zs,
        t,
        ln,
        mn,
        nn
    };

    // Convert back to global coordinates
    RAY_OUT ray_out_global = ConvertRayOutToGlobal(ray_out_local, pslice, transform_func);
    return ray_out_global;
}

double SliceSag(double x, double y, SLICE_PARAMS pslice, TRANSFER_DIST_FUNC transfer_dist_func,
    TRANSFORMATION_FUNC transform_func) {
        // Sag is a special case of ray tracing where the ray is coming in parallel
        // to the z-axis at (x, y)
        RAY_IN ray_in = { x, y, 0.0, 0.0, 0.0, 1.0 };
        RAY_OUT ray_out = SliceRayTrace(ray_in, pslice, transfer_dist_func,
            &NoSurfaceNormal, transform_func, 0);
        return ray_out.zs;
}

void SliceSurfaceNormal(double* ln, double* mn, double* nn, double x, double y, SLICE_PARAMS pslice, TRANSFER_DIST_FUNC transfer_dist_func,
    SURF_NORMAL_FUNC surface_normal_func, TRANSFORMATION_FUNC transform_func, int normalize) {
        RAY_IN ray_in = { x, y, 0.0, 0.0, 0.0, 1.0 };
        RAY_OUT ray_out = SliceRayTrace(ray_in, pslice, transfer_dist_func,
            surface_normal_func, transform_func, normalize);
        *ln = ray_out.ln;
        *mn = ray_out.mn;
        *nn = ray_out.nn;
    }

/* --------------------------------------------------------------------
** Conic solutions (cv != 0, surf_type == 0)
** --------------------------------------------------------------------
*/

// Converts angles in radians to measured from vertex to off-axis distances. Only
// works if cv is non-zero. If cv is close to zero, the plane solutions (below)
// should be used instead, for which no off-axis distances are necessary.
void Conic2DOffAxisDistance(double *x0, double *y0, double cv, double k, double alpha, double beta) {

    // Put a minimum on cv to prevent division by zero
    if (fabs(cv) < 1E-13) cv = 1E-13;

    if (alpha == 0.0) {
        *y0 = 0.0;
    } else {
        double tana = tan(alpha);
        double num = (k - 1.0) + sqrt(4.0 + tana * tana * (3.0 - k));
        double den = tana * tana + (1.0 + k);
        *y0 = (tana / (2.0 * cv)) * (num / den);
    }

    // same as above but for x0 (beta)
    if (beta == 0.0) {
        *x0 = 0.0;
    } else {
        double tanb = tan(beta);
        double num = (k - 1.0) + sqrt(4.0 + tanb * tanb * (3.0 - k));
        double den = tanb * tanb + (1.0 + k);
        *x0 = (tanb / (2.0 * cv)) * (num / den);
    }

    // Make direction of effective rotation consistent with plane. For a reflective
    // surface where alpha and/or beta apply a rotation rather than an OAD, the
    // angles are effectively doubled.
    if (cv <= 0) {*x0 *=-1;}
    else {*y0 *=-1;}
}

// Converts off-axis distances to angles in radians. Only works if cv is non-zero.
double Conic2DOffAxisAngle(double* alpha, double* beta, double cv, double k, double x0, double y0) {
    if (fabs(cv) < 1E-13) return 0;
    double sagx = cv * x0 * x0 / (1 + sqrt(1 - (1 + k) * cv * cv * x0 * x0));
    double denom = 1 / (2*cv) - sagx;
    *alpha = atan(x0 / denom);

    double sagy = cv * y0 * y0 / (1 + sqrt(1 - (1 + k) * cv * cv * y0 * y0));
    denom = 1 / (2*cv) - sagy;
    *beta = atan(y0 / denom);

    if (cv <= 0) {*beta *=-1;}
    else {*alpha *=-1;}

}

void Conic2DTransformation(double coords_out[3], const double coords_in[3], SLICE_PARAMS pslice,
    int direction, int translate) {
    double alpha = ConvertAngle(pslice.alpha) * M_PI / 180.0;
    double beta  = ConvertAngle(pslice.beta)  * M_PI / 180.0;
    double gamma = ConvertAngle(pslice.gamma) * M_PI / 180.0;
    double zp = pslice.zp;
    double syx = pslice.syx;
    double syz = pslice.syz;
    double u = pslice.u;

    double x0, y0;
    Conic2DOffAxisDistance(&x0, &y0, pslice.cv, pslice.k, alpha, beta); 
    double cosg = cos(gamma);
    double sing = sin(gamma);

    double T[4][4] = {
        { 1, 0, 0, -u  },
        { 0, 1, 0,  0  },
        { 0, 0, 1,  zp },
        { 0, 0, 0,  1  }
    };

    double Ry1[4][4] = {
        { 1, 0, 0,  syx },
        { 0, 1, 0,   0  },
        { 0, 0, 1,  syz },
        { 0, 0, 0,   1  }
    };
    double Ry2[4][4] = {
        {  cosg, 0,  sing, 0 },
        {   0  , 1,    0 , 0 },
        { -sing, 0,  cosg, 0 },
        {   0  , 0,    0 , 1 }
    };
    double Ry3[4][4] = {
        { 1, 0, 0, -syx },
        { 0, 1, 0,   0  },
        { 0, 0, 1, -syz },
        { 0, 0, 0,   1  }
    };
    double Ry_temp[4][4];
    double Ry[4][4];
    Mat4Mul(Ry_temp, Ry1, Ry2);
    Mat4Mul(Ry, Ry_temp, Ry3);

    double TOAD[4][4] = {
        { 1, 0, 0, -x0 },
        { 0, 1, 0, -y0 },
        { 0, 0, 1,  0  },
        { 0, 0, 0,  1  }
    };

    // Full transformation matrix
    double Atot_temp[4][4];
    double Atot[4][4];
    Mat4Mul(Atot_temp, T, Ry);
    Mat4Mul(Atot, Atot_temp, TOAD);

    if (direction == -1) {
        double Atot_inv[4][4];
        Mat4AffineInverse(Atot_inv, Atot);
        memcpy(Atot, Atot_inv, sizeof(Atot_inv));
    }

    double w = translate ? 1.0 : 0.0;
    Mat4VecMul(coords_out, Atot, coords_in, w);
    }

double Conic2DTransfer(double xt, double yt, double zt, double l, double m, double n, SLICE_PARAMS pslice) {
    double cv = pslice.cv;
    double k = pslice.k;

    double A = 1 + k*n*n;
    double B = xt*l + yt*m + zt*n*(1 + k) - n/cv;
    double C = xt*xt + yt*yt + zt*zt*(1 + k) - 2*zt/cv;

    double discrim = B*B - A*C;
    if (discrim < 0) {return NAN;}

    int sgn = (cv > 0) ? 1 : -1;
    return C / ( -B + sgn * sqrt(discrim) );
}

void Conic2DSurfaceNormal(double* ln, double* mn, double* nn, double x, double y, SLICE_PARAMS pslice, int normalize) {
    double cv = pslice.cv;
    double k = pslice.k;

    double discrim = 1 - cv*cv*(1+k)*(x*x + y*y);
    if (discrim < 0) {
        *ln = NAN;
        *mn = NAN;
        *nn = NAN;
        return;
    }

    double denom = sqrt(discrim);
    double dervx = cv * x / denom;
    double dervy = cv * y / denom;
    double dervz = -1;

    *ln = dervx;
    *mn = dervy;
    *nn = dervz;

    if (normalize) {
        double norm = sqrt(dervx*dervx + dervy*dervy + dervz*dervz);
        *ln /= norm;
        *mn /= norm;
        *nn /= norm;
    }
}

// Wrapper function to make accessing the partial derivative about x easier. This
// is only used for the critical point calculation. 
static double Conic2DDervX(double x, double y, SLICE_PARAMS pslice) {
    double dervx, junk1, junk2;
    SliceSurfaceNormal(&dervx, &junk1, &junk2, x, y, pslice,
        &Conic2DTransfer, &Conic2DSurfaceNormal, &Conic2DTransformation, 0);
    return dervx;
}

// The critical points are used to find the bounded extrema within a slice.
// Solving for dz/dx = 0 analytically is not feasible, so we can use the secant
// method to find the solution. This should be fast because the derivative is monotonic
// and well-behaved.
void Conic2DCriticalXY(double *xc, double *yc, SLICE_PARAMS pslice) {
    double tol = 1E-13;    // Tolerance for accepting root
    int niter_max = 50;    // Max number of iterations - under normal circumstances
                           // should converge quickly (10 iterations or less)

    double cv = pslice.cv;
    double k = pslice.k;
    double alpha = ConvertAngle(pslice.alpha) * M_PI / 180.0;
    double beta  = ConvertAngle(pslice.beta)  * M_PI / 180.0;
    double gamma = ConvertAngle(pslice.gamma) * M_PI / 180.0;

    // If the curvature is very small this is basically a plane. No critical points.
    if (fabs(cv) < 1E-13) {
        *yc = NAN; *xc = NAN;
        return;
    }
    
    double x0, y0;
    Conic2DOffAxisDistance(&x0, &y0, cv, k, alpha, beta);

    // Perform secant method to find xc; it should usually be close to x0.
    // Use x0 +/- 0.05 * (radius of curvature) as initial guesses.
    double xc0 = x0 + 0.05 * 1/cv;
    double xc1 = x0 - 0.05 * 1/cv;
    double xc2, dervx0, dervx1;
    double err = 1;
    int i = 1;
    while (err < tol && i < niter_max) {
        // Solution along y-direction is guaranteed to always be at y = -y0
        dervx0 = Conic2DDervX(xc0, -y0, pslice);
        dervx1 = Conic2DDervX(xc1, -y0, pslice);
        xc2 = (xc0 * dervx1 - xc1 * dervx0) / (dervx1 - dervx0);
        err = fabs(xc2 - xc1);
        xc0 = xc1; xc1 = xc2;
    }
    // If it was unable to converge for some reason set xc to NAN. The critical
    // point will be ignored in the calculation of the global extrema.
    *xc = (err < tol) ? xc2 : NAN;
    *yc = -y0;
}

/* --------------------------------------------------------------------
** Plane solutions (cv == 0)
** --------------------------------------------------------------------
*/

void PlaneTransformation(double coords_out[3], const double coords_in[3], SLICE_PARAMS pslice,
    int direction, int translate) {
    double alpha = ConvertAngle(pslice.alpha) * M_PI / 180.0;
    double beta  = ConvertAngle(pslice.beta)  * M_PI / 180.0;
    double gamma = ConvertAngle(pslice.gamma) * M_PI / 180.0;

    double zp  = pslice.zp;
    double syx = pslice.syx;
    double syz = pslice.syz;
    double sxy = pslice.sxy;
    double sxz = pslice.sxz;
    double u   = pslice.u;

    double cosa  = cos(alpha);
    double sina  = sin(alpha);
    double cosbg = cos(beta + gamma);
    double sinbg = sin(beta + gamma);

    double T[4][4] = {
        { 1, 0, 0, -u  },
        { 0, 1, 0,  0  },
        { 0, 0, 1,  zp },
        { 0, 0, 0,  1  }
    };
    double Ry1[4][4] = {
        { 1, 0, 0,  syx },
        { 0, 1, 0,   0  },
        { 0, 0, 1,  syz },
        { 0, 0, 0,   1  }
    };
    double Ry2[4][4] = {
        {  cosbg, 0,  sinbg, 0 },
        {   0   , 1,    0  , 0 },
        { -sinbg, 0,  cosbg, 0 },
        {   0   , 0,    0  , 1 }
    };
    double Ry3[4][4] = {
        { 1, 0, 0, -syx },
        { 0, 1, 0,   0  },
        { 0, 0, 1, -syz },
        { 0, 0, 0,   1  }
    };
    double Ry_temp[4][4];
    double Ry[4][4];
    Mat4Mul(Ry_temp, Ry1, Ry2);
    Mat4Mul(Ry, Ry_temp, Ry3);

    double Rx1[4][4] = {
        { 1, 0, 0,  0   },
        { 0, 1, 0,  sxy },
        { 0, 0, 1,  sxz },
        { 0, 0, 0,  1   }
    };
    double Rx2[4][4] = {
        { 1,   0 ,    0 , 0 },
        { 0,  cosa, -sina, 0 },
        { 0,  sina,  cosa, 0 },
        { 0,   0 ,    0 , 1 }
    };
    double Rx3[4][4] = {
        { 1, 0, 0,   0   },
        { 0, 1, 0, -sxy  },
        { 0, 0, 1, -sxz  },
        { 0, 0, 0,   1   }
    };
    double Rx_temp[4][4];
    double Rx[4][4];
    Mat4Mul(Rx_temp, Rx1, Rx2);
    Mat4Mul(Rx, Rx_temp, Rx3);

    double Atot_temp[4][4];
    double Atot[4][4];
    Mat4Mul(Atot_temp, T, Ry);
    Mat4Mul(Atot, Atot_temp, Rx);

    if (direction == -1) {
        double Atot_inv[4][4];
        Mat4AffineInverse(Atot_inv, Atot);
        memcpy(Atot, Atot_inv, sizeof(Atot_inv));
    }

    double w = translate ? 1.0 : 0.0;
    Mat4VecMul(coords_out, Atot, coords_in, w);
}


double PlaneTransfer(double xt, double yt, double zt, double l, double m, double n, SLICE_PARAMS pslice) {
    if (n == 0) { return NAN; }
    return -zt / n;
}

void PlaneSurfaceNormal(double *ln, double *mn, double *nn, double x, double y, SLICE_PARAMS pslice, int normalize) {
    *ln = 0;
    *mn = 0;
    *nn = -1.0;
}

void PlaneCriticalXY(double *xc, double *yc, SLICE_PARAMS pslice) {
     *xc = NAN; *yc = NAN;
}

/* --------------------------------------------------------------------
** Cylinder solutions (cv != 0, surf_type == 1)
** --------------------------------------------------------------------
*/

void CylinderTransformation(double coords_out[3], const double coords_in[3], SLICE_PARAMS pslice,
    int direction, int translate) {
    double alpha = ConvertAngle(pslice.alpha) * M_PI / 180.0;
    double beta  = ConvertAngle(pslice.beta)  * M_PI / 180.0;
    double gamma = ConvertAngle(pslice.gamma) * M_PI / 180.0;

    double zp  = pslice.zp;
    double syx = pslice.syx;
    double syz = pslice.syz;
    double sxy = pslice.sxy;
    double sxz = pslice.sxz;
    double u   = pslice.u;

    double x0, y0;
    Conic2DOffAxisDistance(&x0, &y0, pslice.cv, pslice.k, alpha, beta); 
    double cosg = cos(gamma);
    double sing = sin(gamma);
    double cosa  = cos(alpha);
    double sina  = sin(alpha);

    double T[4][4] = {
        { 1, 0, 0, -u  },
        { 0, 1, 0,  0  },
        { 0, 0, 1,  zp },
        { 0, 0, 0,  1  }
    };
    double Ry1[4][4] = {
        { 1, 0, 0,  syx },
        { 0, 1, 0,   0  },
        { 0, 0, 1,  syz },
        { 0, 0, 0,   1  }
    };
    double Ry2[4][4] = {
        {  cosg, 0,  sing, 0 },
        {   0  , 1,    0 , 0 },
        { -sing, 0,  cosg, 0 },
        {   0  , 0,    0 , 1 }
    };
    double Ry3[4][4] = {
        { 1, 0, 0, -syx },
        { 0, 1, 0,   0  },
        { 0, 0, 1, -syz },
        { 0, 0, 0,   1  }
    };
    double Ry_temp[4][4];
    double Ry[4][4];
    Mat4Mul(Ry_temp, Ry1, Ry2);
    Mat4Mul(Ry, Ry_temp, Ry3);

    double Rx1[4][4] = {
        { 1, 0, 0,  0   },
        { 0, 1, 0,  sxy },
        { 0, 0, 1,  sxz },
        { 0, 0, 0,  1   }
    };
    double Rx2[4][4] = {
        { 1,   0 ,    0 , 0 },
        { 0,  cosa, -sina, 0 },
        { 0,  sina,  cosa, 0 },
        { 0,   0 ,    0 , 1 }
    };
    double Rx3[4][4] = {
        { 1, 0, 0,   0   },
        { 0, 1, 0, -sxy  },
        { 0, 0, 1, -sxz  },
        { 0, 0, 0,   1   }
    };
    double Rx_temp[4][4];
    double Rx[4][4];
    Mat4Mul(Rx_temp, Rx1, Rx2);
    Mat4Mul(Rx, Rx_temp, Rx3);

    double TOAD[4][4] = {
        { 1, 0, 0, -x0 },
        { 0, 1, 0,  0 },
        { 0, 0, 1,  0  },
        { 0, 0, 0,  1  }
    };

    double Atot_temp1[4][4];
    double Atot_temp2[4][4];
    double Atot[4][4];
    Mat4Mul(Atot_temp1, T, Ry);
    Mat4Mul(Atot_temp2, Atot_temp1, Rx);
    Mat4Mul(Atot, Atot_temp2, TOAD);

    if (direction == -1) {
        double Atot_inv[4][4];
        Mat4AffineInverse(Atot_inv, Atot);
        memcpy(Atot, Atot_inv, sizeof(Atot_inv));
    }

    double w = translate ? 1.0 : 0.0;
    Mat4VecMul(coords_out, Atot, coords_in, w);
}

double CylinderTransfer(double xt, double yt, double zt, double l, double m, double n, SLICE_PARAMS pslice) {
    double cv = pslice.cv;
    double k = pslice.k;

    double A = 1 + k*n*n;
    double B = xt*l + zt*n*(1 + k) - n/cv;
    double C = xt*xt + zt*zt*(1 + k) - 2*zt/cv;

    double discrim = B*B - A*C;
    if (discrim < 0) {return NAN;}
    int sgn = (cv > 0) ? 1 : -1;

    return C / ( -B + sgn * sqrt(discrim) );
}

void CylinderSurfaceNormal(double* ln, double* mn, double* nn, double x, double y, SLICE_PARAMS pslice, int normalize) {
    double cv = pslice.cv;
    double k = pslice.k;

    double discrim = 1 - cv*cv*(1+k)*x*x;
    if (discrim < 0) {
        *ln = NAN;
        *mn = NAN;
        *nn = NAN;
        return;
    }
    double denom = sqrt(discrim);
    double dervx = cv * x / denom;
    double dervy = 0;
    double dervz = -1;

    *ln = dervx;
    *mn = dervy;
    *nn = dervz;

    if (normalize) {
        double norm = sqrt(dervx*dervx + dervy*dervy + dervz*dervz);
        *ln /= norm;
        *mn /= norm;
        *nn /= norm;
    }

}

void CylinderCriticalXY(double* xc, double* yc, SLICE_PARAMS pslice) {
    *xc = NAN; *yc = NAN;
}