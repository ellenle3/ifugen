
#ifndef SAG_RAY_SOLNS_H
#define SAG_RAY_SOLNS_H

/**
 * @brief A struct to store the parameters of a single slice.
 */
typedef struct {
    double alpha;   // Off-axis angle along y-axis OR rotation about x-axis
    double beta;    // Off-axis angle along x-axis OR rotation about y-axis
    double gamma;   // Rotation about y-axis (may be combined with beta)
    double cv;      // Curvature = 1/R, where R is the radius of curvature
    double k;       // Conic constant
    double zp;      // Shift along z-axis (piston)
    double syx;     // x-coordinate of axis of rotation about y
    double syz;     // z-coordinate of axis of rotation about y
    double sxy;     // y-coordinate of axis of rotation about x (if applicable)
    double sxz;     // z-coordinate of axis of rotation about x (if applicable)
    double u;       // Row offset along x-axis
} SLICE_PARAMS;


/** @brief Converts an angle to be between -180 and 180 degrees.
*   @param t Angle in degrees.
*   @return The converted angle.
**/
double ConvertAngle(double t);

/**
 * @brief Computes the derivative of the solution to a quadratic:
 *       C / -(B + sgn * sqrt(B*B - A*C))
 *
 * @param sgn Sign of the derivative.
 * @param A Coefficient of the quadratic term.
 * @param B Coefficient of the linear term.
 * @param C Constant term.
 * @param dA Derivative of A.
 * @param dB Derivative of B.
 * @param dC Derivative of C.
 * @return double
 */
double CalcQuadraticDerv(int sgn, double A, double B, double C, double dA, double dB,
    double dC);

// Conic solutions (cv != 0)

/**
 * @brief Computes the off-axis distance from the off-axis angles and curvature.
 * 
 * @param x0 Pointer for off-axis distance in x.
 * @param y0 Pointer for off-axis distance in y.
 * @param cv Curvature of the surface - equal to 1 / R, the radius of curvature.
 * @param alpha Off-axis angle on the y-axis in degrees.
 * @param beta Off-axis angle on the x-axis in degrees.
 */
void Conic2DOffAxisDistance(double *x0, double *y0, double cv, double k, double alpha, double beta);

/**
 * @brief Computes an axially symmetric conic rotated about the global y-axis
 * 
 * @param x x-value to evaluate.
 * @param y y-value to evaluate.
 * @param pslice Parameters that define this slice.
 * @return double Sag at the given point.
 */
double Conic2DSag(double x, double y, SLICE_PARAMS pslice);

/**
 * @brief Computes the ray transfer distance for a rotated conic.
 * 
 * @param xt Starting x-value to transfer.
 * @param yt Starting y-value to transfer.
 * @param l Direction cosine along the x-axis.
 * @param m Direction cosine along the y-axis.
 * @param n Direction cosine along the z-axis.
 * @param pslice Parameters that define this slice.
 * @return double Transfer distance for the given ray to the surface.
 */
double Conic2DTransfer(double xt, double yt, double l, double m, double n, SLICE_PARAMS pslice);

/**
 * @brief Computes the surface normal vectors for a rotated conic. This is the
 * gradient of the sag.
 * 
 * @param ln Pointer for x-component of the normal vector.
 * @param mn Pointer for y-component of the normal vector.
 * @param nn Pointer for z-component of the normal vector.
 * @param x x-value to evaluate.
 * @param y y-value to evaluate.
 * @param pslice Parameters that define this slice.
 * @param normalize If 1, normalizes the vector. If 0, does not normalize.
 */
void Conic2DSurfaceNormal(double *ln, double *mn, double *nn, double x, double y, SLICE_PARAMS pslice, int normalize);

/**
 * @brief Partial derivative about x for the rotated conic sag. This is a wrapper
 * function for Conic3DSurfaceNormal that makes it easier to access the x derivative
 * when computing the critical points.
 * 
 * @param x x-value to evaluate.
 * @param y y-value to evaluate.
 * @param pslice Parameters that define this slice.
 * @return double Partial derivative of the sag with respect to x.
 */
double Conic2DDervX(double x, double y, SLICE_PARAMS pslice);

/**
 * @brief Computes the critical points (dz/dx = 0, dz/dy = 0) for a rotated conic.
 * 
 * @param xc Pointer for the x-coordinate of the critical point.
 * @param yc Pointer for the y-coordinate of the critical point.
 * @param pslice Parameters that define this slice.
 */
void Conic2DCriticalXY(double *xc, double *yc, SLICE_PARAMS pslice);


// Tilted plane solutions (cv = 0)

/**
 * @brief Computes the sag of a tilted plane. Some parameters are not used but are
 * present to match the function signature of the conic sag function.
 * 
 * @param x x-value to evaluate.
 * @param y y-value to evaluate.
 * @param pslice Parameters that define this slice.
 * @return double Sag at the given point.
 */
double TiltedPlaneSag(double x, double y, SLICE_PARAMS pslice);

/**
 * @brief Computes the transfer distance for a tilted plane. Some parameters are
 * not used but are present to match the function signature of the conic transfer
 * distance function.
 * 
 * @param xt Starting x-value to transfer.
 * @param yt Starting y-value to transfer.
 * @param l Direction cosine along the x-axis.
 * @param m Direction cosine along the y-axis.
 * @param n Direction cosine along the z-axis.
 * @param pslice Parameters that define this slice.
 * @return double Transfer distance for the given ray to the surface.
 */
double TiltedPlaneTransfer(double xt, double yt, double l, double m, double n, SLICE_PARAMS pslice);

/**
 * @brief Computes the surface normal vectors for a tilted plane. Some parameters
 * are not used but are present to match the function signature of the conic surface
 * normal function.
 * 
 * @param ln Pointer for x-component of the normal vector.
 * @param mn Pointer for y-component of the normal vector.
 * @param nn Pointer for z-component of the normal vector.
 * @param x x-value to evaluate.
 * @param y y-value to evaluate.
 * @param pslice Parameters that define this slice.
 * @param normalize If 1, normalizes the vector. If 0, does not normalize.
 */
void TiltedPlaneSurfaceNormal(double *ln, double *mn, double *nn, double x, double y, SLICE_PARAMS pslice, int normalize);

/**
 * @brief Sets xc and yc to NAN. This function would normally compute the critical
 * values (dz/dx = 0, dz/dy = 0) for this surface type, but planes do not have
 * critical points so we essentially do nothing here.
 * 
 * @param xc Pointer for the x-coordinate of the critical point.
 * @param yc Pointer for the y-coordinate of the critical point.
 * @param pslice Parameters that define this slice.
 */
void TiltedPlaneCriticalXY(double *xc, double *yc, SLICE_PARAMS pslice);

#endif