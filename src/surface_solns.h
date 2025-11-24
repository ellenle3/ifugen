
#ifndef SURFACE_SOLNS_H
#define SURFACE_SOLNS_H

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

/** Function pointers */                                                                
typedef double (*SAG_FUNC)(double, double, SLICE_PARAMS);
typedef double (*TRANSFER_DIST_FUNC)(double, double, double, double, double, double, SLICE_PARAMS);
typedef void (*SURF_NORMAL_FUNC)(double*, double*, double*, double, double, SLICE_PARAMS, int);
typedef void (*CRITICAL_XY_FUNC)(double*, double*, SLICE_PARAMS);
typedef void (*TRANSFORMATION_FUNC)(double[3], const double[3], SLICE_PARAMS, int, int);

/**
 * @brief A struct to store the input ray parameters.
 */
typedef struct {
    double xt;  // Starting x-coordinate of the ray
    double yt;  // Starting y-coordinate of the ray
    double zt;  // Starting z-coordinate of the ray, usually 0 by convention
    double l;   // Direction cosine along x
    double m;   // Direction cosine along y
    double n;   // Direction cosine along z
} RAY_IN;

/**
 * @brief A structure to store the output ray parameters.
 */
typedef struct {
    double xs;  // x-coordinate at the surface
    double ys;  // y-coordinate at the surface
    double zs;  // z-coordinate at the surface
    double t;   // Transfer distance
    double ln;  // Surface normal along x
    double mn;  // Surface normal along y
    double nn;  // Surface normal along z
} RAY_OUT;

/** 
* @brief Converts an angle to be between -180 and 180 degrees.
*   @param t Angle in degrees.
*   @return The converted angle.
**/
static double ConvertAngle(double t);

/**
 * @brief Multiplies two 4x4 matrices.
 */
static void Mat4Mul(double out[4][4], const double A[4][4], const double B[4][4]);

/**
 * @brief Multiplies a 4x4 matrix with a 4x1 vector.
 * 
 * @param w 1 to apply translation, 0 to ignore translation. The 4x4 matrix represents
 *         an affine transformation in 3D space.
 */
static void Mat4VecMul(double out[3], const double M[4][4], const double v[3], const double w);

/**
 * @brief Computes the inverse of a 4x4 matrix representing an affine transformation.
 */
static void Mat4AffineInverse(double Minv[4][4], const double M[4][4]);

/**
 * @brief Sets surface normal to all NANs for when we don't want to compute them.
 *        Same call signature as other surface normal functions.
 */
static void NoSurfaceNormal(double *ln, double *mn, double *nn, double x, double y, SLICE_PARAMS pslice, int normalize);

/**
 * @brief Converts a ray from global to local coordinates for a given surface.
 * 
 * @param ray_in Input ray parameters in global coordinates.
 * @param pslice Slice parameters.
 * @param transform_func Transformation to apply.
 * @return RAY_IN Ray parameters in the local coordinates of the surface.
 */
RAY_IN ConvertRayInToLocal(RAY_IN ray_in, SLICE_PARAMS pslice, TRANSFORMATION_FUNC transform_func);

RAY_OUT ConvertRayOutToGlobal(RAY_OUT ray_out, SLICE_PARAMS pslice, TRANSFORMATION_FUNC transform_func);

/**
 * @brief Computes the ray trace for the slice.
 * 
 * @param ray_in Input ray parameters in global coordinates.
 * @param transfer_dist_func Transfer distance function to use.
 * @param surface_normal_func Surface normal function to use.
 * @param transform_func Transformation to apply.
 * @return RAY_OUT 
 */
RAY_OUT SliceRayTrace(RAY_IN ray_in, SLICE_PARAMS pslice, TRANSFER_DIST_FUNC transfer_dist_func,
    SURF_NORMAL_FUNC surface_normal_func, TRANSFORMATION_FUNC transform_func, int normalize);

/**
 * @brief Computes the sag of a transformed surface which is a special case where
 * the ray is parallel to the z-axis.
 * 
 * @param x x-coordinate to evaluate.
 * @param y y-coordinate to evaluate.
 */
double SliceSag(double x, double y, SLICE_PARAMS pslice, TRANSFER_DIST_FUNC transfer_dist_func,
    TRANSFORMATION_FUNC transform_func);

/**
 * @brief Computes the surface normal of a transformed surface at (x, y).
 */
void SliceSurfaceNormal(double* ln, double* mn, double* nn, double x, double y, SLICE_PARAMS pslice,
    TRANSFER_DIST_FUNC transfer_dist_func, SURF_NORMAL_FUNC surface_normal_func,
    TRANSFORMATION_FUNC transform_func, int normalize);

// Conic solutions (cv != 0, surf_type == 0)

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
 * @brief Transforms coordinates for a conicoid.
 * 
 * @param coords_out Array to store transformed coordinates.
 * @param coords_in Array of input coordinates.
 * @param pslice Slice parameters.
 * @param direction Direction of the transformation. 1 for forward, -1 for inverse.
 * @param translate If 1, apply the translation component. If 0, ignore translation.
 */
void Conic2DTransformation(double coords_out[3], const double coords_in[3], SLICE_PARAMS pslice,
    int direction, int translate);

/**
 * @brief Computes the ray transfer distance for a rotated conic.
 * 
 * @param xt Starting x-value to transfer.
 * @param yt Starting y-value to transfer.
 * @param zt Starting z-value to transfer. Typically =0 in global coordinates.
 * @param l Direction cosine along the x-axis.
 * @param m Direction cosine along the y-axis.
 * @param n Direction cosine along the z-axis.
 * @return double Transfer distance for the given ray to the surface.
 */
double Conic2DTransfer(double xt, double yt, double zt, double l, double m, double n, SLICE_PARAMS pslice);

/**
 * @brief Computes the surface normal vectors for a rotated conic. This is the
 * gradient of the sag.
 * 
 * @param ln Pointer for x-component of the normal vector.
 * @param mn Pointer for y-component of the normal vector.
 * @param nn Pointer for z-component of the normal vector.
 * @param x x-value to evaluate.
 * @param y y-value to evaluate.
 * @param normalize If 1, normalizes the vector. If 0, does not normalize.
 */
void Conic2DSurfaceNormal(double *ln, double *mn, double *nn, double x, double y, SLICE_PARAMS pslice, int normalize);

/**
 * @brief Partial derivative about x for the rotated conic sag.
 */
static double Conic2DDervX(double x, double y, SLICE_PARAMS pslice);

/**
 * @brief Computes the critical points (dz/dx = 0, dz/dy = 0) for a rotated conic.
 * 
 * @param xc Pointer for the x-coordinate of the critical point.
 * @param yc Pointer for the y-coordinate of the critical point.
 */
void Conic2DCriticalXY(double *xc, double *yc, SLICE_PARAMS pslice);


// Plane solutions (cv = 0)

/**
 * @brief Transforms coordinates for a plane.
 */
void PlaneTransformation(double coords_out[3], const double coords_in[3], SLICE_PARAMS pslice,
    int direction, int translate);

/**
 * @brief Transforms coordinates for a plane.
 */
double PlaneTransfer(double xt, double yt, double zt, double l, double m, double n, SLICE_PARAMS pslice);

/**
 * @brief Computes the surface normal vector for a plane.
 */
void PlaneSurfaceNormal(double *ln, double *mn, double *nn, double x, double y, SLICE_PARAMS pslice, int normalize);

/**
 * @brief Sets xc and yc to NAN because planes do not have critical points.
 */
void PlaneCriticalXY(double *xc, double *yc, SLICE_PARAMS pslice);

// Cylinder solutions (cv != 0, surf_type == 1)

/**
 * @brief Transforms coordinates for a cylinder.
 */
void CylinderTransformation(double coords_out[3], const double coords_in[3], SLICE_PARAMS pslice,
    int direction, int translate);

/**
 * @brief Computes the transfer distance for a cylinder.
 */
double CylinderTransfer(double xt, double yt, double zt, double l, double m, double n, SLICE_PARAMS pslice);

/**
 * @brief Computes the surface normal vector for a cylinder.
 */
void CylinderSurfaceNormal(double* ln, double* mn, double* nn, double x, double y, SLICE_PARAMS pslice, int normalize);

/**
 * @brief Computes the critical points (dz/dx = 0, dz/dy = 0) for a cylinder.
 */
void CylinderCriticalXY(double *xc, double *yc, SLICE_PARAMS pslice);

#endif