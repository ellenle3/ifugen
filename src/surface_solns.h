
#ifndef SAG_RAY_SOLNS_H
#define SAG_RAY_SOLNS_H

/** @brief Converts an angle to be between -180 and 180 degrees.
*   @param t Angle in degrees.
*   @return The converted angle.
**/
double ConvertAngle(double t);


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
void Conic3DOffAxisDistance(double *x0, double *y0, double cv, double alpha, double beta);

/**
 * @brief Computes an axially symmetric conic rotated about the global y-axis
 * 
 * @param x x-value to evaluate.
 * @param y y-value to evaluate.
 * @param cv Curvature of the surface.
 * @param k Conic constant.
 * @param alpha Off-axis angle on the y-axis in degrees.
 * @param beta Off-axis angle on the x-axis in degrees.
 * @param gamma Angle of rotation about the y-axis in degrees.
 * @return double Sag at the given point.
 */
double Conic3DSag(double x, double y, double cv, double k, double alpha, double beta, double gamma);

/**
 * @brief Computes the ray transfer distance for a rotated conic.
 * 
 * @param xt Starting x-value to transfer.
 * @param yt Starting y-value to transfer.
 * @param l Direction cosine along the x-axis.
 * @param m Direction cosine along the y-axis.
 * @param n Direction cosine along the z-axis.
 * @param cv Curvature of the surface.
 * @param k Conic constant.
 * @param alpha Off-axis angle on the y-axis in degrees.
 * @param beta Off-axis angle on the x-axis in degrees.
 * @param gamma Angle of rotation about the y-axis in degrees.
 * @return double Transfer distance for the given ray to the surface.
 */
double Conic3DTransfer(double xt, double yt, double l, double m, double n, double cv, double k, double alpha, double beta, double gamma);

/**
 * @brief Computes the surface normal vectors for a rotated conic. This is the
 * gradient of the sag.
 * 
 * @param ln Pointer for x-component of the normal vector.
 * @param mn Pointer for y-component of the normal vector.
 * @param nn Pointer for z-component of the normal vector.
 * @param x x-value to evaluate.
 * @param y y-value to evaluate.
 * @param cv Curvature of the surface.
 * @param k Conic constant.
 * @param alpha Off-axis angle on the y-axis in degrees.
 * @param beta Off-axis angle on the x-axis in degrees.
 * @param gamma Angle of rotation about the y-axis in degrees.
 * @param normalize If 1, normalizes the vector. If 0, does not normalize.
 */
void Conic3DSurfaceNormal(double *ln, double *mn, double *nn, double x, double y, double cv, double k, double alpha, double beta, double gamma, int normalize);

/**
 * @brief Partial derivative about x for the rotated conic sag. This is a wrapper
 * function for Conic3DSurfaceNormal that makes it easier to access the x derivative
 * when computing the critical points.
 * 
 * @param x x-value to evaluate.
 * @param y y-value to evaluate.
 * @param cv Curvature of the surface.
 * @param k Conic constant.
 * @param alpha Off-axis angle on the y-axis in degrees.
 * @param beta Off-axis angle on the x-axis in degrees.
 * @param gamma Angle of rotation about the y-axis in degrees.
 * @return double Partial derivative of the sag with respect to x.
 */
double Conic3DDervX(double x, double y, double cv, double k, double alpha, double beta, double gamma);

/**
 * @brief Computes the critical points (dz/dx = 0, dz/dy = 0) for a rotated conic.
 * 
 * @param xc Pointer for the x-coordinate of the critical point.
 * @param yc Pointer for the y-coordinate of the critical point.
 * @param cv Curvature of the surface.
 * @param k Conic constant.
 * @param alpha Off-axis angle on the y-axis in degrees.
 * @param beta Off-axis angle on the x-axis in degrees.
 * @param gamma Angle of rotation about the y-axis in degrees.
 */
void Conic3DCriticalXY(double *xc, double *yc, double cv, double k, double alpha, double beta, double gamma);


// Tilted plane solutions (cv = 0)

/**
 * @brief Computes the sag of a tilted plane. Some parameters are not used but are
 * present to match the function signature of the conic sag function.
 * 
 * @param x x-value to evaluate.
 * @param y y-value to evaluate.
 * @param cv Not used.
 * @param k Not used.
 * @param alpha Angle about the global y-axis in degrees.
 * @param beta Angle about the global x-axis in degrees.
 * @param gamma Additional angle about the global x-axis in degrees.
 * @return double Sag at the given point.
 */
double TiltedPlaneSag(double x, double y, double cv, double k, double alpha, double beta, double gamma);

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
 * @param cv Not used
 * @param k Not used.
 * @param alpha Angle about the global y-axis in degrees.
 * @param beta Angle about the global x-axis in degrees.
 * @param gamma Additional angle about the global x-axis in degrees.
 * @return double Transfer distance for the given ray to the surface.
 */
double TiltedPlaneTransfer(double xt, double yt, double l, double m, double n, double cv, double k, double alpha, double beta, double gamma);

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
 * @param cv Not used.
 * @param k Not used.
 * @param alpha Angle about the global y-axis in degrees.
 * @param beta Angle about the global x-axis in degrees.
 * @param gamma Additional angle about the global x-axis in degrees.
 * @param normalize If 1, normalizes the vector. If 0, does not normalize.
 */
void TiltedPlaneSurfaceNormal(double *ln, double *mn, double *nn, double x, double y, double cv, double k, double alpha, double beta, double gamma, int normalize);

/**
 * @brief Sets xc and yc to NAN. This function would normally compute the critical
 * values (dz/dx = 0, dz/dy = 0) for this surface type, but planes do not have
 * critical points so we essentially do nothing here.
 * 
 * @param xc Pointer for the x-coordinate of the critical point.
 * @param yc Pointer for the y-coordinate of the critical point.
 * @param cv Not used.
 * @param k Not used.
 * @param alpha Not used.
 * @param beta Not used.
 * @param gamma Not used.
 */
void TiltedPlaneCriticalXY(double *xc, double *yc, double cv, double k, double alpha, double beta, double gamma);

#endif