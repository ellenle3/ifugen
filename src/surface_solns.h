
#ifndef SAG_RAY_SOLNS_H
#define SAG_RAY_SOLNS_H

/** @brief Converts an angle to be between -180 and 180 degrees.
*   @param t Angle in degrees.
*   @return The converted angle.
**/
double ConvertAngle(double t);


/**  @brief Computes a conic surface rotated about its apex at the origin.
* 
*     [1]  x' = x cos(t) + z sin(t)
*     [2]  y' = y + y0
*     [3]  z' = -x sin(t) + z cos(t)
*     [4]  g  = y'^2 / (r_y + sqrt(r_y^2 - (k_y+1) y'^2))
*     [5]  z' = x'^2 / (r_x + sqrt(r_x^2 - (k_x+1) x'^2)) + g 
* 
*   Substitute 1, 2, 3, and 4 into 5 to get a quadratic equation which can then
*   be solved for z(x,y).
* 
*   @param z Sag of the image slicer.
*   @param x x-value to evaluate. The x-axis is oriented along the slice length.
*   @param y y-value to evaluate. The y-axis is oriented along the slice width.
*   @param c Curvature.
*   @param k Conic constant.
*   @param alpha Off-axis angle about the vertex (how an OAP is defined). Describes
*                a rotation in the direction of the x-axis.
*   @param beta Angle of rotation about the global y-axis in degrees.
*   
*   @return Sag of the conic surface.
**/
void Conic3DOffAxisDistance(double *x0, double *y0, double cv, double alpha, double beta);

double Conic3DSag(double x, double y, double cv, double k, double alpha, double beta, double gamma);

double Conic3DTransfer(double xt, double yt, double l, double m, double n, double cv, double k, double alpha, double beta, double gamma);

void Conic3DSurfaceNormal(double *ln, double *mn, double *nn, double x, double y, double cv, double k, double alpha, double beta, double gamma, int normalize);

double Conic3DDervX(double x, double y, double cv, double k, double alpha, double beta, double gamma);

void Conic3DCriticalXY(double *xc, double *yc, double cv, double k, double alpha, double beta, double gamma);


double TiltedPlaneSag(double x, double y, double cv, double k, double alpha, double beta, double gamma);

void TiltedPlaneCriticalXY(double *xc1, double *xc2, double *yc, double cv, double k, double alpha, double beta, double gamma);

double TiltedPlaneTransfer(double xt, double yt, double l, double m, double n, double cv, double k, double alpha, double beta, double gamma);

void TiltedPlaneSurfaceNormal(double *ln, double *mn, double *nn, double x, double y, double cv, double k, double alpha, double beta, double gamma, int normalize);


#endif