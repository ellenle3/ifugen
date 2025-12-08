#include "surface_solns.h"
#include "slice_param_helpers.h"

#ifndef SLICER_GENERATION_H
#define SLICER_GENERATION_H

/**
 * @brief A struct to store the bounds of a ray. The bounds are computed by the
 * function GetRayBounds.
 */
typedef struct {
    int nc_min;
    int ns_min;
    int nc_max;
    int ns_max;
    int sgnc;
    int sgns;
    double xmin;
    double ymin;
    double xmax;
    double ymax;
} RAY_BOUNDS;

/**
 * @brief Creates a linearly spaced array. This is akin to the numpy linspace
 * function in Python.
 * 
 * @param array Pointer to the array to store the values.
 * @param start Start value, inclusive.
 * @param end End value, exclusive.
 * @param n Number of points.
 */
void linspace(double *array, double start, double end, int n);

/**
 * @brief Get the appropriate functions for the surface type indicated by p.
 * 
 * @param transfer_dist_func Pointer to function that computes the ray transfer distance of a slice.
 * @param critical_xy_func Pointer to function that computes the critical point of a slice.
 * @param surf_normal_func Pointer to function that computes the surface normal of a slice.
 * @param transform_func Pointer to function that applies the transformation matrix of a slice.
 * @param pslice Slice parameters.
 * @param p Image slicer parameters.
 */
void GetSurfaceFuncs(TRANSFER_DIST_FUNC *transfer_dist_func, SURF_NORMAL_FUNC *surf_normal_func,
    CRITICAL_XY_FUNC *critical_xy_func, TRANSFORMATION_FUNC *transform_func, SLICE_PARAMS pslice, GRID_PARAMS_BASIC p);

/**
 * @brief Computes the size of the image slicer.
 * 
 * @param xsize Pointer for size of the image slicer along the x-axis.
 * @param ysize Pointer for size of the image slicer along the y-axis.
 * @param p Image slicer parameters.
 */
void GetSlicerSize(double *xsize, double *ysize, GRID_PARAMS_BASIC p);

/**
 * @brief Computes column, slice indices for a given set of x, y.
 * 
 * @param col_num Pointer for the column number.
 * @param slice_num Pointer for the slice number.
 * @param x x-coordinate.
 * @param y y-coordinate.
 * @param p Image slicer parameters.
 */
void GetSlicerIndex(int *col_num, int *slice_num, double x, double y, GRID_PARAMS_BASIC p, double p_custom[]);

/**
 * @brief Checks whether a point (x, y) is inside a gap.
 * 
 * @param in_xgap Pointer to store whether the point is inside a gap along x.
 * @param in_ygap Pointer to store whether the point is inside a gap along y.
 * @param x x-coordinate to check.
 * @param y y-coordinate to check.
 * @param p Image slicer parameters.
 * @param p_custom Custom slice parameters if applicable.
 */
void IsInsideSlicerGap(int *in_xgap, int *in_ygap, double x, double y, GRID_PARAMS_BASIC p, double p_custom[]);

/**
 * @brief Computes the minimum and maximum u values for the image slicer.
 * 
 * @param umin Pointer to store the minimum u value.
 * @param umax Pointer to store the maximum u value.
 * @param p Image slicer parameters.
 * @param p_custom Custom slice parameters if applicable.
 */
void GetMinMaxU(double *umin, double *umax, GRID_PARAMS_BASIC p, double p_custom[]);

/** @brief Computes the sag of the image slicer.
* 
*   @param x x-coordinate to evaluate.
*   @param y y-coordinate to evaluate.
*   @param p Image slicer parameters.
*   @param p_custom Custom slice parameters if applicable.
* 
*   @return The sag of the image slicer, or NAN if out of bounds.
**/
double ImageSlicerSag(double x, double y, GRID_PARAMS_BASIC p, double p_custom[]);

/**
 * @brief Finds an extremum of a slice within the bounds that it is defined for
 * this image slicer. An initial guess is used to determine which slice to check.
 * 
 * @param x0 Initial guess for the x-coordinate.
 * @param y0 Initial guess for the y-coordinate.
 * @param mode 0 for maximum, 1 for minimum.
 * @param p Image slicer parameters.
 * @param sag_func Function to compute the sag of a slice.
 * @param critical_xy_func Function to compute the critical point of a slice.
 * @return double Bounded maximum or minimum for a slice.
 */
double FindBoundedSliceExtremum(double x0, double y0, int mode, GRID_PARAMS_BASIC p, double p_custom[]);

/**
 * @brief Computes global extrema for the image slicer.
 * 
 * @param zmin Pointer to store the global minimum.
 * @param zmax Pointer to store the global maximum.
 * @param p Image slicer parameters.
 * @param sag_func Function to compute the sag of a slice.
 * @param critical_xy_func Function to compute the critical point of a slice.
 */
void FindSlicerGlobalExtrema(double *zmin, double *zmax, GRID_PARAMS_BASIC p, double p_custom[]);

/**
 * @brief Transfer equation for the entire image slicer. The roots of this function
 * give the transfer distance t.
 * 
 * @param t Transfer distance.
 * @param xt Starting x-coordinate of the ray.
 * @param yt Starting y-coordinate of the ray.
 * @param l Direction cosine along x.
 * @param m Direction cosine along y.
 * @param n Direction cosine along z.
 * @param p Image slicer parameters.
 * @param sag_func Function to compute the sag of a slice.
 * @return double The image slicer sag minus the computed z-coordinate at the surface (zs).
 */
double TransferFunction(double t, RAY_IN ray_in, GRID_PARAMS_BASIC p, double p_custom[]);

/**
 * @brief Computes the bounds of a ray.
 * 
 * @param ray_in Input ray parameters.
 * @param umin Minimum u value for the ray.
 * @param umax Maximum u value for the ray.
 * @param zmin Global minimum of the image slicer.
 * @param zmax Global maximum of the image slicer.
 * @param p Image slicer parameters.
 * @param p_custom Custom slice parameters if applicable.
 * @return RAY_BOUNDS The bounds of the ray.
 */
RAY_BOUNDS GetRayBounds(RAY_IN ray_in, double umin, double umax, double zmin, double zmax,
    GRID_PARAMS_BASIC p, void *p_custom);

int IsRayInBounds(int nc_min, int ns_min, int nc_max, int ns_max, double umax, double umin, GRID_PARAMS_BASIC p);

int IsSectionInBounds(int col_num, int row_num, double umin, double umax, GRID_PARAMS_BASIC p);

void CheckSliceSolution(RAY_OUT *ray_out, double tol, RAY_IN ray_in, int ns_test, int nc_test,
                           GRID_PARAMS_BASIC p, void *p_custom);

void CheckYWallCollision(RAY_OUT *ray_out, RAY_IN ray_in, int ns_test, int nc_test, int sgns,
    GRID_PARAMS_BASIC p, void *p_custom);

void CheckXWallCollision(RAY_OUT *ray_out, RAY_IN ray_in, int ns_test, int nc_test, int sgnc,
    GRID_PARAMS_BASIC p, void *p_custom);

void CalcNextCoords(double *x_next, double *y_next, int *code, RAY_IN ray_in, int sgnc, int sgns, int nc_test, int nr_test,
    double x_test, double y_test, double xmax, double ymax, GRID_PARAMS_BASIC p,
    void *p_custom);

void CalcNumSlicesToCheck(int sgnc, int sgns, int nc_test, int nr_test,
                           double x_test, double y_test,
                           double x_next, double y_next, int code,
                           GRID_PARAMS_BASIC p, void *p_custom,
                           int *n_stocheck, int *nc_new, int *nr_new);

int IsLastSliceInSection(int ns_test, int sgns, GRID_PARAMS_BASIC p);
            
/**
 * @brief Computes the transfered ray and surface normals for the image slicer.
 * 
 * @param ray_out Pointer for the output ray parameters.
 * @param ray_in Input ray parameters.
 * @param zmin Global minimum of the image slicer.
 * @param zmax Global maximum of the image slicer.
 * @param p Image slicer parameters.
 * @param trace_walls If 1, attempts to ray trace walls. If 0, walls are ignored.
 * @param sag_func Function to compute the sag of a slice.
 * @param transfer_dist_func Function to compute the transfer distance for a slice.
 * @param surf_normal_func Function to compute the surface normal of a slice.
 */
void RayTraceSlicer(RAY_OUT *ray_out, RAY_IN ray_in, double zmin, double zmax, double umin, double umax,
     int trace_walls, GRID_PARAMS_BASIC p, double p_custom[]);

void ParaxialRayTraceSlicer(RAY_OUT* ray_out, double* l_out, double* m_out, double* n_out,
    RAY_IN* ray_in, double n1, double n2, int active_x, int active_y,
    GRID_PARAMS_BASIC p, double p_custom[]);

#endif