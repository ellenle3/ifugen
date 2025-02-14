#ifndef SLICER_GENERATION_H
#define SLICER_GENERATION_H

/**
 * @brief A struct to store the parameters of the image slicer.
 */
typedef struct {
    int n_each; // Number of slices in a single row
    int n_rows; // Number of rows
    int n_cols; // Number of columns
    int mode;   // Mode for alternating angles between rows
    int trace_walls; // Flag to trace the walls of the image slicer
    int active_x;    // Flag for central slice in x for an even number of slices
    int active_y;    // Flag for central slice in y for an even number of slices
    double dalpha;   // Difference in off-axis y-angle between rows in degrees
    double dbeta;    // Difference in off-axis x-angle between columns in degrees
    double dgamma;   // Difference in rotation angle between slices in degrees
    double gamma_offset; // Offset for the rotation angle
    double alpha_cen;    // Central off-axis angle along y in degrees
    double beta_cen;     // Central off-axis angle along x in degrees
    double gamma_cen;    // Central rotation angle about global y-axis in degrees
    double dx;           // Slice width
    double dy;           // Slice height
    double gx_width;     // Gap width between columns
    double gx_depth;     // Gap depth betwen column
    double gy_width;     // Gap depth between slices along y
    double gy_depth;     // Gap depth between slices along y
    double cv; // Curvature - equal to 1/R where R is the ROC
    double k;  // Conic constant
} IMAGE_SLICER_PARAMS;

/**
 * @brief A struct to store the input ray parameters.
 */
typedef struct {
    double xt;  // Starting x-coordinate of the ray
    double yt;  // Starting y-coordinate of the ray
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

/** Function pointers */
typedef double (*SAG_FUNC)(double, double, double, double, double, double, double);
typedef double (*TRANSFER_DIST_FUNC)(double, double, double, double, double, double, double, double, double, double);
typedef void (*CRITICAL_XY_FUNC)(double*, double*, double, double, double, double, double);
typedef void (*SURF_NORMAL_FUNC)(double*, double*, double*, double, double, double, double, double, double, double, int);

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
 * @brief Validates image slicer parameters. If parameters are illegal, they
 * are modified to safe values.
 * 
 * @param p Pointer for image slicer parameters.
 * @return int 0 if the parameters are okay. 1 if the parameters were not okay
 *          and were modified.
 */
int ValidateSlicerParams(IMAGE_SLICER_PARAMS *p);

/**
 * @brief Checks whether the image slicer parameters are equal, member by member.
 * 
 * @param p1 First IMAGE_SLICER_PARAMS struct to compare.
 * @param p2 Second IMAGE_SLICER_PARAMS struct to compare.
 * @return int 1 if all members are equal, 0 otherwise.
 */
int IsParametersEqual(IMAGE_SLICER_PARAMS p1, IMAGE_SLICER_PARAMS p2);

/**
 * @brief Get the appropriate functions for the surface type indicated by p.
 * 
 * @param sag_func Pointer to function that computes the sag of a single slice.
 * @param transfer_dist_func Pointer to function that computes the ray transfer distance of a slice.
 * @param critical_xy_func Pointer to function that computes the critical point of a slice.
 * @param surf_normal_func Pointer to function that computes the surface normal of a slice.
 * @param p 
 */
void GetSurfaceFuncs(SAG_FUNC *sag_func, TRANSFER_DIST_FUNC *transfer_dist_func,
CRITICAL_XY_FUNC *critical_xy_func, SURF_NORMAL_FUNC *surf_normal_func, IMAGE_SLICER_PARAMS p);

/**
 * @brief Computes the size of the image slicer.
 * 
 * @param xsize Pointer for size of the image slicer along the x-axis.
 * @param ysize Pointer for size of the image slicer along the y-axis.
 * @param p Image slicer parameters.
 */
void GetSlicerSize(double *xsize, double *ysize, IMAGE_SLICER_PARAMS p);

/**
 * @brief Computes column, slice indices for a given set of x, y.
 * 
 * @param col_num Pointer for the column number.
 * @param slice_num Pointer for the slice number.
 * @param x x-coordinate.
 * @param y y-coordinate.
 * @param p Image slicer parameters.
 */
void GetSlicerIndex(int *col_num, int *slice_num, double x, double y, IMAGE_SLICER_PARAMS p);

void IsInsideSlicerGap(int *in_xgap, int *in_ygap, double x, double y, IMAGE_SLICER_PARAMS p);

void GetParaxialSliceIndex(int *col_num, int *slice_num, IMAGE_SLICER_PARAMS p);

/** @brief Computes the angles alpha and beta for the given slice number. Indexing
*   of the slices starts at 0 from the bottom (negative y-direction) of the image
*   slicer. Rows are indexed the same way.
* 
*   @param alpha Pointer to angle about x.
*   @param beta Pointer to angle about y.
*   @param slice_num See ImageSlicerSag for remaining parameters.
**/
void GetSliceAngles(double* alpha, double* beta, double* gamma, int slice_num, int col_num, IMAGE_SLICER_PARAMS p);

/** @brief Computes the sag of the image slicer.
* 
*   @param z Pointer to sag of the image slicer.
*   @param x x-value to evaluate. The x-axis is oriented along the slice length.
*   @param y y-value to evaluate. The y-axis is oriented along the slice width.
*   @param p Image slicer parameters.
* 
*   @return 0 if success.
**/
double ImageSlicerSag(double x, double y, IMAGE_SLICER_PARAMS p, SAG_FUNC sag_func);

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
double FindBoundedSliceExtremum(double x0, double y0, int mode, IMAGE_SLICER_PARAMS p, SAG_FUNC sag_func, CRITICAL_XY_FUNC critical_xy_func);

/**
 * @brief Computes global extrema for the image slicer.
 * 
 * @param zmin Pointer to store the global minimum.
 * @param zmax Pointer to store the global maximum.
 * @param p Image slicer parameters.
 * @param sag_func Function to compute the sag of a slice.
 * @param critical_xy_func Function to compute the critical point of a slice.
 */
void FindSlicerGlobalExtrema(double *zmin, double *zmax, IMAGE_SLICER_PARAMS p, SAG_FUNC sag_func, CRITICAL_XY_FUNC critical_xy_func);

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
double TransferEquation(double t, double xt, double yt, double l, double m, double n, IMAGE_SLICER_PARAMS p, SAG_FUNC sag_func);

/**
 * @brief Computes the transfered ray and surface normals for the image slicer.
 * 
 * @param ray_out Pointer for the output ray parameters.
 * @param ray_in Input ray parameters.
 * @param zmin Global minimum of the image slicer.
 * @param zmax Global maximum of the image slicer.
 * @param p Image slicer parameters.
 * @param sag_func Function to compute the sag of a slice.
 * @param transfer_dist_func Function to compute the transfer distance for a slice.
 * @param surf_normal_func Function to compute the surface normal of a slice.
 */
void RayTraceSlicer(RAY_OUT *ray_out, RAY_IN ray_in, double zmin, double zmax, IMAGE_SLICER_PARAMS p,
    SAG_FUNC sag_func, TRANSFER_DIST_FUNC transfer_dist_func, SURF_NORMAL_FUNC surf_normal_func);

#endif