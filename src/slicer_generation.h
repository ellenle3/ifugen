#ifndef SLICER_GENERATION_H
#define SLICER_GENERATION_H

/** @struct IMAGE_SLICER_PARAMS
*   @brief Parameters that define the image slicer.
*   @var IMAGE_SLICER_PARAMS::n_each
*       Number of slices per block. n_each > 0
*   @var IMAGE_SLICER_PARAMS::n_rows
*       Number of rows. n_rows > 0
*   @var IMAGE_SLICER_PARAMS::n_cols
*       Number of columns. n_cols > 0
*   @var IMAGE_SLICER_PARAMS::mode
*       Mode for switching angles between blocks.
*   @var IMAGE_SLICER_PARAMS::trace_walls
*       If 0 do not ray trace walls, if 1 attempt to ray trace walls.
*   @var IMAGE_SLICER_PARAMS::dalpha
*       Change in angle between rows in degrees.
*   @var IMAGE_SLICER_PARAMS::dbeta
*       Change in angle between columns in degrees.
*   @var IMAGE_SLICER_PARAMS::dgamma
*       Change in angle between slices in degrees.
*   @var IMAGE_SLICER_PARAMS::alpha_cen
*       If n_each is odd, angle of central row in degrees. If even, mean of angle
*       between two center-most rows.
*   @var IMAGE_SLICER_PARAMS::beta_cen
*       Angle of central column; analagous to alpha_cen for columns.
*   @var IMAGE_SLICER_PARAMS::gamma_cen
*       Angle of central slice; analagous to alpha_cen for slices in each block.
*   @var IMAGE_SLICER_PARAMS::dx
*       Length of each slice. dx > 0
*   @var IMAGE_SLICER_PARAMS::dy
*       Width of each slice. dy > 0
*   @var IMAGE_SLICER_PARAMS::gx_width
*       Width of gaps between slices along x-direction. gx_width >= 0
*   @var IMAGE_SLICER_PARAMS::gx_depth
*       Depth of gaps between slices along x-direction. gx_depth >= 0
*   @var IMAGE_SLICER_PARAMS::gy_width
*       Width of gaps between columns. gy_width >= 0
*   @var IMAGE_SLICER_PARAMS::gy_depth
*       Depth of gaps between columns. gy_depth >= 0
*   @var IMAGE_SLICER_PARAMS::cv
*       Curvature. cv > 0
*   @var IMAGE_SLICER_PARAMS::k
*       Conic constant.
*/
typedef struct {
    int n_each;
    int n_rows;
    int n_cols;      
    int mode;
    int trace_walls;
    int active_x;
    int active_y;
    double dalpha;
    double dbeta;
    double dgamma;
    double alpha_cen;
    double beta_cen;
    double gamma_cen;
    double dx;
    double dy;
    double gx_width;
    double gx_depth;
    double gy_width;
    double gy_depth;
    double cv;
    double k;
} IMAGE_SLICER_PARAMS;

/** Function pointers */
typedef double (*SAG_FUNC)(double, double, double, double, double, double, double);
typedef void (*CRITICAL_XY_FUNC)(double*, double*, double*, double, double, double, double, double);
typedef double (*TRANSFER_DIST_FUNC)(double, double, double, double, double, double, double, double, double, double);
typedef void (*SURF_NORMAL_FUNC)(double*, double*, double*, double, double, double, double, double, double, double);

int CheckSlicerParams(IMAGE_SLICER_PARAMS);

void GetSlicerSize(double *xsize, double *ysize, IMAGE_SLICER_PARAMS p);

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

double FindBoundedSliceExtremum(double x0, double y0, int mode, IMAGE_SLICER_PARAMS p, SAG_FUNC sag_func, CRITICAL_XY_FUNC critical_xy_func);

void FindSlicerGlobalExtrema(double *zmin, double *zmax, IMAGE_SLICER_PARAMS p, SAG_FUNC sag_func, CRITICAL_XY_FUNC critical_xy_func);

double TransferEquation(double t, double xt, double yt, double l, double m, double n, IMAGE_SLICER_PARAMS p, SAG_FUNC sag_func);

void RayTraceSlicer(double* xs, double* ys, double* zs, double* t, double* ln, double* mn, double* nn,
    double xt, double yt, double l, double m, double n, IMAGE_SLICER_PARAMS p,
    SAG_FUNC sag_func, CRITICAL_XY_FUNC critical_xy_func, TRANSFER_DIST_FUNC transfer_dist_func, SURF_NORMAL_FUNC surf_normal_func)

#endif