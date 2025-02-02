#ifndef IFU_HELPERS_H
#define IFU_HELPERS_H

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

/** @struct SUBPUPIL_MIRROR_PARAMS
* 
* 
*/
typedef struct {
    float d_sp;
    double c_sp;
    double k_sp; 
} SUBPUPIL_MIRROR_PARAMS;


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
double Conic3DSag(double x, double y, double c, double k, double alpha, double beta, double gamma);


/** @brief Computes the angles alpha and beta for the given slice number. Indexing
*   of the slices starts at 0 from the bottom (negative y-direction) of the image
*   slicer. Rows are indexed the same way.
* 
*   @param alpha Pointer to angle about x.
*   @param beta Pointer to angle about y.
*   @param slice_num See ImageSlicerSag for remaining parameters.
**/
void SliceAngles (double* alpha, double* beta, double* gamma, int slice_num, int col_num, IMAGE_SLICER_PARAMS p);


/** @brief Computes the sag of the image slicer.
* 
*   @param z Pointer to sag of the image slicer.
*   @param x x-value to evaluate. The x-axis is oriented along the slice length.
*   @param y y-value to evaluate. The y-axis is oriented along the slice width.
*   @param p Image slicer parameters.
* 
*   @return 0 if success.
**/
int ImageSlicerSag(double *z, double x, double y, IMAGE_SLICER_PARAMS p);

#endif