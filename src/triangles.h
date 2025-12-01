#include "slicer_generation.h"
#include "surface_solns.h"

#ifndef TRIANGLES_H
#define TRIANGLES_H

typedef struct {
    double x;
    double y;
    double z;
} POINT3D;

/**
 * @brief Represents a rectangular facet composed of two triangles. The first
 * triangle is defined by points 1, 2, and 3; the second triangle is defined
 * by points 2, 3, and 4.
 */
typedef struct {
    POINT3D p1;
    POINT3D p2;
    POINT3D p3;
    POINT3D p4;
    double code_a;  // See Zemax example for how this 6-digit code is defined
    double code_b;
} FACET;

/**
 * @brief Line defined by two points (w1, z1) and (w2, z2).
 */
typedef struct {
    double w1;
    double z1;
    double w2;
    double z2;
} LINE;

int CalcNumTriangles(IMAGE_SLICER_PARAMS_BASIC p, int Nx, int Ny);

static void CrossoverPoint(double *wc, double* zc, LINE line1, LINE line2);

static int IsPointsEqual(POINT3D p1, POINT3D p2);

static POINT3D CrossProduct(POINT3D u, POINT3D v);

static POINT3D Normalize(POINT3D v);

int IsValidTriangle(POINT3D p1, POINT3D p2, POINT3D p3);

void SetTriListForFacet(double *tri_list, int *num_triangles, FACET facet);

int GetSliceGridIndex(int nc, int ns, IMAGE_SLICER_PARAMS_BASIC p, int Nx, int Ny, int i, int j);

/**
 * @brief Evaluates the sag of every slice in the image slicer on a grid.
 * 
 * @param grid Array containing the sag.
 * @param Nx Number of points to sample along the x-axis.
 * @param Ny Number of points to sample along the y-axis.
 * @param p Image slicer parameters.
 * @param p_custom Custom slice parameters.
 */
void EvalSliceGrid(double slice_grid[], double x_grid[], double y_grid[], int Nx, int Ny, IMAGE_SLICER_PARAMS_BASIC p, double p_custom[]);

/**
 * @brief Generates triangles for the slices.
 * 
 * Each slice is divided into a grid of rectangular facets. Each facet is composed of two triangles.
 * If cv = 0, the parameters Nx and Ny are ignored because a planar slice can always be defined by
 * a single facet.
 * 
 * @param tri_list Triangle list array for a Zemax user-defined object.
 * @param num_triangles Current number of triangles.
 * @param Nx Number of points to sample along the x-axis per slice.
 * @param Ny Number of points to sample along the y-axis per slice.
 * @param p Image slicer parameters.
 * @param p_custom Custom slice parameters.
 */
void MakeSliceTriangles(double *tri_list, int *num_triangles, double slice_grid[], double x_grid[], double y_grid[], 
    int Nx, int Ny, double Z_back, IMAGE_SLICER_PARAMS_BASIC p);

/**
 * @brief Generates triangles for walls along the x-direction.
 */
void MakeXWallTriangles(double *tri_list, int *num_triangles, double slice_grid[], double x_grid[], double y_grid[], 
    int Nx, int Ny, IMAGE_SLICER_PARAMS_BASIC p);

/**
 * @brief Generates triangles for walls along the y-direction.
 */
void MakeYWallTriangles(double *tri_list, int *num_triangles, double slice_grid[], double x_grid[], double y_grid[], 
    int Nx, int Ny, IMAGE_SLICER_PARAMS_BASIC p);

/**
 * @brief Generates triangles for gaps along the x-direction.
 */
void MakeXGapTriangles(double *tri_list, int *num_triangles, double x_grid[], double y_grid[],
    int Nx, int Ny, double Z_back, IMAGE_SLICER_PARAMS_BASIC p);

/**
 * @brief Generates triangles for gaps along the y-direction.
 */
void MakeYGapTriangles(double *tri_list, int *num_triangles, double x_grid[], double y_grid[],
    int Nx, int Ny, double Z_back, IMAGE_SLICER_PARAMS_BASIC p);

void MakeGapBetweenTriangles(double *tri_list, int *num_triangles, double x_grid[], double y_grid[],
    int Nx, int Ny, double Z_back, IMAGE_SLICER_PARAMS_BASIC p);

double CalcZBack(double *tri_list, int *num_triangles, double slice_grid[],
    int Nx, int Ny, double zdiff, IMAGE_SLICER_PARAMS_BASIC p);

/**
 * @brief Generates all triangles for the image slicer.
 */
void MakeAllTrianglesForSlicer(double *tri_list, int *num_triangles, int Nx, int Ny, double zdiff,
    IMAGE_SLICER_PARAMS_BASIC p, double p_custom[]);

/**
 * @brief Converts the triangle list to an STL file.
 * 
 * @param tri_list Triangle list array for a Zemax user-defined object.
 * @param num_triangles Total number of triangles.
 * @param filename Name of the output STL file.
 */
void TrianglesToSTL(double *tri_list, int num_triangles, const char *filename);

#endif