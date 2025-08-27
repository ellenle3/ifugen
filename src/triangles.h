#include "slicer_generation.h"
#include "surface_solns.h"

#ifndef TRIANGLES_H
#define TRIANGLES_H

/**
 * @brief Represents a facet.
 * 
 */
typedef struct {
    double x1;
    double x2;
    double y1;
    double y2;
    double za; // z-coordinate at (x1, y1)
    double zb; // at (x2, y1)
    double zc; // at (x1, y2)
    double zd; // at (x2, y2)
    double code1;
    double code2;
} FACET;

int GetSliceGridIndex(int nc, int ns, int Nx, int Ny, int i, int j);

/**
 * @brief Evaluates the sag of every slice in the image slicer on a grid.
 * 
 * @param grid Array containing the sag.
 * @param Nx Number of points to sample along the x-axis.
 * @param Ny Number of points to sample along the y-axis.
 * @param p Image slicer parameters.
 * @param p_custom Custom slice parameters.
 */
void EvalSliceGrid(double slice_grid[], int Nx, int Ny, IMAGE_SLICER_PARAMS p, double p_custom[]);

void SetTriListForFacet(double *tri_list, int *num_triangles, FACET facet);

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
void MakeSliceTriangles(double *tri_list, int *num_triangles, double slice_grid[], int Nx, int Ny, IMAGE_SLICER_PARAMS p, double p_custom[]);

/**
 * @brief Generates triangles for walls along the x-direction.
 */
void MakeXWallTriangles(double *tri_list, int *num_triangles, double slice_grid[], int Nx, int Ny, IMAGE_SLICER_PARAMS p, double p_custom[]);

/**
 * @brief Generates triangles for walls along the y-direction.
 */
void MakeYWallTriangles(double *tri_list, int *num_triangles, double slice_grid[], int Nx, int Ny, IMAGE_SLICER_PARAMS p, double p_custom[]);

/**
 * @brief Generates triangles for gaps along the x-direction.
 */
void MakeXGapTriangles(double *tri_list, int *num_triangles, IMAGE_SLICER_PARAMS p, double p_custom[]);

/**
 * @brief Generates triangles for gaps along the y-direction.
 */
void MakeYGapTriangles(double *tri_list, int *num_triangles, IMAGE_SLICER_PARAMS p, double p_custom[]);

/**
 * @brief Generates triangles for the outer shell.
 */
void MakeShellTriangles(double *tri_list, int *num_triangles, double slice_grid[], int Nx, int Ny, IMAGE_SLICER_PARAMS p, double p_custom[]);

/**
 * @brief Generates all triangles for the image slicer.
 */
void MakeAllTrianglesForSlicer(double *tri_list, int *num_triangles, int Nx, int Ny, IMAGE_SLICER_PARAMS p, double p_custom[]);

/**
 * @brief Converts the triangle list to an STL file.
 * 
 * @param tri_list Triangle list array for a Zemax user-defined object.
 * @param num_triangles Total number of triangles.
 * @param filename Name of the output STL file.
 */
void TrianglesToSTL(double *tri_list, int num_triangles, const char *filename);

#endif