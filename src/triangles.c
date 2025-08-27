#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "triangles.h"


int GetSliceGridIndex(int nc, int ns, int Nx, int Ny, int i, int j) {
    return (nc * n_slices + ns) * Nx * Ny + i * Ny + j;
}

void EvalSliceGrid(double slice_grid[], int Nx, int Ny, IMAGE_SLICER_PARAMS p, double p_custom[]) {
    int n_slices = p.n_each * p.n_rows;

    SLICE_PARAMS pslice = GetSliceParams(0, 0, p, p_custom);
    SAG_FUNC sag_func;
    TRANSFER_DIST_FUNC transfer_dist_func;
    SURF_NORMAL_FUNC surf_normal_func;
    CRITICAL_XY_FUNC critical_xy_func;
    GetSurfaceFuncs(&sag_func, &transfer_dist_func, &surf_normal_func, &critical_xy_func, pslice, p);

    int xsize; int ysize;
    GetSlicerSize(&xsize, &ysize, p);
    double xstart = -xsize / 2 + pslice.u;
    double ystart = -ysize / 2;
    double xpts[Nx]; double ypts[Ny];
    double z;
    int idx;

    for (int nc = 0; nc < p.n_cols; nc++) {

        for (int ns = 0; ns < n_slices; ns++) {

            linspace(xpts, xstart, xstart + p.dx, Nx);
            linspace(ypts, ystart, ystart + p.dy, Ny);

            for (int i = 0; i < Nx; i++) {
                // need to reset xstart and ystart for u row...
                for (int j = 0; j < Ny; j++) {
                    z = SAG_FUNC(xpts[i], ypts[j], ns, nc, p, p_custom);
                    idx = GetSliceIndex(nc, ns, Nx, Ny, i, j);
                    slice_grid[idx] = z;
                }
            }
            ystart += p.dy + p.gy_width;
        }
    xstart += p.dx + p.gx_width;
    }
}

void SetTriListForFacet(double *tri_list, int *num_triangles, FACET facet) {
    tri_list[(*num_triangles)*10 + 0] = facet.x1;
    tri_list[(*num_triangles)*10 + 1] = facet.y1;
    tri_list[(*num_triangles)*10 + 2] = facet.za;
    tri_list[(*num_triangles)*10 + 3] = facet.x2;
    tri_list[(*num_triangles)*10 + 4] = facet.y1;
    tri_list[(*num_triangles)*10 + 5] = facet.zb;
    tri_list[(*num_triangles)*10 + 6] = facet.x1;
    tri_list[(*num_triangles)*10 + 7] = facet.y2;
    tri_list[(*num_triangles)*10 + 8] = facet.zc;
    tri_list[(*num_triangles)*10 + 9] = facet.code1;
    (*num_triangles)++;

    tri_list[(*num_triangles)*10 + 0] = facet.x2;
    tri_list[(*num_triangles)*10 + 1] = facet.y1;
    tri_list[(*num_triangles)*10 + 2] = facet.zb;
    tri_list[(*num_triangles)*10 + 3] = facet.x1;
    tri_list[(*num_triangles)*10 + 4] = facet.y2;
    tri_list[(*num_triangles)*10 + 5] = facet.zc;
    tri_list[(*num_triangles)*10 + 6] = facet.x2;
    tri_list[(*num_triangles)*10 + 7] = facet.y2;
    tri_list[(*num_triangles)*10 + 8] = facet.zd;
    tri_list[(*num_triangles)*10 + 9] = facet.code2;
    (*num_triangles)++;
}

void MakeSliceTriangles(double *tri_list, int *num_triangles, double slice_grid[], int Nx, int Ny, IMAGE_SLICER_PARAMS p, double p_custom[]) {
    
}

void MakeXWallTriangles(double *tri_list, int *num_triangles, double slice_grid[], int Nx, int Ny, IMAGE_SLICER_PARAMS p, double p_custom[]) {
}

void MakeYWallTriangles(double *tri_list, int *num_triangles, double slice_grid[], int Nx, int Ny, IMAGE_SLICER_PARAMS p, double p_custom[]) {
}

void MakeXGapTriangles(double *tri_list, int *num_triangles, IMAGE_SLICER_PARAMS p, double p_custom[]);

void MakeYGapTriangles(double *tri_list, int *num_triangles, IMAGE_SLICER_PARAMS p, double p_custom[]);

void MakeShellTriangles(double *tri_list, int *num_triangles, double slice_grid[], int Nx, int Ny, IMAGE_SLICER_PARAMS p, double p_custom[]);

void MakeAllTrianglesForSlicer(double *tri_list, int *num_triangles, int Nx, int Ny, IMAGE_SLICER_PARAMS p, double p_custom[]) {
    double *slice_grid = (double *)malloc(Nx * Ny * p.n_each * p.n_rows * sizeof(double));
    if (slice_grid == NULL) {
        fprintf(stderr, "Failed to allocate memory for slice grid\n");
        return;
    }
    MakeSliceTriangles(tri_list, num_triangles, slice_grid, Nx, Ny, p, p_custom);
    MakeXWallTriangles(tri_list, num_triangles, slice_grid, Nx, Ny, p, p_custom);
    MakeYWallTriangles(tri_list, num_triangles, slice_grid, Nx, Ny, p, p_custom);
    MakeXGapTriangles(tri_list, num_triangles, p, p_custom);
    MakeYGapTriangles(tri_list, num_triangles, p, p_custom);
    MakeShellTriangles(tri_list, num_triangles, slice_grid, Nx, Ny, p, p_custom);
    free(slice_grid);
}

void TrianglesToSTL(double *tri_list, int num_triangles, const char *filename);