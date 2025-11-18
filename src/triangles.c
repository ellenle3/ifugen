#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "triangles.h"


void CrossoverPoint(double *wc, double* zc, LINE line1, LINE line2) {

    double m1, b1, m2, b2;
    if (fabs(line1.w2 - line1.w1) < 1e-13) {
        // In case the line is vertical, can still find a solution
        m1 = 1e13;
        b1 = line1.z1 - m1 * line1.w1;
    }
    else {
        m1 = (line1.z2 - line1.z1) / (line1.w2 - line1.w1);
        b1 = line1.z1 - m1 * line1.w1;
    }
    if (fabs(line2.w2 - line2.w1) < 1e-13) {
        m2 = 1e13;
        b2 = line2.z1 - m2 * line2.w1;
    }
    else {
        m2 = (line2.z2 - line2.z1) / (line2.w2 - line2.w1);
        b2 = line2.z1 - m2 * line2.w1;
    }

    // Parallel lines
    if (fabs(m1 - m2) < 1e-13) {
        *wc = NAN;
        *zc = NAN;
        return;
    }

    *wc = (b2 - b1) / (m1 - m2);
    *zc = m1 * (*wc) + b1;

    // Check if the crossover point is within the line segments
    if ((*wc < fmin(line1.w1, line1.w2)) || (*wc > fmax(line1.w1, line1.w2)) ||
        (*wc < fmin(line2.w1, line2.w2)) || (*wc > fmax(line2.w1, line2.w2))) {
        *wc = NAN;
        *zc = NAN;
        return;
    }
}

void SetTriListForTriangle(double *tri_list, int *num_triangles, POINT3D p1, POINT3D p2,
    POINT3D p3, double code) {
    
    tri_list[(*num_triangles)*10 + 0] = p1.x;
    tri_list[(*num_triangles)*10 + 1] = p1.y;
    tri_list[(*num_triangles)*10 + 2] = p1.z;
    tri_list[(*num_triangles)*10 + 3] = p2.x;
    tri_list[(*num_triangles)*10 + 4] = p2.y;
    tri_list[(*num_triangles)*10 + 5] = p2.z;
    tri_list[(*num_triangles)*10 + 6] = p3.x;
    tri_list[(*num_triangles)*10 + 7] = p3.y;
    tri_list[(*num_triangles)*10 + 8] = p3.z;
    tri_list[(*num_triangles)*10 + 9] = code;
    (*num_triangles)++;
}

void SetTriListForFacet(double *tri_list, int *num_triangles, FACET facet) {

    // check if reflective or refractive and update codes accordingly
    // if refractive code_a -=10 codeb-=10

    SetTriListForTriangle(tri_list, num_triangles, facet.p1, facet.p2, facet.p3, facet.code_a);
    SetTriListForTriangle(tri_list, num_triangles, facet.p2, facet.p3, facet.p4, facet.code_b);
}

int GetSliceGridIndex(int nc, int ns, int n_each, int Nx, int Ny, int nx, int ny) {
    return (nc * n_each + ns) * ((Ny + 1) * (Nx + 1)) + ny * (Nx + 1) + nx;
}

void EvalSliceGrid(double slice_grid[], double x_grid[], double y_grid[], int Nx, int Ny, IMAGE_SLICER_PARAMS p, double p_custom[]) {
    double xsize, ysize;
    GetSlicerSize(&xsize, &ysize, p);
    double xstart, ystart;
    double xpts[Nx], ypts[Ny];
    double z;
    int idx;

    SLICE_PARAMS pslice;
    SAG_FUNC sag_func;
    TRANSFER_DIST_FUNC transfer_dist_func;
    SURF_NORMAL_FUNC surf_normal_func;
    CRITICAL_XY_FUNC critical_xy_func;

    for (int nr = 0; nr < p.n_rows; nr++) {

        xstart = -xsize / 2 + GetUForRow(nr, p, p_custom);

        for (int nc = 0; nc < p.n_cols; nc++) {

            ystart = -ysize / 2 + nr * p.n_each * (p.dy + p.gy_width);

            for (int ns = nr * p.n_each; ns < (nr + 1) * p.n_each; ns++) {

                pslice = GetSliceParams(ns, nc, p, p_custom);
                GetSurfaceFuncs(&sag_func, &transfer_dist_func, &surf_normal_func, &critical_xy_func, pslice, p);

                linspace(xpts, xstart, xstart + p.dx, Nx + 1);
                linspace(ypts, ystart, ystart + p.dy, Ny + 1);

                for (int ny = 0; ny < Ny + 1; ny++) {
                    for (int nx = 0; nx < Nx + 1; nx++) {
                        z = sag_func(xpts[nx], ypts[ny], pslice);
                        idx = GetSliceGridIndex(nc, ns, p.n_each, Nx, Ny, nx, ny);
                        slice_grid[idx] = z;
                        x_grid[idx] = xpts[nx];
                        y_grid[idx] = ypts[ny];
                    }
                }
                ystart += p.dy + p.gy_width;
            }
            xstart += p.dx + p.gx_width;
        }
    }
}

void MakeSliceTriangles(double *tri_list, int *num_triangles, double slice_grid[], double x_grid[], double y_grid[], 
    int Nx, int Ny, IMAGE_SLICER_PARAMS p, double p_custom[]) {

    FACET facet;
    POINT3D p1, p2, p3, p4;
    int i1, i2, i3, i4;
    double code_a, code_b;

    for (int nc = 0; nc < p.n_cols; nc++) {

        for (int ns = 0; ns < p.n_each * p.n_rows; ns++) {

            for (int ny = 0; ny < Ny; ny++) {

                for (int nx = 0; nx < Nx; nx++) {

                    code_a = code_b = 010017.0;             // exact 1, CSG 0, all invisible
                    if (nx == 0) { code_a -= 4; }           // side from point 3 to 1 is visible on triangle 1
                    else if (nx == Nx - 1) { code_b -= 4; } // 3 to 1 is visible on triangle 2
                    if (ny == 0) { code_a -= 1; }           // 1 to 2 is visible on triangle 1
                    else if (ny == Ny - 1) { code_b -= 2; } // 2 to 3 is visible on triangle 2
                    
                    i1 = GetSliceGridIndex(nc, ns, p.n_each, Nx, Ny, nx, ny);
                    i2 = GetSliceGridIndex(nc, ns, p.n_each, Nx, Ny, nx + 1, ny);
                    i3 = GetSliceGridIndex(nc, ns, p.n_each, Nx, Ny, nx, ny + 1);
                    i4 = GetSliceGridIndex(nc, ns, p.n_each, Nx, Ny, nx + 1, ny + 1);

                    p1 = (POINT3D){x_grid[i1], y_grid[i1], slice_grid[i1]};
                    p2 = (POINT3D){x_grid[i2], y_grid[i2], slice_grid[i2]};
                    p3 = (POINT3D){x_grid[i3], y_grid[i3], slice_grid[i3]};
                    p4 = (POINT3D){x_grid[i4], y_grid[i4], slice_grid[i4]};

                    facet.p1 = p1;
                    facet.p2 = p2;
                    facet.p3 = p3;
                    facet.p4 = p4;
                    facet.code_a = code_a;
                    facet.code_b = code_b;

                    SetTriListForFacet(tri_list, num_triangles, facet);
                }
            }
        }
    }
}

void MakeXWallTriangles(double *tri_list, int *num_triangles, double slice_grid[], double x_grid[], double y_grid[], 
    int Nx, int Ny, IMAGE_SLICER_PARAMS p, double p_custom[]) {

    FACET facet;
    POINT3D p1, p2, p3, p4, pc;
    LINE line1, line2;
    int i1, i2, i3, i4;
    double x, y1, y2, y3, y4;
    double z1, z2, z3, z4;
    double code_a, code_b;
    double yc, zc;

    for (int nc = 0; nc < p.n_cols - 1; nc++) {

        for (int ns = 0; ns < p.n_each; ns++) {

            for (int ny = 0; ny < Ny; ny++) {

                i1 = GetSliceGridIndex(nc, ns, p.n_each, Nx, Ny, Nx, 0);
                i2 = GetSliceGridIndex(nc, ns, p.n_each, Nx, Ny, Nx, ny + 1);
                i3 = GetSliceGridIndex(nc + 1, ns, p.n_each, Nx, Ny, 0, ny);
                i4 = GetSliceGridIndex(nc + 1, ns, p.n_each, Nx, Ny, 0, ny + 1);

                y1 = y_grid[i1];
                y2 = y_grid[i2];
                x = y_grid[i1];
                z1 = slice_grid[i1];
                z2 = slice_grid[i2];
                if (p.gx_depth == 0) {
                    z3 = slice_grid[i3];
                    z4 = slice_grid[i4];
                    y3 = y_grid[i3];
                    y4 = y_grid[i4];
                }
                else {
                    z3 = z4 = p.gx_depth;
                    y3 = y1; y4 = y2;
                }
                CrossoverPoint(&yc, &zc, (LINE){y1, z1, y2, z2}, (LINE){y3, z3, y4, z4});

                if (!isnan(yc)) {
                    code_a = code_b = 000113.0; // exact 0, CSG 1, left and right invisible
                    if (ny == 0) { code_a -= 3; }           // left is visible
                    else if (ny == Ny - 1) { code_b -= 3; } // right is visible

                    p1 = (POINT3D){x, y1, z1};
                    p2 = (POINT3D){x, y2, z2};
                    p3 = (POINT3D){x, y3, z3};
                    p4 = (POINT3D){x, y4, z4};
                    pc = (POINT3D){x, yc, zc};

                    SetTriListForTriangle(tri_list, num_triangles, p1, pc, p3, code_a);
                    SetTriListForTriangle(tri_list, num_triangles, p2, pc, p4, code_b);
                }
                else {
                    code_a = 000115.0;
                    code_b = 000113.0;                      // exact 0, CSG 1, top and bottom visible
                    if (ny == 0) { code_a -= 4; }           // left is visible
                    else if (ny == Ny - 1) { code_b -= 4; } // right is visible

                    p1 = (POINT3D){x, y1, z1};
                    p2 = (POINT3D){x, y2, z2};
                    p3 = (POINT3D){x, y3, z3};
                    p4 = (POINT3D){x, y4, z4};

                    facet.p1 = p1;
                    facet.p2 = p2;
                    facet.p3 = p3;
                    facet.p4 = p4;
                    facet.code_a = code_a;
                    facet.code_b = code_b;

                    SetTriListForFacet(tri_list, num_triangles, facet);
                }

                if (p.gx_width == 0) continue;

                // Do other side of gap
                x = y_grid[i3];
                y1 = y_grid[i3];
                y2 = y_grid[i4];
                z1 = slice_grid[i3];
                z2 = slice_grid[i4];
                y3 = y1; y4 = y2;
                z3 = z4 = p.gx_depth;
                CrossoverPoint(&yc, &zc, (LINE){y1, z1, y2, z2}, (LINE){y3, z3, y4, z4});

                // Exact copy of above but across the gap, should probably refactor
                if (!isnan(yc)) {
                    code_a = code_b = 000113.0; // exact 0, CSG 1, left and right invisible
                    if (ny == 0) { code_a -= 3; }           // left is visible
                    else if (ny == Ny - 1) { code_b -= 3; } // right is visible

                    p1 = (POINT3D){x, y1, z1};
                    p2 = (POINT3D){x, y2, z2};
                    p3 = (POINT3D){x, y3, z3};
                    p4 = (POINT3D){x, y4, z4};
                    pc = (POINT3D){x, yc, zc};

                    SetTriListForTriangle(tri_list, num_triangles, p1, pc, p3, code_a);
                    SetTriListForTriangle(tri_list, num_triangles, p2, pc, p4, code_b);
                }
                else {
                    code_a = 000115.0;
                    code_b = 000113.0;                      // exact 0, CSG 1, top and bottom visible
                    if (ny == 0) { code_a -= 4; }           // left is visible
                    else if (ny == Ny - 1) { code_b -= 4; } // right is visible

                    p1 = (POINT3D){x, y1, z1};
                    p2 = (POINT3D){x, y2, z2};
                    p3 = (POINT3D){x, y3, z3};
                    p4 = (POINT3D){x, y4, z4};

                    facet.p1 = p1;
                    facet.p2 = p2;
                    facet.p3 = p3;
                    facet.p4 = p4;
                    facet.code_a = code_a;
                    facet.code_b = code_b;

                    SetTriListForFacet(tri_list, num_triangles, facet);
                }
            }
        }
    }
}

void MakeYWallTriangles(double *tri_list, int *num_triangles, double slice_grid[], double x_grid[], double y_grid[], 
    int Nx, int Ny, IMAGE_SLICER_PARAMS p, double p_custom[]) {

    FACET facet;
    POINT3D p1, p2, p3, p4, pc;
    LINE line1, line2;
    int i1, i2, i3, i4;
    double x1, x2, x3, x4, y;
    double z1, z2, z3, z4;
    double code_a, code_b;
    double xc, zc;

    for (int nc = 0; nc < p.n_cols; nc++) {

        for (int ns = 0; ns < p.n_each - 1; ns++) {

            for (int nx = 0; nx < Nx; nx++) {

                i1 = GetSliceGridIndex(nc, ns, p.n_each, Nx, Ny, nx, Ny);
                i2 = GetSliceGridIndex(nc, ns, p.n_each, Nx, Ny, nx + 1, Ny);
                i3 = GetSliceGridIndex(nc, ns + 1, p.n_each, Nx, Ny, nx, 0);
                i4 = GetSliceGridIndex(nc, ns + 1, p.n_each, Nx, Ny, nx + 1, 0);

                x1 = x_grid[i1];
                x2 = x_grid[i2];
                y = y_grid[i1];
                z1 = slice_grid[i1];
                z2 = slice_grid[i2];
                if (p.gy_depth == 0) {  // Connect to next slice directly
                    z3 = slice_grid[i3];
                    z4 = slice_grid[i4];
                    x3 = x_grid[i3];
                    x4 = x_grid[i4];
                }
                else {
                    z3 = z4 = p.gy_depth;  // Connect to the bottom of the gap instead
                    x3 = x1; x4 = x2;
                }
                CrossoverPoint(&xc, &zc, (LINE){x1, z1, x2, z2}, (LINE){x3, z3, x4, z4});

                if (!isnan(xc)) {
                    code_a = code_b = 000113.0; // exact 0, CSG 1, left and right invisible
                    if (nx == 0) { code_a -= 3; }           // left is visible
                    else if (nx == Nx - 1) { code_b -= 3; } // right is visible

                    p1 = (POINT3D){x1, y, z1};
                    p2 = (POINT3D){x2, y, z2};
                    p3 = (POINT3D){x3, y, z3};
                    p4 = (POINT3D){x4, y, z4};
                    pc = (POINT3D){xc, y, zc};

                    SetTriListForTriangle(tri_list, num_triangles, p1, pc, p3, code_a);
                    SetTriListForTriangle(tri_list, num_triangles, p2, pc, p4, code_b);
                }
                else {
                    code_a = 000115.0;
                    code_b = 000113.0;                      // exact 0, CSG 1, top and bottom visible
                    if (nx == 0) { code_a -= 4; }           // left is visible
                    else if (nx == Nx - 1) { code_b -= 4; } // right is visible
                    
                    p1 = (POINT3D){x1, y, z1};
                    p2 = (POINT3D){x2, y, z2};
                    p3 = (POINT3D){x3, y, z3};
                    p4 = (POINT3D){x4, y, z4};

                    facet.p1 = p1;
                    facet.p2 = p2;
                    facet.p3 = p3;
                    facet.p4 = p4;
                    facet.code_a = code_a;
                    facet.code_b = code_b;

                    SetTriListForFacet(tri_list, num_triangles, facet);
                }

                if (p.gy_width == 0) continue;

                // Do other side of gap
                y = y_grid[i3];
                x1 = slice_grid[i3];  // x-values need to be updated to match sampling of slice across the gap
                x2 = slice_grid[i4];
                z1 = slice_grid[i3];
                z2 = slice_grid[i4];
                x3 = x1; x4 = x2;
                z3 = z4 = p.gy_depth;
                CrossoverPoint(&xc, &zc, (LINE){x1, z1, x2, z2}, (LINE){x3, z3, x4, z4});

                // Exact copy of above but across the gap, should probably refactor
                if (!isnan(xc)) {
                    code_a = code_b = 000113.0; // exact 0, CSG 1, left and right invisible
                    if (nx == 0) { code_a -= 3; }           // left is visible
                    else if (nx == Nx - 1) { code_b -= 3; } // right is visible

                    p1 = (POINT3D){x1, y, z1};
                    p2 = (POINT3D){x2, y, z2};
                    p3 = (POINT3D){x3, y, z3};
                    p4 = (POINT3D){x4, y, z4};
                    pc = (POINT3D){xc, y, zc};

                    SetTriListForTriangle(tri_list, num_triangles, p1, pc, p3, code_a);
                    SetTriListForTriangle(tri_list, num_triangles, p2, pc, p4, code_b);
                }
                else {
                    code_a = 000115.0;
                    code_b = 000113.0;                      // exact 0, CSG 1, top and bottom visible
                    if (nx == 0) { code_a -= 4; }           // left is visible
                    else if (nx == Nx - 1) { code_b -= 4; } // right is visible
                    
                    p1 = (POINT3D){x1, y, z1};
                    p2 = (POINT3D){x2, y, z2};
                    p3 = (POINT3D){x3, y, z3};
                    p4 = (POINT3D){x4, y, z4};

                    facet.p1 = p1;
                    facet.p2 = p2;
                    facet.p3 = p3;
                    facet.p4 = p4;
                    facet.code_a = code_a;
                    facet.code_b = code_b;

                    SetTriListForFacet(tri_list, num_triangles, facet);
                }
            }
        }
    }
}

void MakeXGapTriangles(double *tri_list, int *num_triangles, double x_grid[], double y_grid[],
    int Nx, int Ny, IMAGE_SLICER_PARAMS p, double p_custom[]) {

    if (p.gx_width <= 0) {
        return;
    }

    FACET facet;
    POINT3D p1, p2, p3, p4;
    int i1, i2, i3, i4;
    double x1, x2, y1, y2, y3, y4;
    double code_a, code_b;
    
    for (int nc = 0; nc < p.n_cols - 1; nc++) {

        for (int ns = 0; ns < p.n_each; ns++) {

            for (int ny = 0; ny < Ny; ny++) {

                i1 = GetSliceGridIndex(nc, ns, p.n_each, Nx, Ny, Nx, ny);
                i2 = GetSliceGridIndex(nc, ns, p.n_each, Nx, Ny, Nx, ny + 1);
                i3 = GetSliceGridIndex(nc, ns + 1, p.n_each, Nx, Ny, 0, ny);
                i4 = GetSliceGridIndex(nc, ns + 1, p.n_each, Nx, Ny, 0, ny + 1);
                x1 = x_grid[i1];
                x2 = x1 + p.gx_width;
                y1 = y_grid[i1];
                y2 = y_grid[i2];
                y3 = y_grid[i3];
                y4 = y_grid[i4];

                code_a = code_b = 000213.0;             // exact 0, CSG 2, left and right visible
                if (ny == 0) { code_a -= 1; }           // bottom is visible
                else if (ny == Ny - 1) { code_b -= 2; } // top is visible

                p1 = (POINT3D){x1, y1, p.gx_depth};
                p2 = (POINT3D){x2, y1, p.gx_depth};
                p3 = (POINT3D){x1, y3, p.gx_depth};
                p4 = (POINT3D){x2, y4, p.gx_depth};

                facet.p1 = p1;
                facet.p2 = p2;
                facet.p3 = p3;
                facet.p4 = p4;
                facet.code_a = code_a;
                facet.code_b = code_b;

                SetTriListForFacet(tri_list, num_triangles, facet);
            }
        }
    }
}

void MakeYGapTriangles(double *tri_list, int *num_triangles, double x_grid[], double y_grid[],
    int Nx, int Ny, IMAGE_SLICER_PARAMS p, double p_custom[]) {

    if (p.gy_width <= 0) {
        return;
    }

    FACET facet;
    POINT3D p1, p2, p3, p4;
    int i1, i2, i3, i4;;
    double x1, x2, x3, x4, y1, y2;
    double code_a, code_b;

    for (int nc = 0; nc < p.n_cols; nc++) {

        for (int ns = 0; ns < p.n_each - 1; ns++) {

            for (int nx = 0; nx < Nx; nx++) {

                i1 = GetSliceGridIndex(nc, ns, p.n_each, Nx, Ny, nx, Ny);
                i2 = GetSliceGridIndex(nc, ns, p.n_each, Nx, Ny, nx + 1, Ny);
                i3 = GetSliceGridIndex(nc, ns + 1, p.n_each, Nx, Ny, nx, 0);
                i4 = GetSliceGridIndex(nc, ns + 1, p.n_each, Nx, Ny, nx + 1, 0);
                x1 = x_grid[i1];
                x2 = x_grid[i2];
                x3 = x_grid[i3];
                x4 = x_grid[i4];
                y1 = y_grid[i1];
                y2 = y1 + p.gy_width;

                code_a = 000215.0;
                code_b = 000213.0;                      // exact 0, CSG 2, top and bottom visible
                if (nx == 0) { code_a -= 4; }           // left is visible
                else if (nx == Nx - 1) { code_b -= 4; } // right is visible

                p1 = (POINT3D){x1, y1, p.gy_depth};
                p2 = (POINT3D){x2, y1, p.gy_depth};
                p3 = (POINT3D){x3, y2, p.gy_depth};
                p4 = (POINT3D){x4, y2, p.gy_depth};

                facet.p1 = p1;
                facet.p2 = p2;
                facet.p3 = p3;
                facet.p4 = p4;
                facet.code_a = code_a;
                facet.code_b = code_b;

                SetTriListForFacet(tri_list, num_triangles, facet);
            }
        }
    }
}

void MakeGapBetweenTriangles(double *tri_list, int *num_triangles, double x_grid[], double y_grid[],
    int Nx, int Ny, IMAGE_SLICER_PARAMS p, double p_custom[]) {

    if (p.gx_width <= 0 && p.gy_width <= 0) {
        return;
    }
    // if only one is nonzero then still have to make walls

}

void MakeShellTriangles(double *tri_list, int *num_triangles, double slice_grid[], int Nx, int Ny, IMAGE_SLICER_PARAMS p, double p_custom[]) {

}

void MakeAllTrianglesForSlicer(double *tri_list, int *num_triangles, int Nx, int Ny, IMAGE_SLICER_PARAMS p, double p_custom[]) {

    int Npts = (Nx + 1) * (Ny + 1) * p.n_each * p.n_rows * p.n_cols;
    double *slice_grid = (double *)malloc(Npts * sizeof(double));
    double *x_grid = (double *)malloc(Npts * sizeof(double));
    double *y_grid = (double *)malloc(Npts * sizeof(double));
    if (slice_grid == NULL || x_grid == NULL || y_grid == NULL) {
        fprintf(stderr, "Failed to allocate memory for slice grids\n");
        free(slice_grid);
        free(x_grid);
        free(y_grid);
        return;
    }

    EvalSliceGrid(slice_grid, x_grid, y_grid, Nx, Ny, p, p_custom);
    MakeSliceTriangles(tri_list, num_triangles, slice_grid, x_grid, y_grid, Nx, Ny, p, p_custom);
    MakeXWallTriangles(tri_list, num_triangles, slice_grid, x_grid, y_grid, Nx, Ny, p, p_custom);
    MakeYWallTriangles(tri_list, num_triangles, slice_grid, x_grid, y_grid, Nx, Ny, p, p_custom);
    MakeXGapTriangles(tri_list, num_triangles, x_grid, y_grid, Nx, Ny, p, p_custom);
    MakeYGapTriangles(tri_list, num_triangles, x_grid, y_grid, Nx, Ny, p, p_custom);
    MakeShellTriangles(tri_list, num_triangles, slice_grid, Nx, Ny, p, p_custom);

    free(slice_grid);
    free(x_grid);
    free(y_grid);
}

void TrianglesToSTL(double *tri_list, int num_triangles, const char *filename);