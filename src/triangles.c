#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "triangles.h"

int CalcNumTriangles(IMAGE_SLICER_PARAMS_BASIC p, int Nx, int Ny) {

    // For computing how many triangles we need
	int num_slices_total = 2;
	int num_gaps_x = 0, num_gaps_y = 0;
	int num_walls_x = 0, num_walls_y = 0;
	int num_triangles_surface = 2, num_triangles_sides = 2;

    num_slices_total = p.n_each * p.n_rows * p.n_cols;
    if (p.gx_width > 0) {
        num_gaps_x = p.n_cols - 1;
        num_walls_x = num_gaps_x * 2;
    }
    else { num_gaps_x = 0; num_walls_x = 0; }
    if (p.gy_width > 0) {
        num_gaps_y = p.n_each * p.n_rows - 1;
        num_walls_y = num_gaps_y * 2;
    }
    else { num_gaps_y = 0; num_walls_y = 0; }
    // Number of triangles for just the surface
    num_triangles_surface = 2*Nx*Ny*num_slices_total + 2*Nx*num_walls_x + 2*Ny*num_walls_y + 2*(num_gaps_x + num_gaps_y);
    // Need to do 5 more panels - 4 for the sides and one for the back
    // THIS IS WRONG BECAUSE OF HOW SIDE PANELS ARE MADE
    num_triangles_sides = 2 + 2 * (Nx*num_slices_total + num_gaps_x) + 2 * (Ny*num_slices_total + num_gaps_y);
    num_triangles_sides = 10000; // temp for testing
    return num_triangles_surface + num_triangles_sides;
}

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

int IsPointsEqual(POINT3D p1, POINT3D p2) {
    if (fabs(p1.x - p2.x) < 1e-13 &&
        fabs(p1.y - p2.y) < 1e-13 &&
        fabs(p1.z - p2.z) < 1e-13) {
            return 1;
    }
    return 0;
}

int IsValidTriangle(POINT3D p1, POINT3D p2, POINT3D p3) {
    /* Ensure no two points coincide */
    if (IsPointsEqual(p1, p2) || IsPointsEqual(p1, p3) || IsPointsEqual(p2, p3)) {
        return 0;
    }

    /* Compute vectors */
    double v1x = p2.x - p1.x;
    double v1y = p2.y - p1.y;
    double v1z = p2.z - p1.z;

    double v2x = p3.x - p1.x;
    double v2y = p3.y - p1.y;
    double v2z = p3.z - p1.z;

    /* Cross product (v1 × v2) */
    double cx = v1y * v2z - v1z * v2y;
    double cy = v1z * v2x - v1x * v2z;
    double cz = v1x * v2y - v1y * v2x;

    /* Magnitude of cross product (2 * area of triangle) */
    double cross_mag = sqrt(cx*cx + cy*cy + cz*cz);

    /* If magnitude is nearly zero, points are collinear */
    if (cross_mag < 1e-13) {
        return 0;
    }

    return 1;
}


void SetTriListForTriangle(double *tri_list, int *num_triangles, POINT3D p1, POINT3D p2,
    POINT3D p3, double code) {
    if (!IsValidTriangle(p1, p2, p3)) {
        return;
    }
    
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

int GetSliceGridIndex(int nc, int ns, IMAGE_SLICER_PARAMS_BASIC p, int Nx, int Ny, int nx, int ny) {
    return (nc * p.n_rows * p.n_each + ns) * ((Ny + 1) * (Nx + 1)) + ny * (Nx + 1) + nx;
}

void EvalSliceGrid(double slice_grid[], double x_grid[], double y_grid[], int Nx, int Ny, IMAGE_SLICER_PARAMS_BASIC p, double p_custom[]) {
    double xsize, ysize;
    GetSlicerSize(&xsize, &ysize, p);
    double xstart, ystart;
    double xpts[Nx + 1], ypts[Ny + 1];
    double z;
    int idx;

    SLICE_PARAMS pslice;
    TRANSFER_DIST_FUNC transfer_dist_func;
    SURF_NORMAL_FUNC surf_normal_func;
    CRITICAL_XY_FUNC critical_xy_func;
    TRANSFORMATION_FUNC transform_func;

    for (int nr = 0; nr < p.n_rows; nr++) {

        xstart = -xsize / 2 + GetUForRow(nr, p_custom);

        for (int nc = 0; nc < p.n_cols; nc++) {

            ystart = -ysize / 2 + nr * p.n_each * (p.dy + p.gy_width);
            
            for (int ns = nr * p.n_each; ns < (nr + 1) * p.n_each; ns++) {

                pslice = GetSliceParams(ns, nc, p_custom);
                GetSurfaceFuncs(&transfer_dist_func, &surf_normal_func, &critical_xy_func, &transform_func, pslice, p);

                linspace(xpts, xstart, xstart + p.dx, Nx + 1);
                linspace(ypts, ystart, ystart + p.dy, Ny + 1);

                for (int ny = 0; ny < Ny + 1; ny++) {
                    for (int nx = 0; nx < Nx + 1; nx++) {
                        z = SliceSag(xpts[nx], ypts[ny], pslice, transfer_dist_func, transform_func);
                        idx = GetSliceGridIndex(nc, ns, p, Nx, Ny, nx, ny);
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
    int Nx, int Ny, IMAGE_SLICER_PARAMS_BASIC p) {

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
                    
                    i1 = GetSliceGridIndex(nc, ns, p, Nx, Ny, nx, ny);
                    i2 = GetSliceGridIndex(nc, ns, p, Nx, Ny, nx + 1, ny);
                    i3 = GetSliceGridIndex(nc, ns, p, Nx, Ny, nx, ny + 1);
                    i4 = GetSliceGridIndex(nc, ns, p, Nx, Ny, nx + 1, ny + 1);

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
    int Nx, int Ny, IMAGE_SLICER_PARAMS_BASIC p) {

    FACET facet;
    POINT3D p1, p2, p3, p4, pc;
    LINE line1, line2;
    int i1, i2, i3, i4;
    double x, y1, y2, y3, y4;
    double z1, z2, z3, z4;
    double code_a, code_b;
    double yc, zc;

    for (int nc = 0; nc < p.n_cols - 1; nc++) {

        for (int ns = 0; ns < p.n_each * p.n_rows; ns++) {

            for (int ny = 0; ny < Ny; ny++) {

                i1 = GetSliceGridIndex(nc, ns, p, Nx, Ny, Nx, ny);
                i2 = GetSliceGridIndex(nc, ns, p, Nx, Ny, Nx, ny + 1);
                i3 = GetSliceGridIndex(nc + 1, ns, p, Nx, Ny, 0, ny);
                i4 = GetSliceGridIndex(nc + 1, ns, p, Nx, Ny, 0, ny + 1);

                y1 = y_grid[i1];
                y2 = y_grid[i2];
                x = x_grid[i1];
                z1 = slice_grid[i1];
                z2 = slice_grid[i2];
                if (p.gx_width <= 0) {
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

                if (p.gx_width <= 0) continue;

                // Do other side of gap
                x = x_grid[i3];
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
    int Nx, int Ny, IMAGE_SLICER_PARAMS_BASIC p) {

    FACET facet;
    POINT3D p1, p2, p3, p4, pc;
    LINE line1, line2;
    int i1, i2, i3, i4;
    double x1, x2, x3, x4, y;
    double z1, z2, z3, z4;
    double code_a, code_b;
    double xc, zc;

    for (int nc = 0; nc < p.n_cols; nc++) {

        for (int ns = 0; ns < p.n_each * p.n_rows - 1; ns++) {

            for (int nx = 0; nx < Nx; nx++) {

                i1 = GetSliceGridIndex(nc, ns, p, Nx, Ny, nx, Ny);
                i2 = GetSliceGridIndex(nc, ns, p, Nx, Ny, nx + 1, Ny);
                i3 = GetSliceGridIndex(nc, ns + 1, p, Nx, Ny, nx, 0);
                i4 = GetSliceGridIndex(nc, ns + 1, p, Nx, Ny, nx + 1, 0);

                x1 = x_grid[i1];
                x2 = x_grid[i2];
                y = y_grid[i1];
                z1 = slice_grid[i1];
                z2 = slice_grid[i2];
                if (p.gy_width <= 0) {  // Connect to next slice directly
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
                    if (nx == 0) { code_b -= 3; }           // left is visible
                    else if (nx == Nx - 1) { code_a -= 3; } // right is visible

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
                    if (nx == 0) { code_b -= 4; }           // left is visible
                    else if (nx == Nx - 1) { code_a -= 4; } // right is visible
                    
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

                if (p.gy_width <= 0) continue;

                // Do other side of gap
                y = y_grid[i3];
                x1 = x_grid[i3];  // x-values need to be updated to match sampling of slice across the gap
                x2 = x_grid[i4];
                z1 = slice_grid[i3];
                z2 = slice_grid[i4];
                x3 = x1; x4 = x2;
                z3 = z4 = p.gy_depth;
                CrossoverPoint(&xc, &zc, (LINE){x1, z1, x2, z2}, (LINE){x3, z3, x4, z4});

                // Exact copy of above but across the gap, should probably refactor
                if (!isnan(xc)) {
                    code_a = code_b = 000113.0; // exact 0, CSG 1, left and right invisible
                    if (nx == 0) { code_b -= 3; }           // left is visible
                    else if (nx == Nx - 1) { code_a -= 3; } // right is visible

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
                    if (nx == 0) { code_b -= 4; }           // left is visible
                    else if (nx == Nx - 1) { code_a -= 4; } // right is visible
                    
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
    int Nx, int Ny, IMAGE_SLICER_PARAMS_BASIC p) {

    if (p.gx_width <= 0) {
        return;
    }

    FACET facet;
    POINT3D p1, p2, p3, p4;
    int i1, i2, i3, i4;
    double x1, x2, y1, y2, y3, y4;
    double code_a, code_b;
    
    for (int nc = 0; nc < p.n_cols - 1; nc++) {

        for (int ns = 0; ns < p.n_each * p.n_rows; ns++) {

            for (int ny = 0; ny < Ny; ny++) {

                i1 = GetSliceGridIndex(nc, ns, p, Nx, Ny, Nx, ny);
                i2 = GetSliceGridIndex(nc, ns, p, Nx, Ny, Nx, ny + 1);
                i3 = GetSliceGridIndex(nc + 1, ns, p, Nx, Ny, 0, ny);
                i4 = GetSliceGridIndex(nc + 1, ns, p, Nx, Ny, 0, ny + 1);
                x1 = x_grid[i1];
                x2 = x_grid[i3];
                y1 = y_grid[i1];
                y2 = y_grid[i2];
                y3 = y_grid[i3];
                y4 = y_grid[i4];

                code_a = code_b = 000213.0;             // exact 0, CSG 2, left and right visible
                if (ny == 0) { code_a -= 1; }           // bottom is visible
                else if (ny == Ny - 1) { code_b -= 2; } // top is visible

                p1 = (POINT3D){x1, y1, p.gx_depth};
                p2 = (POINT3D){x2, y3, p.gx_depth};
                p3 = (POINT3D){x1, y2, p.gx_depth};
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
    int Nx, int Ny, IMAGE_SLICER_PARAMS_BASIC p) {

    if (p.gy_width <= 0) {
        return;
    }

    FACET facet;
    POINT3D p1, p2, p3, p4;
    int i1, i2, i3, i4;;
    double x1, x2, x3, x4, y1, y2;
    double code_a, code_b;

    for (int nc = 0; nc < p.n_cols; nc++) {

        for (int ns = 0; ns < p.n_each * p.n_rows - 1; ns++) {

            for (int nx = 0; nx < Nx; nx++) {

                i1 = GetSliceGridIndex(nc, ns, p, Nx, Ny, nx, Ny);
                i2 = GetSliceGridIndex(nc, ns, p, Nx, Ny, nx + 1, Ny);
                i3 = GetSliceGridIndex(nc, ns + 1, p, Nx, Ny, nx, 0);
                i4 = GetSliceGridIndex(nc, ns + 1, p, Nx, Ny, nx + 1, 0);
                x1 = x_grid[i1];
                x2 = x_grid[i2];
                x3 = x_grid[i3];
                x4 = x_grid[i4];
                y1 = y_grid[i1];
                y2 = y_grid[i3];

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
    int Nx, int Ny, IMAGE_SLICER_PARAMS_BASIC p) {

    if (p.gx_width <= 0 && p.gy_width <= 0) {
        return;
    }
    // if only one is nonzero then still have to make walls

}

void MakeShellTriangles(double *tri_list, int *num_triangles, double slice_grid[], int Nx, int Ny, IMAGE_SLICER_PARAMS_BASIC p) {
    return;
}

void MakeAllTrianglesForSlicer(double *tri_list, int *num_triangles, int Nx, int Ny, IMAGE_SLICER_PARAMS_BASIC p, double p_custom[]) {

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
    MakeSliceTriangles(tri_list, num_triangles, slice_grid, x_grid, y_grid, Nx, Ny, p);
    MakeXWallTriangles(tri_list, num_triangles, slice_grid, x_grid, y_grid, Nx, Ny, p);
    MakeYWallTriangles(tri_list, num_triangles, slice_grid, x_grid, y_grid, Nx, Ny, p);
    MakeXGapTriangles(tri_list, num_triangles, x_grid, y_grid, Nx, Ny, p);
    MakeYGapTriangles(tri_list, num_triangles, x_grid, y_grid, Nx, Ny, p);
    //MakeShellTriangles(tri_list, num_triangles, slice_grid, Nx, Ny, p);

    free(slice_grid);
    free(x_grid);
    free(y_grid);
}

void TrianglesToSTL(double *tri_list, int num_triangles, const char *filename);