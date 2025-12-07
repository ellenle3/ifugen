#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "triangles.h"

int CalcNumTriangles(IMAGE_SLICER_PARAMS_BASIC p, int Nx, int Ny) {
    // Likely a slight overestimate of the number of triangles required.

    // For computing how many triangles we need
	int num_gaps_x, num_gaps_y;
	int num_walls_x, num_walls_y;
    int num_gaps_between = 0, num_walls_between = 0;

    // Number of triangles for the surface
    int num_slices_total = p.n_each * p.n_rows * p.n_cols;

    // Gap widths are always > 0 for now, even if they are negligibly small.
    num_gaps_x = p.n_cols * p.n_each * p.n_rows + 1;
    num_walls_x = num_gaps_x * 2;
    num_gaps_y = p.n_each * p.n_rows * p.n_cols + 1;
    num_walls_y = num_gaps_y * 2;
    num_gaps_between = (p.n_cols - 1) * (p.n_each * p.n_rows - 1);
    if (p.gx_depth != p.gy_depth) num_walls_between = num_gaps_between * 2;

    int num_facets_surface = Nx * Ny * num_slices_total  // Slice faces
        + Ny * num_walls_x + Nx * num_walls_y            // Walls
        + Ny * num_gaps_x
        + Nx * num_gaps_y
        + num_gaps_between + num_walls_between; 
    ;

    int zextra = (p.gy_depth == p.gx_depth) ? 2 : 3;
    int gy_extra = (p.gy_depth > p.gx_depth) ? 1 : 2;
    int gx_extra = (p.gx_depth > p.gy_depth) ? 1 : 2;
    
    // For the outer shell
    int num_facets_shell = 2 * (Ny * p.n_each * p.n_rows * zextra) // left and right slices
        + 2 * (Nx * p.n_cols * zextra)              // bottom and top slices
        + 2 * (num_gaps_y * gy_extra)   // left and right gaps
        + 2 * (num_gaps_x * gx_extra)   // bottom and top gaps
        + Nx * Ny * num_slices_total
        + Ny * num_gaps_x + Nx * num_gaps_y
        + num_gaps_between;
    ;

    int num_triangles_total = 2 * (num_facets_surface + num_facets_shell);

    // If we ever understimate the number of triangles needed, the DLL will crash.
    // Add another 5% just in case...

    return (int)(num_triangles_total * 1.05);
}

static void CrossoverPoint(double *wc, double* zc, LINE line1, LINE line2) {

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

static int IsPointsEqual(POINT3D p1, POINT3D p2) {
    if (fabs(p1.x - p2.x) < 1e-13 &&
        fabs(p1.y - p2.y) < 1e-13 &&
        fabs(p1.z - p2.z) < 1e-13) {
            return 1;
    }
    return 0;
}

// Used to create normal vector for STL file
static POINT3D CrossProduct(POINT3D u, POINT3D v) {
    POINT3D result;
    result.x = u.y * v.z - u.z * v.y;
    result.y = u.z * v.x - u.x * v.z;
    result.z = u.x * v.y - u.y * v.x;
    return result;
}

static POINT3D Normalize(POINT3D v) {
    double len = sqrt(v.x*v.x + v.y*v.y + v.z*v.z);
    if (len < 1e-13) return (POINT3D){0,0,0};
    return (POINT3D){v.x/len, v.y/len, v.z/len};
}

int IsValidTriangle(POINT3D p1, POINT3D p2, POINT3D p3) {
    /* Ensure no two points coincide */
    if (IsPointsEqual(p1, p2) || IsPointsEqual(p1, p3) || IsPointsEqual(p2, p3)) {
        return 0;
    }

    /* Compute edge vectors */
    POINT3D v1 = {p2.x - p1.x, p2.y - p1.y, p2.z - p1.z};
    POINT3D v2 = {p3.x - p1.x, p3.y - p1.y, p3.z - p1.z};

    /* Compute cross product */
    POINT3D cross = CrossProduct(v1, v2);

    /* Magnitude of cross product (2 * area of triangle) */
    double cross_mag = sqrt(cross.x*cross.x + cross.y*cross.y + cross.z*cross.z);

    /* If magnitude is nearly zero, points are collinear */
    if (cross_mag < 1e-13) {
        return 0;
    }
    return 1;
}

void CalcZBack(double* Z_back, double* Z_max_slice, double slice_grid[], int Nx, int Ny,
    double zdiff, IMAGE_SLICER_PARAMS_BASIC p) {

    double zmax = slice_grid[0];
    int Npts = (Nx + 1) * (Ny + 1) * p.n_each * p.n_rows * p.n_cols;
    for (int i = 1; i < Npts; i++) {
        if (slice_grid[i] > zmax) {
            zmax = slice_grid[i];
        }
    }
    *Z_max_slice = zmax;
    zmax = (p.gx_depth > 0 && zmax < p.gx_depth) ? p.gx_depth : zmax;
    zmax = (p.gy_depth > 0 && zmax < p.gy_depth) ? p.gy_depth : zmax;

    if (zdiff == 0) {
        zdiff = 2e-6; // Small extra margin
    }
    *Z_back = zmax + fabs(zdiff);
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
    double* xpts = (double*)malloc( (Nx + 1) * sizeof(double));
    double* ypts = (double*)malloc( (Ny + 1) * sizeof(double));
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
    free(xpts);
    free(ypts);
}

void MakeSliceTriangles(double *tri_list, int *num_triangles, double slice_grid[], double x_grid[], double y_grid[], 
    int Nx, int Ny, double Z_back, IMAGE_SLICER_PARAMS_BASIC p) {

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

                    // This is for the back of the shell
                    p1 = (POINT3D){x_grid[i1], y_grid[i1], Z_back};
                    p2 = (POINT3D){x_grid[i2], y_grid[i2], Z_back};
                    p3 = (POINT3D){x_grid[i3], y_grid[i3], Z_back};
                    p4 = (POINT3D){x_grid[i4], y_grid[i4], Z_back};

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
    int Nx, int Ny, double Z_back, IMAGE_SLICER_PARAMS_BASIC p) {

    if (p.gx_width <= 0) {
        return;
    }

    FACET facet;
    POINT3D p1, p2, p3, p4;
    int i1, i2, i3, i4;
    double x1, x2, x3, x4, y1, y2, y3, y4;
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


                // This is for the back of the shell
                code_a = 000310.0;                      // exact 0, CSG 1, top and bottom visible
                code_b = 000310.0;

                p1 = (POINT3D){x1, y1, Z_back};
                p2 = (POINT3D){x2, y3, Z_back};
                p3 = (POINT3D){x1, y2, Z_back};
                p4 = (POINT3D){x2, y4, Z_back};


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
    int Nx, int Ny, double Z_back, IMAGE_SLICER_PARAMS_BASIC p) {

    if (p.gy_width <= 0) {
        return;
    }

    FACET facet;
    POINT3D p1, p2, p3, p4;
    int i1, i2, i3, i4;;
    double x1, x2, x3, x4, y1, y2, y3, y4;
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

                // This is for the back of the shell
                code_a = 000310.0;
                code_b = 000310.0;

                p1 = (POINT3D){x1, y1, Z_back};
                p2 = (POINT3D){x2, y1, Z_back};
                p3 = (POINT3D){x3, y2, Z_back};
                p4 = (POINT3D){x4, y2, Z_back};

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
    int Nx, int Ny, double Z_back, IMAGE_SLICER_PARAMS_BASIC p) {

    if (p.gy_width <= 0 || p.gx_width <= 0) {
        return;
    }

    FACET facet;
    POINT3D p1, p2, p3, p4;
    int i1, i2, i3, i4;
    double x1, x2, x3, x4, y1, y2, y3, y4;
    double code_a, code_b;
    
    for (int nc = 0; nc < p.n_cols - 1; nc++) {

        for (int ns = 0; ns < p.n_each * p.n_rows - 1; ns++) {

            i1 = GetSliceGridIndex(nc, ns, p, Nx, Ny, Nx, Ny);
            i2 = GetSliceGridIndex(nc, ns + 1, p, Nx, Ny, Nx, 0);
            i3 = GetSliceGridIndex(nc + 1, ns, p, Nx, Ny, 0, Ny);
            i4 = GetSliceGridIndex(nc + 1, ns + 1, p, Nx, Ny, 0, 0);
            x1 = x_grid[i1];
            x2 = x_grid[i2];
            x3 = x_grid[i3];
            x4 = x_grid[i4];
            y1 = y_grid[i1];
            y2 = y_grid[i2];
            y3 = y_grid[i3];
            y4 = y_grid[i4];

            code_a = code_b = 000210.0;             // exact 0, CSG 2, left and right visible
            
            // First, fill in the x gap
            p1 = (POINT3D){x1, y1, p.gx_depth};
            p2 = (POINT3D){x2, y2, p.gx_depth};
            p3 = (POINT3D){x3, y3, p.gx_depth};
            p4 = (POINT3D){x4, y4, p.gx_depth};

            facet.p1 = p1;
            facet.p2 = p2;
            facet.p3 = p3;
            facet.p4 = p4;
            facet.code_a = code_a;
            facet.code_b = code_b;

            SetTriListForFacet(tri_list, num_triangles, facet);

            // Do the back of the shell as well
            p1 = (POINT3D){x1, y1, Z_back};
            p2 = (POINT3D){x2, y2, Z_back};
            p3 = (POINT3D){x3, y3, Z_back};
            p4 = (POINT3D){x4, y4, Z_back};

            facet.p1 = p1;
            facet.p2 = p2;
            facet.p3 = p3;
            facet.p4 = p4;
            facet.code_a = code_a;
            facet.code_b = code_b;

            SetTriListForFacet(tri_list, num_triangles, facet);

            // Now fill in the two walls on the sides

            if (p.gx_depth != p.gy_depth) {
                p1 = (POINT3D){x1, y1, p.gx_depth};
                p2 = (POINT3D){x2, y2, p.gx_depth};
                p3 = (POINT3D){x1, y1, p.gy_depth};
                p4 = (POINT3D){x2, y2, p.gy_depth};

                facet.p1 = p1;
                facet.p2 = p2;
                facet.p3 = p3;
                facet.p4 = p4;
                facet.code_a = code_a;
                facet.code_b = code_b;

                SetTriListForFacet(tri_list, num_triangles, facet);

                p1 = (POINT3D){x3, y3, p.gx_depth};
                p2 = (POINT3D){x4, y4, p.gx_depth};
                p3 = (POINT3D){x3, y3, p.gy_depth};
                p4 = (POINT3D){x4, y4, p.gy_depth};

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

void MakeShellSideTriangles(double *tri_list, int *num_triangles, double slice_grid[], double x_grid[], double y_grid[], 
    int Nx, int Ny, double Z_back, IMAGE_SLICER_PARAMS_BASIC p) {

    FACET facet;
    POINT3D p1, p2, p3, p4, pc;
    LINE line1, line2;
    int i1, i2, i3, i4;
    double x, y, x1, x2, y1, y2;
    double z1, z2, z3;
    double code_a, code_b;

    // May need to split into 2 or 3 panels depending on gap intersections.
    double zvals[3] = {Z_back, Z_back, Z_back};
    int num_z = 1;
    int isYDeeper = 1;  // deeper is larger in +z

    if (p.gx_width > 0) {
        zvals[num_z-1] = p.gx_depth;
        num_z++;
        if (p.gy_width > 0 && p.gx_depth != p.gy_depth) {
            zvals[num_z-1] = p.gy_depth;
            num_z++;

            // Z_back is always the largest value. Swap the first two if needed
            // to sort from smallest to largest.
            if (zvals[0] > zvals[1]) {
                isYDeeper = 0;
                double temp = zvals[0];
                zvals[0] = zvals[1];
                zvals[1] = temp;
            }
        }
    }
    else if (p.gy_width > 0) {
        zvals[num_z] = p.gy_depth;
        num_z++;
    }

    // Left side : nc = 0
    for (int ns = 0; ns < p.n_each * p.n_rows; ns++) {

        for (int ny = 0; ny < Ny; ny++) {

            i1 = GetSliceGridIndex(0, ns, p, Nx, Ny, 0, ny);
            i2 = GetSliceGridIndex(0, ns, p, Nx, Ny, 0, ny + 1);

            x = x_grid[i1];
            y1 = y_grid[i1];
            y2 = y_grid[i2];
            z1 = slice_grid[i1];
            z2 = slice_grid[i2];

            for (int nz = 0; nz < num_z; nz++) {

                z3 = zvals[nz];

                code_a = 000310.0;
                code_b = 000310.0;                      // exact 0, CSG 3, top and bottom visible
                //if (ny == 0) { code_a -= 4; }           // left is visible
                //else if (ny == Ny - 1) { code_b -= 4; } // right is visible

                p1 = (POINT3D){x, y1, z1};
                p2 = (POINT3D){x, y2, z2};
                p3 = (POINT3D){x, y1, z3};
                p4 = (POINT3D){x, y2, z3};

                facet.p1 = p1;
                facet.p2 = p2;
                facet.p3 = p3;
                facet.p4 = p4;
                facet.code_a = code_a;
                facet.code_b = code_b;

                SetTriListForFacet(tri_list, num_triangles, facet);

                z1 = z2 = z3; // For next iteration
            }
        }
    }

    // Right side : nc = p.n_cols - 1
    for (int ns = 0; ns < p.n_each * p.n_rows; ns++) {

        for (int ny = 0; ny < Ny; ny++) {

            i1 = GetSliceGridIndex(p.n_cols - 1, ns, p, Nx, Ny, Nx, ny);
            i2 = GetSliceGridIndex(p.n_cols - 1, ns, p, Nx, Ny, Nx, ny + 1);

            x = x_grid[i1];
            y1 = y_grid[i1];
            y2 = y_grid[i2];
            z1 = slice_grid[i1];
            z2 = slice_grid[i2];

            for (int nz = 0; nz < num_z; nz++) {

                z3 = zvals[nz];

                code_a = 000310.0;
                code_b = 000310.0;                      // exact 0, CSG 3, top and bottom visible
                //if (ny == 0) { code_a -= 4; }           // left is visible
                //else if (ny == Ny - 1) { code_b -= 4; } // right is visible

                p1 = (POINT3D){x, y1, z1};
                p2 = (POINT3D){x, y2, z2};
                p3 = (POINT3D){x, y1, z3};
                p4 = (POINT3D){x, y2, z3};

                facet.p1 = p1;
                facet.p2 = p2;
                facet.p3 = p3;
                facet.p4 = p4;
                facet.code_a = code_a;
                facet.code_b = code_b;

                SetTriListForFacet(tri_list, num_triangles, facet);

                z1 = z2 = z3; // For next iteration
            }
        }
    }

    // Bottom side : ns = 0
    for (int nc = 0; nc < p.n_cols; nc++) {

        for (int nx = 0; nx < Nx; nx++) {

            i1 = GetSliceGridIndex(nc, 0, p, Nx, Ny, nx, 0);
            i2 = GetSliceGridIndex(nc, 0, p, Nx, Ny, nx + 1, 0);

            y = y_grid[i1];
            x1 = x_grid[i1];
            x2 = x_grid[i2];
            z1 = slice_grid[i1];
            z2 = slice_grid[i2];

            for (int nz = 0; nz < num_z; nz++) {

                z3 = zvals[nz];

                code_a = 000310.0;
                code_b = 000310.0;                      // exact 0, CSG 3, top and bottom visible
                //if (ny == 0) { code_a -= 4; }           // left is visible
                //else if (ny == Ny - 1) { code_b -= 4; } // right is visible

                p1 = (POINT3D){x1, y, z1};
                p2 = (POINT3D){x2, y, z2};
                p3 = (POINT3D){x1, y, z3};
                p4 = (POINT3D){x2, y, z3};

                facet.p1 = p1;
                facet.p2 = p2;
                facet.p3 = p3;
                facet.p4 = p4;
                facet.code_a = code_a;
                facet.code_b = code_b;

                SetTriListForFacet(tri_list, num_triangles, facet);

                z1 = z2 = z3; // For next iteration
            }
        }
    }

    // Top side : ns = p.n_each * p.n_rows - 1
    for (int nc = 0; nc < p.n_cols; nc++) {

        for (int nx = 0; nx < Nx; nx++) {

            i1 = GetSliceGridIndex(nc, p.n_each * p.n_rows - 1, p, Nx, Ny, nx, Ny);
            i2 = GetSliceGridIndex(nc, p.n_each * p.n_rows - 1, p, Nx, Ny, nx + 1, Ny);

            y = y_grid[i1];
            x1 = x_grid[i1];
            x2 = x_grid[i2];
            z1 = slice_grid[i1];
            z2 = slice_grid[i2];

            for (int nz = 0; nz < num_z; nz++) {

                z3 = zvals[nz];

                code_a = 000310.0;
                code_b = 000310.0;                      // exact 0, CSG 3, top and bottom visible
                //if (ny == 0) { code_a -= 4; }           // left is visible
                //else if (ny == Ny - 1) { code_b -= 4; } // right is visible

                p1 = (POINT3D){x1, y, z1};
                p2 = (POINT3D){x2, y, z2};
                p3 = (POINT3D){x1, y, z3};
                p4 = (POINT3D){x2, y, z3};

                facet.p1 = p1;
                facet.p2 = p2;
                facet.p3 = p3;
                facet.p4 = p4;
                facet.code_a = code_a;
                facet.code_b = code_b;

                SetTriListForFacet(tri_list, num_triangles, facet);

                z1 = z2 = z3; // For next iteration
            }
        }
    }

    // Now do sides between gaps, if applicable
    double x3, x4, y3, y4;

    if (p.gy_width > 0) {   

        for (int ns = 0; ns < p.n_each * p.n_rows - 1; ns++) {

            // Left side
            i1 = GetSliceGridIndex(0, ns, p, Nx, Ny, 0, Ny);
            i2 = GetSliceGridIndex(0, ns + 1, p, Nx, Ny, 0, 0);
            // Right side
            i3 = GetSliceGridIndex(p.n_cols - 1, ns, p, Nx, Ny, Nx, Ny);
            i4 = GetSliceGridIndex(p.n_cols - 1, ns + 1, p, Nx, Ny, Nx, 0);

            x1 = x_grid[i1];
            x2 = x_grid[i2];
            x3 = x_grid[i3];
            x4 = x_grid[i4];
            y1 = y_grid[i1];
            y2 = y_grid[i2];
            y3 = y_grid[i3];
            y4 = y_grid[i4];
            z1 = p.gy_depth;

            for (int nz = 0; nz < num_z - 1; nz++) {
                z2 = zvals[nz + 1];
                if (isYDeeper) {
                    z2 = Z_back;
                    nz++; // Skip the next iteration
                }

                code_a = 000310.0;
                code_b = 000310.0;

                p1 = (POINT3D){x1, y1, z1};
                p2 = (POINT3D){x2, y2, z1};
                p3 = (POINT3D){x1, y1, z2};
                p4 = (POINT3D){x2, y2, z2};

                facet.p1 = p1;
                facet.p2 = p2;
                facet.p3 = p3;
                facet.p4 = p4;
                facet.code_a = code_a;
                facet.code_b = code_b;

                SetTriListForFacet(tri_list, num_triangles, facet);

                p1 = (POINT3D){x3, y3, z1};
                p2 = (POINT3D){x4, y4, z1};
                p3 = (POINT3D){x3, y3, z2};
                p4 = (POINT3D){x4, y4, z2};

                facet.p1 = p1;
                facet.p2 = p2;
                facet.p3 = p3;
                facet.p4 = p4;
                facet.code_a = code_a;
                facet.code_b = code_b;

                SetTriListForFacet(tri_list, num_triangles, facet);

                z1 = z2;
            }

        }
    }
    if (p.gx_width > 0) {
        
        for (int nc = 0; nc < p.n_cols - 1; nc++) {

            // Bottom side
            i1 = GetSliceGridIndex(nc, 0, p, Nx, Ny, Nx, 0);
            i2 = GetSliceGridIndex(nc + 1, 0, p, Nx, Ny, 0, 0);
            // Top side
            i3 = GetSliceGridIndex(nc, p.n_rows * p.n_each - 1, p, Nx, Ny, Nx, Ny);
            i4 = GetSliceGridIndex(nc + 1, p.n_rows * p.n_each - 1, p, Nx, Ny, 0, Ny);

            x1 = x_grid[i1];
            x2 = x_grid[i2];
            x3 = x_grid[i3];
            x4 = x_grid[i4];
            y1 = y_grid[i1];
            y2 = y_grid[i2];
            y3 = y_grid[i3];
            y4 = y_grid[i4];
            z1 = p.gx_depth;

            for (int nz = 0; nz < num_z - 1; nz++) {
                z2 = zvals[nz + 1];
                if (!isYDeeper) {
                    z2 = Z_back; 
                    nz++; // Skip the next iteration
                } 

                code_a = 000310.0;
                code_b = 000310.0;

                p1 = (POINT3D){x1, y1, z1};
                p2 = (POINT3D){x2, y2, z1};
                p3 = (POINT3D){x1, y1, z2};
                p4 = (POINT3D){x2, y2, z2};

                facet.p1 = p1;
                facet.p2 = p2;
                facet.p3 = p3;
                facet.p4 = p4;
                facet.code_a = code_a;
                facet.code_b = code_b;

                SetTriListForFacet(tri_list, num_triangles, facet);

                p1 = (POINT3D){x3, y3, z1};
                p2 = (POINT3D){x4, y4, z1};
                p3 = (POINT3D){x3, y3, z2};
                p4 = (POINT3D){x4, y4, z2};

                facet.p1 = p1;
                facet.p2 = p2;
                facet.p3 = p3;
                facet.p4 = p4;
                facet.code_a = code_a;
                facet.code_b = code_b;

                SetTriListForFacet(tri_list, num_triangles, facet);

                z1 = z2;
            }

        }
    }
}

void MakeAllTrianglesForSlicer(double *tri_list, int *num_triangles, int Nx, int Ny, double zdiff,
    IMAGE_SLICER_PARAMS_BASIC p, double p_custom[]) {

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

    double Z_back, Z_max_slice;
    CalcZBack(&Z_back, &Z_max_slice, slice_grid, Nx, Ny, zdiff, p);

    IMAGE_SLICER_PARAMS_BASIC p2 = p;
    if (p.gy_width == 0) {
        fprintf(stderr, "Error: Triangle generation currently requires non-zero gap widths. Setting gy_width to infinitesimal value.\n");
        p2.gy_width = 1e-9;
        p2.gy_depth = Z_max_slice + 1e-9;
    }
    if (p.gx_width == 0) {
        fprintf(stderr, "Error: Triangle generation currently requires non-zero gap widths. Setting gx_width to infinitesimal value.\n");
        p2.gx_width = 1e-9;
        p2.gx_depth = Z_max_slice + 1e-9;
    }

    MakeSliceTriangles(tri_list, num_triangles, slice_grid, x_grid, y_grid, Nx, Ny, Z_back, p2);
    MakeXWallTriangles(tri_list, num_triangles, slice_grid, x_grid, y_grid, Nx, Ny, p2);
    MakeYWallTriangles(tri_list, num_triangles, slice_grid, x_grid, y_grid, Nx, Ny, p2);
    MakeXGapTriangles(tri_list, num_triangles, x_grid, y_grid, Nx, Ny, Z_back, p2);
    MakeYGapTriangles(tri_list, num_triangles, x_grid, y_grid, Nx, Ny, Z_back, p2);
    MakeGapBetweenTriangles(tri_list, num_triangles, x_grid, y_grid, Nx, Ny, Z_back, p2);
    MakeShellSideTriangles(tri_list, num_triangles, slice_grid, x_grid, y_grid, Nx, Ny, Z_back, p2);

    free(slice_grid);
    free(x_grid);
    free(y_grid);
}

void MakeTrianglesRefractive(double* tri_list, int num_triangles) {
    // subtract all codes by ...
    return;
}

void TrianglesToSTL(double *tri_list, int num_triangles, const char *filename) {

    FILE *f = fopen(filename, "w");
    if (!f) {
        perror("Failed to open file");
        return;
    }

    fprintf(f, "solid model\n");

    for (int i = 0; i < num_triangles; i++) {
        POINT3D p1 = {tri_list[i*10+0], tri_list[i*10+1], tri_list[i*10+2]};
        POINT3D p2 = {tri_list[i*10+3], tri_list[i*10+4], tri_list[i*10+5]};
        POINT3D p3 = {tri_list[i*10+6], tri_list[i*10+7], tri_list[i*10+8]};

        // Compute normal
        POINT3D u = {p2.x - p1.x, p2.y - p1.y, p2.z - p1.z};
        POINT3D v = {p3.x - p1.x, p3.y - p1.y, p3.z - p1.z};
        POINT3D normal = Normalize(CrossProduct(u, v));

        fprintf(f, "  facet normal %.12f %.12f %.12f\n", normal.x, normal.y, normal.z);
        fprintf(f, "    outer loop\n");
        fprintf(f, "      vertex %.12f %.12f %.12f\n", p1.x, p1.y, p1.z);
        fprintf(f, "      vertex %.12f %.12f %.12f\n", p2.x, p2.y, p2.z);
        fprintf(f, "      vertex %.12f %.12f %.12f\n", p3.x, p3.y, p3.z);
        fprintf(f, "    endloop\n");
        fprintf(f, "  endfacet\n");
    }

    fprintf(f, "endsolid model\n");
    fclose(f);
}