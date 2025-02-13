#include <stdio.h>
#include <stdlib.h>
#include "slicer_generation.h"
#include "surface_solns.h"

void TestImageSlicerSag(FILE* fptr, IMAGE_SLICER_PARAMS p, SAG_FUNC sag_func);

int main() {
    FILE* fptr = fopen("test_output.txt", "w+");

    if (fptr==NULL) {
        printf("The file could not be opened.");
        return 1;
    }
    
    // Set up the image slicer

    IMAGE_SLICER_PARAMS p = {
        .n_each = 5,
        .n_rows = 3,
        .n_cols = 2,
        .mode = 0,
        .trace_walls = 0,
        .active_x = 0,
        .active_y = 0,
        .dalpha = 8,
        .dbeta = -5,
        .dgamma = 1.5,
        .alpha_cen = -10,
        .beta_cen = 11,
        .gamma_cen = 0,
        .dx = 9,
        .dy = 0.8,
        .gx_width = 0,
        .gx_depth = 0,
        .gy_width = 0,
        .gy_depth = 0,
        .cv = 0.0,
        .k = 0
    };

    SAG_FUNC sag_func;
    CRITICAL_XY_FUNC critical_xy_func;
    TRANSFER_DIST_FUNC transfer_dist_func;
    SURF_NORMAL_FUNC surf_normal_func; 

    if (p.cv == 0) {
            sag_func = &TiltedPlaneSag;
            critical_xy_func = &TiltedPlaneCriticalXY;
            transfer_dist_func = &TiltedPlaneTransfer;
            surf_normal_func = &TiltedPlaneSurfaceNormal;
         }
         else {
            sag_func = &Conic3DSag;
            critical_xy_func = &Conic3DCriticalXY;
            transfer_dist_func = &Conic3DTransfer;
            surf_normal_func = &Conic3DSurfaceNormal;
         }
    
    TestImageSlicerSag(fptr, p, sag_func);

    fclose(fptr);
}


void TestImageSlicerSag(FILE* fptr, IMAGE_SLICER_PARAMS p, SAG_FUNC sag_func) {
    
    // Compute test function over an array of points and save to a TXT file
    double xsize, ysize;
    GetSlicerSize(&xsize, &ysize, p);

    int nx = 20; int ny = 20;
    double xpts[nx]; double ypts[ny];
    linspace(xpts, -xsize/2, xsize/2, nx);
    linspace(ypts, -ysize/2, ysize/2, ny);

    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            double output = ImageSlicerSag(xpts[i], ypts[j], p, sag_func);
            fprintf(fptr, "%.10f %.10f %.10f\n", xpts[i], ypts[j], output);
        }
    }
}