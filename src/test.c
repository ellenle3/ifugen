#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "slicer_generation.h"
#include "surface_solns.h"
#include "custom_slicer_helpers.h"

void TestImageSlicerSag(FILE* fptr, IMAGE_SLICER_PARAMS p, SAG_FUNC sag_func);
void TestGlobalExtrema(FILE* fptr, IMAGE_SLICER_PARAMS p, SAG_FUNC sag_func, CRITICAL_XY_FUNC critical_xy_func);
void TestTransferDistance(FILE* fptr, IMAGE_SLICER_PARAMS p, TRANSFER_DIST_FUNC transfer_dist_func);
void TestRayTrace(FILE* fptr, IMAGE_SLICER_PARAMS p, SAG_FUNC sag_func, TRANSFER_DIST_FUNC transfer_dist_func,
SURF_NORMAL_FUNC surf_normal_func, CRITICAL_XY_FUNC critical_xy_func);
void TestLoadCustomParams(FILE* fptr, int file_num);

int main() {
    FILE* fptr = fopen("test_output.txt", "w+");

    if (fptr==NULL) {
        printf("The file could not be opened.");
        return 1;
    }
    
    // Set up the image slicer

    IMAGE_SLICER_PARAMS p = {
        .custom = 0,
        .cylinder = 0,
        .n_each = 1,
        .n_rows = 1,
        .n_cols = 1,
        .angle_mode = 0,
        .dalpha = 0,
        .dbeta = 0,
        .dgamma = 0,
        .alpha_cen = 0,
        .beta_cen = 0,
        .gamma_cen = 0,
        .gamma_offset = 0,
        .dx = 10,
        .dy = 10,
        .gx_width = 0,
        .gx_depth = 0,
        .gy_width = 0,
        .gy_depth = 0,
        .cv = -0.01,
        .k = -1
    };

    double custom_slice_params[1] = {0};

    SAG_FUNC sag_func;
    CRITICAL_XY_FUNC critical_xy_func;
    TRANSFER_DIST_FUNC transfer_dist_func;
    SURF_NORMAL_FUNC surf_normal_func; 

    GetSurfaceFuncs(&sag_func, &transfer_dist_func, &surf_normal_func, &critical_xy_func, p);
    
    //TestImageSlicerSag(fptr, p, sag_func);
    //TestGlobalExtrema(fptr, p, sag_func, critical_xy_func);
    //TestRayTrace(fptr, p, sag_func, transfer_dist_func, surf_normal_func, critical_xy_func);
    //TestTransferDistance(fptr, p, transfer_dist_func);
    TestLoadCustomParams(fptr, 0);

    fclose(fptr);
}

void TestImageSlicerSag(FILE* fptr, IMAGE_SLICER_PARAMS p, SAG_FUNC sag_func) {
    
    // Compute test function over an array of points and save to a TXT file
    double xsize, ysize;
    double output;
    double custom_slice_params[1] = {0};
    GetSlicerSize(&xsize, &ysize, p);

    int nx = 20; int ny = 20;
    double xpts[nx]; double ypts[ny];
    linspace(xpts, -xsize/2, xsize/2, nx);
    linspace(ypts, -ysize/2, ysize/2, ny);

    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            output = ImageSlicerSag(xpts[i], ypts[j], p, custom_slice_params, sag_func);
            fprintf(fptr, "%.10f %.10f %.10f\n", xpts[i], ypts[j], output);
        }
    }
}

void TestGlobalExtrema(FILE* fptr, IMAGE_SLICER_PARAMS p, SAG_FUNC sag_func, CRITICAL_XY_FUNC critical_xy_func) {
    
    double zmin, zmax;
    double custom_slice_params[1] = {0};
    FindSlicerGlobalExtrema(&zmin, &zmax, p, custom_slice_params, sag_func, critical_xy_func);
    fprintf(fptr, "%.10f %.10f\n", zmin, zmax);

}

void TestTransferDistance(FILE* fptr, IMAGE_SLICER_PARAMS p, TRANSFER_DIST_FUNC transfer_dist_func) {
    int n = 5;
    double t;
    double xpts[n]; double ypts[n];
    double lpts[n]; double mpts[n]; double npts[n];
    linspace(xpts, -3, 3, n);
    linspace(ypts, -3, 3, n);
    linspace(lpts, -0.2, 0.2, n);
    linspace(mpts, -0.2, 0.2, n);
    linspace(npts, -0.1, 1.3, n);

    for (int i1 = 0; i1 < n; i1++) {
        for (int i2 = 0; i2 < n; i2++) {
            for (int i3 = 0; i3 < n; i3++) {
                for (int i4 = 0; i4 < n; i4++) {
                    for (int i5 = 0; i5 < n; i5++) {

                        t = transfer_dist_func(xpts[i1], ypts[i2], lpts[i3], mpts[i4], npts[i5], p.cv, p.k, p.alpha_cen, p.beta_cen, p.gamma_cen);

                        fprintf(fptr, "%.10f %.10f %.10f %.10f %.10f %.10f\n", xpts[i1], ypts[i2], lpts[i3], mpts[i4], npts[i5], t);
                    }
                }
            }
        }
    }
}


void TestRayTrace(FILE* fptr, IMAGE_SLICER_PARAMS p, SAG_FUNC sag_func, TRANSFER_DIST_FUNC transfer_dist_func,
SURF_NORMAL_FUNC surf_normal_func, CRITICAL_XY_FUNC critical_xy_func) {

    RAY_IN ray_in;
    RAY_OUT ray_out;

    double zmin, zmax, norm, l, m, n;
    double custom_slice_params[1] = {0};
    FindSlicerGlobalExtrema(&zmin, &zmax, p, custom_slice_params, sag_func, critical_xy_func);

    int num = 2;
    double xpts[num]; double ypts[num];
    double lpts[num]; double mpts[num]; double npts[num];
    linspace(xpts, -3, 3, num);
    linspace(ypts, -3, 3, num);
    linspace(lpts, -0.2, 0.2, num);
    linspace(mpts, -0.2, 0.2, num);
    linspace(npts, -0.1, 1.3, num);

    for (int i1 = 0; i1 < num; i1++) {
        for (int i2 = 0; i2 < num; i2++) {
            for (int i3 = 0; i3 < num; i3++) {
                for (int i4 = 0; i4 < num; i4++) {
                    for (int i5 = 0; i5 < num; i5++) {
                        l = lpts[i3]; m = mpts[i4]; n = npts[i5];
                        norm = sqrt(l*l+m*m+n*n);
                        l /= norm;
                        m /= norm;
                        n /= norm;

                        ray_in.xt = xpts[i1]; ray_in.yt = ypts[i2];
                        ray_in.l = l; ray_in.m = m; ray_in.n = n;

                        RayTraceSlicer(&ray_out, ray_in, zmin, zmax, 1,
                        p, custom_slice_params, sag_func, transfer_dist_func, surf_normal_func);

                        fprintf(fptr, "%.10f %.10f %.10f %.10f %.10f %.10f %.10f %.10f %.10f\n", xpts[i1], ypts[i2], l, m, n, ray_out.t, ray_out.ln, ray_out.mn, ray_out.nn);
                    }
                }
            }
        }
    }
}

void TestLoadCustomParams(FILE* fptr, int file_num) {
    int array_size = 0;
    double* params = LoadCustomParamsFromFile(file_num, &array_size);
    
    if (params == NULL) {
        printf("Failed to load slice parameters.\n");
        return;
    }

    // Print the first few parameters as a test
    printf("Array Size: %d\n", array_size);
    for (int i = 0; i < array_size; i++) {
        printf("params[%d] = %f\n", i, params[i]);
    }

    free(params);
}