#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "slicer_generation.h"
#include "surface_solns.h"
#include "custom_slicer_helpers.h"

#define MAX_ELEMENTS 6250009

void TestImageSlicerSag(FILE* fptr, IMAGE_SLICER_PARAMS p, double p_custom[]);
void TestGlobalExtrema(FILE* fptr, IMAGE_SLICER_PARAMS p);
void TestRayTrace(FILE* fptr, IMAGE_SLICER_PARAMS p);
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
        .surface_type = 0,
        .n_each = 4,
        .n_rows = 3,
        .n_cols = 1,
        .angle_mode = 1,

        .dalpha = 0,
        .dbeta = 0,
        .dgamma = 2,
        .gamma_offset = 0,

        .dzps = 0,
        .dzp_col = 0,
        .dzp_row = 0,
        .dsyx = 0,
        .dsyz = 0,
        .dsxy = 0,
        .dsxz = 0,
        .du = 0,

        .alpha_cen = 0,
        .beta_cen = 0,
        .gamma_cen = 0,
        .zps_cen = 0,
        .zp_cen = 0,
        .syx_cen = 0,
        .syz_cen = 0,
        .sxy_cen = 0,
        .sxz_cen = 0,
        .u_cen = 0,

        .dx = 12,
        .dy = 1,
        .gx_width = 0,
        .gx_depth = 0,
        .gy_width = 0,
        .gy_depth = 0,
        .cv = -0.01,
        .k = 0
    };

    double p_custom[1] = {0};

    TestImageSlicerSag(fptr, p, p_custom);
    //TestGlobalExtrema(fptr, p);
    //TestRayTrace(fptr, p);

    //double* params = (double *)calloc(MAX_ELEMENTS, sizeof(double));
    //char params_dir[] = "/Users/ellenlee/Documents/Zemax_dll/ifugen/python/";
    //LoadCustomParamsFromFile(params, 1, params_dir, MAX_ELEMENTS);
    //p = MakeSlicerParamsFromCustom(params);

    //TestImageSlicerSag(fptr, p, params);

    fclose(fptr);
}

void TestImageSlicerSag(FILE* fptr, IMAGE_SLICER_PARAMS p, double p_custom[]) {
    
    // Compute test function over an array of points and save to a TXT file
    double xsize, ysize;
    double output;
    GetSlicerSize(&xsize, &ysize, p);

    int nx = 200; int ny = 200;
    double xpts[nx]; double ypts[ny];
    linspace(xpts, -xsize/2 - 3, xsize/2 + 3, nx);
    linspace(ypts, -ysize/2 - 3, ysize/2 + 3, ny);

    // for (int i = 0; i < nx; i++) {
    //     output = ImageSlicerSag(xpts[i], ysize/2, p, p_custom);
    //     fprintf(fptr, "%.10f %.10f %.10f\n", xpts[i], ysize/2, output);
    // }

    int nc, ns;
    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            output = ImageSlicerSag(xpts[i], ypts[j], p, p_custom);
            GetSlicerIndex(&nc, &ns, xpts[i], ypts[j], p, p_custom);
            fprintf(fptr, "%.15f %.15f %.15f %d %d\n", xpts[i], ypts[j], output, nc, ns);
        }
    }
}

void TestGlobalExtrema(FILE* fptr, IMAGE_SLICER_PARAMS p) {
    
    double zmin, zmax;
    double p_custom[1] = {0};
    FindSlicerGlobalExtrema(&zmin, &zmax, p, p_custom);
    fprintf(fptr, "%.10f %.10f\n", zmin, zmax);

}

void TestRayTrace(FILE* fptr, IMAGE_SLICER_PARAMS p) {

    RAY_IN ray_in;
    RAY_OUT ray_out;

    double zmin, zmax, norm, l, m, n;
    double p_custom[1] = {0};
    FindSlicerGlobalExtrema(&zmin, &zmax, p, p_custom);

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

                        //RayTraceSlicer(&ray_out, ray_in, zmin, zmax, 1,
                        //p, p_custom);

                        fprintf(fptr, "%.10f %.10f %.10f %.10f %.10f %.10f %.10f %.10f %.10f\n", xpts[i1], ypts[i2], l, m, n, ray_out.t, ray_out.ln, ray_out.mn, ray_out.nn);
                    }
                }
            }
        }
    }
}

void TestLoadCustomParams(FILE* fptr, int file_num) {
    double* params = (double *)calloc(MAX_ELEMENTS, sizeof(double));
    char params_dir[] = "/Users/ellenlee/Documents/Zemax_dll/ifugen/python/";
    LoadCustomParamsFromFile(params, file_num, params_dir, MAX_ELEMENTS);

    int n_slices_per_col = params[0];
    int n_cols = params[1];
    int array_size = 9 + n_slices_per_col * n_cols * 5;

    // Print the first few parameters as a test
    printf("Array Size: %d\n", array_size);
    for (int i = 0; i < array_size; i++) {
        printf("params[%d] = %f\n", i, params[i]);
    }

    free(params);
}