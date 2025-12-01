#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "slicer_generation.h"
#include "surface_solns.h"
#include "slice_param_helpers.h"
#include "triangles.h"

void TestImageSlicerSag(FILE* fptr, IMAGE_SLICER_PARAMS_BASIC p, double p_custom[]);
void TestGlobalExtrema(FILE* fptr, IMAGE_SLICER_PARAMS_BASIC p, double p_custom[]);
void TestRayTrace(FILE* fptr, IMAGE_SLICER_PARAMS_BASIC p, double p_custom[]);
void TestLoadCustomParams(FILE* fptr, int file_num);

int main() {

    FILE* fptr = fopen("triangles.txt", "w+");

    if (fptr==NULL) {
        printf("The file could not be opened.");
        return 1;
    }
    
    // Set up the image slicer

    IMAGE_SLICER_PARAMS_ANGULAR p = {
        .surface_type = 0,
        .n_each = 2,
        .n_rows = 2,
        .n_cols = 3,
        .angle_mode = 3,

        .dalpha = 0,
        .dbeta = 0,
        .dgamma = 10,
        .gamma_offset = 0,

        .azps = 0,
        .dsyx = 0,
        .dsyz = 0,
        .dsxy = 0,
        .dsxz = 0,
        .du = 4,

        .alpha_cen = 0,
        .beta_cen = 0,
        .gamma_cen = 0,
        .syx_cen = 0,
        .syz_cen = 0,
        .sxy_cen = 0,
        .sxz_cen = 0,
        .u_cen = 0,

        .dx = 6,
        .dy = 1,
        .gx_width = 1,
        .gx_depth = 1,
        .gy_width = 1,
        .gy_depth = 2,
        .cv = -0.03,
        .k = 0
    };

    double* p_custom = (double*)malloc(MAX_ELEMENTS * sizeof(double));
    //LoadCustomParamsFromFile(p_custom, 0, "/Users/ellenlee/Documents/Zemax_dll/ifugen/python/", 10000);
    MakeSliceParamsArrayAngular(p_custom, p);
    IMAGE_SLICER_PARAMS_BASIC p_basic = MakeBasicParamsFromCustom(p_custom);
    int Nx = 2;
    int Ny = 2;
    int Ntotal = CalcNumTriangles(p_basic, Nx, Ny);
    double* tri_list = (double*)malloc(Ntotal * 10 * sizeof(double));
    int num_triangles = 0;
    MakeAllTrianglesForSlicer(tri_list, &num_triangles, Nx, Ny, p_basic, p_custom);

    for (int i = 0; i < num_triangles; i++) {
    int base = i * 10;
    fprintf(fptr, "%f,%f,%f,%f,%f,%f,%f,%f,%f,%d\n",
        tri_list[base+0], tri_list[base+1], tri_list[base+2],
        tri_list[base+3], tri_list[base+4], tri_list[base+5],
        tri_list[base+6], tri_list[base+7], tri_list[base+8],
        (int)tri_list[base+9]);
    }


    free(p_custom);
    free(tri_list);
    fclose(fptr);
    return 0;
}

int main_stdmode() {
    FILE* fptr = fopen("test_output.txt", "w+");

    if (fptr==NULL) {
        printf("The file could not be opened.");
        return 1;
    }
    
    // Set up the image slicer

    IMAGE_SLICER_PARAMS_ANGULAR p = {
        .surface_type = 0,
        .n_each = 4,
        .n_rows = 3,
        .n_cols = 1,
        .angle_mode = 1,

        .dalpha = 0,
        .dbeta = 0,
        .dgamma = 5,
        .gamma_offset = 0,

        .azps = 0,
        .dsyx = 0,
        .dsyz = 0,
        .dsxy = 0,
        .dsxz = 0,
        .du = 0,

        .alpha_cen = 0,
        .beta_cen = 0,
        .gamma_cen = 0,
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

    double* p_custom = (double*)malloc(MAX_ELEMENTS * sizeof(double));
    //LoadCustomParamsFromFile(p_custom, 0, "/Users/ellenlee/Documents/Zemax_dll/ifugen/python/", 10000);
    MakeSliceParamsArrayAngular(p_custom, p);
    IMAGE_SLICER_PARAMS_BASIC p_basic = MakeBasicParamsFromCustom(p_custom);

    SLICE_PARAMS pslice = GetSliceParams(9, 0, p_custom); // Test function call
    //SLICE_PARAMS pslice = GetSliceParamsAngular(2, 0, p); // Test function call
    TestRayTrace(fptr, p_basic, p_custom);
    //TestImageSlicerSag(fptr, p_basic, p_custom);

    int array_size = NUM_BASE_PARAMS + p_basic.n_rows + NUM_PARAMS_PER_SLICE * p_basic.n_rows * p_basic.n_cols * p_basic.n_each;
    FILE *file = fopen("test_output-2.txt", "wb");
    for (int i = 0; i < array_size; i++) {
        fprintf(file, "%.10f\n", p_custom[i]);
    }
    fclose(file);

    free(p_custom);
    //TestGlobalExtrema(fptr, p);

    fclose(fptr);
    return 0;
}

void TestImageSlicerSag(FILE* fptr, IMAGE_SLICER_PARAMS_BASIC p, double p_custom[]) {
    
    // Compute test function over an array of points and save to a TXT file
    double xsize, ysize;
    double output;
    GetSlicerSize(&xsize, &ysize, p);

    int nx = 200; int ny = 200;
    double xpts[nx]; double ypts[ny];
    linspace(xpts, -xsize/2 - 3, xsize/2 + 3, nx);
    linspace(ypts, -ysize/2 - 3, ysize/2 + 3, ny);

    int nc, ns;
    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            output = ImageSlicerSag(xpts[i], ypts[j], p, p_custom);
            GetSlicerIndex(&nc, &ns, xpts[i], ypts[j], p, p_custom);
            fprintf(fptr, "%.15f %.15f %.15f %d %d\n", xpts[i], ypts[j], output, nc, ns);
        }
    }
}

void TestGlobalExtrema(FILE* fptr, IMAGE_SLICER_PARAMS_BASIC p, double p_custom[]) {

    double zmin, zmax;
    FindSlicerGlobalExtrema(&zmin, &zmax, p, p_custom);
    fprintf(fptr, "%.10f %.10f\n", zmin, zmax);

}

void TestRayTrace(FILE* fptr, IMAGE_SLICER_PARAMS_BASIC p, double p_custom[]) {

    RAY_IN ray_in;
    RAY_OUT ray_out;

    double zmin, zmax, umin, umax, norm, l, m, n;
    FindSlicerGlobalExtrema(&zmin, &zmax, p, p_custom);
    GetMinMaxU(&umin, &umax, p, p_custom);

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

                        ray_in.xt = xpts[i1]; ray_in.yt = ypts[i2]; ray_in.zt = 0;
                        ray_in.l = l; ray_in.m = m; ray_in.n = n;

                        RayTraceSlicer(&ray_out, ray_in, zmin, zmax, umin, umax, 1,
                        p, p_custom);

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
    LoadCustomParamsFromFile(params, file_num, params_dir);

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