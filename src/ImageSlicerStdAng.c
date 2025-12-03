#define _USE_MATH_DEFINES
#include <windows.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "slicer_generation.h"
#include "slice_param_helpers.h"
#include "surface_solns.h"
#include "triangles.h"

/*
Ellen Lee
Feb 2025

combine flat and standard types. triangle generation for slices is slightly diffferent
but everything else is the same
*/


int __declspec(dllexport) APIENTRY UserObjectDefinition(double *data, double *tri_list);
int __declspec(dllexport) APIENTRY UserParamNames(char *data);
void SetDataFromSlicerParams(IMAGE_SLICER_PARAMS *p, double *data);
void SetSlicerParamsFromData(IMAGE_SLICER_PARAMS *p, double *data);


BOOL WINAPI DllMain (HANDLE hInst, ULONG ul_reason_for_call, LPVOID lpReserved)
	{
   return TRUE;
   }

/* the data is stored as follows:

For calls to UserObjectDefinition:
data[1]  = code indicating what data the DLL should compute (see code sample below)
data[10] - data[12] = data back from the DLL (see code sample below)

data[100] = the number of data parameters passed to the DLL
data[101+] = the value of parameter 1, 2, etc...

For calls to UserParamNames:
data = the parameter number whose name is requested, stored as an ASCII integer;

*/

/*

This DLL models a half cylinder. The parameters are the length and radius of the cylinder.
The number of facets to use in approximating the curve of the cylinder is also a user definable parameter.

*/

static double *p_custom = NULL;
static IMAGE_SLICER_PARAMS_ANGULAR P_OLD = {
        .surface_type = -1,
        .n_each = -1,
        .n_rows = -1,
        .n_cols = -1,
        .angle_mode = -1,
        .dalpha = -1,
        .dbeta = -1,
        .dgamma = -1,
        .gamma_offset = -1,
        .azps = -1,
        .dsyx = -1,
        .dsyz = -1,
        .dsxy = -1,
        .dsxz = -1,
        .du = -1,
        .alpha_cen = -1,
        .beta_cen = -1,
        .gamma_cen = -1,
        .syx_cen = -1,
        .syz_cen = -1,
        .sxy_cen = -1,
        .sxz_cen = -1,
        .u_cen = -1,
        .dx = -1,
        .dy = -1,
        .cv = -1,
        .k = -1,
        .gx_width = -1,
        .gx_depth = -1,
        .gy_width = -1,
        .gy_depth = -1
    };

BOOL WINAPI DllMain(HANDLE hInst, ULONG ul_reason_for_call, LPVOID lpReserved)
{
	switch (ul_reason_for_call)
	{
      case DLL_PROCESS_ATTACH:

          p_custom = (double*)malloc(MAX_ELEMENTS * sizeof(double));

         if (p_custom == NULL) {
            MessageBoxA(NULL, "Memory allocation failed for custom slice parameters", "Error", MB_OK);
            return FALSE;
         }
         for (int i = 0; i < MAX_ELEMENTS; i++) {
             p_custom[i] = 0.0;
         }
         break;

      case DLL_PROCESS_DETACH:
         free(p_custom);
         break;
	}
	return TRUE;
}

int __declspec(dllexport) APIENTRY UserObjectDefinition(double *data, double *tri_list)
	{
	IMAGE_SLICER_PARAMS_ANGULAR p;
	SetSlicerParamsFromData(&p, data);
	ValidateSlicerParams(&p);
	if ( !IsParametersEqualAngular(p, P_OLD) ) {
            MakeSliceParamsArrayAngular(p_custom, p);
            P_OLD = p;
         };
	IMAGE_SLICER_PARAMS_BASIC p_basic = MakeBasicParamsFromCustom(p_custom);

    int Nx = (int) data[101];
    int Ny = (int) data[102];
	double Zdiff = data[103];

	/* what we do now depends upon what was requested */
	code = (int) data[1];

	switch(code)
		{

		/* basic data */
		case 0:
			/* Compute the total number of triangular facets used to render and trace this object */
			/* put this value in data[10] */

			data[10] = CalcNumTriangles(p_basic, Nx, Ny);
			
			/* is this object a solid? put 1 in data[11], use 0 if shell */
			data[11] = 1;

			/* is this object potentially diffractive? put in data[12] the CSG number that is diffractive, cannot use 0. Use 0 if not diffractive */
			data[12] = 0;

			return 0;

		/* generate triangles */
		case 1:
			{
			/*
			We are being asked to generate the triangle list.
			A triangle consists of 3 triplets of coordinates (x, y, and z for each of the 3 corners),
			plus a 6 digit integer that indicates the properties of the triangle.
			That is 10 numbers for each triangle.
			The triangle data will be stored in the tri_list array

			The 6 digit number is constructed as follows.

			The "visible" flag is in the 1's place. If all 3 sides of a
			triangle are visible, use 0. If the side from point 1 to point
			2 is invisible, use 1. If the side from point 2 to point 3 is
			invisible, use 2. If the side from point 3 to point 1 is invisible,
			use 4. Add these flags if more than one side is invisible. For example,
			if all 3 sides are invisible, use 7. Note invisible simply means the
			side is not drawn on the 3D Layout plot, it has nothing to do with ray
			tracing.

			The "reflective" flag is in the 10's place, use zero for refractive and 1 for reflective.

			The "coat scatter group" is a 2 digit integer in the 1,000 and 100's place.
			Only CSG 0-3 are currently supported, but 2 digits are provided for future
			expansion.

			The "exact group" flag is a 2 digit integer in the 100,000 and 10,000's place.
			The exact group number indicates which set of exact formulas are used to iterate
			to the actual object surface. Use 0 to indicate the flat triangular facet is the
			correct solution. The exact flags are used in the intercept portion of the code that is case 2 below.

			Example: 030212 means the exact group is 3, the CSG is 2, the facet
			is reflective, and the side of the triangle from point 2 to point 3 is not drawn.

			*/

			int num_triangles = 0;
			MakeAllTrianglesForSlicer(tri_list, &num_triangles, Nx, Ny, p_basic, p_custom);
			}
			break;

		/* iterate to exact solution given starting point */
		case 2:
			{
			/*
			We are giving starting data, and must compute the distance to the actual surface.
			We also need to compute the normal vector at the point.

			The data ZEMAX sends is formatted as follows:
			data[2], data[3], data[4] = x, y, z
			data[5], data[6], data[7] = l, m, n
			data[8] = exact code

			*/

			switch((int)data[8])
				{
				default:
					/* no need to iterate, assume exact solution is on flat triangle face */
					break;
				case 1:
					/* slice face */
					{
					int col_num, slice_num;
					double x, y, z, l, m, n;
					x = data[2];
					y = data[3];
					z = data[4];
					l = data[5];
					m = data[6];
					n = data[7];
					RAY_IN ray_in = {x, y, z, l, m, n};

					// Use x and y to figure which slice to use
					GetSlicerIndex(&col_num, &slice_num, x, y, p);
					SLICE_PARAMS pslice = GetSliceParams(slice_num, col_num, p_custom);

					TRANSFER_DIST_FUNC transfer_dist_func;
					SURF_NORMAL_FUNC surf_normal_func;
					CRITICAL_XY_FUNC critical_xy_func;
					TRANSFORMATION_FUNC transform_func;
					GetSurfaceFuncs(&transfer_dist_func, &surf_normal_func, &critical_xy_func, &transform_func, pslice, p);

					RAY_OUT ray_out = SliceRayTrace(ray_in, pslice, transfer_dist_func, surf_normal_func, transform_func, 1);

					if (isnan(ray_out.t)) return 1; // something went wrong... ray missed somehow?

					if (isnan(ray_out.ln) || isnan(ray_out.mn) || isnan(ray_out.nn)) return 1;

					data[10] = ray_out.t;
					data[11] = ray_out.ln;
					data[12] = ray_out.mn;
					data[13] = ray_out.nn;
					return 0;
					}
					break;
				}
			}
			break;

		/* coating data */
		case 3:
			return -1;

		/* safe data */
		case 4:
			/* set safe parameter data values the first time the DLL is initialized */
			data[101] = 4; 		  // Nx
			data[102] = 3; 		  // Ny
			data[103] = 5;		  // Zdiff

			 // slicer params
			data[104] = -0.01;    // cv
			data[105] = 0.0;      // k
			data[106] = (int) 0;  // surface_type
			data[107] = (int) 6;  // n_each
			data[108] = (int) 2;  // n_rows
			data[109] = (int) 1;  // n_cols
			data[110] = (int) 2;  // angle_mode
			data[111] = 4.0;      // dalpha
			data[112] = 0.0;      // dbeta
			data[113] = 0.5;      // dgamma
			data[114] = 0.0;      // gamma_offset
			data[115] = 0.0;      // azps
			data[116] = 0.0;      // dsyx
			data[117] = 0.0;      // dsyz
			data[118] = 0.0;      // du
			data[119] = 0.0;      // alpha_cen
			data[120] = 0.0;      // beta_cen
			data[121] = 0.0;      // gamma_cen
			data[122] = 4.0;      // syx_cen
			data[123] = 0.0;      // syz_cen
			data[124] = 0.0;      // sxy_cen
			data[125] = 0.0;      // sxz_cen
			data[126] = 0.0;      // u_cen
			data[127] = 8.0;      // dx
			data[128] = 0.75;     // dy
			data[129] = 0.05;     // gx_width
			data[130] = 0.0;      // gx_depth
			data[131] = 0.0;      // gy_width
			data[132] = 0.0;      // gy_depth

			SetSlicerParamsFromData(&p, data);
			return 0;
		}

	/* we did not recognize the request */
	return -1;
   }

int __declspec(dllexport) APIENTRY UserParamNames(char *data) {
	/* this function returns the name of the parameter requested */
	int i;
	i = atoi(data);
	strcpy(data,""); // blank means unused

	/* negative numbers or zero mean names for coat/scatter groups */
	if (i ==  0) strcpy(data,"Slicer Face");
	if (i == -1) strcpy(data,"Wall Face");
	if (i == -2) strcpy(data,"Gap Face");
	if (i == -3) strcpy(data,"Shell Face");

	if (i == 1) strcpy(data,"Nx");
	if (i == 2) strcpy(data,"Ny");
	if (i == 3) strcpy(data,"Zdiff");
	if (i == 4) strcpy(data,"cv");
	if (i == 5) strcpy(data,"k");
	if (i == 6) strcpy(data,"surface_type");
	if (i == 7) strcpy(data,"n_each");
	if (i == 8) strcpy(data,"n_rows");
	if (i == 9) strcpy(data,"n_cols");
	if (i == 10) strcpy(data,"angle_mode");
	if (i == 11) strcpy(data,"dalpha");
	if (i == 12) strcpy(data,"dbeta");
	if (i == 13) strcpy(data,"dgamma");
	if (i == 14) strcpy(data,"gamma_offset");
	if (i == 15) strcpy(data,"azps");
	if (i == 16) strcpy(data,"dsyx");
	if (i == 17) strcpy(data,"dsyz");
	if (i == 18) strcpy(data,"dsxy");
	if (i == 19) strcpy(data,"dsxz");
	if (i == 20) strcpy(data,"du");
	if (i == 21) strcpy(data,"alpha_cen");
	if (i == 22) strcpy(data,"beta_cen");
	if (i == 23) strcpy(data,"gamma_cen");
	if (i == 24) strcpy(data,"syx_cen");
	if (i == 25) strcpy(data,"syz_cen");
	if (i == 26) strcpy(data,"sxy_cen");
	if (i == 27) strcpy(data,"sxz_cen");
	if (i == 28) strcpy(data,"u_cen");
	if (i == 29) strcpy(data,"dx");
	if (i == 30) strcpy(data,"dy");
	if (i == 31) strcpy(data,"gx_width");
	if (i == 32) strcpy(data,"gx_depth");
	if (i == 33) strcpy(data,"gy_width");
	if (i == 34) strcpy(data,"gy_depth");
	return 0;
}

// Functions to convert between Zemax data and image slicer params struct
void SetDataFromSlicerParams(IMAGE_SLICER_PARAMS *p, double *data) {
	data[104] = p->cv;
	data[105] = p->k;
	data[106] = (int) p->surface_type;
	data[107] = (int) p->n_each;
	data[108] = (int) p->n_rows;
	data[109] = (int) p->n_cols;
	data[110] = (int) p->angle_mode;
	data[111] = p->dalpha;
	data[112] = p->dbeta;
	data[113] = p->dgamma;
	data[114] = p->gamma_offset;
	data[115] = p->azps;
	data[116] = p->dsyx;
	data[117] = p->dsyz;
	data[118] = p->dsxy;
	data[119] = p->dsxz;
	data[120] = p->du;
	data[121] = p->alpha_cen;
	data[122] = p->beta_cen;
	data[123] = p->gamma_cen;
	data[124] = p->syx_cen;
	data[125] = p->syz_cen;
	data[126] = p->sxy_cen;
	data[127] = p->sxz_cen;
	data[128] = p->u_cen;
	data[129] = p->dx;
	data[130] = p->dy;
	data[131] = p->gx_width;
	data[132] = p->gx_depth;
	data[133] = p->gy_width;
	data[134] = p->gy_depth;
}

void SetSlicerParamsFromData(IMAGE_SLICER_PARAMS *p, double *data) {
	p->cv = data[104];
	p->k = data[105];
	p->surface_type = (int) data[106];
	p->n_each = (int) data[107];
	p->n_rows = (int) data[108];
	p->n_cols = (int) data[109];
	p->angle_mode = (int) data[110];
	p->dalpha = data[111];
	p->dbeta = data[112];
	p->dgamma = data[113];
	p->gamma_offset = data[114];
	p->azps = data[115];
	p->dsyx = data[116];
	p->dsyz = data[117];
	p->dsxy = data[118];
	p->dsxz = data[119];
	p->du = data[120];
	p->alpha_cen = data[121];
	p->beta_cen = data[122];
	p->gamma_cen = data[123];
	p->syx_cen = data[124];
	p->syz_cen = data[125];
	p->sxy_cen = data[126];
	p->sxz_cen = data[127];
	p->u_cen = data[128];
	p->dx = data[129];
	p->dy = data[130];
	p->gx_width = data[131];
	p->gx_depth = data[132];
	p->gy_width = data[133];
	p->gy_depth = data[134];
}