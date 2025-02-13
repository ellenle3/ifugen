#define _USE_MATH_DEFINES
#include <windows.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "slicer_generation.h"
#include "surface_solns.h"

/*
Ellen Lee
Feb 2025
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

int __declspec(dllexport) APIENTRY UserObjectDefinition(double *data, double *tri_list)
	{
	IMAGE_SLICER_PARAMS p;
	SetSlicerParamsFromData(&p, data);

	// For computing how many triangles we need
	int num_slices_total = 2;
	int num_gaps_x = 0, num_gaps_y = 0;
	int num_walls_x = 0, num_walls_y = 0;
	int num_triangles_surface = 2, num_triangles_sides = 2;

	// Validate parameters
	ValidateSlicerParams(&p);

	/* what we do now depends upon what was requested */
	code = (int) data[1];

	switch(code)
		{

		/* basic data */
		case 0:
			/* Compute the total number of triangular facets used to render and trace this object */
			/* put this value in data[10] */

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
			num_triangles_surface = 2*num_slices_total + 2*num_walls_x + 2*num_walls_y + 2*(num_gaps_x + num_gaps_y)
			// Need to do 5 more panels - 4 for the sides and one for the back
			num_triangles_sides = 2 + 2 * (num_slices_total + num_gaps_x) + 2 * (num_slices_total + num_gaps_y)
			data[10] = num_triangles_surface + num_triangles_sides;
			
			/* is this object a solid? put 1 in data[11], use 0 if shell */
			data[11] = 1;

			/* is this object potentially diffractive? put in data[12] the CSG number that is diffractive, cannot use 0. Use 0 if not diffractive */
			data[12] = 0;

			return 0;

		/* generate triangles */
		case 1:
			{
			int num_triangles;
			int col_num = 0, slice_num = 0;
			int i = 0, j = 0;
			double xstart, xend, xstep;
			double ystart, yend, ystep;
			double pt1, pt2, pt3, pt4;
			double code1, code2;

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

			num_triangles = 0;

			xstep = p.dx / (Nx - 1); ystep = p.dy / (Ny - 1);
			// First, do the surfaces of all of the slices
			for (col_num = 0, col_num < p.n_cols, col_num++) {

				for (slice_num = 0, slice_num < p.n_each * p.n_rows, slice_num++) {

                    code1 = code2 = 010017.0;   // exact 1, CSG 0, reflective, all invisible

                    // Determine face visibility
                    if (i == 0) {code1 -= 4;}  		 // side from point 3 to 1 is visible on triangle 1
                    elif (i == Nx - 1) {code2 -= 4;} // 3 to 1 is visible on triangle 2
                    if (j == 0) {code1 -= 1;} 		 // 1 to 2 is visible on triangle 1
                    elif (j == Ny - 1) {code -= 2;}  // 2 to 3 is visible on triangle 2


                    x1 = xstart + xstep*i  // you don't need to recompute these!
                    x2 = xstart + xstep*(i+1) // only needs to be computed if i=0, j=0
                    y1 = ystart + ystep*j  // recycle somehow...
                    y2 = ystart + ystep*(j+1)
                    z1 = ImageSlicerSag(x1, y1, p, TiltedPlaneSag);
                    z2 = ImageSlicerSag(x2, y1, p, TiltedPlaneSag);
                    z3 = ImageSlicerSag(x1, y2, p, TiltedPlaneSag);
                    z4 = ImageSlicerSag(x2, y2, p, TiltedPlaneSag);

                    tri_list[num_triangles*10 + 0] = x1;  // x1
                    tri_list[num_triangles*10 + 1] = y1;  // y1
                    tri_list[num_triangles*10 + 2] = z1;  // z1
                    tri_list[num_triangles*10 + 3] = x2;  // x2
                    tri_list[num_triangles*10 + 4] = y1;  // y2
                    tri_list[num_triangles*10 + 5] = z2;  // z2
                    tri_list[num_triangles*10 + 6] = x1;  // x3
                    tri_list[num_triangles*10 + 7] = y2;  // y3
                    tri_list[num_triangles*10 + 8] = z3;  // z3
                    tri_list[num_triangles*10 + 9] = code1;
                    num_triangles++;

                    tri_list[num_triangles*10 + 0] = x2;  // x1
                    tri_list[num_triangles*10 + 1] = y1;  // y1
                    tri_list[num_triangles*10 + 2] = z2;  // z1
                    tri_list[num_triangles*10 + 3] = x1;  // x2
                    tri_list[num_triangles*10 + 4] = y2;  // y2
                    tri_list[num_triangles*10 + 5] = z3;  // z2
                    tri_list[num_triangles*10 + 6] = x2;  // x3
                    tri_list[num_triangles*10 + 7] = y2;  // y3
                    tri_list[num_triangles*10 + 8] = z4;  // z3
                    tri_list[num_triangles*10 + 9] = code2;
                    num_triangles++;
				}
			}

			// Now do x-walls

			// ...y-walls

			// Do any x-gaps
			for (col_num = 0, col_num < p.n_cols, col_num++) {
			}

			// ...y-gaps

			// if xgap and ygap depths are different these also need walls! Find
			// number of intersection points between gaps. 4 walls (8) triangles
			// for each one... No need if either xgap or ygap size is 0.

			// Finish by doing the side + back panels

			data[10] = num_triangles; /* how many we actually wrote out */
			return 0;
			}

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
					double x, y, z, l, m, n, tp, xt, yt, t, xs, ys Fx, Fy, Fz;
					int col_num, slice_num;
					double alpha, beta, gamma;

					x = data[2];
					y = data[3];
					z = data[4];
					l = data[5];
					m = data[6];
					n = data[7];
					t = 0.0;
					// Use x and y to figure which slice to use
					GetSlicerIndex(&col_num, &slice_num, x, y, p);
					GetSliceAngles(&alpha, &beta, &gamma, slice_num, col_num, p);

					// Our expression for the transfer distance is dependent on
					// knowing (xt, yt, 0) for the ray. Calculate xt and yt from (x, y, z).
					if (fabs(n) < 1E-13) {n = 1E-13;} // avoid divide by zero
					tp = z / n;
					xt = x - tp * l; yt = y - tp * m;
					// Analytic solution for the transfer distance for this slice		
					t = TiltedPlaneTransfer(xt, yt, l, m, n, p.cv, p.k, alpha, beta, gamma);
					if (isnan(t)) return -1; // something went wrong... ray missed somehow?

					// Compute surface normal from xs and ys
					xs = xt + t*l;
					ys = yt + t*m;
					TiltedPlaneSurfaceNormal(&Fx, &Fy, &Fz, xs, ys, p.cv, p.k, alpha, beta, gamma, 1);
					if (isnan(Fx) || isnan(Fy) || isnan(Fz)) return 1;
					
					data[10] = t;
					data[11] = Fx;
					data[12] = Fy;
					data[13] = Fz;
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
			R = 1.0;
			L = 1.0;
			N = 6;
			data[101] = R;
			data[102] = L;
			data[103] = (double) N;
			return 0;

		}

	/* we did not recognize the request */
	return -1;
   }

int __declspec(dllexport) APIENTRY UserParamNames(char *data)
	{
	/* this function returns the name of the parameter requested */
	int i;
	i = atoi(data);
	strcpy(data,""); // blank means unused

	/* negative numbers or zero mean names for coat/scatter groups */
	if (i ==  0) strcpy(data,"Slicer Face");
	if (i == -1) strcpy(data,"Wall Face");
	if (i == -2) strcpy(data,"Gap Face");
	if (i == -3) strcpy(dta,"Outside Face")

	if (i == 1) strcpy(data,"n_each");
	if (i == 2) strcpy(data,"n_rows");
	if (i == 3) strcpy(data,"n_cols");
	if (i == 4) strcpy(data,"mode");
	if (i == 5) strcpy(data,"trace_walls");
	if (i == 6) strcpy(data,"active_x");
	if (i == 7) strcpy(data,"active_y");
	if (i == 8) strcpy(data,"dalpha");
	if (i == 9) strcpy(data,"dbeta");
	if (i == 10) strcpy(data,"dgamma");
	if (i == 11) strcpy(data,"alpha_cen");
	if (i == 12) strcpy(data,"beta_cen");
	if (i == 13) strcpy(data,"gamma_cen");
	if (i == 14) strcpy(data,"dx");
	if (i == 15) strcpy(data,"dy");
	if (i == 16) strcpy(data,"gx_width");
	if (i == 17) strcpy(data,"gx_depth");
	if (i == 18) strcpy(data,"gy_width");
	if (i == 19) strcpy(data,"gy_depth");
	if (i == 20) strcpy(data,"# Facets x");
	if (i == 21) strcpy(data,"# Facets y");
	
	return 0;
	}

// Functions to convert between Zemax data and image slicer params struct
void SetDataFromSlicerParams(IMAGE_SLICER_PARAMS *p, double *data) {
	data[101] = p->n_each;
	data[102] = p->n_rows;
	data[103] = p->n_cols;
	data[104] = p->mode;
	data[105] = p->trace_walls;
	data[106] = p->active_x;
	data[107] = p->active_y;
	data[108] = p->dalpha;
	data[109] = p->dbeta;
	data[110] = p->dgamma;
	data[111] = p->alpha_cen;
	data[112] = p->beta_cen;
	data[113] = p->gamma_cen;
	data[114] = p->dx;
	data[115] = p->dy;
	data[116] = p->gx_width;
	data[117] = p->gx_depth;
	data[118] = p->gy_width;
	data[119] = p->gy_depth;
}

void SetSlicerParamsFromData(IMAGE_SLICER_PARAMS *p, double *data) {
	p->n_each = (int) data[101];
	p->n_rows = (int) data[102];
	p->n_cols = (int) data[103];
	p->mode = (int) data[104];
	p->trace_walls = (int) data[105];
	p->active_x = (int) data[106];
	p->active_y = (int) data[107];
	p->dalpha = data[108];
	p->dbeta = data[109];
	p->dgamma = data[110];
	p->alpha_cen = data[111];
	p->beta_cen = data[112];
	p->gamma_cen = data[113];
	p->dx = data[114];
	p->dy = data[115];
	p->gx_width = data[116];
	p->gx_depth = data[117];
	p->gy_width = data[118];
	p->gy_depth = data[119];
}