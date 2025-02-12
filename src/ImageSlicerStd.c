#include <windows.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

/*
Written by Kenneth E. Moore
January 17, 2003
Modified December 16, 2003 to support the "diffractive flag" on data[12].
Modified July 29, 2005 to support the "safe values" code = 3.
*/

// modified KEM 4-12-2006 to include coating placeholders, safe data moved to code = 4


int __declspec(dllexport) APIENTRY UserObjectDefinition(double *data, double *tri_list);
int __declspec(dllexport) APIENTRY UserParamNames(char *data);

#define PI 3.14159265358979323846

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
	int Nx, Ny, code;
    double R, L;

    R = data[101];
	L = data[102];
    Nx = (int) data[103];
    Nx = (int) data[103];

	/* do some basic error trapping and handling */
	if (R <= 0.0) R = 1.0;
	if (L <= 0.0) L = 1.0;
	if (N < 6) N = 6;

	/* what we do now depends upon what was requested */
	code = (int) data[1];

	switch(code)
		{

		/* basic data */
		case 0:
			/* Compute the total number of triangular facets used to render and trace this object */
			/* we need N on each end cap, 2N along the face, and 2 on the flat back */
			/* put this value in data[10] */
			data[10] = 4*N+2;
			
			/* is this object a solid? put 1 in data[11], use 0 if shell */
			data[11] = 1;

			/* is this object potentially diffractive? put in data[12] the CSG number that is diffractive, cannot use 0. Use 0 if not diffractive */
			data[12] = 0;

			return 0;

		/* generate triangles */
		case 1:
			{
			int num_triangles, t;
			double a1, a2;

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

			/* first do the bottom face 2 triangles */
			/* note we use the visible flag to hide the diagonal line across the face */
			tri_list[num_triangles*10 + 0] = R;   // x1
			tri_list[num_triangles*10 + 1] = 0.0; // y1
			tri_list[num_triangles*10 + 2] = 0.0; // z1
			tri_list[num_triangles*10 + 3] = R;   // x2
			tri_list[num_triangles*10 + 4] = 0.0; // y2
			tri_list[num_triangles*10 + 5] = L;   // z2
			tri_list[num_triangles*10 + 6] = -R;  // x3
			tri_list[num_triangles*10 + 7] = 0.0; // y3
			tri_list[num_triangles*10 + 8] = L;   // z3
			tri_list[num_triangles*10 + 9] = 000104.0;   // exact 0, CSG 1, refractive, 3-1 invisible
			num_triangles++;

			tri_list[num_triangles*10 + 0] = R;   // x1
			tri_list[num_triangles*10 + 1] = 0.0; // y1
			tri_list[num_triangles*10 + 2] = 0.0; // z1
			tri_list[num_triangles*10 + 3] = -R;  // x2
			tri_list[num_triangles*10 + 4] = 0.0; // y2
			tri_list[num_triangles*10 + 5] = 0.0; // z2
			tri_list[num_triangles*10 + 6] = -R;  // x3
			tri_list[num_triangles*10 + 7] = 0.0; // y3
			tri_list[num_triangles*10 + 8] = L;   // z3
			tri_list[num_triangles*10 + 9] = 000104.0;   // exact 0, CSG 1, refractive, 3-1 invisible
			num_triangles++;

			/* now do the end panel at z = 0 */
			for (t = 1; t <= N; t++)
				{
				a1 = PI * (double) (t - 1) / (double) N;
				a2 = PI * (double) (t    ) / (double) N;

				tri_list[num_triangles*10 + 0] = 0.0; // x1
				tri_list[num_triangles*10 + 1] = 0.0; // y1
				tri_list[num_triangles*10 + 2] = 0.0; // z1
				tri_list[num_triangles*10 + 3] = R*cos(a1);  // x2
				tri_list[num_triangles*10 + 4] = R*sin(a1);  // y2
				tri_list[num_triangles*10 + 5] = 0.0; // z2
				tri_list[num_triangles*10 + 6] = R*cos(a2);  // x3
				tri_list[num_triangles*10 + 7] = R*sin(a2);  // y3
				tri_list[num_triangles*10 + 8] = 0.0;  // z3
				if (t == 1)
					{
					tri_list[num_triangles*10 + 9] = 000204.0;   // exact 0, CSG 2, refractive, 3-1 invisible
					}
				if (t > 1 && t < N)
					{
					tri_list[num_triangles*10 + 9] = 000205.0;   // exact 0, CSG 2, refractive, 1-2 and 3-1 invisible
					}
				if (t == N)
					{
					tri_list[num_triangles*10 + 9] = 000201.0;   // exact 0, CSG 2, refractive, 1-2 invisible
					}
				num_triangles++;
				}

			/* now do the end panel at z = L, almost identical to code above */
			for (t = 1; t <= N; t++)
				{
				a1 = PI * (double) (t - 1) / (double) N;
				a2 = PI * (double) (t    ) / (double) N;

				tri_list[num_triangles*10 + 0] = 0.0; // x1
				tri_list[num_triangles*10 + 1] = 0.0; // y1
				tri_list[num_triangles*10 + 2] = L; // z1
				tri_list[num_triangles*10 + 3] = R*cos(a1);  // x2
				tri_list[num_triangles*10 + 4] = R*sin(a1);  // y2
				tri_list[num_triangles*10 + 5] = L; // z2
				tri_list[num_triangles*10 + 6] = R*cos(a2);  // x3
				tri_list[num_triangles*10 + 7] = R*sin(a2);  // y3
				tri_list[num_triangles*10 + 8] = L;  // z3
				if (t == 1)
					{
					tri_list[num_triangles*10 + 9] = 000304.0;   // exact 0, CSG 3, refractive, 3-1 invisible
					}
				if (t > 1 && t < N)
					{
					tri_list[num_triangles*10 + 9] = 000305.0;   // exact 0, CSG 3, refractive, 1-2 and 3-1 invisible
					}
				if (t == N)
					{
					tri_list[num_triangles*10 + 9] = 000301.0;   // exact 0, CSG 3, refractive, 1-2 invisible
					}
				num_triangles++;
				}

			/* finally, do the triangles that make up the cylinder face, note these are in pairs */
			for (t = 1; t <= N; t++)
				{
				a1 = PI * (double) (t - 1) / (double) N;
				a2 = PI * (double) (t    ) / (double) N;

				tri_list[num_triangles*10 + 0] = R*cos(a1); // x1
				tri_list[num_triangles*10 + 1] = R*sin(a1); // y1
				tri_list[num_triangles*10 + 2] = 0.0; // z1
				tri_list[num_triangles*10 + 3] = R*cos(a1);  // x2
				tri_list[num_triangles*10 + 4] = R*sin(a1);  // y2
				tri_list[num_triangles*10 + 5] = L; // z2
				tri_list[num_triangles*10 + 6] = R*cos(a2);  // x3
				tri_list[num_triangles*10 + 7] = R*sin(a2);  // y3
				tri_list[num_triangles*10 + 8] = L;  // z3
				tri_list[num_triangles*10 + 9] = 010004.0;   // exact 1, CSG 0, refractive, 3-1 invisible
				num_triangles++;

				tri_list[num_triangles*10 + 0] = R*cos(a1); // x1
				tri_list[num_triangles*10 + 1] = R*sin(a1); // y1
				tri_list[num_triangles*10 + 2] = 0.0; // z1
				tri_list[num_triangles*10 + 3] = R*cos(a2);  // x2
				tri_list[num_triangles*10 + 4] = R*sin(a2);  // y2
				tri_list[num_triangles*10 + 5] = 0.0; // z2
				tri_list[num_triangles*10 + 6] = R*cos(a2);  // x3
				tri_list[num_triangles*10 + 7] = R*sin(a2);  // y3
				tri_list[num_triangles*10 + 8] = L;  // z3
				tri_list[num_triangles*10 + 9] = 010004.0;   // exact 1, CSG 0, refractive, 3-1 invisible
				num_triangles++;
				}

			data[10] = num_triangles; /* how many we actually wrote out */
			return 0;
			}

		/* iterate to exact solution given starting point */
		case 2:
			{
			/*
			We are giving starting data, and must compute the distance to the actual surface.
			We also need to compute the normal vector at the point.

			This is simple iteration. All surfaces must be described as:
			F(x,y,z) = 0
			We will also need Fx = dF/dx, Fy = dF/dy, and Fz = dF/dz
			For the cylindrical face, we have set the triangles to exact = 1
			F  = R*R - X*X - Y*Y
			Fx = -2X
			Fy = -2Y
			Fz = 0.0

			then compute
			Fp = Fx*l + Fy*m + Fz*n where l, m, and n are the ray cosines

			the propagation distance is just
			delt = -F/Fp;
			keep a running total on the position
			t += delt;
							
			increment the ray coordinates
			x += l*delt;
			y += m*delt;
			z += n*delt;

			repeat until delt is small!

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
					/* cylinder face */
					{
					int loop;
					double x, y, z, l, m, n, t, delt, F, Fx, Fy, Fz, Fp;

					x = data[2];
					y = data[3];
					z = data[4];
					l = data[5];
					m = data[6];
					n = data[7];
					t = 0.0;
					delt = 100.0; /* any big number to start */
					loop = 0; /* loop is to prevent infinite loops, which will HANG ZEMAX! */

					while (fabs(delt) > 1E-10 && loop < 200)
						{
						loop++;

						F = R*R - x*x - y*y;
						Fx = -2.0*x;
						Fy = -2.0*y;
						Fz = 0.0; /* this could obviously be pulled from the loop for this surface but is included for clarity */
						Fp = Fx*l + Fy*m + Fz*n;

						delt = -F/Fp;
						t += delt;
						
						x += l*delt;
						y += m*delt;
						z += n*delt;
						}

					/* we have converged, hopefully */
					if (fabs(delt) > 1E-8)
						{
						/* assume ray misses */
						return -1;
						}

					/* normalize Fx, Fy, and Fz to get the normal vector */
					Fp = sqrt(Fx*Fx+Fy*Fy+Fz*Fz);
					if (Fp == 0.0) Fp = 1.0; /* this is not possible for a real surface, but better safe than sorry! */
					Fx /= Fp;
					Fy /= Fp;
					Fz /= Fp;
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
	if (i ==  0) strcpy(data,"Cylinder Face");
	if (i == -1) strcpy(data,"Bottom Flat");
	if (i == -2) strcpy(data,"Front Face");
	if (i == -3) strcpy(data,"Back Face");


	if (i == 1) strcpy(data,"Radius");
	if (i == 2) strcpy(data,"Length");
	if (i == 3) strcpy(data,"# Facets");
	return 0;
	}

