#define _USE_MATH_DEFINES
#include <windows.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "usersurf.h"
#include "slicer_generation.h"
#include "surface_solns.h"
#include "slice_param_helpers.h"

/*
Ellen Lee
Feb 2025
*/

// Maximum number of characters in any path string
#define MAX_PATH_LENGTH 512

int __declspec(dllexport) APIENTRY UserDefinedSurface5(USER_DATA *UD, FIXED_DATA5 *FD);

/* a generic Snells law refraction routine */
int Refract(double thisn, double nextn, double *l, double *m, double *n, double ln, double mn, double nn);

// Normally having global variables that can change is bad practice, but this is
// necessary for us to store the global extrema without having to recalculate them
// every time we trace rays. Because each analysis window in Zemax gets its own
// copy of the DLL, we shouldn't have to worry about locks or race conditions.
static double ZMIN, ZMAX, UMIN, UMAX;
static int FILE_NUM_OLD = -9999; // Store previous file number to check if it changed

// Keep the custom slice parameters in a global array so we don't need to reload
// the file every time this DLL is called.
static double *p_custom = NULL;
static char PARAMS_DIR[MAX_PATH_LENGTH];


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

         // Get the absolute directory of the current DLL file
         char dll_path[MAX_PATH_LENGTH];
         int ret;

         if (GetModuleFileNameA(hInst, dll_path, sizeof(dll_path)) == 0)
         {
            ret = GetLastError();
            MessageBoxA(NULL, "GetModuleHandle failed. Could not resolve the absolute path of the DLL file.", "Error", MB_OK);
            return FALSE;
         }

         // Find the last backslash to get the directory
         char *last_slash = strrchr(dll_path, '\\');
         if (last_slash) {
            *last_slash = '\0';  // Trim the DLL file name, keeping only the directory
         }
         snprintf(PARAMS_DIR, sizeof(PARAMS_DIR), "%s\\ifugen_params\\", dll_path);
         break;

      case DLL_PROCESS_DETACH:
         free(p_custom);
         break;
	}
	return TRUE;
}

int __declspec(dllexport) APIENTRY UserDefinedSurface5(USER_DATA *UD, FIXED_DATA5 *FD)
	{
   double sag, l_par, m_par, n_par;
   RAY_IN ray_in = {
      .xt = 0,
      .yt = 0,
      .l = 0,
      .m = 0,
      .n = 1
   };
   RAY_OUT ray_out = {
      .xs = NAN,
      .ys = NAN,
      .zs = NAN,
      .t = NAN,
      .ln = NAN,
      .mn = NAN,
      .nn = NAN
   };

   // Read in the custom slice array params
   GRID_PARAMS_BASIC p = MakeBasicParamsFromCustom(p_custom);
   ValidateBasicParams(&p);

   // Store FD params
   int trace_walls = FD->param[0];
   int active_x = FD->param[1];
   int active_y = FD->param[2];
   int file_num = FD->param[3];

   active_x = active_y = 0;

   switch(FD->type)
   	{
      case 0:
      	/* ZEMAX is requesting general information about the surface */
         switch(FD->numb)
         	{
            case 0:
            	/* ZEMAX wants to know the name of the surface */
		         /* do not exceed 12 characters */
		         strcpy(UD->string,"SlicerCustom");
               break;
            case 1:
            	/* ZEMAX wants to know if this surface is rotationally symmetric */
               /* it is not, so return a null string */
            	UD->string[0] = '\0';
               break;
            case 2:
            	/* ZEMAX wants to know if this surface is a gradient index media */
               /* it is not, so return a null string */
            	UD->string[0] = '\0';
            	break;
            }
         break;

      case 1:
      	/* ZEMAX is requesting the names of the parameter columns */
         /* the value FD->numb will indicate which value ZEMAX wants. */
         switch(FD->numb)
         	{
            case 0:
               strcpy(UD->string,"trace_walls");
               break;
            case 1:
               strcpy(UD->string,"active_x");
               break;
            case 2:
               strcpy(UD->string,"active_y");
               break;
            case 3:
               strcpy(UD->string,"file_num");
               break;
            default:
            	UD->string[0] = '\0';
            	break;
            }
      	break;

      case 2:
      	/*
         This case is now deprecated. It used to be necessary when inputs for a surface
         needed to be separately defined as "parameter" data and "extra" data. This is no
         longer the case, as OpticStudio now supports the ability to define all inputs as
         "parameter" data, as long as the FIXED_DATA5 structure is being used by the DLL.
         */
         break;

      case 3:
      	/* ZEMAX wants to know the sag of the surface */
         /* if there is an alternate sag, return it as well */
         /* otherwise, set the alternate sag identical to the sag */
         /* The sag is sag1, alternate is sag2. */
         UD->sag1 = 0.0;
         UD->sag2 = 0.0;

         sag = ImageSlicerSag(UD->x, UD->y, p, p_custom);

         if (isnan(sag)) return 0;    // Out of bounds, keep sag at 0
         else {
            UD->sag1 = sag;
            UD->sag2 = sag;
         }
         break;
      case 4:
      	/* ZEMAX wants a paraxial ray trace to this surface */
         /* x, y, z, and the path are unaffected, at least for this surface type */
         /* for paraxial ray tracing, the return z coordinate should always be zero. */
         /* paraxial surfaces are always planes with the following normals */
        ray_in.xt = (UD->x); ray_in.yt = (UD->y);
        ray_in.l = (UD->l); ray_in.m = (UD->m); ray_in.n = (UD->n);
        ParaxialRayTraceSlicer(&ray_out, &l_par, &m_par, &n_par, &ray_in,
            FD->n1, FD->n2, active_x, active_y, p, p_custom);

        if (isnan(ray_out.t)) return (FD->surf);  // Missed somehow?

        UD->l = l_par,
        UD->m = m_par;
        UD->n = n_par;
        UD->ln = ray_out.ln;
        UD->mn = ray_out.mn;
        UD->nn = ray_out.nn;
         break;

      case 5:
      	/* ZEMAX wants a real ray trace to this surface */
         ray_in.xt = (UD->x); ray_in.yt = (UD->y);
         ray_in.l = (UD->l); ray_in.m = (UD->m); ray_in.n = (UD->n);

         // The first thing this function does is reset the members of ray_out
         // to be all NANs
         RayTraceSlicer(&ray_out, ray_in, ZMIN, ZMAX, UMIN, UMAX, trace_walls,
            p, p_custom);

         // Ray missed if transfer distance or normal vector could not be computed
         if (isnan(ray_out.t) || isnan(ray_out.ln)) return (FD->surf);

         (UD->x) = ray_out.xs; (UD->y) = ray_out.ys; (UD->z) = ray_out.zs;
         (UD->path) = ray_out.t;
         (UD->ln) = ray_out.ln; (UD->mn) = ray_out.mn; (UD->nn) = ray_out.nn;

         if (Refract(FD->n1, FD->n2, &UD->l, &UD->m, &UD->n, UD->ln, UD->mn, UD->nn)) return(-FD->surf);

         break;

      case 6:
      	/* ZEMAX wants the index, dn/dx, dn/dy, and dn/dz at the given x, y, z. */

         /* This is only required for gradient index surfaces, so return dummy values */
         UD->index = FD->n2;
         UD->dndx = 0.0;
         UD->dndy = 0.0;
         UD->dndz = 0.0;
      	break;

      case 7:
      	/* ZEMAX wants the "safe" data. */
         /* this is used by ZEMAX to set the initial values for all parameters and extra data */
         /* when the user first changes to this surface type. */
         /* this is the only time the DLL should modify the data in the FIXED_DATA FD structure */
         FD->param[0] = 0;      // trace_walls
         FD->param[1] = 0;      // active_x
         FD->param[2] = 0;      // active_y
         FD->param[3] = 0;      // file_num
         FD->cv = 0;
         FD->k = 0;

         break;

      case 8:
         /* ZEMAX is calling the DLL for the first time, do any memory or data initialization here. */
          p = MakeBasicParamsFromCustom(p_custom);
          ValidateBasicParams(&p);

         if ( FILE_NUM_OLD != file_num ) {
            // Update custom slice params and global extrema
            LoadCustomParamsFromFile(p_custom, file_num, PARAMS_DIR);
            p = MakeBasicParamsFromCustom(p_custom);
            FindSlicerGlobalExtrema(&ZMIN, &ZMAX, p, p_custom);
            GetMinMaxU(&UMIN, &UMAX, p, p_custom);
            FILE_NUM_OLD = file_num;
         };

         ValidateBasicParams(&p);
         break;

      case 9:
      	/* ZEMAX is calling the DLL for the last time, do any memory release here. */
         break;

      case 10:
         /* Scaling of parameter and extra data values */
         break;
      }
   return 0;
   }

int Refract(double thisn, double nextn, double *l, double *m, double *n, double ln, double mn, double nn)
{
double nr, cosi, cosi2, rad, cosr, gamma;
if (thisn != nextn)
	{
	nr = thisn / nextn;
	cosi = fabs((*l) * ln + (*m) * mn + (*n) * nn);
	cosi2 = cosi * cosi;
	if (cosi2 > 1) cosi2 = 1;
	rad = 1 - ((1 - cosi2) * (nr * nr));
	if (rad < 0) return(-1);
	cosr = sqrt(rad);
	gamma = nr * cosi - cosr;
	(*l) = (nr * (*l)) + (gamma * ln);
	(*m) = (nr * (*m)) + (gamma * mn);
	(*n) = (nr * (*n)) + (gamma * nn);
	}
return 0;
}
