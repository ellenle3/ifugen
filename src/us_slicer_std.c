#define _USE_MATH_DEFINES
#include <windows.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "usersurf.h"
#include "slicer_generation.h"
#include "surface_solns.h"

/*
Ellen Lee
Feb 2025
*/


int __declspec(dllexport) APIENTRY UserDefinedSurface5(USER_DATA *UD, FIXED_DATA5 *FD);

/* a generic Snells law refraction routine */
int Refract(double thisn, double nextn, double *l, double *m, double *n, double ln, double mn, double nn);
void SetFDFromSlicerParams(IMAGE_SLICER_PARAMS *p, FIXED_DATA5 *FD);
void SetSlicerParamsFromFD(IMAGE_SLICER_PARAMS *p, FIXED_DATA5 *FD);

// Normally having global variables that can change is bad practice, but this is
// necessary for us to store the global extrema without having to recalculate them
// every time we trace rays. Because each analysis window in Zemax gets its own
// copy of the DLL, we shouldn't have to worry about locks or race conditions.
static double ZMIN, ZMAX;
// Need to keep track of whether parameters changed
static IMAGE_SLICER_PARAMS P_OLD = {
        .custom = 0,
        .surface_type = -1,
        .n_each = -1,
        .n_rows = -1,
        .n_cols = -1,
        .angle_mode = -1,
        .dalpha = -1,
        .dbeta = -1,
        .dgamma = -1,
        .alpha_cen = -1,
        .beta_cen = -1,
        .gamma_cen = -1,
        .gamma_offset = -1,
        .dx = -1,
        .dy = -1,
        .gx_width = -1,
        .gx_depth = -1,
        .gy_width = -1,
        .gy_depth = -1,
        .cv = -1,
        .k = -1
    };


BOOL WINAPI DllMain (HANDLE hInst, ULONG ul_reason_for_call, LPVOID lpReserved)
	{
   return TRUE;
   }

int __declspec(dllexport) APIENTRY UserDefinedSurface5(USER_DATA *UD, FIXED_DATA5 *FD)
	{
   double power, sag;
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

   IMAGE_SLICER_PARAMS p; SetSlicerParamsFromFD(&p, FD);
   // We will never need the custom_slice_params array in standard mode but we need
   // to pass it as an argument. It will never be accessed because p.custom will
   // always be set to 0.
   double custom_slice_params[1] = {0};
   p.custom = 0;

   // Store some other FD params
   int trace_walls = FD->param[0];
   int active_x = FD->param[1];
   int active_y = FD->param[2];

   switch(FD->type)
   	{
      case 0:
      	/* ZEMAX is requesting general information about the surface */
         switch(FD->numb)
         	{
            case 0:
            	/* ZEMAX wants to know the name of the surface */
		         /* do not exceed 12 characters */
		         strcpy(UD->string,"SlicerStd");
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
               strcpy(UD->string,"surface_type");
               break;
            case 4:
               strcpy(UD->string,"n_each");
               break;
            case 5:
               strcpy(UD->string,"n_rows");
               break;
            case 6:
               strcpy(UD->string,"n_cols");
               break;
            case 7:
            	strcpy(UD->string,"angle_mode");
               break;
            case 8:
            	strcpy(UD->string,"dalpha");
               break;
            case 9:
            	strcpy(UD->string,"dbeta");
               break;
            case 10:
            	strcpy(UD->string,"dgamma");
               break;
            case 11:
               strcpy(UD->string,"gamma_offset");
               break;
            case 12:
            	strcpy(UD->string,"alpha_cen");
               break;
            case 13:
            	strcpy(UD->string,"beta_cen");
               break;
            case 14:
            	strcpy(UD->string,"gamma_cen");
               break;
            case 15:
            	strcpy(UD->string,"dx");
               break;
            case 16:
            	strcpy(UD->string,"dy");
               break;
            case 17:
            	strcpy(UD->string,"gx_width");
               break;
            case 18:
            	strcpy(UD->string,"gx_depth");
               break;
            case 19:
            	strcpy(UD->string,"gy_width");
               break;
            case 20:
            	strcpy(UD->string,"gy_depth");
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

         sag = ImageSlicerSag(UD->x, UD->y, p, custom_slice_params);

         if (isnan(sag)) return 0;    // Out of bounds, keep sag at 0... should I return -1?
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
         
         UD->ln =  0.0;
         UD->mn =  0.0;
         UD->nn = -1.0;
         power = (FD->n2 - FD->n1)*FD->cv;
         // 
         // if ((UD->n) != 0.0)
         // 	{
         //    (UD->l) = (UD->l)/(UD->n);
         //    (UD->m) = (UD->m)/(UD->n);

         //    (UD->l) = (FD->n1*(UD->l) - (UD->x)*power)/(FD->n2);
         //    (UD->m) = (FD->n1*(UD->m) - (UD->y)*power)/(FD->n2);

         //    /* normalize */
         //    (UD->n) = sqrt(1/(1 + (UD->l)*(UD->l) + (UD->m)*(UD->m) ) );
         //    /* de-paraxialize */
         //    (UD->l) = (UD->l)*(UD->n);
         //    (UD->m) = (UD->m)*(UD->n);
         //    }
         // break;

      case 5:
      	/* ZEMAX wants a real ray trace to this surface */
         ray_in.xt = (UD->x); ray_in.yt = (UD->y);
         ray_in.l = (UD->l); ray_in.m = (UD->m); ray_in.n = (UD->n);

         RayTraceSlicer(&ray_out, ray_in, ZMIN, ZMAX, trace_walls,
            p, custom_slice_params);

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
         FD->param[3] = 0;      // surface_type
         FD->param[4] = 5;      // n_each
         FD->param[5] = 2;      // n_rows
         FD->param[6] = 1;      // n_cols
         FD->param[7] = 0;      // angle_mode
         FD->param[8] = 4.0;    // dalpha
         FD->param[9] = 4.0;    // dbeta
         FD->param[10] = 1.0;   // dgamma
         FD->param[11] = 0.0;   // gamma_offset
         FD->param[12] = 0.0;   // alpha_cen
         FD->param[13] = 0.0;   // beta_cen
         FD->param[14] = 0.0;   // gamma_cen
         FD->param[15] = 8.0;   // dx
         FD->param[16] = 1.5;   // dy
         FD->param[17] = 0.0;   // gx_width
         FD->param[18] = 0.0;   // gx_depth
         FD->param[19] = 0.0;   // gy_width
         FD->param[20] = 0.0;   // gy_depth
         FD->cv = -0.01;
         FD->k = 0;

         SetSlicerParamsFromFD(&p, FD);
         break;

      case 8:
         /* ZEMAX is calling the DLL for the first time, do any memory or data initialization here. */

         //if ( ValidateSlicerParams(&p) ) SetFDFromSlicerParams(&p, FD); // prevent illegal values
         ValidateSlicerParams(&p);

         if ( !IsParametersEqual(p, P_OLD) ) {
            // Update global extrema
            FindSlicerGlobalExtrema(&ZMIN, &ZMAX, p, custom_slice_params);
            P_OLD = p;
         };
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

// Functions to convert between Zemax FIXED_DATA and image slicer params struct
void SetFDFromSlicerParams(IMAGE_SLICER_PARAMS *p, FIXED_DATA5 *FD) {
   FD->param[3] =  p->surface_type;
   FD->param[4] =  p->n_each;
   FD->param[5] =  p->n_rows;
   FD->param[6] =  p->n_cols;
   FD->param[7] =  p->angle_mode;
   FD->param[8] =  p->dalpha;
   FD->param[9] =  p->dbeta;
   FD->param[10] = p->dgamma;
   FD->param[11] = p->gamma_offset;
   FD->param[12] = p->alpha_cen;
   FD->param[13] = p->beta_cen;
   FD->param[14] = p->gamma_cen;
   FD->param[15] = p->dx;
   FD->param[16] = p->dy;
   FD->param[17] = p->gx_width;
   FD->param[18] = p->gx_depth;
   FD->param[19] = p->gy_width;
   FD->param[20] = p->gy_depth;
   FD->cv = p->cv;
   FD->k = p->k;
}

void SetSlicerParamsFromFD(IMAGE_SLICER_PARAMS *p, FIXED_DATA5 *FD) {
   p->surface_type =     FD->param[3];
   p->n_each =       FD->param[4];
   p->n_rows =       FD->param[5];
   p->n_cols =       FD->param[6];
   p->angle_mode =   FD->param[7];
   p->dalpha =       FD->param[8];
   p->dbeta =        FD->param[9];
   p->dgamma =       FD->param[10];
   p->gamma_offset = FD->param[11];
   p->alpha_cen =    FD->param[12];
   p->beta_cen =     FD->param[13];
   p->gamma_cen =    FD->param[14];
   p->dx =           FD->param[15];
   p->dy =           FD->param[16];
   p->gx_width =     FD->param[17];
   p->gx_depth =     FD->param[18];
   p->gy_width =     FD->param[19];
   p->gy_depth =     FD->param[20];
   p->cv = FD->cv;
   p->k =  FD->k;
}