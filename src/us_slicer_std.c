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
*/

int __declspec(dllexport) APIENTRY UserDefinedSurface5(USER_DATA *UD, FIXED_DATA5 *FD);

/* a generic Snells law refraction routine */
int Refract(double thisn, double nextn, double *l, double *m, double *n, double ln, double mn, double nn);
void SetSlicerParamsFromFD(IMAGE_SLICER_PARAMS *p, FIXED_DATA5 *FD);

BOOL WINAPI DllMain (HANDLE hInst, ULONG ul_reason_for_call, LPVOID lpReserved)
	{
   return TRUE;
   }

int __declspec(dllexport) APIENTRY UserDefinedSurface5(USER_DATA *UD, FIXED_DATA5 *FD)
	{
   //FILE* fptr;
   int i;
   double p2, alpha, a, b, c, rad, casp, t, zc;
   double power;
   double zmax = 1E13, zmin = -1E13, sag = 0.0;
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

   IMAGE_SLICER_PARAMS p = {
      .n_each = 5,
      .n_rows = 1,
      .n_cols = 1,
      .mode = 0,
      .trace_walls = 0,
      .active_x = 0,
      .active_y = 0,
      .dalpha = 4,
      .dbeta = 4,
      .dgamma = 1,
      .alpha_cen = 0,
      .beta_cen = 0,
      .gamma_cen = 0,
      .dx = 10,
      .dy = 2,
      .gx_width = 0,
      .gx_depth = 0,
      .gy_width = 0,
      .gy_depth = 0,
      .cv = -0.01,
      .k = 0
   };
   SAG_FUNC sag_func = &Conic3DSag;
   TRANSFER_DIST_FUNC transfer_dist_func = &Conic3DTransfer;
   CRITICAL_XY_FUNC critical_xy_func = &Conic3DCriticalXY;
   SURF_NORMAL_FUNC surf_normal_func = &Conic3DSurfaceNormal;

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
            	strcpy(UD->string,"n_each");
               break;
            case 1:
            	strcpy(UD->string,"n_rows");
               break;
            case 2:
            	strcpy(UD->string,"n_cols");
               break;
            case 3:
            	strcpy(UD->string,"mode");
               break;
            case 4:
            	strcpy(UD->string,"trace_walls");
               break;
            case 5:
               strcpy(UD->string,"active_x");
               break;
            case 6:
               strcpy(UD->string,"active_y");
               break;
            case 7:
            	strcpy(UD->string,"dalpha");
               break;
            case 8:
            	strcpy(UD->string,"dbeta");
               break;
            case 9:
            	strcpy(UD->string,"dgamma");
               break;
            case 10:
            	strcpy(UD->string,"alpha_cen");
               break;
            case 11:
            	strcpy(UD->string,"beta_cen");
               break;
            case 12:
            	strcpy(UD->string,"gamma_cen");
               break;
            case 13:
            	strcpy(UD->string,"dx");
               break;
            case 14:
            	strcpy(UD->string,"dy");
               break;
            case 15:
            	strcpy(UD->string,"gx_width");
               break;
            case 16:
            	strcpy(UD->string,"gx_depth");
               break;
            case 17:
            	strcpy(UD->string,"gy_width");
               break;
            case 18:
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

         sag = ImageSlicerSag(UD->x, UD->y, p, sag_func);
         //fptr = fopen("C:\\Projects\\ifugen\\test_output.txt", "a");
         //fprintf(fptr, "%.10f %.10f %.10f\n", UD->x, UD->y, sag);
         //fclose(fptr);

         if (isnan(sag)) return -1;    // Out of bounds
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

         // Use active_x and active_y to get the central slice if there are an
         // odd number of slices

         // Calculate transfer distance based on these parameters

         UD->ln =  0.0;
         UD->mn =  0.0;
         UD->nn = -1.0;
         power = (FD->n2 - FD->n1)*FD->cv;
         if ((UD->n) != 0.0)
         	{
            (UD->l) = (UD->l)/(UD->n);
            (UD->m) = (UD->m)/(UD->n);

            (UD->l) = (FD->n1*(UD->l) - (UD->x)*power)/(FD->n2);
            (UD->m) = (FD->n1*(UD->m) - (UD->y)*power)/(FD->n2);

            /* normalize */
            (UD->n) = sqrt(1/(1 + (UD->l)*(UD->l) + (UD->m)*(UD->m) ) );
            /* de-paraxialize */
            (UD->l) = (UD->l)*(UD->n);
            (UD->m) = (UD->m)*(UD->n);
            }
         break;

      case 5:
      	/* ZEMAX wants a real ray trace to this surface */
         // ray_in.xt = (UD->x); ray_in.yt = (UD->y);
         // ray_in.l = (UD->l); ray_in.m = (UD->m); ray_in.n = (UD->n);

         // // The first thing this function does is reset the members of ray_out
         // // to be all NANs
         // RayTraceSlicer(&ray_out, ray_in, zmin, zmax,
         //    p, sag_func, transfer_dist_func, surf_normal_func);

         // // Ray missed if transfer distance or normal vector could not be computed
         // if (isnan(ray_out.t) || isnan(ray_out.ln)) return (FD->surf);

         // (UD->x) = ray_out.xs; (UD->y) = ray_out.ys; (UD->z) = ray_out.zs;
         // (UD->path) = ray_out.t;
         // (UD->ln) = ray_out.ln; (UD->mn) = ray_out.mn; (UD->nn) = ray_out.nn;

         // if (Refract(FD->n1, FD->n2, &UD->l, &UD->m, &UD->n, UD->ln, UD->mn, UD->nn)) return(-FD->surf);

         /* ZEMAX wants a real ray trace to this surface */
         if (FD->cv == 0.0)
         	{
	         UD->ln =  0.0;
   	      UD->mn =  0.0;
      	   UD->nn = -1.0;
			   if (Refract(FD->n1, FD->n2, &UD->l, &UD->m, &UD->n, UD->ln, UD->mn, UD->nn)) return(-FD->surf);
            return(0);
            }
         /* okay, not a plane. */
			a = (UD->n) * (UD->n) * FD->k + 1;
			b = ((UD->n)/FD->cv) - (UD->x) * (UD->l) - (UD->y) * (UD->m);
			c = (UD->x) * (UD->x) + (UD->y) * (UD->y);
			rad = b * b - a * c;
			if (rad < 0) return(FD->surf);  /* ray missed this surface */
			if (FD->cv > 0) t = c / (b + sqrt(rad));
			else           t = c / (b - sqrt(rad));
			(UD->x) = (UD->l) * t + (UD->x);
			(UD->y) = (UD->m) * t + (UD->y);
			(UD->z) = (UD->n) * t + (UD->z);
			UD->path = t;
			zc = (UD->z) * FD->cv;
			rad = zc * FD->k * (zc * (FD->k + 1) - 2) + 1;
			casp = FD->cv / sqrt(rad);
			UD->ln = (UD->x) * casp;
			UD->mn = (UD->y) * casp;
			UD->nn = ((UD->z) - ((1/FD->cv) - (UD->z) * FD->k)) * casp;
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
         FD->param[0] = 5;      // n_each
         FD->param[1] = 1;      // n_rows
         FD->param[2] = 1;      // n_cols
         FD->param[3] = 0;      // mode
         FD->param[4] = 0;      // trace_walls
         FD->param[5] = 0;      // active_x
         FD->param[6] = 4.0;    // active_y
         FD->param[7] = 4.0;    // dalpha
         FD->param[8] = 4.0;    // dbeta
         FD->param[9] = 1.0;    // dgamma
         FD->param[10] = 0.0;   // alpha_cen
         FD->param[11] = 0.0;   // beta_cen
         FD->param[12] = 0.0;   // gamma_cen
         FD->param[13] = 10.0;  // dx
         FD->param[14] = 2;     // dy
         FD->param[15] = 0.0;   // gx_width
         FD->param[16] = 0.0;   // gx_depth
         FD->param[17] = 0.0;   // gy_width
         FD->param[18] = 0.0;   // gy_depth
         FD->cv = -0.01;
         FD->k = 0;

         SetSlicerParamsFromFD(&p, FD);
         break;

      case 8:
      	/* ZEMAX is calling the DLL for the first time, do any memory or data initialization here. */
         // Initialize the parameter struct
         SetSlicerParamsFromFD(&p, FD);

         // Validate the parameters. Technically we shouldn't change the 
         //CheckSlicerParams(&p);

         // If plane, use different solutions
         if (p.cv == 0) {
            sag_func = &TiltedPlaneSag;
            critical_xy_func = &TiltedPlaneCriticalXY;
            transfer_dist_func = &TiltedPlaneTransfer;
            surf_normal_func = &TiltedPlaneSurfaceNormal;
         }

         // Compute and store the global extrema
         FindSlicerGlobalExtrema(&zmin, &zmax, p, sag_func, critical_xy_func);
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

void SetSlicerParamsFromFD(IMAGE_SLICER_PARAMS *p, FIXED_DATA5 *FD) {
    p->n_each =      FD->param[0];
    p->n_rows =      FD->param[1];
    p->n_cols =      FD->param[2];
    p->mode =        FD->param[3];
    p->trace_walls = FD->param[4];
    p->active_x =    FD->param[5];
    p->active_y =    FD->param[6];
    p->dalpha =      FD->param[7];
    p->dbeta =       FD->param[8];
    p->dgamma =      FD->param[9];
    p->alpha_cen =   FD->param[10];
    p->beta_cen =    FD->param[11];
    p->gamma_cen =   FD->param[12];
    p->dx =          FD->param[13];
    p->dy =          FD->param[14];
    p->gx_width =    FD->param[15];
    p->gx_depth =    FD->param[16];
    p->gy_width =    FD->param[17];
    p->gy_depth =    FD->param[18];
    p->cv = FD->cv;
    p->k =  FD->k;
}