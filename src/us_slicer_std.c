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
static double ZMIN, ZMAX, UMIN, UMAX;
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
        .gamma_offset = -1,
        .dzps = -1,
        .dzp_col = -1,
        .dzp_row = -1,
        .dsyx = -1,
        .dsyz = -1,
        .dsxy = -1,
        .dsxz = -1,
        .du = -1,
        .alpha_cen = -1,
        .beta_cen = -1,
        .gamma_cen = -1,
        .zps_cen = -1,
        .zp_cen = -1,
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


BOOL WINAPI DllMain (HANDLE hInst, ULONG ul_reason_for_call, LPVOID lpReserved)
	{
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

   IMAGE_SLICER_PARAMS p; SetSlicerParamsFromFD(&p, FD);
   // We will never need the p_custom array in standard mode but we need
   // to pass it as an argument. It will never be accessed because p.custom will
   // always be set to 0.
   double p_custom[1] = {0};
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
            	strcpy(UD->string,"dzps");
               break;
            case 13:
            	strcpy(UD->string,"dzp_col");
               break;
            case 14:
            	strcpy(UD->string,"dzp_row");
               break;
            case 15:
            	strcpy(UD->string,"dsyx");
               break;
            case 16:
            	strcpy(UD->string,"dsyz");
               break;
            case 17:
            	strcpy(UD->string,"dsxy");
               break;
            case 18:
            	strcpy(UD->string,"dsxz");
               break;
            case 19:
            	strcpy(UD->string,"du");
               break;
            case 20:
            	strcpy(UD->string,"alpha_cen");
               break;
            case 21:
            	strcpy(UD->string,"beta_cen");
               break;
            case 22:
            	strcpy(UD->string,"gamma_cen");
               break;
            case 23:
            	strcpy(UD->string,"zps_cen");
               break;
            case 24:
            	strcpy(UD->string,"zp_cen");
               break;
            case 25:
            	strcpy(UD->string,"syx_cen");
               break;
            case 26:
            	strcpy(UD->string,"syz_cen");
               break;
            case 27:
            	strcpy(UD->string,"sxy_cen");
               break;
            case 28:
            	strcpy(UD->string,"sxz_cen");
               break;
            case 29:
            	strcpy(UD->string,"u_cen");
               break;
            case 30:
            	strcpy(UD->string,"dx");
               break;
            case 31:
            	strcpy(UD->string,"dy");
               break;
            case 32:
            	strcpy(UD->string,"gx_width");
               break;
            case 33:
            	strcpy(UD->string,"gx_depth");
               break;
            case 34:
            	strcpy(UD->string,"gy_width");
               break;
            case 35:
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

         sag = ImageSlicerSag(UD->x, UD->y, p, p_custom);

         if (isnan(sag)) return 0;    // Out of bounds
         else {
            UD->sag1 = sag;
            UD->sag2 = sag;
         }
         break;
      case 4:
      	/* ZEMAX wants a paraxial ray trace to this surface */
         /* x, y, z, and the path are unaffected, at least for this surface type */
         /* for paraxial ray tracing, the return z coordinate should always be zero. */
         ray_in.xt = (UD->x); ray_in.yt = (UD->y);
         ray_in.l = (UD->l); ray_in.m = (UD->m); ray_in.n = (UD->n);
         ParaxialRayTraceSlicer(&ray_out, ray_in, &l_par, &m_par, &n_par,
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
         FD->param[3] = 0;      // surface_type
         FD->param[4] = 5;      // n_each
         FD->param[5] = 2;      // n_rows
         FD->param[6] = 1;      // n_cols
         FD->param[7] = 0;      // angle_mode
         FD->param[8] = 4.0;    // dalpha
         FD->param[9] = 4.0;    // dbeta
         FD->param[10] = 1.0;   // dgamma
         FD->param[11] = 0.0;   // gamma_offset
         FD->param[12] = 0.0;   // dzps
         FD->param[13] = 0.0;   // dzp_col
         FD->param[14] = 0.0;   // dzp_row
         FD->param[15] = 0.0;   // dsyx
         FD->param[16] = 0.0;   // dsyz
         FD->param[17] = 0.0;   // dsxy
         FD->param[18] = 0.0;   // dsxz
         FD->param[19] = 0.0;   // du
         FD->param[20] = 0.0;   // alpha_cen
         FD->param[21] = 0.0;   // beta_cen
         FD->param[22] = 0.0;   // gamma_cen
         FD->param[23] = 0.0;   // zps_cen
         FD->param[24] = 0.0;   // zp_cen
         FD->param[25] = 0.0;   // syx_cen
         FD->param[26] = 0.0;   // syz_cen
         FD->param[27] = 0.0;   // sxy_cen
         FD->param[28] = 0.0;   // sxz_cen
         FD->param[29] = 0.0;   // u_cen
         FD->param[30] = 8.0;   // dx
         FD->param[31] = 1.5;   // dy
         FD->param[32] = 0.0;   // gx_width
         FD->param[33] = 0.0;   // gx_depth
         FD->param[34] = 0.0;   // gy_width
         FD->param[35] = 0.0;   // gy_depth
         FD->cv = -0.01;
         FD->k = 0;

         SetSlicerParamsFromFD(&p, FD);
         break;

      case 8:
         /* ZEMAX is calling the DLL for the first time, do any memory or data initialization here. */
         ValidateSlicerParams(&p); // prevent illegal values - not allowed to modify FD here...

         if ( !IsParametersEqual(p, P_OLD) ) {
            // Update global extrema
            FindSlicerGlobalExtrema(&ZMIN, &ZMAX, p, p_custom);
            GetMinMaxU(&UMIN, &UMAX, p, p_custom);
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
   FD->param[12] = p->dzps;
   FD->param[13] = p->dzp_col;
   FD->param[14] = p->dzp_row;
   FD->param[15] = p->dsyx;
   FD->param[16] = p->dsyz;
   FD->param[17] = p->dsxy;
   FD->param[18] = p->dsxz;
   FD->param[19] = p->du;
   FD->param[20] = p->alpha_cen;
   FD->param[21] = p->beta_cen;
   FD->param[22] = p->gamma_cen;
   FD->param[23] = p->zps_cen;
   FD->param[24] = p->zp_cen;
   FD->param[25] = p->syx_cen;
   FD->param[26] = p->syz_cen;
   FD->param[27] = p->sxy_cen;
   FD->param[28] = p->sxz_cen;
   FD->param[29] = p->u_cen;
   FD->param[30] = p->dx;
   FD->param[31] = p->dy;
   FD->param[32] = p->gx_width;
   FD->param[33] = p->gx_depth;
   FD->param[34] = p->gy_width;
   FD->param[35] = p->gy_depth;
   FD->cv = p->cv;
   FD->k = p->k;
}

void SetSlicerParamsFromFD(IMAGE_SLICER_PARAMS *p, FIXED_DATA5 *FD) {
   p->surface_type = FD->param[3];
   p->n_each =       FD->param[4];
   p->n_rows =       FD->param[5];
   p->n_cols =       FD->param[6];
   p->angle_mode =   FD->param[7];
   p->dalpha =       FD->param[8];
   p->dbeta =        FD->param[9];
   p->dgamma =       FD->param[10];
   p->gamma_offset = FD->param[11];
   p->dzps =         FD->param[12];
   p->dzp_col =      FD->param[13];
   p->dzp_row =      FD->param[14];
   p->dsyx =         FD->param[15];
   p->dsyz =         FD->param[16];
   p->dsxy =         FD->param[17];
   p->dsxz =         FD->param[18];
   p->du =           FD->param[19];
   p->alpha_cen =    FD->param[20];
   p->beta_cen =     FD->param[21];
   p->gamma_cen =    FD->param[22];
   p->zps_cen =      FD->param[23];
   p->zp_cen =       FD->param[24];
   p->syx_cen =      FD->param[25];
   p->syz_cen =      FD->param[26];
   p->sxy_cen =      FD->param[27];
   p->sxz_cen =      FD->param[28];
   p->u_cen =        FD->param[29];
   p->dx =           FD->param[30];
   p->dy =           FD->param[31];
   p->gx_width =     FD->param[32];
   p->gx_depth =     FD->param[33];
   p->gy_width =     FD->param[34];
   p->gy_depth =     FD->param[35];
   p->cv = FD->cv;
   p->k =  FD->k;
}