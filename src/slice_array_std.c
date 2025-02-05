#include <windows.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "usersurf.h"
#include "ifu_helpers.h"

/*
Ellen Lee
*/

int __declspec(dllexport) APIENTRY UserDefinedSurface5(USER_DATA *UD, FIXED_DATA5 *FD);

/* a generic Snells law refraction routine */
int Refract(double thisn, double nextn, double *l, double *m, double *n, double ln, double mn, double nn);

BOOL WINAPI DllMain (HANDLE hInst, ULONG ul_reason_for_call, LPVOID lpReserved)
	{
   return TRUE;
   }

int __declspec(dllexport) APIENTRY UserDefinedSurface5(USER_DATA *UD, FIXED_DATA5 *FD)
	{
   int i;
   int n_each, n_rows, mode;
   IMAGE_SLICER_PARAMS p;
   double sag, zmax, zmin;

   switch(FD->type)
   	{
      case 0:
      	/* ZEMAX is requesting general information about the surface */
         switch(FD->numb)
         	{
            case 0:
            	/* ZEMAX wants to know the name of the surface */
		         /* do not exceed 12 characters */
		         strcpy(UD->string,"SliceArrStd");
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
            default:
            	UD->string[0] = '\0';
            	break;
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
            	strcpy(UD->string,"trace_gaps");
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
            	strcpy(UD->string,"gax_width");
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
            }
      	break;
      case 2:
      	/* ZEMAX is requesting the names of the extra data columns */
         /* the value FD->numb will indicate which value ZEMAX wants. */
         /* they are all "Unused" for this surface type */
         /* returning a null string indicates that the extradata value is unused. */
         switch(FD->numb)
         	{
            default:
            	UD->string[0] = '\0';
            	break;
            }
      	break;
      case 3:
      	/* ZEMAX wants to know the sag of the surface */
         /* if there is an alternate sag, return it as well */
         /* otherwise, set the alternate sag identical to the sag */
         /* The sag is sag1, alternate is sag2. */
         UD->sag1 = 0.0;
         UD->sag2 = 0.0;

         x = UD->x;
         y = UD->y;
         sag = 0.0;
         ImageSlicerSag(&sag, x, y, p);

         // Set the sag to zero if it isn't defined at the given x, y. This is
         // only used to draw the surface so it shouldn't matter. We will tell Zemax
         // that the ray missed the surface when it is traced.
         if (isnan(sag)) {
            UD->sag1 = 0.0;
            UD->sag2 = 0.0;
         }
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
         p.n_each = 5;
         p.n_rows = 3;
         p.n_cols = 2;
         p.mode = 0;
         p.trace_walls = 0;
         p.active_x = 0;
         p.active_y = 0;
         p.dalpha = 4;
         p.dbeta = 4;
         p.dgamma = 0.2;
         p.alpha_cen = 0;
         p.beta_cen = 0;
         p.gamma_cen = 0;
         p.dx = 10;
         p.dy = 0.5;
         p.gx_width = 0;
         p.gx_depth = 0;
         p.gy_width = 0;
         p.gy_depth = 0;
         break;
      case 8:
      	/* ZEMAX is calling the DLL for the first time, do any memory or data initialization here. */
         // Initialize the parameter struct
         p.n_each =     FD->param[0];
         p.n_rows =     FD->param[1];
         p.n_cols =     FD->param[2];
         p.mode =       FD->param[3];
         p.trace_gaps = FD->param[4];
         p.active_x =   FD->param[5];
         p.active_y =   FD->param[6];
         p.dalpha =     FD->param[7];
         p.dbeta =      FD->param[8];
         p.dgamma =     FD->param[9];
         p.alpha_cen =  FD->param[10];
         p.beta_cen =   FD->param[11];
         p.gamma_cen =  FD->param[12];
         p.dx =         FD->param[13];
         p.dy =         FD->param[14];
         p.gx_width =   FD->param[15];
         p.gx_depth =   FD->param[16];
         p.gy_width =   FD->param[17];
         p.gy_depth =   FD->param[18];
         p.cv = FD->cv;
         p.k = FD->k;
         // Compute and store the global maxima
         //zmax=0; zmin=0;
         break;
      case 9:
      	/* ZEMAX is calling the DLL for the last time, do any memory release here. */
         break;
      case 10;
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