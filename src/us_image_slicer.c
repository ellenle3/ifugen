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

int __declspec(dllexport) APIENTRY UserDefinedSurface5(USER_DATA *UD, FIXED_DATA3 *FD);

/* a generic Snells law refraction routine */
int Refract(double thisn, double nextn, double *l, double *m, double *n, double ln, double mn, double nn);

BOOL WINAPI DllMain (HANDLE hInst, ULONG ul_reason_for_call, LPVOID lpReserved)
	{
   return TRUE;
   }

int __declspec(dllexport) APIENTRY UserDefinedSurface5(USER_DATA *UD, FIXED_DATA3 *FD)
	{
   int i;
   int n_each, n_rows, mode;
   imageSlicerParams p;
   double sag;
   double* zptr;

   switch(FD->type)
   	{
      case 0:
      	/* ZEMAX is requesting general information about the surface */
         switch(FD->numb)
         	{
            case 0:
            	/* ZEMAX wants to know the name of the surface */
		         /* do not exceed 12 characters */
		         strcpy(UD->string,"Image Slicer");
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
            case 1:
            	strcpy(UD->string,"dx");
               break;
            case 2:
            	strcpy(UD->string,"dy");
               break;
            case 3:
            	strcpy(UD->string,"Rx");
               break;
            case 4:
            	strcpy(UD->string,"Ry");
               break;
            case 5:
            	strcpy(UD->string,"kx");
               break;
            case 6:
            	strcpy(UD->string,"ky");
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

         c = FD->c;
         k = FD->k;
         n_each = FD->param[1];
         n_rows = FD->param[2];
         mode = FD->param[3];
         dx = FD->param[4];
         dy = FD->param[5];
         gap_width = FD->param[6];
         gap_depth =FD->param[7];

         x = UD->x;
         y = UD->y;
         sag = 0.0;
         ImageSlicerSag(zptr, x, y, n_each, n_rows, mode, dalpha, dbeta, alpha_cen, beta_cen, dx, dy, c, k, gap_width, gap_depth);
         double sag = *zptr;
         UD->sag1 = sag;
         UD->sag2 = sag;
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
         for (i = 0; i <= FD->max_parameter; i++) FD->param[i] = 0.0;
         for (i = 0; i <= FD->max_extradata; i++) FD->xdata[i] = 0.0;
         break;
      case 8:
      	/* ZEMAX is calling the DLL for the first time, do any memory or data initialization here. */
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