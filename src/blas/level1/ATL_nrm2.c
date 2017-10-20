/*
 *             Automatically Tuned Linear Algebra Software v3.10.3
 *                    (C) Copyright 1999 R. Clint Whaley
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *   1. Redistributions of source code must retain the above copyright
 *      notice, this list of conditions and the following disclaimer.
 *   2. Redistributions in binary form must reproduce the above copyright
 *      notice, this list of conditions, and the following disclaimer in the
 *      documentation and/or other materials provided with the distribution.
 *   3. The name of the ATLAS group or the names of its contributers may
 *      not be used to endorse or promote products derived from this
 *      software without specific written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
 * TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 * PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE ATLAS GROUP OR ITS CONTRIBUTORS
 * BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 *
 */

#include "atlas_misc.h"
#include "atlas_level1.h"
#include <math.h>

void Mjoin(PATL,ssq)(const int N, const TYPE *X, const int incX,
                     TYPE *scal0, TYPE *ssq0)
{
#ifdef TREAL
   int i;
   TYPE t0, ssq=(*ssq0), ax, scal=(*scal0);

   if (N < 1 || incX < 1) return;
   for (i=N; i; i--, X += incX)
   {
      ax = *X;
      if (ax != ATL_rzero)
      {
         ax = Mabs(ax);
         if (scal < ax)
         {
            t0 = scal / ax;
            t0 *= t0;
            ssq = ATL_rone + ssq * t0;
            scal = ax;
         }
         else
         {
            t0 = ax / scal;
            ssq += t0*t0;
         }
      }
   }
   *scal0 = scal;
   *ssq0 = ssq;
#else
   int i;
   const int incx = incX<<1;
   TYPE t0, ax, ssq=(*ssq0), scal=(*scal0);

   if (N < 1 || incX < 1) return;
   else if (incX == 1) Mjoin(PATLU,ssq)(N<<1, X, 1, scal0, ssq0);
   else
   {
      for (i=N; i; i--, X += incx)
      {
         ax = *X;
         if (ax != ATL_rzero)
         {
            ax = Mabs(ax);
            if (scal < ax)
            {
               t0 = scal / ax;
               t0 *= t0;
               scal = ax;
               ssq = ATL_rone + ssq * t0;
            }
            else
            {
               t0 = ax / scal;
               ssq += t0*t0;
            }
         }
         ax = X[1];
         if (ax != ATL_rzero)
         {
            ax = Mabs(ax);
            if (scal < ax)
            {
               t0 = scal / ax;
               t0 *= t0;
               scal = ax;
               ssq = ATL_rone + ssq * t0;
            }
            else
            {
               t0 = ax / scal;
               ssq += t0*t0;
            }
         }
      }
      *scal0 = scal;
      *ssq0 = ssq;
   }
#endif
}

#ifdef TREAL
TYPE Mjoin(PATL,nrm2)
#else
TYPE Mjoin(Mjoin(Mjoin(ATL_,UPR),PRE),nrm2)
#endif
   (const int N, const TYPE *X, const int incX)
{
   TYPE ssq=ATL_rone, scal=ATL_rzero;
   if (N < 1 || incX < 1) return(ATL_rzero);
   #ifdef TREAL
      else if (N == 1) return(Mabs(*X));
   #endif
   Mjoin(PATL,ssq)(N, X, incX, &scal, &ssq);
   return(scal * sqrt(ssq));
}
