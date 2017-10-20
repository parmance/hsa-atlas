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

#ifdef TREAL
static int iamax_1(const int N, const TYPE *X, const int incX)
{
   const int n = (N >> 2)<<2;
   const int nr = N - n;
   const TYPE *x=X, *stX = X + n, *xmax;
   register TYPE pmax=0.0, nmax=0.0, x0, x1, x2, x3;

   if (N < 2) return(0);
   if (n)
   {
      x0 = *x;
      x1 = x[1];
      x2 = x[2];
      x3 = x[3];
      x += 4;
      xmax = x;
      if (n != 4)
      {
         do
         {
            if (x0 <= pmax && x0 >= nmax) goto L10;
            if (x0 > pmax)
            {
               nmax = -x0;
               pmax =  x0;
               xmax = x;
            }
            else /* if (x0 < nmax) */
            {
               pmax = -x0;
               nmax =  x0;
               xmax = x;
            }
L10 :        x0 = *x;
            if (x1 <= pmax && x1 >= nmax) goto L20;
            if (x1 > pmax)
            {
               nmax = -x1;
               pmax =  x1;
               xmax = x+1;
            }
            else /* if (x1 < nmax) */
            {
               pmax = -x1;
               nmax =  x1;
               xmax = x+1;
            }
L20 :        x1 = x[1];
            if (x2 <= pmax && x2 >= nmax) goto L30;
            if (x2 > pmax)
            {
               nmax = -x2;
               pmax =  x2;
               xmax = x+2;
            }
            else /* if (x2 < nmax) */
            {
               pmax = -x2;
               nmax =  x2;
               xmax = x+2;
            }
L30 :        x2 = x[2];
            if (x3 <= pmax && x3 >= nmax) goto L40;
            if (x3 > pmax)
            {
               nmax = -x3;
               pmax =  x3;
               xmax = x+3;
            }
            else /* if (x3 < nmax) */
            {
               pmax = -x3;
               nmax =  x3;
               xmax = x+3;
            }
L40 :        x3 = x[3];
            x += 4;
         }
         while (x != stX);
      }
      if (x0 <= pmax && x0 >= nmax) goto L15;
      if (x0 > pmax)
      {
         nmax = -x0;
         pmax =  x0;
         xmax = x;
      }
      else /* if (x0 < nmax) */
      {
         pmax = -x0;
         nmax =  x0;
         xmax = x;
      }
L15 :
      if (x1 <= pmax && x1 >= nmax) goto L25;
      if (x1 > pmax)
      {
         nmax = -x1;
         pmax =  x1;
         xmax = x+1;
      }
      else /* if (x1 < nmax) */
      {
         pmax = -x1;
         nmax =  x1;
         xmax = x+1;
      }
L25 :
      if (x2 <= pmax && x2 >= nmax) goto L35;
      if (x2 > pmax)
      {
         nmax = -x2;
         pmax =  x2;
         xmax = x+2;
      }
      else /* if (x2 < nmax) */
      {
         pmax = -x2;
         nmax =  x2;
         xmax = x+2;
      }
L35 :
      if (x3 <= pmax && x3 >= nmax) goto L45;
      if (x3 > pmax)
      {
         nmax = -x3;
         pmax =  x3;
         xmax = x+3;
      }
      else /* if (x3 < nmax) */
      {
         pmax = -x3;
         nmax =  x3;
         xmax = x+3;
      }
L45 :
      xmax -= 4;
   }
   else xmax = X+1;
   if (nr)
   {
      stX = x + nr;
      do
      {
         x0 = *x++;
         if (x0 <= pmax && x0 >= nmax) continue;
         if (x0 > pmax)
         {
            nmax = -x0;
            pmax =  x0;
            xmax = x-1;
         }
         else /* if (x0 < nmax) */
         {
            pmax = -x0;
            nmax =  x0;
            xmax = x-1;
         }
      }
      while (x != stX);
   }
   if (pmax == 0.0 || nmax == 0.0) return(0);
   return(xmax-X);
}
#endif

int Mjoin(Mjoin(ATL_i,PRE),amax)(const int N, const TYPE *X, const int incX)
{
#ifdef TREAL
   int i, imax=N;
   register TYPE pmax=0.0, nmax=0.0, x0;

   if (incX == 1) return(iamax_1(N, X, 1));
   if (N < 2) return(0);
   for(i=N; i; i--, X += incX)
   {
      x0 = *X;
      if (x0 <= pmax && x0 >= nmax) continue;
      if (x0 > pmax)
      {
         nmax = -x0;
         pmax =  x0;
         imax = i;
      }
      else   /* if (x0 < nmax) */
      {
         nmax =  x0;
         pmax = -x0;
         imax = i;
      }
   }
#else
   int i, imax=N;
   const int incx = incX<<1;
   register TYPE pmax=0, nmax=0, x0, x1, tmp;
   if (N > 1)
   {
      if (incX == 1)
      {
         for(i=N; i; i--, X += 2)
         {
            x0 = *X;
            tmp = X[1];
            x1 = x0 - tmp;
            x0 += tmp;

            if (x0 >= x1)
            {
               if (x0 <= pmax && x1 >= nmax) continue;
               if (x0 > pmax) { pmax = x0; nmax = -x0; }
               else { pmax = -x1; nmax = x1; }
            }
            else
            {
               if (x1 <= pmax && x0 >= nmax) continue;
               if (x1 > pmax) { pmax = x1; nmax = -x1; }
               else { pmax = -x0; nmax = x0; }
            }
            imax = i;
         }
      }
      else
      {
         for(i=N; i; i--, X += incx)
         {
            x0 = *X;
            tmp = X[1];
            x1 = x0 - tmp;
            x0 += tmp;

            if (x0 >= x1)
            {
               if (x0 <= pmax && x1 >= nmax) continue;
               if (x0 > pmax) { pmax = x0; nmax = -x0; }
               else { pmax = -x1; nmax = x1; }
            }
            else
            {
               if (x1 <= pmax && x0 >= nmax) continue;
               if (x1 > pmax) { pmax = x1; nmax = -x1; }
               else { pmax = -x0; nmax = x0; }
            }
            imax = i;
         }
      }
   }

#endif
   return(N-imax);
}
