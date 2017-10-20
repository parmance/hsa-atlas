/*
 *             Automatically Tuned Linear Algebra Software v3.10.3
 * Copyright (C) 2012 Anthony M. Castaldo
 *
 * Code contributers : Anthony M. Castaldo, Siju Samuel, R. Clint Whaley
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
/*----------------------------------------------------------------------------*/
/* *** WARNING *** WARNING *** WARNING *** WARNING *** WARNING *** WARNING ****/
/* WARNING: THIS CODE IS ONLY LIGHTLY TESTED, AT BEST. I believe I made a     */
/* faithful translation of the LAPACK Fortran code. Tony Castaldo Jan 2012.   */
/*----------------------------------------------------------------------------*/
/*
 * This is the C translation of the standard LAPACK Fortran routine:
 *     SUBROUTINE DLARFGP( N, ALPHA, X, INCX, TAU )
 *
 * ATL_larfgp.c :
 * void ATL_larfgp( const int N, TYPE *ALPHA, TYPE *X, int INCX, TYPE *TAU)
 *
 *  Purpose
 *  =======
 *
 *  Real Precision
 *  --------------
 *
 *  ATL_larfgp generates a real/complex elementary reflector H of order n,
 *  such that:
 *
 *        H * ( alpha ) = ( beta ),   H' * H = I.
 *            (   x   )   (   0  )
 *
 *  where alpha and beta are scalars, beta is POSITIVE,
 *  and x is an (n-1)-element real vector. H is represented in the form
 *
 *        H = I - tau * ( 1 ) * ( 1 v' ) ,                     (Real Precisions)
 *                      ( v )
 *
 *        H = I - tau * ( 1 ) * ( 1 conjugate(v)' ) ,       (Complex Precisions)
 *                      ( v )
 *
 *  where tau is a real/complex scalar and v is a real/complex (n-1)-element
 *  vector.
 *
 *  If the elements of x are all zero, then tau = 0 and H is taken to be
 *  the unit matrix.
 *
 *  Otherwise  1 <= tau <= 2.
 *
 *  The difference between this routine and larfg is that beta is guaranteed to
 *  be non-negative; this also ensures a unique H, and for a QR factorization,
 *  a unique Q for any given matrix A, along with additional mathematical
 *  properties for R that do not necessarily exist with negative numbers on the
 *  diagonal. However, empirical analysis indicates these features come at the
 *  price of some minor losses in precision.
 *
 *
 *  Arguments
 *  =========
 *
 *  N       (input) INTEGER
 *          The order of the elementary reflector.
 *
 *  ALPHA   (input/output)
 *          On entry, the value alpha.
 *          On exit, it is overwritten with the value beta.
 *
 *  X       (input/output)   array pointer, dimension
 *                         (1+(N-2)*abs(INCX))
 *          On entry, the vector x.
 *          On exit, it is overwritten with the vector v.
 *
 *  INCX    (input) INTEGER
 *          The increment between elements of X. INCX > 0.
 *
 *  TAU     (output)
 *          value of tau.
 */
#include "atlas_misc.h"
#include <math.h>
#include "cblas.h"
#include "atlas_lapack.h"
#include "atlas_lamch.h"

/*----------------------------------------------------------------------------*/
/* This is NOT a textbook Householder reflection, because [1 v]^T does not    */
/* have a norm of 1. The purpose seems to be to save the first storage        */
/* location; the '1' is not stored anywhere.                                  */
/*                                                                            */
/* Real and Complex do share some structure, but they are complicated enough  */
/* as it is, so for clarity we just have two codes separated by an IFDEF.     */
/*----------------------------------------------------------------------------*/

void ATL_larfgp(ATL_CINT N, TYPE *ALPHA, TYPE *X, ATL_CINT INCX, TYPE *TAU)
{
   #ifdef TREAL

   TYPE TWO=2.0, ONE=1.0, ZERO=0.0, ALPHR, BETA, BETAp, BIGNUM, SMLNUM, XNORM;
   TYPE TEMP, SAVEALPHA;
   int    j, KNT;

   if (N < 0)
   {
      *TAU = ZERO;
      return;
   }

   XNORM = cblas_nrm2(N-1, X, INCX);      /* Norm2 of X. */
   ALPHR = *ALPHA;                        /* Use ALPHR throughout... */

   if (XNORM == ZERO)
   {
      /* Set H = [+/-1, 0; I], sign is chosen so ALPHA >= 0. */
      if (ALPHR >= ZERO )
      {
         /* When TAU.eq.ZERO, the vector is special-cased to be     */
         /* all zeros in the application routines.  We do not need  */
         /* to clear it. Also, *ALPHA remains as is.                */
         *TAU = ZERO;
      }
      else
      {
         /* However, the application routines rely on explicit zero */
         /* checks when TAU.ne.ZERO, so we must clear X.            */
         *TAU = TWO;
         for (j=0; j<(N-1); j++)
         {
            X[j*INCX] = ZERO;
         }
         *ALPHA = -ALPHR;                 /* Must negate alpha in data. */
      }
   }
   else                                   /* XNORM NE ZERO. */
   {
      BETAp = ATL_lapy2(ALPHR, XNORM);    /* Full norm2=sqrt(a^2+b^2) */
      BETA = BETAp;                       /* Assume ALPHA < 0 */
      if (ALPHR < 0) BETA= -BETAp;        /* Change sign if alpha negative.*/
      SMLNUM = ATL_laSAFMIN;              /* Get smallest safe invertible. */
      KNT = 0;
      if (BETAp < SMLNUM)
      {
         /* XNORM, BETA may be inaccurate; scale X and recompute them. */
         BIGNUM = ONE / SMLNUM;           /* Get large number. */

         while (BETAp < SMLNUM)
         {
            KNT++;
            cblas_scal(N-1, BIGNUM, X, INCX);
            BETA *= BIGNUM;
            BETAp *= BIGNUM;
            ALPHR *= BIGNUM;
         }
         /* New BETA is at most 1, at least SMLNUM. */

         XNORM = cblas_nrm2(N-1, X, INCX);
         BETAp = ATL_lapy2(ALPHR, XNORM);   /* Will always be positive */
         BETA = BETAp;                      /* Assume ALPHA < 0        */
         if (ALPHR > 0) BETA=-BETAp;        /* -SIGN(BETA, ALPHA)      */
      } /* BETAp < SAFMIN */

      SAVEALPHA = ALPHR;
      ALPHR = ALPHR + BETA;
      if (BETA < ZERO)
      {
         BETA = -BETA;
         *TAU = -(ALPHR/BETA);
      }
      else
      {
         ALPHR = XNORM * (XNORM/ALPHR);
         *TAU   = ALPHR/BETA;
         ALPHR = -ALPHR;
      }

      #if defined(SREAL)
      TEMP = fabsf(*TAU);
      #else
      TEMP = fabs(*TAU);
      #endif

      if (TEMP <= SMLNUM)
      {
         /* In the case where the computed TAU ends up being a denormalized  */
         /* number, it loses relative accuracy. This is a BIG problem.       */
         /* Solution: flush TAU to ZERO. This explains the next IF statement.*/
         /*                                                                  */
         /*(Bug report provided by Pat Quillen from MathWorks Jul 29, 2009.) */
         /*(Thanks Pat. Thanks MathWorks.)                                   */
         if (SAVEALPHA >= ZERO)
         {
            *TAU = ZERO;
         } else {
            *TAU = TWO;
            for (j=0; j<(N-1); j++)
            {
               X[j*INCX] = ZERO;
            }

            BETA = -SAVEALPHA;
         }
      } else {
         /* This is the general case. */
         cblas_scal(N-1, ONE / ALPHR, X, INCX);
      }

      for (j=0; j<KNT; j++)
      {
         BETA *= SMLNUM;
      }

      *ALPHA = BETA;

   } /* end XNORM !=0 */
   return; /* END of REAL Precision ATL_larfgp. */

#else /* COMPLEX CASE  of ATL_larfgp. */
   TYPE TWO=2.0, ONE=1.0, ZERO=0.0, BETA, BETAp, BIGNUM, SMLNUM, XNORM,
        ALPHI, ALPHR, TEMP;
   TYPE ONEVAL[2] =  {ATL_rone, ATL_rzero};
   TYPE SAVEALPHA[2] =  {ATL_rzero, ATL_rzero};
   int j, KNT;

   if ( N < 0)
   {
/*
 *    H  =  I
 */
      *(TAU)  = 0.0;
      *(TAU + 1) = 0.0;
      return;
   }

   XNORM = cblas_nrm2(N-1, X, INCX);
   ALPHR = *(ALPHA);          /* Real component. */
   ALPHI = *(ALPHA+1);        /* Imag component. */

   if (XNORM == ZERO )
   {
      /* Set H=[1-alpha/abs(alpha) 0; 0 I], sign chosen so ALPHA >= 0. */
      if ( ALPHI == ZERO )
      {
         if(ALPHR > ZERO)
         {
/*          When TAU.eq.ZERO, the vector is special-cased to be
 *          all zeros in the application routines.  We do not need
 *          to clear it.
 */
            *(TAU)  = 0.0;
            *(TAU + 1) = 0.0;
         }
         else
         {
/*          However, the application routines rely on explicit
 *          zero checks when TAU.ne.ZERO, and we must clear X.
 */
            *(TAU)  = 2.0;
            *(TAU+1)= 0.0;

            for (j=0; j<(N-1); j++)                /* For all entries,        */
            {
               X[((j*INCX) SHIFT)]   = ZERO;       /* .. Clear real part.     */
               X[((j*INCX) SHIFT)+1] = ZERO;       /* .. Clear imaginary part.*/
            }

            *(ALPHA)  = -ALPHR;              /* Negate Alpha. */
            *(ALPHA+1)= -ALPHI;
         }
      }                                     /* ALPHI != ZERO                  */
      else
      {
/*       Only "reflecting" the diagonal entry to be real and non-negative.    */

         XNORM = ATL_lapy2(ALPHR,ALPHI);     /* Get abs value of Alpha. */
         *(TAU)   = ONE-(ALPHR/XNORM);
         *(TAU+1) = -(ALPHI/XNORM);

         for (j=0; j<(N-1); j++)                   /* For all entries,        */
         {
            X[((j*INCX) SHIFT)]   = ZERO;          /* .. Clear real part.     */
            X[((j*INCX) SHIFT)+1] = ZERO;          /* .. Clear imaginary part.*/
         }

         *(ALPHA)   = XNORM;              /* Real Part of alpha */
         *(ALPHA+1) = 0.0;                /* Set Imaginary part to Zero */
      }                                   /* ALPHI = 0 ends */
   }
   else  /* general case, XNORM != zero. */
   {
      BETAp = ATL_lapy3(ALPHR, ALPHI, XNORM);  /* Full norm 2. */
      BETA = BETAp;                 /* Assume ALPHR < 0 */
      if (ALPHR < 0)                /* If AlphR negative, */
         BETA = -BETAp;             /* Make beta match sign of ALPHR  */

      SMLNUM = ATL_laSAFMIN;        /* Get safe inversion number. */
      BIGNUM = ONE / SMLNUM;
      KNT = 0;

      if (BETAp < SMLNUM)           /* If norm2 is too small, */
      {
         /* XNORM, BETA may be inaccurate; scale X and recompute them. */
         while (BETAp < SMLNUM)
         {
            KNT++;
            #ifdef DCPLX
                cblas_zdscal(N-1, BIGNUM, X, INCX);
            #else
                cblas_csscal(N-1, BIGNUM, X, INCX);
            #endif
            BETA *= BIGNUM;
            BETAp *= BIGNUM;
            ALPHI = ALPHI*BIGNUM;
            ALPHR = ALPHR*BIGNUM;
         }

         /* New BETA is at most 1, at least SMLNUM */
         XNORM = cblas_nrm2(N-1, X, INCX);
         *(ALPHA)  = ALPHR;
         *(ALPHA+1)= ALPHI;

         BETAp = ATL_lapy3(ALPHR, ALPHI, XNORM); /* always positive. */
         BETA = BETAp;
         if (ALPHR < 0) BETA= -BETAp;  /* SIGN(BETA, ALPHR) */
      } /* END if BETAp < SMLNUM */

      SAVEALPHA[0] = *(ALPHA);
      SAVEALPHA[1] = *(ALPHA+1);
      *ALPHA = *ALPHA + BETA;

      if (BETA < ZERO)
      {
         BETA = -BETA;

         /* TAU = -ALPHA/BETA (Note TAU is complex but Beta is REAL) */
         *TAU     = -(*ALPHA)/BETA;
         *(TAU+1) = -(*(ALPHA+1))/BETA;
      }
      else
      {
         ALPHR = ALPHI * (ALPHI/ (*ALPHA));
         ALPHR = ALPHR + XNORM * (XNORM/(*ALPHA));
         *TAU     = ALPHR/BETA;
         *(TAU+1) = - (ALPHI/BETA);
         *ALPHA     = -ALPHR;
         *(ALPHA+1) = ALPHI;
      }

      /* Perform complex division before scaling the X vector */
      ATL_ladiv(ONEVAL, ALPHA, ALPHA);  /* ALPHA will have the result */

      TEMP = ATL_lapy2( TAU[0], TAU[1]); /* Get abs(TAU). */
      if (TEMP <= SMLNUM)
      {
         /* In the case where the computed TAU ends up being a             */
         /* denormalized number, it loses relative accuracy. This is a     */
         /* BIG problem. Solution: flush TAU to ZERO (or TWO or whatever   */
         /* makes a nonnegative real number for BETA).                     */

         /* (Bug report provided by Pat Quillen from MathWorks on Jul      */
         /* 29, 2009.) (Thanks Pat. Thanks MathWorks.)                     */

         ALPHR = SAVEALPHA[0];
         ALPHI = SAVEALPHA[1];
         if (ALPHI == ZERO)
         {
            if (ALPHR >= ZERO)
            {
               *(TAU) = ZERO;
               *(TAU+1) = ZERO;
            } else {
               *(TAU) = TWO;
               *(TAU+1) = ZERO;
               for (j=0; j<(N-1); j++)
               {
                  X[((j*INCX) SHIFT)] = ZERO;
                  X[((j*INCX) SHIFT)+1] = ZERO;
               }

               BETA = -ALPHR;             /* Note ALPHI == zero, here. */
            }
         } else {
            XNORM = ATL_lapy2(ALPHR, ALPHI);
            *(TAU) = ONE - (ALPHR/XNORM);
            *(TAU+1) = -(ALPHI/XNORM);
            for (j=0; j<(N-1); j++)
            {
               X[((j*INCX) SHIFT)] = ZERO;
               X[((j*INCX) SHIFT)+1] = ZERO;
            }

            BETA = XNORM;
         }
      } else /* TAU large enough, handle general case */
      {
         cblas_scal(N-1, ALPHA, X, INCX);
      }

      /* If BETA is subnormal, it may lose relative accuracy. */
      for (j=0; j<KNT; j++)
      {
         BETA *= SMLNUM;
      }

     *(ALPHA)   = BETA;          /* Note BETA is just REAL, not complex. */
     *(ALPHA+1) = ZERO;          /* ...                                  */
   }  /* XNORM != 0 ends  */
   return;
#endif /* COMPLEX CASE */
} /* END ATL_larfgp */
