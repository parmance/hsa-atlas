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
#ifdef DIRECTHSA /* DEVTEMP */
#define HSADECLS
#endif
#include "atlas_misc.h"

HSA_FUNCTION
void Mjoin3(PATL,gescal,PHSA_FN)
(const int M, const int N, const SCALAR beta, TYPE *C, const int ldc)
/*
 * C <- beta*C
 */
{
#ifdef TREAL
   if (beta == ATL_rzero) Mjoin3(PATL,gezero,PHSA_FN)(M, N, C, ldc);
   else if (beta == ATL_rone) return;
   else Mjoin3(PATL,gescal_bX,PHSA_FN)(M, N, beta, C, ldc);
#else
   TYPE rbeta = *beta;
   if (beta[1] == ATL_rzero)
   {
      if (rbeta == ATL_rzero) Mjoin(PATL,gezero)(M, N, C, ldc);
      else if (rbeta == ATL_rone) return;
      else Mjoin(PATL,gescal_bXi0)(M, N, beta, C, ldc);
   }
   else Mjoin(PATL,gescal_bX)(M, N, beta, C, ldc);
#endif
}

#ifdef DIRECTHSA

typedef struct Mjoin(PATL,gescal_args_s) {
   const int M;
   const int N;
   const SCALAR beta;
   TYPE *C;
   const int ldc;
} Mjoin(PATL,gescal_args_t);

HSA_KERNEL
void Mjoin(PATL,gescal_kernel)(Mjoin(PATL,gescal_args_t)* args)
{
   const int M = args->M;
   const int N = args->N;
   const SCALAR beta = args->beta;
   TYPE *C = args->C;
   const int ldc = args->ldc;

   Mjoin3(PATL,gescal,PHSA_FN)(M, N, beta, C, ldc);
}

void Mjoin3(PATL,gescal,PHSA)
   (const int M, const int N, const SCALAR beta, TYPE *C, const int ldc)
{
   Mjoin(PATL,gescal_args_t) args = { M, N, beta, C, ldc };
   HSA_LAUNCH(Mjoin(PATL,gescal_kernel), &args);
}

#endif /* DIRECTHSA */
