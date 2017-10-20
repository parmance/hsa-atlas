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
#  define HSADECLS
#endif
#include "atlas_misc.h"

HSA_FUNCTION
void Mjoin3(PATL,gemove,PHSA_FN)
   (const int M, const int N, const SCALAR alpha, const TYPE *A, const int lda,
    TYPE *C, const int ldc)
/*
 * C <- alpha * A
 */
{
#ifdef TREAL
   if (alpha == ATL_rone) Mjoin3(PATL,gecopy,PHSA_FN)(M, N, A, lda, C, ldc);
   else if (alpha == ATL_rzero) Mjoin3(PATL,gezero,PHSA_FN)(M, N, C, ldc);
   else Mjoin3(PATL,gemove_aX,PHSA_FN)(M, N, alpha, A, lda, C, ldc);
#else
   TYPE ralpha = *alpha;

   if (alpha[1] == ATL_rzero)
   {
      if (ralpha == ATL_rone) Mjoin(PATL,gecopy)(M, N, A, lda, C, ldc);
      else if (ralpha == ATL_rzero) Mjoin(PATL,gezero)(M, N, C, ldc);
      else Mjoin(PATL,gemove_aXi0)(M, N, alpha, A, lda, C, ldc);
   }
   else Mjoin(PATL,gemove_aX)(M, N, alpha, A, lda, C, ldc);
#endif
}

#ifdef DIRECTHSA

typedef struct Mjoin(PATL,gemove_args_s) {
   const int M;
   const int N;
   const SCALAR alpha;
   const TYPE *A;
   const int lda;
   TYPE *C;
   const int ldc;
} Mjoin(PATL,gemove_args_t);

HSA_KERNEL
void Mjoin3(PATL,gemove,_kernel)(Mjoin(PATL,gemove_args_t)* args)
{
   const int M = args->M;
   const int N = args->N;
   const SCALAR alpha = args->alpha;
   const TYPE *A = args->A;
   const int lda = args->lda;
   TYPE *C = args->C;
   const int ldc = args->ldc;

   Mjoin3(PATL,gemove,PHSA_FN)(M, N, alpha, A, lda, C, ldc);
}

void Mjoin3(PATL,gemove,PHSA)
(const int M, const int N, const SCALAR alpha, const TYPE *A, const int lda,
 TYPE *C, const int ldc)
{
   Mjoin(PATL,gemove_args_t) args = { M, N, alpha, A, lda, C, ldc };
   HSA_LAUNCH(Mjoin3(PATL,gemove,_kernel), &args);
}

#endif /* DIRECTHSA */
