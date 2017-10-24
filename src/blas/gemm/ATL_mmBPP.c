/*
 *             Automatically Tuned Linear Algebra Software v3.10.3
 *                    (C) Copyright 2007 R. Clint Whaley
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

#ifdef DIRECTHSA
#define ATL_NoICalls
#endif

#include "atlas_misc.h"
#include "atlas_lvl3.h"
#include "atlas_malloc.h"

HSA_FUNCTION
int Mjoin3(PATL,mmBPP,PHSA)
   (MemBlob* memBlob,
    const enum ATLAS_TRANS TA, const enum ATLAS_TRANS TB,
    const int M, const int N, const int K,
    const SCALAR alpha, const TYPE *A, const int lda,
    const TYPE *B, const int ldb, const SCALAR beta,
    TYPE *C, const int ldc0)
/*
 * Copy algorithm, assuming M <= MB && N <= NB, K large (shape: block, panel,
 * panel); copies A and B on-the-fly
 */
{
   void *vC;
   TYPE *pA, *pB, *pC;
   MAT2BLK2 A2blk = ATL_NullFn;
   MAT2BLK2 B2blk = ATL_NullFn;
   NBMM0 NBmm0 = ATL_NullFn;
   NBMM0 NBmm1 = ATL_NullFn;
   size_t ldc, incA, incB;
   int m, n, nblk, k, kr;

   if (M > MB || N > NB)  /* don't handle multiple M/N blocks */
      return 1;
   if (M < MB && M+ATL_mmMU >= MB)
      m = MB;
   else
      m = M;
   if (N < NB && N+ATL_mmNU >= NB)
      n = NB;
   else
      n = N;
   ldc = (((m*sizeof(TYPE)+ATL_Cachelen-1)/ATL_Cachelen)*ATL_Cachelen)
          / sizeof(TYPE);
   vC = Mjoin(ATL_Malloc,PHSA)(memBlob,
                               ATL_Cachelen+ATL_MulBySize(ldc*n+KB*(m+n)));
   if (!vC)
      return -1;
   pC = ATL_AlignPtr(vC);
   pA = pC + ldc*n;
   pB = pA + KB*m;
   if (TA == AtlasNoTrans)
   {
      A2blk = ATL_TargetFn(Mjoin3(PATL,gemoveT,PHSA));
      incA = lda*KB;
   }
   else
   {
      A2blk = ATL_TargetFn(Mjoin3(PATL,gemove,PHSA));
      incA = KB;
   }
   if (TB == AtlasNoTrans)
   {
      B2blk = ATL_TargetFn(Mjoin3(PATL,gemove,PHSA));
      incB = KB;
   }
   else
   {
      B2blk = ATL_TargetFn(Mjoin3(PATL,gemoveT,PHSA));
      incB = ldb*KB;
   }
/*
 * If we are going to multiply zeros to avoid cleanup, zero workspace
 */
   if (m != M || n != N)
      Mjoin3(PATL,zero,PHSA)(ldc*n+KB*(m+n), pC, 1);
/*
 * See what kernel we're calling
 */
   if (m == MB)
   {
      if (n == NB)  /* no cleanup */
      {
         NBmm0 = ATL_TargetFn(NBmm_b0);
         NBmm1 = ATL_TargetFn(NBmm_b1);
      }
      else         /* need to call N-cleanup kernel */
      {
         NBmm0 = ATL_TargetFn(Mjoin3(PATL,pNBmm_b0,PHSA));
         NBmm1 = ATL_TargetFn(Mjoin3(PATL,pNBmm_b1,PHSA));
      }
   }
   else if (n == NB) /* call M-cleanup kernel */
   {
      NBmm0 = ATL_TargetFn(Mjoin3(PATL,pMBmm_b0,PHSA));
      NBmm1 = ATL_TargetFn(Mjoin3(PATL,pMBmm_b1,PHSA));
   }
   else              /* both N & M are cleanup, call general K clean */
   {
      NBmm0 = ATL_TargetFn(Mjoin3(PATL,pKBmm,PHSA));
      NBmm1 = ATL_TargetFn(Mjoin3(PATL,pKBmm,PHSA));
      if (m == M && n == N)  /* must zero pC if not done above */
         Mjoin3(PATL,zero,PHSA)(ldc*n, pC, 1);

   }
   nblk = K / KB;
   kr = K - nblk*KB;
   if (!nblk && kr)
      Mjoin3(PATL,zero,PHSA)(ldc*n, pC, 1);
   Mjoin3(PATL,mmK,PHSA)(M, m, N, n, nblk, kr, (kr && kr+4 >= KB) ? KB : 0,
                         ATL_rone, ATL_rone, ATL_rzero, A, lda, incA, pA, 0,
                         B, ldb, incB, pB, 0, pC, ldc,
                         A2blk, B2blk, NBmm0, NBmm1);
   Mjoin3(PATL,geadd,PHSA)(M, N, alpha, pC, ldc, beta, C, ldc0);
   Mjoin(ATL_Free,PHSA)(memBlob, vC);
   return 0;
}

