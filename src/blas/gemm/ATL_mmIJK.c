/*
 *             Automatically Tuned Linear Algebra Software v3.10.3
 *                    (C) Copyright 1997 R. Clint Whaley
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
#  define ATL_no_icalls
#  define HSADECLS
#endif

#include "atlas_misc.h"
#include "atlas_lvl3.h"
#include "atlas_malloc.h"

#define KBmm Mjoin3(PATL,pKBmm,PHSA_FN)
#define IBNBmm Mjoin3(PATL,IBNBmm,PHSA_FN)
#define NBJBmm Mjoin3(PATL,MBJBmm,PHSA_FN)
#define IBJBmm Mjoin3(PATL,IBJBmm,PHSA_FN)
#define col2blk Mjoin3(PATL,col2blk_a1,PHSA_FN)

#define NBMM_ICALL Mjoin3(PATL,icall_site_NBMM0,HSADECL)
#define PUTBLK_ICALL Mjoin3(PATL,icall_site_PUTBLK,HSADECL)
#define MAT2BLK_ICALL Mjoin3(PATL,icall_site_MAT2BLK,HSADECL)
#include "ATL_indir_call.c"

HSA_FUNCTION
void Mjoin3(PATL,mmIJK2,PHSA_FN)(
   int K, int nMb, int nNb, int nKb, int ib, int jb,
   int kb, const SCALAR alpha, const TYPE *A, int lda,
   TYPE *pA0, int incA, MAT2BLK A2blk, const TYPE *pB0,
   const SCALAR beta, TYPE *C, int ldc, TYPE *pC,
   PUTBLK putblk, NBMM0 NBmm0)
/*
 * Outer three loops for matmul with outer loop over rows of A
 */
{
   int i, j, ldpc;
   const int ZEROC = ((putblk == ATL_NullFn) && SCALAR_IS_ZERO(beta));
   const int incK = ATL_MulByNB(K), incC = ATL_MulByNB(ldc);
   TYPE *pA=pA0, *stA=pA0+ATL_MulByNBNB(nKb);
   const TYPE *pB=pB0;
   const TYPE cubeta = ( (putblk) ? ATL_rzero : beta );
   TYPE *c;

   if (putblk)
   {
      ldpc = NB;
      if (!nKb && kb) Mjoin3(PATL,gezero,PHSA_FN)(MB, NB, pC, MB);
   }
   else ldpc = ldc;
   for (i=nMb; i; i--)    /* loop over full row panels of A */
   {
      if (A)
      {
         MAT2BLK_ICALL(A2blk, K, NB, A, lda, pA, alpha);  /* get 1 row panel of A */
         A += incA;
      }
      if (!putblk) pC = C;
      c = C;
      C += NB;
      for (j=nNb; j; j--)  /* full column panels of B */
      {
         if (nKb)
         {
            NBMM_ICALL(NBmm0, MB, NB, KB, ATL_rone, pA, KB, pB, KB, beta,
                       pC, ldpc);
            pA += NBNB;
            pB += NBNB;
            if (nKb != 1)
            {
               do
               {
                  NBmm(MB, NB, KB, ATL_rone, pA, KB, pB, KB, ATL_rone,
                       pC, ldpc);
                  pA += NBNB;
                  pB += NBNB;
               }
               while (pA != stA);
            }
            if (kb)
            {
               KBmm(MB, NB, kb, ATL_rone, pA, kb, pB, kb, ATL_rone, pC, ldpc);
               pB += kb*NB;
            }
         }
         else
         {
            if (ZEROC) Mjoin3(PATL,gezero,PHSA_FN)(MB, NB, pC, ldpc);
            if (kb)
            {
               KBmm(MB, NB, kb, ATL_rone, pA, kb, pB, kb, cubeta, pC, ldpc);
               pB += kb*NB;
            }
         }
         pA = pA0;
         if (putblk) PUTBLK_ICALL(putblk, NB, NB, pC, c, ldc, beta);
         else pC += incC;
         c += incC;
      }
      if (jb)
      {
         NBJBmm(jb, K, pA, pB, cubeta, pC, ldpc);
         if (putblk) PUTBLK_ICALL(putblk, NB, jb, pC, c, ldc, beta);
      }
      pB = pB0;
      if (!A)
      {
         pA0 += incK;
         pA = pA0;
         stA += incK;
      }
   }
   if (ib)
   {
      c = C;
      /* get last row panel of A */
      if (A) MAT2BLK_ICALL(A2blk, K, ib, A, lda, pA, alpha);
      for (j=nNb; j; j--)  /* full column panels of B */
      {
         if (putblk)
         {
            IBNBmm(ib, K, pA, pB, ATL_rzero, pC, ib);
            PUTBLK_ICALL(putblk, ib, NB, pC, c, ldc, beta);
         }
         else IBNBmm(ib, K, pA, pB, beta, c, ldc);
         pB += incK;
         c += incC;
      }
      if (jb)
      {
         if (putblk)
         {
            IBJBmm(ib, jb, K, pA, pB, ATL_rzero, pC, ib);
            PUTBLK_ICALL(putblk, ib, jb, pC, c, ldc, beta);
         }
         else IBJBmm(ib, jb, K, pA, pB, beta, c, ldc);
      }
   }
}

HSA_FUNCTION
int Mjoin3(PATL,mmIJK,PHSA_FN)(
   MemBlob* memBlob,
   const enum ATLAS_TRANS TA, const enum ATLAS_TRANS TB,
   const int M, const int N0, const int K,
   const SCALAR alpha, const TYPE *A, const int lda0,
   const TYPE *B, const int ldb0, const SCALAR beta,
   TYPE *C, const int ldc0)
{
   size_t incA, incB, incC;
   const size_t lda=lda0, ldb=ldb0, ldc=ldc0;
   const size_t incK = ATL_MulByNB((size_t)K);
   int N = N0;
   int nMb, nNb, nKb, ib, jb, kb, jb2, h, i, j, k, n;
   void *vA=NULL, *vC=NULL;
   TYPE *pA, *pB, *pC;
   MAT2BLK A2blk, B2blk;
   A2blk = B2blk = ATL_NullFn;
   PUTBLK putblk = ATL_NullFn;
   NBMM0 NBmm0 = ATL_NullFn;

   nMb = ATL_DivByNB(M);
   nNb = ATL_DivByNB(N);
   nKb = ATL_DivByNB(K);
   ib = M - ATL_MulByNB(nMb);
   jb = N - ATL_MulByNB(nNb);
   kb = K - ATL_MulByNB(nKb);

/*
 * If K sufficiently large, write to temporary C as safety measure;  otherwise
 * write directly to C
 */
   if (nKb < 12)
   {
      putblk = ATL_NullFn;
      pC = C;
      if ( SCALAR_IS_ONE(beta) ) NBmm0 = ATL_TargetFn(NBmm_b1);
      else if ( SCALAR_IS_ZERO(beta) ) NBmm0 = ATL_TargetFn(NBmm_b0);
      else NBmm0 = ATL_TargetFn(NBmm_bX);
   }
   else
   {
      NBmm0 = ATL_TargetFn(NBmm_b0);
      vC = Mjoin(simple_malloc,PHSA_FN)(memBlob,
                                        ATL_Cachelen + ATL_MulBySize(NBNB));
      if (!vC) return(-1);
      pC = ATL_AlignPtr(vC);
      if ( SCALAR_IS_ONE(beta) )
         putblk = ATL_TargetFn(Mjoin3(PATL,putblk_b1,PHSA_FN));
      else if ( SCALAR_IS_ZERO(beta) )
         putblk = ATL_TargetFn(Mjoin3(PATL,putblk_b0,PHSA_FN));
      else if ( SCALAR_IS_NONE(beta) )
         putblk = ATL_TargetFn(Mjoin3(PATL,putblk_bn1,PHSA_FN));
      else putblk = ATL_TargetFn(Mjoin3(PATL,putblk_bX,PHSA_FN));
   }
/*
 * Special case if we don't need to copy one or more input matrix
 */
   if (K == NB && TB == AtlasNoTrans && ldb == NB && ATL_DataIsMinAligned(B))
   {
      if (lda == NB && TA == AtlasTrans && SCALAR_IS_ONE(alpha) &&
          ATL_DataIsMinAligned(A))
      {
         i = NBNB;
         pA = (TYPE *) A;
         A = NULL;
         A2blk = ATL_NullFn;
         incA = 0;
      }
      else
      {
         vA = Mjoin(simple_malloc,PHSA_FN)(memBlob,
                                           ATL_Cachelen + ATL_MulBySize(incK));
         if (!vA)
         {
            Mjoin(simple_free,PHSA_FN)(memBlob, vC);
            return(-1);
         }
         pA = ATL_AlignPtr(vA);
         if (TA == AtlasNoTrans)
         {
            incA = NB;
            if ( SCALAR_IS_ONE(alpha) )
               A2blk = ATL_TargetFn(Mjoin3(PATL,row2blkT_a1,PHSA_FN));
            else A2blk = ATL_TargetFn(Mjoin3(PATL,row2blkT_aX,PHSA_FN));
         }
         else
         {
            incA = ATL_MulByNB(lda);
            if ( SCALAR_IS_ONE(alpha) )
               A2blk = ATL_TargetFn(Mjoin3(PATL,col2blk_a1,PHSA_FN));
            else A2blk = ATL_TargetFn(Mjoin3(PATL,col2blk_aX,PHSA_FN));
         }
      }
      Mjoin3(PATL,mmIJK2,PHSA_FN)(
         K, nMb, nNb, nKb, ib, jb, kb, alpha, A, lda, pA,
         incA, A2blk, B, beta, C, ldc, pC, putblk, NBmm0);
      if (vA) Mjoin(simple_free,PHSA_FN)(memBlob, vA);
      if (vC) Mjoin(simple_free,PHSA_FN)(memBlob, vC);
      return(0);
   }
   i = ATL_Cachelen + ATL_MulBySize(N*K + incK);
   if (i <= ATL_MaxMalloc) vA = Mjoin(simple_malloc,PHSA_FN)(memBlob, i);
   if (!vA)
   {
      if (TA == AtlasNoTrans && TB == AtlasNoTrans)
      {
         if (vC) Mjoin(simple_free,PHSA_FN)(memBlob, vC);
         return(1);
      }
      if (jb) n = nNb + 1;
      else n = nNb;
      for (j=2; !vA; j++)
      {
         k = n / j;
         if (k < 1) break;
         if (k*j < n) k++;
         h = ATL_Cachelen + ATL_MulBySize((k+1)*incK);
         if (h <= ATL_MaxMalloc) vA = Mjoin(simple_malloc,PHSA_FN)(memBlob, h);
      }
      if (!vA)
      {
         if (vC) Mjoin(simple_free,PHSA_FN)(memBlob, vC);
         return(-1);
      }
      n = ATL_MulByNB(k);
      jb2 = 0;
   }
   else
   {
      jb2 = jb;
      k = nNb;
      n = N;
   }
   pA = ATL_AlignPtr(vA);
   if (TB == AtlasNoTrans)
   {
      incB = ldb*n;
      if ( SCALAR_IS_ONE(alpha) )
         B2blk = ATL_TargetFn(Mjoin3(PATL,col2blk2_a1,PHSA_FN));
      else B2blk = ATL_TargetFn(Mjoin3(PATL,col2blk2_aX,PHSA_FN));
   }
   else
   {
      incB = n;
      if ( SCALAR_IS_ONE(alpha) )
         B2blk = ATL_TargetFn(Mjoin3(PATL,row2blkT2_a1,PHSA_FN));
      else B2blk = ATL_TargetFn(Mjoin3(PATL,row2blkT2_aX,PHSA_FN));
   }
   if (TA == AtlasNoTrans)
   {
      incA = NB;
      A2blk = ATL_TargetFn(Mjoin3(PATL,row2blkT_a1,PHSA_FN));
   }
   else
   {
      incA = ATL_MulByNB(lda);
      A2blk = ATL_TargetFn(Mjoin3(PATL,col2blk_a1,PHSA_FN));
   }
   incC = ldc*n;
   pB = pA + incK;

   do
   {
      if (TB == AtlasNoTrans) MAT2BLK_ICALL(B2blk, K, n, B, ldb, pB, alpha);
      else MAT2BLK_ICALL(B2blk, n, K, B, ldb, pB, alpha);
      Mjoin3(PATL,mmIJK2,PHSA_FN)(
         K, nMb, k, nKb, ib, jb2, kb, alpha, A, lda, pA,
         incA, A2blk, pB, beta, C, ldc, pC, putblk, NBmm0);
      N -= n;
      nNb -= k;
      if (N < n)
      {
         jb2 = jb;
         n = N;
         k = nNb;
      }
      C += incC;
      B += incB;
      if (!putblk) pC = C;
   }
   while (N);

   if (vC) Mjoin(simple_free,PHSA_FN)(memBlob, vC);
   Mjoin(simple_free,PHSA_FN)(memBlob, vA);
   return(0);
}

