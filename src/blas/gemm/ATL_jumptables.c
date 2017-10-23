/*
 *             Automatically Tuned Linear Algebra Software v3.10.3
 *                    (C) Copyright 2017 Parmance
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

#define ATL_no_icalls
#include "atlas_misc.h"
#include "atlas_lvl3.h"

#ifdef DIRECTHSA

#if defined(GEN_JT_MAT2BLK)
HSA_FUNCTION
void Mjoin3(PATL,icall_MAT2BLK,PHSA)(
   MAT2BLK mat2blk, const int M, const int N, const TYPE *A,
   const int lda, TYPE *V, const SCALAR alpha0)
{
#undef DO_CALL
#define DO_CALL(fn)                             \
   fn(M, N, A, lda, V, alpha0);                 \
   break

   switch (mat2blk)
   {
   default:
   case MAT2BLK_null:
      ATL_assert (0);
      break;
   case Mjoin4(PATL,col2blk_a1,PHSA,_fnid):
      DO_CALL(Mjoin3(PATL,col2blk_a1,PHSA));
   case Mjoin4(PATL,col2blk_aX,PHSA,_fnid):
      DO_CALL(Mjoin3(PATL,col2blk_aX,PHSA));
   case Mjoin4(PATL,col2blk2_a1,PHSA,_fnid):
      DO_CALL(Mjoin3(PATL,col2blk2_a1,PHSA));
   case Mjoin4(PATL,col2blk2_aX,PHSA,_fnid):
      DO_CALL(Mjoin3(PATL,col2blk2_aX,PHSA));
   case Mjoin4(PATL,row2blkT_a1,PHSA,_fnid):
      DO_CALL(Mjoin3(PATL,row2blkT_a1,PHSA));
   case Mjoin4(PATL,row2blkT_aX,PHSA,_fnid):
      DO_CALL(Mjoin3(PATL,row2blkT_aX,PHSA));
   case Mjoin4(PATL,row2blkT2_a1,PHSA,_fnid):
      DO_CALL(Mjoin3(PATL,row2blkT2_a1,PHSA));
   case Mjoin4(PATL,row2blkT2_aX,PHSA,_fnid):
      DO_CALL(Mjoin3(PATL,row2blkT2_aX,PHSA));
   }
}

#elif defined(GEN_JT_MAT2BLK2)
HSA_FUNCTION
void Mjoin3(PATL,icall_MAT2BLK2,PHSA)(
   MAT2BLK2 mat2blk2, const int M, const int N, const float alpha,
   const float *A, const int lda, float *C, const int ldc)
{
#undef DO_CALL
#define DO_CALL(fn)                             \
   fn(M, N, alpha, A, lda, C, ldc);             \
   break

   switch (mat2blk2)
   {
   default:
   case MAT2BLK2_null:
      ATL_assert (0);
      break;
   case Mjoin4(PATL,gemoveT,PHSA,_fnid):
      DO_CALL(Mjoin3(PATL,gemoveT,PHSA));
   case Mjoin4(PATL,gemove,PHSA,_fnid):
      DO_CALL(Mjoin3(PATL,gemove,PHSA));
   }
}

#elif defined(GEN_JT_NBMM0)
HSA_FUNCTION
void Mjoin3(PATL,icall_NBMM0,PHSA)(
   NBMM0 NBmm0, const int M, const int N, const int K, const TYPE alpha,
   const TYPE* A, const int lda, const TYPE* B, const int ldb,
   const TYPE beta, TYPE* C, const int ldc) {

#undef DO_CALL
#define DO_CALL(fn)                                  \
   fn(M, N, K, alpha, A, lda, B, ldb, beta, C, ldc); \
   break

   switch (NBmm0)
   {
   default:
   case NBMM0_null:
      ATL_assert (0);
   case Mjoin(NBmm_b0,_fnid):
      DO_CALL(NBmm_b0);
   case Mjoin(NBmm_b1,_fnid):
      DO_CALL(NBmm_b1);
   case Mjoin(NBmm_bX,_fnid):
      DO_CALL(NBmm_bX);
   case Mjoin4(PATL,pNBmm_b0,PHSA,_fnid):
      DO_CALL(Mjoin3(PATL,pNBmm_b0,PHSA));
   case Mjoin4(PATL,pNBmm_b1,PHSA,_fnid):
      DO_CALL(Mjoin3(PATL,pNBmm_b1,PHSA));
   case Mjoin4(PATL,pNBmm_bX,PHSA,_fnid):
      DO_CALL(Mjoin3(PATL,pNBmm_bX,PHSA));
   case Mjoin4(PATL,pMBmm_b0,PHSA,_fnid):
      DO_CALL(Mjoin3(PATL,pMBmm_b0,PHSA));
   case Mjoin4(PATL,pMBmm_b1,PHSA,_fnid):
      DO_CALL(Mjoin3(PATL,pMBmm_b1,PHSA));
   case Mjoin4(PATL,pMBmm_bX,PHSA,_fnid):
      DO_CALL(Mjoin3(PATL,pMBmm_bX,PHSA));
   case Mjoin4(PATL,pKBmm,PHSA,_fnid):
      DO_CALL(Mjoin3(PATL,pKBmm,PHSA));
   }
}

#elif defined(GEN_JT_PUTBLK)
HSA_FUNCTION
void Mjoin3(PATL,icall_PUTBLK,PHSA)(
   PUTBLK putblk, int M, int N, TYPE *V, TYPE *C, int ldc, const SCALAR beta0)
{
#undef DO_CALL
#define DO_CALL(fn)                             \
   fn(M, N, V, C, ldc, beta0);                  \
   break

   switch (putblk)
   {
   default:
   case PUTBLK_null:
      ATL_assert (0);
   case Mjoin4(PATL,putblk_b0,PHSA,_fnid):
      DO_CALL(Mjoin3(PATL,putblk_b0,PHSA));
   case Mjoin4(PATL,putblk_b1,PHSA,_fnid):
      DO_CALL(Mjoin3(PATL,putblk_b1,PHSA));
   case Mjoin4(PATL,putblk_bn1,PHSA,_fnid):
      DO_CALL(Mjoin3(PATL,putblk_bn1,PHSA));
   case Mjoin4(PATL,putblk_bX,PHSA,_fnid):
      DO_CALL(Mjoin3(PATL,putblk_bX,PHSA));
   }
}

#elif defined(GEN_JT_MMINTR)
HSA_FUNCTION
int Mjoin3(PATL,icall_MMINTR,PHSA)(
   MMINTR mmintr, MemBlob* memBlob,
   const enum ATLAS_TRANS TA, const enum ATLAS_TRANS TB,
   const int M, const int N, const int K, const SCALAR alpha,
   const TYPE *A, const int lda, const TYPE *B, const int ldb,
   const SCALAR beta, TYPE *C, const int ldc)
{
#undef DO_CALL
#define DO_CALL(fn)                                                     \
   return fn(memBlob, TA, TB, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc)

   switch (mmintr)
   {
   default:
   case MMINTR_null:
      ATL_assert (0);
      return -1;
   case Mjoin4(PATL,mmIJK,PHSA,_fnid):
      DO_CALL(Mjoin3(PATL,mmIJK,PHSA));
   case Mjoin4(PATL,mmJIK,PHSA,_fnid):
      DO_CALL(Mjoin3(PATL,mmJIK,PHSA));
   case Mjoin4(PATL,mmJITcp,PHSA,_fnid):
      DO_CALL(Mjoin3(PATL,mmJITcp,PHSA));
   case Mjoin4(PATL,NCmmIJK,PHSA,_fnid):
      DO_CALL(Mjoin3(PATL,NCmmIJK,PHSA));
   case Mjoin4(PATL,NCmmJIK,PHSA,_fnid):
      DO_CALL(Mjoin3(PATL,NCmmJIK,PHSA));
#ifdef USERGEMM
   case Mjoin3(PATL,usergemm_wrapper,_fnid):
      DO_CALL(Mjoin(PATL,usergemm_wrapper));
#endif
   }
}

#elif defined(GEN_JT_GEADD)
HSA_FUNCTION
void Mjoin3(PATL,icall_GEADD,PHSA)(
   GEADD geadd, const int M, const int N, const SCALAR alpha, const TYPE *A,
   const int lda, const SCALAR beta, TYPE *C, const int ldc)
{
#undef DO_CALL
#define DO_CALL(fn)                                 \
   fn(M, N, alpha, A, lda, beta, C, ldc);           \
   break

   switch (geadd)
   {
   default:
   case GEADD_null:
      ATL_assert (0);
      break;
   case Mjoin5(PATL,geadd,_a1_b0,PHSA,_fnid):
      DO_CALL(Mjoin4(PATL,geadd,_a1_b0,PHSA));
   case Mjoin5(PATL,geadd,_a1_b1,PHSA,_fnid):
      DO_CALL(Mjoin4(PATL,geadd,_a1_b1,PHSA));
   case Mjoin5(PATL,geadd,_a1_bX,PHSA,_fnid):
      DO_CALL(Mjoin4(PATL,geadd,_a1_bX,PHSA));
   case Mjoin5(PATL,geadd,_aX_b0,PHSA,_fnid):
      DO_CALL(Mjoin4(PATL,geadd,_aX_b0,PHSA));
   case Mjoin5(PATL,geadd,_aX_b1,PHSA,_fnid):
      DO_CALL(Mjoin4(PATL,geadd,_aX_b1,PHSA));
   case Mjoin5(PATL,geadd,_aX_bX,PHSA,_fnid):
      DO_CALL(Mjoin4(PATL,geadd,_aX_bX,PHSA));
   }
}

#endif /* GEN_JT_* */
#endif /* DIRECTHSA */

#undef DO_CALL

