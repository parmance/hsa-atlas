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
#include <stdio.h>
#include <stdarg.h>
#ifdef DIRECTHSA
#define ATL_NoICalls
#endif
#include "atlas_misc.h"
#include "atlas_lvl3.h"
#include "atlas_cacheedge.h"

#include "atlas_malloc.h"

#ifdef SMALL_MM
   #define Cgemm Mjoin(PATL,small_mm)
#elif defined(SMALLK_MM)
   extern int ATL_bigmmOutOfMem;
   #define Cgemm Mjoin(PATL,smallK_mm)
#elif defined (BIG_MM)
   int ATL_bigmmOutOfMem=0;
   #define Cgemm Mjoin(PATL,big_mm)
#elif defined (BIGNORK_MM)
   #define Cgemm Mjoin(PATL,bignork_mm)
#elif defined (CRBIG_MM)
   extern int ATL_bigmmOutOfMem;
   #define Cgemm Mjoin(PATL,Mjoin(UPR,big_mm))
#elif defined(FindingCE) || defined(FindingJITCPCE)
   #define Cgemm Mjoin(PATL,FindCE_mm)
#elif defined (ATLGEMM)
#   define Cgemm Mjoin(PATL,gemm)
#elif defined (ATLGEMM_CPU_REF)
#   define Cgemm Mjoin(PATL,gemm_cpu)
#   define ATLGEMM
#elif defined (ALIASED_GEMM)
#   define Cgemm Mjoin(PATL,aliased_gemm)
#endif

/*
 * This is for a gemm where the matrix C can overlap with A and/or B
*/
#ifdef ALIASED_GEMM
   #define CgemmNN Mjoin(PATL,aliased_gemmNN)
   #define CgemmNT Mjoin(PATL,aliased_gemmNT)
   #define CgemmTN Mjoin(PATL,aliased_gemmTN)
   #define CgemmTT Mjoin(PATL,aliased_gemmTT)
   #ifdef TCPLX
      #define CgemmCN Mjoin(PATL,aliased_gemmCN)
      #define CgemmNC Mjoin(PATL,aliased_gemmNC)
      #define CgemmCT Mjoin(PATL,aliased_gemmCT)
      #define CgemmTC Mjoin(PATL,aliased_gemmTC)
      #define CgemmCC Mjoin(PATL,aliased_gemmCC)
   #endif
/*
 * Otherwise, include routines for doing the various transpose cases of gemm.
 * These are included & declared "static void" to encourage inlining
 */
#else
   #define CgemmNN Mjoin4(PATL,GEMM2,NN,PHSA)
   #define CgemmNT Mjoin4(PATL,GEMM2,NT,PHSA)
   #define CgemmTN Mjoin4(PATL,GEMM2,TN,PHSA)
   #define CgemmTT Mjoin4(PATL,GEMM2,TT,PHSA)
   #ifdef TCPLX
      #define CgemmCN Mjoin(Mjoin(PATL,GEMM2),CN)
      #define CgemmNC Mjoin(Mjoin(PATL,GEMM2),NC)
      #define CgemmCT Mjoin(Mjoin(PATL,GEMM2),CT)
      #define CgemmTC Mjoin(Mjoin(PATL,GEMM2),TC)
      #define CgemmCC Mjoin(Mjoin(PATL,GEMM2),CC)
   #endif
   #define ATL_VOID static void

/*
 * Cases for A is NoTranspose
 */
   #define NoTransA_

   #define NoTransB_
   #define Cgemm__ CgemmNN
   #include "ATL_gemmXX.c"
   #undef Cgemm__
   #undef NoTransB_

   #define TransB_
   #define Cgemm__ CgemmNT
   #include "ATL_gemmXX.c"
   #undef Cgemm__
   #undef TransB_

   #ifdef TCPLX
       #define ConjTransB_
       #define Cgemm__ CgemmNC
       #include "ATL_gemmXX.c"
       #undef Cgemm__
       #undef ConjTransB_
   #endif

   #undef NoTransA_

/*
 * Cases for A is ConjTrans
 */
   #ifdef TCPLX
      #define ConjTransA_

      #define NoTransB_
      #define Cgemm__ CgemmCN
      #include "ATL_gemmXX.c"
      #undef Cgemm__
      #undef  NoTransB_

      #define TransB_
      #define Cgemm__ CgemmCT
      #include "ATL_gemmXX.c"
      #undef Cgemm__
      #undef  TransB_

      #define ConjTransB_
      #define Cgemm__ CgemmCC
      #include "ATL_gemmXX.c"
      #undef Cgemm__
      #undef  ConjTransB_

      #undef  ConjTransA_
   #endif

/*
 * Cases for A is transpose
 */
   #define TransA_

   #define NoTransB_
   #define Cgemm__ CgemmTN
   #include "ATL_gemmXX.c"
   #undef Cgemm__
   #undef NoTransB_

   #define TransB_
   #define Cgemm__ CgemmTT
   #include "ATL_gemmXX.c"
   #undef Cgemm__
   #undef TransB_

   #ifdef TCPLX
       #define ConjTransB_
       #define Cgemm__ CgemmTC
       #include "ATL_gemmXX.c"
       #undef Cgemm__
       #undef ConjTransB_
   #endif
   #undef  TransA_

#endif

typedef struct Mjoin(PATL,gemm_args_s)
{
   MemBlob* memBlob;
   const enum ATLAS_TRANS TA;
   const enum ATLAS_TRANS TB;
   const int M;
   const int N;
   const int K;
   const SCALAR alpha;
   const TYPE *A;
   const int lda;
   const TYPE *B;
   const int ldb;
   const SCALAR beta;
   TYPE *C;
   const int ldc;
} Mjoin(PATL,gemm_args_t);

HSA_KERNEL
static void Mjoin3(Cgemm,_kernel,PHSA)(Mjoin(PATL,gemm_args_t)* args)
/*
 * Actual function to do work.
 */
{
   MemBlob* memBlob = args->memBlob;
   const enum ATLAS_TRANS TA = args->TA;
   const enum ATLAS_TRANS TB = args->TB;
   const int M = args->M;
   const int N = args->N;
   const int K = args->K;
   const SCALAR alpha = args->alpha;
   const TYPE *A = args->A;
   const int lda = args->lda;
   const TYPE *B = args->B;
   const int ldb = args->ldb;
   const SCALAR beta = args->beta;
   TYPE *C = args->C;
   const int ldc = args->ldc;

   if (!M  ||  !N) return;  /* quick return */
   if ( SCALAR_IS_ZERO(alpha) || !K)
   {
      #ifdef TREAL
      if (beta == ATL_rzero)
         Mjoin3(PATL,gezero,PHSA)(M, N, C, ldc);
      else if (beta != ATL_rone)
         Mjoin3(PATL,gescal_bX,PHSA)(M, N, beta, C, ldc);
      #else
         if (beta[1] == ATL_rzero)
         {
            if (*beta == ATL_rzero)
               Mjoin3(PATL,gezero,PHSA)(M, N, C, ldc);
            else if (*beta != ATL_rone)
               Mjoin3(PATL,gescal_bXi0,PHSA)(M, N, beta, C, ldc);
         }
         else
            Mjoin3(PATL,gescal_bX,PHSA)(M, N, beta, C, ldc);
      #endif
      return;
   }
   if (TA == AtlasNoTrans)
   {
      if (TB == AtlasNoTrans)
         CgemmNN(memBlob, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
#ifdef TCPLX
      else if (TB == AtlasConjTrans)
         CgemmNC(memBlob, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
#endif
      else
         CgemmNT(memBlob, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
   }
#ifdef TCPLX
   else if (TA == AtlasConjTrans)
   {
      if (TB == AtlasNoTrans)
         CgemmCN(memBlob, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
      else if (TB == AtlasConjTrans)
         CgemmCC(memBlob, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
      else
         CgemmCT(memBlob, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
   }
#endif
   else
   {
      if (TB == AtlasNoTrans)
         CgemmTN(memBlob, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
#ifdef TCPLX
      else if (TB == AtlasConjTrans)
         CgemmTC(memBlob, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
#endif
      else
         CgemmTT(memBlob, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
   }
}


#if defined(ATLGEMM) && defined(DIRECTHSA)
static char Mjoin(Cgemm,MemPool)[ATL_MaxMalloc];
static MemBlob Mjoin(Cgemm,DynMemBlobData) = {
   &Mjoin(Cgemm,MemPool), sizeof(Mjoin(Cgemm,MemPool)), NULL };
static MemBlob* Mjoin(Cgemm,DynMemBlob) = &Mjoin(Cgemm,DynMemBlobData);
#endif

void Cgemm(const enum ATLAS_TRANS TA, const enum ATLAS_TRANS TB,
           const int M, const int N, const int K, const SCALAR alpha,
           const TYPE *A, const int lda, const TYPE *B, const int ldb,
           const SCALAR beta, TYPE *C, const int ldc)
/*
 * Entry to gemm function. Error checks have been done by interface routine
 */
{
#if defined(ATLGEMM) && defined(DIRECTHSA)
    Mjoin(PATL,gemm_args_t) args = {
       Mjoin(Cgemm,DynMemBlob), TA, TB, M, N, K, alpha, A, lda, B, ldb, beta,
       C, ldc };
    HSA_LAUNCH(Mjoin3(Cgemm,_kernel,PHSA), &args);
#else
    Mjoin(PATL,gemm_args_t) args = {
       globalMemBlob, TA, TB, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc };
    Mjoin3(Cgemm,_kernel,PHSA)(&args);
#endif
}

