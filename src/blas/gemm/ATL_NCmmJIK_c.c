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
#  define ATL_no_icalls /* DEVTEMP */
#  define HSADECLS
#endif

#include "atlas_misc.h"
#include Mstr(Mjoin(Mjoin(atlas_,PRE),NCmm.h))
#include "atlas_NCmm.h"
#include "atlas_lvl3.h"
#include "atlas_malloc.h"

#define NCMM_ICALL Mjoin3(PATL,icall_site_NCMM,HSADECL)
#define GEADD_ICALL Mjoin3(PATL,icall_site_GEADD,HSADECL)
#include "ATL_indir_call.c"

HSA_FUNCTION
int Mjoin3(PATL,NCmmJIK_c,PHSA_FN)
   (MemBlob* memBlob,
    const enum ATLAS_TRANS TA, const enum ATLAS_TRANS TB,
    const int M, const int N, const int K, const SCALAR alpha,
    const TYPE *A, const int lda0, const TYPE *B, const int ldb0,
    const SCALAR beta, TYPE *C, const int ldc0)
/*
 * JIK loop-ordered matmul with no matrix copy
 */
{
   size_t incAk, incAm, incAn, incBk, incBm, incBn;
   const size_t lda=lda0, ldb=ldb0, ldc=ldc0;
   const int Mb = M / MB, Nb = N / NB, Kb = K / KB;
   const int mr = M - Mb*MB, nr = N - Nb*NB, kr = K - Kb*KB;
   #define incCm MB
   const size_t incCn = ldc*NB - M + mr;
   int i, j, k;
   const TYPE *a=A, *b=B;
   TYPE *c=C;
   TYPE btmp;
   void *vp;
   TYPE *cp;
   GEADD geadd = ATL_NullFn;
   NCMM mm_bX, mm_b1, mmcu, mm_fixedKcu;
   mm_bX = mm_b1 = mmcu = mm_fixedKcu = ATL_NullFn;

   if (TA == AtlasNoTrans)
   {
      if (TB == AtlasNoTrans)
      {
         mm_bX = ATL_TargetFn(Mjoin5(NCmm0,NN,0x0x0,_a1_b0,PHSA_FN));
         mm_b1 = ATL_TargetFn(Mjoin5(NCmm0,NN,0x0x0,_a1_b1,PHSA_FN));
         mm_fixedKcu =
            ATL_TargetFn(Mjoin5(NCmm00,Mjoin(0x0x,KB),NN,0x0x0_aX_bX,PHSA_FN));
         mmcu = ATL_TargetFn(Mjoin5(NCmm00,0x0x0,NN,0x0x0_aX_bX,PHSA_FN));
      }
      else
      {
         mm_bX = ATL_TargetFn(Mjoin5(NCmm0,NT,0x0x0,_a1_b0,PHSA_FN));
         mm_b1 = ATL_TargetFn(Mjoin5(NCmm0,NT,0x0x0,_a1_b1,PHSA_FN));
         mm_fixedKcu =
            ATL_TargetFn(Mjoin5(NCmm00,Mjoin(0x0x,KB),NT,0x0x0_aX_bX,PHSA_FN));
         mmcu = ATL_TargetFn(Mjoin5(NCmm00,0x0x0,NT,0x0x0_aX_bX,PHSA_FN));
      }
      incAk = lda * KB;
      incAm = MB - Kb * incAk;
      incAn = -Mb * MB;
   }
   else
   {
      if (TB == AtlasNoTrans)
      {
         mm_bX = ATL_TargetFn(Mjoin5(NCmm0,TN,0x0x0,_a1_b0,PHSA_FN));
         mm_b1 = ATL_TargetFn(Mjoin5(NCmm0,TN,0x0x0,_a1_b1,PHSA_FN));
         mm_fixedKcu =
            ATL_TargetFn(Mjoin5(NCmm00,Mjoin(0x0x,KB),TN,0x0x0_aX_bX,PHSA_FN));
         mmcu = ATL_TargetFn(Mjoin5(NCmm00,0x0x0,TN,0x0x0_aX_bX,PHSA_FN));
      }
      else
      {
         mm_bX = ATL_TargetFn(Mjoin5(NCmm0,TT,0x0x0,_a1_b0,PHSA_FN));
         mm_b1 = ATL_TargetFn(Mjoin5(NCmm0,TT,0x0x0,_a1_b1,PHSA_FN));
         mm_fixedKcu =
            ATL_TargetFn(Mjoin5(NCmm00,Mjoin(0x0x,KB),TT,0x0x0_aX_bX,PHSA_FN));
         mmcu = ATL_TargetFn(Mjoin5(NCmm00,0x0x0,TT,0x0x0_aX_bX,PHSA_FN));
      }
      incAk = KB;
      incAm = lda*MB - Kb*KB;
      incAn = -lda*MB*Mb;
   }
   if (TB == AtlasNoTrans)
   {
      incBk = KB;
      incBm = -KB*Kb;
      incBn = ldb*NB;
   }
   else
   {
      incBk = KB*ldb;
      incBm = -Kb * incBk;
      incBn = NB;
   }

   if (alpha == ATL_rone)
   {
      if (beta == ATL_rzero)
         geadd = ATL_TargetFn(Mjoin5(PATL,geadd,_a1,_b0,PHSA_FN));
      else if (beta == ATL_rone)
         geadd = ATL_TargetFn(Mjoin5(PATL,geadd,_a1,_b1,PHSA_FN));
      else geadd = ATL_TargetFn(Mjoin5(PATL,geadd,_a1,_bX,PHSA_FN));
   }
   else if (beta == ATL_rzero)
      geadd = ATL_TargetFn(Mjoin5(PATL,geadd,_aX,_b0,PHSA_FN));
   else if (beta == ATL_rone)
      geadd = ATL_TargetFn(Mjoin5(PATL,geadd,_aX,_b1,PHSA_FN));
   else geadd = ATL_TargetFn(Mjoin5(PATL,geadd,_aX,_bX,PHSA_FN));
   vp = Mjoin(simple_malloc,PHSA_FN)(memBlob,
                                     ATL_Cachelen + ATL_MulBySize(MB * NB));
   ATL_assert(vp);
   cp = ATL_AlignPtr(vp);
   if (mr || nr || kr) for (j=MB*NB, i=0; i != j; i++) cp[i] = ATL_rzero;

   for (j=Nb; j; j--, a += incAn, b += incBn, c += incCn)
   {
      for (i=Mb; i; i--, a += incAm, b += incBm, c += incCm)
      {
         if (Kb)
         {
            NCMM_ICALL(mm_bX, MB, NB, KB, ATL_rone, a, lda, b, ldb, ATL_rzero,
                       cp, MB);
            a += incAk;  b += incBk;
            for (k=Kb-1; k; k--, a += incAk, b += incBk)
               NCMM_ICALL(mm_b1, MB, NB, KB, ATL_rone, a, lda, b, ldb,
                          ATL_rone, cp, MB);
            if (kr)
               NCMM_ICALL(mmcu, MB, NB, kr, ATL_rone, a, lda, b, ldb,
                          ATL_rone, cp, MB);
         }
         else if (kr)
         {
            /* kill NaN/INF from last time */
            Mjoin3(PATL,zero,PHSA_FN)(MB*NB, cp, 1);
            NCMM_ICALL(mmcu, MB, NB, kr, ATL_rone, a, lda, b, ldb, ATL_rzero,
                       cp, MB);
         }
         GEADD_ICALL(geadd, MB, NB, alpha, cp, MB, beta, c, ldc);
      }
   }
   if (mr && N != nr)
      ATL_assert(Mjoin3(PATL,NCmmIJK,PHSA_FN)(
                    memBlob, TA, TB, mr, N-nr, K, alpha,
                    A+Mb*(incAm+Kb*incAk), lda, B, ldb,
                    beta, C+Mb*MB, ldc) ==0);
   if (nr)
   {
      for (i=Mb; i; i--, a += incAm, b += incBm, c += incCm)
      {
         /* kill NaN and INF from last time */
         Mjoin3(PATL,zero,PHSA_FN)(MB*nr, cp, 1);
         if (Kb)
         {
            NCMM_ICALL(mm_fixedKcu, MB, nr, KB, ATL_rone, a, lda, b, ldb,
                       ATL_rzero,  cp, MB);
            a += incAk;  b += incBk;
            for (k=Kb-1; k; k--, a += incAk, b += incBk)
               NCMM_ICALL(mm_fixedKcu, MB, nr, KB, ATL_rone, a, lda, b, ldb,
                          ATL_rone, cp, MB);
            if (kr)
               NCMM_ICALL(mmcu, MB, nr, kr, ATL_rone, a, lda, b, ldb, ATL_rone,
                          cp, MB);
         }
         else if (kr)
            NCMM_ICALL(mmcu, MB, nr, kr, ATL_rone, a, lda, b, ldb, ATL_rzero,
                       cp, MB);
         GEADD_ICALL(geadd, MB, nr, alpha, cp, MB, beta, c, ldc);
      }
      if (mr)  /* cleanup small mr x nr block of C */
      {
         c = C + Mb*MB + ldc*Nb*NB;
         a = A + Mb*(incAm+Kb*incAk);
         b = B + Nb*( incBn+(Mb*(incBm+Kb*incBk)) );
         /* kill NaN and INF from last time */
         Mjoin3(PATL,zero,PHSA_FN)(MB*nr, cp, 1);
         if (Kb)
         {
            NCMM_ICALL(mm_fixedKcu, mr, nr, KB, ATL_rone, a, lda, b, ldb,
                       ATL_rzero, cp, MB);
            a += incAk;  b += incBk;
            for (k=Kb-1; k; k--, a += incAk, b += incBk)
               NCMM_ICALL(mm_fixedKcu, mr, nr, KB, ATL_rone, a, lda, b, ldb,
                          ATL_rone, cp, MB);
            if (kr)
               NCMM_ICALL(mmcu, mr, nr, kr, ATL_rone, a, lda, b, ldb, ATL_rone,
                          cp, MB);
         }
         else if (kr)
            NCMM_ICALL(mmcu, mr, nr, kr, ATL_rone, a, lda, b, ldb, ATL_rzero,
                       cp, MB);
         GEADD_ICALL(geadd, mr, nr, alpha, cp, MB, beta, c, ldc);
      }
   }
   Mjoin(simple_free,PHSA_FN)(memBlob, vp);
   return(0);
}

#ifdef DIRECTHSA

typedef struct Mjoin(PATL,NCmmJIK_c_args_s) {
   MemBlob* memBlob;
   const enum ATLAS_TRANS TA;
   const enum ATLAS_TRANS TB;
   const int M;
   const int N;
   const int K;
   const SCALAR alpha;
   const TYPE *A;
   const int lda0;
   const TYPE *B;
   const int ldb0;
   const SCALAR beta;
   TYPE *C;
   const int ldc0;
   int retval;
} Mjoin(PATL,NCmmJIK_c_args_t);

HSA_KERNEL
void Mjoin3(PATL,NCmmJIK_c,_kernel)(Mjoin(PATL,NCmmJIK_c_args_t)* args)
{
   MemBlob* memBlob = args->memBlob;
   const enum ATLAS_TRANS TA = args->TA;
   const enum ATLAS_TRANS TB = args->TB;
   const int M = args->M;
   const int N = args->N;
   const int K = args->K;
   const SCALAR alpha = args->alpha;
   const TYPE *A = args->A;
   const int lda0 = args->lda0;
   const TYPE *B = args->B;
   const int ldb0 = args->ldb0;
   const SCALAR beta = args->beta;
   TYPE *C = args->C;
   const int ldc0 = args->ldc0;

   args->retval = Mjoin3(PATL,NCmmJIK_c,PHSA_FN)(
      memBlob, TA, TB, M, N, K, alpha, A, lda0, B, ldb0, beta, C, ldc0);
}

int Mjoin3(PATL,NCmmJIK_c,PHSA)
   (MemBlob* memBlob,
    const enum ATLAS_TRANS TA, const enum ATLAS_TRANS TB,
    const int M, const int N, const int K, const SCALAR alpha,
    const TYPE *A, const int lda0, const TYPE *B, const int ldb0,
    const SCALAR beta, TYPE *C, const int ldc0)
{
   Mjoin(PATL,NCmmJIK_c_args_t) args = {
      memBlob, TA, TB, M, N, K, alpha, A, lda0, B, ldb0, beta, C, ldc0, -1 };
   HSA_LAUNCH(Mjoin3(PATL,NCmmJIK_c,_kernel), &args);
   return args.retval;
}

#endif /* DIRECTHSA */
