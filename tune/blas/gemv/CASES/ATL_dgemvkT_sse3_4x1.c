#include "atlas_asm.h"
#ifndef ATL_SSE3
   #error "This routine requires SSE3!"
#endif
/*
 * Rename for backward compatibility
 */
#ifndef ATL_UGEMVT
   #if BETA0
      #define ATL_UGEMVT ATL_dgemvkT_b0
   #elif defined(BETAX)
      #define ATL_UGEMVT ATL_dgemvkT_bX
   #else
      #ifndef BETA1
         #define BETA1
      #endif
      #define ATL_UGEMVT ATL_dgemvkT_b1
   #endif
#endif
/*
 * ASSUMES: M >= 6, M%2=0, N%4=0, A,X,Y all 16-byte aligned
 */

#define rX0     %xmm0
#define rA00    %xmm1
#define rA01    %xmm2
#define rA02    %xmm3
#define rA03    %xmm4
#define rY0     %xmm5
#define rY1     %xmm6
#define rY2     %xmm7
#define rY3     %xmm8
#define rBETA   %xmm9

#define pA0     %rdx
#define lda     %rcx
#define lda3    %rax
#define pX      %rbx
#define II      %rsi
#define pY      %rbp
#define JJ      %rdi
#define pAn     %r8
#define II0     %r9
#define Mr      %r10 /* 0 or 1 iteration left at end of vector loop? */
/*                      rdi,         rsi,             xmm0
void ATL_UGEMVT(const int M, const int N, const float alpha,
                           rdx            rcx              r8              r9
                const float *A, const int lda, const float *X, const int incX,
                           %xmm1     8(rsp)        16(rsp)
                const float beta, float *Y, const int incY)
*/

        .text
.global ATL_asmdecor(ATL_UGEMVT)
ATL_asmdecor(ATL_UGEMVT):
/*
 * Save callee-saved registers, and load in-memory parameters to registers
 */
   movq %rbx, -8(%rsp)
   movq %rbp, -16(%rsp)

   movq %r8, pX
   movq 8(%rsp), pY
   #ifdef BETAX
      movddup   %xmm1, rBETA
   #endif
/*
 *  Compute the number of remainder iterations after vector loop complete
 */
   xor Mr, Mr                  /* Mr = 0 */
   mov $-1, lda3               /* lda = -1 */
   mov pX, II0                 /* II0 = X */
   lea (pX, II, 8), pX         /* pX pts to end of X */
   bt  $3, pX                  /* CF=1 if pX not aligned */
   cmovc lda3, Mr              /* Mr = (X ends aligned?) 0 : -1 */
   lea (pX, Mr, 8), pX         /* pts to end of aligned space */
   sub pX, II0                 /* II0 = -na*sizeof, na= aligned space elts */
/*
 * Compute derived values
 */
   shl  $3, II                  /* II = M*sizeof */
   movq II, II0
   shl $3, lda                  /* lda *= sizeof */
   lea (lda,lda,2), lda3        /* lda3 = lda*3 */
   shl $3, JJ                   /* JJ = N*sizeof */
   lea (pY, JJ), pY             /* pY += N */
   NLOOP:
      lea (pA0, lda,4), pAn      /* pAn = pA0 + lda*4 */
/*
 *    Peel 1st iteration to zero Y, and then preload all A/B regs
 */
      movaps  (pX,II), rY0
      movaps  rY0, rY1
      mulpd   (pA0), rY0
      movaps  rY1, rY2
      mulpd   (pA0,lda), rY1
      movaps  rY2, rY3
      mulpd   (pA0,lda,2), rY2
      mulpd   (pA0,lda3), rY3

      movaps 16(pA0), rA00
      movaps 16(pA0,lda), rA01
      movaps 16(pA0,lda,2), rA02
      movaps 16(pA0,lda3), rA03
      movaps 16(pX), rX0
      add    $32, pA0
      add    $48, pX
   .align 16
   LOOPM:
      mulpd  rX0, rA00
      addpd  rA00, rY0
      movaps (pA0), rA00

      mulpd  rX0, rA01
      addpd  rA01, rY1
      movaps (pA0,lda), rA01

      mulpd  rX0, rA02
      addpd  rA01, rY2
      movaps (pA0,lda,2), rA02

      mulpd  rX0, rA03
      movaps -16(pX,II), rX0
      addpd  rA03, rY3
      movaps (pA0,lda3), rA03
      add    $16, pA0
      add    $16, II
      jnz    LOOPM
/*
 *    Drain the fetch pipe
 */
      mulpd  rX0, rA00
      addpd  rA00, rY0
      mulpd  rX0, rA01
      addpd  rA01, rY1
      mulpd  rX0, rA02
      addpd  rA01, rY2
      mulpd  rX0, rA03
      addpd  rA03, rY3
      haddpd rY1, rY0
      haddpd rY3, rY2
      #ifdef BETA1
         addps  (pY,JJ), rY0
         addps  16(pY,JJ), rY2
      #elif defined(BETAX)
         movaps (pY,JJ), rA00
         movaps 16(pY,JJ), rA01
         mulps  rBETA, rA00
         mulps  rBETA, rA01
         addps  rA00, rY0
         addps  rA01, rY2
      #endif
      movaps    rY0, (pY,JJ)
      movaps    rY2, 16(pY,JJ)
      movq      pAn, pA0
      movq      II0, II
   add       $16, JJ
   jnz  NLOOP

   movq -8(%rsp), %rbx
   movq -16(%rsp), %rbp
   ret
