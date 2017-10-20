#include "atlas_asm.h"
#ifndef ATL_SSE3
   #error "This routine requires SSE3!"
#endif
/*
 * This kernel assumes (lda%2 = 0) && X%16 == A%16
 */
#define II      %rdi
#define JJ      %rsi
#define pA5     %rdx
#define lda     %rcx
#define pX      %rbx
#define pY      %rbp
#define Mr      %rax  /* 0 or 1 of iteration left at end of vector loop? */
#define mlda    %r8
#define mlda5   %r9
#define lda3    %r10
#define mlda3   %r11
#define lda5    %r12
#define II0     %r13
#define lda7    %r14
#define pAn     %r15

#define rX0     %xmm0
#define rA0     %xmm1
#define rY0     %xmm2
#define rY1     %xmm3
#define rY2     %xmm4
#define rY3     %xmm5
#define rY4     %xmm6
#define rY5     %xmm7
#define rY6     %xmm8
#define rY7     %xmm9
#define rY8     %xmm10
#define rY9     %xmm11
#define rY10    %xmm12
#define rY11    %xmm13
#define rY12    %xmm14
#define rY13    %xmm15

#define MOVAPD  movaps
#define MOVSD movsd
#ifdef BETAX
   #define BETAOFF -72
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
 * Swap M & N for backwards compatibility
 */
   xchg %rdi, %rsi
/*
 * Save callee-saved registers, and load in-memory parameters to registers
 */
   movq %rbx, -8(%rsp)
   movq %rbp, -16(%rsp)
   movq %r12, -24(%rsp)
   movq %r13, -32(%rsp)
   movq %r14, -40(%rsp)
   movq %r15, -48(%rsp)

   movq 8(%rsp), pY             /* pY = Y */
   mov  %r8, pX
#ifdef BETAX
   unpcklpd  %xmm1, %xmm1       /* xmm1 = {beta, beta} */
   movapd %xmm1, BETAOFF(%rsp)  /* save beta for later use */
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
   shl $3, lda                  /* lda *= sizeof */
   lea  (lda,lda,4), mlda5      /* mlda5 = 5*lda */
   mov  mlda5, lda5
   add  mlda5, pA5              /* pA5 = pA0 + 5*lda */
   neg  mlda5                   /* mlda5 = -5*lda */
   mov  lda, mlda               /* mlda = lda */
   lea  (mlda5,lda,2), mlda3    /* mlda3 = -5*lda+2*lda = -3*lda */
   neg  mlda                    /* mlda = -lda */
   lea  (lda,lda,2), lda3       /* lda3 = 3*lda */
   lea  (lda3,lda,4), lda7      /* lda7 = 3*lda+4*lda = 7*lda */

   shl  $3, JJ                  /* JJ = N * sizeof */
   lea  (pY, JJ), pY            /* Y += N */
   neg  JJ                      /* JJ = -N*sizeof */

   NLOOP:
      lea (pA5, lda7,2), pAn
      mov       II0, II
/*
 *    First iteration peeled to handle misalignment and zeroing of Y regs
 */
/*
 *    If A is not 16-byte aligned, peel one iteration to force 16-byte alignment
 */
      bt $3, pA5
      jc  PEEL_SCALAR
/*
 *    This code entered when A already 16-byte aligned, peels 1st 2 scalar iters
 */
      MOVAPD (pX,II), rY0
      MOVAPD rY0, rY1
      mulpd  (pA5,mlda5), rY0
      MOVAPD rY1, rY2
      mulpd  (pA5,mlda,4), rY1
      MOVAPD rY2, rY3
      mulpd  (pA5,mlda3), rY2
      MOVAPD rY3, rY4
      mulpd  (pA5,mlda,2), rY3
      MOVAPD rY4, rY5
      mulpd  (pA5,mlda), rY4
      MOVAPD rY5, rY6
      mulpd  (pA5), rY5
      MOVAPD rY6, rY7
      mulpd  (pA5,lda), rY6
      MOVAPD rY7, rY8
      mulpd  (pA5,lda,2), rY7
      MOVAPD rY8, rY9
      mulpd  (pA5,lda3), rY8
      MOVAPD rY9, rY10
      mulpd  (pA5,lda,4), rY9
      MOVAPD rY10, rY11
      mulpd  (pA5,lda5), rY10
      MOVAPD rY11, rY12
      mulpd  (pA5,lda3,2), rY11
      MOVAPD rY12, rY13
      mulpd  (pA5,lda7), rY12
      mulpd  (pA5,lda,8), rY13
      add    $16, pA5
      MOVAPD 16(pX,II), rX0
      add    $32, II
      MLOOP:
         MOVAPD  (pA5,mlda5), rA0
         mulpd   rX0, rA0
         addpd   rA0, rY0
         MOVAPD  (pA5,mlda,4), rA0
         mulpd   rX0, rA0
         addpd   rA0, rY1
         MOVAPD  (pA5,mlda3), rA0
         mulpd   rX0, rA0
         addpd   rA0, rY2
         MOVAPD  (pA5,mlda,2), rA0
         mulpd   rX0, rA0
         addpd   rA0, rY3
         MOVAPD  (pA5,mlda), rA0
         mulpd   rX0, rA0
         addpd   rA0, rY4
         MOVAPD  (pA5), rA0
         mulpd   rX0, rA0
         addpd   rA0, rY5
         MOVAPD  (pA5,lda), rA0
         mulpd   rX0, rA0
         addpd   rA0, rY6
         MOVAPD  (pA5,lda,2), rA0
         mulpd   rX0, rA0
         addpd   rA0, rY7
         MOVAPD  (pA5,lda3), rA0
         mulpd   rX0, rA0
         addpd   rA0, rY8
         MOVAPD  (pA5,lda,4), rA0
         mulpd   rX0, rA0
         addpd   rA0, rY9
         MOVAPD  (pA5,lda5), rA0
         mulpd   rX0, rA0
         addpd   rA0, rY10
         MOVAPD  (pA5,lda3,2), rA0
         mulpd   rX0, rA0
         addpd   rA0, rY11
         MOVAPD  (pA5,lda7), rA0
         mulpd   rX0, rA0
         addpd   rA0, rY12
         MOVAPD  (pA5,lda,8), rA0
         mulpd   rX0, rA0
         addpd   rA0, rY13
         add    $16, pA5          /* pA5 += 16 */
         MOVAPD (pX,II), rX0
      add  $16,   II              /* II += 16 */
      jnz MLOOP
/*
 *    Last iteration peeled, to stop preload of rX0
 */
      MOVAPD  (pA5,mlda5), rA0
      mulpd   rX0, rA0
      addpd   rA0, rY0
      MOVAPD  (pA5,mlda,4), rA0
      mulpd   rX0, rA0
      addpd   rA0, rY1
      MOVAPD  (pA5,mlda3), rA0
      mulpd   rX0, rA0
      addpd   rA0, rY2
      MOVAPD  (pA5,mlda,2), rA0
      mulpd   rX0, rA0
      addpd   rA0, rY3
      MOVAPD  (pA5,mlda), rA0
      mulpd   rX0, rA0
      addpd   rA0, rY4
      MOVAPD  (pA5), rA0
      mulpd   rX0, rA0
      addpd   rA0, rY5
      MOVAPD  (pA5,lda), rA0
      mulpd   rX0, rA0
      addpd   rA0, rY6
      MOVAPD  (pA5,lda,2), rA0
      mulpd   rX0, rA0
      addpd   rA0, rY7
      MOVAPD  (pA5,lda3), rA0
      mulpd   rX0, rA0
      addpd   rA0, rY8
      MOVAPD  (pA5,lda,4), rA0
      mulpd   rX0, rA0
      addpd   rA0, rY9
      MOVAPD  (pA5,lda5), rA0
      mulpd   rX0, rA0
      addpd   rA0, rY10
      MOVAPD  (pA5,lda3,2), rA0
      mulpd   rX0, rA0
      addpd   rA0, rY11
      MOVAPD  (pA5,lda7), rA0
      mulpd   rX0, rA0
      addpd   rA0, rY12
      MOVAPD  (pA5,lda,8), rA0
      mulpd   rX0, rA0
      addpd   rA0, rY13
/*
 *    Handle N=1 cleanup if necessary
 */
      bt $0, Mr
      jc SCALAR_LAST_ITER
YSUM:
/*
 *    Sum up and apply Y/BETA
 */
      #ifdef BETAX
         MOVAPD BETAOFF(%rsp), rX0
         MOVAPD (pY,JJ), rA0
         mulpd  rX0, rA0
      #endif
      haddpd  rY1, rY0          /* rY0 = {Y1ab, Y0ab} */
      #ifdef BETA1
         addpd (pY,JJ), rY0
      #elif BETAX
         addpd rA0, rY0
         MOVAPD 16(pY,JJ), rA0
         mulpd  rX0, rA0
      #endif
      haddpd  rY3, rY2          /* rY2 = {Y3ab, Y2ab} */
      #ifdef BETA1
         addpd 16(pY,JJ), rY2
      #elif BETAX
         addpd rA0, rY2
         MOVAPD 32(pY,JJ), rA0
         mulpd  rX0, rA0
      #endif
      haddpd  rY5, rY4          /* rY4 = {Y5ab, Y4ab} */
      #ifdef BETA1
         addpd 32(pY,JJ), rY4
      #elif BETAX
         addpd rA0, rY4
         MOVAPD 48(pY,JJ), rA0
         mulpd  rX0, rA0
      #endif
      haddpd  rY7, rY6          /* rY6 = {Y7ab, Y6ab} */
      #ifdef BETA1
         addpd 48(pY,JJ), rY6
      #elif BETAX
         addpd rA0, rY6
         MOVAPD 64(pY,JJ), rA0
         mulpd  rX0, rA0
      #endif
      haddpd  rY9, rY8          /* rY8 = {Y9ab, Y8ab} */
      #ifdef BETA1
         addpd 64(pY,JJ), rY8
      #elif BETAX
         addpd rA0, rY8
         MOVAPD 80(pY,JJ), rA0
         mulpd  rX0, rA0
      #endif
      haddpd  rY11,rY10         /* rY10= {Y11ab, Y10ab} */
      #ifdef BETA1
         addpd 80(pY,JJ), rY10
      #elif BETAX
         addpd rA0, rY10
         MOVAPD 96(pY,JJ), rA0
         mulpd  rX0, rA0
      #endif
      haddpd  rY13,rY12         /* rY12= {Y13ab, Y12ab} */
      #ifdef BETA1
         addpd 96(pY,JJ), rY12
      #elif BETAX
         addpd rA0, rY12
      #endif
      MOVAPD    rY0, (pY,JJ)
         mov    pAn, pA5
      MOVAPD    rY2, 16(pY,JJ)
         mov    II0, II
      MOVAPD    rY4, 32(pY,JJ)
      MOVAPD    rY6, 48(pY,JJ)
   prefetchnta  112(pY,JJ)
      MOVAPD    rY8, 64(pY,JJ)
   prefetchnta  112+64(pY,JJ)
      MOVAPD    rY10, 80(pY,JJ)
      MOVAPD    rY12, 96(pY,JJ)
   add   $112, JJ
   jnz NLOOP

DONE:
   movq -8(%rsp), %rbx
   movq -16(%rsp), %rbp
   movq -24(%rsp), %r12
   movq -32(%rsp), %r13
   movq -40(%rsp), %r14
   movq -48(%rsp), %r15
   ret
/*
 *    Peel of one iteration of X loop to handle misalign of X & A
 */
PEEL_SCALAR:
      MOVSD  (pX,II), rY0
      MOVAPD rY0, rY1
      mulsd  (pA5,mlda5), rY0
      MOVAPD rY1, rY2
      mulsd  (pA5,mlda,4), rY1
      MOVAPD rY2, rY3
      mulsd  (pA5,mlda3), rY2
      MOVAPD rY3, rY4
      mulsd  (pA5,mlda,2), rY3
      MOVAPD rY4, rY5
      mulsd  (pA5,mlda), rY4
      MOVAPD rY5, rY6
      mulsd  (pA5), rY5
      MOVAPD rY6, rY7
      mulsd  (pA5,lda), rY6
      MOVAPD rY7, rY8
      mulsd  (pA5,lda,2), rY7
      MOVAPD rY8, rY9
      mulsd  (pA5,lda3), rY8
      MOVAPD rY9, rY10
      mulsd  (pA5,lda,4), rY9
      MOVAPD rY10, rY11
      mulsd  (pA5,lda5), rY10
      MOVAPD rY11, rY12
      mulsd  (pA5,lda3,2), rY11
      MOVAPD rY12, rY13
      mulsd  (pA5,lda7), rY12
      mulsd  (pA5,lda,8), rY13
         addq  $8, pA5
      MOVAPD 8(pX,II), rX0
      add    $24, II
      jmp    MLOOP
/*
 *    When our X+N was did not end with 16-byte alignment, we must do one
 *    extra scalar iteration
 */
SCALAR_LAST_ITER:
   movsd   (pX,II), rX0
   movsd   16(pA5,mlda5), rA0
   mulsd   rX0, rA0
   addsd   rA0, rY0
   movsd   16(pA5,mlda,4), rA0
   mulsd   rX0, rA0
   addsd   rA0, rY1
   movsd   16(pA5,mlda3), rA0
   mulsd   rX0, rA0
   addsd   rA0, rY2
   movsd   16(pA5,mlda,2), rA0
   mulsd   rX0, rA0
   addsd   rA0, rY3
   movsd   16(pA5,mlda), rA0
   mulsd   rX0, rA0
   addsd   rA0, rY4
   movsd   16(pA5), rA0
   mulsd   rX0, rA0
   addsd   rA0, rY5
   movsd   16(pA5,lda), rA0
   mulsd   rX0, rA0
   addsd   rA0, rY6
   movsd   16(pA5,lda,2), rA0
   mulsd   rX0, rA0
   addsd   rA0, rY7
   movsd   16(pA5,lda3), rA0
   mulsd   rX0, rA0
   addsd   rA0, rY8
   movsd   16(pA5,lda,4), rA0
   mulsd   rX0, rA0
   addsd   rA0, rY9
   movsd   16(pA5,lda5), rA0
   mulsd   rX0, rA0
   addsd   rA0, rY10
   movsd   16(pA5,lda3,2), rA0
   mulsd   rX0, rA0
   addsd   rA0, rY11
   movsd   16(pA5,lda7), rA0
   mulsd   rX0, rA0
   addsd   rA0, rY12
   movsd   16(pA5,lda,8), rA0
   mulsd   rX0, rA0
   addsd   rA0, rY13
   jmp     YSUM
