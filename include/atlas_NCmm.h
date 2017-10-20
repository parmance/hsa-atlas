/* TODO copyright */

#ifndef ATLAS_NCMM_H
#define ATLAS_NCMM_H

#include "atlas_misc.h"
#include Mstr(Mjoin(Mjoin(atlas_,PRE),NCmm.h))

#ifndef MB
   #define MB NB
#endif
#ifndef KB
   #define KB NB
#endif
#define NBnam Mjoin(Mjoin(Mjoin(Mjoin(MB,x),NB),x),KB)
#define NCmm0 Mjoin(Mjoin(PATL,JIK),NBnam)
#define NCmm00 Mjoin(PATL,JIK)
void Mjoin5(NCmm0,NN,0x0x0,_a1_b1,HSADECL)
   (const int M, const int N, const int K, const SCALAR alpha, const TYPE *A,
    const int lda, const TYPE *B, const int ldb,
    const SCALAR beta, TYPE *C, const int ldc);
void Mjoin5(NCmm0,NN,0x0x0,_a1_b0,HSADECL)
   (const int M, const int N, const int K, const SCALAR alpha, const TYPE *A,
    const int lda, const TYPE *B, const int ldb,
    const SCALAR beta, TYPE *C, const int ldc);
void Mjoin5(NCmm0,NN,0x0x0,_aX_bX,HSADECL)
   (const int M, const int N, const int K, const SCALAR alpha, const TYPE *A,
    const int lda, const TYPE *B, const int ldb,
    const SCALAR beta, TYPE *C, const int ldc);
void Mjoin5(NCmm0,NN,0x0x0,_aX_b0,HSADECL)
   (const int M, const int N, const int K, const SCALAR alpha, const TYPE *A,
    const int lda, const TYPE *B, const int ldb,
    const SCALAR beta, TYPE *C, const int ldc);
void Mjoin5(NCmm0,NN,0x0x0,_a1_bX,HSADECL)
   (const int M, const int N, const int K, const SCALAR alpha, const TYPE *A,
    const int lda, const TYPE *B, const int ldb,
    const SCALAR beta, TYPE *C, const int ldc);
void Mjoin5(NCmm00,Mjoin(0x0x,KB),NN,0x0x0_aX_bX,HSADECL)
   (const int M, const int N, const int K, const SCALAR alpha, const TYPE *A,
    const int lda, const TYPE *B, const int ldb,
    const SCALAR beta, TYPE *C, const int ldc);
void Mjoin5(NCmm00,0x0x0,NN,0x0x0_aX_bX,HSADECL)
   (const int M, const int N, const int K, const SCALAR alpha, const TYPE *A,
    const int lda, const TYPE *B, const int ldb,
    const SCALAR beta, TYPE *C, const int ldc);

void Mjoin5(NCmm0,NT,0x0x0,_a1_b1,HSADECL)
   (const int M, const int N, const int K, const SCALAR alpha, const TYPE *A,
    const int lda, const TYPE *B, const int ldb,
    const SCALAR beta, TYPE *C, const int ldc);
void Mjoin5(NCmm0,NT,0x0x0,_a1_b0,HSADECL)
   (const int M, const int N, const int K, const SCALAR alpha, const TYPE *A,
    const int lda, const TYPE *B, const int ldb,
    const SCALAR beta, TYPE *C, const int ldc);
void Mjoin5(NCmm0,NT,0x0x0,_aX_bX,HSADECL)
   (const int M, const int N, const int K, const SCALAR alpha, const TYPE *A,
    const int lda, const TYPE *B, const int ldb,
    const SCALAR beta, TYPE *C, const int ldc);
void Mjoin5(NCmm0,NT,0x0x0,_aX_b0,HSADECL)
   (const int M, const int N, const int K, const SCALAR alpha, const TYPE *A,
    const int lda, const TYPE *B, const int ldb,
    const SCALAR beta, TYPE *C, const int ldc);
void Mjoin5(NCmm0,NT,0x0x0,_a1_bX,HSADECL)
   (const int M, const int N, const int K, const SCALAR alpha, const TYPE *A,
    const int lda, const TYPE *B, const int ldb,
    const SCALAR beta, TYPE *C, const int ldc);
void Mjoin5(NCmm00,Mjoin(0x0x,KB),NT,0x0x0_aX_bX,HSADECL)
   (const int M, const int N, const int K, const SCALAR alpha, const TYPE *A,
    const int lda, const TYPE *B, const int ldb,
    const SCALAR beta, TYPE *C, const int ldc);
void Mjoin5(NCmm00,0x0x0,NT,0x0x0_aX_bX,HSADECL)
   (const int M, const int N, const int K, const SCALAR alpha, const TYPE *A,
    const int lda, const TYPE *B, const int ldb,
    const SCALAR beta, TYPE *C, const int ldc);

void Mjoin5(NCmm0,TN,0x0x0,_a1_b1,HSADECL)
   (const int M, const int N, const int K, const SCALAR alpha, const TYPE *A,
    const int lda, const TYPE *B, const int ldb,
    const SCALAR beta, TYPE *C, const int ldc);
void Mjoin5(NCmm0,TN,0x0x0,_a1_b0,HSADECL)
   (const int M, const int N, const int K, const SCALAR alpha, const TYPE *A,
    const int lda, const TYPE *B, const int ldb,
    const SCALAR beta, TYPE *C, const int ldc);
void Mjoin5(NCmm0,TN,0x0x0,_aX_bX,HSADECL)
   (const int M, const int N, const int K, const SCALAR alpha, const TYPE *A,
    const int lda, const TYPE *B, const int ldb,
    const SCALAR beta, TYPE *C, const int ldc);
void Mjoin5(NCmm0,TN,0x0x0,_aX_b0,HSADECL)
   (const int M, const int N, const int K, const SCALAR alpha, const TYPE *A,
    const int lda, const TYPE *B, const int ldb,
    const SCALAR beta, TYPE *C, const int ldc);
void Mjoin5(NCmm0,TN,0x0x0,_a1_bX,HSADECL)
   (const int M, const int N, const int K, const SCALAR alpha, const TYPE *A,
    const int lda, const TYPE *B, const int ldb,
    const SCALAR beta, TYPE *C, const int ldc);
void Mjoin5(NCmm00,Mjoin(0x0x,KB),TN,0x0x0_aX_bX,HSADECL)
   (const int M, const int N, const int K, const SCALAR alpha, const TYPE *A,
    const int lda, const TYPE *B, const int ldb,
    const SCALAR beta, TYPE *C, const int ldc);
void Mjoin5(NCmm00,0x0x0,TN,0x0x0_aX_bX,HSADECL)
   (const int M, const int N, const int K, const SCALAR alpha, const TYPE *A,
    const int lda, const TYPE *B, const int ldb,
    const SCALAR beta, TYPE *C, const int ldc);

void Mjoin5(NCmm0,TT,0x0x0,_a1_b1,HSADECL)
   (const int M, const int N, const int K, const SCALAR alpha, const TYPE *A,
    const int lda, const TYPE *B, const int ldb,
    const SCALAR beta, TYPE *C, const int ldc);
void Mjoin5(NCmm0,TT,0x0x0,_a1_b0,HSADECL)
   (const int M, const int N, const int K, const SCALAR alpha, const TYPE *A,
    const int lda, const TYPE *B, const int ldb,
    const SCALAR beta, TYPE *C, const int ldc);
void Mjoin5(NCmm0,TT,0x0x0,_aX_bX,HSADECL)
   (const int M, const int N, const int K, const SCALAR alpha, const TYPE *A,
    const int lda, const TYPE *B, const int ldb,
    const SCALAR beta, TYPE *C, const int ldc);
void Mjoin5(NCmm0,TT,0x0x0,_aX_b0,HSADECL)
   (const int M, const int N, const int K, const SCALAR alpha, const TYPE *A,
    const int lda, const TYPE *B, const int ldb,
    const SCALAR beta, TYPE *C, const int ldc);
void Mjoin5(NCmm0,TT,0x0x0,_a1_bX,HSADECL)
   (const int M, const int N, const int K, const SCALAR alpha, const TYPE *A,
    const int lda, const TYPE *B, const int ldb,
    const SCALAR beta, TYPE *C, const int ldc);
void Mjoin5(NCmm00,Mjoin(0x0x,KB),TT,0x0x0_aX_bX,HSADECL)
   (const int M, const int N, const int K, const SCALAR alpha, const TYPE *A,
    const int lda, const TYPE *B, const int ldb,
    const SCALAR beta, TYPE *C, const int ldc);
void Mjoin5(NCmm00,0x0x0,TT,0x0x0_aX_bX,HSADECL)
   (const int M, const int N, const int K, const SCALAR alpha, const TYPE *A,
    const int lda, const TYPE *B, const int ldb,
    const SCALAR beta, TYPE *C, const int ldc);

#ifndef ATL_MaxMMalpha
   #define ATL_MaxMMalpha 3
#endif
#ifndef MB
   #define MB NB
#endif
#ifndef KB
   #define KB NB
#endif

#ifdef ATL_no_icalls
typedef enum NCMM_fnid {
   NCMM_null,
   Mjoin(Mjoin5(NCmm00,0x0x0,NN,0x0x0_aX_bX,HSADECL),_fnid),
   Mjoin(Mjoin5(NCmm00,0x0x0,NT,0x0x0_aX_bX,HSADECL),_fnid),
   Mjoin(Mjoin5(NCmm00,0x0x0,TN,0x0x0_aX_bX,HSADECL),_fnid),
   Mjoin(Mjoin5(NCmm00,0x0x0,TT,0x0x0_aX_bX,HSADECL),_fnid),
   Mjoin(Mjoin5(NCmm00,Mjoin(0x0x,KB),NN,0x0x0_aX_bX,HSADECL),_fnid),
   Mjoin(Mjoin5(NCmm00,Mjoin(0x0x,KB),NT,0x0x0_aX_bX,HSADECL),_fnid),
   Mjoin(Mjoin5(NCmm00,Mjoin(0x0x,KB),TN,0x0x0_aX_bX,HSADECL),_fnid),
   Mjoin(Mjoin5(NCmm00,Mjoin(0x0x,KB),TT,0x0x0_aX_bX,HSADECL),_fnid),
   Mjoin(Mjoin5(NCmm0,NN,0x0x0,_a1_b0,HSADECL),_fnid),
   Mjoin(Mjoin5(NCmm0,NN,0x0x0,_a1_b1,HSADECL),_fnid),
   Mjoin(Mjoin5(NCmm0,NN,0x0x0,_a1_bX,HSADECL),_fnid),
   Mjoin(Mjoin5(NCmm0,NN,0x0x0,_aX_b0,HSADECL),_fnid),
   Mjoin(Mjoin5(NCmm0,NN,0x0x0,_aX_bX,HSADECL),_fnid),
   Mjoin(Mjoin5(NCmm0,NT,0x0x0,_a1_b0,HSADECL),_fnid),
   Mjoin(Mjoin5(NCmm0,NT,0x0x0,_a1_b1,HSADECL),_fnid),
   Mjoin(Mjoin5(NCmm0,NT,0x0x0,_a1_bX,HSADECL),_fnid),
   Mjoin(Mjoin5(NCmm0,NT,0x0x0,_aX_b0,HSADECL),_fnid),
   Mjoin(Mjoin5(NCmm0,NT,0x0x0,_aX_bX,HSADECL),_fnid),
   Mjoin(Mjoin5(NCmm0,TN,0x0x0,_a1_b0,HSADECL),_fnid),
   Mjoin(Mjoin5(NCmm0,TN,0x0x0,_a1_b1,HSADECL),_fnid),
   Mjoin(Mjoin5(NCmm0,TN,0x0x0,_a1_bX,HSADECL),_fnid),
   Mjoin(Mjoin5(NCmm0,TN,0x0x0,_aX_b0,HSADECL),_fnid),
   Mjoin(Mjoin5(NCmm0,TN,0x0x0,_aX_bX,HSADECL),_fnid),
   Mjoin(Mjoin5(NCmm0,TT,0x0x0,_a1_b0,HSADECL),_fnid),
   Mjoin(Mjoin5(NCmm0,TT,0x0x0,_a1_b1,HSADECL),_fnid),
   Mjoin(Mjoin5(NCmm0,TT,0x0x0,_a1_bX,HSADECL),_fnid),
   Mjoin(Mjoin5(NCmm0,TT,0x0x0,_aX_b0,HSADECL),_fnid),
   Mjoin(Mjoin5(NCmm0,TT,0x0x0,_aX_bX,HSADECL),_fnid),
} NCMM;

HSA_FUNCTION
void Mjoin3(PATL,icall_NCMM,HSADECL)(
   NCMM NCmm, const int M, const int N, const int K, const SCALAR alpha,
   const TYPE *A, const int lda, const TYPE *B, const int ldb,
   const SCALAR beta, TYPE *C, const int ldc);

#else

typedef void (*NCMM)(const int M, const int N, const int K, const SCALAR alpha,
                     const TYPE *A, const int lda, const TYPE *B, const int ldb,
                     const SCALAR beta, TYPE *C, const int ldc);
#endif


#endif
