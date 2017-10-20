/* TODO copyright */

#define ATL_no_icalls
#include "atlas_misc.h"
#include "atlas_lvl3.h"

#ifdef DIRECTHSA

#ifndef HSADECLS /* DEVTEMP */
#define DBGLOG_CALLED_ONCE                       \
   do                                            \
   {                                             \
      static unsigned called_once = 0;           \
      if (!called_once)                          \
      {                                          \
         called_once++;                          \
         printf ("icall-emul: %s\n", __func__);  \
      }                                          \
   }                                             \
   while(0)

#define DBGLOG_CALLED_UNKNOWN                           \
   do                                                   \
   {                                                    \
      static unsigned called_unknown = 0;               \
      if (!called_unknown)                              \
      {                                                 \
         called_unknown++;                              \
         printf("%s: Called unknown fn", __func__);     \
      }                                                 \
   }                                                    \
   while(0)
#else
#define DBGLOG_CALLED_ONCE
#define DBGLOG_CALLED_UNKNOWN
#endif

#if defined(GEN_JT_MAT2BLK)
HSA_FUNCTION
void Mjoin3(PATL,icall_MAT2BLK,HSADECL)(
   MAT2BLK mat2blk, const int M, const int N, const TYPE *A,
   const int lda, TYPE *V, const SCALAR alpha0)
{
#undef DO_CALL
#define DO_CALL(fn)                             \
   fn(M, N, A, lda, V, alpha0);                 \
   break

   DBGLOG_CALLED_ONCE;

   switch (mat2blk)
   {
   default:
      DBGLOG_CALLED_UNKNOWN;
   case MAT2BLK_null:
      ATL_assert (0);
      break;
   case Mjoin4(PATL,col2blk_a1,HSADECL,_fnid):
      DO_CALL(Mjoin3(PATL,col2blk_a1,HSADECL));
   case Mjoin4(PATL,col2blk_aX,HSADECL,_fnid):
      DO_CALL(Mjoin3(PATL,col2blk_aX,HSADECL));
   case Mjoin4(PATL,col2blk2_a1,HSADECL,_fnid):
      DO_CALL(Mjoin3(PATL,col2blk2_a1,HSADECL));
   case Mjoin4(PATL,col2blk2_aX,HSADECL,_fnid):
      DO_CALL(Mjoin3(PATL,col2blk2_aX,HSADECL));
   case Mjoin4(PATL,row2blkT_a1,HSADECL,_fnid):
      DO_CALL(Mjoin3(PATL,row2blkT_a1,HSADECL));
   case Mjoin4(PATL,row2blkT_aX,HSADECL,_fnid):
      DO_CALL(Mjoin3(PATL,row2blkT_aX,HSADECL));
   case Mjoin4(PATL,row2blkT2_a1,HSADECL,_fnid):
      DO_CALL(Mjoin3(PATL,row2blkT2_a1,HSADECL));
   case Mjoin4(PATL,row2blkT2_aX,HSADECL,_fnid):
      DO_CALL(Mjoin3(PATL,row2blkT2_aX,HSADECL));
   }
}

#elif defined(GEN_JT_MAT2BLK2)
HSA_FUNCTION
void Mjoin3(PATL,icall_MAT2BLK2,HSADECL)(
   MAT2BLK2 mat2blk2, const int M, const int N, const float alpha,
   const float *A, const int lda, float *C, const int ldc)
{
#undef DO_CALL
#define DO_CALL(fn)                             \
   fn(M, N, alpha, A, lda, C, ldc);             \
   break

   DBGLOG_CALLED_ONCE;

   switch (mat2blk2)
   {
   default:
      DBGLOG_CALLED_UNKNOWN;
   case MAT2BLK2_null:
      ATL_assert (0);
      break;
   case Mjoin4(PATL,gemoveT,HSADECL,_fnid):
      DO_CALL(Mjoin3(PATL,gemoveT,HSADECL));
   case Mjoin4(PATL,gemove,HSADECL,_fnid):
      DO_CALL(Mjoin3(PATL,gemove,HSADECL));
   }
}

#elif defined(GEN_JT_NBMM0)
HSA_FUNCTION
void Mjoin3(PATL,icall_NBMM0,HSADECL)(
   NBMM0 NBmm0, const int M, const int N, const int K, const TYPE alpha,
   const TYPE* A, const int lda, const TYPE* B, const int ldb,
   const TYPE beta, TYPE* C, const int ldc) {

#undef DO_CALL
#define DO_CALL(fn)                                  \
   fn(M, N, K, alpha, A, lda, B, ldb, beta, C, ldc); \
   break

   DBGLOG_CALLED_ONCE;

   switch (NBmm0)
   {
   default:
      DBGLOG_CALLED_UNKNOWN;
   case NBMM0_null:
      ATL_assert (0);
   case Mjoin(NBmm_b0,_fnid):
      DO_CALL(NBmm_b0);
   case Mjoin(NBmm_b1,_fnid):
      DO_CALL(NBmm_b1);
   case Mjoin(NBmm_bX,_fnid):
      DO_CALL(NBmm_bX);
   case Mjoin4(PATL,pNBmm_b0,HSADECL,_fnid):
      DO_CALL(Mjoin3(PATL,pNBmm_b0,HSADECL));
   case Mjoin4(PATL,pNBmm_b1,HSADECL,_fnid):
      DO_CALL(Mjoin3(PATL,pNBmm_b1,HSADECL));
   case Mjoin4(PATL,pNBmm_bX,HSADECL,_fnid):
      DO_CALL(Mjoin3(PATL,pNBmm_bX,HSADECL));
   case Mjoin4(PATL,pMBmm_b0,HSADECL,_fnid):
      DO_CALL(Mjoin3(PATL,pMBmm_b0,HSADECL));
   case Mjoin4(PATL,pMBmm_b1,HSADECL,_fnid):
      DO_CALL(Mjoin3(PATL,pMBmm_b1,HSADECL));
   case Mjoin4(PATL,pMBmm_bX,HSADECL,_fnid):
      DO_CALL(Mjoin3(PATL,pMBmm_bX,HSADECL));
   case Mjoin4(PATL,pKBmm,HSADECL,_fnid):
      DO_CALL(Mjoin3(PATL,pKBmm,HSADECL));
   }
}

#elif defined(GEN_JT_PUTBLK)
HSA_FUNCTION
void Mjoin3(PATL,icall_PUTBLK,HSADECL)(
   PUTBLK putblk, int M, int N, TYPE *V, TYPE *C, int ldc, const SCALAR beta0)
{
#undef DO_CALL
#define DO_CALL(fn)                             \
   fn(M, N, V, C, ldc, beta0);                  \
   break

   DBGLOG_CALLED_ONCE;

   switch (putblk)
   {
   default:
      DBGLOG_CALLED_UNKNOWN;
   case PUTBLK_null:
      ATL_assert (0);
   case Mjoin4(PATL,putblk_b0,HSADECL,_fnid):
      DO_CALL(Mjoin3(PATL,putblk_b0,HSADECL));
   case Mjoin4(PATL,putblk_b1,HSADECL,_fnid):
      DO_CALL(Mjoin3(PATL,putblk_b1,HSADECL));
   case Mjoin4(PATL,putblk_bn1,HSADECL,_fnid):
      DO_CALL(Mjoin3(PATL,putblk_bn1,HSADECL));
   case Mjoin4(PATL,putblk_bX,HSADECL,_fnid):
      DO_CALL(Mjoin3(PATL,putblk_bX,HSADECL));
   }
}

#elif defined(GEN_JT_MMINTR)
HSA_FUNCTION
int Mjoin3(PATL,icall_MMINTR,HSADECL)(
   MMINTR mmintr, MemBlob* memBlob,
   const enum ATLAS_TRANS TA, const enum ATLAS_TRANS TB,
   const int M, const int N, const int K, const SCALAR alpha,
   const TYPE *A, const int lda, const TYPE *B, const int ldb,
   const SCALAR beta, TYPE *C, const int ldc)
{
#undef DO_CALL
#define DO_CALL(fn)                                                     \
   return fn(memBlob, TA, TB, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc)

   DBGLOG_CALLED_ONCE;

   switch (mmintr)
   {
   default:
      DBGLOG_CALLED_UNKNOWN;
   case MMINTR_null:
      ATL_assert (0);
      return -1;
   case Mjoin4(PATL,mmIJK,HSADECL,_fnid):
      DO_CALL(Mjoin3(PATL,mmIJK,HSADECL));
   case Mjoin4(PATL,mmJIK,HSADECL,_fnid):
      DO_CALL(Mjoin3(PATL,mmJIK,HSADECL));
   case Mjoin4(PATL,mmJITcp,HSADECL,_fnid):
      DO_CALL(Mjoin3(PATL,mmJITcp,HSADECL));
   case Mjoin4(PATL,NCmmIJK,HSADECL,_fnid):
      DO_CALL(Mjoin3(PATL,NCmmIJK,HSADECL));
   case Mjoin4(PATL,NCmmJIK,HSADECL,_fnid):
      DO_CALL(Mjoin3(PATL,NCmmJIK,HSADECL));
#ifdef USERGEMM
   case Mjoin3(PATL,usergemm_wrapper,_fnid):
      DO_CALL(Mjoin(PATL,usergemm_wrapper));
#endif
   }
}

#elif defined(GEN_JT_GEADD)
HSA_FUNCTION
void Mjoin3(PATL,icall_GEADD,HSADECL)(
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
      DBGLOG_CALLED_UNKNOWN;
   case GEADD_null:
      ATL_assert (0);
      break;
   case Mjoin5(PATL,geadd,_a1_b0,HSADECL,_fnid):
      DO_CALL(Mjoin4(PATL,geadd,_a1_b0,HSADECL));
   case Mjoin5(PATL,geadd,_a1_b1,HSADECL,_fnid):
      DO_CALL(Mjoin4(PATL,geadd,_a1_b1,HSADECL));
   case Mjoin5(PATL,geadd,_a1_bX,HSADECL,_fnid):
      DO_CALL(Mjoin4(PATL,geadd,_a1_bX,HSADECL));
   case Mjoin5(PATL,geadd,_aX_b0,HSADECL,_fnid):
      DO_CALL(Mjoin4(PATL,geadd,_aX_b0,HSADECL));
   case Mjoin5(PATL,geadd,_aX_b1,HSADECL,_fnid):
      DO_CALL(Mjoin4(PATL,geadd,_aX_b1,HSADECL));
   case Mjoin5(PATL,geadd,_aX_bX,HSADECL,_fnid):
      DO_CALL(Mjoin4(PATL,geadd,_aX_bX,HSADECL));
   }
}

#endif /* GEN_JT_* */
#endif /* DIRECTHSA */

#undef DO_CALL

