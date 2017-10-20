/* TODO copyright */

/* Architecture/target depend indirect call interfaces. For now, only
 * HSA (DirectHSA) target has alternate indirect call
 * implementation. See ATL_jumptables.c for the implementation
 * details.
 *
 * Others have the usual C indirect call via function pointers.
 *
 * This file is meant to be included directly. All functions here are
 * static to encourage inlining.
 */

#ifdef MAT2BLK_ICALL
static void MAT2BLK_ICALL(
   MAT2BLK mat2blk, const int M, const int N, const TYPE *A, const int lda,
   TYPE *V, const SCALAR alpha0)
{
#ifdef DIRECTHSA
   Mjoin3(PATL,icall_MAT2BLK,HSADECL)(mat2blk, M, N, A, lda, V, alpha0);
#else
   mat2blk(M, N, A, lda, V, alpha0);
#endif
}
#endif

#ifdef MAT2BLK2_ICALL
static void MAT2BLK2_ICALL(
   MAT2BLK2 mat2blk2, const int M, const int N, const TYPE alpha,
   const TYPE *A, const int lda, TYPE *C, const int ldc)
{
#ifdef DIRECTHSA
   Mjoin3(PATL,icall_MAT2BLK2,HSADECL)(
      mat2blk2, M, N, alpha, A, lda, C, ldc);
#else
   mat2blk2(M, N, alpha, A, lda, C, ldc);
#endif
}
#endif

#ifdef PUTBLK_ICALL
static void PUTBLK_ICALL(
   PUTBLK putblk, int M, int N, TYPE *V, TYPE *C, int ldc, const SCALAR beta0)
{
#  ifdef DIRECTHSA
   Mjoin3(PATL,icall_PUTBLK,HSADECL)(
      putblk, M, N, V, C, ldc, beta0);
#  else
   putblk(M, N, V, C, ldc, beta0);
#  endif
}
#endif

#ifdef NBMM_ICALL
static void NBMM_ICALL(
   NBMM0 NBmm0, const int M, const int N, const int K, const TYPE alpha,
   const TYPE* A, const int lda, const TYPE* B, const int ldb,
   const TYPE beta, TYPE* C, const int ldc)
{
#ifdef DIRECTHSA
   Mjoin3(PATL,icall_NBMM0,HSADECL)(
      NBmm0, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
#else
   NBmm0(M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
#endif
}
#endif

#ifdef NCMM_ICALL
static void NCMM_ICALL(
   NCMM NCmm, const int M, const int N, const int K, const SCALAR alpha,
   const TYPE *A, const int lda, const TYPE *B, const int ldb,
   const SCALAR beta, TYPE *C, const int ldc)
{
#ifdef DIRECTHSA
   Mjoin3(PATL,icall_NCMM,HSADECL)(
      NCmm, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
#else
   NCmm(M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
#endif
}
#endif

#ifdef MMINTR_ICALL
static int MMINTR_ICALL(
   MMINTR mmintr, MemBlob* memBlob,
   const enum ATLAS_TRANS TA, const enum ATLAS_TRANS TB,
   const int M, const int N, const int K, const SCALAR alpha,
   const TYPE *A, const int lda, const TYPE *B, const int ldb,
   const SCALAR beta, TYPE *C, const int ldc)
{
#ifdef DIRECTHSA
   return Mjoin3(PATL,icall_MMINTR,HSADECL)(
      mmintr, memBlob, TA, TB, M, N, K, alpha, A, lda, B ,ldb, beta, C, ldc);
#else
   return mmintr(
      memBlob, TA, TB, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
#endif
}
#endif

#ifdef GEADD_ICALL
static void GEADD_ICALL(
   GEADD geadd, const int M, const int N, const SCALAR alpha, const TYPE *A,
   const int lda, const SCALAR beta, TYPE *C, const int ldc)
{
#ifdef DIRECTHSA
   Mjoin3(PATL,icall_GEADD,HSADECL)(geadd, M, N, alpha, A, lda, beta, C, ldc);
#else
   geadd(M, N, alpha, A, lda, beta, C, ldc);
#endif
}
#endif

