/* TODO copyright */

#define ATL_no_icalls
#include "atlas_misc.h"
#include Mstr(Mjoin(Mjoin(atlas_,PRE),NCmm.h))
#include "atlas_NCmm.h"

#ifdef DIRECTHSA

#ifndef HSADECLS
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

HSA_FUNCTION
void Mjoin3(PATL,icall_NCMM,HSADECL)(
   NCMM NCmm, const int M, const int N, const int K, const SCALAR alpha,
   const TYPE *A, const int lda, const TYPE *B, const int ldb,
   const SCALAR beta, TYPE *C, const int ldc)
{
#undef DO_CALL
#define DO_CALL(fn)                                         \
   fn(M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);        \
   break

   DBGLOG_CALLED_ONCE;

   switch (NCmm)
   {
   default:
      DBGLOG_CALLED_UNKNOWN;
   case NCMM_null:
      ATL_assert (0);
      return;
   case Mjoin(Mjoin5(NCmm00,0x0x0,NN,0x0x0_aX_bX,HSADECL),_fnid):
      DO_CALL(Mjoin5(NCmm00,0x0x0,NN,0x0x0_aX_bX,HSADECL));
   case Mjoin(Mjoin5(NCmm00,0x0x0,NT,0x0x0_aX_bX,HSADECL),_fnid):
      DO_CALL(Mjoin5(NCmm00,0x0x0,NT,0x0x0_aX_bX,HSADECL));
   case Mjoin(Mjoin5(NCmm00,0x0x0,TN,0x0x0_aX_bX,HSADECL),_fnid):
      DO_CALL(Mjoin5(NCmm00,0x0x0,TN,0x0x0_aX_bX,HSADECL));
   case Mjoin(Mjoin5(NCmm00,0x0x0,TT,0x0x0_aX_bX,HSADECL),_fnid):
      DO_CALL(Mjoin5(NCmm00,0x0x0,TT,0x0x0_aX_bX,HSADECL));
   case Mjoin(Mjoin5(NCmm00,Mjoin(0x0x,KB),NN,0x0x0_aX_bX,HSADECL),_fnid):
      DO_CALL(Mjoin5(NCmm00,Mjoin(0x0x,KB),NN,0x0x0_aX_bX,HSADECL));
   case Mjoin(Mjoin5(NCmm00,Mjoin(0x0x,KB),NT,0x0x0_aX_bX,HSADECL),_fnid):
      DO_CALL(Mjoin5(NCmm00,Mjoin(0x0x,KB),NT,0x0x0_aX_bX,HSADECL));
   case Mjoin(Mjoin5(NCmm00,Mjoin(0x0x,KB),TN,0x0x0_aX_bX,HSADECL),_fnid):
      DO_CALL(Mjoin5(NCmm00,Mjoin(0x0x,KB),TN,0x0x0_aX_bX,HSADECL));
   case Mjoin(Mjoin5(NCmm00,Mjoin(0x0x,KB),TT,0x0x0_aX_bX,HSADECL),_fnid):
      DO_CALL(Mjoin5(NCmm00,Mjoin(0x0x,KB),TT,0x0x0_aX_bX,HSADECL));
   case Mjoin(Mjoin5(NCmm0,NN,0x0x0,_a1_b0,HSADECL),_fnid):
      DO_CALL(Mjoin5(NCmm0,NN,0x0x0,_a1_b0,HSADECL));
   case Mjoin(Mjoin5(NCmm0,NN,0x0x0,_a1_b1,HSADECL),_fnid):
      DO_CALL(Mjoin5(NCmm0,NN,0x0x0,_a1_b1,HSADECL));
   case Mjoin(Mjoin5(NCmm0,NN,0x0x0,_a1_bX,HSADECL),_fnid):
      DO_CALL(Mjoin5(NCmm0,NN,0x0x0,_a1_bX,HSADECL));
   case Mjoin(Mjoin5(NCmm0,NN,0x0x0,_aX_b0,HSADECL),_fnid):
      DO_CALL(Mjoin5(NCmm0,NN,0x0x0,_aX_b0,HSADECL));
   case Mjoin(Mjoin5(NCmm0,NN,0x0x0,_aX_bX,HSADECL),_fnid):
      DO_CALL(Mjoin5(NCmm0,NN,0x0x0,_aX_bX,HSADECL));
   case Mjoin(Mjoin5(NCmm0,NT,0x0x0,_a1_b0,HSADECL),_fnid):
      DO_CALL(Mjoin5(NCmm0,NT,0x0x0,_a1_b0,HSADECL));
   case Mjoin(Mjoin5(NCmm0,NT,0x0x0,_a1_b1,HSADECL),_fnid):
      DO_CALL(Mjoin5(NCmm0,NT,0x0x0,_a1_b1,HSADECL));
   case Mjoin(Mjoin5(NCmm0,NT,0x0x0,_a1_bX,HSADECL),_fnid):
      DO_CALL(Mjoin5(NCmm0,NT,0x0x0,_a1_bX,HSADECL));
   case Mjoin(Mjoin5(NCmm0,NT,0x0x0,_aX_b0,HSADECL),_fnid):
      DO_CALL(Mjoin5(NCmm0,NT,0x0x0,_aX_b0,HSADECL));
   case Mjoin(Mjoin5(NCmm0,NT,0x0x0,_aX_bX,HSADECL),_fnid):
      DO_CALL(Mjoin5(NCmm0,NT,0x0x0,_aX_bX,HSADECL));
   case Mjoin(Mjoin5(NCmm0,TN,0x0x0,_a1_b0,HSADECL),_fnid):
      DO_CALL(Mjoin5(NCmm0,TN,0x0x0,_a1_b0,HSADECL));
   case Mjoin(Mjoin5(NCmm0,TN,0x0x0,_a1_b1,HSADECL),_fnid):
      DO_CALL(Mjoin5(NCmm0,TN,0x0x0,_a1_b1,HSADECL));
   case Mjoin(Mjoin5(NCmm0,TN,0x0x0,_a1_bX,HSADECL),_fnid):
      DO_CALL(Mjoin5(NCmm0,TN,0x0x0,_a1_bX,HSADECL));
   case Mjoin(Mjoin5(NCmm0,TN,0x0x0,_aX_b0,HSADECL),_fnid):
      DO_CALL(Mjoin5(NCmm0,TN,0x0x0,_aX_b0,HSADECL));
   case Mjoin(Mjoin5(NCmm0,TN,0x0x0,_aX_bX,HSADECL),_fnid):
      DO_CALL(Mjoin5(NCmm0,TN,0x0x0,_aX_bX,HSADECL));
   case Mjoin(Mjoin5(NCmm0,TT,0x0x0,_a1_b0,HSADECL),_fnid):
      DO_CALL(Mjoin5(NCmm0,TT,0x0x0,_a1_b0,HSADECL));
   case Mjoin(Mjoin5(NCmm0,TT,0x0x0,_a1_b1,HSADECL),_fnid):
      DO_CALL(Mjoin5(NCmm0,TT,0x0x0,_a1_b1,HSADECL));
   case Mjoin(Mjoin5(NCmm0,TT,0x0x0,_a1_bX,HSADECL),_fnid):
      DO_CALL(Mjoin5(NCmm0,TT,0x0x0,_a1_bX,HSADECL));
   case Mjoin(Mjoin5(NCmm0,TT,0x0x0,_aX_b0,HSADECL),_fnid):
      DO_CALL(Mjoin5(NCmm0,TT,0x0x0,_aX_b0,HSADECL));
   case Mjoin(Mjoin5(NCmm0,TT,0x0x0,_aX_bX,HSADECL),_fnid):
      DO_CALL(Mjoin5(NCmm0,TT,0x0x0,_aX_bX,HSADECL));
   }
}

#endif /* DIRECTHSA */
