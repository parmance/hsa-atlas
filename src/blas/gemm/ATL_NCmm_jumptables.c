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
#define ATL_NoICalls
#include "atlas_misc.h"
#include Mstr(Mjoin(Mjoin(atlas_,PRE),NCmm.h))
#include "atlas_NCmm.h"

#ifdef DIRECTHSA

HSA_FUNCTION
void Mjoin3(PATL,icall_NCMM,PHSA)(
   NCMM NCmm, const int M, const int N, const int K, const SCALAR alpha,
   const TYPE *A, const int lda, const TYPE *B, const int ldb,
   const SCALAR beta, TYPE *C, const int ldc)
{
#undef DO_CALL
#define DO_CALL(fn)                                         \
   fn(M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);        \
   break

   switch (NCmm)
   {
   default:
   case NCMM_null:
      ATL_assert (0);
      return;
   case Mjoin(Mjoin5(NCmm00,0x0x0,NN,0x0x0_aX_bX,PHSA),_fnid):
      DO_CALL(Mjoin5(NCmm00,0x0x0,NN,0x0x0_aX_bX,PHSA));
   case Mjoin(Mjoin5(NCmm00,0x0x0,NT,0x0x0_aX_bX,PHSA),_fnid):
      DO_CALL(Mjoin5(NCmm00,0x0x0,NT,0x0x0_aX_bX,PHSA));
   case Mjoin(Mjoin5(NCmm00,0x0x0,TN,0x0x0_aX_bX,PHSA),_fnid):
      DO_CALL(Mjoin5(NCmm00,0x0x0,TN,0x0x0_aX_bX,PHSA));
   case Mjoin(Mjoin5(NCmm00,0x0x0,TT,0x0x0_aX_bX,PHSA),_fnid):
      DO_CALL(Mjoin5(NCmm00,0x0x0,TT,0x0x0_aX_bX,PHSA));
   case Mjoin(Mjoin5(NCmm00,Mjoin(0x0x,KB),NN,0x0x0_aX_bX,PHSA),_fnid):
      DO_CALL(Mjoin5(NCmm00,Mjoin(0x0x,KB),NN,0x0x0_aX_bX,PHSA));
   case Mjoin(Mjoin5(NCmm00,Mjoin(0x0x,KB),NT,0x0x0_aX_bX,PHSA),_fnid):
      DO_CALL(Mjoin5(NCmm00,Mjoin(0x0x,KB),NT,0x0x0_aX_bX,PHSA));
   case Mjoin(Mjoin5(NCmm00,Mjoin(0x0x,KB),TN,0x0x0_aX_bX,PHSA),_fnid):
      DO_CALL(Mjoin5(NCmm00,Mjoin(0x0x,KB),TN,0x0x0_aX_bX,PHSA));
   case Mjoin(Mjoin5(NCmm00,Mjoin(0x0x,KB),TT,0x0x0_aX_bX,PHSA),_fnid):
      DO_CALL(Mjoin5(NCmm00,Mjoin(0x0x,KB),TT,0x0x0_aX_bX,PHSA));
   case Mjoin(Mjoin5(NCmm0,NN,0x0x0,_a1_b0,PHSA),_fnid):
      DO_CALL(Mjoin5(NCmm0,NN,0x0x0,_a1_b0,PHSA));
   case Mjoin(Mjoin5(NCmm0,NN,0x0x0,_a1_b1,PHSA),_fnid):
      DO_CALL(Mjoin5(NCmm0,NN,0x0x0,_a1_b1,PHSA));
   case Mjoin(Mjoin5(NCmm0,NN,0x0x0,_a1_bX,PHSA),_fnid):
      DO_CALL(Mjoin5(NCmm0,NN,0x0x0,_a1_bX,PHSA));
   case Mjoin(Mjoin5(NCmm0,NN,0x0x0,_aX_b0,PHSA),_fnid):
      DO_CALL(Mjoin5(NCmm0,NN,0x0x0,_aX_b0,PHSA));
   case Mjoin(Mjoin5(NCmm0,NN,0x0x0,_aX_bX,PHSA),_fnid):
      DO_CALL(Mjoin5(NCmm0,NN,0x0x0,_aX_bX,PHSA));
   case Mjoin(Mjoin5(NCmm0,NT,0x0x0,_a1_b0,PHSA),_fnid):
      DO_CALL(Mjoin5(NCmm0,NT,0x0x0,_a1_b0,PHSA));
   case Mjoin(Mjoin5(NCmm0,NT,0x0x0,_a1_b1,PHSA),_fnid):
      DO_CALL(Mjoin5(NCmm0,NT,0x0x0,_a1_b1,PHSA));
   case Mjoin(Mjoin5(NCmm0,NT,0x0x0,_a1_bX,PHSA),_fnid):
      DO_CALL(Mjoin5(NCmm0,NT,0x0x0,_a1_bX,PHSA));
   case Mjoin(Mjoin5(NCmm0,NT,0x0x0,_aX_b0,PHSA),_fnid):
      DO_CALL(Mjoin5(NCmm0,NT,0x0x0,_aX_b0,PHSA));
   case Mjoin(Mjoin5(NCmm0,NT,0x0x0,_aX_bX,PHSA),_fnid):
      DO_CALL(Mjoin5(NCmm0,NT,0x0x0,_aX_bX,PHSA));
   case Mjoin(Mjoin5(NCmm0,TN,0x0x0,_a1_b0,PHSA),_fnid):
      DO_CALL(Mjoin5(NCmm0,TN,0x0x0,_a1_b0,PHSA));
   case Mjoin(Mjoin5(NCmm0,TN,0x0x0,_a1_b1,PHSA),_fnid):
      DO_CALL(Mjoin5(NCmm0,TN,0x0x0,_a1_b1,PHSA));
   case Mjoin(Mjoin5(NCmm0,TN,0x0x0,_a1_bX,PHSA),_fnid):
      DO_CALL(Mjoin5(NCmm0,TN,0x0x0,_a1_bX,PHSA));
   case Mjoin(Mjoin5(NCmm0,TN,0x0x0,_aX_b0,PHSA),_fnid):
      DO_CALL(Mjoin5(NCmm0,TN,0x0x0,_aX_b0,PHSA));
   case Mjoin(Mjoin5(NCmm0,TN,0x0x0,_aX_bX,PHSA),_fnid):
      DO_CALL(Mjoin5(NCmm0,TN,0x0x0,_aX_bX,PHSA));
   case Mjoin(Mjoin5(NCmm0,TT,0x0x0,_a1_b0,PHSA),_fnid):
      DO_CALL(Mjoin5(NCmm0,TT,0x0x0,_a1_b0,PHSA));
   case Mjoin(Mjoin5(NCmm0,TT,0x0x0,_a1_b1,PHSA),_fnid):
      DO_CALL(Mjoin5(NCmm0,TT,0x0x0,_a1_b1,PHSA));
   case Mjoin(Mjoin5(NCmm0,TT,0x0x0,_a1_bX,PHSA),_fnid):
      DO_CALL(Mjoin5(NCmm0,TT,0x0x0,_a1_bX,PHSA));
   case Mjoin(Mjoin5(NCmm0,TT,0x0x0,_aX_b0,PHSA),_fnid):
      DO_CALL(Mjoin5(NCmm0,TT,0x0x0,_aX_b0,PHSA));
   case Mjoin(Mjoin5(NCmm0,TT,0x0x0,_aX_bX,PHSA),_fnid):
      DO_CALL(Mjoin5(NCmm0,TT,0x0x0,_aX_bX,PHSA));
   }
}

#undef DO_CALL
#endif /* DIRECTHSA */
