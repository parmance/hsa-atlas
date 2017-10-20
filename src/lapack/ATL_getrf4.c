#include "atlas_misc.h"
#include "atlas_lapack.h"
#include "atlas_level2.h"
#include "atlas_lamch.h"

#if defined(ATL_AVX) && defined(DREAL)
   #include <immintrin.h>
static int irk1amax(const int M, const TYPE s0, TYPE *A0, const int lda)
/*
 * This routine is used after 1 column of LU has been pivoted, and it
 * merges three steps into one:
 *    A0[1:M] = s0 * A0[1:M]
 *    A1[1:M] = A1[0] * A1[1:M]
 * RETURNS: iamax of A1 *after* update
 */
{
   const int im = Mmin(4,M), m4=(M>>2)<<2, mr = M-im-m4;
   TYPE *A1 = A0 + lda;
   register int i, imax=1;
   const register TYPE s1 = A1[0];
   register TYPE amax = ATL_rzero;
//   printf("M=%d, im=%d, m4=%d, mr=%d\n", M, im, m4, mr);
/*
 * This initial peel will keep vectors aligned if original matrix is
 */
   for (i=1; i < im; i++)
   {
      register TYPE t0 = A0[i], t1 = A1[i];
      t0 *= s0;
      A0[i] = t0;
      t1 -= t0 * s1;
      A1[i] = t1;
      t1 = Mabs(t1);
      if (t1 > amax) { amax = t1; imax = i; }
   }
//         printf("A1[0:3] = %e, %e, %e, %e\n", A1[0], A1[1], A1[2], A1[3]);
   if (m4)
   {
      register __m256d vindx = {4.0, 5.0, 6.0, 7.0};
      register __m256d vimax={(double)imax,(double)imax,(double)imax,
                              (double)imax};
      const register __m256d viinc = {4.0,4.0,4.0,4.0};
      const register __m256d v0 = {s0,s0,s0,s0}, vabs = {-0.0,-0.0,-0.0,-0.0};
      const register __m256d v1 = {s1,s1,s1,s1};
      register __m256d vamax = {amax,amax,amax,amax};
      TYPE VV[4], VI[4];

      for (; i != m4; i += 4)
      {
         register __m256d t0, t1, t2;
         t0 = _mm256_loadu_pd(A0+i);
         t0 = _mm256_mul_pd(t0, v0);
         t1 = _mm256_loadu_pd(A1+i);
         t2 = _mm256_mul_pd(t0, v1);
         _mm256_storeu_pd(A0+i, t0);
         t1 = _mm256_sub_pd(t1, t2);
         _mm256_storeu_pd(A1+i, t1);
//         printf("A1[4:7] = %e, %e, %e, %e\n", A1[4], A1[5], A1[6], A1[7]);
         t1 = _mm256_andnot_pd(vabs,t1);   /* t1 = ABS(t1) */
         t2 = _mm256_cmp_pd(t1, vamax, 14);  /* t2 =  (t1 > vamax) */
         vimax = _mm256_blendv_pd(vimax, vindx, t2);
         vindx = _mm256_add_pd(vindx, viinc);
         vamax = _mm256_blendv_pd(vamax, t1, t2);
      }
/*
 *    Reduce vector vamax/vimax to scalar
 */
      {
         __m128d t0, t1, t2, t3, t4;
         TYPE vv;

         t0 = _mm256_extractf128_pd(vamax, 0);
         t2 = _mm256_extractf128_pd(vamax, 1);
         t4 = _mm_cmp_pd(t2, t0, 14);           /* t4 = (t2 > t0) */
         t0 = _mm_blendv_pd(t0, t2, t4);
         t1 = _mm256_extractf128_pd(vimax, 0);
         t3 = _mm256_extractf128_pd(vimax, 1);
         t1 = _mm_blendv_pd(t1, t3, t4);
         _mm_store_sd(&vv, t0);
         amax = vv;
         _mm_storeh_pd(&vv, t0);
         if (amax < vv)
            t1 = _mm_unpackhi_pd(t1, t1);
         imax = _mm_cvtsd_si32(t1);
      }
   }
/*
 * Finish off any remaining elements with scalar code
 */
   for (; i < M; i++)
   {
      register TYPE t0 = A0[i], t1 = A1[i];
      t0 *= s0;
      A0[i] = t0;
      t1 -= t0 * s1;
      A1[i] = t1;
      t1 = Mabs(t1);
      if (t1 > amax) { amax = t1; imax = i; }
   }
   return(imax);
}
#elif defined(ATL_AVX) && defined(SREAL)
   #include <immintrin.h>
static int irk1amax(const int M, const TYPE s0, TYPE *A0, const int lda)
/*
 * This routine is used after 1 column of LU has been pivoted, and it
 * merges three steps into one:
 *    A0[1:M] = s0 * A0[1:M]
 *    A1[1:M] = A1[0] * A1[1:M]
 * RETURNS: iamax of A1 *after* update
 */
{
   const int im = Mmin(8,M), m8=(M>>3)<<3;
   TYPE *A1 = A0 + lda;
   register int i, imax=1;
   const register TYPE s1 = A1[0];
   register TYPE amax = ATL_rzero;
//   printf("M=%d, im=%d, m4=%d, mr=%d\n", M, im, m4, mr);
/*
 * This initial peel will keep vectors aligned if original matrix is
 */
   for (i=1; i < im; i++)
   {
      register TYPE t0 = A0[i], t1 = A1[i];
      t0 *= s0;
      A0[i] = t0;
      t1 -= t0 * s1;
      A1[i] = t1;
      t1 = Mabs(t1);
      if (t1 > amax) { amax = t1; imax = i; }
   }
//         printf("A1[0:3] = %e, %e, %e, %e\n", A1[0], A1[1], A1[2], A1[3]);
//   printf("imax=%d, amax=%e\n", imax, amax);
   if (m8)
   {
      register __m256 vindx={8.0f,9.0f,10.0f,11.0f,12.0f,13.0f,14.0f,15.0f};
      register __m256 vimax={(float)imax,(float)imax,(float)imax, (float)imax,
                             (float)imax,(float)imax,(float)imax, (float)imax};
      const register __m256 viinc = {8.0f,8.0f,8.0f,8.0f,8.0f,8.0f,8.0f,8.0f};
      const register __m256 v0 = {s0,s0,s0,s0,s0,s0,s0,s0};
      const register __m256 vabs = {-0.0f,-0.0f,-0.0f,-0.0f,
                                    -0.0f,-0.0f,-0.0f,-0.0f};
      const register __m256 v1 = {s1,s1,s1,s1,s1,s1,s1,s1};
      register __m256 vamax = {amax,amax,amax,amax,amax,amax,amax,amax};
      TYPE VV[8], VI[8];

      for (; i != m8; i += 8)
      {
         register __m256 t0, t1, t2;
         t0 = _mm256_loadu_ps(A0+i);
         t0 = _mm256_mul_ps(t0, v0);
         t1 = _mm256_loadu_ps(A1+i);
         t2 = _mm256_mul_ps(t0, v1);
         _mm256_storeu_ps(A0+i, t0);
         t1 = _mm256_sub_ps(t1, t2);
         _mm256_storeu_ps(A1+i, t1);
//         printf("A1[4:7] = %e, %e, %e, %e\n", A1[4], A1[5], A1[6], A1[7]);
         t1 = _mm256_andnot_ps(vabs,t1);     /* t1 = ABS(t1) */
         t2 = _mm256_cmp_ps(t1, vamax, 14);  /* t2 =  (t1 > vamax) */
         vimax = _mm256_blendv_ps(vimax, vindx, t2);
         vindx = _mm256_add_ps(vindx, viinc);
         vamax = _mm256_blendv_ps(vamax, t1, t2);
      }
/*
 *    Reduce vector values vamax/vimax to scalar values amax/imax
 */
      {
         register __m128 t0, t1, t2, t3, t4;
         TYPE vv;

         t0 = _mm256_extractf128_ps(vamax, 0);
         t1 = _mm256_extractf128_ps(vamax, 1);
         t2 = _mm_cmp_ps(t1, t0, 14);  /* t2 =  (t1 > vamax) */
         t0 = _mm_blendv_ps(t0, t1, t2);
         t1 = _mm256_extractf128_ps(vimax, 0);
         t3 = _mm256_extractf128_ps(vimax, 1);
         t1 = _mm_blendv_ps(t1, t3, t2);

         t2 = _mm_movehl_ps(t0, t0);
         t4 = _mm_cmp_ps(t2, t0, 14);  /* t4 =  (t1 > t0) */
         t3 = _mm_movehl_ps(t1, t1);
         t0 = _mm_blendv_ps(t0, t2, t4);
         t1 = _mm_blendv_ps(t1, t3, t4);

         _mm_store_ss(&vv, t0);
         amax = vv;
         t0 = _mm_shuffle_ps(t0, t0, 1);
         _mm_store_ss(&vv, t0);
         if (vv > amax)
         {
            amax = vv;
            t1 = _mm_shuffle_ps(t1, t1, 1);
         }
         imax = _mm_cvtss_si32(t1);
      }
   }
/*
 * Finish off any remaining elements with scalar code
 */
   for (; i < M; i++)
   {
      register TYPE t0 = A0[i], t1 = A1[i];
      t0 *= s0;
      A0[i] = t0;
      t1 -= t0 * s1;
      A1[i] = t1;
      t1 = Mabs(t1);
      if (t1 > amax) { amax = t1; imax = i; }
   }
   return(imax);
}
#else
static int irk1amax(const int M, const TYPE s0, TYPE *A0, const int lda)
/*
 * This routine is used after 1 column of LU has been pivoted, and it
 * merges three steps into one:
 *    A0[1:M] = s0 * A0[1:M]
 *    A1[1:M] = A1[0] * A1[1:M]
 * RETURNS: iamax of A1 *after* update
 */
{
   TYPE *A1 = A0 + lda;
   register int i, imax=1;
   const register TYPE s1 = A1[0];
   register TYPE amax = ATL_rzero;
   for (i=1; i != M; i++)
   {
      register TYPE t0 = A0[i], t1 = A1[i];
      t0 *= s0;
      A0[i] = t0;
      t1 -= t0 * s1;
      A1[i] = t1;
      t1 = Mabs(t1);
      if (t1 > amax) { amax = t1; imax = i; }
   }
   return(imax);
}
#endif

#if defined(ATL_AVX) && defined(DREAL)
static int irk2amax(const int M, const TYPE s0, TYPE *A0, const int lda)
/*
 * This routine is used after 2 columns of LU have been pivoted, and it
 * merges several updates into one pass through memory:
 *   A2[1:M] -= A2[0]*A0[1:M] // GER frm 1st step of LU
 *   A1[2:M] = s0 * A1[2:M]   // scale A1 by previously discovered pivot
 *   A2[2:M] -= A2[1]*A1[2:M] // GER frm 2nd step of LU
 * RETURNS: iamax of A2[2:M] AFTER updates applied
 */
{
   TYPE *A1 = A0 + lda, *A2 = A0 + lda+lda;
   const register TYPE s1 = -A2[0];
   register TYPE s2 = A2[1];
   register TYPE amax = ATL_rzero;
   register int i, imax=2;
   const int im=Mmin(4,M), m4=(M>>2)<<2;
   s2 += A0[1] * s1;
   A2[1] = s2;
   s2 = - s2;
   for (i=2; i != im; i++)
   {
      register TYPE t0 = A0[i], t1 = A1[i], t2 = A2[i];
      t1 *= s0;
      A1[i] = t1;
      t2 += t0*s1;
      t2 += t1*s2;
      A2[i] = t2;
      t2 = Mabs(t2);
      if (t2 > amax) { amax = t2; imax = i; }
   }
   if (m4)
   {
      register __m256d vindx = {4.0, 5.0, 6.0, 7.0};
      register __m256d vimax={(double)imax,(double)imax,(double)imax,
                              (double)imax};
      const register __m256d viinc = {4.0,4.0,4.0,4.0};
      const register __m256d v0 = {s0,s0,s0,s0}, vabs = {-0.0,-0.0,-0.0,-0.0};
      const register __m256d v1 = {s1,s1,s1,s1}, v2 = {s2,s2,s2,s2};
      register __m256d vamax = {amax,amax,amax,amax};
      TYPE VV[4], VI[4];

      for (; i != m4; i += 4)
      {
         register __m256d t0, t1, t2, t3;
         t1 = _mm256_loadu_pd(A1+i);
         t1 = _mm256_mul_pd(t1, v0);
         t0 = _mm256_loadu_pd(A0+i);
         t3 = _mm256_mul_pd(t0, v1);
         t2 = _mm256_loadu_pd(A2+i);
         t2 = _mm256_add_pd(t2, t3);
         _mm256_storeu_pd(A1+i, t1);
         t3 = _mm256_mul_pd(t1, v2);
         t2 = _mm256_add_pd(t2, t3);
         _mm256_storeu_pd(A2+i, t2);
         t2 = _mm256_andnot_pd(vabs,t2);        /* t2 = ABS(t2) */
         t3 = _mm256_cmp_pd(t2, vamax, 14);     /* t3 =  (t2 > vamax) */
         vimax = _mm256_blendv_pd(vimax, vindx, t3);
         vindx = _mm256_add_pd(vindx, viinc);
         vamax = _mm256_blendv_pd(vamax, t2, t3);
      }
/*
 *    Reduce vector vamax/vimax to scalar
 */
      {
         __m128d t0, t1, t2, t3, t4;
         TYPE vv;

         t0 = _mm256_extractf128_pd(vamax, 0);
         t2 = _mm256_extractf128_pd(vamax, 1);
         t4 = _mm_cmp_pd(t2, t0, 14);           /* t4 = (t2 > t0) */
         t0 = _mm_blendv_pd(t0, t2, t4);
         t1 = _mm256_extractf128_pd(vimax, 0);
         t3 = _mm256_extractf128_pd(vimax, 1);
         t1 = _mm_blendv_pd(t1, t3, t4);
         _mm_store_sd(&vv, t0);
         amax = vv;
         _mm_storeh_pd(&vv, t0);
         if (amax < vv)
            t1 = _mm_unpackhi_pd(t1, t1);
         imax = _mm_cvtsd_si32(t1);
      }
   }
   for (; i != M; i++)
   {
      register TYPE t0 = A0[i], t1 = A1[i], t2 = A2[i];
      t1 *= s0;
      A1[i] = t1;
      t2 += t0*s1;
      t2 += t1*s2;
      A2[i] = t2;
      t2 = Mabs(t2);
      if (t2 > amax) { amax = t2; imax = i; }
   }
   return(imax);
}
#elif defined(ATL_AVX) && defined(SREAL)
static int irk2amax(const int M, const TYPE s0, TYPE *A0, const int lda)
/*
 * This routine is used after 2 columns of LU have been pivoted, and it
 * merges several updates into one pass through memory:
 *   A2[1:M] -= A2[0]*A0[1:M] // GER frm 1st step of LU
 *   A1[2:M] = s0 * A1[2:M]   // scale A1 by previously discovered pivot
 *   A2[2:M] -= A2[1]*A1[2:M] // GER frm 2nd step of LU
 * RETURNS: iamax of A2[2:M] AFTER updates applied
 */
{
   TYPE *A1 = A0 + lda, *A2 = A0 + lda+lda;
   const register TYPE s1 = -A2[0];
   register TYPE s2 = A2[1];
   register TYPE amax = ATL_rzero;
   register int i, imax=2;
   const int im=Mmin(8,M), m8=(M>>3)<<3;
   s2 += A0[1] * s1;
   A2[1] = s2;
   s2 = - s2;
   for (i=2; i != im; i++)
   {
      register TYPE t0 = A0[i], t1 = A1[i], t2 = A2[i];
      t1 *= s0;
      A1[i] = t1;
      t2 += t0*s1;
      t2 += t1*s2;
      A2[i] = t2;
      t2 = Mabs(t2);
      if (t2 > amax) { amax = t2; imax = i; }
   }
   if (m8)
   {
      register __m256 vindx = {8.0f,9.0f,10.0f,11.0f,12.0f,13.0f,14.0f,15.0f};
      register __m256 vimax={(float)imax,(float)imax,(float)imax,(float)imax,
                             (float)imax,(float)imax,(float)imax,(float)imax};
      const register __m256 viinc = {8.0,8.0,8.0,8.0,8.0,8.0,8.0,8.0};
      const register __m256 vabs = {-0.0,-0.0,-0.0,-0.0};
      const register __m256 v0 = {s0,s0,s0,s0,s0,s0,s0,s0};
      const register __m256 v1 = {s1,s1,s1,s1,s1,s1,s1,s1};
      const register __m256 v2 = {s2,s2,s2,s2,s2,s2,s2,s2};
      register __m256 vamax = {amax,amax,amax,amax,amax,amax,amax,amax};
      TYPE VV[8], VI[8];

      for (; i != m8; i += 8)
      {
         register __m256 t0, t1, t2, t3;
         t1 = _mm256_loadu_ps(A1+i);
         t1 = _mm256_mul_ps(t1, v0);
         t0 = _mm256_loadu_ps(A0+i);
         t3 = _mm256_mul_ps(t0, v1);
         t2 = _mm256_loadu_ps(A2+i);
         t2 = _mm256_add_ps(t2, t3);
         _mm256_storeu_ps(A1+i, t1);
         t3 = _mm256_mul_ps(t1, v2);
         t2 = _mm256_add_ps(t2, t3);
         _mm256_storeu_ps(A2+i, t2);
         t2 = _mm256_andnot_ps(vabs,t2);        /* t2 = ABS(t2) */
         t3 = _mm256_cmp_ps(t2, vamax, 14);     /* t3 =  (t2 > vamax) */
         vimax = _mm256_blendv_ps(vimax, vindx, t3);
         vindx = _mm256_add_ps(vindx, viinc);
         vamax = _mm256_blendv_ps(vamax, t2, t3);
      }
/*
 *    Reduce vector values vamax/vimax to scalar values amax/imax
 */
      {
         register __m128 t0, t1, t2, t3, t4;
         TYPE vv;

         t0 = _mm256_extractf128_ps(vamax, 0);
         t1 = _mm256_extractf128_ps(vamax, 1);
         t2 = _mm_cmp_ps(t1, t0, 14);  /* t2 =  (t1 > vamax) */
         t0 = _mm_blendv_ps(t0, t1, t2);
         t1 = _mm256_extractf128_ps(vimax, 0);
         t3 = _mm256_extractf128_ps(vimax, 1);
         t1 = _mm_blendv_ps(t1, t3, t2);

         t2 = _mm_movehl_ps(t0, t0);
         t4 = _mm_cmp_ps(t2, t0, 14);  /* t4 =  (t1 > t0) */
         t3 = _mm_movehl_ps(t1, t1);
         t0 = _mm_blendv_ps(t0, t2, t4);
         t1 = _mm_blendv_ps(t1, t3, t4);

         _mm_store_ss(&vv, t0);
         amax = vv;
         t0 = _mm_shuffle_ps(t0, t0, 1);
         _mm_store_ss(&vv, t0);
         if (vv > amax)
         {
            amax = vv;
            t1 = _mm_shuffle_ps(t1, t1, 1);
         }
         imax = _mm_cvtss_si32(t1);
      }
   }
   for (; i != M; i++)
   {
      register TYPE t0 = A0[i], t1 = A1[i], t2 = A2[i];
      t1 *= s0;
      A1[i] = t1;
      t2 += t0*s1;
      t2 += t1*s2;
      A2[i] = t2;
      t2 = Mabs(t2);
      if (t2 > amax) { amax = t2; imax = i; }
   }
   return(imax);
}
#else
static int irk2amax(const int M, const TYPE s0, TYPE *A0, const int lda)
/*
 * This routine is used after 2 columns of LU have been pivoted, and it
 * merges several updates into one pass through memory:
 *   A2[1:M] -= A2[0]*A0[1:M] // GER frm 1st step of LU
 *   A1[2:M] = s0 * A1[2:M]   // scale A1 by previously discovered pivot
 *   A2[2:M] -= A2[1]*A1[2:M] // GER frm 2nd step of LU
 * RETURNS: iamax of A2[2:M] AFTER updates applied
 */
{
   TYPE *A1 = A0 + lda, *A2 = A0 + lda+lda;
   const register TYPE s1 = -A2[0];
   register TYPE s2 = A2[1];
   register TYPE amax = ATL_rzero;
   register int i, imax=2;
   s2 += A0[1] * s1;
   A2[1] = s2;
   s2 = - s2;
   for (i=2; i != M; i++)
   {
      register TYPE t0 = A0[i], t1 = A1[i], t2 = A2[i];
      t1 *= s0;
      A1[i] = t1;
      t2 += t0*s1;
      t2 += t1*s2;
      A2[i] = t2;
      t2 = Mabs(t2);
      if (t2 > amax) { amax = t2; imax = i; }
   }
   return(imax);
}
#endif

#if defined(ATL_AVX) && defined(DREAL)
static int irk3amax(const int M, const TYPE s0, TYPE *A0, const int lda)
/*
 * This routine is used after 3 columns of LU have been pivoted, and it
 * merges several updates into one pass through memory:
 *   A3[1:M] -= A3[0]*A0[1:M] // GER frm 1st step of LU
 *   A3[2:M] -= A3[1]*A1[2:M] // GER frm 2nd step of LU
 *   A2[3:M] = s0 * A2[3:M]   // scale A2 by previously discovered pivot
 *   A3[3:M] -= A3[2]*A2[3:M] // GER frm 3rd step of LU
 * RETURNS: iamax of A3[2:M] AFTER updates applied
 */
{
   TYPE *A1 = A0 + lda, *A2 = A0 + lda+lda, *A3 = A2+lda;
   const register TYPE s1 = -A3[0];
   register TYPE s2 = A3[1], s3 = A3[2];
   register TYPE amax = ATL_rzero;
   register int i, imax=3;
   const int im = Mmin(4,M), m4=(M>>2)<<2;

   s2 += s1 * A0[1];
   A3[1] = s2;
   s2 = -s2;
   s3 += A0[2]*s1;
   s3 += A1[2]*s2;
   A3[2] = s3;
   s3 = -s3;
   for (i=3; i != im; i++)
   {
      register TYPE t0 = A0[i], t1 = A1[i], t2 = A2[i], t3 = A3[i];
      t2 *= s0;
      A2[i] = t2;
      t3 += t0*s1;
      t3 += t1*s2;
      t3 += t2*s3;
      A3[i] = t3;
      t3 = Mabs(t3);
      if (t3 > amax) { amax = t3; imax = i; }
   }
   if (m4)
   {
      register __m256d vindx = {4.0, 5.0, 6.0, 7.0};
      register __m256d vimax={(double)imax,(double)imax,(double)imax,
                              (double)imax};
      const register __m256d viinc = {4.0,4.0,4.0,4.0};
      const register __m256d v0 = {s0,s0,s0,s0}, vabs = {-0.0,-0.0,-0.0,-0.0};
      const register __m256d v1 = {s1,s1,s1,s1}, v2 = {s2,s2,s2,s2};
      const register __m256d v3 = {s3,s3,s3,s3};
      register __m256d vamax = {amax,amax,amax,amax};
      TYPE VV[4], VI[4];
      for (; i != m4; i += 4)
      {
         register __m256d t0, t1, t2, t3, t4;
         t0 = _mm256_loadu_pd(A0+i);
         t1 = _mm256_loadu_pd(A1+i);
         t2 = _mm256_loadu_pd(A2+i);
         t3 = _mm256_loadu_pd(A3+i);

         t2 = _mm256_mul_pd(t2, v0);
         _mm256_storeu_pd(A2+i, t2);
         t4 = _mm256_mul_pd(t0, v1);
         t3 = _mm256_add_pd(t3, t4);
         t4 = _mm256_mul_pd(t1, v2);
         t3 = _mm256_add_pd(t3, t4);
         t4 = _mm256_mul_pd(t2, v3);
         t3 = _mm256_add_pd(t3, t4);
         _mm256_storeu_pd(A3+i, t3);
         t3 = _mm256_andnot_pd(vabs,t3);        /* t3 = ABS(t3) */
         t4 = _mm256_cmp_pd(t3, vamax, 14);  /* t4 =  (t3 > vamax) */
         vimax = _mm256_blendv_pd(vimax, vindx, t4);
         vindx = _mm256_add_pd(vindx, viinc);
         vamax = _mm256_blendv_pd(vamax, t3, t4);
      }
/*
 *    Reduce vector vamax/vimax to scalar
 */
      {
         __m128d t0, t1, t2, t3, t4;
         TYPE vv;

         t0 = _mm256_extractf128_pd(vamax, 0);
         t2 = _mm256_extractf128_pd(vamax, 1);
         t4 = _mm_cmp_pd(t2, t0, 14);           /* t4 = (t2 > t0) */
         t0 = _mm_blendv_pd(t0, t2, t4);
         t1 = _mm256_extractf128_pd(vimax, 0);
         t3 = _mm256_extractf128_pd(vimax, 1);
         t1 = _mm_blendv_pd(t1, t3, t4);
         _mm_store_sd(&vv, t0);
         amax = vv;
         _mm_storeh_pd(&vv, t0);
         if (amax < vv)
            t1 = _mm_unpackhi_pd(t1, t1);
         imax = _mm_cvtsd_si32(t1);
      }
   }
   for (; i != M; i++)
   {
      register TYPE t0 = A0[i], t1 = A1[i], t2 = A2[i], t3 = A3[i];
      t2 *= s0;
      A2[i] = t2;
      t3 += t0*s1;
      t3 += t1*s2;
      t3 += t2*s3;
      A3[i] = t3;
      t3 = Mabs(t3);
      if (t3 > amax) { amax = t3; imax = i; }
   }
   return(imax);
}
#elif defined(ATL_AVX) && defined(SREAL)
static int irk3amax(const int M, const TYPE s0, TYPE *A0, const int lda)
/*
 * This routine is used after 3 columns of LU have been pivoted, and it
 * merges several updates into one pass through memory:
 *   A3[1:M] -= A3[0]*A0[1:M] // GER frm 1st step of LU
 *   A3[2:M] -= A3[1]*A1[2:M] // GER frm 2nd step of LU
 *   A2[3:M] = s0 * A2[3:M]   // scale A2 by previously discovered pivot
 *   A3[3:M] -= A3[2]*A2[3:M] // GER frm 3rd step of LU
 * RETURNS: iamax of A3[2:M] AFTER updates applied
 */
{
   TYPE *A1 = A0 + lda, *A2 = A0 + lda+lda, *A3 = A2+lda;
   const register TYPE s1 = -A3[0];
   register TYPE s2 = A3[1], s3 = A3[2];
   register TYPE amax = ATL_rzero;
   register int i, imax=3;
   const int im = Mmin(8,M), m8=(M>>3)<<3;

   s2 += s1 * A0[1];
   A3[1] = s2;
   s2 = -s2;
   s3 += A0[2]*s1;
   s3 += A1[2]*s2;
   A3[2] = s3;
   s3 = -s3;
/*
 * This initial peel will keep vectors aligned if original matrix is aligned
 */
   for (i=3; i != im; i++)
   {
      register TYPE t0 = A0[i], t1 = A1[i], t2 = A2[i], t3 = A3[i];
      t2 *= s0;
      A2[i] = t2;
      t3 += t0*s1;
      t3 += t1*s2;
      t3 += t2*s3;
      A3[i] = t3;
      t3 = Mabs(t3);
      if (t3 > amax) { amax = t3; imax = i; }
   }
   if (m8)
   {
      register __m256 vindx = {8.0f,9.0f,10.0f,11.0f,12.0f,13.0f,14.0f,15.0f};
      register __m256 vimax={(float)imax,(float)imax,(float)imax, (float)imax};
      const register __m256 viinc = {8.0,8.0,8.0,8.0,8.0,8.0,8.0,8.0};
      const register __m256 v0 = {s0,s0,s0,s0,s0,s0,s0,s0};
      const register __m256 vabs = {-0.0f,-0.0f,-0.0f,-0.0f,
                                    -0.0f,-0.0f,-0.0f,-0.0f};
      const register __m256 v1 = {s1,s1,s1,s1,s1,s1,s1,s1};
      const register __m256 v2 = {s2,s2,s2,s2,s2,s2,s2,s2};
      const register __m256 v3 = {s3,s3,s3,s3,s3,s3,s3,s3};
      register __m256 vamax = {amax,amax,amax,amax,amax,amax,amax,amax};
      TYPE VV[8], VI[8];
      for (; i != m8; i += 8)
      {
         register __m256 t0, t1, t2, t3, t4;
         t0 = _mm256_loadu_ps(A0+i);
         t1 = _mm256_loadu_ps(A1+i);
         t2 = _mm256_loadu_ps(A2+i);
         t3 = _mm256_loadu_ps(A3+i);

         t2 = _mm256_mul_ps(t2, v0);
         _mm256_storeu_ps(A2+i, t2);
         t4 = _mm256_mul_ps(t0, v1);
         t3 = _mm256_add_ps(t3, t4);
         t4 = _mm256_mul_ps(t1, v2);
         t3 = _mm256_add_ps(t3, t4);
         t4 = _mm256_mul_ps(t2, v3);
         t3 = _mm256_add_ps(t3, t4);
         _mm256_storeu_ps(A3+i, t3);
         t3 = _mm256_andnot_ps(vabs,t3);        /* t3 = ABS(t3) */
         t4 = _mm256_cmp_ps(t3, vamax, 14);  /* t4 =  (t3 > vamax) */
         vimax = _mm256_blendv_ps(vimax, vindx, t4);
         vindx = _mm256_add_ps(vindx, viinc);
         vamax = _mm256_blendv_ps(vamax, t3, t4);
      }
/*
 *    Reduce vector values vamax/vimax to scalar values amax/imax
 */
      {
         register __m128 t0, t1, t2, t3, t4;
         TYPE vv;

         t0 = _mm256_extractf128_ps(vamax, 0);
         t1 = _mm256_extractf128_ps(vamax, 1);
         t2 = _mm_cmp_ps(t1, t0, 14);  /* t2 =  (t1 > vamax) */
         t0 = _mm_blendv_ps(t0, t1, t2);
         t1 = _mm256_extractf128_ps(vimax, 0);
         t3 = _mm256_extractf128_ps(vimax, 1);
         t1 = _mm_blendv_ps(t1, t3, t2);

         t2 = _mm_movehl_ps(t0, t0);
         t4 = _mm_cmp_ps(t2, t0, 14);  /* t4 =  (t1 > t0) */
         t3 = _mm_movehl_ps(t1, t1);
         t0 = _mm_blendv_ps(t0, t2, t4);
         t1 = _mm_blendv_ps(t1, t3, t4);

         _mm_store_ss(&vv, t0);
         amax = vv;
         t0 = _mm_shuffle_ps(t0, t0, 1);
         _mm_store_ss(&vv, t0);
         if (vv > amax)
         {
            amax = vv;
            t1 = _mm_shuffle_ps(t1, t1, 1);
         }
         imax = _mm_cvtss_si32(t1);
      }
//      printf("VV=%e, %e, %e, %e\n", VV[0], VV[1], VV[2], VV[3]);
//      printf("VI=%d, %d, %d, %d\n",(int)VI[0],(int)VI[1],(int)VI[2],(int)VI[3]);
   }
   for (; i != M; i++)
   {
      register TYPE t0 = A0[i], t1 = A1[i], t2 = A2[i], t3 = A3[i];
      t2 *= s0;
      A2[i] = t2;
      t3 += t0*s1;
      t3 += t1*s2;
      t3 += t2*s3;
      A3[i] = t3;
      t3 = Mabs(t3);
      if (t3 > amax) { amax = t3; imax = i; }
   }
   return(imax);
}
#else
static int irk3amax(const int M, const TYPE s0, TYPE *A0, const int lda)
/*
 * This routine is used after 3 columns of LU have been pivoted, and it
 * merges several updates into one pass through memory:
 *   A3[1:M] -= A3[0]*A0[1:M] // GER frm 1st step of LU
 *   A3[2:M] -= A3[1]*A1[2:M] // GER frm 2nd step of LU
 *   A2[3:M] = s0 * A2[3:M]   // scale A2 by previously discovered pivot
 *   A3[3:M] -= A3[2]*A2[3:M] // GER frm 3rd step of LU
 * RETURNS: iamax of A3[2:M] AFTER updates applied
 */
{
   TYPE *A1 = A0 + lda, *A2 = A0 + lda+lda, *A3 = A2+lda;
   const register TYPE s1 = -A3[0];
   register TYPE s2 = A3[1], s3 = A3[2];
   register TYPE amax = ATL_rzero;
   register int i, imax=3;

   s2 += s1 * A0[1];
   A3[1] = s2;
   s2 = -s2;
   s3 += A0[2]*s1;
   s3 += A1[2]*s2;
   A3[2] = s3;
   s3 = -s3;
   for (i=3; i != M; i++)
   {
      register TYPE t0 = A0[i], t1 = A1[i], t2 = A2[i], t3 = A3[i];
      t2 *= s0;
      A2[i] = t2;
      t3 += t0*s1;
      t3 += t1*s2;
      t3 += t2*s3;
      A3[i] = t3;
      t3 = Mabs(t3);
      if (t3 > amax) { amax = t3; imax = i; }
   }
   return(imax);
}
#endif

int Mjoin(PATL,getrf2)(const int M, TYPE *A0, const int lda, int *ipiv)
/*
 * Factors 2 columns of LU, using left-looking algorithm with minimal
 * number of reads and writes by using special merged blas routines
 */
{
   int ip, iret=0;
   register TYPE t0;
   *ipiv = ip = cblas_iamax(M, A0, 1);
   t0 = A0[ip];
   if (M >= 2 && Mabs(t0) >= ATL_laSAFMIN)
   {
/*
 *    Apply pivot to both columns
 */
      TYPE *A1=A0+lda;
      register TYPE t1;
      if (ip)
      {
         t1 = A1[ip];
         A0[ip] = *A0;
         A1[ip] = *A1;
         *A0 = t0;
         *A1 = t1;
      }
/*
 *    Now, scale first column by 1/maxval, while doing a rank-1
 *    update on the second column
 */
      ipiv[1] = ip = irk1amax(M, ATL_rone/t0, A0, lda);
      t1 = A1[ip];
      if (t1 == ATL_rzero)
         return(2);
/*
 *    Apply pivot to both columns
 */
      if (ip != 1)
      {
         t0 = A0[ip];
         A0[ip] = A0[1];
         A1[ip] = A1[1];
         A0[1] = t0;
         A1[1] = t1;
      }
/*
 *    Scale last column by 1/maxval
 */
      if (Mabs(t1) >= ATL_laSAFMIN)
         Mjoin(PATL,scal)(M-2, ATL_rone/t1, A1+2, 1);
      else
      {
         register int i;
         for (i=2; i < M; i++)
            A1[i] /= t1;
      }
   }
   else
      return(Mjoin(PATL,getf2)(M, 2, A0, lda, ipiv));
   return(0);
}

int Mjoin(PATL,getrf3)(const int M, TYPE *A0, const int lda, int *ipiv)
/*
 * Factors 3 columns of LU, using left-looking algorithm with minimal
 * number of reads and writes by using special merged blas routines
 */
{
   int ip, iret=0;
   register TYPE t0;
   *ipiv = ip = cblas_iamax(M, A0, 1);
   t0 = A0[ip];
   if (M >= 3 && Mabs(t0) >= ATL_laSAFMIN)
   {
/*
 *    Apply pivot to all 3 columns
 */
      TYPE *A1=A0+lda, *A2=A0+lda+lda;
      register TYPE t1, t2, t3;
      if (ip)
      {
         t1 = A1[ip];
         t2 = A2[ip];
         A0[ip] = *A0;
         A1[ip] = *A1;
         A2[ip] = *A2;
         *A0 = t0;
         *A1 = t1;
         *A2 = t2;
      }
/*
 *    Now, scale remaining first column by 1/maxval, while doing a rank-1
 *    update on the second column
 */
      ipiv[1] = ip = irk1amax(M, ATL_rone/t0, A0, lda);
      t1 = A1[ip];
      if (Mabs(t1) >= ATL_laSAFMIN)
      {
/*
 *       Apply pivot to all 3 columns
 */
         if (ip != 1)
         {
            t0 = A0[ip];
            t2 = A2[ip];
            A0[ip] = A0[1];
            A1[ip] = A1[1];
            A2[ip] = A2[1];
            A0[1] = t0;
            A1[1] = t1;
            A2[1] = t2;
         }
/*
 *       Now, scale remaining second column by 1/maxval, while doing a rank-2
 *       update on the third column
 */
         ipiv[2] = ip = irk2amax(M, ATL_rone/t1, A0, lda);
         t2 = A2[ip];
         if (t2 == ATL_rzero)
            return(3);
/*
 *       Apply pivot to all 3 columns
 */
         if (ip != 2)
         {
            t0 = A0[ip];
            t1 = A1[ip];
            A0[ip] = A0[2];
            A1[ip] = A1[2];
            A2[ip] = A2[2];
            A0[2] = t0;
            A1[2] = t1;
            A2[2] = t2;
         }
/*
 *       Scale last column by 1/maxval
 */
         if (Mabs(t2) >= ATL_laSAFMIN)
            Mjoin(PATL,scal)(M-3, ATL_rone/t2, A2+3, 1);
         else
         {
            register int i;
            for (i=3; i < M; i++)
               A2[i] /= t2;
         }
      }
      else
      {
         iret = Mjoin(PATL,getf2)(M-1, 3, A1+1, lda, ipiv+1);
         if (iret)
            iret++;
         ipiv[1]++;
         ipiv[2]++;
         return(iret);
      }
   }
   else
      return(Mjoin(PATL,getf2)(M, 3, A0, lda, ipiv));
//   printf("IPIV = [%d, %d, %d, %d]\n", ipiv[0], ipiv[1], ipiv[2], ipiv[3]);
   return(0);
}
int Mjoin(PATL,getrf4)(const int M, TYPE *A, const int lda, int *ipiv)
/*
 * Factors 4 columns of LU, using left-looking algorithm with minimal
 * number of reads and writes by using special merged blas routines
 */
{
   int ip, iret=0;
   register TYPE t0;
   *ipiv = ip = cblas_iamax(M, A, 1);
   t0 = A[ip];
   if (M >= 4 && Mabs(t0) >= ATL_laSAFMIN)
   {
/*
 *    Apply pivot to all 4 columns
 */
      TYPE *A1=A+lda, *A2=A+lda+lda, *A3=A2+lda;
      register TYPE t1, t2, t3;
      if (ip)
      {
         t1 = A1[ip];
         t2 = A2[ip];
         t3 = A3[ip];
         A[ip] = *A;
         A1[ip] = *A1;
         A2[ip] = *A2;
         A3[ip] = *A3;
         *A = t0;
         *A1 = t1;
         *A2 = t2;
         *A3 = t3;
      }
/*
 *    Now, scale remaining first column by 1/maxval, while doing a rank-1
 *    update on the second column
 */
      ipiv[1] = ip = irk1amax(M, ATL_rone/t0, A, lda);
      t1 = A1[ip];
      if (Mabs(t1) >= ATL_laSAFMIN)
      {
/*
 *       Apply pivot to all 4 columns
 */
         if (ip != 1)
         {
            t0 = A[ip];
            t2 = A2[ip];
            t3 = A3[ip];
            A[ip] = A[1];
            A1[ip] = A1[1];
            A2[ip] = A2[1];
            A3[ip] = A3[1];
            A[1] = t0;
            A1[1] = t1;
            A2[1] = t2;
            A3[1] = t3;
         }
/*
 *       Now, scale remaining second column by 1/maxval, while doing a rank-2
 *       update on the third column
 */
         ipiv[2] = ip = irk2amax(M, ATL_rone/t1, A, lda);
         t2 = A2[ip];
         if (Mabs(t2) >= ATL_laSAFMIN)
         {
/*
 *          Apply pivot to all 4 columns
 */
            if (ip != 2)
            {
               t0 = A[ip];
               t1 = A1[ip];
               t3 = A3[ip];
               A[ip] = A[2];
               A1[ip] = A1[2];
               A2[ip] = A2[2];
               A3[ip] = A3[2];
               A[2] = t0;
               A1[2] = t1;
               A2[2] = t2;
               A3[2] = t3;
            }
/*
 *          Now, scale remaining third column by 1/maxval, while doing a rank-3
 *          update on the fourth column, and taking amax of it
 */
            ipiv[3] = ip = irk3amax(M, ATL_rone/t2, A, lda);
            t3 = A3[ip];
            if (t3 == ATL_rzero)
               return(4);
/*
 *          Apply pivot to all 4 columns
 */
            if (ip != 3)
            {
               t0 = A[ip];
               t1 = A1[ip];
               t2 = A2[ip];
               A[ip]  = A[3];
               A1[ip] = A1[3];
               A2[ip] = A2[3];
               A3[ip] = A3[3];
               A[3]  = t0;
               A1[3] = t1;
               A2[3] = t2;
               A3[3] = t3;
            }
/*
 *          Scale last column by 1/maxval
 */
            if (Mabs(t3) >= ATL_laSAFMIN)
               Mjoin(PATL,scal)(M-4, ATL_rone/t3, A3+4, 1);
            else if (t3 != ATL_rzero)
            {
               register int i;
               for (i=4; i < M; i++)
                  A3[i] /= t3;
            }
            else   /* pivot of 0 means rank-deficient matrix */
               return(4);
         }
         else
         {
            iret = Mjoin(PATL,getf2)(M-2, 2, A2+2, lda, ipiv+2);
            if (iret)
               iret += 2;
            ipiv[2] += 2;
            ipiv[3] += 2;
            return(iret);
         }
      }
      else
      {
         iret = Mjoin(PATL,getf2)(M-1, 3, A1+1, lda, ipiv+1);
         if (iret)
            iret++;
         ipiv[1]++;
         ipiv[2]++;
         ipiv[3]++;
         return(iret);
      }
   }
   else
      return(Mjoin(PATL,getf2)(M, 4, A, lda, ipiv));
//   printf("IPIV = [%d, %d, %d, %d]\n", ipiv[0], ipiv[1], ipiv[2], ipiv[3]);
   return(0);
}

