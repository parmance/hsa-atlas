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
#include <stdlib.h>
#include <assert.h>
/*
 * Note that this program assumes:
 * (1) pagesize is the size of a page pointed to by a TLB entry
 * (2) Each TLB entry points to only 1 page
 */

FILE *fpout=NULL;
int *tlbinitmem(int n, size_t pglen)
/*
 * allocates and initializes relevant memory with values that won't overflow
 */
{
   int k=0;
   int *mem, *p, *stp;
   size_t len = n * pglen, ilen = pglen / sizeof(int);

   assert(n*ilen > 0);
   mem = p = malloc(len);
   assert(mem != NULL);
   stp = p + ilen * n;
   do
   {
      if (k == 0) *p = 0;
      else if (k == 1) *p = 1;
      else k = *p = -1;
      k++;
      p += ilen;
   }
   while(p != stp);

   return(mem);
}

double tlbaccesstime(int n, size_t pglen, int refs)
/*
 * Finds time to perform (refs) references, which are pglen apart.  We
 * allocate n separate pages
 */
{
   int *m, *mem, *p, *stp, reps;
   size_t len = pglen * n, ilen = pglen / sizeof(int);
   double tim;
   double time00();

   m = mem = tlbinitmem(n, pglen);
   reps = refs / n;
   stp = mem + n * ilen;
   tim = time00();
   do
   {
      p = mem;
      do
      {
         *m += *p;
         p += ilen;
      }
      while (p != stp);
      if (++m == stp) m = mem;
   }
   while(--reps);
   tim = time00() - tim;
   free(mem);
   return(tim);
}

int tlbnentries(size_t refs, int maxN, size_t maxpglen)
{
   int i, j, n=(-1);
   double tim1, tim2, diff, maxdiff=0.0;

   for (j=maxN; j >= 16; j /= 2)
   {
      tim1 = tlbaccesstime(j+1, maxpglen, refs);
      for (i=(-1); i < 6; i++)
      {
         tim2 = tlbaccesstime(j-i, maxpglen, refs);
         diff = tim1 - tim2;
         if (diff > maxdiff) {  maxdiff = diff;  n = j-i;  }
         printf("n=%d, time=%f, diff=%f\n", j-i, tim2, diff);
         if (fpout) fprintf(fpout, "n=%d, time=%f, diff=%f\n", j-i, tim2, diff);
         tim1 = tim2;
      }
   }
   return(n);
}

size_t tlbpglen(size_t refs, int nentries, int maxpglen)
/*
 * Given nentries in TLB, finds the page length
 */
{
   int n = nentries+8;
   double tim1, tim2, diff, maxdiff=0.0;
   size_t i, pglen=0;

   tim1 = tlbaccesstime(n, maxpglen, refs);
   tim1 = tlbaccesstime(n, maxpglen, refs);
   for (i=maxpglen/2; i > 256; i /= 2)
   {
      tim2 = tlbaccesstime(n, i, refs);
      diff = tim1 - tim2;
      if (diff > maxdiff) {  maxdiff = diff; pglen = i * 2;  }
      printf("len=%d, time=%f, diff=%f\n", i, tim2, diff);
      if (fpout) fprintf(fpout, "len=%d, time=%f, diff=%f\n", i, tim2, diff);
      tim1 = tim2;
   }
   return(pglen);
}

int main(int nargs, char *args)
{
   size_t refs=6400000, maxlen=512*1024, pglen;
   int maxn = 128, n;
   FILE *fp;

   fpout = fopen("res/TLBrun", "w");  /* stores run */
   assert (fpout != NULL);
#ifdef NoPageSize
   n = tlbnentries(refs, maxn, maxlen);
   printf("\n\n*****  Number of usable entries in TLB = %d\n\n", n);
   pglen = tlbpglen(refs, n, maxlen);
   printf("\n\n*****  Length of page in bytes         = %d\n\n", pglen);
#else
   pglen = getpagesize();
   printf("\n\n*****  Length of page in bytes         = %d\n\n", pglen);
   n = tlbnentries(refs, maxn, pglen);
   printf("\n\n*****  Number of usable entries in TLB = %d\n\n", n);
#endif
   fp = fopen("res/TLB", "w");
   assert (fp != NULL);
   fprintf(fp, "%d\n", n);
   fprintf(fp, "%d\n", pglen);
   fclose(fp);
   return(0);
}

