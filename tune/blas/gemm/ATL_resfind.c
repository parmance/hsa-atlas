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
double GetAvg(int n, double tolerance, double *mflop)
{
   int i, j;
   double t0, tavg;
/*
 * Sort results, largest first
 */
   for (i=0; i != n; i++)
   {
      for (j=i+1; j < n; j++)
      {
         if (mflop[i] < mflop[j])
         {
            t0 = mflop[i];
            mflop[i] = mflop[j];
            mflop[j] = t0;
         }
      }
   }
/*
 * Not doing tolerance anymore, just take largest mflop rate if doing wall
 * times, or median value if doing CPU
 */

#if 1
   #ifdef WALL
      tavg = mflop[0];
   #else
      tavg = mflop[n/2];
   #endif
#else
/*
 * Throw out result if it is outside tolerance; rerun if two mflop not within
 * tolerance;  this code assumes n == 3
 */
   if (tolerance*mflop[1] < mflop[0])  /* too big a range in results */
   {
      if (tolerance*mflop[2] < mflop[1]) return(-1.0);
      tavg = (mflop[1] + mflop[2]) / 2.0;
   }
   else if (tolerance*mflop[2] < mflop[0]) tavg = (mflop[0] + mflop[1]) / 2.0;
   else tavg = (mflop[0] + mflop[1] + mflop[2]) / 3.0;
#endif

   return(tavg);
}

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
PrintUsage(char *nam)
{
   fprintf(stderr, "USAGE: %s <resfile1> .... <resfileN>\n", nam);
}
#define NTIM 3
#define TOLERANCE 1.2
double GetMflop(char *fnam)
{
   int i;
   double mflop[NTIM], t0;
   FILE *fp;

   fp = fopen(fnam, "r");
   assert(fp);
   for (i=0; i != NTIM; i++)
   {
      assert( fscanf(fp, "%lf", &mflop[i]) );
   }
   fclose(fp);

   t0 = GetAvg(NTIM, TOLERANCE, mflop);
   if (t0 == -1.0)
   {
      fprintf(stderr, "file=%s: rerun with higher reps; variation exceeds tolerence\n", fnam);
      exit(-1);
   }
   return(t0);
}

void ATL_resfind(FILE *fpout, int nfiles, char **files)
{
   char c1, pre, ln[64], *spc="           ";
   int i, j, k;
   int nb, mu, nu, ku, lat, muladd;
   double mflop;

   c1  = files[0][0];
   pre = files[0][1];
   fprintf(fpout, "%catlres = [\n", pre);
   for(i=0; i != nfiles; i++)
   {
      mflop = GetMflop(files[i]);
      assert(files[i][0] == c1  &&  files[i][1] == pre);
      assert(files[i][2] == 'N' && files[i][3] == 'B');

      for (j=4; isdigit(files[i][j]); j++) ln[j-4] = files[i][j];
      ln[j-4] = '\0';
      nb = atoi(ln);
      assert(files[i][j] == '_');

      for (j++, k=j; isdigit(files[i][j]); j++) ln[j-k] = files[i][j];
      ln[j-k] = '\0';
      mu = atoi(ln);
      assert(files[i][j] == 'x');

      for (j++, k=j; isdigit(files[i][j]); j++) ln[j-k] = files[i][j];
      ln[j-k] = '\0';
      nu = atoi(ln);
      assert(files[i][j] == 'x');

      for (j++, k=j; isdigit(files[i][j]); j++) ln[j-k] = files[i][j];
      ln[j-k] = '\0';
      ku = atoi(ln);
      assert(files[i][j] == '_');

      for (j++, k=j; isdigit(files[i][j]); j++) ln[j-k] = files[i][j];
      ln[j-k] = '\0';
      muladd = atoi(ln);
      assert(files[i][j] == '-');

      for (j++, k=j; isdigit(files[i][j]); j++) ln[j-k] = files[i][j];
      ln[j-k] = '\0';
      lat = atoi(ln);
      assert(files[i][j] == '.');

      fprintf(fpout, "%s%2d, %2d, %2d, %2d, %2d, %1d, %.2f;\n",
              spc, nb, mu, nu, ku, lat, muladd, mflop);
   }
   fprintf(fpout, "          ];\n");
}

void GetMesh(FILE *fpout, char pre)
{
   fprintf(fpout, "if (0)\n");
   fprintf(fpout, "   m = 0; n = 0;\n");
   fprintf(fpout, "   for k=1:length(datlres(:,1))\n");
      fprintf(fpout, "   if (datlres(k,2) > m)\n");
         fprintf(fpout, "   m = datlres(k,2);\n");
      fprintf(fpout, "   end\n");
      fprintf(fpout, "   if (datlres(k,3) > n)\n");
         fprintf(fpout, "   n = datlres(k,3);\n");
      fprintf(fpout, "   end\n");
   fprintf(fpout, "   end\n");
   fprintf(fpout, "   \n");
   fprintf(fpout, "   for i=1:m\n");
     fprintf(fpout, "   x(i) = i;\n");
   fprintf(fpout, "   end\n");
   fprintf(fpout, "   for j=1:n\n");
     fprintf(fpout, "   y(j) = j;\n");
   fprintf(fpout, "   end\n");
   fprintf(fpout, "   mf(1:m,1:n) = 0.0;\n");
   fprintf(fpout, "   \n");
   fprintf(fpout, "   ii = 1\n");
   fprintf(fpout, "   jj = 1\n");
   fprintf(fpout, "   for i=1:m\n");
   fprintf(fpout, "      for j=1:n\n");
   fprintf(fpout, "         for k=1:length(datlres(:,1))\n");
   fprintf(fpout, "            ii = datlres(k,2);\n");
   fprintf(fpout, "            jj = datlres(k,3);\n");
   fprintf(fpout, "            if (ii == i && jj == j)\n");
   fprintf(fpout, "               if (datlres(k,7) > mf(i,j))\n");
   fprintf(fpout, "                  mf(i,j) = datlres(k,7);\n");
   fprintf(fpout, "               end\n");
   fprintf(fpout, "            end\n");
   fprintf(fpout, "         end\n");
   fprintf(fpout, "      end\n");
   fprintf(fpout, "      grid \"on\"\n");
   fprintf(fpout, "      mesh(x, y, mf);\n");
   fprintf(fpout, "   end\n");
   fprintf(fpout, "end\n");
}

int main(int nargs, char **args)
{
   assert(nargs > 1);
   ATL_resfind(stdout, nargs-1, args+1);
   GetMesh(stdout, args[1][1]);
   return(0);
}
