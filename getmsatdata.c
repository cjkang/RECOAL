#ifndef GETMSATDATA_INCLUDE
#include "getmsatdata.h"
#endif

#ifdef DMEMDEBUG
#define MEMDEBUG
#include "memdebug.h"
#endif

#define NUMCHROM 2L

extern FILE *simlog, *weightfile;

/***********************************************************
 * setupmsatdata() initializes the microsat data structure */
void setupmsatdata(msatdata **ms, long numpop, long numloci,
   long numind)
{
long i, j, k;

(*ms) = (msatdata *)calloc(1,sizeof(msatdata));

(*ms)->title = (char *)calloc(LINESIZE,sizeof(char));
(*ms)->popnames = (char **)calloc(numpop,sizeof(char *));
(*ms)->popnames[0] = (char *)calloc(numpop*LINESIZE,sizeof(char));
for(i = 1; i < numpop; i++)
   (*ms)->popnames[i] = (*ms)->popnames[0] + i*LINESIZE;

(*ms)->numpop = numpop;
(*ms)->numloci = numloci;

(*ms)->numind = (long *)calloc(numpop,sizeof(long));
for(i = 0; i < numpop; i++) (*ms)->numind[i] = numind;

(*ms)->msats = (long ****)calloc(numpop,sizeof(long ***));
(*ms)->mspace = (double ***)calloc(numpop,sizeof(double **));
(*ms)->indnames = (char ****)calloc(numpop,sizeof(char ***));

(*ms)->msats[0] = (long ***)calloc(numpop*numind,sizeof(long **));
(*ms)->mspace[0] = (double **)calloc(numpop*numind,sizeof(double *));
(*ms)->indnames[0] = (char ***)calloc(numpop,sizeof(char **));
for(i = 1; i < numpop; i++) {
   (*ms)->msats[i] = (*ms)->msats[0] + i*numind;
   (*ms)->mspace[i] = (*ms)->mspace[0] + i*numind;
   (*ms)->indnames[i] = (*ms)->indnames[0] + i*numind;
}

(*ms)->msats[0][0] = 
   (long **)calloc(numpop*numind*numloci,sizeof(long *));
(*ms)->mspace[0][0] = (double *)calloc(numpop*numind*numloci,sizeof(double));
(*ms)->indnames[0][0] = (char **)calloc(numpop*numind*numloci,sizeof(char *));
for(i = 0; i < numpop; i++)
   for(j = 0; j < numind; j++) {
      (*ms)->msats[i][j] = (*ms)->msats[0][0] + i*numind*numloci +
         j*numloci;
      (*ms)->mspace[i][j] = (*ms)->mspace[0][0] + i*numind*numloci +
         j*numloci;
      (*ms)->indnames[i][j] = (*ms)->indnames[0][0] + i*numind*numloci +
         j*numloci;
   }

(*ms)->msats[0][0][0] = 
   (long *)calloc(numpop*numind*numloci*NUMCHROM,sizeof(long));
(*ms)->indnames[0][0][0] = 
   (char *)calloc(numpop*numind*numloci*(NMLNGTH+1),sizeof(char));
for(i = 0; i < numpop; i++)
   for(j = 0; j < numind; j++)
      for(k = 0; k < numloci; k++) {
          (*ms)->msats[i][j][k] = (*ms)->msats[0][0][0] +
             i*numind*numloci*NUMCHROM + j*numloci*NUMCHROM + k*NUMCHROM;
          (*ms)->indnames[i][j][k] = (*ms)->indnames[0][0][0] +
             i*numind*numloci*(NMLNGTH+1) + j*numloci*(NMLNGTH+1) +
             k*(NMLNGTH+1);
      }

(*ms)->steps = (double **)calloc(MICRO_MAXCHANGE,sizeof(double *));
(*ms)->steps[0] = (double *)calloc(MICRO_MAXCHANGE*MICRO_MAXCHANGE,
   sizeof(double));
for(i = 1; i < MICRO_MAXCHANGE; i++)
   (*ms)->steps[i] = (*ms)->steps[0] + i*MICRO_MAXCHANGE;

} /* setupmsatdata */


/****************************************************
 * freemsatdata() frees the microsat data structure */
void freemsatdata(msatdata *ms)
{
free(ms->msats[0][0][0]);
free(ms->msats[0][0]);
free(ms->msats[0]);
free(ms->msats);
free(ms->mspace[0][0]);
free(ms->mspace[0]);
free(ms->mspace);
free(ms->popnames[0]);
free(ms->popnames);
free(ms->title);
free(ms->steps[0]);
free(ms->steps);
free(ms);
} /* freemsatdata */


/******************************************************************
 * calculate_steps() initializes the table of microsat transition *
 * probabilities.                                                 */
void calculate_steps(msatdata *ms)
{
long k, diff;
long stepnum = MICRO_MAXCHANGE;

for (diff = 0; diff < stepnum; diff++) {
   for (k = diff; k < stepnum; k += 2) {
      ms->steps[diff][k-diff] = LOG2 * k +
      logfac((k-diff)/2)+logfac((k+diff)/2);
   }
}

} /* calculate_steps */


#ifndef HAVE_LGAMMA
double lgamma(double z)
{
    const double a[9] = { .9999999999995183,676.5203681218835,
            -1259.139216722289,771.3234287757674,-176.6150291498386,
            12.50734324009056,-.1385710331296526,9.934937113930748e-6,
           1.659470187408462e-7 };
    const double lnsqrt2pi = .9189385332046727;
    double result;

    long j;
    double  tmp;


    /*       Uses Lanczos-type approximation to ln(gamma) for z > 0. */
    /*       Reference:                                              */
    /*       Lanczos, C. 'A precision approximation of the gamma     */
    /*                    function', J. SIAM Numer. Anal., B, 1,     */
    /*                    86-96, 1964.                               */
    /*       Accuracy: About 14 significant digits except for small  */
    /*                    regions in the vicinity of 1 and 2.        */
    /*       Programmer: Alan Miller                                 */
    /*                   CSIRO Division of Mathematics & Statistics  */
    /*       Latest revision - 17 April 1988                         */
    /* translated and modified into C by Peter Beerli 1997           */


    if (z <= 0.) {
        return DBL_MAX; /*this will kill the receiving calculation*/
    }
    result = 0.;
    tmp = z + 7.;
    for (j = 9; j >= 2; --j) {
        result += a[j - 1] / tmp;
        tmp += -1.;
    }
    result += a[0];
    result = log(result) + lnsqrt2pi - (z + 6.5) + (z - .5) *
      log(z + 6.5);
    return result;
} /* lgamma */
#endif


/*****************************************************************
 * logfac() returns a value of log(n!), the first 12 values of n *
 * were calculated with Mathematica, 30 digits of precision.     */
double logfac(long n)
{
   switch (n) {
   case 0: return 0.;
   case 1: return 0.;
   case 2: return 0.693147180559945309417232121458;
   case 3: return 1.791759469228055000812477358381;
   case 4: return 3.1780538303479456196469416013;
   case 5: return 4.78749174278204599424770093452;
   case 6: return 6.5792512120101009950601782929;
   case 7: return 8.52516136106541430016553103635;
   case 8: return 10.60460290274525022841722740072;
   case 9: return 12.80182748008146961120771787457;
   case 10: return 15.10441257307551529522570932925;
   case 11: return 17.50230784587388583928765290722;
   case 12: return 19.98721449566188614951736238706;
   default: return lgamma(n + 1.);
   }
} /* logfac */
