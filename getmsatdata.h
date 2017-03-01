#define GETMSATDATA_INCLUDE

#ifndef RECOMBINE_INCLUDE
#include "recombine.h"
#endif

void setupmsatdata(msatdata **ms, long numpop, long numloci,
   long numind);
/* OBSOLETE */
void freemsatdata(msatdata *ms);

void calculate_steps(msatdata *ms);
double logfac(long n);
#ifndef HAVE_LGAMMA
double lgamma(double z);
#endif
