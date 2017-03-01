#include "rec_modellike.h"
#include "jdrop.h"
#include "world.h"

#ifdef DMEMDBUG
#define MEMDEBUG
#include "memdebug.h"
#endif

#ifdef DMALLOC_FUNC_CHECK
#include "/usr/local/include/dmalloc.h"
#endif

#define GOODUNDERFLOW
/*#undef GOODUNDERFLOW */

/***************************************************************************
 *  ALPHA                                                                  *
 *  version 1.00. (c) Copyright 1986, 1991, 1992 by the University of      *
 *  Washington and Joseph Felsenstein.  Written by Joseph Felsenstein,     *
 *  Mary K. Kuhner and Jon A. Yamato, with some additional grunt work by   *
 *  Sean T. Lamont.  Permission is granted to copy and use this program    *
 *   provided no fee is charged for it and provided that this copyright    *
 *  notice is not removed.                                                 *
 *                                                                         *
 ***************************************************************************/

/************************************************************
 This file contains all the functions used in calculating
 P(D|theta), the model_likelihood.
***********************************************************/

extern FILE *simlog, *outfile;
extern long population, locus, totchains, numtrees;
extern double watttheta, rec0;
extern tree *curtree;
extern boolean **sametree;

/* needed for output of chain results used in rec_estimate() */
extern long slid, slacc, hap, hacc, swap, swacc, numdropped, apps, indecks;

extern void histogram(cutpointrec dataToPlotHistogram);
extern int doublecmp(const void *v1, const void *v2);

extern recnumb *recnum0;

treerec ***sum;
double **reci, **theti;
double **reciH;
extern long numspaces;
extern double theta0;
extern long *spaceinfo;
extern long seq_length;
extern double *recarray;
double temptheta, temprecrates[2], templamda[2];
double *alpha0, *alpha1, *beta0, *beta1;


/**************************************************************
 * MYEXP() is a math::exp() wrapper that deals with underflow *
 * for machines that return NAN.                              */
#ifndef GOODUNDERFLOW
#define MYEXP(a) (((a) < EXPMIN) ? 0.0 : exp(a))
#else
#define MYEXP(a) exp(a)
#endif

/****************************************************************
 * zerocheck returns TRUE if value is zero, and FALSE otherwise */
boolean zerocheck(double value)
{
  if (value == 0.0) return(TRUE);
  else return(FALSE);
} /* zerocheck */

/******************************************************************
 * whatsign returns a 0 if value is zero, -1 if value is negative *
 * and 1 if value is positive.                                    */
long whatsign(double value)
{
  if (zerocheck(value)) return(0);
  if (value < 0.0) return(-1);
  else return(1);
} /* whatsign */

/*******************************************************************
 * makepositive halves the value of change until "value-change" is *
 * greater than 0, then that difference is returned.               */
double makepositive(double value, double change)
{

  while (value - change <= 0.0) change /= 2.0;

  return(value-change);

} /* makepositive */

/******************************************************************
 * changeparamLN assumes that everything is done in log parameter *
 * values!                                                        */
double changeparamLN(double value, double change)
{

  return(value * exp(-1.0 * change));

} /* changeparamLN */

/**********************************************************************
 * model_alloc allocates space for variables that will be used over   *
 * more than one locus in calculating the model_likelihood.           *
 * ONLY CALL ONCE.                                                    */
void model_alloc(option_struct *op, data_fmt *data)
{
  long i, j, k, numloci;

  numloci = getdata_numloci(op,data);

  sum = (treerec ***)calloc(1,numloci * sizeof(treerec **));
  sum[0] = (treerec **)calloc(1,numloci*(1+op->numchains[1]) *
                              sizeof(treerec *));
  for (i = 1; i < numloci; i++)
    sum[i] = sum[0] + i*(1+op->numchains[1]);
  sum[0][0] = (treerec *)calloc(1,numloci*(1+op->numchains[1])*numtrees *
				sizeof(treerec));
  for (i = 0; i < numloci; i++)
    for(j = 0; j < (1+op->numchains[1]); j++)
      sum[i][j] = sum[0][0] + i*(1+op->numchains[1])*numtrees + j*numtrees;
  /* the rest of the sum array is allocated with each tree */
  for (i = 0; i < numloci; i++)
    for (j = 0; j < 1+op->numchains[1]; j++)
      for (k = 0; k < numtrees; k++) {
        sum[i][j][k].kk = NULL;
        sum[i][j][k].kend = NULL;
        sum[i][j][k].eventtype = NULL;
        sum[i][j][k].actives = NULL;
        sum[i][j][k].numcoals = 0;
        sum[i][j][k].numrecombs = 0;
      }

  reci = (double **)calloc(1,numloci * sizeof(double *));
  theti = (double **)calloc(1,numloci * sizeof(double *));
  reci[0] = (double *)
    calloc(1,numloci*(totchains+1) * sizeof(double));
  theti[0] = (double *)
    calloc(1,numloci*(totchains+1) * sizeof(double));
  for (i = 1; i < numloci; i++) {
    reci[i] = reci[0] + i*(totchains+1);
    theti[i] = theti[0] + i*(totchains+1);
  }
}  /* model_alloc */

/*********************************************************
 * count_active_tlist returns the number of active links *
 * within a time-slice of the tree.                      */
double count_active_tlist(option_struct *op, data_fmt *data, tlist *t)
{
  long i;
  double accum;

  accum = 0.0;
  for (i = 0; i < t->numbranch; i++) 
    if (op->fc) accum += count_activefc(op,data,t->branchlist[i]->back);
    else accum += count_active(op,data,t->branchlist[i]->back);

  return (accum);

} /* count_active_tlist */


/***********************************************************
 * scoretree saves the values necessary for evaluating the *
 * model_likelihood for a single tree in the 'sum' array   */
void scoretree(option_struct *op, data_fmt *data, long chain)
{
  tlist *t;
  treerec *trii;
  long i, j, refchain, chaintype, entries;
  double temp;

  refchain = REF_CHAIN(chain);
  chaintype = TYPE_CHAIN(chain);

  temp = 0.0;
  t = curtree->tymelist;
  /* count the tymelist entries */
  entries = 0;
  while (t!=NULL) {
    entries++;
    t = t->succ;
  }
  t = curtree->tymelist;

  /* allocate space for them */
  trii = &sum[locus][refchain][op->numout[chaintype]-1];

#if GROWTHUSED
  if (trii->kk) free(trii->kk);
  trii->kk = (long *)calloc(1,entries*sizeof(long));
  if(trii->kend) free(trii->kend);
  trii->kend = (double *)calloc(1,entries*sizeof(double));
  if (trii->actives) free(trii->actives);
  trii->actives = (double *)calloc(1,entries*sizeof(double));
#endif
  if(trii->eventtype) free(trii->eventtype);
  trii->eventtype = (long *)calloc(1,entries*sizeof(long));
  if (trii->sitescore) free(trii->sitescore);
  /* the +1 is to guarantee a minimum allocation of size 1 */
  trii->sitescore = (long *)calloc(curtree->numrecombs+1,sizeof(long));

  trii->tk = 0.0;
  trii->ts = 0.0;

  for(i=0,j=0,temp=0.0;i<entries;i++) {
    if (t->numbranch == 1) break;
#if GROWTHUSED
    trii->kk[i] = t->numbranch * (t->numbranch - 1);
    trii->kend[i] = t->age;
    trii->actives[i] = count_active_tlist(op,data,t);
#endif
    trii->eventtype[i] = isrecomb(t->eventnode);
    temp += t->numbranch * (t->numbranch - 1) * (t->age - t->eventnode->tyme);
    trii->tk = temp;
    trii->ts += count_active_tlist(op,data,t) * (t->age - t->eventnode->tyme);
    if (trii->eventtype[i]) {
      trii->sitescore[j] = findlink(t->eventnode);
      j++;
    }
    t = t->succ;
  }

  trii->numcoals = curtree->numcoals;
  trii->numrecombs = curtree->numrecombs;
  trii->llike = curtree->likelihood;

#if !ALWAYS_REJECT
  if (temp == 0.0) fprintf(ERRFILE,"WARNING:  Tree has become length zero\n");
#endif
  if (trii->numcoals - trii->numrecombs != getdata_numtips(op,data) - 1) {
    fprintf(ERRFILE,"ERROR:scoretree says: bad nodes in the tree!\n");
    exit(-1);
  }

}  /* scoretree */


double lastintervalcalc(option_struct *op, data_fmt *data, long chain)
{
  treerec *trii;
  long i, refchain, chaintype, numintervals;
  double endtyme, prevendtyme, intsum = 0.0;

  refchain = REF_CHAIN(chain);
  chaintype = TYPE_CHAIN(chain);

  for(i = 0; i < op->numout[chaintype]; i++) {
    trii = &sum[locus][refchain][i];
    numintervals = trii->numcoals + trii->numrecombs;
    prevendtyme = (numintervals > 1) ? 
      trii->kend[numintervals-2] : 0.0;
    endtyme = trii->kend[numintervals-1];
    intsum += endtyme - prevendtyme;
  }

  return(intsum/op->numout[chaintype]);

} /* lastintervalcalc */


void intervalcheck(option_struct *op, data_fmt *data, long chain)
{
  long i, numintervals, refchain, chaintype;
  treerec *trii;

  refchain = REF_CHAIN(chain);
  chaintype = TYPE_CHAIN(chain);

  for(i = 0; i < op->numout[chaintype]; i++) {
    trii = &sum[locus][refchain][i];
    if (trii->numrecombs != 1) continue;
    numintervals = trii->numcoals + trii->numrecombs;
  }

} /* intervalcheck */

/****************************************************************
 * rec_scoreprint prints the contents of the sum array out into *
 * 6 files for use in Mathematica.  Note that this function     *
 * only works for the recombination case.                       */
void rec_scoreprint(option_struct *op, data_fmt *data, long lowcus,
		    long chain, boolean mathematica)
{
  long i, lowc, lowcstart, lowcend, refchain, numloci, chaintype;
#if GROWTHUSED
  long j, numintervals;
#endif
  char openchar, closechar, dlmchar;
  FILE *starts, *ends, *kks, *actives, *numrecombs, *numcoals;

  starts = fopen("starts","w+");
  ends = fopen("ends","w+");
  kks = fopen("kks","w+");
  actives = fopen("actives","w+");
  numrecombs = fopen("numrecombs","w+");
  numcoals = fopen("numcoals","w+");

  refchain = REF_CHAIN(chain);
  chaintype = TYPE_CHAIN(chain);
  numloci = getdata_numloci(op,data);

  if (lowcus < 0) {
    lowcstart = 0;
    lowcend = numloci;
  } else {
    lowcstart = lowcus;
    lowcend = lowcus + 1;
  }

  if(mathematica) {
    openchar = '{';
    closechar = '}';
    dlmchar = ',';
    fprintf(starts,"stv = {\n");
    fprintf(ends,"etv = {\n");
    fprintf(kks,"kkv = {\n");
    fprintf(actives,"activev = {\n");
    fprintf(numrecombs,"numrecv = {\n");
    fprintf(numcoals,"numcoalv = {\n");
  } else {
    openchar = closechar = dlmchar = ' ';
    fprintf(kks,"%ld %ld %ld\n",lowcend-lowcstart,totchains,
	    op->numout[chaintype]);
  }

  for(lowc = lowcstart; lowc < lowcend; lowc++) {
    fprintf(starts,"%c",openchar);
    fprintf(ends,"%c",openchar);
    fprintf(kks,"%c",openchar);
    fprintf(actives,"%c",openchar);
    fprintf(numrecombs,"%c",openchar);
    fprintf(numcoals,"%c",openchar);
    for(i = 0; i < op->numout[chaintype]; i++) {
      fprintf(starts,"%c",openchar);
      fprintf(ends,"%c",openchar);
      fprintf(kks,"%c",openchar);
      fprintf(actives,"%c",openchar);
      fprintf(numrecombs,"%18.16f",sum[lowc][refchain][i].numrecombs);
      fprintf(numcoals,"%ld",sum[lowc][refchain][i].numcoals);
      fprintf(kks,"%18.16f",sum[lowc][refchain][i].tk);
      fprintf(actives,"%18.16f",sum[lowc][refchain][i].ts);
#if GROWTHUSED
      numintervals = sum[lowc][refchain][i].numrecombs +
	sum[lowc][refchain][i].numcoals;
      for(j = 0; j < numintervals; j++) {
	if (!j) fprintf(starts,"0.0");
	else
	  fprintf(starts,"%18.16f",sum[lowc][refchain][i].kend[j-1]);
	fprintf(ends,"%18.16f",sum[lowc][refchain][i].kend[j]);
	fprintf(kks,"%ld",sum[lowc][refchain][i].kk[j]);
	fprintf(actives,"%18.6f",sum[lowc][refchain][i].actives[j]);
	if (j != numintervals - 1) {
	  fprintf(starts,"%c",dlmchar);
	  fprintf(ends,"%c",dlmchar);
	  fprintf(kks,"%c",dlmchar);
	  fprintf(actives,"%c",dlmchar);
	}
      }
#endif
      if (i != op->numout[chaintype]-1) {
	fprintf(starts,"%c%c\n",closechar,dlmchar);
	fprintf(ends,"%c%c\n",closechar,dlmchar);
	fprintf(kks,"%c%c\n",closechar,dlmchar);
	fprintf(actives,"%c%c\n",closechar,dlmchar);
	fprintf(numrecombs,"%c",dlmchar);
	fprintf(numcoals,"%c",dlmchar);
      } else {
	fprintf(starts,"%c",closechar);
	fprintf(ends,"%c",closechar);
	fprintf(kks,"%c",closechar);
	fprintf(actives,"%c",closechar);
      }
    }
    if (lowc != lowcend-1) {
      fprintf(starts,"%c%c",closechar,dlmchar);
      fprintf(ends,"%c%c",closechar,dlmchar);
      fprintf(kks,"%c%c",closechar,dlmchar);
      fprintf(actives,"%c%c",closechar,dlmchar);
      fprintf(numrecombs,"%c%c",closechar,dlmchar);
      fprintf(numcoals,"%c%c",closechar,dlmchar);
    } else {
      fprintf(starts,"%c",closechar);
      fprintf(ends,"%c",closechar);
      fprintf(kks,"%c",closechar);
      fprintf(actives,"%c",closechar);
      fprintf(numrecombs,"%c",closechar);
      fprintf(numcoals,"%c",closechar);
    }
  }

  fprintf(starts,"%c",closechar);
  fprintf(ends,"%c",closechar);
  fprintf(kks,"%c",closechar);
  fprintf(actives,"%c",closechar);
  fprintf(numrecombs,"%c",closechar);
  fprintf(numcoals,"%c",closechar);

  if (mathematica) fprintf(kks,"\n\ntheti = {");
  else fprintf(kks,"\n\n");
  for(lowc = lowcstart; lowc < lowcend; lowc++) {
    fprintf(kks,"%c",openchar);
    for(i = 0; i < totchains; i++) {
      fprintf(kks,"%18.16f",theti[lowc][i]);
      if (i != totchains-1) fprintf(kks,"%c",dlmchar);
    }
    fprintf(kks,"%c",closechar);
    if (lowc != lowcend-1) fprintf(kks,"%c\n",dlmchar);
  }
  fprintf(kks,"%c\n",closechar);

  if (mathematica) fprintf(kks,"\n\nreci = {");
  else fprintf(kks,"\n\n");
  for(lowc = lowcstart; lowc < lowcend; lowc++) {
    fprintf(kks,"%c",openchar);
    for(i = 0; i < totchains; i++) {
      fprintf(kks,"%18.16f",reci[lowc][i]);
      if (i != totchains-1) fprintf(kks,"%c",dlmchar);
    }
    fprintf(kks,"%c",closechar);
    if (lowc != lowcend-1) fprintf(kks,"%c\n",dlmchar);
  }
  fprintf(kks,"%c\n",closechar);

  fclose(kks);
  fclose(ends);
  fclose(starts);
  fclose(actives);
  fclose(numrecombs);
  fclose(numcoals);

} /* rec_scoreprint */


/*****************************************************************
 * rec_readscoretree() reads in a single tree to the "sum" array */
void rec_readscoretree(FILE *numcoals, FILE *numrecombs, FILE *kks,
		       FILE *ends, FILE *actives, option_struct *op,
		       long lowcus, long chain, long trii)
{
  long refchain;
#if GROWTHUSED
  long i, numintervals;
#endif
  treerec *tr;

  refchain = REF_CHAIN(chain);

  tr = &sum[lowcus][refchain][trii];

  fscanf(numrecombs,"%lf",&(tr->numrecombs));
  fscanf(numcoals,"%ld",&(tr->numcoals));

  fscanf(kks,"%lf",&(tr->tk));
  fscanf(actives,"%lf",&(tr->ts));
#if GROWTHUSED
  numintervals = tr->numrecombs + tr->numcoals;

  tr->kk = (long *)calloc(numintervals,sizeof(long));
  tr->kend = (double *)calloc(numintervals,sizeof(double));
  tr->actives = (double *)calloc(numintervals,sizeof(double));

  for(i = 0; i < numintervals; i++) {
    fscanf(ends,"%lf",&(tr->kend[i]));
    fscanf(kks,"%ld",&(tr->kk[i]));
    fscanf(actives,"%lf",&(tr->actives[i]));
  }
#endif


} /* rec_readscoretree */


/****************************************************************
 * rec_scoreread() reads the files created by rec_scoreprint to *
 * recreate the treesummary array, "sum".  rec_scoreprint needs *
 * to be passed "mathematica = FALSE".                          */
void rec_scoreread(option_struct *op, data_fmt *data)
{
  long lowcus, chain, trii, numtriis, numchains, numlowci;
  FILE *ends, *kks, *actives, *numrecombs, *numcoals;

  ends = fopen("ends","r");
  kks = fopen("kks","r");
  actives = fopen("actives","r");
  numrecombs = fopen("numrecombs","r");
  numcoals = fopen("numcoals","r");

  fscanf(kks,"%ld %ld %ld",&numlowci,&numchains,&numtriis);

  op->numout[1] = numtriis;

  chain = numchains-1;
  for(lowcus = 0; lowcus < numlowci; lowcus++)
    for(trii = 0; trii < numtriis; trii++)
      rec_readscoretree(numcoals,numrecombs,kks,ends,
			actives,op,lowcus,chain,trii);

  for(lowcus = 0; lowcus < numlowci; lowcus++)
    for(chain = 0; chain < numchains; chain++)
      fscanf(kks,"%lf",&theti[lowcus][chain]);

  for(lowcus = 0; lowcus < numlowci; lowcus++)
    for(chain = 0; chain < numchains; chain++)
      fscanf(kks,"%lf",&reci[lowcus][chain]);

  fclose(ends);
  fclose(kks);
  fclose(actives);
  fclose(numrecombs);
  fclose(numcoals);

} /* rec_scoreread */


#define LNUMTH 10000
#define LTHMAX 10.0
#define LTHMIN 0.001
#define LNUMNUM 10 /* scaled in # of numbers per line */

/****************************************************************
 * rec_likeprint() prints out a set of ordered pairs for use by *
 * Mathematica to construct a likelihood curve                  */
void rec_likeprint(option_struct *op, data_fmt *data, long lowcus,
		   long chain)
{
  long i, linebreak;
  double llike[LNUMTH], theta;
  FILE *likefile;

  theta = LTHMAX;
  linebreak = LNUMNUM;
  likefile = fopen("likefile","w+");

  if (op->progress)
    printf("\ncalculating the exact curve for printing\n");

  fprintf(likefile,"likecurve = {");
  for(i = 0; i < LNUMTH; i++) {
    llike[i] = model_likelihood(op,data,theta,0.0,lowcus,chain);
    fprintf(likefile,"{%12.6f,%12.6f}",log(theta),llike[i]);
    theta -= (LTHMAX - LTHMIN) / LNUMTH;
    if (i != LNUMTH-1) fprintf(likefile,",");
    if (i == linebreak) {
      fprintf(likefile,"\n");
      linebreak += LNUMNUM;
    }
  }
  fprintf(likefile,"}");

  fclose(likefile);

} /* rec_likeprint */

#define LNUMPTS 1000
#define LRECMAX 10.0
#define LRECMIN 0.0

/*********************************************************************
 * rec_likeprint2() prints out likelihood values at the max in theta *
 * for use by Mathematica to construct a likelihood curve            */
void rec_likeprint2(option_struct *op, data_fmt *data, double thmax, 
		    long lowcus, long chain)
{
  long i, linebreak;
  double llike[LNUMPTS], rec;
  FILE *likefile;

  rec = LRECMAX;
  linebreak = LNUMNUM;
  likefile = fopen("thmaxlikefile","w+");

  if (op->progress)
    printf("\ncalculating the exact curve for printing\n");

  fprintf(likefile,"likecurve = {");
  for(i = 0; i < LNUMPTS; i++) {
    llike[i] = model_likelihood(op,data,thmax,rec,lowcus,chain);
    fprintf(likefile,"{%12.6f,%12.6f}",rec,llike[i]);
    rec -= (LRECMAX - LRECMIN) / LNUMPTS;
    if (i != LNUMPTS-1) fprintf(likefile,",");
    if (i == linebreak) {
      fprintf(likefile,"\n");
      linebreak += LNUMNUM;
    }
  }
  fprintf(likefile,"}");

  fclose(likefile);

} /* rec_likeprint2 */


/************************************************************************
 * stuff returns the exponent of the waiting time for the likelihood of *
 * a coalescence.                                                       *
 *                                                                      *
 * Eqn: kk = (#active lineages) * (#active lineages - 1)                *
 *       s = active site count over all active lineages                 *
 *       t = length of interval                                         *
 *                                                                      *
 * stuff = Sum_over_all_intervals[-t * ((kk/theta)+(rec*s))]            *
 *                                                                      */
double stuff(double theta, double rec, treerec *tr)
{
  double temp;
#if GROWTHUSED
  double starttyme, endtyme;
  long i, numintervals;

  temp = 0.0;

  numintervals = tr->numcoals + tr->numrecombs;

  for(i = 0; i < numintervals; i++) {
    if (i == 0) starttyme = 0;
    else starttyme = tr->kend[i - 1];
    endtyme = tr->kend[i];

    temp += (starttyme - endtyme) * 
      (tr->kk[i]/theta + rec*tr->actives[i]);
  }
#endif

  temp = -tr->tk/theta - rec*tr->ts;

  return(temp);

} /* stuff */
 
/************************************************************************
 * rectreellike returns the Ln(Likelihood) of a tree with recombination */
double rectreellike(option_struct *op, double theta, double rec,
		    double lth, double lrec, long lowcus, long chain, long trii)
{
  long refchain;
  double temp;
  treerec *tr;

  refchain = REF_CHAIN(chain);

  tr = (&sum[lowcus][refchain][trii]);

  if (tr->numrecombs && !rec) return(NEGMAX);

  if (rec)
    temp = tr->numcoals * lth + tr->numrecombs * lrec + 
      stuff(theta, rec, tr);
  else
    temp = tr->numcoals*lth + stuff(theta,rec,tr);

  return(temp);

} /* rectreellike */

/**************************************************************************
 * recchainllike returns the Ln(Likelihood) of a chain with recombination */
double recchainllike(option_struct *op, double theta, double rec,
		     long lowcus, long chain)
{
  long i, chaintype, refchain;
  double *triillike, max, answ, lth, lrec, lth0, lrec0, th0, r0;

  chaintype = TYPE_CHAIN(chain);
  refchain = REF_CHAIN(chain);

  triillike = (double *)calloc(1,op->numout[chaintype] * sizeof(double)); 

  /* program speed ups */
  th0 = theti[lowcus][chain];
  r0 = reci[lowcus][chain];
  lth = log(2.0/theta);
  if (rec) lrec = log(rec);
  else lrec = log(epsilon);
  lth0 = log(2.0/th0);
  if (r0) lrec0 = log(r0);
  else lrec0 = log(epsilon);

  max = NEGMAX;
  for(i = 0; i < op->numout[chaintype]; i++) {
    if (i == 0 || !sametree[lowcus][i]) {
      triillike[i] = rectreellike(op,theta,rec,lth,lrec,lowcus,chain,i) -
	rectreellike(op,th0,r0,lth0,lrec0,lowcus,chain,i);
      if (triillike[i] > max) max = triillike[i]; 
    } else triillike[i] = triillike[i - 1];
  }

  answ = 0.0;
  for(i = 0; i < op->numout[chaintype]; i++) {
    if (sum[lowcus][refchain][i].numrecombs && !rec) continue;
    if (triillike[i] - max > EXPMIN) answ += exp(triillike[i] - max);
  }

  if (answ) answ = log(answ) + max - log(op->numout[chaintype]);
  else answ = NEGMAX;

  free(triillike);

  return(answ);

} /* recchainllike */


/****************************************************************
 * rec_locusllike returns the Ln(Likelihood) over all loci with *
 * recombination.                                               */
double rec_locusllike(option_struct *op, data_fmt *data,
		      double theta, double rec)
{
  long lowcus, lastchain, numloci;
  double answ, temp;

  lastchain = totchains - 1;
  numloci = getdata_numloci(op,data);

  answ = 0.0;
  for(lowcus = 0; lowcus < numloci; lowcus++) {
    temp = recchainllike(op,op->ne_ratio[lowcus]*theta,rec,lowcus,
			 lastchain);
    if (temp != NEGMAX) answ += temp;
    else {answ = NEGMAX; break;}
  }

  return(answ);

} /* rec_locusllike */

/******************************************************************
 * model_likelihood is the driver for the likelihood calculations */
double model_likelihood(option_struct *op, data_fmt *data,
			double theta, double rec, long lowcus, long chain)
{
  if (lowcus == -1) return(rec_locusllike(op,data,theta,rec));
  else return(recchainllike(op,theta,rec,lowcus,chain));

} /* model_likelihood */

/********************************************************************
 * frec_chainllike returns the Ln(Likelihood) ratio of a chain with *
 * recombination.  This is the "fast" version.                      */
double frec_chainllike(option_struct *op, double theta, double rec,
		       double **denom, long lowcus, long chain)
{
  long i, chaintype, refchain;
  double *triillike, max, answ, lth, lrec;

  chaintype = TYPE_CHAIN(chain);
  refchain = REF_CHAIN(chain);

  triillike = (double *)calloc(1,op->numout[chaintype] * sizeof(double));

  /* program speed ups */
  lth = log(2.0/theta);
  if (rec) lrec = log(rec);
  else lrec = log(epsilon);

  max = NEGMAX;
  for(i = 0; i < op->numout[chaintype]; i++) {
    if (i == 0 || !sametree[lowcus][i]) {
      triillike[i] = rectreellike(op,theta,rec,lth,lrec,lowcus,chain,i) -
	denom[lowcus][i];
      if (triillike[i] > max) max = triillike[i];
    } else triillike[i] = triillike[i - 1];
  }

  answ = 0.0;
  for(i = 0; i < op->numout[chaintype]; i++) {
    if (sum[lowcus][refchain][i].numrecombs && !rec) continue;
    if (triillike[i] - max > EXPMIN) answ += exp(triillike[i] - max);
  }

  if (answ) answ = log(answ) + max - log(op->numout[chaintype]);
  else answ = NEGMAX;

  free(triillike);

  return(answ);

} /* frec_chainllike */

/*****************************************************************
 * frec_locusllike returns the Ln(Likelihood) over all loci with *
 * recombination.  This is the "fast" version.                   */
double frec_locusllike(option_struct *op, data_fmt *data,
		       double theta, double rec, double **denom)
{
  long lowcus, lastchain, numloci;
  double answ, temp;

  lastchain = totchains - 1;
  numloci = getdata_numloci(op,data);

  answ = 0.0;
  for(lowcus = 0; lowcus < numloci; lowcus++) {
    temp = frec_chainllike(op,op->ne_ratio[lowcus]*theta,rec,
			   denom,lowcus,lastchain);
    if (temp != NEGMAX) answ += temp;
    else {answ = NEGMAX; break;}
  }

  return(answ);

} /* frec_locusllike */

/*******************************************************************
 * fmodel_likelihood is the driver for the likelihood calculations *
 * This version is speeded up by passing in the denominator dealing*
 * with the LnL under which the chain was run.                     */
double fmodel_likelihood(option_struct *op, data_fmt *data,
			 double theta, double rec, double **denom, long lowcus, long chain)
{
  if (lowcus == -1) return(frec_locusllike(op,data,theta,rec,denom));
  else return(frec_chainllike(op,theta,rec,denom,lowcus,chain));

} /* fmodel_likelihood */

/***************************************************************
 * dstuff is a helper function that calculates either theta or *
 * recombine relevant values.                                  */
double dstuff(double theta, treerec *tr)
{
  double answ;

#if GROWTHUSED
  long i, numintervals;
  double starttyme, endtyme;

  i = 0;
  answ = 0.0;

  numintervals = tr->numcoals + tr->numrecombs;

  for(i = 0; i < numintervals; i++) {
    if (i == 0) starttyme = 0;
    else starttyme = tr->kend[i - 1];
    endtyme = tr->kend[i];

    if (!theta) answ += (starttyme - endtyme) * tr->actives[i];
    else answ += (endtyme - starttyme) * tr->kk[i] / (theta * theta);
  }
#endif

  if (!theta) answ = -1.0 * tr->ts;
  else answ = tr->tk/(theta*theta);

  return (answ);

} /* dstuff */

/**********************************************************************
 * rec_thetalderiv returns the ratio of the 1st to the 2nd derivative *
 * in theta.                                                          *
 *                                                                    *
 * CHANGED foR B-F-G-S to only do 1st derivative!                     *
 *                                                                    *
 * Eqn: kk = (#active lineages) * (#active lineages - 1)              *
 *       c = #coalescences present in a given tree                    *
 *       s = active site count over all active lineages               *
 *       t = length of interval                                       *
 *       L = Prior Likelihood under theta and rec                     *
 *    SQ() = square of the arguments given                            *
 *                                                                    *
 * dstuff = Sum_over_all_intervals[t * kk/(theta*theta)]              *
 *                                                                    *
 * 1st derivative = sum_over_all_trees[L * (dstuff - (c/theta))]      *
 *                                                                    *
 * 2nd derivative = sum_over_all_trees of                             *
 *    L * [SQ(dstuff - (c/theta)) + (c/SQ(theta)) - (2/theta)*dstuff] *
 *                                                                    */
double rec_thetalderiv(option_struct *op, double theta, double rec,
		       double **llike0, long lowcus, long chain, double *dfn, long *dfnplus)
{
  long trii, refchain, chaintype, numc;
  double ds, ds2, fn, lth, lrec, temp;
  treerec *tr;

  refchain = REF_CHAIN(chain);
  chaintype = TYPE_CHAIN(chain);
  lth = log(2.0/theta);
  lrec = log(rec);

  *dfn = 0.0;

  for(trii = 0; trii < op->numout[chaintype]; trii++) {
    tr = (&sum[lowcus][refchain][trii]);
    if (tr->numrecombs && !rec) continue;
    ds = dstuff(theta, tr);
    numc = tr->numcoals;
    ds2 = ds - numc/theta;
    fn = rectreellike(op,theta,rec,lth,lrec,lowcus,chain,trii) -
      llike0[lowcus][trii];

    temp = exp(fn) * ds2;
    (*dfn) += exp(fn) * ds2;
  }

  *dfnplus = whatsign(*dfn);

  if (*dfn == 0.0) *dfn = epsilon;

  *dfn = log(fabs(*dfn)) - log(op->numout[chaintype]);

  return((*dfn));

} /* rec_thetalderiv */

/********************************************************************
 * rec_reclderiv returns the ratio of the 1st to the 2nd derivative *
 * in recombination rate.                                           *
 *                                                                  *
 * CHANGED foR B-F-G-S to only do 1st derivative!                   *
 *                                                                  *
 * Eqn: kk = (#active lineages) * (#active lineages - 1)            *
 *       r = #recombinations present in a given tree                *
 *       s = active site count over all active lineages             *
 *       t = length of interval                                     *
 *       L = Prior Likelihood under theta and rec                   *
 *    SQ() = square of the arguments given                          *
 *                                                                  *
 * dstuff = Sum_over_all_intervals[-t * s]                          *
 *                                                                  *
 * 1st derivative = sum_over_all_trees[L * (dstuff + (r/rec))]      *
 *                                                                  *
 * 2nd derivative = sum_over_all_trees of                           *
 *    L * [ SQ(dstuff + (r/rec)) - (r/SQ(rec))]                     *
 *                                                                  */
double rec_reclderiv(option_struct *op, double theta, double rec,
		     double **llike0, long lowcus, long chain, double *dfn, long *dfnplus)
{
  long trii, refchain, chaintype, numr;
  double ds, ds2, fn, lth, lrec, temp;
  treerec *tr;

  refchain = REF_CHAIN(chain);
  chaintype = TYPE_CHAIN(chain);
  lth = log(2.0/theta);
  lrec = log(rec);

  *dfn = 0.0;

  for(trii = 0; trii < op->numout[chaintype]; trii++) {
    tr = (&sum[lowcus][refchain][trii]);
    if (tr->numrecombs && !rec) continue;
    numr = tr->numrecombs;
    ds = dstuff(0.0, tr);
    ds2 = ds + numr/rec;
    fn = rectreellike(op,theta,rec,lth,lrec,lowcus,chain,trii) -
      llike0[lowcus][trii];

    temp = exp(fn) * ds2;
    (*dfn) += exp(fn) * ds2;
  }

  *dfnplus = whatsign(*dfn);

  if (*dfn == 0.0) *dfn = epsilon;

  *dfn = log(fabs(*dfn)) - log(op->numout[chaintype]);

  return((*dfn));

} /* rec_reclderiv */


/******************************************************
 * rec_locus_thetalderiv computes the multi-locus 1st *
 * derivatives wrt. theta                             */
double rec_locus_thetalderiv(option_struct *op, data_fmt *data,
			     double theta, double rec, double **llike0, double *dfx, long *dfxplus)
{
  long lowcus, lfxplus, lastchain, numloci;
  double rtheta, fntheta, dfntheta, dummy, temp1, temp2;

  lastchain = totchains - 1;
  numloci = getdata_numloci(op,data);
  temp1 = 0.0;
  temp2 = 0.0;

  for (lowcus = 0; lowcus < numloci; lowcus++) {
    rtheta = theta * op->ne_ratio[lowcus];
    fntheta =
      fmodel_likelihood(op,data,rtheta,rec,llike0,lowcus,lastchain);
    dummy =
      rec_thetalderiv(op,rtheta,rec,llike0,lowcus,lastchain,&dfntheta,
		      &lfxplus);

    temp1 += fntheta;
    if (lfxplus)
      temp2 += op->ne_ratio[lowcus] * lfxplus * exp(dfntheta - fntheta);
  }

  *dfxplus = whatsign(temp2);

  if (!*dfxplus) *dfx = 0.0;
  else *dfx = temp1 + log(fabs(temp2));

  return(*dfx);

} /* rec_locus_thetalderiv */


double rec_locus_reclderiv(option_struct *op, data_fmt *data, double theta,
			   double rec, double **llike0, double *dfx, long *dfxplus)
{
  long lowcus, lfxplus, lastchain, numloci;
  double rtheta, fnrec, dfnrec, dummy, temp1, temp2;

  lastchain = totchains - 1;
  numloci = getdata_numloci(op,data);
  temp1 = 0.0;
  temp2 = 0.0;

  for (lowcus = 0; lowcus < numloci; lowcus++) {
    rtheta = theta * op->ne_ratio[lowcus];
    fnrec =
      fmodel_likelihood(op,data,rtheta,rec,llike0,lowcus,lastchain);
    dummy = rec_reclderiv(op,rtheta,rec,llike0,lowcus,lastchain,&dfnrec,
			  &lfxplus);

    temp1 += fnrec;
    if (lfxplus)
      temp2 += lfxplus * exp(dfnrec - fnrec);
  }

  *dfxplus = whatsign(temp2);

  if (!*dfxplus) *dfx = 0.0;
  else *dfx = temp1 + log(fabs(temp2));

  return(*dfx);

} /* rec_locus_reclderiv */


/************************************************************************
 * NRrec_thetalderiv returns the ratio of the 1st to the 2nd derivative *
 * in theta.                                                            *
 *                                                                      *
 * Eqn: kk = (#active lineages) * (#active lineages - 1)                *
 *       c = #coalescences present in a given tree                      *
 *       s = active site count over all active lineages                 *
 *       t = length of interval                                         *
 *       L = Prior Likelihood under theta and rec                       *
 *    SQ() = square of the arguments given                              *
 *                                                                      *
 * dstuff = Sum_over_all_intervals[t * kk/(theta*theta)]                *
 *                                                                      *
 * 1st derivative = sum_over_all_trees[L * (dstuff - (c/theta))]        *
 *                                                                      *
 * 2nd derivative = sum_over_all_trees of                               *
 *    L * [SQ(dstuff - (c/theta)) + (c/SQ(theta)) - (2/theta)*dstuff]   *
 *                                                                      */
double NRrec_thetalderiv(option_struct *op, double theta, double rec,
			 double **llike0, long lowcus, long chain, double *dfn, double *ddfn,
			 long *dfnplus, long *ddfnplus)
{
  long trii, refchain, chaintype, numc;
  double ds, ds2, fn, lth, lrec;
  treerec *tr;

  refchain = REF_CHAIN(chain);
  chaintype = TYPE_CHAIN(chain);
  lth = log(2.0/theta);
  lrec = log(rec);

  *dfn = 0.0;
  *ddfn = 0.0;

  for(trii = 0; trii < op->numout[chaintype]; trii++) {
    tr = (&sum[lowcus][refchain][trii]);
    if (tr->numrecombs && !rec) continue;
    ds = dstuff(theta, tr);
    numc = tr->numcoals;
    ds2 = ds - numc/theta;
    fn = rectreellike(op,theta,rec,lth,lrec,lowcus,chain,trii) -
      llike0[lowcus][trii];

    (*dfn) += exp(fn) * ds2;
    (*ddfn) += exp(fn) * (ds2 * ds2 + numc/(theta*theta) -
			  (2.0/theta)*ds);
  }

  *dfnplus = whatsign(*dfn);
  *ddfnplus = whatsign(*ddfn);

  if (*dfn == 0.0) *dfn = epsilon;
  if (*ddfn == 0.0) *ddfn = epsilon;

  *dfn = log(fabs(*dfn)) - log(op->numout[chaintype]);
  *ddfn = log(fabs(*ddfn)) - log(op->numout[chaintype]);

  return((*dfn) - (*ddfn));

} /* NRrec_thetalderiv */


/**********************************************************************
 * NRrec_reclderiv returns the ratio of the 1st to the 2nd derivative *
 * in recombination rate.                                             *
 *                                                                    *
 * Eqn: kk = (#active lineages) * (#active lineages - 1)              *
 *       r = #recombinations present in a given tree                  *
 *       s = active site count over all active lineages               *
 *       t = length of interval                                       *
 *       L = Prior Likelihood under theta and rec                     *
 *    SQ() = square of the arguments given                            *
 *                                                                    *
 * dstuff = Sum_over_all_intervals[-t * s]                            *
 *                                                                    *
 * 1st derivative = sum_over_all_trees[L * (dstuff + (r/rec))]        *
 *                                                                    *
 * 2nd derivative = sum_over_all_trees of                             *
 *    L * [ SQ(dstuff + (r/rec)) - (r/SQ(rec))]                       *
 *                                                                    */
double NRrec_reclderiv(option_struct *op, double theta, double rec,
		       double **llike0, long lowcus, long chain, double *dfn, double *ddfn,
		       long *dfnplus, long *ddfnplus)
{
  long trii, refchain, chaintype, numr;
  double ds, ds2, fn, lth, lrec;
  treerec *tr;

  refchain = REF_CHAIN(chain);
  chaintype = TYPE_CHAIN(chain);
  lth = log(2.0/theta);
  lrec = log(rec);

  *dfn = 0.0;
  *ddfn = 0.0;

  for(trii = 0; trii < op->numout[chaintype]; trii++) {
    tr = (&sum[lowcus][refchain][trii]);
    if (tr->numrecombs && !rec) continue;
    numr = tr->numrecombs;
    ds = dstuff(0.0, tr);
    ds2 = ds + numr/rec;
    fn = rectreellike(op,theta,rec,lth,lrec,lowcus,chain,trii) -
      llike0[lowcus][trii];

    (*dfn) += exp(fn) * ds2;
    (*ddfn) += exp(fn) * (ds2 * ds2 - numr/(rec*rec));
  }

  *dfnplus = whatsign(*dfn);
  *ddfnplus = whatsign(*ddfn);

  if (*dfn == 0.0) *dfn = epsilon;
  if (*ddfn == 0.0) *ddfn = epsilon;

  *dfn = log(fabs(*dfn)) - log(op->numout[chaintype]);
  *ddfn = log(fabs(*ddfn)) - log(op->numout[chaintype]);

  return((*dfn) - (*ddfn));

} /* NRrec_reclderiv */


/*****************************************************************
 * NRrec_partiallderiv returns the mixed derivative in theta and *
 * recombination rate.                                           *
 *                                                               *
 * Eqn: kk = (#active lineages) * (#active lineages - 1)         *
 *       r = #recombinations present in a given tree             *
 *       c = #coalescences present in a given tree               *
 *       s = active site count over all active lineages          *
 *       t = length of interval                                  *
 *       L = Prior Likelihood under theta and rec                *
 *    SQ() = square of the arguments given                       *
 *                                                               *
 * dtheta = Sum_over_all_intervals[t * kk/(theta*theta)]         *
 * drec = Sum_over_all_intervals[-t * s]                         *
 *                                                               *
 * mixed derivative = sum_over_all_trees of                      *
 *    L * [(drec + (r/rec) * (dtheta - c/theta)]                 *
 *                                                               */
double NRrec_partiallderiv(option_struct *op, double theta,
			   double rec, double **llike0, long lowcus, long chain,
			   double *ddfn, long *ddfnplus)
{
  long trii, refchain, chaintype, numr, numc;
  double dtheta, drec, fn, lth, lrec;
  treerec *tr;

  refchain = REF_CHAIN(chain);
  chaintype = TYPE_CHAIN(chain);
  lth = log(2.0/theta);
  lrec = log(rec);

  *ddfn = 0.0;

  for(trii = 0; trii < op->numout[chaintype]; trii++) {
    tr = (&sum[lowcus][refchain][trii]);
    if (tr->numrecombs && !rec) continue;
    numr = tr->numrecombs;
    numc = tr->numcoals;
    drec = dstuff(0.0, tr);
    dtheta = dstuff(theta, tr);
    fn = rectreellike(op,theta,rec,lth,lrec,lowcus,chain,trii) -
      llike0[lowcus][trii];

    (*ddfn) += exp(fn) * (drec + numr/rec) * (dtheta - numc/theta);
  }

  *ddfnplus = whatsign(*ddfn);

  if (*ddfn == 0.0) *ddfn = epsilon;

  *ddfn = log(fabs(*ddfn)) - log(op->numout[chaintype]);

  return(*ddfn);

} /* NRrec_partiallderiv */

/****************************************************************
 * NRrec_locus_thetalderiv computes the multi-locus 1st and 2nd *
 * derivatives wrt. theta                                       */
double NRrec_locus_thetalderiv(option_struct *op, data_fmt *data,
			       double theta, double rec, double **llike0, double *fx,
			       double *dfx, long *fxplus, long *dfxplus)
{
  long lowcus, lfxplus, ldfxplus, lastchain, numloci;
  double rtheta, fntheta, dfntheta, ddfntheta, dummy, temp, temp1, 
    temp2, temp3, temp4;

  lastchain = totchains - 1;
  numloci = getdata_numloci(op,data);
  temp1 = 0.0;
  temp2 = 0.0;
  temp3 = 0.0;
  temp4 = 0.0;

  for (lowcus = 0; lowcus < numloci; lowcus++) {
    rtheta = theta * op->ne_ratio[lowcus];
    fntheta = fmodel_likelihood(op,data,rtheta,rec,llike0,lowcus,lastchain);
    dummy = NRrec_thetalderiv(op,rtheta,rec,llike0,lowcus,lastchain,&dfntheta,
			      &ddfntheta,&lfxplus,&ldfxplus);

    temp1 += fntheta;
    if (lfxplus) {
      temp2 += op->ne_ratio[lowcus] * lfxplus * exp(dfntheta - fntheta);
      temp4 -= op->ne_ratio[lowcus] * op->ne_ratio[lowcus] * 
	exp(2 * (dfntheta - fntheta));
    }
    if (ldfxplus) temp3 += op->ne_ratio[lowcus] * 
		    ldfxplus * exp(ddfntheta - fntheta);
  }

  temp = (temp2*temp2+temp3+temp4);

  *fxplus = whatsign(temp2);
  *dfxplus = whatsign(temp);

  if (!*fxplus) *fx = 0.0;
  else *fx = temp1 + log(fabs(temp2));
  if (!*dfxplus) *dfx = 0.0;
  else *dfx = temp1 + log(fabs(temp));

  return(*fx - *dfx);

} /* NRrec_locus_thetalderiv */

double NRrec_locus_reclderiv(option_struct *op, data_fmt *data, double theta,
			     double rec, double **llike0, double *fx, double *dfx, long *fxplus,
			     long *dfxplus)
{
  long lowcus, lfxplus, ldfxplus, lastchain, numloci;
  double rtheta, fnrec, dfnrec, ddfnrec, dummy, temp, temp1, temp2, 
    temp3, temp4;

  lastchain = totchains - 1;
  numloci = getdata_numloci(op,data);
  temp1 = 0.0;
  temp2 = 0.0;
  temp3 = 0.0;
  temp4 = 0.0;

  for (lowcus = 0; lowcus < numloci; lowcus++) {
    rtheta = theta * op->ne_ratio[lowcus];
    fnrec = fmodel_likelihood(op,data,rtheta,rec,llike0,lowcus,lastchain);
    dummy = NRrec_reclderiv(op,rtheta,rec,llike0,lowcus,lastchain,&dfnrec,
			    &ddfnrec,&lfxplus,&ldfxplus);

    temp1 += fnrec;
    if (lfxplus) { 
      temp2 += lfxplus * exp(dfnrec - fnrec); 
      temp4 -= exp(2 * (dfnrec - fnrec));
    }
    if (ldfxplus) temp3 += ldfxplus * exp(ddfnrec - fnrec);
  }

  temp = (temp2*temp2+temp3+temp4);

  *fxplus = whatsign(temp2);
  *dfxplus = whatsign(temp);

  if (!*fxplus) *fx = 0.0;
  else *fx = temp1 + log(fabs(temp2));
  if (!*dfxplus) *dfx = 0.0;
  else *dfx = temp1 + log(fabs(temp));

  return(*fx - *dfx);

} /* NRrec_locus_reclderiv */

double NRrec_locus_partiallderiv(option_struct *op, data_fmt *data,
				 double theta, double rec, double **llike0, double *df, long *dfplus)
{
  long lowcus, recplus, thetaplus, pplus, lastchain, ldummy, numloci;
  double rtheta, fn, fnrec, fntheta, fnp, dummy, temp1, temp2, temp3,
    temp4, temp5;

  lastchain = totchains - 1;
  numloci = getdata_numloci(op,data);
  temp1 = 0.0;
  temp2 = 0.0;
  temp3 = 0.0;
  temp4 = 0.0;
  temp5 = 0.0;

  for (lowcus = 0; lowcus < numloci; lowcus++) {
    rtheta = theta * op->ne_ratio[lowcus];
    fn = fmodel_likelihood(op,data,rtheta,rec,llike0,lowcus,lastchain);
    dummy = NRrec_reclderiv(op,rtheta,rec,llike0,lowcus,lastchain,&fnrec,
			    &dummy,&recplus,&ldummy);
    dummy = NRrec_thetalderiv(op,rtheta,rec,llike0,lowcus,lastchain,&fntheta,
			      &dummy,&thetaplus,&ldummy);
    dummy = NRrec_partiallderiv(op,rtheta,rec,llike0,lowcus,lastchain,&fnp,&pplus);

    temp1 += fn;
    temp2 += recplus * exp(fnrec - fn);
    temp3 += op->ne_ratio[lowcus] * thetaplus * exp(fntheta - fn);
    temp4 += op->ne_ratio[lowcus] * pplus * exp(fnp - fn);
    temp5 -= op->ne_ratio[lowcus] * recplus * thetaplus *
      exp(fnrec+fntheta - 2*fn);
  }

  temp2 = temp2*temp3+temp4+temp5;

  *dfplus = whatsign(temp2);
  if (!*dfplus) *df = 0.0;
  else *df = temp1 + log(fabs(temp2));

  return(*df);
} /* NRrec_locus_partiallderiv */

#if 0
/**************************************************************
 * init_matrix fills the derivative matrices for matrixpoint  *
 *                                                            *
 * 1st derivative vector looks like:                          *
 *           0                             1                  *
 *----------------------------------------------------------- *
 *|1st-derivative-wrt-recombination 1st-derivative-wrt-theta| *
 *----------------------------------------------------------- *
 *                                                            *
 * 2nd derivative matrix looks like:                          *
 *             0                             1                *
 *   -------------------------------------------------------- *
 * 0 |2nd-derivative-wrt-recombination  mixed-derivative    | *
 * 1 |   mixed-derivative          2nd-derivative-wrt-theta | *
 *   -------------------------------------------------------- */
void init_matrix(option_struct *op, data_fmt *data, double *fx, double **dfx,
		 double rec, double theta, double **llike0, long *fxplus,
		 long **dfxplus, long lowcus, long chain, boolean thetaonly)
{
  double dummy;

  if (lowcus == -1) {
    dummy = rec_locus_thetalderiv(op,data,theta,rec,llike0,&fx[1],&dfx[1][1],
				  &fxplus[1],&dfxplus[1][1]);
    if (!thetaonly) {
      dummy = rec_locus_reclderiv(op,data,theta,rec,llike0,&fx[0],&dfx[0][0],
				  &fxplus[0],&dfxplus[0][0]);
      dummy = rec_locus_partiallderiv(op,data,theta,rec,llike0,&dfx[0][1],
				      &dfxplus[0][1]);
    }
  } else {
    dummy = rec_thetalderiv(op,theta,rec,llike0,lowcus,chain,&fx[1],&dfx[1][1],
			    &fxplus[1],&dfxplus[1][1]);
    if (!thetaonly) {
      dummy = rec_reclderiv(op,theta,rec,llike0,lowcus,chain,&fx[0],&dfx[0][0],
			    &fxplus[0],&dfxplus[0][0]);
      dummy = rec_partiallderiv(op,theta,rec,llike0,lowcus,chain,&dfx[0][1],
				&dfxplus[0][1]);
    }
  }

  dfx[1][0] = dfx[0][1];
  dfxplus[1][0] = dfxplus[0][1];

} /* init_matrix */

/****************************************************************
 * to_LNparameters performs the change of variables from normal *
 * parameters to LN(parameters).   Remember that all params are *
 * currently stored as LNs with sign bits!                      *
 * WARNING: this code assumes that legal parameter values can   *
 * only be positive numbers                                     */
void to_LNparameters(double theta, double rec, double *fx, double **dfx,
		     long *fxplus, long **dfxplus)
{

  /* first change the 1st derivatives */
  fx[0] += log(rec);
  fx[1] += log(theta);

  /* now the 2nd derivatives */
  dfx[0][0] = fxplus[0] * exp(fx[0]) + 
    dfxplus[0][0] * exp(dfx[0][0]) * rec * rec;
  dfxplus[0][0] = whatsign(dfx[0][0]);
  dfx[0][0] = log(fabs(dfx[0][0]));
  dfx[1][1] = fxplus[1] * exp(fx[1]) + 
    dfxplus[1][1] * exp(dfx[1][1]) * theta * theta;
  dfxplus[1][1] = whatsign(dfx[1][1]);
  dfx[1][1] = log(fabs(dfx[1][1]));

  /* now the mixed derivatives */
  dfx[0][1] += log(rec) + log(theta);

} /* to_LNparameters */

/****************************************************************
 * matrix_denom returns TRUE if it succeeds in its calculations *
 * and FALSE otherwise.  Ln(result of the calculation) is       *
 * returned in denom with the appropiate sign in denomsign.     *
 *                                                              *
 * EQN: ddtheta = 2nd derivative in theta                       *
 *      ddrec = 2nd derivative in recombination rate            *
 *      ddreta = mixed derivative in parameters                 *
 *                                                              *
 *      denom = ddrec*ddtheta - ddreta*ddreta                   *
 *                                                              */
boolean matrix_denom(double **dfx, long **dfxplus, double *denom,
		     long *denomsign)
{
  double temp1, temp2, max;

  temp1 = dfx[0][0] + dfx[1][1];

  if (dfxplus[0][1] && dfxplus[1][0]) {
    temp2 = dfx[0][1] + dfx[1][0];
    if (temp1 > temp2) max = temp1;
    else max = temp2;
    *denom = 0.0;
    if (temp1 - max > EXPMIN)
      *denom += dfxplus[0][0] * dfxplus[1][1] * exp(temp1 - max);
    if (temp2 - max > EXPMIN)
      *denom -= dfxplus[0][1] * dfxplus[1][0] * exp(temp2 - max);
    *denomsign = whatsign(*denom);
    *denom = log(fabs(*denom)) + max;
  } else {
    *denomsign = dfxplus[0][0] * dfxplus[1][1];
    *denom = temp1;
  }

  /* check that dfx is invertable */
  if(!*denom) return(FALSE);
  else return(TRUE);

} /* matrix_denom */

/****************************************************************
 * check_curvature returns TRUE if we are within the attraction *
 * zone of a maximum, otherwise it returns FALSE.               *
 *                                                              *
 * EQN: dtheta = 1st derivative in theta                        *
 *      ddtheta = 2nd derivative in theta                       *
 *      drec = 1st derivative in recombination rate             *
 *      ddrec = 2nd derivative in recombination rate            *
 *      ddreta = mixed derivative in parameters                 *
 *      SQ() = square of the arguments                          *
 *                                                              *
 * numerator =                                                  *
 *   ddtheta*SQ(drec) - 2*ddreta*dtheta*drec + ddrec*SQ(dtheta) *
 *                                                              *
 * result = numerator / (ddrec*ddtheta - SQ(ddreta))            *
 *                                                              *
 * if result is negative then we're within the maximum's zone.  *
 *                                                              */
boolean check_curvature(double *fx, double **dfx, long *fxplus,
			long **dfxplus)
{
  double denom, max, temp1, temp2, temp3, temp4, result;
  long denomsign;

  /* initializations for lint */
  temp3 = temp4 = 0.0;

  if (!matrix_denom(dfx,dfxplus,&denom,&denomsign)) {
    fprintf(ERRFILE,"WARNING--infinite curvature encountered!");
    return(FALSE);
  }

  max = NEGMAX;

  temp1 = dfx[1][1] + 2 * fx[0] - denom;
  if (temp1 > max) max = temp1;
  temp2 = dfx[0][0] + 2 * fx[1] - denom;
  if (temp2 > max) max = temp2;
  if (dfxplus[1][0]) {
    temp3 = dfx[1][0] + fx[0] + fx[1] - denom;
    if (temp3 > max) max = temp3;
  }
  if (dfxplus[0][1]) {
    temp4 = dfx[0][1] + fx[0] + fx[1] - denom;
    if (temp4 > max) max = temp4;
  }

  result = 0.0;

  if (temp1 - max > EXPMIN)
    /* the sign of the 1st derivative is irrelevant since it is squared */
    result += dfxplus[1][1] * denomsign * exp(temp1 - max);
  if (temp2 - max > EXPMIN)
    /* the sign of the 1st derivative is irrelevant since it is squared */
    result += dfxplus[0][0] * denomsign * exp(temp2 - max);
  if (dfxplus[1][0] && (temp3 - max > EXPMIN))
    result -= dfxplus[1][0] * fxplus[0] * fxplus[1] * denomsign *
      exp(temp3 - max);
  if (dfxplus[0][1] && (temp4 - max > EXPMIN))
    result -= dfxplus[0][1] * fxplus[0] * fxplus[1] * denomsign *
      exp(temp4 - max);

  /* note that we have actually calculated our desired value divided
     by Exp(max), but this should not affect the sign */
  if (result < 0) return(TRUE);
  else return(FALSE);

} /* check_curvature */

/***************************************************************
 * calc_change calculates the vector used in Newton-Raphson to *
 * adjust the parameters being maximized.  The vector is       *
 * returned in answ.                                           *
 *                                                             *
 * EQN: dtheta = 1st derivative in theta                       *
 *      ddtheta = 2nd derivative in theta                      *
 *      drec = 1st derivative in recombination rate            *
 *      ddrec = 2nd derivative in recombination rate           *
 *      ddreta = mixed derivative in parameters                *
 *      SQ() = square of the arguments                         *
 *                                                             *
 * denom = ddtheta*ddrec - SQ(ddreta)                          *
 * rec change = (ddtheta*drec - ddreta*dtheta) / denom         *
 * theta change = (ddrec*dtheta - ddreta*drec) / denom         *
 *                                                             */
boolean calc_change(double *fx, double **dfx, long *fxplus, 
		    long **dfxplus, double *answ)
{
  double denom, temp1, temp2;
  long denomsign;

  if (!matrix_denom(dfx,dfxplus,&denom,&denomsign)) return(FALSE);

  temp1 = dfxplus[1][1] * fxplus[0] * denomsign *
    exp(dfx[1][1] + fx[0] - denom);
  if (dfxplus[0][1]) {
    temp2 = dfxplus[0][1] * fxplus[1] * denomsign *
      exp(dfx[0][1] + fx[1] - denom);
    answ[0] = temp1 - temp2;
  } else answ[0] = temp1;

  temp1 = dfxplus[0][0] * fxplus[1] * denomsign *
    exp(dfx[0][0] + fx[1] - denom);
  if (dfxplus[1][0]) {
    temp2 = dfxplus[1][0] * fxplus[0] * denomsign *
      exp(dfx[1][0] + fx[0] - denom);
    answ[1] = temp1 - temp2;
  } else answ[1] = temp1;

  return(TRUE);

} /* calc_change */

#define HALF_MAX 100 
/******************************************************************
 * halfback attempts to find a maximum by halving one or more of  *
 * its parameters until the likelihood increases.                 *
 *                                                                *
 * pass vflag = 1 for halving on rec, vflag = 2 for theta         *
 * and vflag = 0 for both                                         *
 *                                                                *
 * halfback returns TRUE if it succeeds in finding a higher value.*/
boolean halfback(option_struct *op, data_fmt *data, double theta, double rec,
		 double thetachange, double recchange, double **llike0, double oldlike,
		 double *newtheta, double *newrec, long vflag, long lowcus, long chain,
		 boolean thetaonly)
{
  long numloop;
  double newlike;
  boolean succeeded;

  /* don't do any halving back in "rec" if thetaonly optimization is
     being done */
  if ((vflag == 1 || vflag == 0) && thetaonly) return(FALSE);

  numloop = 0;
  succeeded = FALSE;
  while(!succeeded) {
    numloop++;
    if (vflag == 1) recchange /= 2.0;
    if (vflag == 2) thetachange /= 2.0;
    if (vflag == 0) {recchange /= 2.0; thetachange /= 2.0;}
#if 1
    if (!thetaonly) *newrec = makepositive(rec,recchange);
    *newtheta = makepositive(theta,thetachange);
#endif
#if 0
    if (!thetaonly) *newrec = changeparamLN(rec,recchange);
    *newtheta = changeparamLN(theta,thetachange);
#endif
    newlike = fmodel_likelihood(op,data,*newtheta,*newrec,llike0,lowcus,chain);
    if(newlike >= oldlike) succeeded = TRUE;
    if(numloop > HALF_MAX) break;
  }

  return(succeeded);

} /* halfback */


/***********************************************************************
 * non_nr_step sets the change in the parameters if a non-NR algorithm *
 * step is indicated.                                                  */
void non_nr_step(option_struct *op, data_fmt *data, double theta, double rec,
		 double **llike0, long lowcus, long chain, double *recchange,
		 double *thetachange, double *pchange, long *fxplus)
{
  double tempth, temprec, templike;

  *recchange = -1 * fxplus[0] * rec/2.0;
  *thetachange = -1 * fxplus[1] * theta/2.0;
  temprec = makepositive(rec,*recchange);
  tempth = makepositive(theta,*thetachange);
  templike = fmodel_likelihood(op,data,tempth,temprec,llike0,lowcus,chain);

  temprec = makepositive(rec,-1.0 * pchange[0]);
  tempth = makepositive(theta,-1.0 * pchange[1]);
  if (templike < fmodel_likelihood(op,data,tempth,temprec,llike0,lowcus,chain)) {
    *recchange = -1.0 * pchange[0];
    *thetachange = -1.0 * pchange[1];
    templike = llike;
  }

  temprec = makepositive(rec, fxplus[0] * rec/2.0);
  tempth = makepositive(theta, fxplus[1] * theta/2.0);
  llike = fmodel_likelihood(op,data,tempth,temprec,llike0,lowcus,chain);
  if (templike < llike) {
    *recchange = fxplus[0]*rec/2.0;
    *thetachange = fxplus[1]*theta/2.0;
  }

  temprec = makepositive(rec, pchange[0]/2.0);
  tempth = makepositive(theta, pchange[1]/2.0);
  llike = fmodel_likelihood(op,data,tempth,temprec,llike0,lowcus,chain);
  if (templike < llike) {
    *recchange = pchange[0];
    *thetachange = pchange[1];
  }

} /* non_nr_step */
#endif


/******************************************************************
 * found_recomb() returns TRUE if the passed chain sampled a tree *
 * containing any recombinations, FALSE otherwise.                */
boolean found_recomb(option_struct *op, data_fmt *data, long lowcus, long chain)
{
  long i, j, chaintype, refchain, startlowcus, endlowcus;

  chaintype = TYPE_CHAIN(chain);
  refchain = REF_CHAIN(chain);

  if (lowcus == -1) {startlowcus = 0; endlowcus = getdata_numloci(op,data);}
  else {startlowcus = lowcus; endlowcus = lowcus + 1;}

  for (i = startlowcus; i < endlowcus; i++)
    for (j = 0; j < op->numout[chaintype]; j++) 
      if (sum[i][refchain][j].numrecombs) return(TRUE);

  return(FALSE);
} /* found_recomb */


#if 0
/*******************************************************************
 * matrixpoint_alloc allocates memory for the local arrays used in *
 * matrixpoint.                                                    */
void matrixpoint_alloc(option_struct *op, data_fmt *data, long chain, double **fx,
		       long **fxplus, double ***dfx, long ***dfxplus, double **fchange,
		       long **fplus, double ***llike0)
{
  long i, chaintype, numloci;

  chaintype = TYPE_CHAIN(chain);
  numloci = getdata_numloci(op,data);

  (*fxplus) = (long *)calloc(1,NUMPARAMETERS * sizeof(long));
  (*dfxplus) = (long **)calloc(1,NUMPARAMETERS * sizeof(long *));
  (*fplus) = (long *)calloc(1,NUMPARAMETERS * sizeof(long));
  (*dfxplus)[0] = (long *)calloc(1,NUMPARAMETERS*NUMPARAMETERS * sizeof(long));
  (*fx) = (double *)calloc(1,NUMPARAMETERS * sizeof(double));
  (*fchange) = (double *)calloc(1,NUMPARAMETERS * sizeof(double));
  (*dfx) = (double **)calloc(1,NUMPARAMETERS * sizeof(double *));
  (*dfx)[0] = (double *)calloc(1,NUMPARAMETERS*NUMPARAMETERS * sizeof(double));
  for(i = 1; i < NUMPARAMETERS; i++) {
    (*dfxplus)[i] = (*dfxplus)[0] + i * NUMPARAMETERS;
    (*dfx)[i] = (*dfx)[0] + i * NUMPARAMETERS;
  }
  (*llike0) = (double **)calloc(1,numloci*sizeof(double *));
  (*llike0)[0] = (double *)calloc(1,numloci*op->numout[chaintype]*sizeof(double));
  for(i = 1; i < numloci; i++) (*llike0)[i] = (*llike0)[0] + i*op->numout[chaintype];

} /* matrixpoint_alloc */
#endif


/**********************************************************************
 * precalc_llike0 calculates a common denominator used in much of the *
 * prior likelihood and derivative calculations.  This is done for    *
 * speed reasons.                                                     */
void precalc_llike0(option_struct *op, data_fmt *data, long lowcus, long chain, double **llike0)
{
  long i, j, chaintype, numloci;
  double th0, r0, lth0, lrec0;

  chaintype = TYPE_CHAIN(chain);
  numloci = getdata_numloci(op,data);

  if (lowcus == -1)
    for(i = 0; i < numloci; i++) {
      th0 = theti[i][chain];
      r0 = reci[i][chain];
      lth0 = log(2.0/th0);
      lrec0 = log(r0);
      for(j = 0; j < op->numout[chaintype]; j++)
	llike0[i][j] = rectreellike(op,th0,r0,lth0,lrec0,i,chain,j);
    }
  else {
    th0 = theti[lowcus][chain];
    r0 = reci[lowcus][chain];
    lth0 = log(2.0/th0);
    lrec0 = log(r0);
    for (j = 0; j < op->numout[chaintype]; j++)
      llike0[lowcus][j] = rectreellike(op,th0,r0,lth0,lrec0,lowcus,chain,j); 
  }

} /* precalc_llike0 */


#if 0
/**********************************************************************
 * matrixpoint runs the Newton-Raphson maximizer for the parameters.  *
 * Maximized answers are returned in thetaresult and recresult.  Pass *
 * lowcus = -1 for multi-locus estimation.                            */
void matrixpoint(option_struct *op, data_fmt *data, long lowcus, long chain,
		 double *thetaresult, double *recresult)
{
  long i, report, samecount, *fxplus, **dfxplus, *pplus;
  double theta, newtheta, rec, newrec, oldlike, newlike, templike,
    thetachange, recchange, *pchange, *fmat, **dfmat, ltolerance, **llike0;
  boolean succeeded, thetaonly, firstreport;

  firstreport = TRUE;

  matrixpoint_alloc(op,data,chain,&fmat,&fxplus,&dfmat,&dfxplus,
		    &pchange,&pplus,&llike0);

  theta = *thetaresult;
  if (found_recomb(op,data,lowcus,chain)) {
    thetaonly = FALSE;
    rec = *recresult;
    if (rec == 0.0) rec = EPSILON;
  }
  else {thetaonly = TRUE; rec = EPSILON;}

  /* speedups involving pre-calculating the LnL wrt. th0 and rec0 */
  precalc_llike0(op,data,lowcus,chain,llike0);

  /* point estimate of growth & theta simultaneously */
  oldlike = fmodel_likelihood(op,data,theta,rec,llike0,lowcus,chain);
  i = 0;
  report = 1;
  samecount = 0;
  ltolerance = log(epsilon1);

  /* solve by modified Newton Raphson */
  while (1) {

    if (!thetaonly) {
      init_matrix(op,data,fmat,dfmat,rec,theta,llike0,fxplus,
		  dfxplus,lowcus,chain,thetaonly);

      /* if within epsilon of the maximum, we're done! */
      if (fmat[0] < ltolerance && fmat[1] < ltolerance) {
	if (!firstreport) printf(" steps\n");
	break;
      }

      while (!dfxplus[0][0] || !dfxplus[1][1]) {
	if (!dfxplus[0][0]) rec += epsilon;
	if (!dfxplus[1][1]) theta += epsilon;
	init_matrix(op,data,fmat,dfmat,rec,theta,llike0,fxplus,
		    dfxplus,lowcus,chain,thetaonly);
      }
    } else {
      if (lowcus == -1)
	thetachange = rec_locus_thetalderiv(op,data,theta,rec,
					    llike0,&fmat[1],&dfmat[1][1],&fxplus[1],&dfxplus[1][1]);
      else 
	thetachange = rec_thetalderiv(op,theta,rec,llike0,
				      lowcus,chain,&fmat[1],&dfmat[1][1],&fxplus[1],&dfxplus[1][1]);

      /* just check theta for within epsilon of maximum */
      if (fmat[1] < ltolerance) break;

      /* we're only dealing with theta changes, so our Newton proposal
	 involves *adding* the 1st/abs(2nd); still need to initialize
	 the other change though */
      recchange = 0.0;
      newrec = EPSILON;
      thetachange = -1 * fxplus[1] * exp(thetachange);
    }

#if 0
    /* change to LN(parameters) to make NR behave better */
    to_LNparameters(theta,rec,fmat,dfmat,fxplus,dfxplus);
#endif

    if (!thetaonly) {
      while(!calc_change(fmat,dfmat,fxplus,dfxplus,pchange)) {
	fprintf(ERRFILE,"\nERROR:OH NO, no inversion possible!!!!!\n");
	rec += epsilon;
	theta += epsilon;
	init_matrix(op,data,fmat,dfmat,rec,theta,llike0,fxplus,
		    dfxplus,lowcus,chain,thetaonly);
      }
      
      if (check_curvature(fmat,dfmat,fxplus,dfxplus)) {
	recchange = pchange[0];
	thetachange = pchange[1];
      } else {
	non_nr_step(op,data,theta,rec,llike0,lowcus,chain,
		    &recchange,&thetachange,pchange,fxplus);
      }
    }

#if 1
    if (!thetaonly) newrec = makepositive(rec,recchange);
    newtheta = makepositive(theta,thetachange);
#endif
#if 0
    if (!thetaonly) newrec = changeparamLN(rec,recchange);
    newtheta = changeparamLN(theta,thetachange);
#endif

    newlike = fmodel_likelihood(op,data,newtheta,newrec,llike0,lowcus,chain);

    templike = newlike;
    succeeded = TRUE;
    while (succeeded) {
      succeeded = FALSE;
      recchange *= 2.0;
      thetachange *= 2.0;
#if 1
      if (!thetaonly) newrec = makepositive(rec,recchange);
      newtheta = makepositive(theta,thetachange);
#endif
#if 0
      if (!thetaonly) newrec = changeparamLN(rec,recchange);
      newtheta = changeparamLN(theta,thetachange);
#endif
      newlike = 
	fmodel_likelihood(op,data,newtheta,newrec,llike0,lowcus,chain);
      if (newlike > templike) {
	templike = newlike;
	succeeded = TRUE;
	continue;
      }
      /* else */
      recchange /= 2.0;
      thetachange /= 2.0;
#if 1
      if (!thetaonly) newrec = makepositive(rec,recchange);
      newtheta = makepositive(theta,thetachange);
#endif
#if 0
      if (!thetaonly) newrec = changeparamLN(rec,recchange);
      newtheta = changeparamLN(theta,thetachange);
#endif
      newlike = 
	fmodel_likelihood(op,data,newtheta,newrec,llike0,lowcus,chain);
    }

    if(newlike < oldlike) {
      /* in case we overshoot the maximum, don't jump so far...*/
      if(!halfback(op,data,theta,rec,thetachange,recchange,llike0,oldlike,
		   &newtheta,&newrec,1L,lowcus,chain,thetaonly))
	if(!halfback(op,data,theta,rec,thetachange,recchange,llike0,oldlike,
		     &newtheta,&newrec,2L,lowcus,chain,thetaonly))
	  if(!halfback(op,data,theta,rec,thetachange,recchange,llike0,oldlike,
		       &newtheta,&newrec,0L,lowcus,chain,thetaonly))
	    /* we expect to come here if we are doing theta-only optimization */
	    if (!thetaonly)
	      fprintf(ERRFILE,"ERROR:failure in matrixpoint halving!!\n");
      newlike = fmodel_likelihood(op,data,newtheta,newrec,llike0,lowcus,chain);
    }

    if (newrec < epsilon) {
      newrec = EPSILON;
      thetaonly = TRUE;
    }

    if (newlike == oldlike) {
      samecount++;
      if (samecount == REPEAT_TOLERANCE) {
#if 0
	if (!firstreport) printf(" steps\n");
	fprintf(ERRFILE,"claimed to find something flat in matrixpoint\n");
#endif
	break;
      }
    } else samecount = 0;

    oldlike = newlike;
    if (!thetaonly) rec = newrec;
    theta = newtheta;
    i++;
#if 0
    if (i == report) {
      if (!firstreport) for(j = 0; j < STEPPRINT; j++) printf("\b");
      else printf("made ");
      printf("%*ld",(int)STEPPRINT,i);
      report += 1;
      firstreport = FALSE;
    }
#endif
    if (i > STEPMAX) {
#if 0
      if (!firstreport) printf(" steps\n");
#endif
      fprintf(ERRFILE,"ERROR:Matrixpoint failure, wah!!!!\n");
      break;
    }
  }

  *recresult = rec;
  *thetaresult = theta;

  free(pchange);
  free(pplus);
  free(fmat);
  free(fxplus);
  free(dfmat[0]);
  free(dfmat);
  free(dfxplus[0]);
  free(dfxplus);
  free(llike0[0]);
  free(llike0);

} /* matrixpoint */
#endif


#define NUMPARAM 2
#define STARTNORM 1e20
#define BADLIKE  -100000000.0
#define LOGPARAM 1

/**************************************************************
 * reset_hess() sets the passed matrix to an Identity matrix. */
void reset_hess (double **hess, long n)
{
  long pop;
  memset (hess[0], 0, sizeof (double) * n * n);
  for (pop = 1; pop < n; pop++)
    {
      hess[pop] = hess[0] + pop * n;
      hess[pop][pop] = 1.0;
    }
  hess[0][0] = 1.0;
} /* reset_hess */


void calc_loci_param(option_struct *op, data_fmt *data, double *oparam,
		     double *nparam, double lambda, double *dv)
{
  long i;

  for(i = 0; i < NUMPARAM; i++) {
#if LOGPARAM
    nparam[i] = oparam[i] * exp(MAX(-100,MIN(-lambda * dv[i],100)));
#endif
#if !LOGPARAM
    nparam[i] = oparam[i] + -lambda * dv[i];
#endif
    if (nparam[i] < epsilon) nparam[i] = epsilon;
  }

} /* calc_loci_param */


/* psi() is a scaled likelihood calculator, using lambda */
double psi(option_struct *op, data_fmt *data, double *param,
	   double **denom, long lowcus, long chain, double lambda, double *dv)
{
  double *nparam, llike;

  nparam = (double *)calloc(NUMPARAM,sizeof(double));

  calc_loci_param(op,data,param,nparam,lambda,dv);
  llike = fmodel_likelihood(op,data,nparam[1],nparam[0],denom,lowcus,chain);
  llike *= -1.0;

  free(nparam);

  return(llike);

} /* psi */


/* Code by Peter Beerli */
/* line searcher
   finds the maximum in a direction
   this should be replaced with something more efficient. */
#define PP 0.61803399
#define QQ 0.38196601
#define MOVE3(a,b,c,d) (a)=(b);(b)=(c);(c)=(d)

/* Code by Peter Beerli, adapted by Jon Yamato */
double calc_line (option_struct *op, data_fmt *data, double *param,
		  double **denom, double *dv, long lowcus, long chain)
{
  /* a < b < c AND psia > psib < psic */
  double d, psib, psic, a = -10, b = 0.1, c = 10;
  d = c;
  if (fabs (c - b) > fabs (b - a))
    {
      c = b + QQ * (c - b);
    }
  else
    {
      c = b;
      b = b - QQ * (b - a);
    }
  psib = psi(op,data,param,denom,lowcus,chain,b,dv);
  psic = psi(op,data,param,denom,lowcus,chain,c,dv);
  while (fabs (d - a) > EPSILON * (fabs (b) + fabs (c)))
    {
      if (psic < psib)
        {
          MOVE3 (a, b, c, PP * b + QQ * d);
          psib = psic;
          psic = psi(op,data,param,denom,lowcus,chain,c,dv);
        }
      else
        {
          MOVE3 (d, c, b, PP * c + QQ * a);
          psic = psib;
          psib = psi(op,data,param,denom,lowcus,chain,b,dv);
        }
    }
  if (psib < psic)
    {
      return b;
    }
  else
    {
      return c;
    }
} /* calc_line */


/********************************************************
 * calcnorm() calculates the norm of the passed vector. *
 * Code by Peter Beerli                                 */
double calcnorm (double *d, long size)
{
  int i;
  double total = 0.;
  for (i = 0; i < (int) size; i++)
    {
      total += d[i] * d[i];
    }
  return sqrt (total);
} /* calcnorm */


/*********************************************************************
 * calc_hessian() calculates the derived 2nd derivative matrix for a *
 * new set of parameters.                                            *
 * Code by Peter Beerli                                              *
 *                                                                   *
 * "delta" is the vector of (new_parameter - old_parameter)          *
 * "gamma" is the vector of (new_1stderiv - old_1stderiv)            */
void calc_hessian (double **hess, long n, double *delta, double *gama)
{
  double **dd, *temp, t;
  long i, j, k;
  double dtg;
  temp = (double *) calloc (n, sizeof (double));
  dd = (double **) calloc (n, sizeof (double *));
  dd[0] = (double *) calloc (n * n, sizeof (double));
  dtg = delta[0] * gama[0];
  for (i = 1; i < n; i++)
    {
      dd[i] = dd[0] + n * i;
      dtg += delta[i] * gama[i];
    }
  if (dtg != 0.0)
    dtg = 1. / dtg;
  else
    {
      reset_hess (hess, n);
      return;
    }
  for (i = 0; i < n; i++)
    {
      for (j = 0; j < n; j++)
        {
          temp[i] += gama[j] * hess[j][i];
        }
    }
  t = 0.0;
  for (i = 0; i < n; i++)
    t += temp[i] * gama[i];
  t = (1.0 + t * dtg) * dtg;
  for (i = 0; i < n; i++)
    {
      for (j = 0; j < n; j++)
        {
          for (k = 0; k < n; k++)
            {
              dd[i][j] += delta[i] * gama[k] * hess[k][j];
            }
        }
    }
  for (i = 0; i < n; i++)
    {
      temp[i] = 0.0;
      for (j = 0; j < n; j++)
        {
          temp[i] += hess[i][j] * gama[j];
        }
    }
  for (i = 0; i < n; i++)
    {
      for (j = 0; j < n; j++)
        {
          dd[i][j] += temp[i] * delta[j];
          dd[i][j] *= dtg;
          hess[i][j] += delta[i] * delta[j] * t - dd[i][j];
        }
    }
  free (temp);
  free (dd[0]);
  free (dd);
} /* calc_hessian */


/***********************************************************************
 * calc_dirv() calculates the new direction vector for use in the BFGS *
 * search.                                                             *
 * Code by Peter Beerli                                                *
 *                                                                     *
 * "gxv" is the vector of 1st derivatives                              *
 * "hess" is the hessian built off of the same 1st derivatives         */
void calc_dirv (double *dv, double **hess, double *gxv, long n)
{
  long i, j;
  memset (dv, 0, sizeof (double) * n);
  for (i = 0; i < n; i++)
    {
      for (j = 0; j < n; j++)
        {
          dv[i] += hess[i][j] * gxv[j];
        }
    }
} /* calc_dirv */


/*******************************************************************
 * grad2loggrad() transforms the 1st derivatives to log(parameter) *
 * 1st derivatives.                                                */
void grad2loggrad (double *param, double *d, double PGC, long nn)
{
  long i;
  for (i = 0; i < nn; i++)
    {
      d[i] = -param[i] * d[i] / PGC;
      /* to log derivatives */
      /* the division by PGC
	 (the underflow protected uncorrected likelihood)
	 is made here, because
	 in multiple loci where I have to sum over loci
	 before "taking logs"
	 I use the derivatives of L with PGC=1 instead of Log(L)
	 where I use PGC = nr->PGC, the minus is for minimizing
	 the function -L
	 instead of maximizing L
      */
    }
} /* grad2loggrad */


/*********************************************************************
 * calc_grad() fills the gradient vector with the derivatives of the *
 * log_likelihood functions.                                         */
void calc_grad(option_struct *op, data_fmt *data, double *param,
	       double **denom, long lowcus, long chain, double *grad)
{
  long i, *gplus;
  double th, rec;
  boolean multilocus;

  gplus = (long *)calloc(NUMPARAM,sizeof(long));

  th = param[1];
  rec = param[0];
  multilocus = ((lowcus == -1) ? TRUE : FALSE);

  /* calculate the simple, non-log_likelihood, derivatives in log form */
  grad[0] = ((multilocus) ? 
	     rec_locus_reclderiv(op,data,th,rec,denom,&(grad[0]),&(gplus[0])) :
	     rec_reclderiv(op,th,rec,denom,lowcus,chain,&(grad[0]),&(gplus[0])));
  grad[1] = ((multilocus) ?
	     rec_locus_thetalderiv(op,data,th,rec,denom,&(grad[1]),&(gplus[1])) :
	     rec_thetalderiv(op,th,rec,denom,lowcus,chain,&(grad[1]),&(gplus[1])));

  /* now divide by likelihood function, or since in logs, subtract
     log_likelihood */
  for(i = 0; i <NUMPARAM; i++)
    grad[i] -= frec_chainllike(op,th,rec,denom,lowcus,chain);

  /* now de-log the deriviatives */
  for(i = 0; i <NUMPARAM; i++)
    grad[i] = exp(grad[i]) * gplus[i];

  /* now transform to log-parameters */
#if LOGPARAM
  grad2loggrad(param,grad,1.0,NUMPARAM);
#endif

  free(gplus);

} /* calc_grad */

#define TRIALPRINT 5

/********************************************************************
 * broydenpoint() runs a broyen-fletcher-goldfarb-shanno search for *
 * the conjoint maxima of the likelihood curve in theta & recrate.  *
 *                                                                  *
 * Original code by Peter Beerli.  Modified for recombine by JAY.   */
void broydenpoint(option_struct *op, data_fmt *data, long lowcus, long chain,
		  double *thetaresult, double *recresult, double **llike0)
{
  long i, numloci, chaintype, trials;
  double *deltap, *deltad, *dirv, *nderiv, *oderiv, *nparam, *oparam,
    **hess, nllike, ollike, lambda, norm, norm20;

  numloci = getdata_numloci(op,data);
  chaintype = TYPE_CHAIN(chain);

  deltap = (double *)calloc(NUMPARAM,sizeof(double));
  deltad = (double *)calloc(NUMPARAM,sizeof(double));
  oparam = (double *)calloc(NUMPARAM,sizeof(double));
  nparam = (double *)calloc(NUMPARAM,sizeof(double));
  nderiv = (double *)calloc(NUMPARAM,sizeof(double));
  oderiv = (double *)calloc(NUMPARAM,sizeof(double));
  dirv = (double *)calloc(NUMPARAM,sizeof(double));

  hess = (double **)calloc(NUMPARAM,sizeof(double *));
  hess[0] = (double *)calloc(NUMPARAM*NUMPARAM,sizeof(double));
  for(i = 1; i < NUMPARAM; i++)
    hess[i] = hess[0] + i*NUMPARAM;

  reset_hess(hess,NUMPARAM);
  oparam[0] = *recresult;
  if (oparam[0] < epsilon) oparam[0] = epsilon;
  oparam[1] = *thetaresult;
  ollike = fmodel_likelihood(op,data,oparam[1],oparam[0],llike0,lowcus,chain);
  calc_grad(op,data,oparam,llike0,lowcus,chain,oderiv);
  memcpy(nderiv,oderiv,NUMPARAM*sizeof(double));
  memcpy(nparam,oparam,NUMPARAM*sizeof(double));
  memcpy(dirv,nderiv,NUMPARAM*sizeof(double));
  norm20 = STARTNORM; 

  for(trials = 0; trials < STEPMAX; trials++) {
#ifdef MAC
    eventloop();
#endif
    nllike = NEGMAX;
    if (ollike < BADLIKE) lambda = 1;
    else lambda = calc_line(op,data,oparam,llike0,dirv,lowcus,chain);
    while(nllike < ollike && fabs(lambda) > EPSILON ) {
      calc_loci_param(op,data,oparam,nparam,lambda,dirv);
      nllike =
	fmodel_likelihood(op,data,nparam[1],nparam[0],llike0,lowcus,chain);
      lambda /= -2.0;
    }

    /* are we done yet? */
    norm = calcnorm(nderiv,NUMPARAM);
    if ((trials+1)%20 == 0) {
      if (fabs(norm-norm20) < EPSILON) break;
      norm20 = norm;
    }
    if (norm < epsilon1) break;
    else if (ollike >= nllike) lambda = 0; /* zero is a flag value */
    else ollike = nllike;

    if (fabs(lambda) <= EPSILON) {
      reset_hess(hess,NUMPARAM);
      memcpy(dirv,nderiv,NUMPARAM*sizeof(double));
      continue;
    }

    /* not done yet */
    memcpy(oderiv,nderiv,NUMPARAM*sizeof(double));
    calc_grad(op,data,nparam,llike0,lowcus,chain,nderiv);
    for(i = 0; i < NUMPARAM; i++) {
#if LOGPARAM
      deltap[i] = log(nparam[i]) - log(oparam[i]);
#endif
#if !LOGPARAM
      deltap[i] = nparam[i] - oparam[i];
#endif
      deltad[i] = nderiv[i] - oderiv[i];
    }
    calc_hessian(hess,NUMPARAM,deltap,deltad);
    calc_dirv(dirv,hess,nderiv,NUMPARAM);
    memcpy(oparam,nparam,NUMPARAM*sizeof(double));

#if 0
    if (op->progress) {
      if (trials)
	for(i = 0; i < TRIALPRINT; i++) printf("\b");
      else printf("maximized curve in: ");
      printf("%*ld",(int)TRIALPRINT,trials+1);
      fflush(stdout);
    }
#endif

  }

#if 0
  if (op->progress)
    printf(" trials");
#endif

  (*recresult) = nparam[0];
  (*thetaresult) = nparam[1];

  free(deltap);
  free(deltad);
  free(oparam);
  free(nparam);
  free(nderiv);
  free(oderiv);
  free(dirv);
  free(hess[0]);
  free(hess);

} /* broydenpoint */


/********************************************************************
 * This function does a simple grid search for a better likelihood, *
 * returning TRUE if it finds one and FALSE otherwise.              */
#define GRIDSIDE 5

boolean gridcheck(option_struct *op, data_fmt *data, double test_theta, double 
		  test_r, long lowcus, long chain)
{
  double **grid, **llike0;
  long i, j, numloci, chaintype;
  double temp_theta, temp_r;

  /* if the user has forced one (or both) variables constant,
     then abort this check */
  if (op->holding != 0) return(FALSE);

  /* if the estimate of rec-rate has been arbitrarily truncated
     by the maximizer code, abort this check */
  if (test_r == EPSILON) return(FALSE);

  /* precalculate parts of the likelihood */
  chaintype = TYPE_CHAIN(chain);
  numloci = getdata_numloci(op,data);

  llike0 = (double **)calloc(1,numloci*sizeof(double *));
  llike0[0] = (double *)
    calloc(1,numloci*op->numout[chaintype]*sizeof(double));
  for(i = 1; i < numloci; i++)
    llike0[i] = llike0[0] + i*op->numout[chaintype];

  precalc_llike0(op,data,lowcus,chain,llike0);

  /* make the grid */
  grid = (double **)calloc(GRIDSIDE, sizeof(double *));
  grid[0] = (double *)calloc(GRIDSIDE*GRIDSIDE, sizeof(double));
  for (i = 1; i < GRIDSIDE; i++)
    grid[i] = grid[0] + i * GRIDSIDE;

  /* note:  change the following if you change GRIDSIDE! */

  temp_theta = 0.9 * test_theta;
  for (i = 0; i < GRIDSIDE; i++) {
    temp_r = 0.9 * test_r;
    for (j = 0; j < GRIDSIDE; j++) {
      grid[i][j] = fmodel_likelihood(op,data,temp_theta,temp_r,llike0,
				     lowcus,chain);
      temp_r = temp_r + 0.05 * test_r;
    }
    temp_theta = temp_theta + 0.05 * test_theta;
  }
 
  for (i = 0; i < GRIDSIDE; i++) {
    for (j = 0; j < GRIDSIDE; j++) {
      if (i == 2 && j == 2) continue;
      if (grid[i][j] > grid[2][2]+epsilon1) {
	fprintf(ERRFILE,"\nfound something! real max: %f, fake max: %f\n",
		grid[i][j],grid[2][2]);
	free(grid[0]);
	free(grid);
	free(llike0[0]);
	free(llike0);
	return(TRUE);
      }
    }
  }

  free(grid[0]);
  free(grid);
  free(llike0[0]);
  free(llike0);
  return(FALSE);

} /* gridcheck */


/**********************************************************************
 * find_thetamax() finds the maximum of theta for the passed rec-rate */
double find_thetamax(option_struct *op, data_fmt *data, long lowcus,
		     long chain, double startrec, double starttheta, double **llike0,
		     double *llike)
{
  long chaintype, numloci, dfxplus, ddfxplus, repcount;
  double ltolerance, theta, rec, oldlike, newlike, change, dfx,
    ddfx, newtheta, newchange;
  boolean multilocus;

  multilocus = ((lowcus == -1) ? TRUE : FALSE);
  chaintype = TYPE_CHAIN(chain);
  ltolerance = log(epsilon1);
  numloci = getdata_numloci(op,data);
  theta = starttheta;
  rec = ((startrec > epsilon1) ? startrec : epsilon1); 
  repcount = 0;

  oldlike = fmodel_likelihood(op,data,theta,rec,llike0,lowcus,chain);

  while(1) {
    change = ((multilocus) ?
	      NRrec_locus_thetalderiv(op,data,theta,rec,llike0,&dfx,
				      &ddfx,&dfxplus,&ddfxplus) :
	      NRrec_thetalderiv(op,theta,rec,llike0,lowcus,chain,
				&dfx,&ddfx,&dfxplus,&ddfxplus));

    if (dfx < ltolerance) break;

    change = dfxplus*exp(change);
    newtheta = theta + change;
    if (newtheta < epsilon1) newtheta = epsilon1;
    newlike = fmodel_likelihood(op,data,newtheta,rec,llike0,
				lowcus,chain);

    for(newchange = -0.5*change; newlike < oldlike; newchange *= -0.5) {
      newtheta = theta + newchange;
      if (newtheta < epsilon1) newtheta = epsilon1;
      newlike = fmodel_likelihood(op,data,newtheta,rec,llike0,
				  lowcus,chain);
    }

    if (newlike == oldlike) {
      repcount++;
      if (repcount == REPEAT_TOLERANCE) break;
    } else repcount = 0;

    oldlike = newlike;
    theta = newtheta;
  }

  *llike = oldlike;
  return(theta);

} /* find_thetamax */


/********************************************************************
 * find_recmax() finds the maximum of rec-rate for the passed theta */
double find_recmax(option_struct *op, data_fmt *data, long lowcus,
		   long chain, double startrec, double starttheta, double **llike0,
		   double *llike)
{
  long chaintype, numloci, dfxplus, ddfxplus, repcount;
  double ltolerance, theta, rec, oldlike, newlike, change, dfx,
    ddfx, newrec, newchange;
  boolean multilocus;

  multilocus = ((lowcus == -1) ? TRUE : FALSE);
  chaintype = TYPE_CHAIN(chain);
  ltolerance = log(epsilon1);
  numloci = getdata_numloci(op,data);
  theta = starttheta;
  rec = ((startrec > epsilon1) ? startrec : epsilon1); 
  repcount = 0;

  oldlike = fmodel_likelihood(op,data,theta,rec,llike0,lowcus,chain);

  while(1) {
    change = ((multilocus) ?
	      NRrec_locus_reclderiv(op,data,theta,rec,llike0,&dfx,
				    &ddfx,&dfxplus,&ddfxplus) :
	      NRrec_reclderiv(op,theta,rec,llike0,lowcus,chain,
			      &dfx,&ddfx,&dfxplus,&ddfxplus));

    if (dfx < ltolerance) break;

    change = dfxplus*exp(change);
    newrec = rec + change;
    if (newrec < epsilon1) newrec = epsilon1;
    newlike = fmodel_likelihood(op,data,theta,newrec,llike0,
				lowcus,chain);

    for(newchange = -0.5*change; newlike < oldlike; newchange *= -0.5) {
      newrec = rec + newchange;
      if (newrec < epsilon1) newrec = epsilon1;
      newlike = fmodel_likelihood(op,data,theta,newrec,llike0,
				  lowcus,chain);
    }

    if (newlike == oldlike) {
      repcount++;
      if (repcount == REPEAT_TOLERANCE) break;
    } else repcount = 0;

    oldlike = newlike;
    rec = newrec;
  }

  *llike = oldlike;
  return(rec);

} /* find_recmax */


/******************************************************************
 * bisectfindth() is a helper function for profile_estimate() and *
 * performs a simple bisection line search in theta.              *
 * "startlike" must be less than "endlike".                       */
double bisectfindth(option_struct *op, data_fmt *data, long lowcus,
		    long chain, double start, double startlike, double end,
		    double endlike, double donelike, double **llike0, double rec)
{
  long lastchain;
  double mid, midlike, dummy;

  lastchain = totchains-1;

  mid = (start+end)/2.0;
  dummy = find_recmax(op,data,lowcus,lastchain,rec,mid,llike0,&midlike);

  while(!(donelike-epsilon1 < midlike && midlike < donelike+epsilon1)) {
    if (donelike > midlike) {
      start = mid;
      startlike = midlike;
    } else {
      end = mid;
      endlike = midlike;
    }
    mid = (start+end)/2.0;
    dummy = find_recmax(op,data,lowcus,lastchain,rec,mid,llike0,&midlike);
  } 

  return(mid);

} /* bisectfindth */


/*******************************************************************
 * bisectfindrec() is a helper function for profile_estimate() and *
 * performs a simple bisection line search in rec-rate.            *
 * "startlike" must be less than "endlike".                        */
double bisectfindrec(option_struct *op, data_fmt *data, long lowcus,
		     long chain, double start, double startlike, double end,
		     double endlike, double donelike, double **llike0, double th)
{
  long lastchain;
  double mid, midlike, dummy;

  lastchain = totchains-1;

  mid = (start+end)/2.0;
  dummy = find_thetamax(op,data,lowcus,lastchain,mid,th,llike0,&midlike);

  while(!(donelike-epsilon1 < midlike && midlike < donelike+epsilon1)) {
    if (donelike > midlike) { 
      start = mid;
      startlike = midlike;
    } else {
      end = mid;
      endlike = midlike;
    }
    mid = (start+end)/2.0;
    dummy = find_thetamax(op,data,lowcus,lastchain,mid,th,llike0,&midlike);
  } 

  return(mid);

} /* bisectfindrec */


/****************************************************************
 * dobounds() calculates the upper and lower bounds for a given *
 * likelihood value in "theta" and "rec-rate".                  */
void dobounds(option_struct *op, data_fmt *data, long lowcus,
	      double bestth, double bestrec, double **llike0, double donelike,
	      double maxlike, double *thlb, double *thub, double *reclb,
	      double *recub)
{
  long lastchain;
  double start, end, startlike, endlike = 0, dummy;

  lastchain = totchains - 1;

  /* using simple modified-Newton-Raphson to find max's */
  /* don't calculate confidence intervals for user fixed parameters */
  if (op->holding != 1) {
    start = bestth*0.001;
    dummy = find_recmax(op,data,lowcus,lastchain,bestrec,
			start,llike0,&startlike);
    end = bestth;
    endlike = maxlike;
    *thlb = bisectfindth(op,data,lowcus,lastchain,start,startlike,end,
			 endlike, donelike,llike0,bestrec);
  } else
    *thlb = bestth;

  if (op->holding != 2) {
    start = bestrec;
    do {
      start *= 0.001;
      dummy = find_thetamax(op,data,lowcus,lastchain,start,
			    bestth,llike0,&startlike);
    } while(start > EPSILON && startlike >= donelike);
    if (start <= EPSILON) {
      *reclb = 0.0;
    } else {
      end = bestrec;
      *reclb = bisectfindrec(op,data,lowcus,lastchain,start,startlike,
			     end,endlike,donelike,llike0,bestth);
    }
  } else
    *reclb = bestrec;

  if (op->holding != 1) {
    start = bestth*1000.0;
    dummy = find_recmax(op,data,lowcus,lastchain,bestrec,
			start,llike0,&startlike);
    end = bestth;
    endlike = maxlike;
    *thub = bisectfindth(op,data,lowcus,lastchain,start,startlike,end,
			 endlike,donelike,llike0,bestrec);
  } else
    *thub = bestth;

  if (op->holding != 2) {
    start = MAX(bestrec*1000.0,1.0);
    dummy = find_thetamax(op,data,lowcus,lastchain,start,
			  bestth,llike0,&startlike);
    end = bestrec;
    *recub = bisectfindrec(op,data,lowcus,lastchain,start,startlike,end,
			   endlike,donelike,llike0,bestth);
  } else
    *recub = bestrec;

} /* dobounds */


/****************************************************************
 * profile_estimate() is the main driver for profile likelihood *
 * calculation and output.                                      */
void profile_estimate(option_struct *op, data_fmt *data, long lowcus,
		      double bestth, double bestrec, double **llike0)
{
  long lastchain;
  double maxlike, donelike, thlb, thub, reclb, recub;

  lastchain = totchains - 1;
  maxlike = fmodel_likelihood(op,data,bestth,bestrec,llike0,lowcus,lastchain);

  printf("\ncalculating approximte 75%% confidence intervals\n");
  donelike = maxlike - (0.575364/2.0);
  dobounds(op,data,lowcus,bestth,bestrec,llike0,donelike,maxlike,&thlb,&thub,
	   &reclb,&recub);
  fprintf(outfile,"75%%\n");
  if (op->holding != 1) fprintf(outfile,"theta    | %e to %e\n",thlb,thub);
  else fprintf(outfile,"theta    | none calculated\n");
  if (op->holding != 2) fprintf(outfile,"rec-rate | %e to %e\n\n",reclb,recub);
  else fprintf(outfile,"rec-rate | none calculated\n");

  printf("\ncalculating approximte 95%% confidence intervals\n");
  donelike = maxlike - (5.99147/2.0);
  dobounds(op,data,lowcus,bestth,bestrec,llike0,donelike,maxlike,&thlb,&thub,
	   &reclb,&recub);
  fprintf(outfile,"95%%\n");
  if (op->holding != 1) fprintf(outfile,"theta    | %e to %e\n",thlb,thub);
  else fprintf(outfile,"theta    | none calculated\n");
  if (op->holding != 2) fprintf(outfile,"rec-rate | %e to %e\n\n",reclb,recub);
  else fprintf(outfile,"rec-rate | none calculated\n");

  printf("\ncalculating approximte 99%% confidence intervals\n");
  donelike = maxlike - (9.21034/2.0);
  dobounds(op,data,lowcus,bestth,bestrec,llike0,donelike,maxlike,&thlb,&thub,
	   &reclb,&recub);
  fprintf(outfile,"99%%\n");
  if (op->holding != 1) fprintf(outfile,"theta    | %e to %e\n",thlb,thub);
  else fprintf(outfile,"theta    | none calculated\n");
  if (op->holding != 2) fprintf(outfile,"rec-rate | %e to %e\n\n",reclb,recub);
  else fprintf(outfile,"rec-rate | none calculated\n");

  /* pass the bounding box information to the likelihood curve printer,
     currently using the 99% confidence interval measures */
  op->thlb = thlb;
  op->thub = thub;
  op->reclb = reclb;
  op->recub = recub;
   
} /* profile_estimate */


/**********************************************************************
 * rec_estimate is the main calling function for the point estimation *
 * procedures.  Pass lowcus = -1 for multi-locus estimation.          */
void rec_estimate(option_struct *op, data_fmt *data, long lowcus,
		  long chain, boolean locusend)
{
  double th, rec, like_at_max, **llike0, dummy;
  long i, numloci, chaintype;
#if DISTRIBUTION
  long nrecs[1000], refchain;
#endif

  numloci = getdata_numloci(op,data);

  chaintype = TYPE_CHAIN(chain);
  llike0 = (double **)calloc(1,numloci*sizeof(double *));
  llike0[0] = (double *)calloc(1,numloci*op->numout[chaintype]*sizeof(double));
  for(i = 1; i < numloci; i++) 
    llike0[i] = llike0[0] + i*op->numout[chaintype];

  precalc_llike0(op,data,lowcus,chain,llike0);

  if (lowcus == -1) {
    th = theti[numloci-1][chain];
    rec = reci[numloci-1][chain];
  } else {
    th = theti[lowcus][chain];
    rec = reci[lowcus][chain];
  }

#if !ALWAYS_REJECT
  if (lowcus == -1 && op->progress) printf("Beginning final estimator:\n");

#ifdef MAC
  eventloop();
#endif
  if (op->holding == 1) 
    rec = find_recmax(op,data,lowcus,chain,rec,th,llike0,&dummy);
  else if (op->holding == 2)
    th = find_thetamax(op,data,lowcus,chain,rec,th,llike0,&dummy);
  else
    broydenpoint(op,data,lowcus,chain,&th,&rec,llike0);
#ifdef MAC
  eventloop();
#endif

  if (op->progress) {
    if (lowcus != -1) printf("\nTheta=%10.5f, r=%10.5f, ",th,rec);
    else printf("\nFinal estimates: Theta=%10.5f, r=%10.5f\n",th,rec);
#if DISTRIBUTION
    fprintf(outfile,"\nchain %ld estimates are: ", apps+1);
    if (op->holding != 1) fprintf(outfile,"th = %e", th);
    if (op->holding != 2) fprintf(outfile,", rec = %e",rec);
    fprintf(outfile," accepted %ld/%ld trees", slacc,slid-NUMBURNIN);
    if (numdropped) fprintf(outfile,", dropped %ld", numdropped);
    if (swap) fprintf(outfile, ", swapped %ld/%ld\n", swacc,swap);
    else fprintf(outfile,"\n");
    if (op->haplotyping)
      fprintf(outfile,", %ld/%ld haplotype changes\n", hacc,hap);
#endif
#if !DISTRIBUTION
    fprintf(outfile,"  %ld      ",apps+1);
    if (op->holding != 1) fprintf(outfile,"%8.6f    ", th);
    else fprintf(outfile,"    const   ");
    if (op->holding != 2) fprintf(outfile,"%8.6f ",rec);
    else fprintf(outfile,"    const   ");
    fprintf(outfile,"accepted %ld/%ld trees", slacc,slid-NUMBURNIN);
    if (numdropped) fprintf(outfile,", dropped %ld", numdropped);
    if (swap) fprintf(outfile, ", swapped %ld/%ld", swacc,swap);
    if (op->haplotyping)
      fprintf(outfile,", %ld/%ld haplotype changes", hacc,hap);
    fprintf(outfile,"\n");
#endif

  }
#endif

  if (chain == totchains-1) {
    /*   rec_likeprint(op, data, lowcus, chain); */
    /* debug DEBUG warning WARNING */
    op->print_recbythmaxcurve = TRUE;
    if (op->print_recbythmaxcurve) 
      rec_likeprint2(op, data, th, lowcus, chain);
#ifdef MAC
    eventloop();
#endif
    if (op->mhmcsave) rec_scoreprint(op,data,lowcus,chain,FALSE);
#if DISTRIBUTION
    refchain = REF_CHAIN(chain);
    fprintf(outfile,"\nSampled # of recombinations by tree\n");
    for(i = 0; i < 40; i++) nrecs[i] = 0;
    for(i = 0; i < numtrees; i++) 
      nrecs[(long)(sum[locus][refchain][i].numrecombs)]++;
    for(i = 0; i < 40; i++) fprintf(outfile,"%ld recs occured %ld times\n",
				    i,nrecs[i]);
    fprintf(outfile,"\n");
#endif
  }

  if (lowcus != -1) {
    /* set the values for the next chain to use */
    if (th > 0 && rec > 0) {
      if (op->holding != 1) theti[lowcus][chain+1] = th;
      else theti[lowcus][chain+1] = theti[lowcus][chain];
      if (op->holding != 2) reci[lowcus][chain+1] = rec;
      else reci[lowcus][chain+1] = reci[lowcus][chain];
    } else {
      fprintf(ERRFILE,"WARNING, point estimation of parameters failed!\n");
      fprintf(ERRFILE,"using previous chain's parameter estimates\n");
      theti[locus][chain+1] = theti[lowcus][chain];
      reci[locus][chain+1] = reci[lowcus][chain];
    }
  }

  if (locusend) {
    if (op->holding == 1) {
      fprintf(outfile,"\n\nTheta was held constant at %e\n\n",
	      theti[lowcus][chain]);
    }
    if (op->holding == 2) {
      fprintf(outfile,"\n\nRec-rate was held constant at %e\n\n",
	      reci[lowcus][chain]);
    }
    like_at_max = model_likelihood(op,data,th,rec,lowcus,chain);
    if (gridcheck(op,data,th,rec,lowcus,chain)) { 
      fprintf(outfile,"WARNING, Maximizer failed:  you may wish to try");
      fprintf(outfile," again with a different random seed.\n");
      /* kludge needed for table printing in rec_liketable() */
      op->thlb = th/10.0;
      op->thub = 10.0*th;
      op->reclb = rec/10.0;
      op->recub = 10.0*rec;
    } else {
      if (lowcus == -1) {
	fprintf(outfile,"-----------------------------------------------\n");
	fprintf(outfile,"\n\nOver %ld loci, best estimate of parameters:\n",
		numloci);
	if (!MENU && op->holding == 2)
	  fprintf(outfile,"   Theta = %f Recombination = *****\n",th);
	else
	  fprintf(outfile,"   Theta = %f Recombination = %f\n",th,rec);
	fprintf(outfile,"Log(Likelihood) of the estimate %f",like_at_max);
      } else {
	fprintf(outfile,"\nPoint estimate of parameters for locus %ld\n",
		lowcus+1);
	if (!MENU && op->holding == 2)
	  fprintf(outfile,"Theta = %f Recombination = *****\n",th);
	else
	  fprintf(outfile,"Theta = %f Recombination = %f\n",th,rec);
	fprintf(outfile,"Log(Likelihood) at maximum %f\n\n",
		like_at_max);
      }
      if (op->profile) {
	fprintf(outfile,"\n\n-------------------------------------------\n");
	fprintf(outfile,"Approximate Confidence Intervals\n");
	fprintf(outfile,"-------------------------------------------\n");
	profile_estimate(op,data,lowcus,th,rec,llike0);
      } else {
	/* kludge needed for table printing in rec_liketable() */
	op->thlb = th/10.0;
	op->thub = 10.0*th;
	op->reclb = rec/10.0;
	op->recub = 10.0*rec;
      }
    } 
#if !ALWAYS_REJECT
    if (chain == totchains-1)
      rec_liketable(op, data, th, rec, chain, lowcus, 
		    (boolean)(MENU && op->progress));
#ifdef MAC
    eventloop();
#endif
    fprintf(outfile,"\n\n");
    print_locusplot(op,data,lowcus);
    if (chain == totchains-1)
      rec_linkgraph(op, data, sum, chain, lowcus);
#endif
  }

  free(llike0[0]);
  free(llike0);

}  /* rec_estimate */

#define MINTHETA 0.0001
#define MAXTHETA 10.0

/****************************************************************
 * rec_liketable() prints the likelihood curve in rec and theta *
 * for the chain "chain".                                       */
void rec_liketable(option_struct *op, data_fmt *data, double th,
		   double rec, long chain, long lowcus, boolean to_screen)
{
  long i, j, chaintype, garbage, numloci;
  double ltheta, *theta, *rrs, *llike, **dummy, **llike0;
  double lowtheta, lowrec, hightheta, highrec, fract;
  FILE *location;
  boolean done, logtheta;
  char inchar;
  long numfix, numrec;

  chaintype = TYPE_CHAIN(chain);
  numloci = getdata_numloci(op,data);

  llike0 = (double **)calloc(1,numloci*sizeof(double *));
  llike0[0] = (double *)
    calloc(1,numloci*op->numout[chaintype]*sizeof(double));
  for(i = 1; i < numloci; i++) 
    llike0[i] = llike0[0] + i*op->numout[chaintype];

  precalc_llike0(op,data,lowcus,chain,llike0);

  /* if the user holds theta constant, then only show
     single dimensional curve in rec-rate */
  if (op->holding == 1) numfix = 0;
  else numfix = 9;

  /* if the user holds rec-rate constant, then only show
     single dimensional curve in theta */
  if (op->holding == 2) {
    numrec = 1;
    rrs = (double *)calloc(numrec,sizeof(double));
    rrs[0] = rec;
  } else {
    numrec = 6;
    rrs = (double *)calloc(numrec,sizeof(double));
    rrs[0] = 0.0;
    rrs[1] = op->reclb;
    rrs[2] = 0.25 * (3*op->reclb+op->recub);
    rrs[3] = rec;
    rrs[4] = 0.25 * (op->reclb+3*op->recub);
    rrs[5] = op->recub;
    qsort((void *)rrs,6L,sizeof(double),doublecmp);
  }

  garbage = 0; /* just in case */
  dummy = NULL;

  llike = (double *)calloc(numtrees,sizeof(double));
  theta = (double *)calloc(numfix+1,sizeof(double));

  theta[numfix] = th;

  fract = (op->thub-op->thlb)/(double)(numfix-1);
  for (i = 0; i < numfix; i++) {
    theta[i] = op->thlb + ((double)i * fract);
  }

  done = FALSE;
  logtheta = TRUE;

  do {

    if (to_screen) location=stdout;
    else {
      location=outfile;
      if (op->progress) fprintf(stdout,"Table printed to file.\n");
    }

    /* now do the actual likelihood calculations */

    fprintf(location,
	    "\n  ============================================================================\n");
    fprintf(location,"             ln(Likelihood) for various values of Theta and rec\n");
    fprintf(location,
	    "\n  ============================================================================\n");
    fprintf(location,"                                   rec\n");
    fprintf(location,"\n   Theta  ||  ");
    for(i = 0; i < numrec; i++) {
      fprintf(location," %9.7f",rrs[i]);
      if (rrs[i]==rec) fprintf(location,"*");
      else fprintf(location," ");
    }
    fprintf(location,
	    "\n  ============================================================================\n");
    for(i = 0; i < numfix; i++) {
      fprintf(location,"%10.4f ||",theta[i]);
      for(j = 0; j < numrec; j++) {
	ltheta = 
	  fmodel_likelihood(op,data,theta[i],rrs[j],llike0,lowcus,chain);
	if (ltheta <= -1000.0) fprintf(location,"     ----  ");
	else fprintf(location," %10.4f",ltheta);
      }
      fprintf(location,"\n");
      if (theta[numfix]>theta[i] && theta[numfix]<=theta[i+1]) {
        fprintf(location,"%10.4f*|| ",theta[numfix]);
        for(j = 0; j < numrec; j++) {
	  ltheta = 
	    fmodel_likelihood(op,data,theta[numfix],rrs[j],llike0,lowcus,chain);
	  if (ltheta <= -1000.0) fprintf(location,"    ----  ");
	  else fprintf(location,"%10.4f",ltheta);
	  if (rrs[j]==rec)fprintf(location,"*");
	  else fprintf(location," ");
        }
        fprintf(location,"\n");
      }
    }
    fprintf(location,
	    "  ============================================================================\n");

    fprintf(location,"Dashed lines indicate lnL values less than -1000;\n");
    fprintf(location,"Asterisks indicate the MLEs of Theta and rec\n");

    if (to_screen) {
      fprintf(location,"\nAccept this table? [y/n]\n");
      scanf("%c%*[^\n]",&inchar);
      getchar();
      if(inchar=='y' || inchar=='Y') to_screen=FALSE;
      else {
	/* change Theta values in table */
	fprintf(location,"\nChange Theta values? [y/n]\n");
	scanf("%c%*[^\n]",&inchar);
	getchar();
	if(inchar=='y' || inchar=='Y') {
	  fprintf(location,"\nUse logarithmic scale? [y/n]\n");
	  scanf("%c%*[^\n]",&inchar);
	  getchar();
	  if (inchar=='y' || inchar=='Y') logtheta = TRUE;
	  else logtheta = FALSE;
	  do {
	    fprintf(location,"Lower bound on Theta?\n");
	    scanf("%lf%*[^\n]",&lowtheta);
	    getchar();
	    if (lowtheta <= 0.0) fprintf(location,"Illegal value\n");
	  } while (lowtheta <= 0.0);
	  if (logtheta) lowtheta = log(lowtheta);
	  do {
	    fprintf(location,"Upper bound on Theta?\n");
	    scanf("%lf%*[^\n]",&hightheta);
	    getchar();
	    if (hightheta <= lowtheta) fprintf(location,"Illegal value\n");
	  } while (hightheta <= lowtheta);
	  if (logtheta) hightheta = log(hightheta);
	  fract = (hightheta-lowtheta)/(double)(numfix-1);
	  for (i = 0; i < numfix; i++) {
	    theta[i] = lowtheta + ((double)i * fract);
	    if (logtheta) theta[i] = exp(theta[i]);
	  }
	}
	/* change rec values in table */
	fprintf(location,"\nChange rec values? [y/n]\n");
	scanf("%c%*[^\n]",&inchar);
	getchar();
	if(inchar=='y' || inchar=='Y') {
	  fprintf(location,"Lower bound on rec?\n");
	  scanf("%lf%*[^\n]",&lowrec);
	  getchar();
	  do {
	    fprintf(location,"Upper bound on rec?\n");
	    scanf("%lf%*[^\n]",&highrec);
	    getchar();
	    if (highrec <= lowrec) fprintf(location,"Illegal value\n");
	  } while (highrec <= lowrec);
	  fract = (highrec - lowrec)/3.0;
	  for (i = 0; i < 4; i++) rrs[i] = lowrec + ((double)i * fract);
	  for (i = 0; i < 4; i++) {
	    if (rrs[i] > 0.0) {
	      for (j = 4; j > i; j--) {
		rrs[j] = rrs[j-1];
	      }
	      rrs[i] = 0.0;
	      break;
	    }
	  }
	  for (i = 0; i < 5; i++) {
	    if (rrs[i] > rec) {
	      for (j = 5; j > i; j--) {
		rrs[j] = rrs[j-1];
	      }
	      rrs[i] = rec;
	      break;
	    }
	  }
	}
      }
    } else done=TRUE;

  } while (!done);

  free(llike);
  free(llike0[0]);
  free(llike0);
  free(theta);
  free(rrs);

} /* rec_liketable */


/*************************************************************
 * printhistplot() prints an ascii histogram to outfile.     *
 *                                                           *
 * It uses the scaled number of "sites" to bound the Y-axis, *
 * and HIST_MAXWD to bound the X-axis of the histogram.      *
 *                                                           *
 * XOUT and YOUT set the number of columns or lines between  *
 * tick marks on the respective axis (end point included).   *
 *                                                           *
 * LEFT_MARGIN is the size of the LEFT_MARGIN in characters. */
void printhistplot(option_struct *op, data_fmt *data, long *bars,
		   long htmult, long wdmult)
{
  long i, line, xtick, ytick, yscale, numticks;

  yscale = countsites(op,data) / htmult;

  /* print the histogram, one line at a time. */

  /* start with the x-axis labels */
  numticks = (HIST_MAXWD+2)/XOUT;
  for(i = 0; i < LEFT_MARGIN+1; i++) fprintf(outfile," ");
  for(i = 1; i < numticks+1; i++) fprintf(outfile,"%5ld",i*XOUT*wdmult);
  fprintf(outfile,"\n");

  /* now the x-axis */
  for(i = 0; i < LEFT_MARGIN; i++) fprintf(outfile," ");
  fprintf(outfile,"+");
  for(i = 0, xtick = 1; i < HIST_MAXWD+2; i++) {
    if (xtick == XOUT) {fprintf(outfile,"+"); xtick = 1;}
    else {fprintf(outfile,"-"); xtick++;}
  }
  fprintf(outfile,"\n");

  /* now the histogram body */
  ytick = yscale;
  for(line = 0; line < yscale; line++) {
    /* print the y-axis and label (if any) for this line */
    if(yscale-line == ytick) {
      fprintf(outfile,"%5ld+",line);
      ytick -= YOUT;
    } else {
      for(i = 0; i < LEFT_MARGIN; i++) fprintf(outfile," ");
      fprintf(outfile,"|");
    }
    /* print the bar */
    for(i = 0; i < HIST_MAXWD; i++)
      fprintf(outfile,"%s",((i < bars[line]) ? "*" : " "));
    fprintf(outfile,"\n");
  }

} /* printhistplot */


/*************************************************************
 * rechistout() prints an ascii histogram to outfile.        *
 *                                                           *
 * It uses the two constants HIST_MAXHT and HIST_MAXWD to set*
 * the boundaries of the histogram.                          */
void rechistout(option_struct *op, data_fmt *data, cutpointrec *cutp)
{
  long i, j, site, sbins[HIST_MAXWD+1], wdmult, htmult, *temp, numsites,
    maxbinht, minbinht, line, ytick, xtick;

  numsites = countsites(op,data);
  temp = (long *)calloc(numsites,sizeof(long));

  for(i = 0; i < HIST_MAXWD+1; i++) sbins[i] = 0;

  for(i = 0; i < cutp->numrecombs; i++) temp[cutp->sitescore[i]]++;

  for(wdmult = 1; numsites/wdmult > HIST_MAXWD; wdmult++)
    ;

  minbinht = RECOMB_MAX;
  maxbinht = 0;
  for(i = 0, site = 0; i < HIST_MAXWD && site < numsites; i++) {
    sbins[i] = temp[site++];
    for(j = 1; j < wdmult && site < numsites; j++) {
      sbins[i] += temp[site++];
    }
    if (sbins[i] > maxbinht) maxbinht = sbins[i];
    if (sbins[i] && (sbins[i] < minbinht)) minbinht = sbins[i];
  }

  minbinht--; /* leave at least 1 element in all bins with elements */
  for(i = 0; i < HIST_MAXWD; i++) if (sbins[i]) sbins[i] -= minbinht;
  maxbinht -= minbinht;

  for(htmult = 1; maxbinht/htmult > HIST_MAXHT; htmult++) 
    ;

  for(i = 0; i < HIST_MAXWD; i++) sbins[i] /= htmult;

  /* print the histogram header */
  fprintf(outfile,"\n\n");
  fprintf(outfile,"Histogram of recombination events along the sequence\n");
  fprintf(outfile,"\n");
  fprintf(outfile,"Scaling: x-axis = %ld (sites per tick)\n",wdmult);
  fprintf(outfile,"         y-axis = %ld (events per tick)\n",htmult);
  fprintf(outfile,"         y-axis start at %ld events\n",minbinht);

  /* print the histogram, one line at a time */
  ytick = HIST_MAXHT;
  for(line = 0; line < HIST_MAXHT; line++) {
    if (HIST_MAXHT-line == ytick) {
      fprintf(outfile,"%5ld+",(HIST_MAXHT-line)*htmult+minbinht);
      ytick -= YOUT;
    } else {
      for(i = 0; i < LEFT_MARGIN; i++) fprintf(outfile," ");
      fprintf(outfile,"|");
    }
    for(i = 0; i < HIST_MAXWD;i++)
      fprintf(outfile,"%s",((HIST_MAXHT-line <= sbins[i]) ? "*" : " "));
    fprintf(outfile,"\n");
  }

  /* print the x-axis, then the associated labels */
  for(i = 0; i < LEFT_MARGIN; i++) fprintf(outfile," ");
  fprintf(outfile,"+");
  for(i = 0, xtick = 1, j = 0; i < HIST_MAXHT+2; i++) {
    if (xtick == XOUT) {fprintf(outfile,"+"); xtick = 1; j++;}
    else {fprintf(outfile,"-"); xtick++;}
  }
  fprintf(outfile,"\n");
  for(i = 0; i < LEFT_MARGIN; i++) fprintf(outfile," ");
  fprintf(outfile," ");
  for(i = 1; i < j+1; i++) fprintf(outfile,"%5ld",i*XOUT*wdmult);
  fprintf(outfile,"\n\n");

  /* debug DEBUG warning WARNING */
  maxbinht += minbinht;
  for(htmult = 1; maxbinht/htmult > HIST_MAXWD; htmult++) ;
  for(i = 0; i < numsites; i++) temp[i] /= htmult;

  /*printhistplot(op,data,temp,1L,htmult); graph usually too big!! */
  free(temp);

} /* rechistout */


/*********************************************************************
 * rec_linkgraph() outputs the histogram of recombination cut points *
 * to the outfile and for use by Mathematica into a mathlinkfile.   */
void rec_linkgraph(option_struct *op, data_fmt *data, treerec ***triis,
		   long chain, long lowcus)
{
  long i, chaintype, trii, numr;
  treerec *tr;
  cutpointrec cutp;

  if(lowcus == -1) return;

  chaintype = TYPE_CHAIN(chain);

  cutp.numrecombs = 0;
  for(trii = 0; trii < op->numout[chaintype]; trii++) {
    cutp.numrecombs += triis[lowcus][chaintype][trii].numrecombs;
  }

  if(cutp.numrecombs == 0) {
    fprintf(outfile,"\n\nNo recombinations present");
    fprintf(outfile," in sampled trees.\n\n");
    return;
  }

  cutp.sitescore = (long *)calloc(cutp.numrecombs,sizeof(long));

  for(trii = 0, i = 0; trii < op->numout[chaintype]; trii++) {
    tr = &triis[lowcus][chaintype][trii];
    for(numr = 0; numr < tr->numrecombs; numr++)
      cutp.sitescore[i++] = tr->sitescore[numr];
  }

  /* output to outfile */
  rechistout(op,data,&cutp);

  /* output to mathlinkfile */

  free(cutp.sitescore);

} /* rec_linkgraph */


/******************************************************************
 * runlike() produces output to be processed by Mathematica using *
 * something like Peter's mixing notebook script.                 */
void runlike(option_struct *op, data_fmt *data, treerec ***triis,
	     long lowcus)
{
  long trii, chain, chaintype, refchain;
  double llike, th, rec, lth, lrec;
  FILE *rlike;
  treerec *tr;

  rlike = fopen("runlike","w+");

  for(chain = 0; chain < totchains; chain++) {
    chaintype = TYPE_CHAIN(chain);
    refchain = REF_CHAIN(chain);
    th = theti[lowcus][chain];
    rec = reci[lowcus][chain];
    lth = log(th);
    lrec = log(rec);
    for(trii = 0; trii < op->numout[chaintype]; trii++) {
      tr = &triis[lowcus][refchain][trii];
      llike = rectreellike(op,th,rec,lth,lrec,lowcus,chain,trii);
      fprintf(rlike," %12.6f %12.6f",tr->llike,llike);
    }
    fprintf(rlike,"\n\n");
  }

  fclose(rlike);

} /* runlike */

/************************************************************
 * Computes ln(P(G|params)) for use in chain heating        */
double coalprob(option_struct *op, data_fmt *data, tree *tr, double th, double r)
{
  tlist *t;
  double tk=0, ts=0;
  double tyme, prob;

  for (t = tr->tymelist; t->succ != NULL; t=t->succ) {
    tyme = t->age - t->eventnode->tyme;
    tk += t->numbranch * (t->numbranch - 1.0) * tyme;
    ts += count_active_tlist(op,data,t) * tyme;
  }
  if (r)
    prob = -tk/th - r*ts + (log(2.0/th))*tr->numcoals + log(r)*tr->numrecombs;
  else
    prob = -tk/th + (log(2.0/th))*tr->numcoals;

  return(prob);

} /* coalprob */



//updaterecnum - need more work!!!
void updaterecnum(long start, long end, double utyme, recnumb *recnum){
    recnumb *cur;
    // branch fragment should be larger than or same as recnum fragment!
    for (cur = recnum; cur != NULL; cur=cur->next){
	if (cur->beg >= start && cur->end <= end) {
	    cur->swt += utyme;
	    if (cur->end < end) cur->last_weight += utyme;
	}
    }
    return;
}

void addrecs(option_struct *op, data_fmt *data, node *p, tlist *t, double tyme, recnumb *recnum){
    node *q; 
    q = p;
/*   if (!q->top) q = q->back; */
    recnumb *cur;
    if (istip(q)) {
	for (cur=recnum; cur != NULL; cur = cur->next){
	    cur->swt = cur->swt + tyme;
	  if (cur->next != NULL)
	      cur->last_weight = cur->last_weight + tyme;
	}
	return;
    }
  

  if (isrecomb(q)) updaterec(op,data,q,tyme,recnum);
  else updatecoal(op,data,q,tyme, recnum);
  return;
}

void updatecoal(option_struct *op, data_fmt *data, node *p, double tyme, recnumb *recnum)
{
  long newstart, newend;

  finds(p,&newstart,&newend);
  
  if (newstart > newend) printf("wrong site info\n");

  
  updaterecnum(newstart,newend,tyme, recnum);
  return;
} 


void updaterec(option_struct *op, data_fmt *data, node *p, double tyme, recnumb *recnum)
{
  updaterecnum(p->coal[1],p->coal[2*p->coal[0]],tyme,recnum);
  return;

}


void addrecnum(recnumb *recnum, long site){
  recnumb *cur,*temp;

  for (cur = recnum; cur != NULL; cur = cur->next){
    if (cur->beg <= site && cur->end > site){
	temp = (struct recnumb *)malloc(sizeof(struct recnumb));
	temp->beg = site + 1;
	temp->end = cur->end;
	temp->swt = cur->swt;
	temp->numrec = cur->numrec;
	temp->next = cur->next;

	cur->end = site;
	cur->numrec = 1;
	cur->next = temp;
	return;
    }
    if (cur->end == site){
	cur->numrec ++;
	return;
    }
  }
}


struct sumtree *treescore(option_struct *op, data_fmt *data,tree *target){
  tlist *t;
  recnumb *recnums,*temp,*curtemp,*cur;
  long i, j,  entries;
  double temptyme,tyme;
  struct sumtree *tempsum;

  tempsum = (struct sumtree *)malloc(sizeof(struct sumtree));

  temp = (struct recnumb *)malloc(sizeof(struct recnumb));
  recnums = temp;
  cur = recnum0;
  temp->beg = cur->beg;
  temp->end = cur->end;
  temp->next = NULL;
  temp->numrec = 0;
  temp->swt = 0;
  while (1)
    {
    if (cur->next == NULL) break;
    cur = cur->next;
    curtemp = (struct recnumb *)malloc(sizeof(struct recnumb));
    temp->next = curtemp;
    temp = temp->next;
    temp->beg = cur->beg;
    temp->end = cur->end;
    temp->next = NULL;
    temp->numrec = 0;
    temp->swt = 0;
    }

  temptyme = 0.0;
  t = target->tymelist;
  /* count the tymelist entries */
  entries = 0;
  while (t!=NULL) {
    entries++;
    t = t->succ;
  }
 
  t = curtree->tymelist;
  
  for(i=0,j=0,temptyme =0.0;i<entries;i++) {
    if (t->numbranch == 1) break;
    tyme = t->age - t->eventnode->tyme;
    temptyme += t->numbranch * (t->numbranch - 1) * tyme;
    tempsum->wlinks += temptyme;
    addrec(op,data,t,tyme, recnums);

    t = t->succ;
  }

  tempsum->numcoal = curtree->numcoals;
  tempsum->recinfo = recnums;  
  
  return tempsum;
}



/**********************************************************************
 * rec_estimate is the main calling function for the point estimation *
 * procedures.  Pass lowcus = -1 for multi-locus estimation.          */
void rec_estimateH(option_struct *op, data_fmt *data, long lowcus, long chain, boolean locusend)
{
  double th, rec, like_at_max, **llike0, dummy;
  long i, numloci, chaintype;
#if DISTRIBUTION
  long nrecs[1000], refchain;
#endif

  numloci = getdata_numloci(op,data);

  chaintype = TYPE_CHAIN(chain);
  llike0 = (double **)calloc(1,numloci*sizeof(double *));
  llike0[0] = (double *)calloc(1,numloci*op->numout[chaintype]*sizeof(double));
  for(i = 1; i < numloci; i++) 
    llike0[i] = llike0[0] + i*op->numout[chaintype];

  precalc_llike0(op,data,lowcus,chain,llike0);

  if (lowcus == -1) {
    th = theti[numloci-1][chain];
    rec = reci[numloci-1][chain];
  } else {
    th = theti[lowcus][chain];
    rec = reci[lowcus][chain];
  }

#if !ALWAYS_REJECT
  if (lowcus == -1 && op->progress) printf("Beginning final estimator:\n");

#ifdef MAC
  eventloop();
#endif
  if (op->holding == 1) 
    rec = find_recmax(op,data,lowcus,chain,rec,th,llike0,&dummy);
  else if (op->holding == 2)
    th = find_thetamax(op,data,lowcus,chain,rec,th,llike0,&dummy);
  else
    broydenpoint(op,data,lowcus,chain,&th,&rec,llike0);
#ifdef MAC
  eventloop();
#endif

  if (op->progress) {
    if (lowcus != -1) printf("\nTheta=%10.5f, r=%10.5f, ",th,rec);
    else printf("\nFinal estimates: Theta=%10.5f, r=%10.5f\n",th,rec);
#if DISTRIBUTION
    fprintf(outfile,"\nchain %ld estimates are: ", apps+1);
    if (op->holding != 1) fprintf(outfile,"th = %e", th);
    if (op->holding != 2) fprintf(outfile,", rec = %e",rec);
    fprintf(outfile," accepted %ld/%ld trees", slacc,slid-NUMBURNIN);
    if (numdropped) fprintf(outfile,", dropped %ld", numdropped);
    if (swap) fprintf(outfile, ", swapped %ld/%ld\n", swacc,swap);
    else fprintf(outfile,"\n");
    if (op->haplotyping)
      fprintf(outfile,", %ld/%ld haplotype changes\n", hacc,hap);
#endif
#if !DISTRIBUTION
    fprintf(outfile,"  %ld      ",apps+1);
    if (op->holding != 1) fprintf(outfile,"%8.6f    ", th);
    else fprintf(outfile,"    const   ");
    if (op->holding != 2) fprintf(outfile,"%8.6f ",rec);
    else fprintf(outfile,"    const   ");
    fprintf(outfile,"accepted %ld/%ld trees", slacc,slid-NUMBURNIN);
    if (numdropped) fprintf(outfile,", dropped %ld", numdropped);
    if (swap) fprintf(outfile, ", swapped %ld/%ld", swacc,swap);
    if (op->haplotyping)
      fprintf(outfile,", %ld/%ld haplotype changes", hacc,hap);
    fprintf(outfile,"\n");
#endif

  }
#endif

  if (chain == totchains-1) {
    /*   rec_likeprint(op, data, lowcus, chain); */
    /* debug DEBUG warning WARNING */
    op->print_recbythmaxcurve = TRUE;
    if (op->print_recbythmaxcurve) 
      rec_likeprint2(op, data, th, lowcus, chain);
#ifdef MAC
    eventloop();
#endif
    if (op->mhmcsave) rec_scoreprint(op,data,lowcus,chain,FALSE);
#if DISTRIBUTION
    refchain = REF_CHAIN(chain);
    fprintf(outfile,"\nSampled # of recombinations by tree\n");
    for(i = 0; i < 40; i++) nrecs[i] = 0;
    for(i = 0; i < numtrees; i++) 
      nrecs[(long)(sum[locus][refchain][i].numrecombs)]++;
    for(i = 0; i < 40; i++) fprintf(outfile,"%ld recs occured %ld times\n",
				    i,nrecs[i]);
    fprintf(outfile,"\n");
#endif
  }

  if (lowcus != -1) {
    /* set the values for the next chain to use */
    if (th > 0 && rec > 0) {
      if (op->holding != 1) theti[lowcus][chain+1] = th;
      else theti[lowcus][chain+1] = theti[lowcus][chain];
      if (op->holding != 2) reci[lowcus][chain+1] = rec;
      else reci[lowcus][chain+1] = reci[lowcus][chain];
    } else {
      fprintf(ERRFILE,"WARNING, point estimation of parameters failed!\n");
      fprintf(ERRFILE,"using previous chain's parameter estimates\n");
      theti[locus][chain+1] = theti[lowcus][chain];
      reci[locus][chain+1] = reci[lowcus][chain];
    }
  }

  if (locusend) {
    if (op->holding == 1) {
      fprintf(outfile,"\n\nTheta was held constant at %e\n\n",
	      theti[lowcus][chain]);
    }
    if (op->holding == 2) {
      fprintf(outfile,"\n\nRec-rate was held constant at %e\n\n",
	      reci[lowcus][chain]);
    }
    like_at_max = model_likelihood(op,data,th,rec,lowcus,chain);
    if (gridcheck(op,data,th,rec,lowcus,chain)) { 
      fprintf(outfile,"WARNING, Maximizer failed:  you may wish to try");
      fprintf(outfile," again with a different random seed.\n");
      /* kludge needed for table printing in rec_liketable() */
      op->thlb = th/10.0;
      op->thub = 10.0*th;
      op->reclb = rec/10.0;
      op->recub = 10.0*rec;
    } else {
      if (lowcus == -1) {
	fprintf(outfile,"-----------------------------------------------\n");
	fprintf(outfile,"\n\nOver %ld loci, best estimate of parameters:\n",
		numloci);
	if (!MENU && op->holding == 2)
	  fprintf(outfile,"   Theta = %f Recombination = *****\n",th);
	else
	  fprintf(outfile,"   Theta = %f Recombination = %f\n",th,rec);
	fprintf(outfile,"Log(Likelihood) of the estimate %f",like_at_max);
      } else {
	fprintf(outfile,"\nPoint estimate of parameters for locus %ld\n",
		lowcus+1);
	if (!MENU && op->holding == 2)
	  fprintf(outfile,"Theta = %f Recombination = *****\n",th);
	else
	  fprintf(outfile,"Theta = %f Recombination = %f\n",th,rec);
	fprintf(outfile,"Log(Likelihood) at maximum %f\n\n",
		like_at_max);
      }
      if (op->profile) {
	fprintf(outfile,"\n\n-------------------------------------------\n");
	fprintf(outfile,"Approximate Confidence Intervals\n");
	fprintf(outfile,"-------------------------------------------\n");
	profile_estimate(op,data,lowcus,th,rec,llike0);
      } else {
	/* kludge needed for table printing in rec_liketable() */
	op->thlb = th/10.0;
	op->thub = 10.0*th;
	op->reclb = rec/10.0;
	op->recub = 10.0*rec;
      }
    } 
#if !ALWAYS_REJECT
    if (chain == totchains-1)
      rec_liketable(op, data, th, rec, chain, lowcus, 
		    (boolean)(MENU && op->progress));
#ifdef MAC
    eventloop();
#endif
    fprintf(outfile,"\n\n");
    print_locusplot(op,data,lowcus);
    if (chain == totchains-1)
      rec_linkgraph(op, data, sum, chain, lowcus);
#endif
  }

  free(llike0[0]);
  free(llike0);

}  /* rec_estimate */

double coalprobH(option_struct *op, data_fmt *data, tree *tr, double th, recnumb *recs, Rs *recrates)
{
  tlist *t;
  double tk=0, wts=0;
  double tyme, prob;
  recnumb *cur;
  Rs *rec;

  for (t = tr->tymelist; t->succ != NULL; t=t->succ) {
    tyme = t->age - t->eventnode->tyme;
    tk += t->numbranch * (t->numbranch - 1.0) * tyme;
    wts += count_weighted_tlist(op,data,t) * tyme;
  }
  
  prob = -tk/th - wts + (log(2.0/th)*tr->numcoals);
  
  rec =recrates;
  for(cur = recs; cur->next !=NULL; cur = cur->next){
      prob = prob + log(rec->recrates) * cur->numrec;
      rec = rec->next;
  }
  return prob;
}

 
	    
void finds(node *p, long *start, long *end){
  *start = p->coal[1];
  *end = p->coal[2*p->coal[0]]; 
  return;
}

  
struct recnumb *scorerecs(option_struct *op, data_fmt *data, tree *target)
{
  tlist *t;
  recnumb *recnums;
  double temptyme,tyme;  

  recnums = (struct recnumb *)malloc(sizeof(struct recnumb));

  recnums->beg = 0;
  recnums->end = countsites(op,data)-1;
  recnums->next = NULL;
  recnums->numrec = 0;
  recnums->swt = 0.0;
  recnums->last_weight = 0.0;

  temptyme = 0.0;
  t = target->tymelist;

  while (1){
    tyme = t->succ->eventnode->tyme - t->eventnode->tyme;
    if (isrecomb(t->eventnode)) updaterecnumbs(t->eventnode,recnums); // insert a new recnumb struct 
    addrec(op,data,t,tyme,recnums);
    t = t->succ;
    if (t->succ == NULL) break;
  }
  return recnums;

} 

void updaterecnumbs(node *p, recnumb *recnumbs){
  recnumb *newfragment,*cur;
  node *q;
  long recsite;
  //find recombination site

  q = findunique(p);
  
  if (q->next->recstart > q->next->next->recstart){
    recsite = q->next->recstart - 1;
  }
  else{
    recsite = q->next->next->recstart - 1;
  }

  //update recnumb

  for (cur = recnumbs; cur != NULL;cur=cur->next){
    if (cur->end == recsite){ // recombination at same site no update 
      cur->numrec++;
      return; 
    }
    if (cur->beg <= recsite && cur->end > recsite){ //new recombination site
      newfragment = (struct recnumb *)malloc(sizeof(struct recnumb));
      newfragment->beg = recsite+1;
      newfragment->end = cur->end; 
      newfragment->swt = cur->swt;
      newfragment->numrec = cur->numrec;
      newfragment->next = cur->next;
      newfragment->last_weight = cur->last_weight;
      cur->end = recsite;
      cur->numrec = 1;
      cur->next = newfragment;
      cur->last_weight = cur->swt;
      return;
    }
  }
  printf("error at updaterecnumbs\n");
  return;
}

void addrec(option_struct *op, data_fmt *data, tlist *t, double tyme, recnumb *recnum){
  long i; 
  for(i = 0;i<t->numbranch;i++){
    if (t->branchlist[i]->coal[0] != 0) addrecs(op,data,t->branchlist[i], t, tyme, recnum);
  }
}

   
void EMestimate(option_struct *op, long chain, double *parameters){
  double oldlikelihood = 0, newlikelihood;
  int count = 0;
  long numsamples, chaintype;

  chaintype = TYPE_CHAIN(chain);
  numsamples = op->numout[chaintype];

  while(1){
    //*thetav = getthetav(op, recsv, thetav, chain);
    getrecs(op, parameters, chain);
    newlikelihood = treellikelihood(op, parameters, chain);
    
    if ((newlikelihood < oldlikelihood) && oldlikelihood != 0) printf("error at EM\n");
    if (newlikelihood == oldlikelihood) count++;
    oldlikelihood = newlikelihood;
    if (count == 5) return;
    //return;
  }
}
 
double getthetav(option_struct *op, double *parameters, long chain){
  double sumbranch = 0, sumcoal = 0,*weightarray, maxweight = 0, minweight = 1,c, newtheta = 0;
  long i, numsamples, refchain, chaintype;
  treerec *tr = NULL;

  chaintype = TYPE_CHAIN(chain);
  refchain = REF_CHAIN(chain);
  numsamples = op->numout[chaintype];

  weightarray = (double *)calloc(numsamples,sizeof(double));

  for (i=0 ; i<numsamples; i++){    
    tr = (&sum[0][refchain][i]);
    weightarray[i] = gettrweight(tr, parameters) ;
    if (weightarray[i] < minweight) minweight = weightarray[i];
    if (weightarray[i] > maxweight) maxweight = weightarray[i];
  }
  
  c = (maxweight + minweight)/2.0;

  for (i=0 ; i < numsamples; i++){
    sumbranch += tr->tk * weightarray[i]/c;
    sumcoal += tr->numcoals * weightarray[i]/c;
  }
  parameters[numspaces] = sumbranch/sumcoal;
  free(weightarray);
  return newtheta;
}

void getrecs(option_struct *op, double *parameters, long chain){
  double sumsites = 0, sumrecs = 0, *weightarray, maxweight = 0, minweight = 1,c;
  long i,j,k, refchain, chaintype, numsamples;
  treerec *tr;

  chaintype = TYPE_CHAIN(chain);
  refchain = REF_CHAIN(chain);
  numsamples = op->numout[chaintype];

  weightarray = (double *)calloc(numsamples,sizeof(double));

  for (i = 0; i<numspaces; i++){
    sumsites = sumrecs = 0;
    maxweight = 0;
    minweight = 1;
    for (j = 0; j<numsamples; j++){
      tr = (&sum[0][refchain][j]);
      weightarray[j] = gettrweight(tr, parameters);
      if (weightarray[j] < minweight) minweight = weightarray[j];
      if (weightarray[j] > maxweight) maxweight = weightarray[j];
    
    }
    c = (maxweight + minweight)/2.0;
    for( k = 0; k<numsamples; k++){
      tr = (&sum[0][refchain][k]);
      sumsites += tr->sumsite[i] * weightarray[k]/c;
      sumrecs += tr->rec[i] * weightarray[k]/c;
    }
    parameters[i] = sumrecs / sumsites;
    printf("pos %ld - %ld rec %lf\n",spaceinfo[i], spaceinfo[i+1], parameters[i]);
  }
  free(weightarray);
}


double treellikelihood(option_struct *op, double *parameters, long chain){
  long numsamples, i, j;
  long chaintype,refchain;
  double likelihood = 0;
  treerec *tr;


  chaintype = TYPE_CHAIN(chain);
  refchain = REF_CHAIN(chain);
  numsamples = op->numout[chaintype];

  for (i = 0; i<numsamples; i++){
    tr = (&sum[0][refchain][i]);   
    likelihood = pow((theta0/(parameters[numspaces])),tr->numcoals) * exp(tr->tk * (1.0/(parameters[numspaces]) - 1.0/theta0));

    for (j = 0; j<numspaces; j++){
      likelihood = likelihood * pow(parameters[j]/recarray[spaceinfo[j]], tr->rec[j]) * exp(tr->sumsite[j] * (parameters[j] - recarray[spaceinfo[j]]));
    }
  }
  likelihood = log(likelihood);
  return likelihood;
}


double gettrweight(treerec *tr, double *parameters){
  double weight;
  long j;
  
  weight = pow((theta0/(parameters[numspaces])),tr->numcoals) * exp(tr->tk * (1.0/(parameters[numspaces]) - 1.0/theta0));
  
  for (j = 0; j < numspaces; j++){
    weight = weight * pow(parameters[j]/recarray[spaceinfo[j]], tr->rec[j]) * exp(tr->sumsite[j] * (parameters[j] - recarray[spaceinfo[j]]));
  }
  return weight;
}


void rec_estimatev(option_struct *op, data_fmt *data, long lowcus, long chain, boolean locusend){
  double *tempp = NULL;//r1,r2...theta
  //theta0 = tempp[numspaces];
  
  free(tempp);

  return;
}

void scoretreeH(option_struct *op, data_fmt *data, long chain)
{
  tlist *t;
  treerec *trii;
  long i, j, refchain, chaintype, entries;
  double temptyme,tyme;

  refchain = REF_CHAIN(chain);
  chaintype = TYPE_CHAIN(chain);

  temptyme = 0.0;
  t = curtree->tymelist;
  /* count the tymelist entries */
  entries = 0;
  while (t!=NULL) {
    entries++;
    t = t->succ;
  }
  t = curtree->tymelist;

  /* allocate space for them */
  trii = &sum[locus][refchain][op->numout[chaintype]-1];
  if (trii->eventtype) free(trii->eventtype);
  trii->eventtype = (long *)calloc(1,entries*sizeof(long));
  if (trii->sitescore) free(trii->sitescore);
  if (trii->rec) free(trii->rec);
  if (trii->sumsite) free(trii->sumsite);
  /* the +1 is to guarantee a minimum allocation of size 1 */
  trii->sitescore = (long *)calloc(curtree->numrecombs+1,sizeof(long));

  trii->tk = 0.0;
  trii->ts = 0.0;

  for(i=0,j=0,temptyme =0.0;i<entries;i++) {
    if (t->numbranch == 1) break;
    tyme = t->age - t->eventnode->tyme;
    trii->eventtype[i] = isrecomb(t->eventnode);
    trii->tk += t->numbranch * (t->numbranch - 1) * tyme;
    trii->ts += count_active_tlist(op,data,t) * tyme;

    if (trii->eventtype[i]) {
      trii->sitescore[j] = findlink(t->eventnode);
      j++;
    }
    t = t->succ;
  }
  
  //addrecevent(trii);


  trii->numcoals = curtree->numcoals;
  trii->numrecombs = curtree->numrecombs;
  trii->llike = curtree->likelihood;

  trii->recs = NULL;


  trii->rec = (double *)calloc(seq_length,sizeof(double));
  trii->sumsite = (double *)calloc(seq_length,sizeof(double));

  memcpy(trii->rec, curtree->numrec_array,seq_length*sizeof(double));
  memcpy(trii->sumsite,curtree->weight_array,seq_length*sizeof(double));

  //check rec
  
/*   for (i = 0; i<seq_length; i++){ */
/*     //printf("x %ld %ld %ld\n",i,trii->rec[i],curtree->numrec_array[i]); */
/*     if (trii->rec[i] != curtree->numrec_array[i]) printf("x %ld %ld %ld\n",i,trii->rec[i],curtree->numrec_array[i]); */
/*     //if (trii->sumsite[i] != curtree->weight_array[i]) printf("y %ld\n",i); */
/*   } */

/*   trii->rec = (long *)calloc(seq_length,sizeof(long)); */
/*   trii->sumsite = (double *)calloc(seq_length,sizeof(double)); */

/*   memcpy(trii->rec, curtree->numrec_array,seq_length*sizeof(long)); */
/*   memcpy(trii->sumsite,curtree->numrec_array,seq_length*sizeof(double)) */
  //transformtree(data, trii);

} 

void addrecevent(treerec *trii){
  tlist *t;
  node *curnode, *pnode1, *pnode2;
  long i, pos;
  double *temprec,*tempsumsite;


/*   trii->sumsite = (double *)calloc(seq_length,sizeof(double)); */
/*   trii->rec = (long *)calloc(seq_length,sizeof(long)); */

  temprec = (double *)calloc(seq_length,sizeof(double));
  tempsumsite = (double *)calloc(seq_length,sizeof(double));

  for (i = 0; i<seq_length; i++){
    temprec[i] = 0;
    tempsumsite[i] = 0.0;
  }


  for (t = curtree->tymelist->succ; t != NULL; t = t->succ){
    curnode = findunique(t->eventnode);
    if (iscoal(curnode)){
      pnode1 = curnode->next->back;
      pnode2 = curnode->next->next->back;
      
      for (i = 1; i <= pnode1->coal[0]; i++){
	for (pos = pnode1->coal[1]; pos < pnode1->coal[2 * i]; pos++){
	  tempsumsite[pos] += pnode1->length;
	}
      }     
      for (i = 1; i <= pnode2->coal[0]; i++){
	for (pos = pnode2->coal[1]; pos < pnode2->coal[2 * i]; pos++){
	  tempsumsite[pos] += pnode2->length;
	}
      }
    }

    else{
      pnode1 = curnode->back;  
      for (i = 1; i <= pnode1->coal[0]; i++){
	for (pos = pnode1->coal[1]; pos < pnode1->coal[2 * i]; pos++){
	  tempsumsite[pos] += pnode1->length;
	}
      }
      
      if (curnode->next->recstart == 0) temprec[curnode->next->recend]++;
      else temprec[curnode->next->recstart - 1]++;
    }
  }
  trii->sumsite = tempsumsite;
  trii->rec = temprec;
}
      


void getparameters(option_struct *op, long chain, double *weight0array){

  //oldlikelihood = getlogrelativelikelihood(op,chain,temptheta,temprecrates,templamda,weight0array);

  //estimate theta
  if (op->holdings[0] ==0)  temptheta = gettheta(op,chain,temptheta,temprecrates,templamda,weight0array); 
  //newlikelihood = getlogrelativelikelihood(op,chain,temptheta,temprecrates,templamda,weight0array);
  //if (newlikelihood < oldlikelihood) printf("theta\n");
  //oldlikelihood = newlikelihood;

  //estimate recrate
  if (op->holdings[1] == 0) temprecrates[0] = getrecrate(op,chain,temptheta,temprecrates,templamda,weight0array,0);

  //newlikelihood = getlogrelativelikelihood(op,chain,temptheta,temprecrates,templamda,weight0array);
  //if (newlikelihood < oldlikelihood) printf("rec1\n");
  //oldlikelihood = newlikelihood;

  if (op->holdings[2] == 0) temprecrates[1] = getrecrate(op,chain,temptheta,temprecrates,templamda,weight0array,1);

  // newlikelihood = getlogrelativelikelihood(op,chain,temptheta,temprecrates,templamda,weight0array);
  //if (newlikelihood < oldlikelihood) printf("rec2\n");
  //oldlikelihood = newlikelihood;

  //estimate lamd

  //if (chain > 0){
  if (op->holdings[3] == 0) templamda[0] = getlamda(op,chain,temptheta,temprecrates,templamda,weight0array,0);

  //newlikelihood = getlogrelativelikelihood(op,chain,temptheta,temprecrates,templamda,weight0array);
  ///if (newlikelihood < oldlikelihood) printf("lamda1\n");
  //oldlikelihood = newlikelihood;

  if (op->holdings[4] == 0)  templamda[1] = getlamda(op,chain,temptheta,temprecrates,templamda,weight0array,1);

  //newlikelihood = getlogrelativelikelihood(op,chain,temptheta,temprecrates,templamda,weight0array);
  //if (newlikelihood < oldlikelihood) printf("lamda2\n");
  //oldlikelihood = newlikelihood;


  //}
  //printf("rec_c %lf, rec_h %lf, lamda c %lf, lamda h %lf\n",temprecrates[0],temprecrates[1],templamda[0],templamda[1]);
  return;
}

  

void parameter_estimation_EM(option_struct *op, long chain, double theta, double *recrates, double *lamda){
  double oldlikelihood =0, newlikelihood = 0;
  double *weight0array;
  long numsamples, refchain, chaintype,i; 
  double temp;
  
  chaintype = TYPE_CHAIN(chain);
  refchain = REF_CHAIN(chain);
  numsamples = op->numout[chaintype];
  
  temptheta = theta;
  temprecrates[0] = recrates[0];
  temprecrates[1] = recrates[1];
  templamda[0] = lamda[0];
  templamda[1] = lamda[1];

  alpha0 = (double *)calloc(1, seq_length * sizeof(double));
  alpha1 = (double *)calloc(1, seq_length * sizeof(double));
  beta0 = (double *)calloc(1, seq_length * sizeof(double));
  beta1 = (double *)calloc(1, seq_length * sizeof(double));


  //checktrees(op, chain);

  //faketrees(op,chain);

  weight0array = (double *)calloc(1, op->numout[chaintype] * sizeof(double));
  
  getlogweight0array(op,chain,theta,recrates,lamda,weight0array);
  oldlikelihood = getlogrelativelikelihood(op,chain,temptheta,temprecrates,templamda,weight0array);
  

  for (i = 0; i < 20; i++){
    //printf(" %lf %lf\n",oldlikelihood,newlikelihood);
    getparameters(op,chain,weight0array);
    newlikelihood = getlogrelativelikelihood(op,chain,temptheta,temprecrates,templamda,weight0array);

    if (oldlikelihood > newlikelihood){
      //printf("error  %lf %lf\n",oldlikelihood,newlikelihood);
    break;
    }


    if ((oldlikelihood / newlikelihood) > 0.99) break;
    else oldlikelihood = newlikelihood;
    
  }
  
  if (temprecrates[0] > temprecrates[1]){
    temp = temprecrates[0];
    temprecrates[0] = temprecrates[1];
    temprecrates[1] = temp;
    temp = templamda[0];
    templamda[0] = templamda[1];
    templamda[1] = temp;
  }

  printf("\n theta %lf ", temptheta);
  printf("recrates %lf, %lf ", temprecrates[0], temprecrates[1]);
  printf("lamda %lf, %lf\n", templamda[0], templamda[1]);
  
  //drawlikelihoodcurve(op,chain,temptheta,temprecrates,templamda,weight0array);
  /* if (chain == 5) { */
/*     likesurface_rec1(op,chain,weight0array); */
/*     likesurface_lamda1(op,chain,weight0array); */
/*     //likesurface_rc_rh(op,chain, weight0array); */
/*     //likesurface_rec(op,chain,weight0array); */
/*     //likesurface_theta(op,chain, weight0array); */
/*   } */
  theta0 = temptheta;
  recrates[0] = temprecrates[0];
  recrates[1] = temprecrates[1];
  lamda[0] = templamda[0];
  lamda[1] = templamda[1];

  findoptimal_temp(op,chain, theta0, recrates, lamda, weight0array);
  hotstateprob(op,chain, theta0, recrates,lamda, weight0array);
  free(weight0array);
  if (chain != 13){
    if (lamda[0] == 1.0) lamda[0] = 0.9999;
    if (lamda[1] == 1.0) lamda[1] = 0.9999;
    if (recrates[0] == 0.0) recrates[0] = 0.0001;
  }
  free(alpha0);
  free(alpha1);
  free(beta0);
  free(beta1);

  return;
}



double getloghmmlikelihood(double *recrates, double *lamda, treerec *tr){
  double logalpha0 = 1;
  double logalpha1 = 1;
  double temp0,temp1,mean,likelihood;
  long i;
  double ini[2];

  ini[0] = (1 - lamda[1]) / (1 - lamda[0] + 1 - lamda[1]);
  ini[1] = 1 - ini[0];


  logalpha0 = log(ini[0]) + log(recrates[0]) * tr->rec[0] - (tr->sumsite[0] * recrates[0]);
  logalpha1 = log(ini[1]) + log(recrates[1]) * tr->rec[0] - (tr->sumsite[0] * recrates[1]);

  
  for (i = 1; i<(seq_length - 1); i++){
    
    mean = (logalpha0 + logalpha1) / 2.0;

    temp0 = log(exp(logalpha0 + log(lamda[0]) - mean) + exp(logalpha1 + log(1 - lamda[1]) - mean)) + (log(recrates[0]) * tr->rec[i] - (tr->sumsite[i] * recrates[0])) + mean;

    temp1 = log(exp(logalpha0 + log(1 - lamda[0]) - mean) + exp(logalpha1 + log(lamda[1]) - mean)) + (log(recrates[1]) * tr->rec[i] - (tr->sumsite[i] * recrates[1])) + mean;

    logalpha0 = temp0;
    logalpha1 = temp1;   
    
  }
  mean = (logalpha0 + logalpha1) / 2.0;
  likelihood = log(exp(logalpha0 - mean) + exp(logalpha1 - mean)) + mean;  
  return likelihood;
}


double getlogthetaweight(double theta, treerec *tr){
  double weight = 0;

  weight =  log(2/theta) * tr->numcoals - (tr->tk/theta);
  
  return weight;
}


double gettheta(option_struct *op, long chain, double theta, double *recrates, double *lamda, double *logweight0array){
  double expectedcoal = 0;
  double expectedsum = 0;
  long i,numsamples, refchain, chaintype;
  treerec *tr;
  double *weight, max = 0;
  double temptk =0, tempcoal = 0;
  

  chaintype = TYPE_CHAIN(chain);
  refchain = REF_CHAIN(chain);
  numsamples = op->numout[chaintype];

  weight = (double *)calloc(numsamples,sizeof(double));

  
  for (i = 0; i < numsamples; i++){
    tr = (&sum[0][refchain][i]);
    weight[i] = getlogthetaweight(theta, tr) + getloghmmlikelihood(recrates, lamda, tr) - logweight0array[i];
  }
  
  for (i = 0; i < numsamples; i++){
    if (max < weight[i]) max = weight[i];
  }

  for (i = 0; i < numsamples; i++){
    tr = (&sum[0][refchain][i]);
    temptk += tr->tk;
    tempcoal += tr->numcoals;
    
    expectedcoal +=  tr->numcoals * exp(weight[i] - max);
    expectedsum += tr->tk  * exp(weight[i] - max);
  }
  free(weight);
  return expectedsum /  expectedcoal;
}

double getlamda(option_struct *op, long chain, double theta, double *recrates, double *lamda, double *weight0array, int state){
  int otherstate;
  long i,numsamples, refchain, chaintype;
  double *values;
  treerec *tr;
  double ss = 0.0;
  double sd = 0.0;
  double *weight, *num, *den, max = 1;

  chaintype = TYPE_CHAIN(chain);
  refchain = REF_CHAIN(chain);
  numsamples = op->numout[chaintype];

  values = (double *)calloc(2,sizeof(double));
  num = (double *)calloc(numsamples,sizeof(double));
  den = (double *)calloc(numsamples,sizeof(double));

  weight = (double *)calloc(numsamples,sizeof(double));

  if (state == 0) otherstate = 1;
  else otherstate = 0;

  for (i = 0; i < numsamples; i++){
    tr = (&sum[0][refchain][i]);

    getalphas(tr, recrates, lamda);
    getbetas(tr, recrates, lamda);

    weight[i] = getlogthetaweight(theta, tr)  - weight0array[i];
    expectedtransitionrate(recrates, lamda, tr, state ,values,weight[i]);
    num[i] = values[0];
    den[i] = values[1];
    //printf("lamda %d, %lf %lf\n",state, values[0],values[1]);
  }

  for (i = 0; i<numsamples; i++){
    if (max < num[i]) max = num[i];
    if (max < den[i]) max = den[i];
  }

  //weight[i] = getlogthetaweight(theta, tr) + getloghmmlikelihood(recrates, lamda, tr) - logweight0array[i]
  for (i = 0; i<numsamples; i++){
    tr = (&sum[0][refchain][i]);
    ss += exp(num[i] - max);
    sd += exp(den[i] - max);
    //sd += exp(den[i] - max) + exp(weight[i] + getloghmmlikelihood(recrates,lamda,tr) + log((1-lamda[otherstate])/(2-lamda[0]-lamda[1])) - max) - exp(weight[i] + getloghmmlikelihood(recrates,lamda,tr) + log(1.0/2.0) - max);
  }
  free(weight);
  free(num);
  free(den);
  free(values);

  return ss/sd;
}

double getrecrate(option_struct *op, long chain, double theta, double *recrates, double *lamda, double *weight0array, int state){
  double ser = 0, sel = 0;
  double *values;
  long i,numsamples, refchain, chaintype;
  double *weight, *num, *den, max = 0;
  treerec *tr;
  
  chaintype = TYPE_CHAIN(chain);
  refchain = REF_CHAIN(chain);
  numsamples = op->numout[chaintype];

  values = (double *)calloc(2,sizeof(double));  
  weight = (double *)calloc(numsamples,sizeof(double));
  num = (double *)calloc(numsamples,sizeof(double));
  den = (double *)calloc(numsamples,sizeof(double));
  

  for (i = 0; i < numsamples; i++){
    tr = (&sum[0][refchain][i]);

    getalphas(tr, recrates, lamda);
    getbetas(tr, recrates, lamda);

    weight[i] = getlogthetaweight(theta, tr) - weight0array[i];
    expectedrecrate(recrates, lamda, tr, state, values, weight[i]);
    num[i] = values[0];
    den[i] = values[1];
    //printf("val %lf, %lf, ratio %lf\n",values[0],values[1],ser/sel);
  }
  
  for (i = 0; i<numsamples; i++){
    if (max < num[i]) max = num[i];
    if (max < den[i]) max = den[i];
  }

  for (i = 0; i<numsamples; i++){
    ser += exp(num[i] - max);
    sel += exp(den[i] - max);
  }
  free(weight);
  free(num);
  free(den);
  free(values);
  return ser/sel;
}
 
void getlogweight0array(option_struct *op, long chain, double theta_0, double *recrates_0, double *lamda_0, double *weight0array){
  long i,numsamples, refchain, chaintype;
  treerec *tr;

  chaintype = TYPE_CHAIN(chain);
  refchain = REF_CHAIN(chain);
  numsamples = op->numout[chaintype];
  
  for(i = 0; i< numsamples; i++){
    tr = (&sum[0][refchain][i]);
    weight0array[i] = getlogthetaweight(theta_0, tr) + getloghmmlikelihood(recrates_0, lamda_0, tr);
  }
}

double getlogrelativelikelihood(option_struct *op, long chain, double theta, double *recrates, double *lamda, double *weight0array){
  double rlikelihood = 0;
  long i,numsamples, refchain, chaintype;
  treerec *tr;

  chaintype = TYPE_CHAIN(chain);
  refchain = REF_CHAIN(chain);
  numsamples = op->numout[chaintype];

  for (i = 0; i < numsamples; i++){
    tr = (&sum[0][refchain][i]);
    rlikelihood += exp(getlogthetaweight(theta, tr) + getloghmmlikelihood(recrates, lamda, tr) - weight0array[i]);
  }
  return log(rlikelihood / (numsamples));
}



/* void getalphas(treerec *tr, double *recrates, double *lamda){ */
/*   double temp0,temp1; */
/*   long i; */
/*   double ini[2]; */

/*   ini[0] = (1 - lamda[1]) / (1 - lamda[0] + 1 - lamda[1]); */
/*   ini[1] = 1 - ini[0]; */

/*   temp0 = ini[0] * pow(recrates[0],tr->rec[0]) * exp(-tr->sumsite[0] * recrates[0]); */
/*   temp1 = ini[1] * pow(recrates[1],tr->rec[0]) * exp(-tr->sumsite[0] * recrates[1]); */

/*   alpha0[0] = temp0 / (temp0 + temp1); */
/*   alpha1[0] = 1 - alpha0[0]; */

/*   for (i = 1; i<(seq_length - 1); i++){ */
/*     temp0 = (alpha0[i - 1] * lamda[0] + alpha1[i - 1] * (1 - lamda[1])) * pow(recrates[0],tr->rec[i]) * exp(-tr->sumsite[i] * recrates[0]); */
/*     temp1 = (alpha0[i - 1] * (1 - lamda[0]) + alpha1[i - 1] * lamda[1]) * pow(recrates[1],tr->rec[i]) * exp(-tr->sumsite[i] * recrates[1]); */
    
/*     alpha0[i] = temp0 / (temp0 + temp1); */
/*     alpha1[i] = 1.0 - alpha0[i]; */
/*     //printf("%ld alpha0 %lf\n",i,alpha0[i]); */
/*   } */
/*   return; */
/* } */

/* void getbetas(treerec *tr, double *recrates, double *lamda){ */
/*   long i; */
/*   double temp0,temp1; */

/*   beta0[seq_length - 1] = 1; */
/*   beta1[seq_length - 1] = 1; */

/*   for (i = seq_length - 2; i >= 0; i--){ */
/*     temp0 = beta0[i + 1] * lamda[0] * pow(recrates[0],tr->rec[i + 1]) * exp(-tr->sumsite[i + 1] * recrates[0])  */
/*       + beta1[i + 1] * (1 - lamda[0]) * pow(recrates[1],tr->rec[i + 1]) * exp(-tr->sumsite[i + 1] * recrates[1]); */

/*     temp1 = beta0[i + 1] * (1 - lamda[1]) * pow(recrates[0],tr->rec[i + 1]) * exp(-tr->sumsite[i + 1] * recrates[0])  */
/*       + beta1[i + 1] * lamda[1] * pow(recrates[1],tr->rec[i + 1]) * exp(-tr->sumsite[i + 1] * recrates[1]); */

/*     beta0[i] = temp0 / (temp0 + temp1); */
/*     beta1[i] = 1.0 - beta0[i]; */
/*   } */
/* return; */
/* } */


void expectedrecrate(double *recrates, double *lamda, treerec *tr, int state, double *values, double weight){
  long i;
  double *alphaarray, *betaarray;
  double *otheralpha, *otherbeta;
  

  if (state == 0){
    alphaarray = alpha0;
    betaarray = beta0;
    otheralpha = alpha1;
    otherbeta = beta1;
  }
  else{
    alphaarray = alpha1;
    betaarray = beta1;
    otheralpha = alpha0;
    otherbeta = beta0;
  }
  values[0] = values[1] = 0;

  for (i = 0; i<(seq_length - 1); i++){
    values[0] += exp(alphaarray[i] + betaarray[i] + weight) * tr->rec[i];
    values[1] += exp(alphaarray[i] + betaarray[i] + weight) * tr->sumsite[i];
  }
  values[0] = log(values[0]);
  values[1] = log(values[1]);
  

  return ;
}

void expectedtransitionrate(double *recrates, double *lamda, treerec *tr, int state, double *values, double weight){
  long i;
  double *alphaarray,*otheralpha, *betaarray, *otherbeta, probst, probstst;
  int otherstate;
  double *temp0, *temp1;
  double hmm_likelihood;

  temp0 = (double *)calloc(seq_length,sizeof(double));
  temp1 = (double *)calloc(seq_length,sizeof(double));

  hmm_likelihood = getloghmmlikelihood(recrates, lamda, tr);

  if (state == 0){
    otherstate = 1;
    alphaarray = alpha0;
    otheralpha = alpha1;
    betaarray = beta0;
    otherbeta = beta1;
  }
  else{
    otherstate = 0;
    alphaarray = alpha1;
    otheralpha = alpha0;
    betaarray = beta1;
    otherbeta = beta0;
  }
  values[0] = values[1] = 0;


  for (i = 0; i<(seq_length - 1); i++){

    probstst = alphaarray[i] + log(lamda[state]) + log(recrates[state]) * tr->rec[i+1] - (tr->sumsite[i+1] * recrates[state]) + betaarray[i + 1] + weight;
    probst =   alphaarray[i] + betaarray[i] + weight;


    //printf("%ld  c->c  %lf   c :%lf\n",i,alphaarray[i],betaarray[i]);
    //temp0[i] = probstst;
    //temp1[i] = probst;
    values[0] += exp(probstst);
    values[1] += exp(probst);

  }

/*   minvalue0 = temp0[0]; */
/*   minvalue1 = temp1[0]; */
/*   for (i = 0; i <(seq_length - 1); i++){ */
/*     if (minvalue0 > temp0[i]) minvalue0 = temp0[i]; */
/*     if (minvalue1 > temp1[i]) minvalue1 = temp1[i]; */
/*   } */

/*   for (i = 0; i<(seq_length - 1); i++){ */
/*     values[0] += exp(temp0[i] - minvalue0); */
/*     values[1] += exp(temp1[i] - minvalue1); */
/*   } */
/*   values[0] = log(values[0]) + minvalue0; */
/*   values[1] = log(values[1]) + minvalue1; */
  
  free(temp0);
  free(temp1);  
  
  values[0] = log(values[0]) + hmm_likelihood;
  values[1] = log(values[1]) + hmm_likelihood;
  
  //printf(" %lf, %lf\n",log(values[0]+values[1]),log(test));
  return;
}

void getalphas(treerec *tr, double *recrates, double *lamda){
  double mean;
  long i;
  double ini[2];

  ini[0] = (1 - lamda[1]) / (1 - lamda[0] + 1 - lamda[1]);
  ini[1] = 1 - ini[0];

  alpha0[0] = log(ini[0]) + log(recrates[0]) * tr->rec[0] - (tr->sumsite[0] * recrates[0]);
  alpha1[0] = log(ini[1]) + log(recrates[1]) * tr->rec[0] - (tr->sumsite[0] * recrates[1]);

  for (i = 1; i<(seq_length - 1); i++){

    mean = (alpha0[i-1] + alpha1[i-1]) / 2.0;
    
    alpha0[i] = log(exp(alpha0[i-1] - mean) * lamda[0]  + exp(alpha1[i-1] - mean) * (1 - lamda[1])) + log(recrates[0]) * tr->rec[i] - (tr->sumsite[i] * recrates[0]) + mean;

    alpha1[i] = log(exp(alpha0[i-1] - mean) * (1 - lamda[0])  + exp(alpha1[i-1] - mean) * lamda[1]) + log(recrates[1]) * tr->rec[i] - (tr->sumsite[i] * recrates[1]) + mean;					
					

  }
  return;
}

void getbetas(treerec *tr, double *recrates, double *lamda){
  long i;
  double mean;

  beta0[seq_length - 1] = log(1);
  beta1[seq_length - 1] = log(1);

  for (i = seq_length - 2; i >= 0; i--){
    mean = (beta0[i+1] + beta1[i+1]) / 2.0;
    
    beta0[i] = log(lamda[0] * pow(recrates[0],tr->rec[i+1]) * exp(- tr->sumsite[i+1] * recrates[0]) * exp(beta0[i+1] - mean) +
		   (1 - lamda[0]) * pow(recrates[1],tr->rec[i+1]) * exp(- tr->sumsite[i+1] * recrates[1]) * exp(beta1[i+1] - mean)) + mean;

    beta1[i] = log(lamda[1] * pow(recrates[1],tr->rec[i+1]) * exp(- tr->sumsite[i+1] * recrates[1]) * exp(beta1[i+1] - mean) +
		   (1 - lamda[1]) * pow(recrates[0],tr->rec[i+1]) * exp(- tr->sumsite[i+1] * recrates[0]) * exp(beta0[i+1] - mean)) + mean;

  }
  return;
}


void checktrees(option_struct *op, long chain){
  long i,j,numsamples, refchain, chaintype;
  long recinhot = 0,recincold = 0;
  double siteinhot = 0, siteincold= 0;
  treerec *tr;

  chaintype = TYPE_CHAIN(chain);
  refchain = REF_CHAIN(chain);
  numsamples = op->numout[chaintype];
  
  for(i = 0; i< numsamples; i++){
    tr = (&sum[0][refchain][i]);
    
    for (j = 0; j<seq_length; j++){
      if (j < 450 || j >550){
	recincold += tr->rec[j];
	siteincold += tr->sumsite[j];
      }
      else{
	recinhot += tr->rec[j];
	siteinhot += tr->sumsite[j];
      }
    }
  }
  printf("rate... %lf\n",(recincold + recinhot) / (double)(siteinhot + siteincold));

  return;

}

void faketrees(option_struct *op, long chain)
{
  treerec *tr;
  long i, j, refchain, chaintype,numsamples;
  double depth = 0.1;
  long site;

  refchain = REF_CHAIN(chain);
  chaintype = TYPE_CHAIN(chain);
  numsamples = op->numout[chaintype];

  for (i = 0; i < numsamples; i++){
    tr = (&sum[0][refchain][i]);
    for (j = 0; j < seq_length; j++){
      tr->sumsite[j] = depth;
      tr->rec[j] = 0;
    }
    for (j = 0; j < 90; j++){
      site = (long)(randum() * seq_length);
      tr->rec[site]++;
    }
    for (j = 0; j <50 ;j++){
      site = (long)(randum() * 50);
      tr->rec[site + 450]++;
    }
  }

  return;
}




void findoptimal_temp(option_struct *op, long chain, double theta,double *recrates, double *lamda, double *weight0array){
  double *delta0,*delta1;
  int *psi0,*psi1;
  int *temppath0,*temppath1,*optimalpath;
  double *weighttree;
  double ini,*temp0,temp,ave0,ave1;
  long i,j,numsamples, refchain, chaintype;
  treerec *tr;

  FILE *bestpath;

  chaintype = TYPE_CHAIN(chain);
  refchain = REF_CHAIN(chain);
  numsamples = op->numout[chaintype];

  bestpath = fopen("pathresult","w");

  delta0 = (double *)calloc(1, seq_length * sizeof(double));
  psi0 = (int *)calloc(1, seq_length * sizeof(int));
  delta1 = (double *)calloc(1, seq_length * sizeof(double));
  psi1 = (int *)calloc(1, seq_length * sizeof(int));
  weighttree = (double *)calloc(1, numsamples * sizeof(double));
  temppath0 =  (int *)calloc(1, seq_length * sizeof(int));
  temppath1 =  (int *)calloc(1, seq_length * sizeof(int));
  temp0 = (double *)calloc(1, numsamples * sizeof(double));

  for (i = 0; i<seq_length;i++) delta0[i] = delta1[i] = 0;


  for (i = 0; i<numsamples;i++){
    tr = (&sum[0][refchain][i]);
    weighttree[i] = getlogthetaweight(theta, tr) - weight0array[i];
  }

  ini = (1 - lamda[1]) / (1 - lamda[0] + 1 - lamda[1]);
  
  temp = 0;
  for (i = 0; i<numsamples; i++){    
    tr = (&sum[0][refchain][i]);
    temp0[i] = log(recrates[0]) * tr->rec[0] - (tr->sumsite[0] * recrates[0]) + weighttree[i];
    temp += temp0[i];
  }
  ave0 = temp/(double)numsamples;
 
  temp = 0;
  for (i = 0; i<numsamples; i++) temp += exp(temp0[i] - ave0);

  delta0[0] = ini * temp;
  

  temp = 0;
  for (i = 0; i<numsamples; i++){    
    tr = (&sum[0][refchain][i]);
    temp0[i] = log(recrates[1]) * tr->rec[0] - (tr->sumsite[0] * recrates[1]) + weighttree[i];
    temp += temp0[i];
  }
  ave1 = temp/(double)numsamples;
 
  temp = 0;
  for (i = 0; i<numsamples; i++) temp += exp(temp0[i] - ave1);

  delta1[0] = (1 - ini) * temp;

  
  delta0[0] = delta0[0] /(delta0[0] + delta1[0] * exp(ave1 - ave0));
  delta1[0] = 1 - delta0[0];

  psi0[0] = 0;
  psi1[0] = 1;


  for (j = 1; j < seq_length; j++){

    temp = 0;
    for (i = 0; i<numsamples; i++){    
      tr = (&sum[0][refchain][i]);
      temp0[i] = log(recrates[0]) * tr->rec[j] - (tr->sumsite[j] * recrates[0]) + weighttree[i];
      temp += temp0[i];
    }
    ave0 = temp/(double)numsamples;
    
    temp = 0;
    for (i = 0; i<numsamples; i++) temp += exp(temp0[i] - ave0);

    if ((delta0[j - 1] * lamda[0]) > delta1[j - 1] * (1 - lamda[1])){
      delta0[j] = delta0[j - 1] * lamda[0] * temp;
      memcpy(temppath0,psi0,seq_length);
    }
    else{
      //printf("%ld, 0->  %lf, %lf\n ",j, (delta0[j - 1] * lamda[0]) , delta1[j - 1] * (1 - lamda[1]));
      //printf("1->0");
      delta0[j] = delta1[j - 1] * (1 - lamda[1]) * temp; 
      memcpy(temppath0,psi1,seq_length*sizeof(int));
    }
    
    temppath0[j] = 0;


    temp = 0;
    for (i = 0; i<numsamples; i++){    
      tr = (&sum[0][refchain][i]);
      temp0[i] = log(recrates[1]) * tr->rec[j] - (tr->sumsite[j] * recrates[1]) + weighttree[i];
      temp += temp0[i];
    }
    ave1 = temp/(double)numsamples;
    
    temp = 0;
    for (i = 0; i<numsamples; i++) temp += exp(temp0[i] - ave1);


    if ((delta1[j - 1] * lamda[1]) > delta0[j - 1] * (1 - lamda[0])){
      delta1[j] = delta1[j - 1] * lamda[1] * temp;
      memcpy(temppath1,psi1,seq_length*sizeof(int));
    }
    else{
      //printf("%ld, 1->  %lf, %lf\n ",j, (delta1[j - 1] * lamda[1]) , delta0[j - 1] * (1 - lamda[0]));
      delta1[j] = delta0[j - 1] * (1 - lamda[0]) * temp; 
      memcpy(temppath1,psi0,seq_length*sizeof(int));
    }
    
    temppath1[j] = 1;
  
    
    delta0[j] = delta0[j] /(delta0[j] + delta1[j] * exp(ave1 - ave0));
    delta1[j] = 1 - delta0[j];


    memcpy(psi0,temppath0,seq_length*sizeof(int));
    memcpy(psi1,temppath1,seq_length*sizeof(int));

    //printf("%ld %d -> 0, %d ->1  %lf %lf \n",j,psi0[j-1],psi1[j-1],delta0[j],delta1[j]);
  }

  if(delta0[seq_length - 1] > delta1[seq_length - 1]) optimalpath = psi0;
  else optimalpath = psi1;


  for (i = 0 ; i < seq_length; i++){
    if(optimalpath[i] == 0) fprintf(bestpath,"%lf\n ", recrates[0]);
    else fprintf(bestpath,"%lf\n ",recrates[1]);
    
/*     if (i < 450 || i > 550) printf("0.1\n"); */
/*     else printf("0.5\n"); */
  }
  fclose(bestpath);
  free(delta0);
  free(delta1);
  free(psi0);
  free(psi1);
  free(weighttree);
  free(temppath0);
  free(temppath1);
  free(temp0);
  return;
}

/* void getparameters2(option_struct *op, long chain, double *weight0array){ */
  
/*   temptheta = getthetawithS(op,chain,temptheta,temprecrates,weight0array);  */

/*   temprecrates[0] = getrecrateswithS(op,chain,temptheta,temprecrates,weight0array,0); */
/*   temprecrates[1] = getrecrateswithS(op,chain,temptheta,temprecrates,weight1array,1); */
  
/*   return; */
/* } */



/* double getrecratewithS(option_struct *op, long chain, double theta, double *recrates,double *weight0array, int state){ */
/*   double ser = 0, sel = 0; */
/*   double *recnum, *sites; */
/*   double *values; */
/*   long i,numsamples, refchain, chaintype; */
/*   double *weight, *num, *den, max = 0; */
/*   treerec *tr; */
  
/*   chaintype = TYPE_CHAIN(chain); */
/*   refchain = REF_CHAIN(chain); */
/*   numsamples = op->numout[chaintype]; */

/*   values = (double *)calloc(2,sizeof(double));   */
/*   weight = (double *)calloc(numsamples,sizeof(double)); */
/*   num = (double *)calloc(numsamples,sizeof(double)); */
/*   den = (double *)calloc(numsamples,sizeof(double)); */
  

/*   for (i = 0; i < numsamples; i++){ */
/*     tr = (&sum[0][refchain][i]); */

/*     getalphas(tr, recrates, lamda); */
/*     getbetas(tr, recrates, lamda); */

/*     weight[i] = getlogthetaweight(theta, tr) - weight0array[i]; */
/*     values[0] = values[1] = 0; */
/*     for (j = 0; j<seq_length; j++){ */
/*       if (givenarray[j] == state){ */
/* 	values[0] += tr->recs[j]; */
/* 	values[1] += tr->sumsites[j]; */
/*       } */
/*       else  weight[i] += log(recrates[givenarray[j]]) * tr->recs[j] - tr->sumsites[j] * recrates[givenarray[j]]; */
/*     } */

/*     expectedrecrate(recrates, lamda, tr, state, values); */
/*     num[i] = log(values[0]) + weight[i]; */
/*     den[i] = log(values[1]) + weight[i]; */
/*     //printf("val %lf, %lf, ratio %lf\n",values[0],values[1],ser/sel); */
/*   } */
  
/*   for (i = 0; i<numsamples; i++){ */
/*     if (max < num[i]) max = num[i]; */
/*     if (max < den[i]) max = den[i]; */
/*   } */

/*   for (i = 0; i<numsamples; i++){ */
/*     ser += exp(num[i] - max); */
/*     sel += exp(den[i] - max); */
/*   } */
/*   free(weight); */
/*   free(num); */
/*   free(den); */
/*   free(values); */
/*   return ser/sel; */
/* } */

 
/* void getlogweight0arraywithS(option_struct *op, long chain, double theta_0, double *recrates_0, double *weight0array){ */
/*   long i,j,numsamples, refchain, chaintype; */
/*   treerec *tr; */

/*   chaintype = TYPE_CHAIN(chain); */
/*   refchain = REF_CHAIN(chain); */
/*   numsamples = op->numout[chaintype]; */
  
/*   for(i = 0; i< numsamples; i++){ */
/*     tr = (&sum[0][refchain][i]); */
/*     weight0array[i] = getlogthetaweight(theta_0, tr); */
/*     for (j = 0; j<seq_length; j++){ */
/*       weight0array[i] += log(recrates_0[givenrecarray[j]]) * tr->recs[j] - recrates_0[giverecrarray[j]] * tr->sumsites[i]; */
/*     } */
/*   } */
/*   return; */
/* } */


/* double gettheta(option_struct *op, long chain, double theta, double *recrates, double *lamda, double *logweight0array){ */
/*   double expectedcoal = 0; */
/*   double expectedsum = 0; */
/*   long i,numsamples, refchain, chaintype; */
/*   treerec *tr; */
/*   double *weight, max = 0; */

/*   chaintype = TYPE_CHAIN(chain); */
/*   refchain = REF_CHAIN(chain); */
/*   numsamples = op->numout[chaintype]; */

/*   weight = (double *)calloc(numsamples,sizeof(double)); */

  
/*   for (i = 0; i < numsamples; i++){ */
/*     tr = (&sum[0][refchain][i]); */
/*     weight[i] = getlogthetaweight(theta, tr) - logweight0array[i]; */
/*     for (j = 0; j < seq_length; j++) weight0array[i] += log(recrates[givenrecarray[j]]) * tr->recs[j] - recrates[giverecrarray[j]] * tr->sumsites[i]; */
    
/*   } */
  
/*   for (i = 0; i < numsamples; i++){ */
/*     if (max < weight[i]) max = weight[i]; */
/*   } */

/*   for (i = 0; i < numsamples; i++){ */
/*     tr = (&sum[0][refchain][i]); */
/*     expectedcoal +=  tr->numcoals * exp(weight[i] - max); */
/*     expectedsum += tr->tk  * exp(weight[i] - max); */
/*   } */
/*   free(weight); */
/*   return expectedsum /  expectedcoal; */
/* } */

void drawlikelihoodcurve(option_struct *op, long chain, double temp_theta, double *temp_recrates, double *temp_lamda, double *weight0array){
  double likelihood;
  double temp_hotrec, temp_hotlamda;
  long i,j;
  FILE *likesurface;

  likesurface = fopen("surf","w");
  temp_hotrec = temp_recrates[1];
  temp_hotlamda = temp_lamda[1];

  for (i = 1; i < 500; i++){
    for (j = 1; j < 200; j++){
      temp_recrates[1] = 2.5 + i * 0.01;
      temp_lamda[1] = 0.8 + j * 0.001;
      likelihood = getlogrelativelikelihood(op,chain,temp_theta,temp_recrates,temp_lamda,weight0array);
      printf("%lf\t %lf\t %lf\n",temp_recrates[1],temp_lamda[1],likelihood);
      fprintf(likesurface,"%lf\t%lf\t%lf\n",temp_recrates[1],temp_lamda[1],likelihood);
    }
  }
  return;
}


void hotstateprob(option_struct *op, long chain, double theta, double *recrates, double *lamda, double *weight0array){
  double *rlikelihood0, *rlikelihood1, probtheta;
  long i,pos,numsamples, refchain, chaintype;
  treerec *tr;
  FILE *hotprob;

  hotprob = fopen("hotprob","w");

  chaintype = TYPE_CHAIN(chain);
  refchain = REF_CHAIN(chain);
  numsamples = op->numout[chaintype];

  rlikelihood0 = (double *)calloc(seq_length,sizeof(double));
  rlikelihood1 = (double *)calloc(seq_length,sizeof(double));

  for (pos = 0; pos < seq_length; pos++){
    rlikelihood0[pos] = rlikelihood1[pos] = 0;
  }

  for (i = 0; i < numsamples; i++){
    tr = (&sum[0][refchain][i]);
    getalphas(tr,recrates,lamda);
    getbetas(tr,recrates,lamda);
    probtheta = getlogthetaweight(theta,tr);
    //rlikelihood = exp(probtheta + getloghmmlikelihood(recrates, lamda, tr) - weight0array[i]);
    for (pos = 0; pos < seq_length; pos++){
      rlikelihood1[pos] += exp(alpha1[pos] + beta1[pos] + probtheta - weight0array[i]);
      rlikelihood0[pos] += exp(alpha0[pos] + beta0[pos] + probtheta - weight0array[i]);
    }
  }
  for (pos = 0; pos<seq_length;pos++){
    fprintf(hotprob,"%lf\n",rlikelihood1[pos]/(rlikelihood1[pos] + rlikelihood0[pos]));
  }
  free(rlikelihood0);
  free(rlikelihood1);
  fclose(hotprob);
	    
  return;
}


/* double biserctfindtheta(option_struct *op, data_fmt *data, long chain, double start, double end, double donelike, double theta, double *recrates, double *lamda, double *weight0array){ //donelike - likelihood at CI point */
  
/*   double startlike, endlike, mid, midlike; */
/*   double *dummytheta, *dummyrecrates, *dummylamda; */
  
/*   dummytheta = (double)calloc(1,sizeof(double)); */
/*   dummyrecrates = (double)calloc(2,sizeof(double)); */
/*   dummylamda = (double)calloc(2,sizeof(double)); */

/*   dummytheta[0] = thetea; */
/*   dummyrecrates[0] = recrates[0]; */
/*   dummyrecrates[1] = recrates[1]; */
/*   dummylamda[0] = lamda[0]; */
/*   dummylamda[1] = lamda[1]; */
  
/*   mid = (start + end) * 0.5; */
  
/*   midlike = getmax(op,data,chain,dummytheta, dummyrecrates, dummylamda, weight0array,0); */


/*   //newlikelihood = getlogrelativelikelihood(op,chain,temptheta,temprecrates,templamda,weight0array); */
/*   while(!(donelike-epsilon1 < midlike && midlike < donelike+epsilon1)) { */
/*     if (donelike > midlike) {  */
/*       start = mid; */
/*       startlike = midlike; */
/*     } else { */
/*       end = mid; */
/*       endlike = midlike; */
/*     } */
/*     mid = (start+end)/2.0; */
/*     midlike = getmax(op,data,chain,dummytheta, dummyrecrates, dummylamda, weight0array,1); */
/*   }  */
  
/*   return mid; */
/* } */
/* double biserctfindrec0(option_struct *op, data_fmt *data, long chain, double start, double end, double donelike, double theta, double *recrates, double *lamda, double *weight0array){ //donelike - likelihood at CI point */
  
/*   double startlike, endlike, mid, midlike, *dummytheta, *dummyrecrates, *dummylamda; */
  
/*   dummytheta = (double)calloc(1,sizeof(double)); */
/*   dummyrecrates = (double)calloc(2,sizeof(double)); */
/*   dummylamda = (double)calloc(2,sizeof(double)); */

/*   dummytheta[0] = thetea; */
/*   dummyrecrates[0] = recrates[0]; */
/*   dummyrecrates[1] = recrates[1]; */
/*   dummylamda[0] = lamda[0]; */
/*   dummylamda[1] = lamda[1]; */
  
/*   mid = (start + end) * 0.5; */
  
/*   midlike = getmax(op,data,chain,dummytheta, dummyrecrates, dummylamda, weight0array,1); */


/*   //newlikelihood = getlogrelativelikelihood(op,chain,temptheta,temprecrates,templamda,weight0array); */
/*   while(!(donelike-epsilon1 < midlike && midlike < donelike+epsilon1)) { */
/*     if (donelike > midlike) {  */
/*       start = mid; */
/*       startlike = midlike; */
/*     } else { */
/*       end = mid; */
/*       endlike = midlike; */
/*     } */
/*     mid = (start+end)/2.0; */
/*     midlike = getmax(op,data,chain,dummytheta, dummyrecrates, dummylamda, weight0array,1); */
/*   }  */
  
/*   return mid; */
/* } */
/* double biserctfindrec1(option_struct *op, data_fmt *data, long chain, double start, double end, double donelike, double theta, double *recrates, double *lamda, double *weight0array){ //donelike - likelihood at CI point */
  
/*   double startlike, endlike, mid, midlike, *dummytheta, *dummyrecrates, *dummylamda; */
  
/*   dummytheta = (double)calloc(1,sizeof(double)); */
/*   dummyrecrates = (double)calloc(2,sizeof(double)); */
/*   dummylamda = (double)calloc(2,sizeof(double)); */

/*   dummytheta[0] = thetea; */
/*   dummyrecrates[0] = recrates[0]; */
/*   dummyrecrates[1] = recrates[1]; */
/*   dummylamda[0] = lamda[0]; */
/*   dummylamda[1] = lamda[1]; */
  
/*   mid = (start + end) * 0.5; */
  
/*   midlike = getmax(op,data,chain,dummytheta, dummyrecrates, dummylamda, weight0array,1); */


/*   //newlikelihood = getlogrelativelikelihood(op,chain,temptheta,temprecrates,templamda,weight0array); */
/*   while(!(donelike-epsilon1 < midlike && midlike < donelike+epsilon1)) { */
/*     if (donelike > midlike) {  */
/*       start = mid; */
/*       startlike = midlike; */
/*     } else { */
/*       end = mid; */
/*       endlike = midlike; */
/*     } */
/*     mid = (start+end)/2.0; */
/*     midlike = getmax(op,data,chain,dummytheta, dummyrecrates, dummylamda, weight0array,1); */
/*   }  */
  
/*   return mid; */
/* } */

/* double biserctfindlamda0(option_struct *op, data_fmt *data, long chain, double start, double end, double donelike, double theta, double *recrates, double *lamda, double *weight0array){ //donelike - likelihood at CI point */
  
/*   double startlike, endlike, mid, midlike, *dummytheta, *dummyrecrates, *dummylamda; */
  
/*   dummytheta = (double)calloc(1,sizeof(double)); */
/*   dummyrecrates = (double)calloc(2,sizeof(double)); */
/*   dummylamda = (double)calloc(2,sizeof(double)); */

/*   dummytheta[0] = thetea; */
/*   dummyrecrates[0] = recrates[0]; */
/*   dummyrecrates[1] = recrates[1]; */
/*   dummylamda[0] = lamda[0]; */
/*   dummylamda[1] = lamda[1]; */
  
/*   mid = (start + end) * 0.5; */
  
/*   midlike = getmax(op,data,chain,dummytheta, dummyrecrates, dummylamda, weight0array,1); */


/*   //newlikelihood = getlogrelativelikelihood(op,chain,temptheta,temprecrates,templamda,weight0array); */
/*   while(!(donelike-epsilon1 < midlike && midlike < donelike+epsilon1)) { */
/*     if (donelike > midlike) {  */
/*       start = mid; */
/*       startlike = midlike; */
/*     } else { */
/*       end = mid; */
/*       endlike = midlike; */
/*     } */
/*     mid = (start+end)/2.0; */
/*     midlike = getmax(op,data,chain,dummytheta, dummyrecrates, dummylamda, weight0array,1); */
/*   }  */
  
/*   return mid; */
/* } */
/* double biserctfindlamda1(option_struct *op, data_fmt *data, long chain, double start, double end, double donelike, double theta, double *recrates, double *lamda, double *weight0array){ //donelike - likelihood at CI point */
  
/*   double startlike, endlike, mid, midlike, *dummytheta, *dummyrecrates, *dummylamda; */
  
/*   dummytheta = (double)calloc(1,sizeof(double)); */
/*   dummyrecrates = (double)calloc(2,sizeof(double)); */
/*   dummylamda = (double)calloc(2,sizeof(double)); */

/*   dummytheta[0] = thetea; */
/*   dummyrecrates[0] = recrates[0]; */
/*   dummyrecrates[1] = recrates[1]; */
/*   dummylamda[0] = lamda[0]; */
/*   dummylamda[1] = lamda[1]; */
  
/*   mid = (start + end) * 0.5; */
  
/*   midlike = getmax(op,data,chain,dummytheta, dummyrecrates, dummylamda, weight0array,1); */


/*   //newlikelihood = getlogrelativelikelihood(op,chain,temptheta,temprecrates,templamda,weight0array); */
/*   while(!(donelike-epsilon1 < midlike && midlike < donelike+epsilon1)) { */
/*     if (donelike > midlike) {  */
/*       start = mid; */
/*       startlike = midlike; */
/*     } else { */
/*       end = mid; */
/*       endlike = midlike; */
/*     } */
/*     mid = (start+end)/2.0; */
/*     midlike = getmax(op,data,chain,dummytheta, dummyrecrates, dummylamda, weight0array,1); */
/*   }  */
  
/*   return mid; */
/* } */

/* double getmax(option_struct *op, data_fmt *data, long chain, double *temptheta, double *temprecrates, double *templamda, double *weight0array, int option){ */
/*   double oldlikelihood, newlikelihood; */

/*   oldlikelihood = getlogrelativelikelihood(op,chain,temptheta,temprecrates,templamda,weight0array); */
/*   while (1){ */
/*     if (option != 0) temptheta[0] = gettheta(op,chain,temptheta[0],temprecrates,templamda,weight0array);  */
/*     if (option != 1) temprecrates[0] = getrecrate(op,chain,temptheta[0],temprecrates,templamda,weight0array,0); */
/*     if (option != 2) temprecrates[1] = getrecrate(op,chain,temptheta[0],temprecrates,templamda,weight0array,1); */
/*     if (option != 3) templamda[0] = getlamda(op,chain,temptheta[0],temprecrates,templamda,weight0array,0); */
/*     if (option != 4) templamda[1] = getlamda(op,chain,temptheta[0],temprecrates,templamda,weight0array,1); */
  
/*     newlikelihood = getlogrelativelikelihood(op,chain,temptheta[0],temprecrates,templamda,weight0array); */
/*     if ((oldlikelihood / newlikelihood) > 0.99) break; */
/*   } */

/*   return newlikelihood; */
/* } */

void likesurface_rc_rh(option_struct *op, long chain, double *weight0array){
  double Rtotal, Rcold,Rhot,newlikelihood,lamda0, lamda1,recrates0,recrates1,tmt;
  long i,j,k;
  FILE *surface, *s2;

  surface=fopen("rcrh","w");
  s2 = fopen("lclh","w");
  Rtotal = ((1 - templamda[1]) * temprecrates[0] + (1 - templamda[0]) * temprecrates[1]) / (2 - templamda[0] - templamda[1]);

  tmt = temptheta;
  lamda0 = templamda[0];
  lamda1 = templamda[1];
  recrates0 = temprecrates[0];
  recrates1 = temprecrates[1];

  for (i = 0; i < 20; i++){
    templamda[0] = lamda0 - ((1.0 - lamda0) * 3.0)  + (1.0 - lamda0) / 5.0 * i;
    for (j = 0; j < 20; j++){
      templamda[1] = lamda1 - ((1.0 - lamda1) * 3.0)  + (1.0 - lamda1) / 5.0 * j;
      if (temprecrates[0] == 0) temprecrates[0] = 0.01;
      temprecrates[0] = recrates0;
      temprecrates[1] = recrates1;
      for (k = 0; k < 5; k ++){ 
	temptheta = gettheta(op,chain, temptheta, temprecrates,templamda,weight0array);
	temprecrates[0] = getrecrate(op,chain,temptheta,temprecrates,templamda,weight0array,0);
	temprecrates[1] = getrecrate(op,chain,temptheta,temprecrates,templamda,weight0array,1);
      }
      
      Rcold =  (1 - templamda[1]) * temprecrates[0] / (2 - templamda[0] - templamda[1]);
      Rhot =  (1 - templamda[0]) * temprecrates[1] / (2 - templamda[0] - templamda[1]);

      newlikelihood = getlogrelativelikelihood(op,chain,temptheta,temprecrates,templamda,weight0array);

      fprintf(surface,"%lf\t%lf\t%lf\n",Rcold,Rhot,newlikelihood);
      fprintf(s2, "%lf\t%lf\t%lf\n",templamda[0],templamda[1],newlikelihood);
      printf("%lf, %lf  %lf,%lf %lf\n",templamda[0],temprecrates[0],templamda[1],temprecrates[1],newlikelihood);
    }
  }
  templamda[0] = lamda0;
  templamda[1] = lamda1;
  temprecrates[0] = recrates0;
  temprecrates[1] = recrates1;

  fclose(s2);
  fclose(surface);
  return;
}
  
void likesurface_theta(option_struct *op, long chain, double *weight0array){
  double Rtotal,newlikelihood,tmp,lamda0, lamda1,recrates0,recrates1;
  long i,k;
  FILE *surface;

  surface=fopen("thetar","w");

  Rtotal = ((1 - templamda[1]) * temprecrates[0] + (1 - templamda[0]) * temprecrates[1]) / (2 - templamda[0] - templamda[1]);
  tmp = temptheta;
  lamda0 = templamda[0];
  lamda1 = templamda[1];
  recrates0 = temprecrates[0];
  recrates1 = temprecrates[1];


  for (i = 0; i < 100; i++){
    temptheta = 0.001 * i;
    templamda[0] = lamda0;
    templamda[1] = lamda1;
    temprecrates[0] = recrates0;
    temprecrates[1] = recrates1;
    for (k = 0; k < 10; k ++){
      temprecrates[0] = getrecrate(op,chain,temptheta,temprecrates,templamda,weight0array,0);
      if (temprecrates[0] == 0) temprecrates[0] = 0.000001;
      temprecrates[1] = getrecrate(op,chain,temptheta,temprecrates,templamda,weight0array,1);
      if (temprecrates[1] == 0) temprecrates[1] = 0.000001;
      templamda[0] = getlamda(op,chain,temptheta,temprecrates,templamda,weight0array,0);
      templamda[1] = getlamda(op,chain,temptheta,temprecrates,templamda,weight0array,1);
    }
  
    Rtotal = ((1 - templamda[1]) * temprecrates[0] + (1 - templamda[0]) * temprecrates[1]) / (2 - templamda[0] - templamda[1]);    
    newlikelihood = getlogrelativelikelihood(op,chain,temptheta,temprecrates,templamda,weight0array);
    
    fprintf(surface,"%lf\t%lf\t%lf\n",temptheta, Rtotal,newlikelihood);
    printf("%lf, %lf, %lf\n",temptheta,Rtotal,newlikelihood);
    
  }
  temptheta = tmp;
  fclose(surface);
  return;
}

void likesurface_rec(option_struct *op, long chain, double *weight0array){
  double Rtotal,newlikelihood,rates0,rates1,tmp,lamda0,lamda1;
  long i,j,k;
  FILE *surface;

  surface=fopen("recch","w");
  Rtotal = ((1 - templamda[1]) * temprecrates[0] + (1 - templamda[0]) * temprecrates[1]) / (2 - templamda[0] - templamda[1]);
  tmp = temptheta;
  lamda0 = templamda[0];
  lamda1 = templamda[1];
  rates0 = temprecrates[0];
  rates1 = temprecrates[1];
  
  for (i = 0; i < 20; i++){
    temprecrates[0] = rates0 * 0.75 + 0.025 * rates0 * i;
    for (j = 0; j < 20; j++){
      temprecrates[1] = rates1 * 0.75 + 0.025 * rates0 * j;
      temptheta = tmp;
      templamda[0] = lamda0;
      templamda[1] = lamda1;
      for (k = 0; k < 5; k ++){
	temptheta = gettheta(op,chain, temptheta, temprecrates,templamda,weight0array);
	templamda[0] = getlamda(op,chain,temptheta,temprecrates,templamda,weight0array,0);
	templamda[1] = getlamda(op,chain,temptheta,temprecrates,templamda,weight0array,1);
      }
      

      newlikelihood = getlogrelativelikelihood(op,chain,temptheta,temprecrates,templamda,weight0array);

      fprintf(surface,"%lf\t%lf\t%lf\n",temprecrates[0],temprecrates[1],newlikelihood);
    }
  }
  fclose(surface);
  return;
}

void likesurface_lamda1(option_struct *op, long chain, double *weight0array){
  double Rtotal,newlikelihood,rates0,rates1,tmp,lamda0,lamda1;
  long i,k;
  FILE *surface;

  surface=fopen("lamda1","w");
  Rtotal = ((1 - templamda[1]) * temprecrates[0] + (1 - templamda[0]) * temprecrates[1]) / (2 - templamda[0] - templamda[1]);
  tmp = temptheta;
  lamda0 = templamda[0];
  lamda1 = templamda[1];
  rates0 = temprecrates[0];
  rates1 = temprecrates[1];
  
  for (i = 1; i < 10; i++){
    templamda[1] = 0.1 * i;
    temptheta = tmp;
    templamda[0] = lamda0;
    temprecrates[0] = rates0;
    temprecrates[1] = rates1;
    for (k = 0; k < 10; k ++){
      temptheta = gettheta(op,chain, temptheta, temprecrates,templamda,weight0array);
      temprecrates[0] = getrecrate(op,chain,temptheta,temprecrates,templamda,weight0array,0);
      temprecrates[1] = getrecrate(op,chain,temptheta,temprecrates,templamda,weight0array,1);
      templamda[0] = getlamda(op,chain,temptheta,temprecrates,templamda,weight0array,0);
    }
    newlikelihood = getlogrelativelikelihood(op,chain,temptheta,temprecrates,templamda,weight0array);
    
    fprintf(surface,"%lf\t%lf\n",templamda[1],newlikelihood);
  }
   for (i = 0; i < 20; i++){
    templamda[1] = 0.8 + 0.01 * i;
    temptheta = tmp;
    templamda[0] = lamda0;
    temprecrates[0] = rates0;
    temprecrates[1] = rates1;
    for (k = 0; k < 10; k ++){
      temptheta = gettheta(op,chain, temptheta, temprecrates,templamda,weight0array);
      temprecrates[0] = getrecrate(op,chain,temptheta,temprecrates,templamda,weight0array,0);
      temprecrates[1] = getrecrate(op,chain,temptheta,temprecrates,templamda,weight0array,1);
      templamda[0] = getlamda(op,chain,temptheta,temprecrates,templamda,weight0array,0);
    }
    newlikelihood = getlogrelativelikelihood(op,chain,temptheta,temprecrates,templamda,weight0array);
    
    fprintf(surface,"%lf\t%lf\n",templamda[1],newlikelihood);
  } 

   for (i = 1; i < 10; i++){
    templamda[1] = 0.99 + 0.001 * i;
    temptheta = tmp;
    templamda[0] = lamda0;
    temprecrates[0] = rates0;
    temprecrates[1] = rates1;
    for (k = 0; k < 10; k ++){
      temptheta = gettheta(op,chain, temptheta, temprecrates,templamda,weight0array);
      temprecrates[0] = getrecrate(op,chain,temptheta,temprecrates,templamda,weight0array,0);
      temprecrates[1] = getrecrate(op,chain,temptheta,temprecrates,templamda,weight0array,1);
      templamda[0] = getlamda(op,chain,temptheta,temprecrates,templamda,weight0array,0);
    }
    newlikelihood = getlogrelativelikelihood(op,chain,temptheta,temprecrates,templamda,weight0array);
    
    fprintf(surface,"%lf\t%lf\n",templamda[1],newlikelihood);
  } 
  fclose(surface);
  return;
}

void likesurface_rec1(option_struct *op, long chain, double *weight0array){
  double Rtotal,newlikelihood,rates0,rates1,tmp,lamda0,lamda1;
  long i,k;
  FILE *surface;

  surface=fopen("rec1","w");
  Rtotal = ((1 - templamda[1]) * temprecrates[0] + (1 - templamda[0]) * temprecrates[1]) / (2 - templamda[0] - templamda[1]);
  tmp = temptheta;
  lamda0 = templamda[0];
  lamda1 = templamda[1];
  rates0 = temprecrates[0];
  rates1 = temprecrates[1];
  
  for (i = 1; i < 20; i++){
    temprecrates[1] = 1 * i;
    temptheta = tmp;
    templamda[0] = lamda0;
    temprecrates[0] = rates0;
    templamda[1] = lamda1;
    for (k = 0; k < 10; k ++){
      temptheta = gettheta(op,chain, temptheta, temprecrates,templamda,weight0array);
      temprecrates[0] = getrecrate(op,chain,temptheta,temprecrates,templamda,weight0array,0);
      if (temprecrates[0] > temprecrates[1] ) temprecrates[0] = temprecrates[1];
      templamda[1] = getlamda(op,chain,temptheta,temprecrates,templamda,weight0array,1);
      templamda[0] = getlamda(op,chain,temptheta,temprecrates,templamda,weight0array,0);
    }
    newlikelihood = getlogrelativelikelihood(op,chain,temptheta,temprecrates,templamda,weight0array);
    
    fprintf(surface,"%lf\t%lf\n",temprecrates[1],newlikelihood);
  }
   for (i = 1; i < 20; i++){
    temprecrates[1] = 2.5 + 0.25 * i;
    temptheta = tmp;
    templamda[0] = lamda0;
    temprecrates[0] = rates0;
    templamda[1] = lamda1;
    for (k = 0; k < 10; k ++){
      temptheta = gettheta(op,chain, temptheta, temprecrates,templamda,weight0array);
      temprecrates[0] = getrecrate(op,chain,temptheta,temprecrates,templamda,weight0array,0);
      templamda[1] = getlamda(op,chain,temptheta,temprecrates,templamda,weight0array,1);
      templamda[0] = getlamda(op,chain,temptheta,temprecrates,templamda,weight0array,0);
    }
    newlikelihood = getlogrelativelikelihood(op,chain,temptheta,temprecrates,templamda,weight0array);
    
    fprintf(surface,"%lf\t%lf\n",temprecrates[1],newlikelihood);
  } 
   
  fclose(surface);
  return;
}
