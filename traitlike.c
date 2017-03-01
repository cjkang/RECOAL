#ifndef TRAITLIKE
#include "traitlike.h"
#endif

/* This file contains the data likelihood for a disease
trait model.  We assume we know mutrait, the ratio of the
forward and back mutation rates for the trait locus; and
traitratio, the ratio of the forward mutation rate at the
trait locus to the single-nucleotide mutation rate.  The
input tree is scaled to the single-nucleotide rate. 
We also assume we know pd, the population prevalence of
the disease.  Traitarray is an array of doubles which
correspond to the possible locations of the trait.  */

/* Each top nodelet must provide an array double z[2] for the
disease trait likelihoods.  At the tips these will contain
either [1,0] healthy or [0,1] diseased.  */

#define HEALTHY 1
#define DISEASED 0

extern boolean istip(node *p);
extern boolean iscoal(node *p);
extern boolean inrange(long *ranges, long location);
extern void finddnasubtrees(option_struct *op, data_fmt *data, tlist *t,
   long *ranges);
extern void printhistplot(option_struct *op, data_fmt *data,
   long *bars, long htmult, long wdmult);
extern long sitetomarker(option_struct *op, data_fmt *data, long site);
extern long countsites(option_struct *op, data_fmt *data);

extern FILE *outfile;
extern long population, locus;

void traitsiteplot(option_struct *op, data_fmt *data, double *llike)
{
long numsites, site, *traitllike, wdmult, htmult;
double mostneg, mostpos;

numsites = countsites(op,data);

for(site = 0, mostneg = POSMAX, mostpos = NEGMAX; site < numsites; site++) {
   mostneg = ((llike[site] < mostneg) ? llike[site] : mostneg);
   mostpos = ((llike[site] > mostpos) ? llike[site] : mostpos);
}

mostpos -= mostneg;
for(site = 0; site < numsites; site++) llike[site] -= mostneg;

for(wdmult = 1; mostpos/wdmult > HIST_MAXWD; wdmult++) ;
for(htmult = 1; numsites/htmult > HIST_MAXHT; htmult++) ;

traitllike = (long *)calloc(numsites,sizeof(long));
for(site = 0; site < numsites; site++)
   traitllike[site] = ((llike[site] - (long)llike[site] >= 0.5) ?
                (long)(llike[site]+1.0) : (long)llike[site]);

/* print the histogram header */
fprintf(outfile,"\n\nHistogram of Ln(likelihood) along sequence\n\n");
fprintf(outfile,"Scaling: x-axis = %ld (likelihood units per tick)\n",wdmult);
fprintf(outfile,"         y-axis = %ld (sites per tick)\n",htmult);
fprintf(outfile,"         x-axis start at %f\n",mostneg);

printhistplot(op,data,traitllike,htmult,wdmult);

free(traitllike);

} /* traitsiteplot */
   

void copytraits(node *p, node *q)
{
  long i;

  for (i = 0; i < NUMTRAIT; i++) q->z[i] = p->z[i];
 
} /* copytraits */


double traiteval (option_struct *op, data_fmt *data, tree *tr,
  double mu, double nu, double pd, long location)
{
  double likel;

  traitsmooth(op, data, tr, mu, nu, location);

  likel = tr->root->back->z[0] * pd +
          tr->root->back->z[1] * (1 - pd);

/* NO NO NO!  likel = log(likel); */

  return(likel); 

} /* traiteval */

void traitsmooth (option_struct *op, data_fmt *data, tree *tr,
  double mu, double nu, long location)
{
tlist *t;

for(t = tr->tymelist->succ; t != NULL; t = t->succ)
   if (iscoal(t->eventnode))
      traitview(op,data,t->eventnode,mu,nu,location);

#if 0
  if (istip(p)) return;
  if (!p->next->top && inrange(p->next->back->ranges,location))
     traitsmooth (op, data, p->next->back, mu, nu, location);
  if (!p->next->next->top && inrange(p->next->next->back->ranges,location)) 
     traitsmooth(op, data, p->next->next->back, mu, nu, location);

  if (iscoal(p))
    traitview(op, data, p, mu, nu, location);
#endif

} /* traitsmooth */

double probHH(double tyme, double mu, double nu) 
{
  return(1.0 - probHD(tyme, mu, nu));
} /* probHH */

double probHD(double tyme, double mu, double nu)
{
  double prob;

  prob = (mu/(mu+nu))*(1-exp((-mu-nu)*tyme));
  return(prob);

} /* probHD */

double probDH(double tyme, double mu, double nu)
{
  double prob;

  prob = (nu/(mu+nu))*(1-exp((-mu-nu)*tyme));
  return(prob);

} /* probDH */

double probDD(double tyme, double mu, double nu)
{
  return(1.0 - probDH(tyme, mu, nu));
} /* probDD */

void traitview (option_struct *op, data_fmt *data, node *p, double mu,
  double nu, long location)
{
  node *q, *r;
  double qtyme = 0.0, rtyme = 0.0;

  if (inrange(p->next->back->ranges,location)) {
     q = findcoal(op,data,p->next,&qtyme);
     qtyme = vtol(op,data,qtyme);
  } else q = NULL;
  if (inrange(p->next->next->back->ranges,location)) {
     r = findcoal(op,data,p->next->next,&rtyme);
     rtyme = vtol(op,data,rtyme);
  } else r = NULL;
  
  p->z[0] = 1.0;
  p->z[1] = 1.0;

  if (q) {
     p->z[0] = q->z[0] * probHH(qtyme,mu,nu)
             + q->z[1] * probDH(qtyme,mu,nu);
     p->z[1] = q->z[0] * probHD(qtyme,mu,nu)
             + q->z[1] * probDD(qtyme,mu,nu);
  }
  if (r) {
     p->z[0] *= r->z[0] * probHH(rtyme,mu,nu)
             + r->z[1] * probDH(rtyme,mu,nu);
     p->z[1] *= r->z[0] * probHD(rtyme,mu,nu)
             + r->z[1] * probDD(rtyme,mu,nu);
  }

#if 0
  q = p->next->back;
  r = p->next->next->back;

  if (q->top) {
    if (inrange(q->ranges,location)) {
       qtyme = p->tyme - q->tyme;
       p->z[0] = q->z[0] * probHH(qtyme,mu,nu)
               + q->z[1] * probDH(qtyme,mu,nu);
       p->z[1] = q->z[0] * probHD(qtyme,mu,nu)
               + q->z[1] * probDD(qtyme,mu,nu);
    } 
  }
  if (r->top) {
    if (inrange(r->ranges,location)) {
       rtyme = p->tyme - r->tyme;
       p->z[0] *= r->z[0] * probHH(rtyme,mu,nu)
               + r->z[1] * probDH(rtyme,mu,nu);
       p->z[1] *= r->z[0] * probHD(rtyme,mu,nu)
               + r->z[1] * probDD(rtyme,mu,nu);
    }
  }
#endif

} /* traitview */

void traitlike(option_struct *op, data_fmt *data, tree *tr, long numsites, 
  double mutrait, double traitratio, double pd, double *traitarray)
{
  long location, subtree, *subtree_ranges;
  double mu, nu, like;

  mu = traitratio;
  nu = traitratio/mutrait;

  subtree_ranges = 
     (long *)calloc(2*tr->numrecombs+4,sizeof(long));
  subtree_ranges[0] = 1;

  finddnasubtrees(op,data,tr->tymelist,subtree_ranges);

  for (subtree = 0; subtree < subtree_ranges[0]; subtree++) {
      like = traiteval(op,data,tr,mu,nu,pd,subtree_ranges[2*subtree+1]);
      for(location = subtree_ranges[2*subtree+1];
          location <= subtree_ranges[2*subtree+2]; location++)
         traitarray[location] += like;
  }

  free(subtree_ranges);

} /* traitlike */

void traitprint(long numsites, double *traitarray, FILE *out, long numout)
{
  long i;

  fprintf(out,"Trait ln-likelihoods by site:\n");
  for (i = 0; i < numsites; i++) {
    fprintf(out,"%f ",log(traitarray[i]/(double)numout));
    if (((i+1)/10)*10==i+1) fprintf(out,"\n");
  }
  fprintf(out,"\n");
} /* traitprint */

void traitread(tree *tr, long numseq)
{
  long i, j, n;
  char str[NMLNGTH], gch;
  boolean found;
  FILE *traitfile;

  traitfile = fopen("traitfile","r");

  for(j = 0; j < numseq; j++) {
    /* read a name */
    gch = getc(traitfile);
    if (gch == '\n')
      gch = ' ';
    for (i = 0; i < NMLNGTH; i++)
      str[i] = ' ';
    n = 0;
    do {
      if (gch == '_')
        gch = ' ';
      str[n] = gch;
      if (eoln(traitfile)) {
        fscanf(traitfile, "%*[^\n]");
        getc(traitfile);
      }
      gch = getc(traitfile);
      if (gch == '\n')
        gch = ' ';
      n++;
    } while (gch != ':' && gch != ',' && gch != ')' && n < NMLNGTH);
    n = 1;
    do {
      found = true;
      for (i = 0; i < NMLNGTH; i++)
        found = (found && str[i] == tr->nodep[n]->nayme[i]);
      if (!found)
        n++;
    } while (!(n > numseq || found));
    if (n > numseq) {
      printf("In TRAIT DATA:  Cannot find sequence: ");
      for (i = 0; i < NMLNGTH; i++)
        putchar(str[i]);
      putchar('\n');
      return;
    }
    if (gch == ' ') do {
      gch = getc(traitfile);
    } while (gch == ' ' || gch == '\n');
    if (gch == 'H') {
      tr->nodep[n]->z[0] = 1.0;
      tr->nodep[n]->z[1] = 0.0;
    } else {
      tr->nodep[n]->z[0] = 0.0;
      tr->nodep[n]->z[1] = 1.0;
    }
    do {
      gch = getc(traitfile);
    } while (gch != '\n');
  }
} /* traitread */


void traitresult(option_struct *op, data_fmt *data, double *traitarray,
   FILE *out, long numout)
/* evaluate how many sites were preferred, and whether
the truth was among them */
{
long i, winners, truth, lineout, numsites;
double biggest;
FILE *truetrait, *mathout;
dnadata *dna;

dna = data->dnaptr;

numsites = countsites(op,data);

truetrait = fopen("truetrait","r");
fscanf(truetrait,"%ld",&truth);

biggest = NEGMAX;
winners = 0;
for (i = 0; i < numsites; i++) {
  traitarray[i]=log(traitarray[i]) - log((double)numout);
  if (traitarray[i] > biggest) biggest = traitarray[i];
}
for (i = 0; i < numsites; i++) {
  if (traitarray[i] >= biggest-DF2) {
     winners++; 
  }
}
fprintf(out,"Mapper report:\n");
fprintf(out,"Number of sites within confidence interval: %ld\n",
  winners);
if (traitarray[truth] >= biggest-DF2)
  fprintf(out,"True site (%ld) included.\n",truth);
else fprintf(out,"True site (%ld) excluded.\n",truth);
fprintf(out,"\n");

mathout = fopen("mathout","w+");
fprintf(mathout,"traitlike = {\n");
lineout = 5;
for (i = 0; i < numsites; i++) {
   fprintf(mathout,"%f",traitarray[i]);
   if (i != numsites-1) fprintf(mathout,", ");
   if (i == lineout) {fprintf(out,"\n"); lineout += 5;}
}
fprintf(mathout,"}\n");
fclose(mathout);

fclose(truetrait);

} /* traitresult */

