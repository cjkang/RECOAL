#define REC_MODELLIKE_INCLUDE

#ifndef RECOMBINE_INCLUDE
#include "recombine.h"
#endif

#ifndef TREEREC_H
#include "cutpointrec.h"
#endif




double lastintervalcalc(option_struct *op, data_fmt *data, long chain);
void intervalcheck(option_struct *op, data_fmt *data, long chain);

long whatsign(double value);
boolean zerocheck(double value);
double makepositive(double value, double change);
double changeparamLN(double value, double change);
void model_alloc(option_struct *op, data_fmt *data);
void scoretree(option_struct *op, data_fmt *data, long chain);
double count_active_tlist(option_struct *op, data_fmt *data, tlist *t);

void rec_scoreprint(option_struct *op, data_fmt *data, long lowcus,
   long chain, boolean mathematica);
void rec_readscoretree(FILE *numcoals, FILE *numrecombs, FILE *kks,
   FILE *ends, FILE *actives, option_struct *op,
   long lowcus, long chain, long trii);
void rec_scoreread(option_struct *op, data_fmt *data);
void rec_likeprint(option_struct *op, data_fmt *data, long lowcus,
   long chain);
void rec_likeprint2(option_struct *op, data_fmt *data, double thmax,
   long lowcus, long chain);

double stuff(double theta, double rec, treerec *tr);
double rectreellike(option_struct *op, double theta, double rec,
   double lth, double lrec, long lowcus, long chain, long trii);
double recchainllike(option_struct *op, double theta, double rec, 
   long lowcus, long chain);
double rec_locusllike(option_struct *op, data_fmt *data,
   double theta, double rec);
double model_likelihood(option_struct *op, data_fmt *data,
   double theta, double rec, long lowcus, long chain);
double frec_chainllike(option_struct *op, double theta, double rec,
   double **denom, long lowcus, long chain);
double frec_locusllike(option_struct *op, data_fmt *data,
   double theta, double rec, double **denom);
double fmodel_likelihood(option_struct *op, data_fmt *data,
   double theta, double rec, double **denom, long lowcus, long chain);

double dstuff(double theta, treerec *tr);
double rec_thetalderiv(option_struct *op, double theta, double rec,
   double **llike0, long lowcus, long chain, double *dfn, long *dfnplus);
double rec_reclderiv(option_struct *op, double theta, double rec,
   double **llike0, long lowcus, long chain, double *dfn, long *dfnplus);
double rec_locus_thetalderiv(option_struct *op, data_fmt *data,
   double theta, double rec, double **llike0,double *dfx, long *dfxplus);
double rec_locus_reclderiv(option_struct *op, data_fmt *data, double theta,
   double rec, double **llike0, double *dfx, long *dfxplus);

double NRrec_thetalderiv(option_struct *op, double theta, double rec,
   double **llike0, long lowcus, long chain, double *dfn, double *ddfn,
   long *dfnplus, long *ddfnplus);
double NRrec_reclderiv(option_struct *op, double theta, double rec,
   double **llike0, long lowcus, long chain, double *dfn, double *ddfn,
   long *dfnplus, long *ddfnplus);
double NRrec_partiallderiv(option_struct *op, double theta, double rec,
   double **llike0, long lowcus, long chain, double *ddfn, long *ddfnplus);
double NRrec_locus_thetalderiv(option_struct *op, data_fmt *data,
   double theta, double rec, double **llike0, double *fx,
   double *dfx, long *fxplus, long *dfxplus);
double NRrec_locus_reclderiv(option_struct *op, data_fmt *data, double theta,
   double rec, double **llike0, double *fx, double *dfx, long *fxplus,
   long *dfxplus);
double NRrec_locus_partiallderiv(option_struct *op, data_fmt *data,
   double theta, double rec, double **llike0, double *df, long *dfplus);


void init_matrix(option_struct *op, data_fmt *data, double *fx, double **dfx, double rec,
   double theta, double **llike0, long *fxplus, long **dfxplus,
   long lowcus, long chain, boolean thetaonly);
void to_LNparameters(double theta, double rec, double *fx, double **dfx,
 long *fxplus, long **dfxplus);
boolean matrix_denom(double **dfx, long **dfxplus, double *denom,
   long *denomsign);
boolean check_curvature(double *fx, double **dfx, long *fxplus,
   long **dfxplus);
boolean calc_change(double *fx, double **dfx, long *fxplus,
   long **dfxplus, double *answ);
boolean halfback(option_struct *op, data_fmt *data, double theta, double rec,
   double thetachange, double recchange, double **llike0, double oldlike,
   double *newtheta, double *newrec, long vflag, long lowcus, long chain,
   boolean thetaonly);
void non_nr_step(option_struct *op, data_fmt *data, double theta, double rec, 
   double **llike0, long lowcus, long chain, double *recchange,
   double *thetachange, double *pchange, long *fxplus);
boolean found_recomb(option_struct *op, data_fmt *data, long lowcus, long chain);
void matrixpoint_alloc(option_struct *op, data_fmt *data, long chain, double **fx,
   long **fxplus, double ***dfx, long ***dfxplus, double **fchange,
   long **fplus, double ***llike0);
void precalc_llike0(option_struct *op, data_fmt *data, long lowcus, long chain, 
   double **llike0);
void matrixpoint(option_struct *op, data_fmt *data, long lowcus, long chain,
   double *thetaresult, double *recresult);

void reset_hess (double **hess, long n);
void calc_loci_param(option_struct *op, data_fmt *data, double *oparam,
   double *nparam, double lambda, double *dv);
double psi(option_struct *op, data_fmt *data, double *param,
   double **denom, long lowcus, long chain, double lambda, double *dv);
double calc_line(option_struct *op, data_fmt *data, double *param,
   double **denom, double *dv, long lowcus, long chain);
double calcnorm (double *d, long size);
void calc_hessian (double **hess, long n, double *delta, double *gama);
void calc_dirv (double *dv, double **hess, double *gxv, long n);
void grad2loggrad (double *param, double *d, double PGC, long nn);
void calc_grad(option_struct *op, data_fmt *data, double *param,
   double **denom, long lowcus, long chain, double *grad);
void broydenpoint(option_struct *op, data_fmt *data, long lowcus,
   long chain, double *thetaresult, double *recresult, double **llike0);

void printhistplot(option_struct *op, data_fmt *data, long *bars,
   long htmult, long wdmult);
void rechistout(option_struct *op, data_fmt *data, cutpointrec *cutp);
double find_thetamax(option_struct *op, data_fmt *data, long lowcus,
   long chain, double startrec, double starttheta, double **llike0, 
   double *llike);
double find_recmax(option_struct *op, data_fmt *data, long lowcus,
   long chain, double startrec, double starttheta, double **llike0,
   double *llike);
double bisectfindth(option_struct *op, data_fmt *data, long lowcus,
   long chain, double start, double startlike, double end,
   double endlike, double donelike, double **llike0, double rec);
double bisectfindrec(option_struct *op, data_fmt *data, long lowcus,
   long chain, double start, double startlike, double end,
   double endlike, double donelike, double **llike0, double th);
void dobounds(option_struct *op, data_fmt *data, long lowcus,
   double bestth, double bestrec, double **llike0, double donelike,
   double maxlike, double *thlb, double *thub, double *reclb,
   double *recub);
void profile_estimate(option_struct *op, data_fmt *data, long lowcus,
   double bestth, double bestrec, double **llike0);
void rec_estimate(option_struct *op, data_fmt *data, long lowcus,
   long chain, boolean locusend);
void rec_liketable(option_struct *op, data_fmt *data, double th, double rec,
   long chain, long lowcus, boolean to_screen);
void rec_linkgraph(option_struct *op, data_fmt *data, treerec ***triis,
  long chain, long lowcus);

boolean gridcheck(option_struct *op, data_fmt *data, double test_theta,
   double test_r, long lowcus, long chain);
void runlike(option_struct *op, data_fmt *data, treerec ***triis,
   long lowcus);

double coalprob(option_struct *op, data_fmt *data, tree *tr, double 
   theta, double r);



void scoretreeH(option_struct *op, data_fmt *data, long chain);
void updaterecnum(long start, long end, double tyme,  recnumb *recnum);
void addrec(option_struct *op, data_fmt *data, tlist *t, double tyme, recnumb *recnum);
void addrecs(option_struct *op, data_fmt *data, node *p, tlist *t, double tyme, recnumb *recnum);
void updatecoal(option_struct *op, data_fmt *data, node *p,double tyme,  recnumb *recnum);
void updaterec(option_struct *op, data_fmt *data, node *p, double tyme, recnumb *recnum);
void addrecfc(option_struct *op, data_fmt *data, node *p, double tyme, recnumb *recnum);
void addrecnum(recnumb *recnum, long site);

recnumb *scorerecs(option_struct *op, data_fmt *data, tree *target);
void updaterecnumbs(node *p, recnumb *recnumbs);
void finds(node *p, long *start, long *end);
void EMestimate(option_struct *op, long chain, double *parameters);
double getthetav(option_struct *op, double *parameters, long chain);
void getrecs(option_struct *op, double *parameters, long chain);
double treellikelihood(option_struct *op, double *parameters, long chain);
double gettrweight(treerec *tr, double *parameters);
void rec_estimatev(option_struct *op, data_fmt *data, long lowcus, long chain, boolean locusend);

void addrecevent(treerec *trii);

double sumxi(double *recrates, double *lamda, treerec *tr, double *alphaarray, double *betaarray, double hmml, int state1, int state2);
double expectedrec(double *recrates, double *lamda, treerec *tr, int state);
double expectedbranchlength(double *recrates, double *lamda, treerec *tr, int state);
double expectedtransitions(double *recrates, double *lamda, treerec *tr, int state0, int state1);
void getalphaarray(double *recrates, double *lamda, treerec *tr, int state, double *alphaarray);
void getbetaarray(double *recrates, double *lamda, treerec *tr, int state, double *betaarray);
double gethmmlikelihood(double *recrates, double *lamda, treerec *tr);
double getthetaweight(double theta, treerec *tr);
double gettheta(option_struct *op, long chain, double theta, double *recrates, double *lamda, double *weight0array);
double getlamda(option_struct *op, long chain, double theta, double *recrates, double *lamda, double *weight0array, int state);
double getrecrate(option_struct *op, long chain, double theta, double *recrates, double *lamda, double *weight0array, int state);
void getlogweight0array(option_struct *op, long chain, double theta_0, double *recrates_0, double *lamda_0, double *weight0array);
void getparameters(option_struct *op, long chain, double *weight0array);
double getlogrelativelikelihood(option_struct *op, long chain, double theta, double *recrates, double *lamda, double *weight0array);
void parameter_estimation_EM(option_struct *op, long chain, double theta, double *recrates, double *lamda);

void getalphas(treerec *tr, double *recrates, double *lamda);
void getbetas(treerec *tr, double *recrates, double *lamda);
void expectedrecrate(double *recrates, double *lamda, treerec *tr, int state, double *values, double weight);
void expectedtransitionrate(double *recrates, double *lamda, treerec *tr, int state, double *values, double weight);
void findoptimal_temp(option_struct *op, long chain, double theta,double *recrates, double *lamda, double *weight0array);

void hotstateprob(option_struct *op, long chain, double theta, double *recrates, double *lamda, double *weight0array);
