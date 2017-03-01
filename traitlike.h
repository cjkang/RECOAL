#include "recombine.h"

#define TRAITLIKE

void traitsiteplot(option_struct *op, data_fmt *data, double *llike);
void copytraits(node *p, node *q);
double traiteval (option_struct *op, data_fmt *data, tree *tr,
  double mu, double nu, double pd, long location);
void traitsmooth (option_struct *op, data_fmt *data, tree *tr,
  double mu, double nu, long location);
double probHH(double tyme, double mu, double nu);
double probHD(double tyme, double mu, double nu);
double probDH(double tyme, double mu, double nu);
double probHH(double tyme, double mu, double nu);
double probDD(double tyme, double mu, double nu);
void traitview (option_struct *op, data_fmt *data, node *p, double mu,
  double nu, long location);
void traitlike(option_struct *op, data_fmt *data, tree *tr, long numsites, 
  double mutrait, double traitratio, double pd, double *traitarray);
void traitprint(long numsites, double *traitarray, FILE *outfile, long numout);
void traitread(tree *tr, long numseq);
void traitresult(option_struct *op, data_fmt *data, double *traitarray,
  FILE *outfile, long numout);
