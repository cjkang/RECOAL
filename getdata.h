#define GETDATA_INCLUDE

#ifndef RECOMBINE_INCLUDE
#include "recombine.h"
#endif

#define OUTLINESIZE 50 /* the maximum number of characters allowed
                          for line-output used in printdnadata */

char uppercase(char ch);
void setupdata(data_fmt *data, char datatype, long numpop, long numloci,
 long numind, long *numsites, char *title);
void freedata(data_fmt *data);
long getdata_nummarkers(option_struct *op, data_fmt *data);
long getdata_numtips(option_struct *op, data_fmt *data);
long getdata_numseq(option_struct *op, data_fmt *data);
long getdata_numloci(option_struct *op, data_fmt *data);
long getdata_numpop(option_struct *op, data_fmt *data);
long getdata_sitecount(option_struct *op, data_fmt *data, long marker);
long getdata_markersite(option_struct *op, data_fmt *data, long marker);
double getdata_space(option_struct *op, data_fmt *data, long target);

void read_spacefile(FILE *file, option_struct *op, data_fmt *data);
int longcmp (const void *v1, const void *v2);
int doublecmp (const void *v1, const void *v2);
void read_flipfile(option_struct *op, data_fmt *data, tree *tr, FILE *file);
void pruneflips(option_struct *op, data_fmt *data, tree *tr);


void read_data(FILE *infile, data_fmt *data, char datatype, long pop);
void setpopstuff(data_fmt *data, long pop, long numind,
   char *poptitle, char datatype);
void read_popdata(FILE *infile, data_fmt *data, long pop,
   option_struct *op);
void read_indname(FILE *infile, data_fmt *data, long pop, long lowcus,
   long ind, long namelngth, char datatype);

void read_microalleles(FILE *infile, msatdata *data, long pop,
   long ind);

long read_ind_seq(FILE *infile, dnadata *data, option_struct *op,
   long lowcus, long pop, long ind, long baseread);
void finish_read_seq(FILE *infile, data_fmt *data, option_struct *op,
   long pop, long baseread);

void setupdnadata(dnadata **dna, long *sites, long numseq, long numloci,
  long numpop);
void freednadata(dnadata *dna);
void printdnadata(option_struct *op, data_fmt *data, FILE *out);
void initinvartips(option_struct *op, data_fmt *data, long categs,
  tree *curtree);
void makednavalues(option_struct *op, data_fmt *data, long categs,
  tree *curtree);
void empiricaldnafreqs(option_struct *op, data_fmt *data, tree *curtree);
void getbasednafreqs(dnadata *dna, option_struct *op, double
  locus_ttratio, FILE *outfile);
void inputdnaweights(long numchars, dnadata *dna, option_struct *op);

long hapdist(option_struct *op,data_fmt *hap1, data_fmt *hap2, creature *critter);
void readtruehaps(option_struct *op, data_fmt *hap, long numpop, long
  numloci, long numind, long *numsites);
void copyhaps(option_struct *op, data_fmt *oldhap, data_fmt *newhap,
  long numpop, long numloci, long numind, long *numsites);
