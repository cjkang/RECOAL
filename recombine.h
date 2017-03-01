


#define MEMDEBUG

#ifndef RECOMBINE_INCLUDE
#define RECOMBINE_INCLUDE
#endif

#define VERSION "0.40"

#ifdef  GNUDOS
#define DJGPP
#define DOS
#endif /* GNUDOS */

#ifdef THINK_C
#define MAC
#endif /* THINK_C */

#ifdef __MWERKS__
#undef HAVE_LGAMMA
#ifdef macintosh
#define MAC
#include <ansi_prefix.mac.h> /*fixes time problems*/
/* #include "unix.mach.h" */
#else /* macintosh */
#define DOS
#endif /* windows */
#endif /* MWERKS */


#ifdef GENERATINGPOWERPC
#define MAC
#endif

#ifdef GENERATING68K
#define MAC
#endif

#ifdef GENERATING68881
#define MAC
#endif

#ifdef SystemSevenFiveOrLater
#define MAC
#endif

#ifdef SystemSevenOrLater
#define MAC
#endif

#ifdef SystemSixOrLater
#define MAC
#endif

#ifdef __CMS_OPEN
#define CMS
#define EBCDIC true
#define INFILE "infile data"
#define OUTFILE "outfile data"
#define TREEFILE "treefile data"
#define FONTFILE "fontfile data"
#define PLOTFILE "plotfile data"
#define INTREE "intree data"
#define OUTTREE "outtree data"
#else
#define EBCDIC false
#define INFILE "infile"
#define OUTFILE "outfile"
#define TREEFILE "treefile"
#define FONTFILE "fontfile" /* on unix this might be /usr/local/lib/fontfile */
#define PLOTFILE "plotfile"
#define INTREE "intree"
#define OUTTREE "outtree"
#endif /* CMS_OPEN */

#ifdef L_ctermid            /* try and detect for sysV or V7. */
#define SYSTEM_FIVE
#endif /* L_ctermid */

#ifdef sequent
#define SYSTEM_FIVE
#endif /* sequent */

#ifndef MAC
#ifndef SYSTEM_FIVE
# include<stdlib.h>
# if defined(_STDLIB_H_) || defined(_H_STDLIB) || defined(H_SCCSID) || defined(unix)
# define UNIX
# define MACHINE_TYPE "BSD Unix C"
# endif /* defined.....*/
#endif /* SYSTEM_FIVE */
#endif /* MAC */

#ifdef __STDIO_LOADED
#define VMS
#define MACHINE_TYPE "VAX/VMS C"
#define printf vax_printf_is_broken
#define fprintf vax_fprintf_is_broken
void vax_printf_is_broken(const char *fmt,...);
void vax_fprintf_is_broken(FILE *fp,const char *fmt,...);
void vax_tweak_fmt(char *);
#endif /* __STDIO_LOADED */

#ifdef __WATCOMC__
#define QUICKC
#define WATCOM
#define DOS
#endif /* __WATCOMC__ */
/* watcom-c has graphics library calls that are almost identical to    *
 * quick-c, so the "QUICKC" symbol name stays.                         */

#ifdef _QC
#define MACHINE_TYPE "MS-DOS / Quick C"
#define QUICKC
#include "graph.h"
#define DOS
#endif /* _QC */

#ifdef _DOS_MODE
#define MACHINE_TYPE "MS-DOS /Microsoft C "
#define DOS           /* DOS is  always defined if  on a dos machine */
#define MSC           /* MSC is defined for microsoft C              */
#endif /* _DOS_MODE */

#ifdef __MSDOS__      /* TURBO c compiler, ONLY (no other DOS C compilers) */
#define DOS
#define TURBOC
#include<stdlib.h>
#include<graphics.h>
#endif /* __MSDOS__ */

#ifdef DJGPP          /* DJ's gnu  C/C++ port */
#include<graphics.h>
#endif

#ifndef MACHINE_TYPE
#define MACHINE_TYPE "ANSI C"
#endif

#ifdef DOS
#define MALLOCRETURN void 
#else
#define MALLOCRETURN void
#endif /* DOS */

#ifdef VMS
#define signed /* signed doesn't exist in VMS */
#endif /* VMS */

/* default screen types */
#ifdef DOS
#define IBMCRT true
#define ANSICRT false
#else
#ifdef MAC
#define IBMCRT false 
#define ANSICRT false
#else
#define IBMCRT false 
#define ANSICRT true 
#endif /* MAC */
#endif /* DOS */

#ifdef DJGPP
#undef MALLOCRETURN
#define MALLOCRETURN void
#endif /* DJGPP */


/* includes: */
#ifdef UNIX
#include<strings.h>
#include<string.h>
#else
#include<string.h>
#endif /* UNIX */

#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>

#ifdef MAC
#include "mac_interface.h"
#endif /* MAC */

#include "constants.h"

#ifdef __MWERKS__
#undef DBL_EPSILON
#undef DBL_MAX
#include <float.h>
#endif /* __MWERKS */

#define FClose(file) if (file) fclose(file) ; file=NULL

typedef unsigned char boolean;

/* warning WARNING debug DEBUG */
#define true    1
#define false   0

#ifndef TRUE
#define TRUE 1
#endif

#ifndef FALSE
#define FALSE 0
#endif

#define SETBITS 32

#ifndef DBL_EPSILON
#include <float.h>
#ifndef DBL_MAX
#define DBL_MAX ((double)1.7976931348623157e308)
#define DBL_EPSILON 2.2204460492503131e-16
#endif
#endif

typedef struct recnumb{
    long beg;
    long end;
    struct recnumb *next;
    int numrec;
    double swt;
    double last_weight;    
} recnumb;

typedef struct treerec {
  long *kk;
  double *kend;
  long *eventtype;
  double *actives;
  long numcoals;
  double numrecombs; /* double to allow fatal attraction avoidance */
  long *sitescore;
  boolean hastings_adjust;
  double tk, ts;
  struct recnumb *recs;
  double llike;

  double *rec;
  double *sumsite;

} treerec;



#ifdef MAC
MALLOCRETURN    *mymalloc(long);
#else
MALLOCRETURN    *mymalloc();
#endif

typedef unsigned int    base;
typedef double          sitelike[4];
typedef sitelike        *ratelike;
typedef ratelike        *phenotype;
typedef double          *contribarr;
typedef char            naym[NMLNGTH];
typedef long            longer[3];
typedef char            **sequence;


/* holds the cumulative data likelihoods in the tree structure */
typedef struct datalike_fmt {
    double **a; /* by loci by MICRO_ALLELEMAX */
    double ****s; /* by site by category by slice by base */
} datalike_fmt;

/* used in the tree structure */
typedef struct _node {
    struct _node *next, *back;
    char type;
    long id;
    long number;
  datalike_fmt *x;
    char *nayme;
    boolean top;
    double v, tyme, length;
    boolean updated;
    long futileflag;

    int update;

/* sequence specific stuff */
    long *oldbackranges;
    long *ranges;
    boolean *members;
    long memberstrue;

/* recombination specific stuff */
    long recstart, recend;
    long *coal;

/* migration specific stuff */
    long pop, actualpop;
    double lxmax;

/* disease trait likelihood specific stuff */
    double *z;        
  
  struct brlist *branch;

} node;
typedef struct seglist{
    struct seglist *prev, *next;
    long start,end;
    long numsam;
} seglist;

/* this is a list describing the tree as a series of horizontal
   slices of the tree.  Each slice starts at a node, and extends
   down until the next cut is encountered. */
typedef struct tlist {
  struct tlist *prev, *succ;
  node *eventnode;   /* points to a random tip for the tips */
  double age;        /* tyme from tree top                  */
  long numbranch;
  node **branchlist;
  struct seglist *segments, *oldsegments;
  int update;
  
} tlist;

typedef struct creature {
    long numhaplotypes;
    node **haplotypes;
    long numflipsites;
    long *flipsites;
} creature;

typedef struct tree {
  
    node   **nodep, **oldnodep;
    creature *creatures;
    double likelihood;
    double coalprob;
    node *root;
    tlist *tymelist;
    long numcoals;
    long numrecombs;


  double *dlikelihood;
  struct recnumb *recsummary;
  struct path *recpath;
  double *recrates, *weight_array;
  double *numrec_array; 

  } tree;

typedef struct reclist {
    struct reclist *prev, *succ;
    node *recnode;
    long startsite, endsite;
    boolean iscross;
} reclist;

typedef struct rlrec {
    double *val;
} rlrec;

typedef struct valrec {
    double rat_xi, rat_xv, zz, z1, y1, ww1, zz1, ww2, zz2, z1zz, z1yy,
	xiz1, xiy1_xv, ww1zz1, vv1zz1, ww2zz2, vv2zz2;
} valrec;

typedef struct dnadata {
    char *title;          /* datafile title, if any */
    char **popnames;      /* population names, if any */
    char ****indnames;    /* individual sequence names */
    long numpop;          /* number of populations */
    long numloci;         /* number of loci */
    long *numseq;         /* number of individuals per population */
    long *sites;          /* number of base pairs per locus */
    long **sitecount;     /* number of sites per locus per site */
    long **markersite;    /* map of markers into sites: per locus per marker */

    int *segranges0;      /* 4 special auxiliary arrays to speed up some */
    int *segranges1;      /* brlist, hidden-passages, code; */
    int *segranges2;      /* all referenced by site */
    int *segranges3;

    char ****seqs;        /* DNA sequences */
    double ***sspace;     /* spacing between DNA sites per pop per locus
			     per site */
    char sdlm;            /* spacing file delimiter for site/gap pairs */
    double freqa, freqc, freqg, freqt, freqr, freqy, freqar, freqcy,
	freqgr, freqty;    /* base frequencies */
    double xi;            /* transition pool prob */
    double xv;            /* transversion pool prob */
    double fracchange;    /* probability of change */
    double ttratio;       /* transition/transversion ratio */
    long *dnaweight;      /* arbitrary dna sitewise weighting factor */
} dnadata;

typedef struct msatdata {
    char *title;          /* datafile title, if any */
    char **popnames;      /* population names, if any */
    char ****indnames;    /* individual tip names */
    char dlm;             /* input file delimiter for diploids */
    long numpop;          /* number of populations */
    long numloci;         /* number of loci        */
    long *numind;         /* number of individuals per population */
    long ****msats;       /* micro satellite counts */
    double ***mspace;     /* spacing after a particular msat locus */
    double **steps;       /* table of transition probabilities from
			     state "x" to state "y". */
} msatdata;

typedef struct data_fmt {
    dnadata *dnaptr;
    msatdata *msptr;
    long *siteptr;   /* site aliasing array */
} data_fmt;

typedef struct option_struct {
  boolean ctgry, watt, printdata, usertree, progress,
    treeprint, interleaved, ibmpc, ansi, autocorr,
    freqsfrom, plump, weights, spacing, newdata, same_ne,
    interactive, mhmcsave, datatypeset, panel, map, fc,
    haplotyping, norecsnp, profile, print_recbythmaxcurve,hetero;
  char datatype;

  int hotspot;
  int holdings[5];
    long steps[NUM_TYPE_CHAINS], numchains[NUM_TYPE_CHAINS],
	increm[NUM_TYPE_CHAINS], numout[NUM_TYPE_CHAINS], holding, categs;
    long *numpanel;       /* number of individuals in SNP panel per
			     population */
    double *ne_ratio, *rate, *probcat, lambda;
    double mutrait, traitratio, pd; /* for mapper */
    long hapdrop; /* for haplotyping:  strat 0 = fliphap, 1 = single flipdrop,
		     2 = double flipdrop */
    double happrob;  /* probability of doing haplotype rearrangement */

/*  The next set of variables are for the "full" SNP option */
    boolean full;
    double chance_seen;

/* these variables handle different temperature chains */
    long *temperature, ctemp, numtempchains;
  int userrec;
/* these variables are for use in communicating the upper and lower
   bounds in theta and rec-rate from the confidence interval calculator
   to the likelihood curve printer */
    double thlb, thub, reclb, recub;

} option_struct;

typedef struct Rs{
    long beg;
    long end;
    float recrates;
    struct Rs *next;
} Rs;

typedef struct sumtree{
    int numcoal;
    float wlinks;
    struct recnumb *recinfo;
} sumtree;
 
typedef struct path{
    long beg,end;
    int state;
    struct path *next;
} path;
typedef struct trackset{
    struct track *hottrack;
    struct track *coldtrack;
    struct trackset *prev,*next;
} trackset;
typedef struct track{
    int state;
    long position;
    struct track *hotnext;
    struct track *coldnext;
    struct track *prev;
} track;

typedef struct state_sequence{
  int state;
  long position;
  struct state_sequence *prev;
  int numlink;
} state_sequence;


typedef struct treeinfo{
  long numrec_hot, numrec_cold;
  double weight_hot, weight_cold;
  long numcoal;
  double branch;
  long HH,HC,CC,CH;
} treeinfo;
  
typedef struct nodeinfo{
  struct nodeinfo *rightnode,*leftnode, *back; //rightone - smaller node id, leftone - bigger node id
  double rightlength, leftlength;
  long nodenum;
  long nodeid;
  //new one
  double toplength;
  int leftnum, rightnum;

  int istip; //tip = 1, nontip = 0
  double A,C,G,T;
  int updated;
  int isroot;
  //boolean SNP;
  //int nucleotide; // A = 0, C = 1, G = 2, T = 3
} nodeinfo;

typedef struct treeshape{
  struct nodeinfo *rootnode;
  double A,C,G,T,U;
  struct nodeinfo_new **nodep;
  double *probs;
} treeshape;

#define FILEP(A,B) (MENU ? (A) : (B))
#define ERRFILE FILEP(stderr,simlog)
#define REF_CHAIN(A) (((A) < op->numchains[0]) ? 0 : (A)+1-op->numchains[0])
#define TYPE_CHAIN(A) (((A) < op->numchains[0]) ? 0 : 1)
#define MAX(A,B) (((A) > (B)) ? (A) : (B))
#define MIN(A,B) (((A) < (B)) ? (A) : (B))
#define EXP(A) (((A) > EXPMAX) ? printf("\n\nTOO BIG\n") : exp(A))
#define NUM_CHROM(A) (((A) == 's') ? 1 : 2)
#define BOOLPRINT(A) ((A) ? "true" : "false")

#ifndef GETDATA_INCLUDE
#include "getdata.h"
#endif

void setupoption_struct(option_struct **op);
char lowercase(char c);
void openfile(FILE **fp, char *filename, char *mode, char *application,
	      char *perm);
boolean isrecomb(node *p);
boolean iscoal(node *p);
boolean branchsub(tlist *t, node *oldbranch, node *newbranch);
void insertaftertymelist(tlist *t, node *p);
void subtymelist(tlist *t, node *branchtop, boolean both);
void printtymelist(tlist *t);
node *findtop(node *p);
node *findunique(node *p);
node *otherdtr(node *p);
node *otherparent(node *p);
long findlink(node *p);
void free_z(option_struct *op, node *p);
void allocate_z(option_struct *op, data_fmt *data, node *p);
void free_x(option_struct *op, node *p);
void allocate_x(option_struct *op, data_fmt *data, node *p);
void VarMalloc(dnadata *dna, node *p, boolean allokate);
void init_coal_alloc(long **coal, long numcoalpairs);
void coal_Malloc(node *p, boolean allokate, long numcoalpairs);
void init_ranges_alloc(long **ranges, long numrangepairs);
void ranges_Malloc(node *p, boolean allokate, long numrangepairs);
void meld_adjacent_ranges(long *cranges, long newelem);
void addrange(long **newranges, long newstart, long newend);
boolean inrange(long *ranges, long site);
void subrangefc(long **newranges, long substart, long subend);
void printrange(long *ranges);
void newnode(node **p);
void freenode(node *p);
void newtymenode(tlist **t);
void freetymenode(tlist *t);
void freetymelist(tlist *t);
void hookup(node *p, node *q);
void atr(node *p);
void probatr(node *p);
tlist *gettymenode(tree *tr, long target);
double vtol(option_struct *op, data_fmt *data, double v);
void ltov(option_struct *op, data_fmt *data, node *p);
double findcoal_ltov(option_struct *op, data_fmt *data, double value);
void joinnode(option_struct *op, data_fmt *data, double length,
	      node *p, node *q);
void readparmfile(option_struct *op);
boolean whichopbool(option_struct *op, long i);
void rec_parmfilewrite(option_struct *op, long numloci, long numpop);
void readseedfile(void);
void print_menuheader(option_struct *op);
void print_startmenu(option_struct *op, boolean writeout);
void print_datamenu(option_struct *op);
void print_searchmenu(option_struct *op);
void print_menuend(void);
void initoptions(option_struct *op);
void workmenu1(option_struct *op, char ch, boolean *menu1,
	       boolean *writeout);
void workmenu2(option_struct *op, char ch, boolean *menu1,
	       long *numloci, long *numpop);
void checkparmfile(option_struct *op);
void getoptions(option_struct *op);
void firstinit(option_struct *op, data_fmt *data);
void popinit(option_struct *op, data_fmt *data);
void locusinit(option_struct *op, data_fmt *data);
void inputcategories(option_struct *op, data_fmt *data);
void end_of_population_free(option_struct *op, data_fmt *data);
void end_of_locus_free(option_struct *op, data_fmt *data);
void freetree(option_struct *op, data_fmt *data, tree *target);
void makesiteptr(option_struct *op, data_fmt *data);

node *allocate_nodelet(long num, char type);
void allocate_root(option_struct *op, data_fmt *data, tree *tr);
void allocate_tip(option_struct *op, data_fmt *data, tree *tr,
		  long num);
void allocate_interior(option_struct *op, data_fmt *data, tree *tr, long num);
void init_creature(creature *cr, long numhaplotypes);
void treesetup(option_struct *op, data_fmt *data);

void read_line(FILE *source, char *line);
void read_header(option_struct *op, long *numpop, long *numloci,
		 long **numsites, char *title, char *dlm);
void read_popheader(long *numind, char *poptitle);

void getinput(option_struct *op, data_fmt *data);
void orient(tree *tr, option_struct *op, data_fmt *data, node *p);
void plumptree(option_struct *op, data_fmt *data, double th0);
void finishsetup(option_struct *op, data_fmt *data, node *p);
void initbranchlist(option_struct *op, data_fmt *data);
void inittable(option_struct *op, dnadata *dna);
void initweightrat(option_struct *op, data_fmt *data);
void treeout(node *p, long s, FILE **usefile);
boolean nuview_usebranch(data_fmt *data, node *p, long site);
double prob_micro(msatdata *ms, double t, long diff);
void calcrange(option_struct *op, data_fmt *data, node *p, node *q,
	       node *r, long whichtree, long start, long finish, long numcategs);
void nuview(option_struct *op, data_fmt *data, node *p, long indexsite,
	    long startsite, long endsite);
void nuview_micro(option_struct *op, data_fmt *data, node *p, long indexsite,
		  long startsite, long endsite);
node *findcoal(option_struct *op, data_fmt *data, node *p, double *v);
void localsmooth(option_struct *op, data_fmt *data, node *p,
		 long indexsite, long startsite, long endsite);
void snpsmooth(option_struct *op, data_fmt *data, tree *tr,
	       long indexsite, long startsite, long endsite);
void localeval(option_struct *op, data_fmt *data, node *p, boolean first);
double micro_evaluate(option_struct *op, msatdata *ms, tree *tr,
		      long start, long end);
void seekch(char c);
void getch(char *c);
void processlength(node *p);
void addelement(option_struct *op, data_fmt *data, node *p, long *nextnode);
void treeread(option_struct *op, data_fmt *data);
void finddnasubtrees(option_struct *op, data_fmt *data, tlist *tstart,
		     long *sranges);
void findsubtrees_node(option_struct *op, data_fmt *data,
		       tlist *tlast, tree *tr, long **coal);
void drop_findsubtrees(option_struct *op, data_fmt *data, tree *tr, 
		       node *p);
void findsubtrees(option_struct *op, data_fmt *data, tlist *tstart,
		  long *sranges);
void findsubtreemarkers(option_struct *op, data_fmt *data, long tstart,
			long tend, long *mstart, long *mend);
long markertosite(option_struct *op, data_fmt *data, long marker);
long sitetorightmarker(option_struct *op, data_fmt *data, long site);
long sitetomarker(option_struct *op, data_fmt *data, long site);
boolean sameranges(long *range1, long *range2);
void copyranges(node *source, node *target);
void copycoal(node *source, node *target);
void copylikes(option_struct *op, data_fmt *data, node *source,
	       node *target);
void copynode(option_struct *op, data_fmt *data, node *source,
	      node *target);
void addnode(option_struct *op, data_fmt *data, node *p, node *q);
node *getrecnodelet(node *copy, node *orig);
void make_tree_copy(option_struct *op, data_fmt *data, tree *source,
		    tree *target);
void copycreature(option_struct *op, data_fmt *data, creature *source,
		  creature *target, tree *tr);
void copycreatures(option_struct *op, data_fmt *data, tree *source,
		   tree *target);
tree *copytree(option_struct *op, data_fmt *data, tree *source);
void constructtree(option_struct *op, data_fmt *data, double branchlength);
boolean rearrange(option_struct *op, data_fmt *data, long whichchain);
void temptreeswap(option_struct *op, data_fmt *data, boolean *changed,
		  long chain);
void maketree(option_struct *op, data_fmt *data);
void freenodelet(option_struct *op, data_fmt *data, node *p);
void finalfree(option_struct *op, data_fmt *data);
double randum(void);
double lengthof(node *p);
void fixlength(option_struct *op, data_fmt *data, node *p);
double watterson(option_struct *op, data_fmt *data);
void invardatacheck(option_struct *op, data_fmt *data);
void buildtymelist(tree *tr, option_struct *op, data_fmt *data, node *p);
double eval_calcrange(option_struct *op, data_fmt *data, tree *tr,
		      boolean first, long start, long finish, long numcategs);
double snp_eval_calcrange(option_struct *op, data_fmt *data, tree *tr,
			  boolean first, long start, long finish, long numcategs);
double panel_eval_calcrange(option_struct *op, data_fmt *data, tree *tr,
			    boolean first, long start, long finish, long numcategs);
double evaluate(option_struct *op, data_fmt *data, tree *tr,
		double llike);
double snp_evaluate(option_struct *op, data_fmt *data, tree *tr,
		    double llike, long start, long finish);
boolean hasrec_chain(option_struct *op, long chain);
void chainendcheck(option_struct *op, data_fmt *data, tree *tr,
		   long chain, boolean locusend);
double min_theta_calc(option_struct *op, data_fmt *data, double th);
long boolcheck(char ch);
long countbranches(option_struct *op, data_fmt *data, tree *tr);
int main(int argc, char *argv[]);
boolean booleancheck(option_struct *op, char *var, char *value);
boolean numbercheck(option_struct *op, char *var, char *value);
boolean sitecompare(dnadata *dna, long site1, long site2);
boolean testratio(option_struct *op, data_fmt *data, tree *oldtree,tree *newtree, char ratiotype);

boolean istip(node *p);
void makeinvarvalues(dnadata *dna, long numcategs, tree *tr);
void getnodata(option_struct *op, data_fmt *data);
void rec_outtree(node *p, boolean first, FILE **usefile);

int eof(FILE *f);
int eoln(FILE *f);
float findrate(long pos);
void setupinterval();
boolean rearrangeH(option_struct *op, data_fmt *data,long whichtchain, long chain);
boolean rearrangeHMM(option_struct *op, data_fmt *data, long whichtchain);
struct path *createnewpath(struct path *curpath);
struct path *copypath(struct path *old);
struct path *findblock(struct path *paths, long pos);
struct Rs *updateR(struct path *statepath);
struct Rs *convertR(struct path *statepath);
struct path *rearrangestruct(struct recnumb *currec);
float pGS(struct recnumb *currec, struct Rs *struc);
void sumstruct(struct path *pathstruct, long chain);
int findstate(struct path *curpath, long site);
struct path *viterbi(struct recnumb *currec);
struct path *randomviterbi(struct recnumb *currec);
struct path *randomviterbi_1(struct recnumb *currec);
void constructseglist(option_struct *op, data_fmt *data, tree *tr); 
void freesequence(state_sequence *target, int key);
void calcinitialtree(option_struct *op, data_fmt *data, tree *inittree, long *ranges);
void calcratio(option_struct *op, data_fmt *data, tree *oldtree, tree *newtree, long *ranges);
double computedlikelihood(option_struct *op, data_fmt *data, long i, treeshape *subtree, double weight);
treeshape *findsubs(node *tipnode, long start);
void climbright(node *localroot, node *leftnode, long start, nodeinfo *rootinfo);
void climbup(node *localroot, long start, nodeinfo *curinfo);
void calclikelihood(option_struct *op, data_fmt *data, long position, nodeinfo *localroot);
void checkbase(option_struct *op, data_fmt *data, long position, nodeinfo *tipnode);
void calcrange_modified(option_struct *op, data_fmt *data, nodeinfo *localroot);
void nuview_modified(option_struct *op, data_fmt *data, long position, nodeinfo *localroot);
double dlikelihood(option_struct *op, data_fmt *data, long position, nodeinfo *rootnode);
void empiricaldnafreqs_modified(option_struct *op,data_fmt *data);
char checkSNP(option_struct *op, data_fmt *data, long position, double freq);

void constructseqinfo(option_struct *op, data_fmt *data, double freq);

void freesubtree(treeshape *target);

void freenodeinfo(nodeinfo *rootnode);
void backward(double *cold, double *hot, recnumb *target);
struct path *randompathsampling(struct recnumb *currec);
void convertrecarray(struct path *source);
void randompathsampling_temp(double *oldrecarray, double *newrecarray,double *weight, int *recnums, long start, long end);
void backward_temp(double *cold, double *hot, double *weight, int *recnum, long start, long end, char lastone);
void scorerecs_temp(tree *target);
void addtime(node *target, double *weights);

void realrandompathsampling_temp(double *oldrecarray, double *newrecarray, double *weight, int *recnums, long start, long end);
void realbackward_temp(double *cold, double *hot, double *weight, int *recnum, long start, long end, char lastone);
void whereisrec(int *recnum);
void temptreeswapH(option_struct *op, data_fmt *data, boolean *changed,long chain);
void generateuniformS();
boolean structureheating(tree *target, double *newS, double *oldS, double temperature);
boolean simulatedtemperingswap(double prob, double C1, double C2, double temperature1, double temperature2);
double treelikelihood(tree *target);
double probS(double *recstructure);
double probGS(tree *target, double *rates);
void transformtorecarray(data_fmt *data);
void transformtorecspace(data_fmt *data);
void transformtree(data_fmt *data,treerec *target);
void initrecarray(double r0);

void randompathgenerator();
double forward(tree *curtree);
void generaterecstructure(tree *target);
double logratiowithunirec(tree *oldtree, tree *newtree);
double hmmlikelihood(tree *target);
double treepartlikelihood(tree *newtree, tree *oldtree);
void generaterecstructure_temp(tree *target, double *rec, double *tran);
void getbetaratio_temp(tree *target, double *betaarray0, double *betaarray1,double *rec, double *tran);
double hmmlikelihood_temp(tree *target, double *rec, double *tran);
double hmmlikelihood_two(tree *target, int *otherrecarray, double *otherweightarray, double temperature);
double hmmlikelihood_with_temp(tree *target, double temperature);

void generaterecstructure_with_temp(tree *target, double temperature);

double treeS(long pos, int *recs, double *weights);

boolean hastingswithS(tree *oldtree, tree *newtree);

void stronghotspotputter(tree *target);

boolean randomhottestratio(tree *oldtree, tree *newtree);
double *getalpha1(tree *tr);
double *getbeta1(tree *tr);
double flatrecombinationrateratio(tree *oldtree, tree *newtree);
double hmmlikelihood_two_with_bigger_rec(tree *target, int *otherrecarray, double *otherweightarray);
void findhotspots(option_struct *op, data_fmt *data);
void temptreeswap_rec(option_struct *op, data_fmt *data, boolean *changed,long chain);
void randomhotspotputter();
void generaterecstructure_with_big_lamda(tree *target, double temperature);
void getbetaratio_with_temp(tree *target, double *betaarray0, double *betaarray1, double temperature);
void getbetaratio_with_temp(tree *target, double *betaarray0, double *betaarray1, double temperature);
void getinirate(double *recs, double *weights);
void getalphabeta(double *recs, double *weights, double *recombrates, double *lamdas, double *alpha0, double *alpha1, double *beta0, double *beta1);


void addnewhap(option_struct *op, data_fmt *data);
void maketree_S(option_struct *op, data_fmt *data);
void simhaplotype(option_struct *op, data_fmt *data);
void computenonnewtipnodes(option_struct *op, data_fmt *data, long site, nodeinfo *localroot);
void computenewtipprob(option_struct *op, data_fmt *data, long position, double length, nodeinfo *tipnode, nodeinfo *localroot);
void nuview_single(option_struct *op, data_fmt *data, long position, double length, nodeinfo *branchnode, nodeinfo *localroot);
void calcrange_single(option_struct *op, data_fmt *data, nodeinfo *branchnode, nodeinfo *localroot);
void computenewtipnodes(option_struct *op, data_fmt *data, long site, nodeinfo *localroot);
void simhaplotype_old_2(option_struct *op, data_fmt *data);
void calclikelihood_type(option_struct *op, data_fmt *data, long position, nodeinfo *localroot, int seqtype);
void nuview_modified_type(option_struct *op, data_fmt *data, long position, nodeinfo *localroot, int seqtype);
void checkbase_type(nodeinfo *tipnode, int seqtype);
void simhaplotype_old_2(option_struct *op, data_fmt *data);
void addprobs(treeshape *result);