#include "recombine.h" 
#include "getdata.h"
#ifdef DMEMDEBUG
#define MEMDEBUG
#include "memdebug.h"
#endif

#ifdef DMALLOC_FUNC_CHECK
#include "/usr/local/include/dmalloc.h"
#endif


/* Haplotype simulator based on Recombine-H*/

/* Recombine-H version 0.1 modified from Recombine 0.4 with permission from M. Kuhner*/

/* Written by Chul Joo Kang */

/* This program was modified to estimate heterogineity of recombination rate ("recombination hotspot") */

/* new FC routine was implemented */ 

/* Modified some rutines to improve running speed. */

//simulator option
boolean sim_mode = TRUE;
long numMCMC,minMCMC;
int numrep,repinter;
int old_hap, new_hap, new_tip = 0;
double growth;
short int **newsimdata;
short int **olddata;
short int *newhaps;
double esnp;
boolean menuoption;
data_fmt truehaps;
data_fmt starthaps;
int filtered = 0;

/* these definitions will be extern in other files*/
extern double **reci, **theti;
extern treerec ***sum;
extern unsigned int baseA, baseC, baseG, baseT;

extern long countcoalbranches(tree *tr);
extern void model_alloc(option_struct *op, data_fmt *data);
extern void scoretree(option_struct *op, data_fmt *data, long chain);
extern void rec_estimatev(option_struct *op, data_fmt *data, long lowcus,
						  long chain, boolean locusend);
extern void parameter_estimation_EM(option_struct *op, long chain, double theta, double *recrates, double *lamda);
extern boolean twiddle(option_struct *op, data_fmt *data);
extern boolean fliphap(option_struct *op, data_fmt *data, tree *tr);
extern boolean flipdrop(option_struct *op, data_fmt *data, tree *tr,
						long numdrop);
extern void add_one_recomb(option_struct *op, data_fmt *data, tree *tr,
						   long *siteptr, long chain);
extern void addfractrecomb(option_struct *op, data_fmt *data, 
						   long chain, treerec ***treesum);
extern boolean makedrop(option_struct *op, data_fmt *data);
extern void setupdnadata(dnadata **dna, long *numsites, long numseq,
						 long numloci, long numpop);
extern void freednadata(dnadata *dna);
extern void empiricaldnafreqs(option_struct *op, data_fmt *data,
							  tree *curtree);
extern void inputdnaweights(long numchars, dnadata *dna, option_struct *op);
extern void remove_branch(tlist *start, node *target);
extern void rename_branch(tlist *start, node *source, node *target);
extern void rec_scoreread(option_struct *op, data_fmt *data);
extern void unflag(node *p);
extern void traverse_unflag(tree *tr, node *p);
extern void read_spacefile(FILE *file, option_struct *op, data_fmt *data);
extern void printdnadata(option_struct *op, data_fmt *data, FILE *out);

extern long countsites(option_struct *op, data_fmt *data);
extern long getdata_sitecount(option_struct *op, data_fmt *data,
							  long marker);
extern double getnumlinks(option_struct *op, data_fmt *data, long start,
						  long finish);

/* functions for haplotyping */
extern void read_flipfile(option_struct *op, data_fmt *data, tree *tr,
						  FILE *file);
extern void pruneflips(option_struct *op, data_fmt *data, tree *tr);
extern void copyhaps(option_struct *op, data_fmt *oldhap, data_fmt *newhap,
					 long numpop, long numloci, long numind, long *numsites);

/* functions for heating */
extern double coalprob(option_struct *op, data_fmt *data, tree *tr, double
					   theta, double r);
extern int longcmp(const void *v1, const void *v2);

/* functions from traitlike.c for trait likelihoods */
extern void copytraits(node *source, node *target);
extern void traitlike(option_struct *op, data_fmt *data, tree *tr, 
					  long numsites, double mutrait, double traitratio, double pd,
					  double *traitarray);
extern void traitread(tree *tr, long numseq);
extern void traitprint(long numsites, double *traitarray, FILE *outfile, long numout);
extern void traitsiteplot(option_struct *op, data_fmt *data, double *llike);
extern void traitresult(option_struct *op, data_fmt *data,
						double *traitarray, FILE *outfile, long numout);

extern boolean twiddleH(option_struct *op, data_fmt *data);
extern boolean makedropH(option_struct *op, data_fmt *data);
extern void scoretreeH(option_struct *op, data_fmt *data, long chain);
extern seglist *copyseglist(seglist *original);
extern void freesegment(seglist *target);
extern boolean makedrop_S_new(option_struct *op, data_fmt *data);


extern double Hastingsratio(tree *oldtree, tree *newtree, long site);
//test

double sumrech = 0, sumsiteh = 0;
double sumrecc = 0, sumsitec = 0;

extern recnumb *scorerecs(option_struct *op, data_fmt *data, tree *target);
extern void freeR(Rs *target);

long part1 = 0, part2 = 0,part3 = 0, part4= 0,part5 = 0,part6 = 0, part7 = 0,part8 = 0, part9= 0,part10 = 0;
extern void checkrecnum(recnumb *target);
long tcount = 0, hotcount = 0;
double firstevent = 0, secondevent = 0.0, thirdevent = 0.0, one = 0.0, two = 0.0, three = 0.0;

tree **temptrees;/* the current tree for each temperature-chain */
tree *curtree;   /* the current tree */
long locus;      /* which locus are we currently working on? */
long population; /* which population are we currently working on? */
longer seed;     /* array used to hold randum seed */
long numdropped = 0;               /* number of trees dropped from
									consideration */
long notip = -99;

//test
double *newarray, constant[3];
double Acc = 0 ,Rej = 0;
double likechange1 = 0,likechange2 = 0,likechange3 = 0,likechange4 = 0;
double inhot = 0, outhot = 0;


long sumcoal = 0, sumrec1 = 0, sumrec2 = 0;
int israndom;
double *gpartweight;
long *gpartrec;
long marker1pos = -99,marker2pos = -99, marker3pos, marker4pos;
double recrates0,recrates1;
long sumhot[1000];
long targetsite;
long hott = 0, culd = 0;
extern long HtH,CtC,HtC,CtH;
struct treeinfo **treesummary;
extern void checkpath(path *target);
node *markednode = NULL;
double *calling = NULL;
//for non-HMM
double *rec_space;
long numspaces;
long num_markers;
long *spaceinfo;

double templike = 0;

double totalbranch_s = 0;
int expsnpnum  = 0;

/* next are a largish set of temporary variables used to hold various
 values before their structures are setup.  Mostly data variables. */
double freqa, freqc, freqg, freqt; /* dna freqs */
double locus_ttratio;              /* dna transition/transversion ratio */

FILE *infile, *outfile, *treefile, *bestree, *seedfile, *simlog,
*parmfile, *intree, *weightfile, *spacefile;
long rootnum, apps, col, numtrees, totchains;
boolean  **sametree;
double clearseed, theta0, rec0, branch0;
long		*category;
contribarr	*contribution;
node		*freenodes;
rlrec		*alogf;
char            gch;
long            nodectr;
double	 	sumweightrat, watttheta;
double	 	*weightrat;
valrec   	*tbl;
long		slid, slacc, hap, hacc, swap, swacc, indecks;
int S_O, S_H, S_N;

long windowsize, windowpos,slidesize,windowbeg,windowend;

Rs *R0;
recnumb *recnum0;
long *tempstructs[4];
struct sumtree **tempevaltrees;
struct path *newpath;
long nsites;
double lamda[2], curlamda[2];
double recrates[2], currates[2];
extern struct sumtree *treescore(option_struct *op, data_fmt *data,tree *target);
long startposition,endposition;
double MAG = 1.5;

double cutofffreq = 0;
char *dna_sequence;
double *recarray;
long seq_length, seq_num;

int *newgrec, *newGrec, *oldgrec, *oldGrec, *unchangedGrec;
double *newgweight, *newGweight, *oldgweight, *oldGweight, *unchagedGweight;

long numburnin_current;

double *tiplikelist;

double coldrecs = 0, hotrecs = 0, coldrectyme = 0, hotrectyme= 0;
double freqw1 =0, freqw2 = 0;
/* the following are for reading in parameters (readparmfile),
 and also writing them back out again (rec_parmfilewrite) */
char *booltokens[NUMBOOL] = {"interleaved","printdata","progress",
	"print-trees","freqs-from-data","categories",
	"watterson", "usertree", "autocorrelation",
	"newdata","same-ne","interactive","mhmcsave",
	"panel","map","final-coalescence","full-snp",
	"haplotyping","profile-like","norecsnp","userrates","menu"},
*numbertokens[NUMNUMBER] = {"ttratio","short-chains",
	"short-steps","short-inc","long-chains",
	"long-inc","long-steps","rec-rate",
	"holding","mu-ratio","trait-ratio",
	"dis-freq","hapdrop","heating","recC","recH","lamdaC","lamdaH"};

/* the following are for managing x array recycling */
long numx;
datalike_fmt *sparex[XARRAYSIZE];

void openfile(FILE **fp, char *filename, char *mode, char *application,
			  char *perm)
{
	FILE *of;
	char file[LINESIZE];
	
	strcpy(file,filename);
	while (1){
		of = fopen(file,mode);
		if (of)
			break;
		else {
			switch (*mode){
				case 'r':
					fprintf(stdout,"%s:  can't read %s\n",application,file);
					file[0] = '\0';
					while (file[0] =='\0'){
						fprintf(stdout,"Please enter a new filename>");
						//gets(file);
						fgets(file, LINESIZE, stdin);
					}
					break;
				case 'w':
					fprintf(stdout,"%s: can't write %s\n",application,file);
					file[0] = '\0';
					while (file[0] =='\0'){
						fprintf(stdout,"Please enter a new filename>");
						//gets(file);
						fgets(file, LINESIZE, stdin);
					}
					break;
			}
		}
	}
	*fp=of;
	if (perm != NULL)
		strcpy(perm,file);
} /* openfile */


char lowercase(char c)
{
	return((char)tolower((int) c));
} /* lowercase */


double treeprobs(tree *target){
	tlist *t;
	int i;
	double tk=0;
	double tyme, prob;
	
	for (t = target->tymelist; t->succ != NULL; t=t->succ) {
		tyme = t->age - t->eventnode->tyme;
		tk += t->numbranch * (t->numbranch - 1.0) * tyme;
	}
	
	prob = -(tk/theta0) + (log(2.0/theta0))*target->numcoals;
	
	for (i = 0; i<seq_length; i++){
		prob += log(recarray[i]) * (target->numrec_array[i]) 
		- recarray[i] * (target->weight_array[i]);
	}
	
	return prob;
}




double randum(void)
/* Mary's version--faster but needs 32 bits.  Loops have been unrolled
 for speed. */
{
	long newseed0, newseed1, newseed2;
	
	newseed0 = 1549*seed[0];
	newseed1 = newseed0/2048;
	newseed0 &= 2047;
	newseed2 = newseed1/2048;
	newseed1 &= 2047;
	newseed1 += 1549*seed[1] 
    + 812*seed[0];
	newseed2 += newseed1/2048;
	newseed1 &= 2047;
	newseed2 += 1549*seed[2]
    + 812*seed[1];
	
	seed[0] = newseed0;
	seed[1] = newseed1;
	seed[2] = newseed2 & 1023;
	return (((seed[0]/2048.0 + seed[1])/2048.0 + seed[2])/1024.0);
}  /* randum */


/****************************************************
 * setupoption_struct() creates a new option_struct */
void setupoption_struct(option_struct **op)
{
	
	(*op) = (option_struct *)calloc(1,sizeof(option_struct));	
	(*op)->rate = (double *)calloc(1,sizeof(double));
	(*op)->probcat = (double *)calloc(1,sizeof(double));
	(*op)->temperature = NULL;
	(*op)->numpanel = NULL;
	
} /* setupoption_struct */


/********************************************************************
 * istip() returns TRUE if the nodelet is a tip and FALSE otherwise */
boolean istip (node *p)
{
	return(p->type == 't');
} /* istip */


/***********************************************************************
 * iscoal returns TRUE if the nodelet is a coalescent, FALSE otherwise */
boolean iscoal (node *p)
{
	return(p->type == 'c');
} /* iscoal */


/************************************************************************
 * isrecomb returns TRUE if the nodelet is recombinant, FALSE otherwise */
boolean isrecomb (node *p)
{
	return(p->type == 'r');
} /* isrecomb */


/******************************************************************
 * branchsub() substitutes one branch for another in one interval *
 * in a tymelist.                                                 */
boolean branchsub(tlist *t, node *oldbranch, node *newbranch)
{
	long i;
	boolean found;
	
	found = FALSE;
	for(i = 0; i < t->numbranch; i++)
		if (t->branchlist[i] == oldbranch) {
			t->branchlist[i] = newbranch;
			found = TRUE;
		}
	
	return(found);
} /* branchsub */


/************************************************************************
 * insertaftertymelist adds the node "p" to a tymelist after, inserting *
 *    within the tymeslice pointed to by "t" using bisection.           */
void insertaftertymelist(tlist *t, node *p)
{
	long i, j, temp;
	tlist *s, *prevtime;
	/*   boolean b1, b2; */
	
	newtymenode(&s);
	
	//copy seglist
	//s->segments = copyseglist(t->segments);
	
	s->update = 0;
	
	s->eventnode = findtop(p);
	s->age = t->age;
	t->age = p->tyme;
	s->prev = t;
	s->succ = t->succ;
	t->succ = s;
	
	if (s->succ != NULL) s->succ->prev = s;
	
	/* now deal with the branchlist for the new tymelist entry;
     the "if" setting "temp" is a workaround for not being able to
     allocate zero size chunks of memory by malloc-debug/gcc */
	if (t->numbranch-1 > 0) temp = t->numbranch - 1;
	else temp = t->numbranch;
	s->branchlist = (node **)calloc(1,temp*sizeof(node *));
	s->numbranch = 0;
	
	prevtime = t;
	while (prevtime->update == 1){
		prevtime = prevtime->prev;
	}
	
	s->numbranch = prevtime->numbranch;
	
	i = 0;
	if (iscoal(s->eventnode) && s->eventnode->coal != NULL) if (s->eventnode->coal[0] == 0)   i = 1;
	
	
	if (iscoal(s->eventnode)){
		s->numbranch--;
		if (s->numbranch == 0) s->numbranch++;
		free(s->branchlist);
		s->branchlist = (node **)malloc(s->numbranch * sizeof(node *));
		s->branchlist[0] = s->eventnode;
		
		if (i != 1){
			for (j = 1 ; i < prevtime->numbranch; i++){
				if (prevtime->branchlist[i] != s->eventnode->next->back && prevtime->branchlist[i] != s->eventnode->next->next->back){
					s->branchlist[j] = prevtime->branchlist[i];
					j ++;
				}
			}
		}
		else{ //FC case
			s->numbranch --;
			if (s->numbranch == 0){
				s->numbranch++; // final root
				s->branchlist[0] = s->eventnode;
			}
			else{
				for (j = 0, i = 0; i < prevtime->numbranch; i++){
					if (prevtime->branchlist[i] != s->eventnode->next->back && prevtime->branchlist[i] != s->eventnode->next->next->back){
						s->branchlist[j] = prevtime->branchlist[i];
						j ++;
					}
				}
				s->branchlist[s->numbranch] = s->eventnode;
			}
		}
		
		return;
		
	}
	else{
		s->numbranch++;
		free(s->branchlist);
		s->branchlist = (node **)malloc(s->numbranch * sizeof(node *));
		s->branchlist[0] = findunique(s->eventnode)->next;
		s->branchlist[1] = findunique(s->eventnode)->next->next;
		for (j = 2; i < prevtime->numbranch; i++){
			//printf("num %d, %d\n",i,j);
			if (prevtime->branchlist[i] != findunique(s->eventnode)->back){
				s->branchlist[j] = prevtime->branchlist[i];
				j ++;
			}
		}
		return;
	}
	
} /* insertaftertymelist */


/******************************************************************
 * subtymelist() completely removes the tymelist entry passed in. *
 * the passed in branch should be the exact branch that wants to  *
 * be removed in the case of a recombinant node, or the branch    *
 * that will "replace" the missing branch in a coalescent node.   */
void subtymelist(tlist *t, node *branch, boolean all)
{
	
	if (isrecomb(t->eventnode)) {
		remove_branch(t,branch);
		if (all) remove_branch(t,otherparent(branch));
		else rename_branch(t,findunique(t->eventnode)->back,otherparent(branch));
	} else 
		if (all) remove_branch(t,t->eventnode);
		else rename_branch(t,branch,t->eventnode);
	
	if (t->prev) {t->prev->age = t->age; t->prev->succ = t->succ;}
	else {fprintf(ERRFILE,"ERROR:tried to remove first tymelist entry!\n");}
	if (t->succ) t->succ->prev = t->prev;
	
	freetymenode(t);
	
} /* subtymelist */


/****************************************************************
 * printtymelist() prints the entire tymelist: a debug function */
void printtymelist(tlist *t)
{
	long i;
	
	fprintf(ERRFILE,"TYMELIST BEGINS\n");
	while (t != NULL) {
		fprintf(ERRFILE,"%3ld age%8.6f--",
				t->eventnode->number, t->age);
		for (i = 0; i < t->numbranch; i++) {
			fprintf(ERRFILE," %3ld", t->branchlist[i]->number);
			if (t->branchlist[i]->top) fprintf(ERRFILE,"t ");
		}
		fprintf(ERRFILE,"\n");
		t = t->succ;
	}
	fprintf(ERRFILE,"TYMELIST ENDS\n");
}  /* printtymelist */


/*****************************************************************
 * lengthof() returns the length of the branch "rootward" of the *
 * passed node                                                   */
double lengthof(node *p)
{
	return fabs(p->tyme - p->back->tyme);
}  /* lengthof */


/******************************************************************
 * fixlength() sets the length and v values for a mother-daughter *
 * pair of nodelets.                                              */
void fixlength(option_struct *op, data_fmt *data, node *p)
{
	double newlength;
	
	newlength = lengthof(p);
	p->length = newlength;
	ltov(op,data,p);
	p->back->length = newlength;
	ltov(op,data,p->back);
	
} /* fixlength */


/************************************************
 * findtop() finds a "top" nodelet in the node. */
node *findtop(node *p)
{
	
	if (!istip(p))
		while (!p->top) p = p->next;
	
	return(p);
	
} /* findtop */


/********************************************************************
 * findunique() returns the unique nodelet of a node (i.e. "bottom" *
 * of a recombinant node or "top" of a coalescent node).            */
node *findunique(node *p)
{
	
	if (istip(p)) return(p);
	
	if (isrecomb(p))
		while((p)->top) p = p->next;
	else
		while (!p->top) p = p->next;
	
	return(p);
	
}  /* findunique */


/*********************************************************************
 * otherdtr() finds the other daughter nodelet of a node, it assumes *
 * that a daughter nodelet was passed to it.                         */
node *otherdtr(node *p)
{
	node *q;
	
	if (isrecomb(p)) {
		fprintf(ERRFILE,"ERROR:otherdtr() passed a recombinant node!\n");
		return(NULL);
	}
	
	if (p->top) {
		fprintf(ERRFILE,"ERROR:otherdtr() passed a top nodelet!\n");
		return(NULL);
	}
	
	for(q = p->next; q->top; q = q->next)
		;
	
	return(q);
	
} /* otherdtr */


/*********************************************************************
 * findlink() takes a recombinant node and returns the # of the site *
 * to the left of the cut link; or FLAGLONG in error states.         */
long findlink(node *p)
{
	node *q;
	
	if (!isrecomb(p)) {
		fprintf(ERRFILE,"ERROR:findlink passed non-recombinant node\n");
		return(FLAGLONG);
	}
	
	q = findunique(p);
	if (q->next->recstart) return(q->next->next->recend);
	else return(q->next->recend);
	
} /* findlink */


/**********************************************************************
 * otherparent() finds the other parent nodelet of a node, it assumes *
 * that a parent nodelet was passed to it.                            */
node *otherparent(node *p)
{
	node *q;
	
	if (!isrecomb(p)) {
		fprintf(ERRFILE,"ERROR:otherparent() passed a non-recombinant node!\n");
		return(NULL);
	}
	
	if (!p->top) {
		fprintf(ERRFILE,"ERROR:otherparent() passed a non-top nodelet!\n");
		return(NULL);
	}
	
	for(q = p->next; !q->top; q = q->next)
		;
	
	return(q);
	
} /* otherdtr */


/*********************************************
 * free_z() frees the "z" field of a nodelet */
void free_z(option_struct *op, node *p)
{
	if (p->z == NULL) return;
	
	free(p->z);
	p->z = NULL;
	
} /* free_z */


/***************************************************************
 * allocate_z() allocates space for the "z" field of a nodelet */
void allocate_z(option_struct *op, data_fmt *data, node *p)
{
	if (p->z) free(p->z);
	
	p->z = (double *)calloc(NUMTRAIT,sizeof(double));
	
} /* allocate_z */


/*********************************************
 * free_x() frees the "x" field of a nodelet */
void free_x(option_struct *op, node *p)
{
	
	if (p->x == NULL) return;
	
	if (numx < XARRAYSIZE) {
		sparex[numx] = p->x;
		numx++;
	} else {
		switch(op->datatype) {
			case 'a':
				break;
			case 'b':
			case 'm':
				free(p->x->a);
				break;
			case 'n':
			case 's':
				free(p->x->s[0][0][0]);
				free(p->x->s[0][0]);
				free(p->x->s[0]);
				free(p->x->s);
				break;
			default:
				fprintf(ERRFILE,"ERROR:free_x: can't get here!\n");
				exit(-1);
		}
		free(p->x);
	}
	
	p->x = NULL;
	
} /* free_x */


/************************************************************
 * init_coal_alloc() callocs and initializes a range field. */
void init_coal_alloc(long **coal, long numcoalpairs)
{
	
	if (*coal != NULL) free(*coal);
	
	(*coal) = (long *)calloc(numcoalpairs * 2 + 2,sizeof(long));
	(*coal)[0] = numcoalpairs;
	(*coal)[numcoalpairs*2+1] = FLAGLONG;
	
} /* init_coal_alloc */


/****************************************************************
 * coal_Malloc() callocs or frees the "coal" field of a nodelet */
void coal_Malloc(node *p, boolean allokate, long numcoalpairs)
{
	
	if (allokate){// && p->top) {
		init_coal_alloc(&(p->coal),numcoalpairs);
	} else {
		if (p->coal != NULL) {
			free(p->coal);
			p->coal = NULL;
		}
	}
	
} /* coal_Malloc */


/**************************************************************
 * init_ranges_alloc() callocs and initializes a range field. */
void init_ranges_alloc(long **ranges, long numrangepairs)
{
	
	if (*ranges != NULL) free(*ranges);
	
	(*ranges) = (long *)calloc(numrangepairs * 2 + 2,sizeof(long));
	(*ranges)[0] = numrangepairs;
	(*ranges)[numrangepairs*2+1] = FLAGLONG;
	
} /* init_ranges_alloc */


/********************************************************************
 * ranges_Malloc() callocs or frees the "ranges" field of a nodelet */
void ranges_Malloc(node *p, boolean allokate, long numrangepairs)
{
	
	if (allokate){// && p->top) {
		init_ranges_alloc(&(p->ranges),numrangepairs);
	} else {
		if (p->ranges != NULL) {
			free(p->ranges);
			p->ranges = NULL;
		}
	}
	
} /* ranges_Malloc */


/****************************************************************
 * meld_adjacent_ranges() concatenates elements adjacent to the *
 * passed position if it is appropiate to do so.                */
void meld_adjacent_ranges(long *cranges, long newelem)
{
	long numpairs, nummove, prevelem, nextelem, lastelem;
	
	prevelem = newelem - 1;
	nextelem = newelem + 2;
	lastelem = 2*cranges[0]+1;
	
	numpairs = 0;
	if (cranges[newelem] == cranges[prevelem]+1) numpairs++;
	if (cranges[newelem+1]+1 == cranges[nextelem]) {
		if (numpairs) numpairs++;
		else numpairs = -1L;
	}
	if (numpairs == 2 && prevelem == 0) numpairs = -1L;
	if (numpairs == 2 && nextelem == lastelem) numpairs--;
	
	if (numpairs == 1 && prevelem != 0) {
		nummove = 2*cranges[0] - newelem;
		cranges[prevelem] = cranges[newelem+1];
		memmove(&cranges[newelem],&cranges[nextelem],nummove*sizeof(long));
		cranges[0]--;
	}
	if (numpairs == -1 && nextelem != lastelem) {
		nummove = 2*cranges[0] - newelem - 2;
		cranges[newelem+1] = cranges[nextelem+1];
		memmove(&cranges[nextelem],&cranges[nextelem+2],nummove*sizeof(long));
		cranges[0]--;
	}
	if (numpairs == 2) {
		nummove = 2*cranges[0] - newelem - 2;
		cranges[newelem-1] = cranges[nextelem+1];
		memmove(&cranges[newelem],&cranges[nextelem+2],nummove*sizeof(long));
		cranges[0] -= 2;
	}
	
} /* meld_adjacent_ranges */


/****************************************************
 * addrange() adds a new element to the range field */
void addrange(long **nnewranges, long newstart, long newend)
{
	long i, first, numpairs, nummove, *newranges;
	
	if (newstart > newend) ;
	//fprintf(ERRFILE,"ERROR:addrange--start > end\n");
	
	newranges = (*nnewranges);
	numpairs = newranges[0]*2;
	
	/* was there nothing in the initial range list? */
	if(!numpairs) {
		newranges[0]++;
		newranges = (long *)realloc(newranges,4*sizeof(long));
		(*nnewranges) = newranges;
		newranges[1] = newstart;
		newranges[2] = newend;
		newranges[3] = FLAGLONG;
		return;
	}
	
	/* do I start after everything else ends? */
	if (newstart > newranges[numpairs]) {
		if (newstart == newranges[numpairs]+1) {
			newranges[numpairs] = newend;
		} else {
			newranges[0]++; numpairs += 2;
			newranges = (long *)realloc(newranges,(numpairs+2)*sizeof(long));
			(*nnewranges) = newranges;
			newranges[numpairs-1] = newstart;
			newranges[numpairs] = newend;
			newranges[numpairs+1] = FLAGLONG;
		}
		return;
	}
	
	/* do I end before anything begins? */
	if (newend < newranges[1]) {
		if (newend == newranges[1]-1) {
			newranges[1] = newstart;
		} else {
			newranges[0]++;
			newranges = 
			(long *)realloc(newranges,(newranges[0]*2+2)*sizeof(long));
			(*nnewranges) = newranges;
			memmove(&newranges[3],&newranges[1],(numpairs+1)*sizeof(long));
			newranges[1] = newstart;
			newranges[2] = newend;
		}
		return;
	}
	
	/* am I contained in another interval altogether? OR
     do I overlap one or more intervals?  OR
     am I between two already existing intervals? */
	for(i = 1, numpairs = 0, first = 0; newranges[i] != FLAGLONG; i+=2) {
		if(newstart >= newranges[i] && newend <= newranges[i+1]) return;
		if (i != 1) {
			if (newstart > newranges[i-1] && newend < newranges[i]) {
				first = i;
				break;
			}
		}
		if((newranges[i] >= newstart && newranges[i] <= newend) ||
		   (newranges[i+1] >= newstart && newranges[i+1] <= newend)) {
			numpairs++;
			if (!first) first = i;
		}
	}
	
	/* do I exist between two intervals? */
	if (!numpairs) {
		nummove = 2*newranges[0]+2 - first;
		newranges[0]++;
		newranges = 
		(long *)realloc(newranges,(newranges[0]*2+2)*sizeof(long));
		(*nnewranges) = newranges;
		memmove(&newranges[first+2],&newranges[first],nummove*sizeof(long));
		newranges[first] = newstart;
		newranges[first+1] = newend;
		meld_adjacent_ranges(newranges,first);
		return;
	}
	
	/* I must overlap one or more intervals */
	if (numpairs == 1) {
		newranges[first] = 
		(newranges[first] < newstart) ? newranges[first] : newstart;
		newranges[first+1] = 
		(newranges[first+1] > newend) ? newranges[first+1] : newend;
	} else {
		nummove = 2*newranges[0]+2 - first - 2*numpairs;
		newranges[0] -= numpairs-1;
		newranges[first] = 
		(newranges[first] < newstart) ? newranges[first] : newstart;
		newranges[first+1] =
		(newranges[first+1+2*(numpairs-1)] > newend) ? 
		newranges[first+1+2*(numpairs-1)] : newend;
		memmove(&newranges[first+2],&newranges[first+2*numpairs],
				nummove*sizeof(long));
	}
	meld_adjacent_ranges(newranges,first);
	
} /* addrange */


/******************************************************************
 * inrange() checks to see if a site is active on a given branch */
boolean inrange(long *ranges, long site)
{
	long i;
	
	for(i = 1; ranges[i] != FLAGLONG; i+=2) {
		if (ranges[i] > site) return(FALSE);
		if (ranges[i+1] >= site) return(TRUE);
	}
	
	return(FALSE);
	
} /* inrange */


/************************************************************
 * subrangefc() causes the passed range to be set to "dead" */
void subrangefc(long **newranges, long substart, long subend)
{
	long i, *nranges, numpairs, first, nummove;
	
	nranges = (*newranges);
	
	/* was there nothing in the initial range list? */
	/* do I start after everything else ends? */
	/* do I end before anything begins? */
	if (!nranges[0] || substart > nranges[2*nranges[0]] ||
		subend < nranges[1]) return;
	
	/* am I contained in another interval altogether? OR
     do I overlap one or more intervals?  OR
     am I between two already existing intervals? */
	for(i = 1, numpairs = 0, first = 0; nranges[i] != FLAGLONG; i+=2) {
		if(nranges[i] > subend) break;
		if (i != 1) {
			if (substart > nranges[i-1] && subend < nranges[i])
				return;
		}
		
		/*    if the deletion starts in this range OR
		 the deletion ends in this range OR
		 the deletion overlaps this range THEN */
		if((substart >= nranges[i] && substart <= nranges[i+1]) ||
		   (subend >= nranges[i] && subend <= nranges[i+1]) ||
		   (substart <= nranges[i] && subend >= nranges[i+1])) {
			numpairs++;
			if (!first) first = i;
		}
	}
	
	/* am I contained within a single interval? */
	if (numpairs == 1) {
		/*    should I just remove the whole entry? */
		if ((substart <= nranges[first]) && (subend >= nranges[first+1])) {
			nummove = 2*nranges[0] - first;
			memmove(&nranges[first],&nranges[first+2],nummove*sizeof(long));
			nranges[0]--;
			return;
		}
		/*    or chop off the beginning? */
		if (substart <= nranges[first]) {
			nranges[first] = subend+1;
			return;
		}
		/*    or chop off the end? */
		if (subend >= nranges[first+1]) {
			nranges[first+1] = substart-1;
			return;
		}
		/*    then split it in two! */
		nranges[0]++;
		nranges = 
		(long *)realloc(nranges,(2*nranges[0]+2)*sizeof(long));
		(*newranges) = nranges;
		nummove = 2*nranges[0] - (first+1);
		memmove(&nranges[first+3],&nranges[first+1],nummove*sizeof(long));
		nranges[first+1] = substart-1;
		nranges[first+2] = subend+1;
		return;
	}
	
	/* I must overlap one or more intervals */
	if (substart > nranges[first]) nranges[first+1] = substart-1;
	if (subend < nranges[first+1+2*(numpairs-1)]) 
		nranges[first+2*(numpairs-1)] = subend + 1;
	else {
		/*   remove the last affected rangepair */
		nummove = 2*nranges[0] - (first+2*(numpairs-1));
		memmove(&nranges[first+2*(numpairs-1)],&nranges[first+2*(numpairs)],
				nummove*sizeof(long));
		nranges[0]--;
	}
	/*  remove the middle affected rangepairs, if any */
	nummove = 2*nranges[0]+2 - (first+2*(numpairs-1));
	memmove(&nranges[first+2],&nranges[first+2*(numpairs-1)],
			nummove*sizeof(long));
	nranges[0] -= numpairs-2;
	/*  remove the first affected rangepair, if necessary */
	if (substart <= nranges[first]) {
		nummove = 2*nranges[0] - first;
		memmove(&nranges[first],&nranges[first+2],nummove*sizeof(long));
		nranges[0]--;
	}
	
} /* subrangefc */


void printrange(long *ranges)
{
	long i;
	
	printf("\n%ld range pairs:  ",ranges[0]);
	for(i = 1; ranges[i] != FLAGLONG; i+=2)
		printf("%ld to %ld ",ranges[i],ranges[i+1]);
	printf("\n");
	
} /* printrange */


/* "newnode" & "freenode" are paired memory managers for tree nodes.
 They maintain a linked list of unused nodes which should be faster
 to use than "calloc" & "free" (at least for recombination).
 They can only be used for internal nodes, not tips. */
void newnode(node **p)
{
	long i;
	node *q;
	
	if (freenodes == NULL) { /* need a new node */
		(*p) = allocate_nodelet(3L,'i');
		(*p)->number = nodectr;
		for(q = (*p)->next; q != (*p); q = q->next) q->number = nodectr;
		nodectr++;
		return;
	} else { /* recycle an old node */
		q=freenodes;
		freenodes=freenodes->back;
		for(i = 0; i < 3; i++) {
			q->ranges = NULL;
			//ranges_Malloc(q,FALSE,0L);
			coal_Malloc(q,FALSE,0L);
			q->updated = FALSE;
			q->futileflag = 0L;
			q->recstart = q->recend = -1L;
			q->back = NULL;
			q = q->next;
			q->branch = NULL;
		}
		*p = q;
		return;
	}
}  /* newnode */

void freenode(node *p)
{
	p->back = freenodes;
	freenodes = p;
}  /* freenode */
/* END of treenode allocation functions */

void newtymenode(tlist **t)
{
	*t = (tlist *)malloc(sizeof(tlist));
	(*t)->prev = NULL;
	(*t)->succ = NULL;
	(*t)->segments = NULL;
	
}  /* newtymenode */

void freetymenode(tlist *t)
{
	if(t->branchlist!=NULL) free(t->branchlist);
	if (t->segments != NULL) freesegment(t->segments);
	
	free(t);
}  /* freetymenode*/

void freetymelist(tlist *t)
{
	if (t->succ != NULL) freetymelist(t->succ);
	freetymenode(t);
	
} /* freetymelist */

void hookup(node *p, node *q)
{
	p->back = q;
	q->back = p;
}  /* hookup */

void atr(node *p) 
/* "atr" prints out a text representation of a tree.  Pass 
 curtree->root->back for normal results.  This is a debugging
 function which is not normally called. */
{
	long i;
	node *q;
	
	if (p->back == curtree->root) {
		fprintf(ERRFILE,"next node is root\n");
		fprintf(ERRFILE,"Node %4ld length %12.6f tyme %10.8f",
				p->back->number, lengthof(p), p->back->tyme);
		fprintf(ERRFILE," --> %4ld\n",p->number);
	}
	fprintf(ERRFILE,"Node %4ld update %ld length %10.8f tyme %10.8f -->",
			p->number, (long)p->updated, lengthof(p), p->tyme);
	if (!istip(p))
		for(i = 0, q = p; i < 3; i++, q = q->next)
			if (!q->top) fprintf(ERRFILE,"%4ld\n",q->back->number);
	fprintf(ERRFILE,"\n");
	if (p->top && p->back->top) fprintf(ERRFILE,"TWO TOPS HERE!!!!");
	if (!istip(p)) {
		if (!p->next->top) atr(p->next->back);
		if (!p->next->next->top) atr(p->next->next->back);
	}
	else fprintf(ERRFILE,"\n");
} /* atr */

void probatr(node *p) 
/* "probatr" checks the internal representation of the tree.  Pass 
 curtree->root->back for normal results.  This is a debugging
 function which is not normally called. */
{
	long i;
	
	if (p->back == curtree->root) {
	}
	
	if(p->top) {
		if(p->ranges[0] < 0)
			fprintf(ERRFILE,"node %ld has less then 0 ranges %ld %ld\n",
					p->number,indecks,apps);
		if(p->ranges[0] > 2)
			fprintf(ERRFILE,"node %ld has %ld ranges %ld %ld\n",
					p->number,p->ranges[0],indecks,apps);
		if(p->ranges[2*p->ranges[0]+1] != FLAGLONG)
			fprintf(ERRFILE,"node %ld has misformed end of ranges %ld %ld\n",
					p->number,indecks,apps);
		for(i = 1; p->ranges[i] != FLAGLONG; i+=2) {
			if(p->ranges[i] > p->ranges[i+1])
				fprintf(ERRFILE,"node %ld has begins > ends %ld %ld\n",
						p->number,indecks,apps);
			if(p->ranges[i] < 0 || p->ranges[i+1] < 0)
				fprintf(ERRFILE,"node %ld has negative starts or ends %ld %ld\n",
						p->number,indecks,apps);
		}
	}
	
	if (p->top && p->back->top) fprintf(ERRFILE,"TWO TOPS HERE!!!!");
	if (!istip(p)) {
		if (!p->next->top) probatr(p->next->back);
		if (!p->next->next->top) probatr(p->next->next->back);
	}
	
	
} /* probatr */


/***************************************************************
 * gettymenode() returns a pointer to the tymelist entry whose *
 * 'eventnode' has the number of 'target'.                     */
tlist *gettymenode(tree *tr, long target)
{
	boolean found;
	tlist *t;
	
	if (target == tr->root->number) return (NULL);
	if (tr->nodep[target]->type == 't') return(tr->tymelist); 
	
	for(t = tr->tymelist, found = FALSE; t != NULL; t = t->succ)
		if (t->eventnode->number == target) {found = TRUE; break;}
	
	if (!found) {
		fprintf(ERRFILE,"ERROR:In gettymenode, failed to find node %12ld ", target);
		fprintf(ERRFILE,"CATASTROPHIC ERROR\n");
		exit(-1);
	}
	
	return(t);
	
}  /* gettymenode */


/***************************************************************
 * vtol() takes a "v" value for a branchlength and returns the *
 * branchlength.                                               */
double vtol(option_struct *op, data_fmt *data, double v)
{
	
	switch(op->datatype) {
		case 'a':
			break;
		case 'b':
		case 'm':
			return(v);
			break;
		case 'n':
		case 's':
			return(data->dnaptr->fracchange * -1.0 * v);
			break;
		default:
			fprintf(ERRFILE,"ERROR:vtol, can't get here!\n");
			exit(-1);
	}
	
	return(-1.0);
	
} /* vtol */


/******************************************************************
 * ltov() recalculates the proper "v" value of a branch, from the *
 * tymes at either end of the branch.                             */
void ltov(option_struct *op, data_fmt *data, node *p)
{
	
	switch(op->datatype) {
		case 'a':
			break;
		case 'b':
		case 'm':
			p->v = lengthof(p);
			break;
		case 'n':
		case 's':
			p->v = -1.0 * lengthof(p) / data->dnaptr->fracchange;
			break;
		default:
			fprintf(ERRFILE,"ERROR:ltov, can't get here!\n");
			exit(-1);
	}
	p->back->v = p->v;
	
}  /* ltov */


/***************************************************************
 * findcoal_ltov() is a duplicate ltov specially for findcoal. */
double findcoal_ltov(option_struct *op, data_fmt *data, double value)
{
	
	switch(op->datatype) {
		case 'a':
			break;
		case 'b':
		case 'm':
			return(value);
			break;
		case 'n':
		case 's':
			return(-1.0 * value / data->dnaptr->fracchange);
			break;
		default:
			fprintf(ERRFILE,"ERROR:ltov, can't get here!\n");
			exit(-1);
	}
	
	return(-1.0);
	
} /* findcoal_ltov */


/* boolcheck(), booleancheck(), numbercheck(), and readparmfile() 
 are used in reading the parameter file "parmfile" */
long boolcheck(char ch)
{
	ch = toupper((int)ch);
	if (ch == 'F' || ch == 'N') return 0;
	if (ch == 'T' || ch == 'Y') return 1;
	return -1;
} /* boolcheck */

boolean booleancheck(option_struct *op, char *var, char *value)
{
	long i, j, check;
	
	check = boolcheck(value[0]);
	if(check == -1) return FALSE;
	
	for(i = 0; i < NUMBOOL; i++) {
		if(!strcmp(var,booltokens[i])) {
			if(i == 0) op->interleaved = (boolean)(check);
			if(i == 1) op->printdata = (boolean)(check);
			if(i == 2) op->progress = (boolean)(check);
			if(i == 3) op->treeprint = (boolean)(check);
			if(i == 4) {
				op->freqsfrom = (boolean)(check);
				if(!op->freqsfrom) {
					strtok(value,":");
					freqa = (double)atof(strtok(NULL,";"));
					freqc = (double)atof(strtok(NULL,";"));
					freqg = (double)atof(strtok(NULL,";"));
					freqt = (double)atof(strtok(NULL,";"));
				}
			}
			if(i == 5) {
				op->ctgry = (boolean)(check);
				if(op->ctgry) {
					strtok(value,":");
					op->categs = (long)atof(strtok(NULL,";"));
					op->rate =
					(double *)realloc(op->rate,op->categs*sizeof(double));
					op->probcat =
					(double *)realloc(op->probcat,op->categs*sizeof(double));
					for(j = 0; j < op->categs; j++) {
						op->rate[j] = (double)atof(strtok(NULL,";"));
						op->probcat[j] = (double)atof(strtok(NULL,";"));
					}
				}
			}
			if(i == 6) {
				op->watt = (boolean)(check);
				if (!op->watt) {
					strtok(value,":");
					theta0 = (double)atof(strtok(NULL,";"));
				}
			}
			if(i == 7) op->usertree = (boolean)(check);
			if(i == 8) {
				op->autocorr = (boolean)(check);
				if (op->autocorr) {
					strtok(value,":");
					op->lambda = 1.0 / (double)atof(strtok(NULL,";"));
				}
			}
			if(i == 9) op->newdata = (boolean)(check);
			if(i == 10) {
				op->same_ne = (boolean)(check);
				if (!op->same_ne) {
					strtok(value,":");
					check = atol(strtok(NULL,";"));
					op->ne_ratio = (double *)calloc(check,sizeof(double));
					for(j = 0; j < check; j++)
						op->ne_ratio[j] = (double)atof(strtok(NULL,";"));         
				}
			}
			if(i == 11) op->interactive = (boolean)(check);
			if(i == 12) op->mhmcsave = (boolean)(check);
			if(i == 13) {
				op->panel = (boolean)(check);
				if (op->panel) {
					strtok(value,":");
					check = atol(strtok(NULL,";"));
					op->numpanel = (long *)calloc(check,sizeof(long));
					for (j = 0; j < check; j++) 
						op->numpanel[j] = atol(strtok(NULL,";"));
				}
			}
			if(i == 14) op->map = (boolean)(check);
			if(i == 15) op->fc = (boolean)(check);
			if(i == 16) {
				op->full = (boolean)(check);
				if (op->full) {
					strtok(value,":");
					op->chance_seen = (double)atof(strtok(NULL,";"));
				}
			}
			if(i == 17) {
				op->haplotyping = (boolean)(check);
				if (op->haplotyping) {
					strtok(value,":");
					op->happrob = (double)atof(strtok(NULL,";"));
				}
			}
			if(i == 18) op->profile = (boolean)(check);
			if(i == 19) op->norecsnp = (boolean)(check);
			if (i == 20){
				if ((boolean)(check)) op->userrec = 1;
				else op->userrec = 0;
			}
			if (i == 21) menuoption = (boolean)(check);
			return TRUE;
		}
	}
	return FALSE;
} /* booleancheck */

boolean numbercheck(option_struct *op, char *var, char *value)
{
	long i, j;
	
	for(i = 0; i < NUMNUMBER; i++) {
		if(!strcmp(var,numbertokens[i])) {
			if(i == 0) locus_ttratio = (double)atof(value);
			if(i == 1) op->numchains[0] = atol(value);
			if(i == 2) op->steps[0] = atol(value);
			if(i == 3) op->increm[0] = atol(value);
			if(i == 4) op->numchains[1] = atol(value);
			if(i == 5) op->increm[1] = atol(value);
			if(i == 6) op->steps[1] = atol(value);
			if(i == 7) rec0 = (double)atof(value);
			if(i == 8) {
				value[0] = toupper((int)value[0]);
				switch(value[0]) {
					case 'N':
						op->holding = 0;
						break;
					case 'T':
						op->holding = 1;
						break;
					case 'R':
						op->holding = 2;
						break;
					default:
						fprintf(ERRFILE,"unknown setting for holding: %s,",
								value);
						fprintf(ERRFILE,"  no holding assumed\n");
						op->holding = 0;
						break;
				}
			}
			if(i == 9) op->mutrait = (double)atof(value);
			if(i == 10) op->traitratio = (double)atof(value);
			if(i == 11) op->pd = (double)atof(value);
			if(i == 12) op->hapdrop = atol(value);
			if(i == 13) {
				op->numtempchains = atol(strtok(value,";"));
				if (!op->temperature) free(op->temperature);
				op->temperature = (long *)calloc(op->numtempchains,sizeof(long));
				for(j = 0; j < op->numtempchains; j++)
					op->temperature[j] = atof(strtok(NULL,";"));     
				qsort((void *)(op->temperature),op->numtempchains,sizeof(long),longcmp);
			}
			if(i == 14) recrates[0] = atof(value);
			if(i == 15) recrates[1] = atof(value);
			if(i == 16) lamda[0] = atof(value);
			if(i == 17) lamda[1] = atof(value);
			return TRUE;
		}
	}
	return FALSE;
	
} /* numbercheck */

void readparmfile(option_struct *op)
{
	char fileline[LINESIZE],parmvar[LINESIZE],varvalue[LINESIZE];
	
#ifdef MAC
	char parmfilename[100];
	
	strcpy(parmfilename,"parmfile");
#endif
	
	parmfile = fopen("parmfile","r");
	
	if(parmfile) {
		while(fgets(fileline, LINESIZE, parmfile) != NULL) {
			if(fileline[0] == '#') continue;
			if(!strncmp(fileline,"end",3)) break;
			strcpy(parmvar,strtok(fileline,"="));
			strcpy(varvalue,strtok(NULL,"\n"));
			/* now to process... */
			if(!booleancheck(op,parmvar,varvalue))
				if(!numbercheck(op,parmvar,varvalue)) {
					fprintf(ERRFILE,
							"Inappropiate entry in parmfile: %s\n", fileline);
				}
		}
	} else
		if (!MENU) {
			fprintf(simlog,"Parameter file (parmfile) missing\n");
			exit(-1);
		}
	
	FClose(parmfile);
	
#ifdef MAC
	fixmacfile(parmfilename);
#endif
	
} /* readparmfile */
/* END parameter file read */


boolean whichopbool(option_struct *op, long i)
{
	if (i == 0) return(op->interleaved);
	if (i == 1) return(op->printdata);
	if (i == 2) return(op->progress);
	if (i == 3) return(op->treeprint);
	if (i == 4) return(op->freqsfrom);
	if (i == 5) return(op->ctgry);
	if (i == 6) return(op->watt);
	if (i == 7) return(op->usertree);
	if (i == 8) return(op->autocorr);
	if (i == 9) return(op->newdata);
	if (i == 10) return(op->same_ne);
	if (i == 11) return(op->interactive);
	if (i == 12) return(op->mhmcsave);
	if (i == 13) return(op->panel);
	if (i == 14) return(op->map);
	if (i == 15) return(op->fc);
	if (i == 16) return(op->full);
	if (i == 17) return(op->haplotyping);
	if (i == 18) return(op->profile);
	if (i == 19) return(op->norecsnp);
	
	fprintf(ERRFILE,"Warning:  whichopbool failed to identify option\n");
	return(FALSE);
} /* whichopbool */


void rec_parmfilewrite(option_struct *op, long numloci, long numpop)
{
	long i, j;
	boolean b;
	char parmfilename[100];
	
	
	openfile(&parmfile,"parmfile","w+",NULL,parmfilename);
	
	for(i = 0; i < NUMBOOL; i++) {
		b = whichopbool(op,i);
		fprintf(parmfile,"%s=%s",booltokens[i],BOOLPRINT(b));
		if (i == 4 && !b) { /* freqsfrom */
			fprintf(parmfile,":%f;%f;%f;%f;",
					freqa, freqc, freqg, freqt);
		}
		if (i == 5 && b) { /* category */
			fprintf(parmfile,":%ld;",op->categs);
			for(j = 0; j < op->categs; j++)
				fprintf(parmfile,"%f;%f;",op->rate[j],op->probcat[j]);
		}
		if (i == 6 && !b) { /* starting theta */
			fprintf(parmfile,":%f",theta0);
		}
		if (i == 8 && !b) { /* autocorrelation */
			fprintf(parmfile,":%f",1.0/op->lambda);
		}
		if (i == 10 && !b) { /* same_ne */
			fprintf(parmfile,":%ld;",numloci);
			for(j = 0; j < numloci; j++)
				fprintf(parmfile,"%f;",op->ne_ratio[j]);
		}
		if (i == 13 && b) { /* numpanels */
			fprintf(parmfile,":%ld;",numpop);
			for(j = 0; j < numpop; j++)
				fprintf(parmfile,"%ld;",op->numpanel[j]);
		}
		if (i == 16 && b) { /* full snp */
			fprintf(parmfile,":%f;",op->chance_seen);
		}
		fprintf(parmfile,"\n");
	}
	
	for(i = 0; i < NUMNUMBER; i++) {
		fprintf(parmfile,"%s=",numbertokens[i]);
		if (i == 0) fprintf(parmfile,"%f",locus_ttratio);
		if (i == 1) fprintf(parmfile,"%ld",op->numchains[0]);
		if (i == 2) fprintf(parmfile,"%ld",op->steps[0]);
		if (i == 3) fprintf(parmfile,"%ld",op->increm[0]);
		if (i == 4) fprintf(parmfile,"%ld",op->numchains[1]);
		if (i == 5) fprintf(parmfile,"%ld",op->increm[1]);
		if (i == 6) fprintf(parmfile,"%ld",op->steps[1]);
		if (i == 7) fprintf(parmfile,"%f",rec0);
		if (i == 8) {
			if (!op->holding) fprintf(parmfile,"none");
			if (op->holding == 1) fprintf(parmfile,"theta");
			if (op->holding == 2) fprintf(parmfile,"recombination");
		}
		if (i == 9) fprintf(parmfile,"%f",op->mutrait);
		if (i == 10) fprintf(parmfile,"%f",op->traitratio);
		if (i == 11) fprintf(parmfile,"%f",op->pd);
		if (i == 12) fprintf(parmfile,"%ld",op->hapdrop);
		if (i == 13) {
			fprintf(parmfile,"%ld;",op->numtempchains);
			for(j = 0; j < op->numtempchains; j++)
				fprintf(parmfile,"%ld;",op->temperature[j]);
		}
		fprintf(parmfile,"\n");
	}
	
	fprintf(parmfile,"end\n");
	
	FClose(parmfile);
	
#ifdef MAC
	fixmacfile(parmfilename);
#endif
	
} /* rec_parmfilewrite */


void readseedfile(void)
{
	long inseed;
	
	seedfile = fopen("seedfile","r");
	
	if (!MENU) {
		if (!seedfile) {
			fprintf(ERRFILE,"\nseedfile not present!\n");
			exit(-1);
		}
		fscanf(seedfile, "%ld%*[^\n]", &inseed);
		getc(seedfile);
	} else {
		if (seedfile) {
			fscanf(seedfile, "%ld%*[^\n]", &inseed);
			getc(seedfile);
		} else {
			printf("Random number seed (must be odd)?\n");
			scanf("%ld%*[^\n]", &inseed);
			getchar();
		}
	}
	fclose(seedfile);
	//seedfile = fopen("seedfile","w");
	//fprintf(seedfile,"%ld\n",inseed + 4);
	//fclose(seedfile);
	
	seed[0] = inseed & 2047;
	inseed /= 2048;
	seed[1] = inseed & 2047;
	inseed /= 2048;
	seed[2] = inseed & 1023;
	
} /* readseedfile */


void print_menuheader(option_struct *op)
{
	printf("\n%s", op->ansi ? "\033[2J\033[H" : "\n");
	printf("Metropolis-Hastings Markov Chain Monte Carlo");
	printf(" method, version %3.2f\n\n",VERSION_NUM);
} /* print_menuheader */


void print_startmenu(option_struct *op, boolean writeout)
{
	printf("STARTUP MENU\n");
	printf("  #               Goto Data/Search Menus\n");
	printf("  O         Save current options to file?  %s\n",
		   writeout ? "Yes" : "No");
	printf("  N          Use trees from previous run?  %s\n",
		   op->newdata ? "No" : "Yes");
	if (op->newdata) 
		printf("  E        Echo the data at start of run?  %s\n",
			   op->printdata ? "Yes" : "No");
	else op->printdata = FALSE;
	/*   printf("  S            Save MHMC output to files?  %s\n", */
	/* 	 op->mhmcsave ? "Yes" : "No"); */
	printf("  P Print indications of progress of run?  %s\n",
		   op->progress ? "Yes" : "No");
#if 0 /* DEBUG debug WARNING warning--no tree writer yet! */
	printf("  G                Print out genealogies?  %s\n",
		   op->treeprint ? "Yes" : "No");
#endif
	printf("  U      Use user tree in file \"intree\" ?  %s\n",
		   op->usertree ? "Yes" : "No");
	printf("  V          Number of temperatures used?  %ld\n",
		   op->numtempchains);
	/*   printf("  H                     Infer haplotypes?  %s\n", */
	/* 	 op->haplotyping ? "Yes" : "No"); */
	/*   printf("  L       Calculate confidence intervals?  %s\n", */
	/* 	 op->profile ? "Yes" : "No"); */
	if (op->haplotyping) {
		printf("  A Strategy for haplotype rearrangement?  ");
		if(op->hapdrop==0)printf("No resimulation\n");
		else if (op->hapdrop==1)printf("Single resimulation\n");
		else printf("Double resimulation\n");
		printf("  F Frequency of haplotype rearrangement?  %f\n",
			   op->happrob);
	}
	/*   printf("  M     Map trait information onto trees?  %s\n", */
	/* 	 op->map ? "Yes" : "No"); */
	if (op->map) {
		printf("  R        Relative mutation rate of trait?  %f\n",
			   op->mutrait);
		printf("  T     Forward versus back trait mutation?  %f\n",
			   op->traitratio);
		printf("  D                     Frequency of trait?  %f\n",
			   op->pd);
	}
	
} /* print_startmenu */
void print_datamenu(option_struct *op)
{
	printf("DATA MENU\n");
	printf("  #                     Goto Startup Menu\n");
	/*   printf("  D                             Datatype:"); */
	op->datatype = 's';
	switch(op->datatype) {
		case 'a':
			printf("  Allelelic Markers\n");
			break;
		case 'b':
			printf("  Brownian-Motion Microsats\n");
		case 'm':
			if (op->datatype != 'b') printf("  Microsats\n");
			break;
		case 'n':
			printf("  SNPs\n");
			printf("  P              SNPs derived from panel?  %s\n",
				   op->panel ? "Yes" : "No");
#if 0
			printf("  A   SNPs are all of the variable sites?  %s\n",
				   op->full ? "No" : "Yes");
			if (op->full) {
				printf("  B           Probability of observation?  %6.4f\n",
					   op->chance_seen);
			}
#endif
			printf("  N              SNPs with recombination?  %s\n",
				   op->norecsnp ? "No" : "Yes");
		case 's':
			/*     if (op->datatype != 'n') printf("  Sequence\n"); */
			printf("  I          Input sequences interleaved?  %s\n",
				   op->interleaved ? "Yes" : "No, sequential");
			printf("  T        Transition/transversion ratio:  %8.4f\n",
				   locus_ttratio);
			printf("  F       Use empirical base frequencies?  %s\n",
				   op->freqsfrom ? "Yes" : "No");
			/*     printf("  C   One category of substitution rates?"); */
			/*     if (!op->ctgry || op->categs == 1) */
			/*       printf("  Yes\n"); */
			/*     else { */
			/*       printf("  %ld categories\n", op->categs); */
			/*       if (op->datatype != 'n') { */
			/* 	printf("  R   Rates at adjacent sites correlated?"); */
			/* 	if (!op->autocorr) printf("  No, they are independent\n"); */
			/* 	else printf("  Yes, mean block length =%6.1f\n", 1.0 / op->lambda); */
			/*       } */
			/*     } */
			/*     printf("  V    Poplulation size equal among loci?  %s\n", */
			/* 	   op->same_ne ? "Yes" : "No"); */
			break;
		default:
			printf("  UNKNOWN\n");
			break;
	}
	
} /* print_datamenu */


void print_searchmenu(option_struct *op)
{
	int i,j = 0;
	printf("Recombination hotspots Menu\n");
	printf("  A   Starting recombination rate on cold region:  %lf\n", recrates[0]);
	printf("  B       Starting recombination rate on hotspot:  %lf\n", recrates[1]);
	printf("  J        Transition paramter from cold to cold:  %lf\n",lamda[0]);
	printf("  K          Transition paramter from hot to hot:  %lf\n",lamda[1]);
	//  printf("  1   
	
	printf("  U         Use user recombination rates pattern? ");
	if (op->userrec == 0) printf(" No\n");
	else printf(" Yes\n");
	
	printf("SEARCH MENU\n");
	/*   printf("  Q      Lots of recombinations expected?  %s\n", */
	/* 	 op->fc ? "Yes" : "No"); */
	printf("  H                Hold parameters fixed?");
	/*   if (op->holding == 0) printf("  No\n"); */
	/*   else if (op->holding == 1) printf("  Yes, theta fixed\n"); */
	/*   else if (op->holding == 2) printf("  Yes, rec-rate fixed\n"); */
	/*   else printf("  Unknown option!!!\n"); */
	for (i = 0; i<5; i++){
		if (op->holdings[i] == 1) j = 1;
	}
	if ( j == 0) printf(" None parameter fixed\n");
	else{
		printf("  Yes");
		if (op->holdings[0] == 1) printf(" theta");
		if (op->holdings[1] == 1) printf(" rec_cold");
		if (op->holdings[2] == 1) printf(" rec_hot");
		if (op->holdings[3] == 1) printf(" transition_cold");
		if (op->holdings[4] == 1) printf(" transition_hot");
		
		printf(" fixed\n");
	}
	j = 0;
	
	printf("  W      Starting theta equals Watterson?  %s",
		   op->watt ? "Yes\n" : "No");
	if (!op->watt) printf(", initial theta = %6.4f\n", theta0);
	printf("  S               Number of short chains?  %6ld\n",
		   op->numchains[0]);
	if (op->numchains[0] > 0) {
		printf("  1             Short sampling increment?  %6ld\n",
			   op->increm[0]);
		printf("  2   Number of steps along short chains?  %6ld\n",
			   op->steps[0]);
	}
	printf("  L                Number of long chains?  %6ld\n",
		   op->numchains[1]);
	if (op->numchains[1] > 0) {
		printf("  3              Long sampling increment?  %6ld\n",
			   op->increm[1]);
		printf("  4    Number of steps along long chains?  %6ld\n",
			   op->steps[1]);
	}
	
} /* print_searchmenu */


void print_menuend(void)
{
	printf("\nAre these settings correct? (type Y or the letter for");
	printf(" one to change)\n");
} /* print_menuend */


void initoptions(option_struct *op)
{
	long i;
	
	/* first some multiple rate-categories code stuff */
	op->ctgry = FALSE;
	op->rate[0] = 1.0;
	op->probcat[0] = 1.0;
	op->categs = 1;
	op->lambda = 1.0;
	op->autocorr = FALSE;  /* FALSE if categs == 1 */
	/* end categories code stuff */
	
	op->interactive = TRUE;
	op->newdata = TRUE;
	op->mhmcsave = FALSE;
	op->interleaved = TRUE;
	op->printdata = FALSE;
	op->progress = TRUE;
	op->treeprint = FALSE;
	op->profile = TRUE;
	
	weightfile = fopen("weightfile","r");
	op->weights = (weightfile) ? TRUE : FALSE;
	
	spacefile = fopen("spacefile","r");
	op->spacing = (spacefile) ? TRUE : FALSE;
	
	op->map = FALSE;
	op->mutrait = 1.0;
	op->traitratio = 1.0;
	op->pd = 0.1;
	
	op->panel = FALSE;
	op->numpanel = NULL;
	
	op->haplotyping = FALSE;
	op->hapdrop = 0;
	op->happrob = 0.2;
	
	op->full = FALSE;
	op->chance_seen = 1.0;
	op->norecsnp = FALSE;
	
	op->numtempchains = 1;
	op->temperature = (long *)calloc(op->numtempchains,sizeof(long));
	for(i = 0; i < op->numtempchains; i++) op->temperature[i] = 1;
	op->ctemp = 1;
	op->userrec = 0;
	
	op->fc = TRUE;
	op->freqsfrom = TRUE;
	op->watt = FALSE;
	op->usertree = FALSE;
	op->plump = TRUE;
	op->holding = 0;
	op->same_ne = TRUE;
	op->ne_ratio = NULL;
	op->numchains[0] = 5;
	op->increm[0] = 20;
	op->steps[0] = 10000;
	op->numchains[1] = 1;
	op->increm[1] = 20;
	op->steps[1] = 50000;
	op->datatype = 's';
	op->datatypeset = FALSE;
	
	op->holdings[0] = op->holdings[1] = op->holdings[2] = op->holdings[3] = op->holdings[4] = 0;
	
	op->print_recbythmaxcurve = TRUE;
	op->thlb = FLAGLONG;
	op->thub = FLAGLONG;
	op->reclb = FLAGLONG;
	op->recub = FLAGLONG;
	
} /* initoptions */


void workmenu1(option_struct *op, char ch, boolean *menu1, boolean *writeout)
{
	long i;
	char input[LINESIZE];
	
	if (strchr("#FLAHONESPGUMVRTD",ch) != NULL) {
		switch(ch) {
			case '#':
				*menu1 = !(*menu1);
				break;
			case 'V':
				printf("Number of different temperatures?\n");
				scanf("%ld",&(op->numtempchains));
				//op->numtempchains++;
				scanf("%*[^\n]");
				if (!op->temperature) free(op->temperature);
				op->temperature = (long *)calloc(op->numtempchains,sizeof(long));
				do {
					printf("Temperature of coldest chain? (must be positive)\n");
					scanf("%ld",&(op->temperature[0]));
					scanf("%*[^\n]");
				} while (op->temperature[0] <= 0.0);
				for (i = 1; i < op->numtempchains; i++) {
					do {
						printf("Temperature of temperature-chain %ld?\n",i+1);
						scanf("%ld",&(op->temperature[i]));
						scanf("%*[^\n]");
					} while(op->temperature[i] <= op->temperature[0]);
				}
				break;
			case 'F':
				printf("Probability of proposing a change in haplotype?\n");
				scanf("%lf",&(op->happrob));
				scanf("%*[^\n]");
				break;
			case 'L':
				op->profile = !op->profile;
				break;
			case 'H':
				op->haplotyping = !op->haplotyping;
				break;
			case 'O':
				*writeout = !(*writeout);
				break;
			case 'N':
				op->newdata = !op->newdata;
				break;
			case 'E':
				op->printdata = !op->printdata;
				break;
			case 'P':
				op->progress = !op->progress;
				break;
			case 'G':
				op->treeprint = !op->treeprint;
				break;
			case 'U':
				op->usertree = !op->usertree;
				break;
			case 'S':
				op->mhmcsave = !op->mhmcsave;
				break;
			case 'M':
				op->map = !op->map;
				break;
			case 'R':
				do {
					printf("Relative mutation rate of trait?\n");
					//gets(input);
					fgets(input, LINESIZE, stdin);	
					op->mutrait = atof(input);
				} while (op->mutrait <= 0.0);
				break;
			case 'T':
				do {
					printf("Ratio of forward to back trait mutation?\n");
					//gets(input);
					fgets(input, LINESIZE, stdin);
					op->traitratio = atof(input);
				} while (op->traitratio <= 0.0);
				break;
			case 'D':
				do {
					printf("Frequency of trait?\n");
					//gets(input);
					fgets(input, LINESIZE, stdin);
					op->pd = atof(input);
				} while (op->pd <= 0.0 || op->pd >= 1.0);
				break;
			case 'A':
				do {
					printf("Number of drops while resimulating (0-2)?\n");
					//gets(input);
					
					fgets(input, LINESIZE, stdin);
					op->hapdrop = atol(input);
				} while (op->hapdrop != 0 && op->hapdrop != 1 && op->hapdrop != 2); 
			default:
				fprintf(stderr,"ERROR:Impossible option %c detected!\n",ch);
				break;
		}
	} else fprintf(stderr,"%c is not an option here.\n",ch);
	
} /* workmenu1 */


void workmenu2(option_struct *op, char ch, boolean *menu1, long *numloci,
			   long *numpop)
{
	long i;
	double probsum;
	char input[LINESIZE];
	boolean done;
	
	if(strchr("#NQPDITFCRVHWZS12L34ABJKU",ch) != NULL) {
		switch(ch) {
			case '#':
				*menu1 = !(*menu1);
				break;
			case 'N':
				op->norecsnp = !(op->norecsnp);
				if (op->norecsnp) {
					if (rec0 || op->holding != 2) {
						printf("\nThis SNP model requires recombination");
						printf(" rate to be fixed at zero.\n");
						printf("Fix starting rec-rate to zero first.\n");
						printf("----SNP model unchanged.----\n");
						op->norecsnp = FALSE;
					}
				}
				printf("press <enter> or <return> to continue\n");
				fgets(input,LINESIZE,stdin);
				break;
			case 'Q':
				op->fc = !op->fc;
				break;
			case 'D':
				/* WARNING--most datatypes disabled!!! */
#if 0
				if (op->datatype == 'a') {op->datatype = 'b'; break;}
				if (op->datatype == 'b') {op->datatype = 'm'; break;}
				if (op->datatype == 'm') {op->datatype = 'n'; break;}
				if (op->datatype == 'n') {op->datatype = 's'; break;}
				if (op->datatype == 's') {op->datatype = 'a'; break;}
#endif
				op->datatypeset = TRUE;
				if (op->datatype == 'n') {op->datatype = 's'; break;}
				if (op->datatype == 's') {
					op->datatype = 'n'; 
					if (op->autocorr) {
						printf("Can't use autocorrelation with SNPs\n");
						printf("----autocorrelation disabled-----\n");
						printf("press <enter> or <return> to continue\n");
						fgets(input,LINESIZE,stdin);
						op->autocorr = FALSE;
					}
					break;
				}
				printf("ERROR:Impossible Datatype %c present\n",op->datatype);
				break;
#if 0
			case 'A': /* deliberate fall through to case 'B' */
				op->full = !op->full;
			case 'B':
				if (!op->full) break;
				do {
					printf("Percent chance that a variable site was observed");
					printf(" during sequencing?\n");
					scanf("%lf",&(op->chance_seen));
					scanf("%*[^\n]");
					if (op->chance_seen <= 0.0 || op->chance_seen > 1.0)
						printf("   the percentage must be between 0 and 1\n");
				} while (op->chance_seen <= 0.0 || op->chance_seen > 1.0);
				break;
#endif
				
			case 'I':
				op->interleaved = !op->interleaved;
				break;
			case 'T':
				do {
					printf("Transition/transversion ratio?\n");
					//gets(input);
					fgets(input, LINESIZE, stdin);
					locus_ttratio = atof(input);
					if (locus_ttratio < 0.5)
						printf("TTratio cannot be less than 0.5\n");
				} while (locus_ttratio < 0.5);
				break;
			case 'F':
				op->freqsfrom = !op->freqsfrom;
				if (!op->freqsfrom) {
					printf("Base frequencies for A, C, G, T/U (use blanks");
					printf(" to separate)?\n");
					scanf("%lf%lf%lf%lf", &freqa, &freqc, &freqg, &freqt);
					scanf("%*[^\n]");
				}
				break;
			case 'P':
				op->panel = !op->panel;
				if (op->panel) {
					printf("Number of populations?\n");
					//gets(input);
					
					fgets(input, LINESIZE, stdin);
					*numpop = atol(input);
					op->numpanel = (long *)calloc(*numpop,sizeof(long));
					for(i = 0; i < *numpop; i++) {
						printf("Number of panel haplotypes for population");
						printf(" %ld?\n",i+1);
						//gets(input);
						
						fgets(input, LINESIZE, stdin);
						op->numpanel[i] = atol(input);
					}
				} else
					if (op->numpanel) {
						free(op->numpanel);
						op->numpanel = NULL;
					}
				break;
			case 'C':
				op->ctgry = !op->ctgry;
				if (!op->ctgry) op->autocorr = FALSE;
				if (op->ctgry) {
					do {
						printf("Number of categories ?");
						//gets(input);
						
						fgets(input, LINESIZE, stdin);
						op->categs = atoi(input); 
					} while (op->categs < 1);
					free(op->rate);
					free(op->probcat);
					printf("Rate for each category? (use a space to");
					printf(" separate)\n");
					op->rate = (double *)calloc(op->categs,sizeof(double));
					op->probcat = (double *)calloc(op->categs,sizeof(double));
					for (i = 0; i < op->categs; i++)
						scanf("%lf*[^\n]", &(op->rate[i]));
					getchar();
					do {
						printf("Probability for each category?");
						printf(" (use a space to separate)\n");
						for (i = 0; i < op->categs; i++)
							scanf("%lf", &(op->probcat[i]));
						scanf("%*[^\n]");
						getchar();
						done = TRUE;
						probsum = 0.0;
						for (i = 0; i < op->categs; i++)
							probsum += op->probcat[i];
						if (fabs(1.0 - probsum) > 0.001) {
							done = FALSE;
							printf("Probabilities must add up to");
							printf(" 1.0, plus or minus 0.001.\n");
						}
					} while (!done);
				}
				break;
			case 'R':
				if (op->datatype == 'n') {
					printf("Can't use autocorrelation with SNPs\n");
					printf("----autocorrelation disabled-----\n");
					printf("press <enter> or <return> to continue\n");
					fgets(input,LINESIZE,stdin);
					op->autocorr = FALSE;
					break;
				}
				op->autocorr = !op->autocorr;
				if (op->autocorr) {
					do {
						printf("Mean block length of sites having the same ");
						printf("rate (greater than 1)?\n");
						scanf("%lf%*[^\n]", &(op->lambda));
						getchar();
					} while (op->lambda <= 1.0);
					op->lambda = 1.0 / op->lambda;
				}
				break;
			case 'V':
				op->same_ne = !op->same_ne;
				if (!op->same_ne) {
					printf("Enter number of loci: ");
					scanf("%ld",numloci);
					if (*numloci == 1) op->same_ne = !op->same_ne;
					if (op->ne_ratio) free(op->ne_ratio);
					op->ne_ratio = (double *)calloc(*numloci,sizeof(double));
					printf("\nEnter relative population size/mutation rate");
					printf(" for each locus in input order:\n");
					for(i = 0; i < *numloci; i++) {
						scanf("%lf",&(op->ne_ratio[i]));
						if (op->ne_ratio[i] <= 0) {
							printf("\nratios must be positive, please reenter\n");
							i--;
						}
					}
				}
				break;
				
			case 'A':
				printf("New recombiantion rate on cold region?\n");
				
				//gets(input);
				
				fgets(input, LINESIZE, stdin);
				recrates[0] = atof(input);
				
				break;
				
			case 'B':
				printf("New recombiantion rate on hotspot?\n");
				
				//gets(input);
				
				fgets(input, LINESIZE, stdin);
				recrates[1] = atof(input);
				
				break;   
				
			case 'K':
				printf("New transition rate from hot to hot?\n");
				
				//gets(input);
				
				fgets(input, LINESIZE, stdin);
				lamda[1] = atof(input);
				
				break;     
				
			case 'J':
				printf("New transition rate from cold to cold?\n");
				
				//gets(input);
				
				fgets(input, LINESIZE, stdin);
				lamda[0] = atof(input);
				
				break;
				
			case 'U':
				if (op->userrec ==0)    op->userrec = 1;
				else op->userrec = 0;
				
				break;
				
			case 'H':
				printf("Which parameter? (theta = 0, rec_cold = 1, rec_hot = 2, transition_cold = 3, transition_hot = 4\n)");
				//gets(input);
				
				fgets(input, LINESIZE, stdin);
				i = atoi(input);
				if (op->holdings[i] == 0)  op->holdings[i] = 1;
				else{
					if (op->holdings[i] == 1)  op->holdings[i] = 0;
				}
				
				break;
			case 'W':
				op->watt = !op->watt;
				if (!op->watt) {
					do {
						printf("Initial theta estimate?\n");
						//gets(input);
						
						fgets(input, LINESIZE, stdin);
						theta0 = atof(input);
					} while (theta0 <= 0.0);
				}
				break;
			case 'Z':
				printf("What recombination rate?\n");
				do {
					//gets(input);
					
					fgets(input, LINESIZE, stdin);
					rec0 = atof(input);
					if (rec0 < 0.0)
						printf("recombination rate must be non-negative\n");
				} while (rec0 < 0.0);
				break;
			case 'S':
				do {
					printf("How many Short Chains?\n");
					//gets(input);
					
					fgets(input, LINESIZE, stdin);
					op->numchains[0] = atoi(input);
					if (op->numchains[0] < 0)
						printf("Must be non-negative\n");
				} while (op->numchains[0] < 0); 
				break;
			case '1':
				done = FALSE;
				while (!done) {
					printf("How often to sample trees?\n");
					//gets(input);
					
					fgets(input, LINESIZE, stdin);
					op->increm[0] = atoi(input);
					if (op->increm[0] > 0) done = TRUE;
					else printf("Must be a positive integer\n");
				}
				break;
			case '2':
				done = FALSE;
				while (!done) {
					printf("How many short steps?\n");
					//gets(input);
					
					fgets(input, LINESIZE, stdin);
					op->steps[0] = atoi(input);
					if (op->steps[0] > 0) done = TRUE;
					else printf("Must be a positive integer\n");
				}
				break;
			case 'L':
				do {
					printf("How many Long Chains?\n");
					//gets(input);
					fgets(input, LINESIZE, stdin);
					op->numchains[1] = atoi(input);
					if (op->numchains[1] < 1)
						printf("Must be a positive integer\n");
				} while (op->numchains[1] < 1);
				break;
			case '3':
				done = FALSE;
				while (!done) {
					printf("How often to sample trees?\n");
					//gets(input);
					fgets(input, LINESIZE, stdin);
					op->increm[1] = atoi(input);
					if (op->increm[1] > 0) done = TRUE;
					else printf("Must be a positive integer\n");
				}
				break;
			case '4':
				done = FALSE;
				while (!done) {
					printf("How many long steps?\n");
					//gets(input);
					fgets(input, LINESIZE, stdin);
					op->steps[1] = atoi(input);
					if (op->steps[1] > 0) done = TRUE;
					else printf("Must be a positive integer\n");
				}
				break;
			default:
				fprintf(stderr,"ERROR:Impossible option %c detected!\n",ch);
				break;
		}
	} else fprintf(stderr,"%c is not an option here.\n",ch);
	
} /* workmenu2 */


/****************************************************************
 * checkparmfile() checks the parmfile for internal consistency */
void checkparmfile(option_struct *op)
{
	boolean error = FALSE;
	
	if (op->norecsnp) {
		if (rec0) {
			error = TRUE;
		}
		if (op->holding != 2) {
			error = TRUE;
		}
	}
	
	if (error) {
		exit(-1);
	}
	
} /* checkparmfile */



void getoptions(option_struct *op)
/* interactively set options using a very basic menu */
{
	long numloci, numpop;
	boolean writeout, done, menu1;
	char ch, input[LINESIZE];
	
	/* default initializations */
	initoptions(op);  
	writeout = FALSE;
	locus_ttratio = 2.0;
	numloci = 1;
	numpop = 1;
	rec0 = 0.1;
	/* end defaults */
	
	readparmfile(op);
	checkparmfile(op);
	readseedfile();
	
	if (sim_mode){
		return;
	}
	
	if (!menuoption) return;
	
	menu1 = TRUE;
	do {
		print_menuheader(op);
		if (menu1) print_startmenu(op,writeout);
		else {print_datamenu(op); print_searchmenu(op);}
		print_menuend();
		//gets(input);
		fgets(input, LINESIZE, stdin);
		ch = toupper((int)input[0]);
		done = (ch == 'Y');
		if (!done) {
			if (menu1) workmenu1(op,ch,&menu1,&writeout);
			else workmenu2(op,ch,&menu1,&numloci,&numpop);
		}
	} while (!done);
	
	if (writeout) rec_parmfilewrite(op,numloci,numpop);
	
}  /* getoptions */


/*********************************************************************
 * firstinit() handles initialization for things that are recorded   *
 * over multiple loci/populations, and therefore are allocated once. */
void firstinit(option_struct *op, data_fmt *data)
{
	long i, numloci;
	
	numloci = getdata_numloci(op,data);
	
	totchains = op->numchains[0] + op->numchains[1];
	
	numtrees = MAX(op->steps[0]/op->increm[0],op->steps[1]/op->increm[1]);
	
	if (op->same_ne) {
		op->ne_ratio = (double *)calloc(1,numloci*sizeof(double));
		for(i = 0; i < numloci; i++) op->ne_ratio[i] = 1.0;
	}
	
	sametree = (boolean **)calloc(1,numloci * sizeof(boolean *));
	//sametree[0] = (boolean *)calloc(1,numloci*numtrees * sizeof(boolean));
	//for(i = 1; i < numloci; i++)
	//  sametree[i] = sametree[0] + i*numtrees;
	
	//model_alloc(op,data);
	
	freenodes = NULL;
	
}  /* firstinit */


/********************************************************************
 * popinit() initializes things that are specific to one population */ 
void popinit(option_struct *op, data_fmt *data)
{
	
	switch(op->datatype) {
		case 'a':
			break;
		case 'b':
		case 'm':
			break;
		case 'n':
		case 's':
			rootnum = 0;
			break;
		default:
			fprintf(ERRFILE,"ERROR:popinit, can't get here!\n");
			exit(-1);
	}
	
} /* popinit */


/***********************************************************************
 * makeinvarvalues() sets up the tip likelikelihoods for the invariant *
 * sites tacked on the end of the sequence for SNP data.               */
void makeinvarvalues(dnadata *dna, long numcategs, tree *tr)
{
	long site, ind, categ, slice;
	
	site = dna->sites[locus];
	slice = 0L; /* we only set slice 0 here; the remainder, if any,
				 are set elsewhere */
	
	for(ind = 1; ind <= dna->numseq[population]; ind++)
		for(categ = 0; categ < numcategs; categ++) {
			tr->nodep[ind]->x->s[site][categ][slice][baseA] = 1.0; 
			tr->nodep[ind]->x->s[site][categ][slice][baseC] = 0.0; 
			tr->nodep[ind]->x->s[site][categ][slice][baseG] = 0.0; 
			tr->nodep[ind]->x->s[site][categ][slice][baseT] = 0.0; 
			
			tr->nodep[ind]->x->s[site+1][categ][slice][baseA] = 0.0; 
			tr->nodep[ind]->x->s[site+1][categ][slice][baseC] = 1.0; 
			tr->nodep[ind]->x->s[site+1][categ][slice][baseG] = 0.0; 
			tr->nodep[ind]->x->s[site+1][categ][slice][baseT] = 0.0; 
			
			tr->nodep[ind]->x->s[site+2][categ][slice][baseA] = 0.0; 
			tr->nodep[ind]->x->s[site+2][categ][slice][baseC] = 0.0; 
			tr->nodep[ind]->x->s[site+2][categ][slice][baseG] = 1.0; 
			tr->nodep[ind]->x->s[site+2][categ][slice][baseT] = 0.0; 
			
			tr->nodep[ind]->x->s[site+3][categ][slice][baseA] = 0.0; 
			tr->nodep[ind]->x->s[site+3][categ][slice][baseC] = 0.0; 
			tr->nodep[ind]->x->s[site+3][categ][slice][baseG] = 0.0; 
			tr->nodep[ind]->x->s[site+3][categ][slice][baseT] = 1.0; 
		}
} /* makeinvarvalues */


/*****************************************************************
 * locusinit() initializes things that are specific to one locus */ 
void locusinit(option_struct *op, data_fmt *data)
{
	long i, nummarkers, numseq, numsites;
	dnadata *dna;
	FILE *flipfile;
	
	dna = data->dnaptr;
	
	if (op->spacing) read_spacefile(spacefile,op,data);
	nummarkers = getdata_nummarkers(op,data);
	/* hack done for speed-up */
	dna->sitecount[locus][0] = getdata_space(op,data,0L);
	dna->markersite[locus][0] = dna->sspace[0][0][getdata_nummarkers(op,data)];
	for(i = 1; i < nummarkers; i++) {
		dna->sitecount[locus][i] = dna->sitecount[locus][i-1] +
		getdata_space(op,data,i);
		dna->markersite[locus][i] = dna->markersite[locus][i-1] +
		dna->sspace[0][0][i-1];
	}
	numsites = countsites(op,data);
	dna->segranges0 = (int *)calloc(numsites,sizeof(int));
	dna->segranges1 = (int *)calloc(numsites,sizeof(int));
	dna->segranges2 = (int *)calloc(numsites,sizeof(int));
	dna->segranges3 = (int *)calloc(numsites,sizeof(int));
	for(i = 0; i < numsites; i++) {
		dna->segranges0[i] = 0;
		dna->segranges1[i] = 1;
		dna->segranges2[i] = 2;
		dna->segranges3[i] = -1;
	}
	
	treesetup(op,data);
	
	numseq = getdata_numseq(op,data);
	
	switch(op->datatype) {
		case 'a':
			break;
		case 'b':
		case 'm':
			break;
		case 'n':
			initinvartips(op,data,op->categs,curtree);
			if (!op->panel) invardatacheck(op,data);
		case 's':
			alogf = (rlrec *)calloc(1,sizeof(rlrec));
			alogf->val = (double *)calloc(nummarkers+1,sizeof(double));
			contribution = 
			(contribarr *)calloc(nummarkers+1,sizeof(contribarr));
			contribution[0] = 
			(double *)calloc((nummarkers+1)*op->categs,sizeof(double));
			for(i=1;i<nummarkers+1;i++)
				contribution[i] = contribution[0] + i*op->categs;
			/* now read the categories from the infile, if applicable */
			inputcategories(op,data);
			/* setup fractional likelihoods at tree tips */
			//makednavalues(op,data,op->categs,curtree);
			/* frequencies must be calculated after makednavalues() called */
			if (!op->freqsfrom) {
				data->dnaptr->freqa = freqa;
				data->dnaptr->freqc = freqc;
				data->dnaptr->freqg = freqg;
				data->dnaptr->freqt = freqt;
			} else empiricaldnafreqs_modified(op,data);// else empiricaldnafreqs(op,data,curtree);
			getbasednafreqs(data->dnaptr,op,locus_ttratio,outfile);
			/* table can only be initialized after frequencies are set */
			inittable(op,dna);
			makesiteptr(op,data);
			initweightrat(op,data);
			if (op->haplotyping) {
				flipfile = fopen("flipfile","r");
				read_flipfile(op,data,curtree,flipfile);
				fclose(flipfile);
				pruneflips(op,data,curtree);
			}
			break;
		default:
			fprintf(ERRFILE,"ERROR:locusinit, can't get here!\n");
			exit(-1);
	}
	
	
	if ((op->increm[0] < 0) || (op->increm[1] < 0)) {
		fprintf(ERRFILE,"Error in input sampling increment,");
		fprintf(ERRFILE," increment set to 10\n");
		if (op->increm[0] < 0)
			op->increm[0] = 10;
		if (op->increm[1] < 0)
			op->increm[1] = 10;
	}
	/* stuff for recycling x arrays */
	numx = 0;
	for (i=0;i<XARRAYSIZE;i++)
		sparex[i]=NULL;
	
} /* locusinit */


void inputcategories(option_struct *op, data_fmt *data)
{
	/*  char ch; */
	long i, extranum, nummarkers;
	
	nummarkers = getdata_nummarkers(op,data);
	
	category = (long *)calloc(nummarkers,sizeof(long));
	
	for (i = 0; i < nummarkers; i++) category[i] = 1;
	extranum = 0;
	/* DEBUG debug WARNING warning--categories input code commented out!
     while (!(eoln(infile))) {
     ch = getc(infile);
     if (ch == '\n')
     ch = ' ';
     ch = isupper(ch) ? ch : toupper((int)ch);
     if (ch == 'C')
     extranum++;
     else if (ch != ' ') {
     printf("BAD OPTION CHARACTER: %c\n", ch);
     exit(-1);
     }
     }
     fscanf(infile, "%*[^\n]");
     getc(infile);
     for (i = 1; i <= extranum; i++) {
     ch = getc(infile);
     if (ch == '\n')
     ch = ' ';
     ch = isupper(ch) ? ch : toupper((int)ch);
     if (ch != 'W'){
     printf("ERROR: INCORRECT AUXILIARY OPTIONS LINE WHICH STARTS WITH %c\n",
     ch);
     exit(-1);
     }
     }
	 */
	if (op->categs <= 1)
		return;
}  /* inputcategories */


/*******************************************************************
 * allocate_nodelet() allocates and initializes a set of nodelets, *
 * which are set to be non-tops of the passed "type".  "x" and all *
 * other datatype specific fields are NULLed and must be allocated *
 * elsewhere.                                                      */
node *allocate_nodelet(long num, char type)
{
	boolean isfirst=TRUE;
	long j;
	node *p, *q = NULL, *pfirst=NULL;
	
	for (j = 0; j < num; j++) {
		p = (node *)calloc(1,sizeof(node));
		
		p->top = FALSE;
		p->updated = FALSE;
		p->number = FLAGLONG;
		p->type = type;
		//p->id = unique_id++;
		p->next = q;
		p->back = NULL;
		p->nayme = NULL;
		p->x = NULL;
		p->v = p->tyme = p->length = 0.0;
		p->z = NULL;
		
		p->recstart = p->recend = -1L;
		p->ranges = NULL;
		p->coal = NULL;
		p->members = NULL;
		p->memberstrue = FALSE;
		p->futileflag = 0L;
		
		p->pop = p->actualpop = -1;
		p->lxmax = 0.0;
		
		if(isfirst){
			isfirst=FALSE;
			pfirst = p;
		}
		
		q = p;
	}
	
	pfirst->next = q;
	
	return q;
	
} /* allocate_nodelet */


/*************************************************************
 * allocate_root() allocates and initializes a root for "tr" */
void allocate_root(option_struct *op, data_fmt *data, tree *tr)
{
	long i, j, numslice, nummarkers;
	dnadata *dna;
	msatdata *ms;
	
	tr->nodep[ROOTNUM] = allocate_nodelet(1,'t');
	tr->nodep[ROOTNUM]->number = ROOTNUM;
	//allocate_x(op,data,tr->nodep[ROOTNUM]);
	if (op->map) {
		tr->nodep[ROOTNUM]->z = NULL;
		allocate_z(op,data,tr->nodep[ROOTNUM]);
	}
	ranges_Malloc(tr->nodep[ROOTNUM],FALSE,0L);
	coal_Malloc(tr->nodep[ROOTNUM],FALSE,0L);
	tr->nodep[ROOTNUM]->nayme = (char *)calloc(NMLNGTH,sizeof(char));
	strncpy(tr->nodep[ROOTNUM]->nayme,"ROOT",4);
	
	/* guarantee that the root node contributes nothing to the likelihood
     of a tree (since its supposed to be at the end of a theoretically
     infinite root branch) */
	switch(op->datatype) {
		case 'a':
			break;
		case 'b':
		case 'm':
			ms = data->msptr;
			for(i = 0; i < ms->numloci; i++)
				for(j = 0; j < MICRO_ALLELEMAX; j++)
					tr->nodep[ROOTNUM]->x->a[i][j] = 1.0;
			break;
		case 'n':
		case 's':
			dna = data->dnaptr;
			numslice = (op->panel) ? NUMSLICE : 1L;
			nummarkers = getdata_nummarkers(op,data);
			/*     for (i = 0; i < nummarkers; i++) { */
			/*       for (j = 0; j < op->categs; j++) { */
			/* 	for (k = 0; k < numslice; k++) { */
			/* 	  tr->nodep[ROOTNUM]->x->s[i][j][k][baseA] = 1.0; */
			/* 	  tr->nodep[ROOTNUM]->x->s[i][j][k][baseC] = 1.0; */
			/* 	  tr->nodep[ROOTNUM]->x->s[i][j][k][baseG] = 1.0; */
			/* 	  tr->nodep[ROOTNUM]->x->s[i][j][k][baseT] = 1.0; */
			/* 	} */
			/*       } */
			/*     } */
			break;
		default:
			fprintf(ERRFILE,"ERROR:allocate_root: can't get here!\n");
			exit(-1);
	}
	
} /* allocate_root */


/***********************************************************
 * allocate_tip() allocates and initializes a tip nodelet. */
void allocate_tip(option_struct *op, data_fmt *data, tree *tr,
				  long num)
{
	
	tr->nodep[num] = allocate_nodelet(1,'t');
	tr->nodep[num]->number = num;
	tr->nodep[num]->top = TRUE;
	tr->nodep[num]->tyme = 0.0;
	tr->nodep[num]->nayme = (char *)calloc((NMLNGTH+1),sizeof(char));
	//allocate_x(op,data,tr->nodep[num]);
	if (op->map) {
		tr->nodep[num]->z = NULL;
		allocate_z(op,data,tr->nodep[num]);
	}
	ranges_Malloc(tr->nodep[num],TRUE,1L);
	if (op->fc) coal_Malloc(tr->nodep[num],TRUE,1L);
	
	switch(op->datatype) {
		case 'a':
			break;
		case 'b':
		case 'm':
			break;
		case 'n':
		case 's':
			tr->nodep[num]->ranges[1] = 0L;
			tr->nodep[num]->ranges[2] = countsites(op,data)-1;
			if (op->fc) {
				tr->nodep[num]->coal[1] = 0L;
				tr->nodep[num]->coal[2] = countsites(op,data)-1;
			}
			break;
		default:
			fprintf(ERRFILE,"ERROR:allocate_tip: can't get here!\n");
			exit(-1);
	}
	
} /* allocate_tip */


/*******************************************************************
 * allocate_interior() allocates and initializes an interior node. */
void allocate_interior(option_struct *op, data_fmt *data, tree *tr, long num)
{
	
	if (freenodes == NULL) tr->nodep[num] = allocate_nodelet(3,'c');
	else {tr->nodep[num] = freenodes; freenodes = freenodes->back;}
	tr->nodep[num]->number = tr->nodep[num]->next->number =
    tr->nodep[num]->next->next->number = num;
	
} /* allocate_interior */


/***********************************************************************
 * init_creature() initializes an already allocated creature structure */
void init_creature(creature *cr, long numhaplotypes)
{
	
	cr->numhaplotypes = numhaplotypes;
	cr->numflipsites = 0;
	cr->haplotypes = (node **)calloc(numhaplotypes,sizeof(node *));
	cr->flipsites = NULL;
	
} /* init_creature */


/***************************************************************
 * treesetup() allocates and initializes the global "curtree". *
 * Curtree's tymelist will be handled later.                   */
void treesetup(option_struct *op, data_fmt *data)
{
	long i, j, whichtip, numtips, numnodes, numcreatures;
	node *p;
	
	numtips = getdata_numtips(op,data);
	numnodes = 2 * numtips;
	numcreatures = numtips/NUMHAPLOTYPES;
	curtree = (tree *)calloc(1,sizeof(tree));
	curtree->nodep = (node **)calloc(numnodes,sizeof(node *));
	
	allocate_root(op,data,curtree);
	
	for(i = 1; i <= numtips; i++) allocate_tip(op,data,curtree,i);
	for(i = numtips+1; i < numnodes; i++) allocate_interior(op,data,curtree,i);
	for(p = freenodes; p != NULL; p = p->back, i++)
		p->number = p->next->number = p->next->next->number = i;
	nodectr = i;
	curtree->likelihood = NEGMAX;
	curtree->coalprob = NEGMAX;
	curtree->numcoals = numtips-1;
	curtree->numrecombs = 0;
	
	if (op->haplotyping) {
		curtree->creatures = 
		(creature *)calloc(numcreatures,sizeof(creature));
		for(i = 0, whichtip = 1; i < numcreatures; i++) {
			init_creature(&curtree->creatures[i],NUMHAPLOTYPES);
			for(j = 0; j < curtree->creatures[i].numhaplotypes; j++)
				curtree->creatures[i].haplotypes[j] = curtree->nodep[whichtip++];
		}
		if (whichtip-1 > numtips)
			fprintf(ERRFILE,"\ntreesetup--whichtip too big!\n");
	} else curtree->creatures = NULL;
	
	newtymenode(&curtree->tymelist);
	curtree->tymelist->branchlist = (node **)calloc(numtips,sizeof(node *));
	curtree->tymelist->eventnode = curtree->nodep[1];
	
} /* treesetup */


/***************************************************************
 * end_of_population_free frees the working tree, curtree; and *
 * the freenode list plus the sparex array.                    */
void end_of_population_free(option_struct *op, data_fmt *data)
{
	dnadata *dna;
	
	dna = data->dnaptr;
	
	free(op->ne_ratio);
	
	switch(op->datatype) {
		case 'a':
			break;
		case 'b':
		case 'm':
			break;
		case 'n':
		case 's':
			break;
		default:
			fprintf(ERRFILE,"ERROR:end_of_population_free, can't be here!\n");
			exit(-1);
	}
	
} /* end_of_population_free */


/***************************************************
 * end_of_locus_free() frees locus specific stuff. */
void end_of_locus_free(option_struct *op, data_fmt *data)
{
	long i;
	node *p, *q;
	
	free(data->siteptr);
	
	switch(op->datatype) {
		case 'a':
			break;
		case 'b':
		case 'm':
			break;
		case 'n':
		case 's':
			free(category);  
			free(alogf->val);
			free(alogf);
			/* free the working arrays */
			free(weightrat);
			free(tbl);
			free(contribution[0]);
			free(contribution);
#if 0 /* pre-heating */
			freetree(op,data,curtree);
#endif
			for(i = 0; i < op->numtempchains; i++)
				freetree(op,data,temptrees[i]);
			for(p = freenodes; p != NULL; p = q) {
				q = p->back;
				freenodelet(op,data,p->next->next);
				freenodelet(op,data,p->next);
				freenodelet(op,data,p);
			}
			freenodes = NULL;
			/* the sparex arrays must be done after the tree and freenodes
			 have been freed (or anything else that uses free_x to
			 free it's xarrays). */
			for(i = 0; i < numx; i++) {
				free(sparex[i]->s[0][0][0]);
				free(sparex[i]->s[0][0]);
				free(sparex[i]->s[0]);
				free(sparex[i]->s);
				free(sparex[i]);
				sparex[i] = NULL;
			}
			numx = 0;
			break;
		default:
			fprintf(ERRFILE,"ERROR:end_of_locus_free, can't be here!\n");
			exit(-1);
	}
	
} /* end_of_locus_free */


boolean sitecompare(dnadata *dna, long site1, long site2)
{
	long i;
	
	for(i = 0; i < dna->numseq[population]; i++)
		if(dna->seqs[population][locus][i][site1] != 
		   dna->seqs[population][locus][i][site2])
			return FALSE;
	
	return TRUE;
} /* sitecompare */


/*****************************************************************
 * makesiteptr creates the siteptr array:                        *
 * if (!siteptr[site]) there is no alias, otherwise potentially  *
 * use siteptr[site]-1 for the alias site;                       *
 * if (siteptr[site]-1 < 0) then the alias is currently unusable *
 * due to recombination.                                         */
void makesiteptr(option_struct *op, data_fmt *data)
{
	long whichsite, lastsite, i, *ptrarray;
	boolean found;
	dnadata *dna;
	
	dna = data->dnaptr;
	
	if (op->datatype == 'n') lastsite = dna->sites[locus]+NUMINVAR;
	else lastsite = dna->sites[locus];
	
	ptrarray = (long *)calloc(lastsite,sizeof(long));
	ptrarray[0] = 0;
	
#if ALIASING
	for(whichsite = 1; whichsite < lastsite; whichsite++) {
		if (op->datatype == 'n')
			if (whichsite >= dna->sites[locus]) {
				ptrarray[whichsite] = 0;
				continue;
			}
		found = FALSE;
		for(i = whichsite - 1; i >= 0; i--) {
			if(sitecompare(dna,i,whichsite)) {
				ptrarray[whichsite] = i+1;
				found = TRUE;
				break;
			}
		}
		if (found) continue;
		ptrarray[whichsite] = 0;
	}
#endif
	
#if !ALIASING
	for(i = 0; i < dna->sites[locus]; i++) ptrarray[i] = 0;
#endif
	
	data->siteptr = ptrarray;
	
} /* makesiteptr */


void read_line(FILE *source, char *line)
{
	char *p;
	char *templine = 0;
	fgets(line,LINESIZE,source);
	while(line[0] == '\n' || line[0] == '#') fgets(line,LINESIZE,source);
	if((p = (char *)strchr(line,'\n')) != NULL) *p = '\0';
	
	if (*line == '\0')
		while(isspace(*line)){
			*line++;
		}
	*templine = *line;
	if(*line == '\0') read_line(source,line);
	
} /* read_line */


void read_header(option_struct *op, long *numpop, long *numloci,
				 long **numsites, char *title, char *dlm)
{
	long i;
	char *input, *tok, temp;
	
	input = (char *)calloc(LINESIZE,sizeof(char));
	
	read_line(infile,input);
	
	switch(lowercase(input[0])) {
		case 'a':
			break;
		case 'b':
		case 'm':
			sscanf(input,"%c%ld%ld%c%[^\n]",&temp,numpop,numloci,dlm,title);
			if (op->datatypeset)
				if (op->datatype != temp) {
					fprintf(ERRFILE,"\nYou chose a datatype of %c,",op->datatype);
					fprintf(ERRFILE,"but your infile shows a datatype of %c\n",temp);
					exit(-1);
				}
			op->datatype = temp;
			numsites = NULL;
			break;
		case 'n':
			if (op->autocorr) {
				printf("\nSNP data type is incompatible with autocorrelation");
				printf("\n----autocorrelation disabled-----\n");
				printf("press <enter> or <return> to continue\n");
				fgets(input,LINESIZE,stdin);
				fprintf(outfile,"\nSNP data type is incompatible with");
				fprintf(outfile," autocorrelation");
				fprintf(outfile,"\n----autocorrelation disabled-----\n");
				op->autocorr = FALSE;
			}
		case 's':
			sscanf(input,"%c%ld%ld%[^\n]",&temp,numpop,numloci,title);
			if (op->datatypeset)
				if (op->datatype != temp) {
					fprintf(ERRFILE,"\nYou chose a datatype of \"%c\",",op->datatype);
					fprintf(ERRFILE," but your infile shows a datatype");
					fprintf(ERRFILE," of \"%c\".\n",temp);
					exit(-1);
				}
			op->datatype = temp;
			read_line(infile,input);
			(*numsites) = (long *)calloc(*numloci,sizeof(long));
			for(i = 0; i < *numloci; i++) {
				while(isspace((int)*input)) (*input)++;
				if(i == 0) tok = strtok(input," ");
				else tok = strtok(NULL," ");
				(*numsites)[i] = atol(tok);
			}
			dlm = NULL;
			break;
		default:
			fprintf(ERRFILE,"WARNING--datatype not specified in infile\n");
			switch(op->datatype) {
				case 'a':
					break;
				case 'b':
					break;
				case 'm':
					sscanf(input,"%ld%ld%c%[^\n]",numpop,numloci,dlm,title);
					break;
				case 'n':
				case 's':
					sscanf(input,"%ld%ld%[^\n]",numpop,numloci,title);
					read_line(infile,input);
					(*numsites) = (long *)calloc(*numloci,sizeof(long));
					for(i = 0; i < *numloci; i++) {
						while(isspace((int)*input)) (*input)++;
						if(locus == 0) tok = strtok(input," ");
						else tok = strtok(NULL," ");
						(*numsites)[i] = atol(tok);
					}
					break;
				default:
					fprintf(ERRFILE,"ERROR:unknown datatype %c\n",op->datatype);
					fprintf(ERRFILE,"only the types a, m, s");
					fprintf(ERRFILE,"(electrophoretic alleles,\n");
					fprintf(ERRFILE,"microsatellite data, sequence data)");
					fprintf(ERRFILE,"are allowed.\n");
					exit(-1L);
					break;
			}
			break;
	}
	
	free(input);
	
} /* read_header */


void read_popheader(long *numind, char *poptitle)
{
	char *input;
	
	input = (char *)calloc(LINESIZE,sizeof(char));
	
	read_line(infile,input);
	sscanf(input,"%ld%[^\n]",numind,poptitle);
	
	free(input);
	
} /* read_popheader */


void getinput(option_struct *op, data_fmt *data)
{
	long i, numpop, numloci, numind, *numsites;
	char *title, *poptitle, dlm;
	
	title = (char *)calloc(LINESIZE,sizeof(char));
	poptitle = (char *)calloc(LINESIZE,sizeof(char));
	
	/* setup data structures and read in the first population's
     worth of data */
	read_header(op,&numpop,&numloci,&numsites,title,&dlm);
	read_popheader(&numind,poptitle);
	setupdata(data,op->datatype,numpop,numloci,numind,numsites,title);
	setpopstuff(data,0L,numind,poptitle,op->datatype);
	read_popdata(infile,data,0L,op);
	/* "true haplotype" code */
	if (TRUEHAP) {
		readtruehaps(op,&truehaps,numpop,numloci,numind,numsites);
		copyhaps(op,data,&starthaps,numpop,numloci,numind,numsites);
	}
	
	/* now read in all other populations' data */
	for(i = 1; i < numpop; i++) {
		read_popheader(&numind,poptitle);
		setpopstuff(data,i,numind,poptitle,op->datatype);
		read_popdata(infile,data,i,op);
	}
	
	
	if (op->datatype == 's' || op->datatype == 'n') {
		if (op->weights)
			inputdnaweights(data->dnaptr->sites[locus],data->dnaptr,op);
		free(numsites);
	}
	
	free(poptitle);
	free(title);
	
}  /* getinput */

void getsimoptions(option_struct *op, data_fmt *data){ // read options and reference haplotype data
	FILE *simoption, *recfile, *reffile;
	double rectemp;
	int numind,i;
	long numsites;
	long *seqsize;
	seqsize = (long *)calloc(1,sizeof(long));
	simoption = fopen("simoption","r");
	// number of reference sequence
	fscanf(simoption,"%d",&numind);
	
	// sequence length
	fscanf(simoption,"%ld",&numsites);
	seqsize[0] = numsites;
	
	// number of new haplotype
	fscanf(simoption,"%d",&new_hap);
	newhaps = (short int *)calloc(new_hap,sizeof(short int));
	
	// theta value
	fscanf(simoption,"%lf",&theta0);
	
	// growth rate
	fscanf(simoption,"%lf",&growth);
	
	// recrate - 4Ner ? 
	fscanf(simoption,"%lf",&rectemp);
	
	
	recarray = (double *)calloc(numsites + 1, sizeof(double));
	
	
	if (rectemp == -1){
		recfile = fopen("recrates","r");
		for (i = 0; i < numsites; i++){
			fscanf(recfile,"%lf",&rectemp);
			recarray[i] = rectemp/theta0;
		}
		fclose(recfile);
	}
	else{
		if (rectemp == 0){
			rectemp = 0.0000000000001;
		}
		for (i = 0; i < numsites; i++){
			recarray[i] = rectemp/theta0;
		}
	}
	
	// number of minimum MCMC chain
	fscanf(simoption,"%ld",&minMCMC);
	
	// number of replicates
	fscanf(simoption,"%d",&numrep);
	
	// number of intervals
	fscanf(simoption,"%d",&repinter);
	
	// filtered ?
	fscanf(simoption,"%d",&filtered);
	// cutoff ?
	if (filtered == 1){
		fscanf(simoption,"%lf",&cutofffreq);
		
		calling = (double *)calloc(numind,sizeof(double)); // probability of calling SNPs given i variant
		calling[0] = 0;
		for (i = 1; i <= cutofffreq; i++){
			calling[i] = 0;
		}
		for (i = cutofffreq + 1; i < numind; i++){
			calling[i] = 1;
		}
	}
	
	if (filtered == -1){
		recfile = fopen("probs","r");
		calling = (double *)calloc(numind,sizeof(double));
		for (i = 0; i < numind; i++){
			fscanf(recfile,"%lf",&rectemp);
			calling[i] = rectemp;
		}
		calling[0] = 0;
		fclose(recfile);
	}
	
	numMCMC = minMCMC + numrep * repinter + 1;
	
	// read reference haplotypes
	reffile = fopen("refhap","r");
	
	/* setup data structures and read in the first population's
     worth of data */
	setupdata(data,'s',1,1,numind,seqsize,"");
	setpopstuff(data,0L,numind,"",op->datatype);
	read_popdata(reffile,data,0L,op);	
	fclose(reffile);
	
	
	
	return;
	
}





void getnodata(option_struct *op, data_fmt *data)
{
	long i, numpop, numloci, numchains, numind, *numsites;
	char *title, *poptitle, dlm;
	FILE *kks;
	
	/* debug DEBUG warning WARNING--bogus default settings */
	op->datatype = 's';
	numpop = 1;
	numind = 1;
	title = "no data run";
	dlm = ' ';
	poptitle = "first pop";
	
	kks = fopen("kks","r");
	fscanf(kks,"%ld %ld",&numloci,&numchains);
	numsites = (long *)calloc(numloci,sizeof(long));
	
	/* debug DEBUG warning WARNING--bogus default settings */
	for(i = 0; i < numloci; i++) numsites[i] = 1;
	
	setupdata(data,op->datatype,numpop,numloci,numind,numsites,title);
	setpopstuff(data,0L,numind,poptitle,op->datatype);
	
	fseek(kks,0,SEEK_SET);
	fclose(kks);
	
} /* getnodata */


double watterson(option_struct *op, data_fmt *data)
{
	/* estimate theta using method of Watterson */
	long i, j, k, kn, numprivate, numsites;
	boolean varies, private;
	double watter;
	dnadata *dna;
	
	dna = data->dnaptr;
	numprivate = kn = 0;
	
	for (i = 0; i < dna->sites[locus]; i++) {
		varies = FALSE;
		private = TRUE;
		for (j = 1; j < dna->numseq[population]; j++) {
			if (dna->seqs[population][locus][j][i] != 
				dna->seqs[population][locus][0][i]) {
				if (varies) {
					for(k = 1; k < dna->numseq[population]; k++) {
						if (dna->seqs[population][locus][k][i] !=
							dna->seqs[population][locus][j][i])
							private = FALSE;
					}
				}
				varies = TRUE;
			}
		}
		if (varies) kn++;
		if (private && varies) numprivate++;
	}
	watter = 0.0;
	if (kn > 0) {
		for (i = 1; i < dna->numseq[population]; i++)
			watter += 1.0 / i;
		numsites = countsites(op,data);
		watter = kn / (numsites * watter);
		fprintf(outfile,"\nThere are %ld variable sites in the data\n",kn);
		
		if (kn == 1) {
			op->holding = 2;
			fprintf(outfile,"\nWARNING:  There is only 1 variable site in");
			fprintf(outfile," this data set.\nRecombination-rate cannot");
			fprintf(outfile," be accurately estimated in this case.\n");
			fprintf(outfile,"The \"estimate\" of rec-rate has been fixed");
			fprintf(outfile," to its starting value.\n\n");
			
			printf("\nWARNING:  There is only 1 variable site in");
			printf(" this data set.\nRecombination-rate cannot");
			printf(" be accurately estimated in this case.\n");
			printf("The \"estimate\" of rec-rate has been fixed");
			printf(" to its starting value.\n\n");
		} else {
			if (kn-numprivate < 2) {
				op->holding = 2;
				fprintf(outfile,"\nWARNING:  Too many of the variable sites ");
				fprintf(outfile,"are uninformative.\nRecombination-rate ");
				fprintf(outfile,"cannot be accurately estimated in this ");
				fprintf(outfile,"case.\nThe \"estimate\" of rec-rate has ");
				fprintf(outfile,"been fixed to its starting value.\n\n");
				
				printf("\nWARNING:  Too many of the variable sites ");
				printf("are uninformative.\nRecombination-rate ");
				printf("cannot be accurately estimated in this ");
				printf("case.\nThe \"estimate\" of rec-rate has ");
				printf("been fixed to its starting value.\n\n");
			}
		}
		
		return watter;
	}
	fprintf(outfile, "Warning:  There are no variable sites");
	fprintf(outfile, " in this data set.\n\n");
	if (MENU) printf("Warning:  There are no variable sites in this data set.\n");
	else {
		fprintf(simlog, "Warning:  There are no variable sites");
		fprintf(simlog, " in this data set.\n\n");
	}
	exit(-1);
}  /* watterson */


/**********************************************************************
 * invardatacheck() is called to check non-panel SNP data for         *
 * invariant sites, which are illegal under that model.  If it finds  *
 * any, it changes that site to a datatype of "unknown"; "x" is used. */
void invardatacheck(option_struct *op, data_fmt *data)
{
	long marker, seq, nummarkers, numseq;
	char **dna;
	boolean varies;
	
	nummarkers = getdata_nummarkers(op,data);
	numseq = getdata_numseq(op,data);
	
	dna = data->dnaptr->seqs[population][locus];
	
	for(marker = 0; marker < nummarkers; marker++) {
		varies = FALSE;
		for(seq = 1; seq < numseq; seq++) {
			if(dna[0][marker] != dna[seq][marker]) {
				varies = TRUE;
				break;
			}
		}
		if (!varies)
			for(seq = 0; seq < numseq; seq++) dna[seq][marker] = 'X';
	}
	
} /* invardatacheck */


/*********************************************************************
 * buildtymelist() allocates space and initializes all the fields of *
 * a tymelist except for failure to initialize the branchlists.      *
 * Those are handled by finishsetup().                               *
 *                                                                   *
 * WARNING buildtymelist() currently assumes a non-recombinant tree! */
void buildtymelist(tree *tr, option_struct *op, data_fmt *data, node *p)
{
	long numentries;
	tlist *t, *u;
	
	t = tr->tymelist;
	
	newtymenode(&u);
	
	numentries = getdata_numtips(op,data);
	u->branchlist = (node **)calloc(numentries,sizeof(node *));
	
	u->eventnode = p;
	
	while (t != NULL) {
		if (u->eventnode->tyme < t->eventnode->tyme) {
			u->prev = t->prev;
			t->prev = u;
			u->succ = t;
			u->prev->succ = u;
			break;
		}
		if (t->succ != NULL)
			t = t->succ;
		else {
			t->succ = u;
			u->prev = t;
			u->succ = NULL;
			break;
		}
	}
	
} /* buildtymelist */


/* WARNING warning orient assumes a non-recombinant tree */
void orient(tree *tr, option_struct *op, data_fmt *data, node *p)
{
	long lastsite;
	
	lastsite = countsites(op,data)-1;
	
	if (istip(p)) {
		return;
	}
	
	curtree->nodep[p->number] = p;  /* insure that curtree->nodep points
									 to nodes with info */
	
	/* since p is a top nodelet, it needs to actually store
     likelihood information, x is a NULL pointer
     in all other non-tip nodelets */
	p->top = TRUE;
	//allocate_x(op,data,p);
	if (op->map) {
		p->z = NULL;
		allocate_z(op,data,p);
	}
	ranges_Malloc(p,TRUE,1L);
	if (op->fc) { 
		coal_Malloc(p,TRUE,1L); 
		p->ranges[1] = p->coal[1] = 0;
		p->ranges[2] = p->coal[2] = lastsite;
	} else {
		p->ranges[1] = 0;
		p->ranges[2] = lastsite;
	}
	p->next->top = FALSE;
	free_x(op,p->next);
	ranges_Malloc(p->next,FALSE,0L);
	coal_Malloc(p->next,FALSE,0L);
	p->next->next->top = FALSE;
	free_x(op,p->next->next);
	ranges_Malloc(p->next->next,FALSE,0L);
	coal_Malloc(p->next->next,FALSE,0L);
	
	orient(tr,op,data,p->next->back);
	orient(tr,op,data,p->next->next->back);
	
	p->tyme = p->next->length + p->next->back->tyme;
	p->next->tyme = p->tyme;
	p->next->next->tyme = p->tyme;
	if (p->number == curtree->root->back->number) {
		p->back->top = FALSE;
		p->back->tyme = ROOTLENGTH;
	}
	
	buildtymelist(tr,op,data,p);
	
}  /* orient */

void plumptree(option_struct *op, data_fmt *data, double th0)
/* change an input tree into a perfect coalescent tree for theta0 */
/* WARNING--doesn't work on a tree with recombinations!!!!! */
{
	tlist *t;
	long k;
	double tyme;
	
	/* we assume that tips are tyme 0 already */
	t = curtree->tymelist->succ;
	
	k = getdata_numtips(op,data);
	tyme = 0.0;
	do {
		tyme += th0/(k*(k-1));
		if (!istip(t->eventnode)) {
			t->eventnode->tyme = tyme;
			t->eventnode->next->tyme = tyme;
			t->eventnode->next->next->tyme = tyme;
		}
		k--;
		t = t->succ;
	} while (t != NULL);
	
	/* now you must call finishsetup and initbranchlist or a disaster
     will happen! */
} /* plumptree */

void finishsetup(option_struct *op, data_fmt *data, node *p)
{
	if (istip(p)) {
		ltov(op,data,p);
		return;
	}
	ltov(op,data,p);
	finishsetup(op,data,p->next->back);
	finishsetup(op,data,p->next->next->back);
	return;
} /* finishsetup */

void initbranchlist(option_struct *op, data_fmt *data)
{
	tlist *t;
	node *p, *q;
	long i, j, k, n, numtips;
	
	t = curtree->tymelist;
	numtips = n = getdata_numtips(op,data);
	t->numbranch = numtips;
	for(i = 0; i < numtips; i++) t->branchlist[i] = curtree->nodep[i+1];
	t->age = t->succ->eventnode->tyme;
	t = t->succ;
	/* WARNING assumes initial tree is not recombinant! */
	for (i = 0; i < numtips-1; i++) {
		/* for each interior node, do... */
		n--;
		t->numbranch = n;
		if (n == 1)
			t->age = t->eventnode->tyme + ROOTLENGTH;
		else
			t->age = t->succ->eventnode->tyme;
		p = t->eventnode->next->back;
		q = t->eventnode->next->next->back;
		k = 0;
		for (j = 0; j < t->prev->numbranch ; j++) {
			/* for the number of branches above the coalescent node, do...*/
			if (t->prev->branchlist[j] != p && t->prev->branchlist[j] != q) {
				t->branchlist[k] = t->prev->branchlist[j];
				k++;
			}
		}
		t->branchlist[t->numbranch - 1] = t->eventnode;
		t = t->succ;
	}
}  /* initbranchlist */

void inittable(option_struct *op, dnadata *dna)
{
	long i;
	tbl = (valrec *)calloc(1,op->categs*sizeof(valrec));
	/* Define a lookup table. Precompute values and store them in a table */
	for (i = 0; i < op->categs; i++) {
		tbl[i].rat_xi = op->rate[i] * dna->xi;
		tbl[i].rat_xv = op->rate[i] * dna->xv;
	}
}  /* inittable */

void initweightrat(option_struct *op, data_fmt *data)
{
	/* WARNING non-DNA data is not yet able to do this! */
	long i, nummarkers;
	
	nummarkers = getdata_nummarkers(op,data);
	weightrat = (double *)calloc(nummarkers,sizeof(double));
	sumweightrat = 0.0;
	for (i = 0; i < nummarkers; i++) {
		weightrat[i] = data->dnaptr->dnaweight[i] * op->rate[category[i]-1];
		sumweightrat += weightrat[i];
	}
}  /* initweightrat */

void rec_outtree(node *p, boolean first, FILE **usefile)
/*  write out a recombinant tree in KYB format; call with
 curtree.root->back.  Calling program should append a
 semicolon.  "first" is used to avoid an extra comma
 in the tree list and should be TRUE when this is called. */
{
	long i;
	
	if(istip(p)) {
		if(!first)fprintf(*usefile,",");
		else first=FALSE;
		fprintf(*usefile,"\"");
		for(i=0;i<NMLNGTH;i++) {
			if(p->nayme[i]==' ') break;
			else putc(p->nayme[i],*usefile);
		}
		fprintf(*usefile,"\":%f",p->length);
		return;
	} 
	if(isrecomb(p)) {
		if(!first)fprintf(*usefile,",");
		else first=FALSE;
		fprintf(*usefile,"%ld",p->number);
		fprintf(*usefile,"{%ld-%ld}",p->recstart,p->recend);
		fprintf(*usefile,":%f",p->length);
		p->updated = !p->updated;
		if (p->updated) {
			if (p->next->top) rec_outtree(p->next->back, first, usefile);
			else rec_outtree(p->next->next->back, first, usefile);
		}
		return;
	} else { /* p is coalescent */
		if(!first)fprintf(*usefile,",");
		else first=FALSE;
		fprintf(*usefile,"*%ld",p->number);
		fprintf(*usefile,":%f",p->length);
		rec_outtree(p->next->back, first, usefile);
		rec_outtree(p->next->next->back, first, usefile);
		return;
	}
} /* rec_outtree */

void treeout(node *p, long s, FILE **usefile)
{
	/* write out file with representation of final tree */
	long i, n, w;
	char c;
	double x;
	
	if (istip(p)) {
		n = 0;
		for (i = 1; i <= NMLNGTH; i++) {
			if (p->nayme[i - 1] != ' ')
				n = i;
		}
		for (i = 0; i < n; i++) {
			c = p->nayme[i];
			if (c == ' ')
				c = '_';
			putc(c, *usefile);
		}
		col += n;
	} else {
		putc('(', *usefile);
		col++;
		treeout(p->next->back, s, usefile);
		putc(',', *usefile);
		col++;
		if (col > 45) {
			putc('\n', *usefile);
			col = 0;
		}
		treeout(p->next->next->back, s, usefile);
		putc(')', *usefile);
		col++;
	}
	if (p->v >= 1.0)
		x = -1.0;
	else
		x = lengthof(p);
	if (x > 0.0)
		w = (long)(0.4343 * log(x));
	else if (x == 0.0)
		w = 0;
	else
		w = (long)(0.4343 * log(-x)) + 1;
	if (w < 0)
		w = 0;
	if (p == curtree->root->back)
		putc(';', *usefile);
	else {
		fprintf(*usefile, ":%*.10f", (int)(w + 7), x);
		col += w + 8;
	}
}  /* treeout */


double snp_eval_calcrange(option_struct *op, data_fmt *data, tree *tr,
						  boolean first, long start, long finish, long numcategs)
{
	double sumterm, sumterm2, lterm, lterm2, temp;
	long i, j, k, slice, snpstart, snpfinish, nummarkers, numinvar, *siteptr;
	contribarr term, term2, clai;
	node *p;
	double *x1;
	dnadata *dna;
	
	dna = data->dnaptr;
	
	term = (double *)calloc(numcategs,sizeof(double));
	term2 = (double *)calloc(numcategs,sizeof(double));
	clai = (double *)calloc(numcategs,sizeof(double));
	
	slice = 0;
	temp = 0.0;
	p = tr->root->back;
	siteptr = data->siteptr;
	
	nummarkers = getdata_nummarkers(op, data);
	findsubtreemarkers(op,data,start,finish,&snpstart,&snpfinish);
	
	if (snpstart == FLAGLONG) numinvar = finish-start+1;
	else numinvar = (finish-start+1) - (snpfinish-snpstart+1);
	
	sumterm2 = 0.0;
	for (j = 0; j < numcategs; j++) { /* compute P(Invar) */
		term2[j] = 0.0;
		for(k = nummarkers; k < nummarkers+NUMINVAR; k++) {
			x1 = p->x->s[k][j][slice];
			term2[j] += dna->freqa * x1[baseA] + 
			dna->freqc * x1[baseC] + dna->freqg * x1[baseG] +
			dna->freqt * x1[baseT];
		}
		sumterm2 += op->probcat[j] * term2[j];
	}
	lterm2 = log(sumterm2);
	for (j = 0; j < numcategs; j++) {
		clai[j] = term2[j] / sumterm2;
	}
	memcpy(contribution[nummarkers], clai, numcategs*sizeof(double));
	if (!op->autocorr)
		alogf->val[nummarkers] = lterm2;
	if (!op->norecsnp) temp += numinvar * lterm2;
	
	if (snpstart != FLAGLONG) { /* that is, if there are SNPs */
		for (i = snpstart; i <= snpfinish; i++) {
			for (j = 0; j < numcategs; j++) {
				x1 = p->x->s[i][j][slice];
				term[j] = dna->freqa * x1[baseA] + dna->freqc * x1[baseC] +
				dna->freqg * x1[baseG] + dna->freqt * x1[baseT];
				if(term[j] == 0) {
					fprintf(ERRFILE,"Encountered tree incompatible with data\n");
					if(first) {
						fprintf(ERRFILE,"starting tree needs to be legal\n");
						exit(-1);
					}
					curtree->likelihood = NEGMAX;
					return(-1);
				}
			}
			sumterm = 0.0;
			for (j = 0; j < numcategs; j++) {
				sumterm += op->probcat[j] * term[j];
			}
			lterm = log(sumterm);
			if (op->norecsnp && !op->full) lterm -= log(1.0-sumterm2);
			for (j = 0; j < numcategs; j++) {
				clai[j] = ((op->norecsnp) ? 
						   (term[j] / (1.0-term2[j])) / (sumterm / (1.0-sumterm2))
						   : (term[j] / sumterm));
			}
			memcpy(contribution[i], clai, numcategs*sizeof(double));
			if (!op->autocorr)
				alogf->val[i] = lterm;
			temp += dna->dnaweight[i] * lterm;
		}
		
#if 0 /* This is probably dead DEAD dead--warning WARNING DEBUG debug */
		if (op->datatype == 'n' && op->full) {
			temp += log(op->chance_seen);
			numunobs = getnumlinks(op,data,start,finish)-(finish-start);
			if (numunobs < 0 || (long)numunobs != numunobs) {
				fprintf(ERRFILE,"non integral SNP distances with Full model\n");
				exit(-1);
			}
			temp += numunobs * log(1.0 - (op->chance_seen*(1.0-sumterm2)));
		}
#endif
	}
	
	free(term);
	free(term2);
	free(clai);
	
	return(temp);
	
} /* snp_eval_calcrange */


double panel_eval_calcrange(option_struct *op, data_fmt *data, tree *tr,
							boolean first, long start, long finish, long numcategs)
{
	double sumterm, sumterm2, lterm, lterm2, temp, invarterm;
	long i, j, k, numunobs, nummarkers, numinvar, snpstart, snpfinish, 
    *siteptr;
	contribarr term, term2, clai;
	node *p;
	double **x1;
	dnadata *dna;
	
	dna = data->dnaptr;
	
	term = (double *)calloc(numcategs,sizeof(double));
	term2 = (double *)calloc(numcategs,sizeof(double));
	clai = (double *)calloc(numcategs,sizeof(double));
	
	temp = 0.0;
	p = tr->root->back;
	siteptr = data->siteptr;
	
	/* compute P(INVAR) for each category, used as P(SNP) = 1.0 - P(INVAR) */
	nummarkers = getdata_nummarkers(op,data);
	findsubtreemarkers(op,data,start,finish,&snpstart,&snpfinish);
	
	if (snpstart == FLAGLONG) numinvar = finish-start+1;
	else numinvar = (finish-start+1) - (snpfinish-snpstart+1);
	
	sumterm2 = 0.0;
	for (j = 0; j < numcategs; j++) {
		x1 = p->x->s[nummarkers][j];
		term2[j] = 0.0;
		for (k = 1; k < NUMSLICE; k++) {
			term2[j] += dna->freqa * x1[k][baseA] +
			dna->freqc * x1[k][baseC] + dna->freqg * x1[k][baseG] +
			dna->freqt * x1[k][baseT];
		}
		sumterm2 += op->probcat[j] * term2[j];
	}
	lterm2 = log(sumterm2);
	for (j = 0; j < numcategs; j++) {
		clai[j] = term2[j] / sumterm2;
	}
	memcpy(contribution[nummarkers], clai, numcategs*sizeof(double));
	if (!op->autocorr)
		alogf->val[nummarkers] = lterm2;
	if (!op->norecsnp) temp += numinvar * lterm2;
	
	if (snpstart != FLAGLONG) {
		for (i = snpstart; i <= snpfinish; i++) {
			for (j = 0; j < numcategs; j++) {
				/* compute P(data) for this category */
				x1 = p->x->s[i][j];
				term[j] = dna->freqa * x1[0][baseA] + dna->freqc * x1[0][baseC] +
				dna->freqg * x1[0][baseG] + dna->freqt * x1[0][baseT];
				invarterm = 0.0;
				/* compute P(data && panel invariant) for this category */
				for (k = 1; k < NUMSLICE; k++) {
					invarterm += dna->freqa * x1[k][baseA] + dna->freqc *
					x1[k][baseC] + dna->freqg * x1[k][baseG] + dna->freqt *
					x1[k][baseT];
				}
				term[j] -= invarterm;
				if(term[j] == 0) {
					fprintf(ERRFILE,"Encountered tree incompatible with data\n");
					if(first) {
						fprintf(ERRFILE,"starting tree needs to be legal\n");
						exit(-1);
					}
					curtree->likelihood = NEGMAX;
					return(-1);
				}
			}
			sumterm = 0.0;
			for (j = 0; j < numcategs; j++) {
				sumterm += op->probcat[j] * term[j];
			}
			lterm = log(sumterm);
			if (op->norecsnp && !op->full) lterm -= log(1.0-sumterm2);
			for (j = 0; j < numcategs; j++) {
				clai[j] = ((op->norecsnp) ?
						   (term[j] / (1.0 - term2[j])) / (sumterm / (1.0 - sumterm2))
						   : (term[j] / sumterm));
			}
			memcpy(contribution[i], clai, numcategs*sizeof(double));
			if (!op->autocorr)
				alogf->val[i] = lterm;
			temp += dna->dnaweight[i] * lterm;
		}
		
		if (op->full) {
			temp += log(op->chance_seen);
			numunobs = (long)(getnumlinks(op,data,start,finish)-(finish-start));
			if (numunobs < 0 || (long)numunobs != numunobs) {
				fprintf(ERRFILE,"non integral SNP distances with Full model\n");
				exit(-1);
			}
			temp += numunobs * log(1.0 - (op->chance_seen*(1.0-sumterm2)));
		}
	}
	
	free(term);
	free(term2);
	free(clai);
	
	return(temp);
	
} /* panel_eval_calcrange */


double eval_calcrange(option_struct *op, data_fmt *data, tree *tr,
					  boolean first, long start, long finish, long numcategs)
{
	double sumterm, lterm, temp;
	long i, j, slice, nummarkers, *siteptr;
	contribarr term, term2, clai;
	node *p;
	double *x1;
	dnadata *dna;
	
	dna = data->dnaptr;
	
	term = (double *)calloc(numcategs,sizeof(double));
	term2 = (double *)calloc(numcategs,sizeof(double));
	clai = (double *)calloc(numcategs,sizeof(double));
	
	slice = 0;
	temp = 0.0;
	p = tr->root->back;
	siteptr = data->siteptr;
	nummarkers = getdata_nummarkers(op,data);
	
	for (i = start; i <= finish; i++) {
		for (j = 0; j < numcategs; j++) {
			x1 = p->x->s[i][j][slice];
			term[j] = dna->freqa * x1[baseA] + dna->freqc * x1[baseC] +
			dna->freqg * x1[baseG] + dna->freqt * x1[baseT];
			if(term[j] == 0) {
				fprintf(ERRFILE,"Encountered tree incompatible with data\n");
				if(first) {
					fprintf(ERRFILE,"starting tree needs to be legal\n");
					exit(-1);
				}
				curtree->likelihood = NEGMAX;
				return(-1);
			}
		}
		sumterm = 0.0;
		for (j = 0; j < numcategs; j++) {
			sumterm += op->probcat[j] * term[j];
		}
		lterm = log(sumterm);
		for (j = 0; j < numcategs; j++) {
			clai[j] = term[j] / sumterm;
		}
		memcpy(contribution[i], clai, numcategs*sizeof(double));
		if (!op->autocorr)
			alogf->val[i] = lterm;
		temp += dna->dnaweight[i] * lterm;
	}
	
	free(term);
	free(term2);
	free(clai);
	
	return(temp);
	
} /* eval_calcrange */


double evaluate(option_struct *op, data_fmt *data, tree *tr,
				double llike)
{
	double sum2, sumc, one_minus_lambda;
	contribarr like, nulike, clai;
	long i, j, k, nummarkers;
	
	if (op->categs > 1) {
		nummarkers = getdata_nummarkers(op,data);
		
		like = (double *)calloc(op->categs,sizeof(double));
		nulike = (double *)calloc(op->categs,sizeof(double));
		
		one_minus_lambda = 1.0 - op->lambda;
		
		for (j = 0; j < op->categs; j++) like[j] = 1.0;
		for (i = 0; i < nummarkers; i++) {
			sumc = 0.0;
			for (k = 1; k <= op->categs; k++)
				sumc += op->probcat[k - 1] * like[k - 1];
			sumc *= op->lambda;
			clai = contribution[i];
			for (j = 0; j < op->categs; j++)
				nulike[j] = (one_minus_lambda * like[j] + sumc) * clai[j];
			memcpy(like, nulike, op->categs*sizeof(double));
		}
		sum2 = 0.0;
		for (i = 0; i < op->categs; i++)
			sum2 += op->probcat[i] * like[i];
		llike += log(sum2);
		tr->likelihood = llike;
		
		free(like);
		free(nulike);
		
		return(llike);
	} else {
		tr->likelihood = llike;
		return(llike);
	}
	
}  /* evaluate */


double snp_evaluate(option_struct *op, data_fmt *data, tree *tr,
					double llike, long start, long finish)
{
	double sum2, sumc;
	contribarr like, nulike, clai;
	long i, j, snpstart, snpfinish, nummarkers, numinvar;
	
	if (op->categs > 1) {
		nummarkers = getdata_nummarkers(op,data);
		findsubtreemarkers(op,data,start,finish,&snpstart,&snpfinish);
		
		like = (double *)calloc(op->categs,sizeof(double));
		nulike = (double *)calloc(op->categs,sizeof(double));
		
		for (j = 0; j < op->categs; j++) like[j] = 1.0;
		
		if (snpstart == FLAGLONG) numinvar = finish-start+1;
		else numinvar = (finish-start+1) - (snpfinish-snpstart+1);
		
		if (snpstart != FLAGLONG) {
			for (i = snpstart; i <= snpfinish; i++) {
				sumc = 0.0;
				for (j = 0; j < op->categs; j++)
					sumc += op->probcat[j] * like[j];
				clai = contribution[i];
				for (j = 0; j < op->categs; j++)
					nulike[j] = sumc * clai[j];
				memcpy(like, nulike, op->categs*sizeof(double));
			}
		}
		
		for (i = 0; i < numinvar; i++) { /* invariant sites */
			sumc = 0.0;
			for (j = 0; j < op->categs; j++)
				sumc += op->probcat[j] * like[j];
			clai = contribution[nummarkers];
			for (j = 0; j < op->categs; j++)
				nulike[j] = sumc * clai[j];
			memcpy(like, nulike, op->categs*sizeof(double));
		}
		
		sum2 = 0.0;
		for (j = 0; j < op->categs; j++)
			sum2 += op->probcat[j] * like[j];
		llike += log(sum2);
		
		free(like);
		free(nulike);
		
		tr->likelihood = llike;
		return(llike);
	} else {
		tr->likelihood = llike;
		return(llike);
	}
	
}  /* snp_evaluate */


boolean nuview_usebranch(data_fmt *data, node *p, long site)
{
	
	if(!p->top) return(FALSE);
	
	return(inrange(p->ranges,site));
	
} /* nuview_usebranch */


double prob_micro(msatdata *ms, double t, long diff)
{
	double **steps = ms->steps;
	long stepnum = MICRO_MAXCHANGE;
	long k;
	double newsum = 0.0, oldsum=0.0;
	double logt = log(t);
	
	if(diff>=stepnum)
		return newsum;
	for (k = diff; k < diff + stepnum; k+=2){
		newsum += exp(-t + logt * k - steps[diff][k-diff]);
		if(oldsum-newsum < DBL_EPSILON)
			break;
		oldsum = newsum;
	}
	return newsum;
	
} /* prob_micro */


/********************************************************
 * findcoal() takes a non-top coalescent nodelet and    *
 * returns the closest tipwards top coalescent nodelet. *
 * It also returns the "v" value between the two nodes. *
 * In the case of failure it returns a NULL pointer     *
 * with v set to FLAGLONG.                              *
 *                                                      *
 * findcoal assumes that branchlengths are additive.    */
node *findcoal(option_struct *op, data_fmt *data, node *p, double *v)
{
	double total;
	node *q;
	
	if(!iscoal(p) || p->top) {*v = FLAGLONG; return(NULL);}
	
	for (q = p->back, total = 0; !iscoal(q) && !istip(q); q = findunique(q)->back)
		total += lengthof(q);
	
	total += lengthof(q);
	*v = findcoal_ltov(op,data,total);
	
	return(q);
	
} /* findcoal */


void nuview_micro(option_struct *op, data_fmt *data, node *p,
				  long indexsite, long startsite, long endsite)
{
	long i, a, s, smax, margin, diff;
	double *xx1, *xx2, *xx3, vv1 = 0.0, vv2 = 0.0, pija1s, pija2s,
    lx1, lx2, x3m;
	node *q, *r;
	boolean qactive, ractive;
	
	smax = MICRO_ALLELEMAX;
	margin = MICRO_MAXCHANGE;
	x3m = NEGMAX;
	
	xx1 = NULL;
	xx2 = NULL;
	vv1 = 0.0;
	vv2 = 0.0;
	xx3 = (double *)calloc(smax,sizeof(double));
	
	q = findcoal(op,data,p->next,&vv1);
	r = findcoal(op,data,p->next->next,&vv2);
	
	/* are these sites active and on an upwards branch?  */
	qactive = nuview_usebranch(data,p->next,indexsite);
	ractive = nuview_usebranch(data,p->next->next->back,indexsite);
	
	for(i = startsite; i <= endsite; i++) {
		if (qactive) {
			xx1 = q->x->a[i];
			lx1 = q->lxmax;
		} else {
			lx1 = 0.0;
		}
		
		if (ractive) {
			xx2 = r->x->a[i];
			lx2 = r->lxmax;
		} else {
			lx2 = 0.0;
		}
		
		for(s = 0; s < smax; s++) {
			if (qactive) pija1s = 0.0;
			else pija1s = 1.0;
			if (ractive) pija2s = 0.0;
			else pija2s = 1.0;
			for(a = MAX(0,s-margin); a < s + margin && a < smax; a++) {
				diff = labs(s-a);
				if(xx1[a] > 0 && qactive) {
					pija1s += prob_micro(data->msptr,vv1,diff) * xx1[a];
				}
				if(xx2[a] > 0 && ractive) {
					pija2s += prob_micro(data->msptr,vv2,diff) * xx2[a];
				}
			}
			xx3[s] = pija1s * pija2s;
			if(xx3[s] > x3m) x3m = xx3[s];
		}
		if(x3m == 0.0) p->lxmax = NEGMAX;
		else {
			for(s = 0; s < smax; s++) xx3[s] /= x3m;
			p->lxmax = log(x3m) + lx1 + lx2;
		}
		memcpy(p->x->a[i],xx3,smax*sizeof(double));
	}
	
	free(xx3);
	
} /* nuview_micro */


void nuview(option_struct *op, data_fmt *data, node *p, long indexsite,
			long startsite, long endsite)
{
	long i;
	double lw1 = 0.0, lw2 = 0.0;
	node *q, *r;
	boolean toobig;
	
	/* set "q" and "r" for use in calcrange(), either to a valid nodelet with
     datalikelihoods, or NULL if there is no such nodelet. */
	if (nuview_usebranch(data,p->next->back,indexsite)) {
		q = findcoal(op,data,p->next,&lw1);
		toobig = (exp(lw1) <= 0.0);
	} else {q = NULL; toobig = TRUE;}
	
	if (toobig) {
		for (i = 0; i < op->categs; i++) {
			tbl[i].ww1 = tbl[i].zz1 = 0.0;
			tbl[i].ww1zz1 = tbl[i].vv1zz1 = 0.0;
		}
	} else {
		for (i = 0; i < op->categs; i++) {
			tbl[i].ww1 = exp(tbl[i].rat_xi * lw1);
			tbl[i].zz1 = exp(tbl[i].rat_xv * lw1);
			tbl[i].ww1zz1 = tbl[i].ww1 * tbl[i].zz1;
			tbl[i].vv1zz1 = (1.0 - tbl[i].ww1) * tbl[i].zz1;
		}
	}
	
	if (nuview_usebranch(data,p->next->next->back,indexsite)) {
		r = findcoal(op,data,p->next->next,&lw2);
		toobig = (exp(lw2) <= 0.0);
	} else {r = NULL; toobig = TRUE;}
	
	if (toobig) {
		for (i = 0; i < op->categs; i++) {
			tbl[i].ww2 = tbl[i].zz2 = 0.0;
			tbl[i].ww2zz2 = tbl[i].vv2zz2 = 0.0;
		}
	} else {
		for (i = 0; i < op->categs; i++) {
			tbl[i].ww2 = exp(tbl[i].rat_xi * lw2);
			tbl[i].zz2 = exp(tbl[i].rat_xv * lw2);
			tbl[i].ww2zz2 = tbl[i].ww2 * tbl[i].zz2;
			tbl[i].vv2zz2 = (1.0 - tbl[i].ww2) * tbl[i].zz2;
		}
	}
	
	calcrange(op,data,p,q,r,indexsite,startsite,endsite,op->categs);
	
}  /* nuview */


void calcrange(option_struct *op, data_fmt *data, node *p, node * q,
			   node *r, long whichtree, long start, long finish, long numcategs)
/* whichtree indicates which tree should be used to calculate the
 likelihood of this set of sites:  pass the number of a site which
 has the desired tree.  (This is helpful in calculations for the
 SNP invariant sites.) */
{
	long i, j, k, numslice, *siteptr;
	double yy1, yy2, ww1zz1, vv1zz1, ww2zz2, vv2zz2, vv1zz1_sumr1,
    vv2zz2_sumr2, vv1zz1_sumy1, vv2zz2_sumy2, sum1, sum2,
    sumr1, sumr2, sumy1, sumy2;
	boolean qactive, ractive;
	dnadata *dna;
	double **allones, **xx1, **xx2, **xx3;
	long istart, ifinish;
	
	// DEBUG
	double temptotal;
	
	dna = data->dnaptr;
	numslice = (op->panel) ? NUMSLICE : 1L;
	siteptr = data->siteptr;
	
	if (op->datatype == 'n') {
		if (start == FLAGLONG) { /* extra sites */
			istart = getdata_nummarkers(op,data);
			ifinish = istart+NUMINVAR-1;
		} else {
			findsubtreemarkers(op,data,start,finish,&istart,&ifinish);
			if(istart==FLAGLONG) return;
		}
	} else {
		istart = start;
		ifinish = finish;
	}
	
	/* allocate some working space--this should probably be moved
     upstream at some point to save time! */
	
	allones = (double **)calloc(numslice,sizeof(double *));
	xx3 = (double **)calloc(numslice,sizeof(double *));
	allones[0] = (double *)calloc(numslice*4,sizeof(double));
	xx3[0] = (double *)calloc(numslice*4,sizeof(double));
	for (i = 0; i < numslice; i++) {
		allones[i] = allones[0] + i * 4;
		xx3[i] = xx3[0] + i * 4;
		for (j = 0; j < 4; j++)
			allones[i][j] = 1.0;
	}
	
	/* are these sites active and on an upwards branch?  */
	qactive = (q) ? TRUE : FALSE;
	ractive = (r) ? TRUE : FALSE;
	
	
	for(i = istart; i <= ifinish; i++) {
		
		if(siteptr[i] > 0) { /* is an alias site available? */
			memcpy(p->x->s[i][0][0], p->x->s[siteptr[i]-1][0][0],
				   op->categs*numslice*4L*sizeof(double));
			// DEBUG
			temptotal = p->x->s[i][0][0][baseA] + p->x->s[i][0][0][baseC] + p->x->s[i][0][0][baseG] + p->x->s[i][0][0][baseT];
			if (temptotal == 0.0) 
				printf("bad likelihood in memcpy!");
			//
			continue;  /* go to next site */
		}
		
		for (j = 0; j < numcategs; j++) {
			memcpy(xx3[0], allones[0], numslice * 4L * sizeof(double));
			if (qactive) {
				ww1zz1 = tbl[j].ww1zz1;
				vv1zz1 = tbl[j].vv1zz1;
				yy1 = 1.0 - tbl[j].zz1;
				xx1 = q->x->s[i][j];
				
				for (k=0; k<numslice; k++) {
					sum1 = yy1 * (dna->freqa * xx1[k][baseA] + 
								  dna->freqc * xx1[k][baseC] +
								  dna->freqg * xx1[k][baseG] + dna->freqt * xx1[k][baseT]);
					sumr1 = dna->freqar * xx1[k][baseA] + dna->freqgr * xx1[k][baseG];
					sumy1 = dna->freqcy * xx1[k][baseC] + dna->freqty * xx1[k][baseT];
					vv1zz1_sumr1 = vv1zz1 * sumr1;
					vv1zz1_sumy1 = vv1zz1 * sumy1;
					xx3[k][baseA] *= 
					(sum1 + ww1zz1 * xx1[k][baseA] + vv1zz1_sumr1);
					xx3[k][baseC] *=
					(sum1 + ww1zz1 * xx1[k][baseC] + vv1zz1_sumy1);
					xx3[k][baseG] *=
					(sum1 + ww1zz1 * xx1[k][baseG] + vv1zz1_sumr1);
					xx3[k][baseT] *=
					(sum1 + ww1zz1 * xx1[k][baseT] + vv1zz1_sumy1);
					// DEBUG
					temptotal = xx3[k][baseA] + xx3[k][baseC] + xx3[k][baseG] + xx3[k][baseT];
					if (temptotal == 0.0) 
						printf("bad likelihood in q! - 4185\n");
					//
				}
			}
			
			if (ractive) {
				ww2zz2 = tbl[j].ww2zz2;
				vv2zz2 = tbl[j].vv2zz2;
				yy2 = 1.0 - tbl[j].zz2;
				xx2 = r->x->s[i][j];
				
				for (k=0; k<numslice; k++) {
					sum2 = yy2 * (dna->freqa * xx2[k][baseA] + 
								  dna->freqc * xx2[k][baseC] +
								  dna->freqg * xx2[k][baseG] + dna->freqt * xx2[k][baseT]);
					sumr2 = dna->freqar * xx2[k][baseA] + dna->freqgr * xx2[k][baseG];
					sumy2 = dna->freqcy * xx2[k][baseC] + dna->freqty * xx2[k][baseT];
					vv2zz2_sumr2 = vv2zz2 * sumr2;
					vv2zz2_sumy2 = vv2zz2 * sumy2;
					xx3[k][baseA] *= 
					(sum2 + ww2zz2 * xx2[k][baseA] + vv2zz2_sumr2);
					xx3[k][baseC] *= 
					(sum2 + ww2zz2 * xx2[k][baseC] + vv2zz2_sumy2);
					xx3[k][baseG] *= 
					(sum2 + ww2zz2 * xx2[k][baseG] + vv2zz2_sumr2);
					xx3[k][baseT] *= 
					(sum2 + ww2zz2 * xx2[k][baseT] + vv2zz2_sumy2);
					// DEBUG
					temptotal = xx3[k][baseA] + xx3[k][baseC] + xx3[k][baseG] + xx3[k][baseT];
					if (temptotal == 0.0) 
						printf("bad likelihood in r! - 4215\n");
					//
				}
			}
			
			memcpy(p->x->s[i][j][0], xx3[0], numslice * 4L * sizeof(double));
		}
	}
	
	free(allones[0]);
	free(allones);
	free(xx3[0]);
	free(xx3);
	
} /* calcrange */


void localsmooth(option_struct *op, data_fmt *data, node *p,
				 long indexsite, long startsite, long endsite)
{
	/* never mind if it's a tip */
	if (istip(p)) return;
	
	/* recurse left */
	if (!p->next->top) 
		if (!p->next->back->updated)
			if (inrange(p->next->back->ranges,indexsite))
				localsmooth(op,data,p->next->back,indexsite,startsite,endsite);
	
	/* recurse right */
	if (!p->next->next->top)
		if (!p->next->next->back->updated)
			if (inrange(p->next->next->back->ranges,indexsite))
				localsmooth(op,data,p->next->next->back,indexsite,startsite,
							endsite);
	
	/* update likelihoods */
	if (!p->updated && iscoal(p))
		switch(op->datatype) {
			case 'a':
				break;
			case 'b':
			case 'm':
				nuview_micro(op,data,p,indexsite,startsite,endsite);
				break;
			case 'n':
			case 's':
				nuview(op,data,p,indexsite,startsite,endsite);
				break;
			default:
				fprintf(ERRFILE,"localsmooth:unknown datatype %c\n",op->datatype);
				fprintf(ERRFILE,"only the types a, m, s");
				fprintf(ERRFILE,"(electrophoretic alleles,\n");
				fprintf(ERRFILE,"microsatellite data, sequence data)");
				fprintf(ERRFILE,"are allowed.\n");
				exit(-1);
		}
	
} /* localsmooth */


void snpsmooth (option_struct *op, data_fmt *data, tree *tr,
				long indexsite, long startsite, long endsite)
{
	tlist *t;
	
	for(t = tr->tymelist->succ; t != NULL; t = t->succ)
		if (iscoal(t->eventnode))
			nuview(op,data,t->eventnode,indexsite,startsite,endsite);
	
} /* snpsmooth */


void localeval(option_struct *op, data_fmt *data, node *p, boolean first)
{
	long subtree, *subtree_ranges, substart, subend, lastmarker, *siteptr;
	double dummy, llike;
	dnadata *dna;
	msatdata *ms;
	
	dna = data->dnaptr;
	ms = data->msptr;
	
	subtree_ranges = (long *)calloc(2*curtree->numrecombs+4,sizeof(long));
	subtree_ranges[0] = 1;
	
	siteptr = data->siteptr;
	
	findsubtrees(op,data,curtree->tymelist,subtree_ranges);
	
	llike = 0.0;
	for(subtree = 0; subtree < subtree_ranges[0]; subtree++) {
		substart = subtree_ranges[2*subtree+1];
		subend = subtree_ranges[2*subtree+2];
		/* do the regular sites */
		localsmooth(op,data,p,substart,substart,subend);
		switch(op->datatype) {
			case 'a':
				break;
			case 'b':
			case 'm':
				/* WARNING this throws away llike, which must be wrong! */
				llike += micro_evaluate(op,ms,curtree,substart,subend);
				break;
			case 'n':
				lastmarker = getdata_nummarkers(op,data)-1;
				snpsmooth(op,data,curtree,substart,FLAGLONG,FLAGLONG);
				if (op->panel)
					llike += panel_eval_calcrange(op,data,curtree,first,substart,
												  subend,op->categs);
				else
					llike += snp_eval_calcrange(op,data,curtree,first,substart,
												subend,op->categs);
				dummy = snp_evaluate(op,data,curtree,llike,substart,subend);
				break;
			case 's':
				llike += eval_calcrange(op,data,curtree,first,substart,subend,
										op->categs);
				dummy = evaluate(op,data,curtree,llike);
				break;
			default:
				fprintf(ERRFILE,"\nlocaleval--impossible datatype %c\n",
						op->datatype);
				break;
		}
	}
	
	traverse_unflag(curtree,curtree->nodep[1]);
	free(subtree_ranges);
	
} /* localeval */


/********************************************************************
 * micro_evaluate() calculates the tree data likelihood at the root *
 * assuming that the tree has been properly "smoothed" first.       */
double micro_evaluate(option_struct *op, msatdata *ms, tree *tr, 
					  long start, long end)
{
	long site, a;
	double term;
	node *nn;
	
	nn = tr->root->back;
	term = 0.0;
	
	for(site = start; site <= end; site++)
		for(a = 0; a < MICRO_ALLELEMAX-1; a++)
			term += nn->x->a[site][a];
	
	if(term == 0.0) return(NEGMAX);
	else return(log(term) + nn->lxmax);
	
} /* micro_evaluate */


/******************************************************************
 * countbranches returns the total number of branches in the tree *
 *    not including the root.                                     */
long countbranches(option_struct *op, data_fmt *data, tree *tr)
{
	
	return(tr->numrecombs + tr->numcoals * 2);
	
} /* countbranches */

/******************************************************************
 * testratio() returns TRUE to accept a tree, FALSE to reject it. *
 * ratiotype = "d" if called to decide a "drop"                   *
 *             "t" if called to decide a "twiddle"                *
 *             "f" if called to decide a "flip"                   */
boolean testratio(option_struct *op, data_fmt *data, tree *oldtree, 
				  tree *newtree, char ratiotype)
{
	double test, x, numoldbranches, numnewbranches;
	
	//test
	//;
	//return TRUE;
	
	if (israndom == 1) return TRUE;
	if (op->ctemp == 98) return hastingswithS(oldtree,newtree);
	if (op->ctemp == 99) return randomhottestratio(oldtree,newtree);
	if (op->ctemp == 100) return flatrecombinationrateratio(oldtree,newtree);
	
	if(newtree->likelihood == NEGMAX)
		return FALSE;  
	
	test = (newtree->likelihood - oldtree->likelihood);  
	
	//test 
	
	switch((int)ratiotype) {
		case 'd' :
			numoldbranches = (double)oldtree->numcoals * 2.0 + (double)oldtree->numrecombs;
			numnewbranches = (double)curtree->numcoals * 2.0 + (double)curtree->numrecombs;
			test += log(numoldbranches/numnewbranches);
			
			break;
		case 't' :
			//test += log((double)oldtree->numrecombs/(double)newtree->numrecombs);
			break;
		case 'f' :
			break;
		default :
			fprintf(ERRFILE,"ERROR:testratio encounters unknown type--no Hastings\n");
			break;
	}
	// heating
	// P(D|G')^1/T  P(G'|lamda)^1/T P(G|lamda)
	// ----------   --------------   ------
	// P(D|G)^1/T   P(G|lamda)^1/T  P(G'|lamda)
	
	if (op->ctemp != 1){
		/*     test += hmmlikelihood(newtree) - hmmlikelihood(oldtree); */
		/*     test += hmmlikelihood_with_temp(oldtree, temperature) - hmmlikelihood_with_temp(newtree,temperature); */
		/*     test += hmmlikelihood_two(newtree,oldgrec, oldgweight, temperature) - hmmlikelihood_two(oldtree,newgrec, newgweight,temperature); */
		test = test * (1.0/op->ctemp);
		
		//test += hmmlikelihood(newtree) - hmmlikelihood(oldtree);
		//if (op->ctemp == 101) test += hmmlikelihood_two_with_bigger_rec(newtree,oldgrec, oldgweight) - hmmlikelihood_two_with_bigger_rec(oldtree,newgrec,newgweight);
		//else  test +=  hmmlikelihood_two(newtree,oldgrec, oldgweight, op->ctemp) - hmmlikelihood_two(oldtree,newgrec, newgweight,op->ctemp);
	}
	
	if (test >= 0.0)    return TRUE;
	else {
		x = log(randum());
		if (x <= test)     return TRUE;
		
		else       return FALSE;
	}
} /* testratio */


void seekch(char c) /* use only in reading file intree! */
{
	if (gch == c)
		return;
	do {
		if (eoln(intree)) {
			fscanf(intree, "%*[^\n]");
			getc(intree);
		}
		gch = getc(intree);
		if (gch == '\n')
			gch = ' ';
	} while (gch != c);
}  /* seekch */

void getch(char *c) /* use only in reading file intree! */
{
	/* get next nonblank character */
	do {
		if (eoln(intree)) {
			fscanf(intree, "%*[^\n]");
			getc(intree);
		}
		*c = getc(intree);
		if (*c == '\n')
			*c = ' ';
	} while (*c == ' ');
}  /* getch */

void processlength(node *p)
{
	long digit;
	double valyew, divisor;
	boolean pointread;
	
	pointread = FALSE;
	valyew = 0.0;
	divisor = 1.0;
	getch(&gch);
	digit = gch - '0';
	while (((unsigned long)digit <= 9) || gch == '.'){
		if (gch == '.')
			pointread = TRUE;
		else {
			valyew = valyew * 10.0 + digit;
			if (pointread)
				divisor *= 10.0;
		}
		getch(&gch);
		digit = gch - '0';
	}
	p->length = valyew / divisor;
	p->back->length = p->length;
}  /* processlength */

void addelement(option_struct *op, data_fmt *data, node *p, long *nextnode)
{
	node *q;
	long i, n;
	boolean found;
	char str[NMLNGTH];
	
	getch(&gch);
	if (gch == '(') {
		(*nextnode)++;
		newnode(&q);
		q = curtree->nodep[(*nextnode)];
		hookup(p, q);
		addelement(op,data,q->next,nextnode);
		seekch(',');
		addelement(op,data,q->next->next, nextnode);
		seekch(')');
		getch(&gch);
	} else {
		for (i = 0; i < NMLNGTH; i++)
			str[i] = ' ';
		n = 1;
		do {
			if (gch == '_')
				gch = ' ';
			str[n - 1] = gch;
			if (eoln(intree)) {
				fscanf(intree, "%*[^\n]");
				getc(intree);
			}
			gch = getc(intree);
			if (gch == '\n')
				gch = ' ';
			n++;
		} while (gch != ':' && gch != ',' && gch != ')' && n <= NMLNGTH);
		n = 1;
		do {
			found = TRUE;
			for (i = 0; i < NMLNGTH; i++)
				found = (found && str[i] == curtree->nodep[n]->nayme[i]);
			if (!found)
				n++;
		} while (!(n > getdata_numseq(op,data) || found));
		if (n > getdata_numseq(op,data)) {
			printf("Cannot find sequence: ");
			for (i = 0; i < NMLNGTH; i++)
				putchar(str[i]);
			putchar('\n');
		}
		hookup(curtree->nodep[n], p);
	}
	if (gch == ':')
		processlength(p);
}  /* addelement */

void treeread(option_struct *op, data_fmt *data)
/* WARNING:  this only works with sequence data, not
 microsatellites, not panel SNPs, not allozymes! */
{
	long nextnode;
	node *p;
	
	curtree->root = curtree->nodep[rootnum];
	getch(&gch);
	if (gch == '(') {
		nextnode = getdata_numtips(op,data) + 1;
		p = curtree->nodep[nextnode];
		addelement(op,data,p, &nextnode);
		seekch(',');
		addelement(op,data,p->next, &nextnode);
		hookup(p->next->next, curtree->nodep[rootnum]);
		p->next->next->length = ROOTLENGTH;
		curtree->nodep[rootnum]->length = p->next->next->length;
		ltov(op,data,curtree->nodep[rootnum]);
	}
	fscanf(intree, "%*[^\n]");
	getc(intree);
}  /* treeread */


/* DEBUG debug WARNING warning --
 actually still used in traitlike.c:traitlike() */
void finddnasubtrees(option_struct *op, data_fmt *data,
					 tlist *tstart, long *sranges)
{
	long i, j, numsites, *temp;
	tlist *t;
	
	numsites = countsites(op,data);
	
	temp = (long *)calloc(numsites,sizeof(long));
	
	for(t = tstart, i = 1; t != NULL; t = t->succ)
		if(isrecomb(t->eventnode))
			temp[findlink(t->eventnode)] = 1;
	
	for(i = 0, j = 2, sranges[1] = 0; i < numsites; i++) {
		if(temp[i] == 0) continue;
		sranges[0]++;
		sranges[j] = i;
		sranges[j+1] = i+1;
		j+=2;
	}
	sranges[j] = numsites-1;
	sranges[j+1] = FLAGLONG;
	
	if (j+1 > 2*sranges[0]+3)
		fprintf(ERRFILE,"ERROR:finddnasubtree: j calculated wrong\n");
	
	free(temp);
	
} /* finddnasubtrees */


void findsubtrees_node(option_struct *op, data_fmt *data,
					   tlist *tlast, tree *tr, long **coal)
{
	long i, j, numsites, *temp;
	tlist *t;
	
	numsites = countsites(op,data);
	
	temp = (long *)calloc(numsites,sizeof(long));
	for(t = tr->tymelist, i = 1; t != tlast->succ; t = t->succ)
		if(isrecomb(t->eventnode)) {
			j = findlink(t->eventnode);
			if (temp[j]) continue;
			temp[j] = 1;
			i++;
		}
	
	init_coal_alloc(coal,i);
	
	for(i = 0, j = 2; i < numsites; i++) {
		if (temp[i] == 0) continue;
		(*coal)[j] = i;
		(*coal)[j+1] = i+1;
		j += 2;
	}
	(*coal)[j] = numsites-1;
	
	free(temp);
	
} /* findsubtrees_node */


void drop_findsubtrees(option_struct *op, data_fmt *data, tree *tr,
					   node *p)
{
	node *q;
	tlist *t;
	
	t = gettymenode(tr,p->number);
	
	findsubtrees_node(op,data,t,tr,&(p->coal));
	
	if (isrecomb(p)) {
		q = otherdtr(p);
		findsubtrees_node(op,data,t,tr,&(q->coal));
	}
	
} /* drop_findsubtrees */


void findsubtrees(option_struct *op, data_fmt *data, tlist *tstart,
				  long *sranges)
{
	long i, j, numsites, *temp;
	tlist *t;
	
	numsites = countsites(op,data);
	
	temp = (long *)calloc(numsites,sizeof(long));
	
	for(t = tstart; t != NULL; t = t->succ)
		if(isrecomb(t->eventnode))
			temp[findlink(t->eventnode)] = 1;
	
	for(i = 0, j = 2, sranges[1] = 0; i < numsites; i++) {
		if(temp[i] == 0) continue;
		sranges[0]++;
		sranges[j] = i;
		sranges[j+1] = i+1;
		j+=2;
	}
	sranges[j] = numsites-1;
	sranges[j+1] = FLAGLONG;
	
	if (j+1 > 2*sranges[0]+3)
		fprintf(ERRFILE,"ERROR:findsubtree: j calculated wrong\n");
	
	free(temp);
	
} /* findsubtrees */


/***************************************************************
 * markertosite() returns the exact "site" that corresponds to *
 * the passed marker position.                                 */
long markertosite(option_struct *op, data_fmt *data, long marker)
{
	
	return(getdata_markersite(op,data,marker));
	
} /* markertosite */


/**********************************************************************
 * sitetorightmarker() returns the closest marker to the right of the *
 * given site, if that marker exists.  If that marker doesn't exist,  *
 * FLAGLONG is returned.  If the site exactly corresponds to a marker *
 * then the corresponding marker will be returned.                    */
long sitetorightmarker(option_struct *op, data_fmt *data, long site)
{
	long marker, nummarkers;
	double msite;
	
	nummarkers = getdata_nummarkers(op,data);
	for(marker = 0; marker < nummarkers; marker++) {
		msite = markertosite(op,data,marker);
		if (msite >= site) return(marker);
	}
	
	return(FLAGLONG);
	
} /* sitetorightmarker */


/**********************************************************************
 * sitetomarker() returns the closest marker to the left of the given *
 * site, if that marker exists.  If that marker doesn't exist, the    *
 * closest marker on the right is returned.                           */
long sitetomarker(option_struct *op, data_fmt *data, long site)
{
	long marker, nummarkers;
	double sitecounts;
	
	nummarkers = getdata_nummarkers(op,data);
	for(marker = 0, sitecounts = 0.0; marker < nummarkers; marker++) {
		sitecounts = data->dnaptr->sitecount[locus][marker];
		if (site < sitecounts) return(marker);
	}
	
	fprintf(ERRFILE,"ERROR--failure to find marker for site %ld\n\n",
			site);
	return(FLAGLONG);
	
	
} /* sitetomarker */


/**********************************************************
 * findsubtreemarkers() gets passed the beginning and end *
 * of a subtree (in psuedosites) and returns the position *
 * of the first and last markers found in that subtree.   *
 * Return FLAGLONG in both marker positions if there are  *
 * no markers present in the subtree.                     */
void findsubtreemarkers(option_struct *op, data_fmt *data, long pstart,
						long pend, long *mstart, long *mend)
{
	long nummarkers, marker, site = 0;
	
	nummarkers = getdata_nummarkers(op,data);
	
	for(marker = 0; marker < nummarkers; marker++) {
		site = markertosite(op,data,marker);
		if (site >=  pstart) break;
	}
	
	if (site > pend) {
		(*mstart) = (*mend) = FLAGLONG;
		return;
	}
	
	(*mstart) = marker;
	
	for(;marker < nummarkers; marker++) {
		site = markertosite(op,data,marker);
		if (site > pend) break;
	}
	
	(*mend) = marker-1;
	
} /* findsubtreemarkers */


/**********************************************************
 * sameranges() returns TRUE if the ranges are identical, *
 * FALSE otherwise.                                       */
boolean sameranges(long *range1, long *range2)
{
	long i;
	boolean result = TRUE;
	
	if (range1[0] == range2[0] && range1[0] == 0) return result;
	
	if (range1[0] != range2[0]){
		result = FALSE;
		return result;
	}
	
	
	for(i = 1; i <= range1[0]; i++){
		if(range1[2 * i] != range2[2 * i] || range1[2 * i - 1] != range2[2 * i - 1]) {
			result = FALSE;
			return result;
		}
	}
	
	/*   for(i = 0; range1[i] != FLAGLONG && range2[i] != FLAGLONG; i++) */
	/*     if(range1[i] != range2[i]) return(FALSE); */
	
	/*   if (range1[i] != range2[i]) { */
	/*     printf("ERROR:sameranges bailed on FLAGLONG but not done yet!!! %ld %ld\n", */
	/* 	   indecks,apps); */
	/*     return(FALSE); */
	/*   } */
	
	return result;
	
} /* sameranges */


/*************************************************
 * copycoal() copies the coal array of a nodelet */
void copycoal(node *source, node *target)
{
	
	if (source->coal != NULL) {
		coal_Malloc(target,TRUE,source->coal[0]);
		memcpy(target->coal,source->coal,(source->coal[0]*2+2)*sizeof(long));
	} else coal_Malloc(target,FALSE,0L);
	
} /* copycoal */


/*****************************************************
 * copyranges() copies the ranges array of a nodelet */
void copyranges(node *source, node *target)
{
	
	if (source->ranges != NULL) {
		ranges_Malloc(target,TRUE,source->ranges[0]);
		memcpy(target->ranges,source->ranges,(source->ranges[0]*2+2)*sizeof(long));
	} else ranges_Malloc(target,FALSE,0L);
	
} /* copyranges */


/***************************************************************
 * This function copies the likelihood, "x" and "z", arrays of *
 * a nodelet.                                                  */
void copylikes(option_struct *op, data_fmt *data, node *source,
			   node *target)
{
	long nummarkers = 0, numslice;
	
	numslice = (op->panel) ? NUMSLICE : 1L;
	
	if (op->map) {
		if (source->z != NULL) {
			allocate_z(op,data,target);
			memcpy(target->z,source->z,NUMTRAIT*sizeof(double));
		} else free_z(op,target);
	}
	
	if (source->x != NULL) {
		//allocate_x(op,data,target);
		switch (op->datatype) {
			case 'a':
				break;
			case 'b':
			case 'm':
				memcpy(target->x->a[0],source->x->a[0],
					   getdata_numloci(op,data)*MICRO_ALLELEMAX*sizeof(double));
				break;
			case 'n':
				nummarkers += NUMINVAR;
			case 's':
				nummarkers += getdata_nummarkers(op,data);
				memcpy(target->x->s[0][0][0],source->x->s[0][0][0],
					   nummarkers*op->categs*numslice*4L*sizeof(double));
				break;
			default:
				fprintf(ERRFILE,"\ncopylikes--impossible datatype %c\n",
						op->datatype);
				exit(-1);
				break;
		}
	} else free_x(op,target);
	
} /* copylikes */


/************************************************************
 * copynode copies the "source" node onto the "target" node */
void copynode(option_struct *op, data_fmt *data, node *source,
			  node *target)
{
	long i, nodecount;
	
	if (istip(source)) nodecount = 1;
	else nodecount = 3;
	
	for (i = 0; i < nodecount; i++) {
		target->type = source->type;
		/* NEVER! target->next := source->next; */
		target->back = source->back;
		target->top = source->top;
		/* but NOT target->number := source->number; */
		//copylikes(op,data,source,target);
		copyranges(source,target);
		
		if (op->fc) copycoal(source,target);
		//if (istip(source)) memcpy(target->nayme,source->nayme,sizeof(source->nayme));
		target->v = source->v;
		target->lxmax = source->lxmax;
		target->tyme = source->tyme;
		target->length = source->length;
		target->updated = source->updated;
		target->recstart = source->recstart;
		target->recend = source->recend;
		source = source->next;
		target = target->next;
	}
}  /* copynode */

/**********************************************************************
 * addnode hooks nodelet "p" to node "q" by free back ptrs in p and q */
void addnode(option_struct *op, data_fmt *data, node *p, node *q)
{
	
	if (p->top) while(q->top && q->back) q = q->next;
	else while(!q->top && q->back) q = q->next;
	
	hookup(p,q);
	fixlength(op,data,p);
	
} /* addnode */


/****************************************************************
 * getrcnodelet finds which nodelet in node "copy" corresponds *
 * to the exact nodelet pointed to by "orig".                   */
node *getrecnodelet (node *copy, node *orig)
{
	long i;
	node *p;
	
	p = copy;
	for(i = 0; i < 3; i++) {
		if ((p->recstart == orig->recstart) && (p->recend == orig->recend))
			return(p);
		p = p->next;
	}
	
	fprintf(ERRFILE,"ERROR:getrecnodelet failed to match\n");
	return(NULL);
	
} /* getrecnodelet */


/**********************************************************
 * make_tree_copy copies tree "source" into tree "target" */
void make_tree_copy(option_struct *op, data_fmt *data, tree *source,
					tree *target)
{
	long i, nodenumber, *tnum, numno, numtips;
	tlist *t, *addhere, *new;
	node *p, *q, *r;
	dnadata *dna;
	
	dna = data->dnaptr;
	numtips = getdata_numtips(op,data);
	numno = numtips + source->numcoals + source->numrecombs + 1;
	tnum = (long *)calloc(1,numno*sizeof(long));
	
	/* make the tips */
	target->nodep = (node **)calloc(1,numno*sizeof(node *));
	for (i = 0; i < numtips+1; i++) {
		allocate_tip(op,data,target,source->nodep[i]->number);
		strcpy(target->nodep[i]->nayme,source->nodep[i]->nayme);
		if (i) {
			target->nodep[i]->top = TRUE;
			target->nodep[i]->tyme = 0.0;
		} else {
			target->nodep[i]->top = FALSE;
			target->nodep[i]->tyme = ROOTLENGTH;
		}
		target->nodep[i]->updated = TRUE;
		copyranges(source->nodep[i],target->nodep[i]);
		if(op->map) copytraits(source->nodep[i],target->nodep[i]);
		if(op->fc) copycoal(source->nodep[i],target->nodep[i]);
		//copylikes(op,data,source->nodep[i],target->nodep[i]);
		tnum[i] = i;
	}
	
	/* now construct the first tymelist entry */
	newtymenode(&target->tymelist);
	target->tymelist->numbranch = numtips;
	target->tymelist->branchlist =
    (node **)calloc(numtips,sizeof(node *));
	for(i = 0; i < numtips; i++)
		target->tymelist->branchlist[i] = target->nodep[i+1];
	target->tymelist->eventnode = target->nodep[1];
	target->tymelist->age = source->root->tyme;
	
	/* now do the rest of the tree and tymelist */
	for (t = source->tymelist->succ; t != NULL; t = t->succ) {
		nodenumber = t->eventnode->number;
		newnode(&p);
		p->number = nodenumber;
		p->next->number = nodenumber;
		p->next->next->number = nodenumber;
		target->nodep[nodenumber] = p;
		copynode(op,data,t->eventnode,p);
		p->back = NULL;
		p->next->back = NULL;
		p->next->next->back = NULL;
		tnum[t->eventnode->number] = nodenumber;
	}
	
	addhere = target->tymelist;
	for (t = source->tymelist->succ; t != NULL; t = t->succ) {
		p = findunique(target->nodep[tnum[t->eventnode->number]]);
		q = findunique(t->eventnode);
		if (isrecomb(p)) {
			r = target->nodep[tnum[q->back->number]];
			if (isrecomb(r)) r = getrecnodelet(r,q->back);
			addnode(op,data,p,r);
		} else {
			r = target->nodep[tnum[q->next->back->number]];
			if (isrecomb(r)) r = getrecnodelet(r,q->next->back);
			addnode(op,data,p->next,r);
			r = target->nodep[tnum[q->next->next->back->number]];
			if (isrecomb(r)) r = getrecnodelet(r,q->next->next->back);
			addnode(op,data,p->next->next,r);
		}
		addhere->update = 0;
		insertaftertymelist(addhere,p);
		addhere = addhere->succ;
	}
	new = target->tymelist;
	for (t = source->tymelist; t != NULL; t = t->succ){
		new->segments = copyseglist(t->segments);
		new = new->succ;
	}
	addnode(op,data,target->nodep[0],
			target->nodep[tnum[source->root->back->number]]);
	target->root = target->nodep[0];
	
	for (t = target->tymelist->succ; t != NULL; t = t->succ){
		if (iscoal(t->eventnode)){
			if (t->eventnode->coal[0] == 0) t->eventnode->back = target->root;
		}
	}
	
	free(tnum);
	
} /* make_tree_copy */


/********************************************************************
 * copycreature() copies a single creature to the "target" pointer. *
 * The passed tree should be the tree associated with the "target"  *
 * creature.                                                        */
void copycreature(option_struct *op, data_fmt *data, creature *source,
				  creature *target, tree *tr)
{
	long i;
	
	target->numflipsites = source->numflipsites;
	target->flipsites = (long *)calloc(target->numflipsites,sizeof(long));
	memcpy(target->flipsites,source->flipsites,
		   (target->numflipsites)*sizeof(long));
	
	target->numhaplotypes = source->numhaplotypes;
	target->haplotypes = (node **)calloc(target->numhaplotypes,sizeof(node *));
	for(i = 0; i < target->numhaplotypes; i++) {
		target->haplotypes[i] = tr->nodep[source->haplotypes[i]->number];
	}
	
} /* copycreature */


/********************************************************
 * copycreatures() copies the creatures array of a tree */
void copycreatures(option_struct *op, data_fmt *data, tree *source,
				   tree *target)
{
	long cr, numcreatures;
	
	if (!op->haplotyping) return;
	
	numcreatures = getdata_numtips(op,data)/NUMHAPLOTYPES;
	target->creatures = (creature *)calloc(numcreatures,sizeof(creature));
	
	for(cr = 0; cr < numcreatures; cr++) {
		copycreature(op,data,&(source->creatures[cr]),
					 &(target->creatures[cr]),target);
	}
	
} /* copycreatures */


/**********************************************************************
 * copytree makes a copy from scratch of the "source" tree, returning *
 * a pointer to the copy that it makes.                               */
tree *copytree(option_struct *op, data_fmt *data, tree *source)
{
	long numsites;
	
	tree *target;
	
	numsites = countsites(op,data);
	target = (tree *)calloc(1,sizeof(tree));
	
	target->likelihood = source->likelihood;
	target->coalprob = source->coalprob;
	target->numcoals = source->numcoals;
	target->numrecombs = source->numrecombs;
	make_tree_copy(op,data,source,target);
	copycreatures(op,data,source,target);
	
	target->dlikelihood = (double *)calloc(numsites,sizeof(double));
	memcpy(target->dlikelihood, source->dlikelihood, numsites*sizeof(double));
	target->recrates = (double *)calloc(numsites,sizeof(double));
	target->weight_array= (double *)calloc(numsites,sizeof(double));
	target->numrec_array = (double *)calloc(numsites,sizeof(double));
	if (source->recrates != NULL) memcpy(target->recrates, source->recrates, numsites*sizeof(double));
	if (source->weight_array != NULL) memcpy(target->weight_array, source->weight_array, numsites*sizeof(double));
	if (source->numrec_array != NULL) memcpy(target->numrec_array, source->numrec_array, numsites*sizeof(double));
	
	return(target);
} /* copytree */


/***************************************************
 * freetree frees the "target" tree from memory */
void freetree(option_struct *op, data_fmt *data, tree *target)
{
	long i;
	tlist *t;
	
	if (op->haplotyping) {
		for(i = 0; i < getdata_numtips(op,data)/NUMHAPLOTYPES; i++) {
			free(target->creatures[i].haplotypes);
			free(target->creatures[i].flipsites);
		}
		free(target->creatures);
	}
	
	t = target->tymelist->succ;
	while(t != NULL) {
		freenodelet(op,data,t->eventnode->next->next);
		freenodelet(op,data,t->eventnode->next);
		freenodelet(op,data,t->eventnode);
		t = t->succ;
	}
	
	t = target->tymelist->succ;
	while(t != NULL){
		if (t->segments != NULL){
			freesegment(t->segments);
			t->segments = NULL;
			t = t->succ;
		}
	}
	
	freetymelist(target->tymelist);
	
	for (i = 0; i < getdata_numtips(op,data)+1; i++) {
		freenodelet(op,data,target->nodep[i]);
	}
	free(target->nodep);
	
	free(target->dlikelihood);
	free(target->recrates);
	free(target->weight_array);
	free(target->numrec_array);
	
	free(target);
	
} /* freetree */

/* joinnode and constructtree are used for constructing a rather bad
 starting tree if the user doesn't provide one */
void joinnode(option_struct *op, data_fmt *data, double length,
			  node *p, node *q)
{
	hookup(p,q);
	p->length = length;
	q->length = length;
	ltov(op,data,p);
} /* joinnode */

void constructtree(option_struct *op, data_fmt *data, double branch)
{
	long i, j, numtips, nextnode;
	double height;
	node *p, *q;
	double branchlength;
	
	branchlength = theta0/getdata_numtips(op,data);
	numtips = getdata_numtips(op,data);
	
	curtree->root = curtree->nodep[rootnum];
	nextnode = numtips+1;
	p = curtree->root;
	q = curtree->nodep[nextnode];
	
	p->back = q;
	q->back = p;
	p->length = ROOTLENGTH;
	q->length = ROOTLENGTH;
	ltov(op,data,p);
	
	height = (numtips - 1) * branchlength;
	p->tyme = ROOTLENGTH + height;
	for (i = 1; i < numtips; i++) {
		p = curtree->nodep[i];
		q = curtree->nodep[nextnode]->next;
		joinnode(op,data,height,p,q);
		q = q->next;
		if (i != numtips-1) {
			nextnode++;
			p = curtree->nodep[nextnode];
			joinnode(op,data,branchlength,p,q);
			height -= branchlength;
		} else {
			p = curtree->nodep[numtips];
			joinnode(op,data,height,p,q);
		}
		for (j = 0; j < 3; j++)
			q->tyme = height;
	}
	
} /* constructtree */
/* End bad starting tree construction */


/******************************************************************
 * min_theta_calc() calculates the minimum value which we wish to *
 * let the mutation-rate parameter, theta, attain before final    *
 * estimation.                                                    */
double min_theta_calc(option_struct *op, data_fmt *data, double th)
{
	double value;
	
	value = THETAMIN;
	
	return(value);
	
} /* min_theta_calc */


/***********************************************************
 * hasrec_() returns TRUE if the passed chain has any *
 * recombinations currently sampled, FALSE otherwise.      */
boolean hasrec_chain(option_struct *op, long chain)
{
	long i, chaintype, refchain;
	treerec *trii;
	
	refchain = REF_CHAIN(chain);
	chaintype = TYPE_CHAIN(chain);
	
	for(i = 0; i < op->numout[chaintype]; i++) {
		trii = &sum[locus][refchain][i];
		if (trii->numrecombs > 0) return(TRUE);
	}
	
	return(FALSE);
	
} /* hasrec_chain */


/***********************************************************************
 * chainendcheck() runs at the end of a chain and makes sure the chain *
 * had the following desired properties:                               *
 *                                                                     *
 *    --at least 1 sampled tree contains 1+ recombinations             *
 *      FAILURE = call addfractrecomb()                                */
void chainendcheck(option_struct *op, data_fmt *data, tree *tr,
				   long chain, boolean locusend)
{
	
	if(!locusend && !hasrec_chain(op,chain) && op->holding != 2) {
		if (apps == totchains - 1) {
			fprintf(outfile,"WARNING--forced a recombination");
			fprintf(outfile," which may have affected final result\n");
		}
		addfractrecomb(op,data,chain,sum);
	}
	
} /* chainendcheck */


/******************************************************************
 * rearrange() is the driver for basic rearrangement of a tree,   *
 * returning TRUE if it successfully rearranged a tree, and FALSE *
 * otherwise.                                                     */
boolean rearrange(option_struct *op, data_fmt *data, long whichtchain)
{
	double chance;
	boolean accepted;
	
	curtree = temptrees[whichtchain];
	chance = randum();
	accepted = FALSE;
	op->ctemp = op->temperature[whichtchain];
	
	if (chance < TWIDDLE_PROB) accepted = twiddle(op,data);
	else {
		chance -= TWIDDLE_PROB;
		if (chance < op->happrob && op->haplotyping) {
			if (op->hapdrop == 0) accepted = fliphap(op,data,curtree);
			else if (op->hapdrop == 1) accepted = flipdrop(op,data,curtree,1L);
			else accepted = flipdrop(op,data,curtree,2L);
			hap++;
			if (accepted && whichtchain == 0) hacc++;
		}
		else {
			accepted = makedrop(op,data);
			slid++;
			if (accepted && whichtchain == 0) slacc++;
		}
	}
#ifdef MAC
	eventloop();
#endif
	
	temptrees[whichtchain] = curtree;
	
	return(accepted);
	
} /* rearrange */


/**********************************************************************
 * temptreeswap() checks to see if a swap of trees between different  *
 * temperatures is called for, executing the swap if so.              *
 *                                                                    *
 * it accepts a proposed swap using                                   *
 *    [ P(Gh|Tl) * P(D|Gh) ] ^ 1/l * [ P(Gl|Th) * P(D|Gl) ] ^ 1/h     *
 *    -----------------------------------------------------------     *
 *    [ P(Gh|Th) * P(D|Gh) ] ^ 1/h * [ P(Gl|Tl) * P(D|Gl) ] ^ 1/l     *
 *                                                                    *
 * where Gh = the tree generated at high temp                         *
 *       Gl = the tree generated at low temp                          *
 *       Tl = the parameter values at the low temp                    *
 *       Th = the parameter values at the high temp                   */
void temptreeswap(option_struct *op, data_fmt *data, boolean *changed,
				  long chain)
{
	long chl, chh;
	double chance, num, denom, templ, temph;
	boolean picked;
	tree *trh, *trl;
	
	if (op->numtempchains < 2) return;
	
	/* first pick which two chains to swap */
	
	picked = FALSE;
	while(!picked) {
		chl = (long)(randum()*op->numtempchains);
		chh = (long)(randum()*op->numtempchains);
		if (chl == op->numtempchains || chh == op->numtempchains) continue;
		if (chh - chl > 0) {
			if (chh - chl == 1) picked = TRUE;
		} else {
			if (chl - chh == 1) picked = TRUE;
		}
	}
	
	trh = temptrees[chh];
	temph = op->temperature[chh];
	trl = temptrees[chl];
	templ = op->temperature[chl];
	
	num = (trh->coalprob + trh->likelihood)*(1.0/templ) +
    (trl->coalprob + trl->likelihood)*(1.0/temph);
	
	denom = (trh->coalprob + trh->likelihood)*(1.0/temph) +
    (trl->coalprob + trl->likelihood)*(1.0/templ);
	
	chance = num - denom;
	
	swap++;
	if (chance >= log(randum())) {
		temptrees[chh] = trl;
		temptrees[chl] = trh;
		changed[chh] = TRUE;
		changed[chl] = TRUE;
		swacc++;
	}
	
} /* temptreeswap */

void printnewhap(){
	FILE *newhapfile;
	int site, ind;
	newhapfile = fopen("simhap","w");
	for (ind = 1; ind <= new_hap; ind++){
		for (site = 0; site < seq_length; site++){
			if (newsimdata[ind][site]== 1)
				fprintf(newhapfile,"A");
			if (newsimdata[ind][site] == 2)
				fprintf(newhapfile,"C");
			if (newsimdata[ind][site] == 3)
				fprintf(newhapfile,"G");
			if (newsimdata[ind][site] == 4)
				fprintf(newhapfile,"T");
		}
		fprintf(newhapfile,"\n");
	}
	fclose(newhapfile);
	return;
}


void copyoldhap(data_fmt *data){
	int ind, site;
	char ****y;
	
	
	y = data->dnaptr->seqs;
	
	for (ind = 0; ind < old_hap; ind ++){
		for (site = 0; site < seq_length; site++){
			
			switch (y[0][0][ind][site]){
				case 'A':
					olddata[ind][site] = 1;
					break;
				case 'C':
					olddata[ind][site] = 2;
					break;
				case 'G':
					olddata[ind][site] = 3;
					break;
				case 'T':
					olddata[ind][site] = 4;
					break;
				case 'N':
					olddata[ind][site] = 0;
			}
		}
	}
	
	return;
}


void marksubtree(treeshape *subtree){ // let's use isroot part
	nodeinfo *cur;
	cur = subtree->rootnode->leftnode;
	while(1){
		cur->isroot = -1;
		if (cur->istip == 1){
			if (cur->nodenum != new_tip + old_hap)
				printf("error at mark subtree\n");
			return;
		}
		cur = cur->leftnode;
	}
	return;
}



#define PRINTNUM     8 /* the number of digits used for printing
the number of steps done */

void maketree(option_struct *op, data_fmt *data)
{
	long tempchain, incrprog, metout, progout, chaintype;
	double bestlike;
	boolean runend, *changetree, first;
	char chainlit[2][6] = {"Short","Long"};
	dnadata *dna;
	double *traitarray;
	long hapscore;
	tree *bestlikelihoodtree = NULL;
	
	FILE *resultfile;
	//tets
	//path *temp1,*temp2;
	
	traitarray = NULL;
	
	changetree = (boolean *)calloc(op->numtempchains,sizeof(boolean));
	
	if (op->map)
		traitarray = (double *)calloc(countsites(op,data),sizeof(double));
	
	chainlit[0][5] = chainlit[1][5] = '\0';
	
	/* WARNING DEBUG assumes DNA data (and shouldn't)! */
	dna = data->dnaptr;
	
	getc(infile);
	if(op->haplotyping) {
		fprintf(outfile,"Haplotypes being inferred");
		if (op->hapdrop == 0) fprintf(outfile," without resimulation\n");
		if (op->hapdrop == 1) fprintf(outfile," with single resimulation\n");
		if (op->hapdrop == 2) fprintf(outfile," with double resimulation\n");
	}
	fprintf(outfile,"Watterson estimate of theta is %12.8f\n", watttheta);
#if !DISTRIBUTION
	fprintf(outfile,"\n\nIntermediate chains (only)\n");
	fprintf(outfile,"chain#     theta     rec-rate\n");
	fprintf(outfile,"-----------------------------\n");
#endif
	if (op->usertree)
		treeread(op,data);
	else {
#if !ALWAYS_REJECT
		branch0 = watttheta/getdata_numtips(op,data);
		constructtree(op,data,branch0);
#else
		constructtree(op,data,0.0);
#endif
	}
	orient(curtree,op,data,curtree->root->back);
#if !ALWAYS_REJECT
	if (op->plump) plumptree(op,data,watttheta);
#endif
	finishsetup(op,data,curtree->root->back);
	initbranchlist(op,data);
#if !ALWAYS_REJECT
	//localeval(op,data,curtree->root->back,TRUE);
	// to be corrected !
	curtree->coalprob = coalprob(op,data,curtree,theta0,rec0);
#endif
	bestlike = NEGMAX;
	theti[locus][0] = theta0;
	reci[locus][0] = rec0;
	runend = FALSE;
	
	//construct seglists for each time list
	constructseglist(op,data,curtree);
	calcinitialtree(op,data,curtree, curtree->nodep[1]->coal);
	
	calcratio(op,data,curtree,curtree,curtree->tymelist->eventnode->coal);
	
	
	/* WARNING DEBUG MARY */
	if (TRUEHAP) {
		hapscore = hapdist(op,&truehaps,data,curtree->creatures);
		printf("Distance from truth%ld\n",hapscore);
		fprintf(outfile,"Distance from truth%ld\n",hapscore);
		hapscore = hapdist(op,&starthaps,data,curtree->creatures);
		fprintf(outfile,"Distance from start%ld\n",hapscore);
	}
	
	/*for test */
	op->hetero = TRUE;
	
	/* for test */
	
	
	if(op->map) traitread(curtree, getdata_numseq(op,data));
	/**********************************/
	/* Begin Hastings-Metropolis loop */
	/**********************************/
	
	
	//getinitialrecs(op,data);
	//testinitrees(op,data);
	
	// generate initial tree..
	findhotspots(op,data);
	//gets(apps);
	
	for (apps = 0; apps < totchains; apps++) {
		
		if (apps >= op->numchains[0]) chaintype = 1;
		else chaintype = 0;
		if (op->progress) {
			printf("%s chain %ld ",chainlit[chaintype],
				   ((chaintype) ? (apps + 1 - op->numchains[0]) : apps + 1));
			fflush(stdout);
		}
		numdropped = 0;
		//metout = op->increm[chaintype] - 1 + NUMBURNIN;
		metout = op->increm[chaintype] - 1 + numburnin_current;
		
		/* print a "." to stdout every 10% of the chain(s) finished */
		incrprog = (long)(op->steps[chaintype] / 10.0);
		//progout = incrprog - 1 + NUMBURNIN;
		progout = incrprog - 1 + numburnin_current;
		op->numout[chaintype] = 0;
		slacc = hacc = swacc = 0;
		slid = hap = swap = 0;
		
		/* if apps == 0, start all temperature chains with the same
		 initial tree */
		if (apps == 0) {
			for(tempchain = 0; tempchain < op->numtempchains; tempchain++) {
				temptrees[tempchain] = copytree(op,data,curtree);
			}
			
			/* curtree is no longer needed as a base tree!, it's just a ptr */
			freetree(op,data,curtree);
		}
		else{
			for(tempchain = 0; tempchain < op->numtempchains; tempchain++) {
				temptrees[tempchain] = copytree(op,data,bestlikelihoodtree);
			}
			freetree(op,data,bestlikelihoodtree);
			bestlikelihoodtree = NULL;
			bestlike = NEGMAX;
		}
		
		for(tempchain = 0; tempchain < op->numtempchains; tempchain++)
			changetree[tempchain] = FALSE;
		first = TRUE;
		
		curtree = temptrees[0];
		for (indecks=0; indecks < op->steps[chaintype]+numburnin_current; indecks++) {
			
			col = 0; /* column number, used in treeout */
			
			/* only count acceptances after burn-in */
			//if (indecks == NUMBURNIN)
			if (indecks == numburnin_current)
				hacc = slacc = 0;
			nsites = countsites(op,data);
			
			for (tempchain = 0; tempchain < op->numtempchains; tempchain ++){
				if (apps == 0 && indecks == numburnin_current - INIBURNIN && curtree->likelihood < bestlikelihoodtree->likelihood && curtree != bestlikelihoodtree){
					//freetree(op,data,curtree);
					temptrees[tempchain] = copytree(op,data,bestlikelihoodtree);
				}
				if (rearrangeH(op,data,tempchain,indecks)) changetree[tempchain] = TRUE;
			}
			
			
			temptreeswap_rec(op,data,changetree,apps);
			
			
			//printf(" %ld, %ld, %ld",t1accep,t2accep,t3accep);
			//printf(", swapped %ld/%ld\n",swacc,swap);
			if (apps == totchains - 1) /* end of run? */
				runend = TRUE;
			
			if (temptrees[0]->likelihood > bestlike) {
				if (bestlikelihoodtree != NULL) freetree(op,data,bestlikelihoodtree);
				bestlikelihoodtree = copytree(op,data,curtree);
				bestlike = temptrees[0]->likelihood;
				//printf("  like %lf\n",bestlike);
				if (ONEBESTREE) {
					FClose(bestree);
					bestree = fopen("bestree","w+");
					fprintf(bestree, "Chain #%2ld (%s) Step:%8ld\n",apps+1,
							chainlit[chaintype], indecks+1);
					/* debug DEBUG warning WARNING--no rectreeout routine */
					/*rec_outtree(temptrees[0]->root->back,TRUE, &bestree);*/
					fprintf(bestree, "; [%12.10f]\n", temptrees[0]->likelihood);
					bestlike = temptrees[0]->likelihood;
				}
				else {
					fprintf(bestree, "Chain #%2ld (%s) Step:%8ld\n",apps+1,
							chainlit[chaintype], indecks+1);
					/* debug DEBUG warning WARNING--no rectreeout routine */
					/*rec_outtree(temptrees[0]->root->back,TRUE, &bestree); */
					fprintf(bestree, "; [%12.10f]\n", temptrees[0]->likelihood);
					bestlike = temptrees[0]->likelihood;
				}
			}
			if (indecks == metout) {
				if (op->numout[chaintype] == 0)
					sametree[locus][0] = FALSE;
				else
					sametree[locus][op->numout[chaintype]] = !changetree[0];
				changetree[0] = FALSE;
				op->numout[chaintype]++;
				/* set curtree to the coldest tree for scoretree() scoring */
				curtree = temptrees[0];
				//scoretree(op,data,apps);
				scoretreeH(op,data,apps);
				if(op->map) traitlike(op,data,curtree,countsites(op,data),op->mutrait,op->traitratio, op->pd, traitarray);
				metout += op->increm[chaintype];
				if (op->treeprint) {
					fprintf(treefile,"\nlocus = %ld, chain = %ld, tree = %ld\n",
							locus,apps,indecks);
					/* debug DEBUG warning WARNING--no rectreeout routine */
					/*rec_outtree(temptrees[0]->root->back,TRUE, &treefile);*/
					fprintf(treefile,";\n");
				}
			}
			
			if (op->progress) {
				//if (!first) for(i = 0; i < PRINTNUM; i++) printf("\b");
				//else printf(" examined: ");
				printf("%*ld",(int)PRINTNUM,indecks-numburnin_current+1);
				fflush(stdout);
				first = FALSE;
			}
		}
		if(runend) {
			if(op->map) {
				traitprint(countsites(op,data),traitarray,
						   outfile,op->numout[chaintype]);
#if TRUTHKNOWN
				traitresult(op,data,traitarray,outfile,op->numout[chaintype]);
#endif
			}
		}
		printf(" trees");
		chainendcheck(op,data,curtree,apps,runend);
		//rec_estimatev(op,data,locus,apps,runend);
		//test
		//printf("ave like %lf\n",templike/10000.0);
		
		parameter_estimation_EM(op,apps,theta0, recrates,lamda);
		
		if (TRUEHAP) {
			hapscore = hapdist(op,&truehaps,data,curtree->creatures);
			printf("Distance from truth: %ld\n",hapscore);
			fprintf(outfile,"Distance from truth: %ld\n",hapscore);
			hapscore = hapdist(op,&starthaps,data,curtree->creatures);
			printf("Distance from start: %ld\n",hapscore);
			fprintf(outfile,"Distance from start: %ld\n",hapscore);
		}
		
		if (op->map && runend) traitsiteplot(op,data,traitarray);
		//theta0 = theti[locus][apps+1];
		//rec0 = reci[locus][apps+1];
		/*     if (min_theta_calc(op,data,theta0) > theta0 &&  */
		/* 	op->holding != 1) { */
		/*       if (apps == totchains - 2) { */
		/* 	fprintf(outfile,"WARNING--lower bound on estimate of"); */
		/* 	fprintf(outfile," theta may have affected final result\n"); */
		/*       } */
		/*       theta0 = min_theta_calc(op,data,theta0); */
		/*       theti[locus][apps+1] = theta0; */
		/*     } */
		if (op->progress) {
			printf("accepted %ld/%ld trees",slacc,slid-numburnin_current);
			if (numdropped) printf(", dropped %ld",numdropped);
			if (swap) printf(", swapped %ld/%ld",swacc,swap);
			if (op->haplotyping)
				printf("and %ld/%ld haplotype changes\n",hacc,hap);
			else printf("\n");
		}
		if ((apps == totchains-1) && numdropped) {
			fprintf(outfile,"%ld trees were dropped from",numdropped);
			fprintf(outfile," the final chain\n");
		}
		
	}
	//printbest(bestlikelihoodtree);
	resultfile = fopen("hmmresult","w");
	fprintf(resultfile,"theta %lf\n",theta0);
	fprintf(resultfile,"rec1 %lf lamda %lf\n",recrates[0], lamda[0]);
	fprintf(resultfile,"rec2 %lf lamda %lf\n",recrates[1], lamda[1]);
	fclose(resultfile);
	
	if(slacc == 0) {
		fprintf(outfile,"WARNING--no proposed trees ever accepted\n");
		fprintf(ERRFILE,"WARNING--no proposed trees ever accepted\n");
	}
	free(changetree);
}  /* maketree */

nodeinfo *findnewtip(nodeinfo *root){
	nodeinfo *cur;
	cur = root;
	while(1){
		cur = cur->leftnode;
		if (cur->istip == 1){
			if (cur->isroot == -1){
				return cur;
			}
			printf("problem at findnewtip");
		}
	}
}



void maketree_S(option_struct *op, data_fmt *data)
{
	long tempchain, incrprog, metout, progout, chaintype, i;
	tree *temptree;
	double bestlike;
	boolean runend, *changetree, first;
	char chainlit[2][6] = {"Short","Long"};
	dnadata *dna;
	double *traitarray;
	long hapscore;
	tree *bestlikelihoodtree = NULL;
	
	//tets
	//path *temp1,*temp2;
	
	traitarray = NULL;
	
	changetree = (boolean *)calloc(op->numtempchains,sizeof(boolean));
	
	if (op->map)
		traitarray = (double *)calloc(countsites(op,data),sizeof(double));
	
	chainlit[0][5] = chainlit[1][5] = '\0';
	
	/* WARNING DEBUG assumes DNA data (and shouldn't)! */
	dna = data->dnaptr;
	
	getc(infile);
	branch0 = theta0/getdata_numtips(op,data);
	constructtree(op,data,theta0/getdata_numtips(op,data));
	
	orient(curtree,op,data,curtree->root->back);
#if !ALWAYS_REJECT
	if (op->plump) plumptree(op,data,theta0);
#endif
	finishsetup(op,data,curtree->root->back);
	initbranchlist(op,data);
#if !ALWAYS_REJECT
	//localeval(op,data,curtree->root->back,TRUE);
	// to be corrected !
	curtree->coalprob = coalprob(op,data,curtree,theta0,rec0);
#endif
	
	
	op->numtempchains = 1;
	
	// copy old hap
	
	
	//construct seglists for each time list
	constructseglist(op,data,curtree);
	calcinitialtree(op,data,curtree, curtree->nodep[1]->coal);
	
	calcratio(op,data,curtree,curtree,curtree->tymelist->eventnode->coal);
	
	
	/* WARNING DEBUG MARY */
	if (TRUEHAP) {
		hapscore = hapdist(op,&truehaps,data,curtree->creatures);
		printf("Distance from truth%ld\n",hapscore);
		fprintf(outfile,"Distance from truth%ld\n",hapscore);
		hapscore = hapdist(op,&starthaps,data,curtree->creatures);
		fprintf(outfile,"Distance from start%ld\n",hapscore);
	}
	
	/*for test */
	op->hetero = TRUE;
	
	/* for test */
	
	
	
	for (apps = 0; apps < 1; apps++) {
		
		if (apps >= op->numchains[0]) chaintype = 1;
		else chaintype = 0;
		if (op->progress) {
			printf("Sampling genealogies for reference haplotypes\n");
			printf("MCMC chain\t Posterior\n");
			fflush(stdout);
		}
		numdropped = 0;
		//metout = op->increm[chaintype] - 1 + NUMBURNIN;
		metout = op->increm[chaintype] - 1 + numburnin_current;
		
		/* print a "." to stdout every 10% of the chain(s) finished */
		incrprog = (long)(op->steps[chaintype] / 10.0);
		//progout = incrprog - 1 + NUMBURNIN;
		progout = incrprog - 1 + numburnin_current;
		op->numout[chaintype] = 0;
		slacc = hacc = swacc = 0;
		slid = hap = swap = 0;
		
		op->numtempchains = 1;
		/* if apps == 0, start all temperature chains with the same
		 initial tree */
		if (apps == 0) {
			for(tempchain = 0; tempchain < op->numtempchains; tempchain++) {
				temptrees[tempchain] = copytree(op,data,curtree);
			}
			
			/* curtree is no longer needed as a base tree!, it's just a ptr */
			freetree(op,data,curtree);
		}
		else{
			for(tempchain = 0; tempchain < op->numtempchains; tempchain++) {
				temptrees[tempchain] = copytree(op,data,bestlikelihoodtree);
			}
			freetree(op,data,bestlikelihoodtree);
			bestlikelihoodtree = NULL;
			bestlike = NEGMAX;
		}
		
		for(tempchain = 0; tempchain < op->numtempchains; tempchain++)
			changetree[tempchain] = FALSE;
		first = TRUE;
		
		curtree = temptrees[0];
		for (indecks=0; indecks < numMCMC; indecks++){
			
			
			col = 0; /* column number, used in treeout */
			
			/* only count acceptances after burn-in */
			//if (indecks == NUMBURNIN)
			if (indecks == numburnin_current)
				hacc = slacc = 0;
			nsites = countsites(op,data);
			
			for (tempchain = 0; tempchain < op->numtempchains; tempchain ++){
				if (rearrangeH(op,data,tempchain,indecks)) changetree[tempchain] = TRUE;
			}
			
			
			temptreeswap_rec(op,data,changetree,apps);
    		
			
			//printf(" %ld, %ld, %ld",t1accep,t2accep,t3accep);
			//printf(", swapped %ld/%ld\n",swacc,swap);
			if (apps == totchains - 1) /* end of run? */
				runend = TRUE;
    		
			
			
			if (indecks == metout) {
				/*
				 if (op->numout[chaintype] == 0)
				 sametree[locus][0] = FALSE;
				 else
				 sametree[locus][op->numout[chaintype]] = !changetree[0];
				 changetree[0] = FALSE;
				 op->numout[chaintype]++;
				 */
				/* set curtree to the coldest tree for scoretree() scoring */
				curtree = temptrees[0];
				//scoretree(op,data,apps);
				//scoretreeH(op,data,apps);
				if(op->map) traitlike(op,data,curtree,countsites(op,data),op->mutrait,op->traitratio, op->pd, traitarray);
				metout += op->increm[chaintype];
				if (op->treeprint) {
					fprintf(treefile,"\nlocus = %ld, chain = %ld, tree = %ld\n",
							locus,apps,indecks);
					/* debug DEBUG warning WARNING--no rectreeout routine */
					/*rec_outtree(temptrees[0]->root->back,TRUE, &treefile);*/
					fprintf(treefile,";\n");
				}
			}
			
			if (op->progress) {
				//if (!first) for(i = 0; i < PRINTNUM; i++) printf("\b");
				//else printf(": ");
				//printf("%*ld",(int)PRINTNUM,indecks-numburnin_current+1);
				//fflush(stdout);
				//first = FALSE;
				if ((indecks)%1000 == 0){
					scorerecs_temp(temptrees[0]);
					temptrees[0]->coalprob = treeprobs(temptrees[0]);
					printf("%ld\t%lf\n",indecks,temptrees[0]->likelihood);
					//branchlength();
				}
			}
			if (indecks > minMCMC && indecks % repinter == 0){
				curtree = temptrees[0];
				
				if (filtered != 0){ // introduce new mutations
					
					simhaplotype_old_2(op,data);
				}
				temptree = copytree(op,data,curtree);
				printf("Simulation of new haplotypes\n");
				makedrop_S_new(op,data);	
				S_O = S_H = S_N = 0;
				simhaplotype(op,data);
				//printf("observed SNP\t%d\thidden SNP\t%d\tnew SNP\t%d\n",S_O,S_H,S_N);
				freetree(op,data,curtree);
				curtree = copytree(op,data,temptree);
				temptrees[0] = curtree;
				
				copyoldhap(data); // back to old haplotype
				constructseqinfo(op, data, 0);
				//curtree = copytree(op,data,temptree);
				//temptrees[0] = curtree;
				
			}
			
			
		}
		if(runend) {
			if(op->map) {
				traitprint(countsites(op,data),traitarray,
						   outfile,op->numout[chaintype]);
#if TRUTHKNOWN
				traitresult(op,data,traitarray,outfile,op->numout[chaintype]);
#endif
			}
		}
		printf(" trees");
		chainendcheck(op,data,curtree,apps,runend);
		//rec_estimatev(op,data,locus,apps,runend);
		//test
		//printf("ave like %lf\n",templike/10000.0);
		
		//parameter_estimation_EM(op,apps,theta0, recrates,lamda);
		
		if (TRUEHAP) {
			hapscore = hapdist(op,&truehaps,data,curtree->creatures);
			printf("Distance from truth: %ld\n",hapscore);
			fprintf(outfile,"Distance from truth: %ld\n",hapscore);
			hapscore = hapdist(op,&starthaps,data,curtree->creatures);
			printf("Distance from start: %ld\n",hapscore);
			fprintf(outfile,"Distance from start: %ld\n",hapscore);
		}
		
		if (op->map && runend) traitsiteplot(op,data,traitarray);
		//theta0 = theti[locus][apps+1];
		//rec0 = reci[locus][apps+1];
		/*     if (min_theta_calc(op,data,theta0) > theta0 &&  */
		/* 	op->holding != 1) { */
		/*       if (apps == totchains - 2) { */
		/* 	fprintf(outfile,"WARNING--lower bound on estimate of"); */
		/* 	fprintf(outfile," theta may have affected final result\n"); */
		/*       } */
		/*       theta0 = min_theta_calc(op,data,theta0); */
		/*       theti[locus][apps+1] = theta0; */
		/*     } */
		if (op->progress) {
			printf(" accepted %ld/%ld trees",slacc,slid-numburnin_current);
			if (numdropped) printf(", dropped %ld",numdropped);
			if (swap) printf(", swapped %ld/%ld",swacc,swap);
			if (op->haplotyping)
				printf("and %ld/%ld haplotype changes\n",hacc,hap);
			else printf("\n");
		}
		if ((apps == totchains-1) && numdropped) {
			fprintf(outfile,"%ld trees were dropped from",numdropped);
			fprintf(outfile," the final chain\n");
		}
		
	}
	
	
	// copy recombination array
	
	for (i = 0; i < seq_length; i++)
		newgrec[i] = curtree->numrec_array[i];
	
	
	printf("\n");
	//printf("New haplotype simulation ");
	//for (new_tip = 1; new_tip <= new_hap; new_tip++){
	
	
	
	//printnewhap();
	//printbest(bestlikelihoodtree);
	free(changetree);
}  /* maketree */


/************************************************************
 * freenodelet frees all fields of a nodelet and the actual *
 * nodelet too.                                             */
void freenodelet(option_struct *op, data_fmt *data, node *p)
{
	free_x(op,p);
	if(op->map)free_z(op,p);
	if (p->nayme) free(p->nayme);
	ranges_Malloc(p,FALSE,0L);
	coal_Malloc(p,FALSE,0L);
	if (p->oldbackranges != NULL){
		free(p->oldbackranges);
	}
	free(p);
} /* freenodelet */

/*******************************************************************
 * finalfree() frees                sum, reci, sametree,           *
 * theti,           the options structure, and the data structure. */
void finalfree(option_struct *op, data_fmt *data)
{
	long i, j, k, chtype, numloci;
	
	numloci = 0;
	
	numloci = getdata_numloci(op,data);
	/*
	 free(reci[0]);
	 free(reci);
	 free(theti[0]);
	 free(theti);
	 free(sametree[0]);
	 free(sametree);
	 */	
	for(i = 0; i < numloci; i++)
		for(j = 0; j < 1+op->numchains[1]; j++) {
			if (j) chtype = 1;
			else chtype = 0;
			for(k = 0; k < op->numout[chtype]; k++) {
#if GROWTHUSED
				free(sum[i][j][k].kk);
				free(sum[i][j][k].kend);
				free(sum[i][j][k].actives);
#endif
				free(sum[i][j][k].eventtype);				
				
				
				free(sum[i][j][k].sitescore);
				free(sum[i][j][k].rec);
				free(sum[i][j][k].sumsite);
			}
		}
	/*
	 free(sum[0][0]);
	 free(sum[0]);
	 free(sum);
	 */	
	freedata(data);
	
	free(op->probcat);
	free(op->rate);
	free(op->temperature);
	if (!op->same_ne) free(op->ne_ratio);
	if (op->panel) free(op->numpanel);
	free(op);
	
} /* finalfree */

void expsnp(){
	double sumh = 0;
	int i;	
	
	for (i = 1; i < old_hap; i++){
		sumh = sumh + 1.0/(double)(i);
	}
	sumh = sumh * seq_length * theta0;
	expsnpnum = (int)(sumh - num_markers);
	if (expsnpnum < 0){
		expsnpnum = 0;
	}
	return;
}
void expsnp_freq(){ // using freq
	double sum1 = 0;
	double sum2 = 0;
	int i;
	int cutoffnum = 0;
	
	cutoffnum = (int)(cutofffreq * old_hap);
	
	for (i = 1; i <= cutoffnum; i++){
		sum1 = sum1 + 1/(double)(i);
	}
	for (i = cutoffnum + 1; i <= old_hap; i++){
		sum2 = sum2 + 1/(double)(i);
	}
	freqw1 = sum1/(sum1 + sum2);
	freqw2 = 1 - freqw1;
	expsnpnum = (int)(sum1 * num_markers / sum2);
	
}



int main(int argc, char *argv[])
{  /* Recombine */
	long i, numloci, numpop;
	data_fmt data;
	option_struct *options;
	int ind;
	char infilename[100], outfilename[100], logfilename[100];	
#ifdef MAC
	char spafilename[100], seedfilename[100];
	
	argv[0] = "Recombine";
#endif
	
	
	/* Open various filenames. */
	
	openfile(&infile,"refhap","r",argv[0],infilename);
	openfile(&simlog,"simlog","w",argv[0],logfilename);
	openfile(&outfile,"debug","w",argv[0],outfilename);
	
	/* warning WARNING debug DEBUG--initialization handled wrong */
	population = 0;
	
	setupoption_struct(&options);
	options->ibmpc = IBMCRT;
	options->ansi = ANSICRT;
	getoptions(options);
	
	for (i = 1; i <= 1000; i++)
		clearseed = randum();
	
	
	//if (options->newdata) getinput(options,&data);
	//else getnodata(options,&data);
	getsimoptions(options,&data);	
	firstinit(options,&data);
	
	
	//constant = (double **)calloc(options->numtempchains,sizeof(double));
	
	seq_length = countsites(options,&data);
	
	old_hap = getdata_numseq(options,&data);
	
	
	olddata = (short int  **)calloc((old_hap + 2), sizeof(short int *));
	for (ind = 0; ind < old_hap + 2; ind++){
		olddata[ind] = (short int  *)calloc(seq_length, sizeof(short int));
	}
	copyoldhap(&data);
	dna_sequence = (char *)calloc(seq_length,sizeof(char));
	
	
	newgrec = (int *)calloc(seq_length, sizeof(int));
	oldgrec = (int *)calloc(seq_length, sizeof(int));
	
	newgweight = (double *)calloc(seq_length, sizeof(double));
	oldgweight = (double *)calloc(seq_length, sizeof(double));
	
	gpartweight = (double *)calloc(seq_length,sizeof(double));
	gpartrec = (long *)calloc(seq_length,sizeof(long));
	
	
	
	tiplikelist = (double *)calloc(getdata_numtips(options,&data),sizeof(double));
	
	if (options->newdata) {
		temptrees = (tree **)calloc(options->numtempchains,sizeof(tree *));
		numpop = getdata_numpop(options,&data);
		for(population = 0; population < numpop; population++) {
			popinit(options,&data);
			numloci = getdata_numloci(options,&data);
			for (locus = 0; locus < numloci; locus++) {
				//if (options->progress) printf("Locus %ld\n",locus+1);
				
				
				locusinit(options,&data);
				
				seq_num = getdata_numseq(options, &data);
				constructseqinfo(options, &data,1);
				
				expsnp_freq();
				maketree_S(options, &data);
				end_of_locus_free(options,&data);
			}
			locus--;
			//if (numloci > 1) rec_estimate(options,&data,-1L,totchains-1,TRUE);
			end_of_population_free(options,&data);
		}
		free(temptrees);
	} else {
		rec_scoreread(options,&data);
		numloci = getdata_numloci(options,&data);
		if (numloci > 1) rec_estimatev(options,&data,-1L,totchains-1,TRUE);
		else rec_estimatev(options,&data,0L,totchains-1,TRUE);
	}
	
	finalfree(options,&data);
	
	FClose(infile);
	FClose(outfile);
	FClose(treefile);
	FClose(bestree);
	//FClose(seedfile);
	FClose(simlog);
	FClose(spacefile);
	
#ifdef MAC
	strcpy(seedfilename,"seedfile");
	strcpy(spafilename,"spacefile");
	fixmacfile(outfilename);
	fixmacfile(infilename);
	fixmacfile(bestrfilename);
	fixmacfile(logfilename);
	fixmacfile(intrfilename);
	fixmacfile(seedfilename);
	fixmacfile(trfilename);
	fixmacfile(spafilename);
#endif
	
	printf("PROGRAM DONE\n");
	exit(0);
	
}  /* recombine */

int eof(FILE *f)
{
	register int ch;
	
	if (feof(f))
		return 1;
	if (f == stdin)
		return 0;
	ch = getc(f);
	if (ch == EOF)
		return 1;
	ungetc(ch, f);
	return 0;
} /* eof */

int eoln(FILE *f)
{
	register int ch;
	
	ch = getc(f);
	if (ch == EOF)
		return 1;
	ungetc(ch, f);
	return (ch == '\n');
} /* eoln */


boolean rearrangeH(option_struct *op, data_fmt *data, long whichtchain, long chain){
	long i;
	boolean accepted;
	
	curtree = temptrees[whichtchain];
	accepted = FALSE;
	
	op->ctemp = op->temperature[whichtchain];
	
	if (randum() < 0.5) recarray[seq_length] = 1;
	
	for (i = 0; i < seq_length; i++){
		newgrec[i] = 0;
		newgweight[i] = 0;
	}
	
	
	//if (op->givenS != 1){} // given structure
	
	if (!sim_mode){
		// generate S 
		op->hotspot= 1;
		if (op->hotspot == 1){
			if (marker1pos == -99 && marker2pos == -99){
				if (op->ctemp == 99) randomhotspotputter();
				else{
					if (op->ctemp == 98) stronghotspotputter(curtree);
					else {
						if (op->ctemp == 101) {
							generaterecstructure_with_big_lamda(curtree, 1);
						}
						else generaterecstructure_with_temp(curtree,op->ctemp);
					}
				}
			}
		}
		else{
			if (op->hotspot == 0){
				for (i = 0; i < seq_length; i++) recarray[i] = recrates[0];
			}
		}
		
	}
	//generate G
	
	accepted = makedropH(op,data);
	if (accepted) Acc++;
	else Rej++;
	
	slid++;
	if (accepted && whichtchain == 0) slacc++;
	
	temptrees[whichtchain] = curtree;
	
	return accepted;
	
	
	
}




struct path *createnewpath(struct path *curpath){
	struct path *last;
	struct path *new;
	long start,length;
	int laststate,i;
	
	while (1){
		
		newpath = copypath(curpath);
		
		length = endposition - startposition;
		start = startposition;
		// find last state
		if (start != 0){
			laststate = findstate(newpath,start - 1);
			last = findblock(newpath,start - 1);
		}
		else{
			laststate = 0; //
			last = newpath;
			last->beg = 0;
			last->state = laststate;
			start = 1;
		}
		
		for (i = startposition ; i  <= endposition;i++){
			if (laststate == 0){
				if (lamda[0]<randum()){
					last->end = i;
				}
				else{
					last->end =  i - 1;
					new = (struct path *)malloc(sizeof(struct path));
					new->state = 1;
					new->beg = i;
					new->end = new->beg;
					last->next = new;
					last = new;
					laststate = 1;
					
				}
			}
			else{
				if (lamda[1]<randum()){
					last->end = i;
				}
				else{
					last->end =  i - 1;
					new = (struct path *)malloc(sizeof(struct path));
					new->state = 0;
					new->beg =  i;
					new->end = new->beg;
					last->next = new;
					last = new;
					laststate = 0;
				}
			}
		}
		
		//link to origianl path
		if (curpath->end != last->end){
			new = findblock(curpath, last->end);
			
			if (last->state == new->state){ // same state - merge them
				last->end = new->end;
				last->next = new->next;
				return newpath;
			}
			else free(newpath);
		}
		else{
			last->next = NULL;
			return newpath;
		}
	}
}

struct path *copypath(struct path *old){
	struct path *cur, *new, *temp = NULL;
	
	cur = old;
	newpath = (struct path *)malloc(sizeof(struct path));
	newpath->end = cur->end;
	newpath->state = cur->state;
	newpath->beg = cur->beg;
	
	if (cur->next == NULL) {
		newpath->next = NULL;
		return newpath;
	}
	new = (struct path *)malloc(sizeof(struct path));
	newpath->next = new;
	cur = cur->next;
	
	while (1){
		new->end = cur->end;
		new->state = cur->state;
		new->beg = cur->beg;
		if (cur->next == NULL) {
			new->next = NULL;
			return newpath;
		}
		temp = (struct path *)malloc(sizeof(struct path));
		new->next = temp;
		new = temp;
		cur = cur->next;
	}
}

struct path *findblock(struct path *paths, long pos){
	struct path *cur;
	cur = paths;
	while (1){
		if (cur->beg <= pos && cur->end >=pos) return cur;
		cur = cur->next;
	}
}

void sumstruct(struct path *pathstruct, long chain){
	int key;
	struct path *cur;
	cur = pathstruct;
	tempstructs[chain] = 0;
	while (1){
		tempstructs[chain][cur->state] += cur->end - cur->beg;
		key = cur->state;
		cur = cur->next;
		if (cur == NULL) return;
		tempstructs[chain][key + 2]++;
	}
}

int findstate(struct path *curpath, long site){
	struct path *cur;
	cur = curpath;
	while (1){
		if (cur->beg <= site && cur->end >=site) return cur->state;
		cur = cur->next;
	}
}


void freesequence(state_sequence *target, int key){
	
	if (target->prev != NULL){
		if (key == 0 && target->prev->numlink != 2) freesequence(target->prev, 0);
		else{
			if (key == 1) freesequence(target->prev, 1);
		}
	}
	free(target);
	return;
}


//initialize seglist - only for the nonrecombinant tree
void constructseglist(option_struct *op, data_fmt *data, tree *tr){
	tlist *cur;
	seglist *segments;
	int tips;
	
	tips = getdata_numtips(op,data);
	for (cur = tr->tymelist; cur != NULL; cur = cur->succ){
		segments = (struct seglist *)malloc(sizeof(struct seglist));
		segments->start = 0;
		segments->end  = countsites(op,data) - 1;
		segments->prev = segments->next = NULL;
		segments->numsam = tips;
		tips--;
		cur->segments = segments;
	}
	return;
}




void calcinitialtree(option_struct *op, data_fmt *data, tree *inittree, long *ranges){
	long i,j, start, end, numsites;
	double  *new_dlikelihood;
	seglist *cursegment;
	treeshape *subtree;
	tlist *curtlist;
	double sumdepth = 0;
	
	curtlist = inittree->tymelist;
	
	while (1){
		cursegment = curtlist->segments;
		if (curtlist->succ == NULL) break;
		curtlist = curtlist->succ;
	}
	
	cursegment = curtlist->segments;
	
	numsites = getdata_nummarkers(op,data);
	new_dlikelihood = (double *)calloc(numsites,sizeof(double));
	
	while (1){
		if (cursegment->start <= ranges[ranges[0] * 2 - 1]) break;
		cursegment = cursegment->next;
	}
	
	for (i = 0; i<seq_length; i++){
		if (dna_sequence[i] != 'S')
			sumdepth = sumdepth + 1.0;
	}
	
	
	totalbranch_s = 0;
	
	subtree = findsubs(curtree->nodep[1], 0);
	addprobs(subtree);
	
	subtree->A = subtree->C = subtree->G = subtree->T = subtree->U = 0;
	new_dlikelihood[0] = computedlikelihood(op, data, 0, subtree, 1.0/sumdepth);
	
	for (i = 1; i <= ranges[0]; i++){
		start = ranges[i * 2 - 1];
		end = ranges[i * 2];
		for (j = start; j <= end; j ++){
			if (j > cursegment->end){
				freesubtree(subtree);
				subtree = findsubs(curtree->nodep[1], j);
				subtree->A = subtree->C = subtree->G = subtree->T = subtree->U = 0;
				addprobs(subtree);
				while (1){
					if (cursegment->start <= j) break;
					cursegment = cursegment->next;
				}
			}
			
			new_dlikelihood[j] = computedlikelihood(op, data, j, subtree,1.0/sumdepth);
		}
	}
	
	freesubtree(subtree);
	inittree->dlikelihood = new_dlikelihood;
	inittree->likelihood = 0;
	for (i = 0; i < numsites; i++) inittree->likelihood += new_dlikelihood[i];
	
	//printf("initial %lf\n",inittree->likelihood);
	
	return;
}


// calculate the data likelihood of the changed regions.
void calcratio(option_struct *op, data_fmt *data, tree *oldtree, tree *newtree, long *ranges){
	long i,j, start, end, numsites;
	double *new_dlikelihood;
	treeshape *subtree = NULL;
	double tempvalue = 0;
	char *olddna = NULL;
	double sumweight = 0,siteweight;
	numsites = seq_length;
	new_dlikelihood = (double *)calloc(numsites,sizeof(double));
	
	for (i = 0; i < seq_length; i++){ 
		new_dlikelihood[i] = oldtree->dlikelihood[i];
		if (curtree->weight_array != NULL)
			sumweight = sumweight + curtree->weight_array[i];
		else
			sumweight = 1;
	}
	// ranges - rearranged region
	for (i = 1; i <= ranges[0]; i++){
		
		start = ranges[i * 2 - 1];
		end = ranges[i * 2];
		for (j = start; j <= end; j ++){
			if (subtree == NULL){
				// need to contruct sub tree
				subtree = findsubs(curtree->nodep[1], j);
				addprobs(subtree);
				// need to calculate data likelihood of non variant sites.
				subtree->A = subtree->C = subtree->G = subtree->T = subtree->U = 0;
			}
			// calculate data likelihood
			if (curtree->weight_array != NULL){
				siteweight = curtree->weight_array[j]/sumweight;
			}
			else{
				siteweight = 1;
			}
			new_dlikelihood[j] = computedlikelihood(op, data, j, subtree, siteweight);
			if (newtree->numrec_array != NULL){
				if (newtree->numrec_array[j] > 0){// recombination -> erase current subtree.
					freesubtree(subtree);
					subtree = NULL;
				}
			}
		}
	}
	
	if (olddna != NULL){
		for (i = 0; i < seq_length; i++)
			dna_sequence[i] = olddna[i];
		free(olddna);
	}
	
	if (subtree != NULL) freesubtree(subtree);
	free(newtree->dlikelihood);
	
	newtree->dlikelihood = new_dlikelihood;
	newtree->likelihood = oldtree->likelihood = 0;
	
	tempvalue = 0;
	
	//log(P(D|G)) = sum(log(P(D site| G site)))
	
	for (i = 0; i < numsites; i++){
		if (isinf(oldtree->dlikelihood[i])){
			//printf("inf %ld\n",i); // underflow
			oldtree->dlikelihood[i] = -99999;
		}
		oldtree->likelihood += oldtree->dlikelihood[i];
		newtree->likelihood += new_dlikelihood[i];
	}
	
	return;
}

void branchlength(){
	double sumbr = 0;
	tlist *curtime;
	for (curtime = curtree->tymelist; curtime->numbranch > 1; curtime = curtime->succ){
		sumbr = sumbr + curtime->numbranch * (curtime->age - curtime->eventnode->tyme);
		
		
	}
	printf("sum %lf\n",sumbr);
}

// calculate P(D site | G site)
double dlikelihood_type(option_struct *op, data_fmt *data, long position, nodeinfo *rootnode, int seqtype){
	double temp;
	dnadata *dna;
	
	dna = data->dnaptr;
	calclikelihood_type(op, data, position, rootnode, seqtype);
	
	//printf("%lf, %lf, %lf, %lf\n",rootnode->A, rootnode->C, rootnode->G, rootnode->T);
	temp =  dna->freqa * rootnode->A + dna->freqc * rootnode->C +
    dna->freqg * rootnode->G + dna->freqt * rootnode->T;
	temp = log(temp);
	return temp;
}
// calculate data likelihood of a site given subtree.
double computedlikelihood(option_struct *op, data_fmt *data, long i, treeshape *subtree,double weight){
	double newdlikelihood = 0.0;
	// variant or not?
	double nonlike,snplike;
	int j;
	
	if (filtered == 0){
		switch(dna_sequence[i]){
			case 'S':
				// variant - need to calculate
				newdlikelihood = dlikelihood(op, data, i, subtree->rootnode);
				break;
			case 'A':
				// non variant : if it is not calculated, do it. or return pre-calculated one.
				if (subtree->A == 0 ) subtree->A = dlikelihood(op, data, i, subtree->rootnode);
				newdlikelihood = subtree->A;
				break;
			case 'C':
				if (subtree->C == 0 ) subtree->C = dlikelihood(op, data, i, subtree->rootnode);
				newdlikelihood = subtree->C;
				break;
			case 'G':
				if (subtree->G == 0 ) subtree->G = dlikelihood(op, data, i, subtree->rootnode);
				newdlikelihood = subtree->G;
				break;
			case 'T':
				if (subtree->T == 0 ) subtree->T = dlikelihood(op, data, i, subtree->rootnode);
				newdlikelihood = subtree->T;
				break;
		}
	}
	else{
		if(dna_sequence[i] == 'S')
			newdlikelihood = dlikelihood(op, data, i, subtree->rootnode);
		else{
			//alpha =  expsnpnum / (double)(seq_length - num_markers) * weight;
			//alpha = 1.0 - num_markers/(expsnpnum+num_markers+0.0);
			//hiddenp = expsnpnum/(seq_length - num_markers + 0.0); // probability of hidden SNP
			nonlike = subtree->probs[0];
			snplike = 0;
			
			for (j = 1; j < old_hap; j++){
				snplike = snplike +  calling[j] * subtree->probs[j];
			}
			
			newdlikelihood = nonlike + snplike;
			newdlikelihood = log(newdlikelihood);
			//printf("old %lf, adjusted %lf %lf\n",log(nonlike), newdlikelihood,snplike);
		}
	}
	return newdlikelihood;
}

// find sub G - start with a tip and go down to the local root.
treeshape *findsubs(node *tipnode, long start){
	nodeinfo *curinfo, *leftinfo;
	node *curnode, *leftnode, *rightnode, *oldnode;
	treeshape *result;

	
	curnode = tipnode->back;
	leftnode = oldnode = tipnode;
	leftinfo = (nodeinfo *)calloc(1,sizeof(nodeinfo));
	leftinfo->nodeid = tipnode->id;
	leftinfo->nodenum = tipnode->number;
	leftinfo->istip = 1;
	
	
	
	
	while (1){
		if (iscoal(curnode)){
			// construct node just for coalescent node NO node for recombinant node.
			curnode = findunique(curnode);
			// is this node for coalesent of target site?
			if (inrange(curnode->next->back->coal, start) && inrange(curnode->next->next->back->coal,start)){
				// if so, make new node and "right" node
				if (curnode->next->back == oldnode) rightnode = curnode->next->next->back;
				else rightnode = curnode->next->back;
				
				curnode = findunique(curnode);
				curinfo = (nodeinfo *)calloc(1, sizeof(nodeinfo));
				
				// link current node to old "left" node
				curinfo->leftnode = leftinfo;
				curinfo->rightnode = NULL;
				curinfo->leftlength = curnode->tyme - leftnode->tyme;
				curinfo->rightlength = 0;
				curinfo->nodeid = curnode->id;
				curinfo->nodenum = curnode->number;
				curinfo->istip = 0;
				
				leftinfo->back = curinfo;
				
				// contruct a sub sub tree from "right" node
				climbright(curnode, oldnode, start, curinfo);
				
				// found local fc node or root node - end of the sub tree construction
				if (!inrange(curnode->coal,start) || curnode->back == curtree->root){
					result = (treeshape *)calloc(1, sizeof(treeshape));
					result->rootnode = curinfo; 
					result->probs = NULL;
					
					
					return result;
				}
				// if not, go down one more time.
				// current node became new "left" node
				leftinfo = curinfo;
				leftnode = curnode;
			}
			// we need this oldnode to find out which node is "right" node.
			oldnode = curnode;
			// go down
			curnode = curnode->back;
		}
		// recombinant node. - which recombinant node has target site?
		else{
			if (isrecomb(curnode)){
				curnode = findunique(curnode);
				if (curnode->next->recstart == 0){
					if (start <= curnode->next->recend) {
						// site on curnode->next->back
						oldnode = curnode->next;
						curnode = curnode->next->back;
					}
					else {
						// site on curnode->next->next->back
						oldnode = curnode->next->next;
						curnode = curnode->next->next->back;
					}
				}
				else{
					if (start <= curnode->next->next->recend) {
						// site on curnode->next->next->back
						oldnode = curnode->next->next;
						curnode = curnode->next->next->back;
					}
					else {
						// site on curnode->next->back
						oldnode = curnode->next;
						curnode = curnode->next->back;
					}
				}
			}
		}
	}
}

// construct sub sub tree. - contruct tree above local rootnode.
void climbright(node *localroot, node *leftnode, long start, nodeinfo *rootinfo){
	node *curnode, *root;
	nodeinfo *curinfo;

	
	// L   R <- let's construct this side
	//  \ /
	//  root
	// find "right" node
	root = findunique(localroot);
	if (findunique(root->next->back) == findunique(leftnode)) curnode = root->next->next->back;
	else curnode = root->next->back;
	
	while (1){
		// local root node should be the coalescent node
		if (isrecomb(curnode)){
			curnode = findunique(curnode);
			curnode = curnode->back;
		}
		else{
			curnode = findunique(curnode);
			if (istip(curnode)) break;
			// local root node should be the coalescent node of the coalescent of target site.
			if (inrange(curnode->next->back->coal, start) && inrange(curnode->next->next->back->coal, start)){
				break;
			}
			else{
				if(inrange(curnode->next->back->coal,start)){
					curnode = curnode->next->back;
				}
				else{
					curnode = curnode->next->next->back;
				}
			}
		}
	}
	
	// OK, we found "right" node
	curinfo = (nodeinfo *)calloc(1,sizeof(nodeinfo));
	rootinfo->rightnode = curinfo;
	rootinfo->rightlength = localroot->tyme - curnode->tyme;
	curinfo->nodeid = curnode->id;
	curinfo->nodenum = curnode->number;
	curinfo->back = rootinfo;
	
	if(istip(curnode)){
		curinfo->rightnode = curinfo->leftnode = NULL;
		curinfo->istip = 1;
	}
	// let's construct the sub sub tree which has the "right" node as a root.
	else climbup(curnode, start, curinfo);
	
	return;
}

// construct sub sub tree above local root.
void climbup(node *localroot, long start, nodeinfo *curinfo){ // curnode : always coal node
	node *node1, *node2, *curnode;
	nodeinfo *node1info, *node2info;
	
	curnode = findunique(localroot);
	node1 = curnode->next->back;
	node2 = curnode->next->next->back;
	
	
	while (1){
		if (istip(node1)) break;
		if (isrecomb(node1)) node1 = findunique(node1)->back;
		else{
			if (!(inrange(node1->next->back->coal, start)) || !(inrange(node1->next->next->back->coal,start))){
				if (inrange(node1->next->back->coal,start)) node1 = node1->next->back;
				else node1 = node1->next->next->back;
			}
			else break;
		}
	}
	
	node1info = (nodeinfo *)calloc(1,sizeof(nodeinfo));
	node1info->back = curinfo;
	node1info->nodeid = node1->id;
	node1info->nodenum = node1->number;
	curinfo->rightlength = curnode->tyme - node1->tyme;
	curinfo->rightnode = node1info;
	
	if(istip(node1)){
		node1info->istip = 1;
	}
	else{
		climbup(node1, start, node1info);
	}
	
	while (1){
		if (istip(node2)) break;
		if (isrecomb(node2)) node2 = findunique(node2)->back;
		else{
			if (!inrange(node2->next->back->coal,start) || !inrange(node2->next->next->back->coal,start)){
				if (inrange(node2->next->back->coal,start)) node2 = node2->next->back;
				else node2 = node2->next->next->back;
			}
			else break;
		}
	}
	
	node2info = (nodeinfo *)calloc(1,sizeof(nodeinfo));
	node2info->back = curinfo;
	node2info->nodeid = node2->id;
	node2info->nodenum = node2->number;
	curinfo->leftlength = curnode->tyme - node2->tyme;
	curinfo->leftnode = node2info;
	
	if(istip(node2)){
		node2info->istip = 1;
	}
	else{
		climbup(node2, start, node2info);
	}
	
	return;
}


void nuview_single(option_struct *op, data_fmt *data, long position, double length, nodeinfo *branchnode, nodeinfo *localroot){
	
	int i;
	double lw1 = 0.0;
	boolean toobig;
	
	
	lw1 = -1.0 * length / data->dnaptr->fracchange;
	toobig = (exp(lw1) <= 0.0);
	
	if (toobig) {
		for (i = 0; i < op->categs; i++) {
			tbl[i].ww1 = tbl[i].zz1 = 0.0;
			tbl[i].ww1zz1 = tbl[i].vv1zz1 = 0.0;
		}
	}
	else {
		for (i = 0; i < op->categs; i++) {
			tbl[i].ww1 = exp(tbl[i].rat_xi * lw1);
			tbl[i].zz1 = exp(tbl[i].rat_xv * lw1);
			tbl[i].ww1zz1 = tbl[i].ww1 * tbl[i].zz1;
			tbl[i].vv1zz1 = (1.0 - tbl[i].ww1) * tbl[i].zz1;
		}
	}
	
	
	// calculrate dlikelihood - right now only one rate.(using tbl[0])
	calcrange_single(op,data,branchnode, localroot);
	return;
}


void calcrange_single(option_struct *op, data_fmt *data, nodeinfo *branchnode, nodeinfo *localroot){
	
	double yy1, ww1zz1, vv1zz1, vv1zz1_sumr1, vv1zz1_sumy1, sum1, sumr1, sumy1;
	dnadata *dna;
	
	
	dna = data->dnaptr;
	
	ww1zz1 = tbl[0].ww1zz1;
	vv1zz1 = tbl[0].vv1zz1;
	yy1 = 1.0 - tbl[0].zz1;
	sum1 = yy1 * (dna->freqa * branchnode->A +
				  dna->freqc *  branchnode->C +
				  dna->freqg *  branchnode->G +
				  dna->freqt *  branchnode->T);
	sumr1 = dna->freqar *  branchnode->A + dna->freqgr *  branchnode->G;
	sumy1 = dna->freqcy *  branchnode->C + dna->freqty *  branchnode->T;
	
	vv1zz1_sumr1 = vv1zz1 * sumr1;
	vv1zz1_sumy1 = vv1zz1 * sumy1;
	
	localroot->A =  (sum1 + ww1zz1 * branchnode->A + vv1zz1_sumr1);
	localroot->C =  (sum1 + ww1zz1 *  branchnode->C + vv1zz1_sumy1);
	localroot->G =  (sum1 + ww1zz1 *  branchnode->G + vv1zz1_sumr1);
	localroot->T =  (sum1 + ww1zz1 *  branchnode->T + vv1zz1_sumy1);
	
}

void calclikelihood_type(option_struct *op, data_fmt *data, long position, nodeinfo *localroot, int seqtype){
	if (localroot->leftnode != NULL){
		calclikelihood_type(op,data,position, localroot->leftnode,seqtype);
	}
	if (localroot->rightnode != NULL){
		calclikelihood_type(op,data,position, localroot->rightnode,seqtype);
	}
	nuview_modified_type(op,data,position,localroot,seqtype);
}


void nuview_modified_type(option_struct *op, data_fmt *data, long position, nodeinfo *localroot, int seqtype){
	
	int i;
	double lw1 = 0.0, lw2 = 0.0;
	boolean toobig;
	
	// is it tip node?
	if (localroot->istip == 1){
		checkbase_type(localroot, seqtype);
		return;
	}
	
	// from "left" branch
	lw1 = -1.0 * localroot->leftlength / data->dnaptr->fracchange;
	toobig = (exp(lw1) <= 0.0);
	
	if (toobig) {
		for (i = 0; i < op->categs; i++) {
			tbl[i].ww1 = tbl[i].zz1 = 0.0;
			tbl[i].ww1zz1 = tbl[i].vv1zz1 = 0.0;
		}
	}
	else {
		for (i = 0; i < op->categs; i++) {
			tbl[i].ww1 = exp(tbl[i].rat_xi * lw1);
			tbl[i].zz1 = exp(tbl[i].rat_xv * lw1);
			tbl[i].ww1zz1 = tbl[i].ww1 * tbl[i].zz1;
			tbl[i].vv1zz1 = (1.0 - tbl[i].ww1) * tbl[i].zz1;
		}
	}
	
	//  from "right" branch
	lw2 = -1.0 * localroot->rightlength / data->dnaptr->fracchange;
	toobig = (exp(lw2) <= 0.0);
	
	if (toobig) {
		for (i = 0; i < op->categs; i++) {
			tbl[i].ww2 = tbl[i].zz2 = 0.0;
			tbl[i].ww2zz2 = tbl[i].vv2zz2 = 0.0;
		}
	}
	else {
		for (i = 0; i < op->categs; i++) {
			tbl[i].ww2 = exp(tbl[i].rat_xi * lw2);
			tbl[i].zz2 = exp(tbl[i].rat_xv * lw2);
			tbl[i].ww2zz2 = tbl[i].ww2 * tbl[i].zz2;
			tbl[i].vv2zz2 = (1.0 - tbl[i].ww2) * tbl[i].zz2;
		}
	}
	
	// calculrate dlikelihood - right now only one rate.(using tbl[0])
	calcrange_modified(op,data,localroot);
	return;
}


void calclikelihood(option_struct *op, data_fmt *data, long position, nodeinfo *localroot){
	if (localroot->leftnode != NULL){
		calclikelihood(op,data,position, localroot->leftnode);
	}
	if (localroot->rightnode != NULL){
		calclikelihood(op,data,position, localroot->rightnode);
	}
	nuview_modified(op,data,position,localroot);
}


void nuview_modified(option_struct *op, data_fmt *data, long position, nodeinfo *localroot){
	
	int i;
	double lw1 = 0.0, lw2 = 0.0;
	boolean toobig;
	
	// is it tip node?
	if (localroot->istip == 1){
		checkbase(op, data, position, localroot);
		return;
	}
	
	// from "left" branch
	lw1 = -1.0 * localroot->leftlength / data->dnaptr->fracchange;
	toobig = (exp(lw1) <= 0.0);
	
	if (toobig) {
		for (i = 0; i < op->categs; i++) {
			tbl[i].ww1 = tbl[i].zz1 = 0.0;
			tbl[i].ww1zz1 = tbl[i].vv1zz1 = 0.0;
		}
	}
	else {
		for (i = 0; i < op->categs; i++) {
			tbl[i].ww1 = exp(tbl[i].rat_xi * lw1);
			tbl[i].zz1 = exp(tbl[i].rat_xv * lw1);
			tbl[i].ww1zz1 = tbl[i].ww1 * tbl[i].zz1;
			tbl[i].vv1zz1 = (1.0 - tbl[i].ww1) * tbl[i].zz1;
		}
	}
	
	//  from "right" branch
	lw2 = -1.0 * localroot->rightlength / data->dnaptr->fracchange;
	toobig = (exp(lw2) <= 0.0);
	
	if (toobig) {
		for (i = 0; i < op->categs; i++) {
			tbl[i].ww2 = tbl[i].zz2 = 0.0;
			tbl[i].ww2zz2 = tbl[i].vv2zz2 = 0.0;
		}
	}
	else {
		for (i = 0; i < op->categs; i++) {
			tbl[i].ww2 = exp(tbl[i].rat_xi * lw2);
			tbl[i].zz2 = exp(tbl[i].rat_xv * lw2);
			tbl[i].ww2zz2 = tbl[i].ww2 * tbl[i].zz2;
			tbl[i].vv2zz2 = (1.0 - tbl[i].ww2) * tbl[i].zz2;
		}
	}
	
	// calculrate dlikelihood - right now only one rate.(using tbl[0])
	calcrange_modified(op,data,localroot);
	return;
}

void calcrange_modified(option_struct *op, data_fmt *data, nodeinfo *localroot){
	
	double yy1, yy2, ww1zz1, vv1zz1, ww2zz2, vv2zz2, vv1zz1_sumr1,
    vv2zz2_sumr2, vv1zz1_sumy1, vv2zz2_sumy2, sum1, sum2,
    sumr1, sumr2, sumy1, sumy2;
	dnadata *dna;
	
	
	dna = data->dnaptr;
	
	ww1zz1 = tbl[0].ww1zz1;
	vv1zz1 = tbl[0].vv1zz1;
	yy1 = 1.0 - tbl[0].zz1;
	sum1 = yy1 * (dna->freqa * localroot->leftnode->A +
				  dna->freqc * localroot->leftnode->C +
				  dna->freqg * localroot->leftnode->G +
				  dna->freqt * localroot->leftnode->T);
	sumr1 = dna->freqar * localroot->leftnode->A + dna->freqgr * localroot->leftnode->G;
	sumy1 = dna->freqcy * localroot->leftnode->C + dna->freqty * localroot->leftnode->T;
	
	vv1zz1_sumr1 = vv1zz1 * sumr1;
	vv1zz1_sumy1 = vv1zz1 * sumy1;
	
	localroot->A =  (sum1 + ww1zz1 * localroot->leftnode->A + vv1zz1_sumr1);
	localroot->C =  (sum1 + ww1zz1 * localroot->leftnode->C + vv1zz1_sumy1);
	localroot->G =  (sum1 + ww1zz1 * localroot->leftnode->G + vv1zz1_sumr1);
	localroot->T =  (sum1 + ww1zz1 * localroot->leftnode->T + vv1zz1_sumy1);
	
	ww2zz2 = tbl[0].ww2zz2;
	vv2zz2 = tbl[0].vv2zz2;
	yy2 = 1.0 - tbl[0].zz2;
	
	sum2 = yy2 * (dna->freqa * localroot->rightnode->A +
				  dna->freqc * localroot->rightnode->C +
				  dna->freqg * localroot->rightnode->G +
				  dna->freqt * localroot->rightnode->T);
	sumr2 = dna->freqar * localroot->rightnode->A + dna->freqgr * localroot->rightnode->G;
	sumy2 = dna->freqcy * localroot->rightnode->C + dna->freqty * localroot->rightnode->T;
	
	vv2zz2_sumr2 = vv2zz2 * sumr2;
	vv2zz2_sumy2 = vv2zz2 * sumy2;
	localroot->A *=  (sum2 + ww2zz2 * localroot->rightnode->A + vv2zz2_sumr2);
	localroot->C *=  (sum2 + ww2zz2 * localroot->rightnode->C + vv2zz2_sumy2);
	localroot->G *=  (sum2 + ww2zz2 * localroot->rightnode->G + vv2zz2_sumr2);
	localroot->T *=  (sum2 + ww2zz2 * localroot->rightnode->T + vv2zz2_sumy2);
	
}


// calculate P(D site | G site)
double dlikelihood(option_struct *op, data_fmt *data, long position, nodeinfo *rootnode){
	double temp;
	dnadata *dna;
	
	dna = data->dnaptr;
	calclikelihood(op, data, position, rootnode);
	
	//printf("%lf, %lf, %lf, %lf\n",rootnode->A, rootnode->C, rootnode->G, rootnode->T);
	temp =  dna->freqa * rootnode->A + dna->freqc * rootnode->C +
    dna->freqg * rootnode->G + dna->freqt * rootnode->T;
	temp = log(temp);
	return temp;
}




void empiricaldnafreqs_modified(option_struct *op,data_fmt *data){
	long i, k;
	long A = 0, C = 0, G = 0, T = 0;
	double temp, suma, sumc, sumg, sumt, w;
	dnadata *dna;
	
	dna = data->dnaptr;
	
	dna->freqa = 0.25;
	dna->freqc = 0.25;
	dna->freqg = 0.25;
	dna->freqt = 0.25;
	suma = 0.0;
	sumc = 0.0;
	sumg = 0.0;
	sumt = 0.0;
	
	
	for (k = 0; k < getdata_nummarkers(op,data); k++) {
		w = dna->dnaweight[k];
		for (i = 1; i <= getdata_numseq(op,data); i++) {
			A = C = G = T = 0;
			strcpy(curtree->nodep[i]->nayme,data->dnaptr->indnames[population][locus][i-1]);
			
			switch (dna->seqs[population][locus][i-1][k]) {
					
				case 'A':
					A = 1;
					break;
				case 'C':
					C = 1;
					break;
				case 'G':
					G = 1;
					break;
				case 'T':
					T = 1;
					break;
				case 'U':
					T = 1;
					break;
				case 'M':
					A = C = 1;
					break;
				case 'R':
					A = G = 1;
					break;
				case 'W':
					A = T = 1;
					break;
				case 'S':
					G = C = 1;
					break;
				case 'Y':
					C = T = 1;
					break;
				case 'K':
					G = T = 1;
					break;
				case 'B':
					C = G = T = 1;
					break;
				case 'D':
					A = G = T = 1;
					break;
				case 'H':
					A = C = T = 1;
					break;
				case 'V':
					A = C = G = 1;
					break;
				case 'N':
					A = C = G = T = 1;
					break;
				case 'X':
					A = C = G = T = 1;
					break;
				case '?':
					A = C = G = T = 1;
					break;
				case 'O':
					A = C = G = T = 1;
					break;
				case '-':
					A = C = G = T = 1;
					break;
				default:
					fprintf(ERRFILE,"**ERROR in setting up tips, probable error");
					fprintf(ERRFILE," in input data file.\nInterleaved/sequential");
					fprintf(ERRFILE," status may be incorrectly set.\n\n");
					exit(-1);
					break;
			}
			
			temp = dna->freqa * A;
			temp += dna->freqc * C;
			temp += dna->freqg * G;
			temp += dna->freqt * T;
			suma += w * dna->freqa * A / temp;
			sumc += w * dna->freqc * C /  temp;
			sumg += w * dna->freqg * G / temp;
			sumt += w * dna->freqt * T / temp;
		}
	}
	
	temp = suma + sumc + sumg + sumt;
	dna->freqa = suma / temp;
	dna->freqc = sumc / temp;
	dna->freqg = sumg / temp;
	dna->freqt = sumt / temp;
	return;
}

void constructseqinfo(option_struct *op, data_fmt *data, double freq){
	long x;
	
	num_markers = 0;
	for (x = 0; x < seq_length; x++) {
		dna_sequence[x] = checkSNP(op,data,x, freq);
		if (dna_sequence[x] == 'S')num_markers ++;
	}
	/*   numspaces = num_markers + 1; */
	/*   if (dna_sequence[0] == 'S') numspaces--; */
	/*   if (dna_sequence[numsites - 1] == 'S') numspaces--; */
	
	/*   y = 1; */
	
}

char checkSNP_old(option_struct *op, data_fmt *data, long position){ // check SNP
	long i;
	char ****y;
	char base1 = 'X', base2 = 'X';
	
	y = data->dnaptr->seqs;
	for (i = 1; i <= getdata_numseq(op,data); i++) {
		base2 = y[0][0][i-1][position];
		if (base1 == 'X') base1 = base2;
		if (base1 != base2) return 'S';
	}
	return base1;
}

void checkbase(option_struct *op, data_fmt *data, long position, nodeinfo *tipnode){
	int temp;
	char ****y;
	
	//tipnode->A = tipnode->C = tipnode->G = tipnode->T = 1;
	//return;
	
	y = data->dnaptr->seqs;
	
	tipnode->A = tipnode->C = tipnode->G = tipnode->T = 0;
	
	if (tipnode->nodenum == notip){ // ignoring this tip
		tipnode->A = tipnode->C = tipnode->G = tipnode->T = 1;
		return;
	}
	if (sim_mode && tipnode->nodenum > old_hap){ // newtips
		temp = newhaps[tipnode->nodenum - old_hap];
		switch (temp){ 
			case 1:
				tipnode->A = 1;
				break;
			case 2:
				tipnode->C = 1;
				break;
			case 3:
				tipnode->G = 1;
				break;
			case 4:
				tipnode->T = 1;
				break;
			case 0:
				tipnode->A = tipnode->C = tipnode->G  = tipnode->T = 1;
				break;
		}
		return;
	}
	
	if (sim_mode && tipnode->nodenum <= old_hap){ // oldtips (simulated)
		switch (olddata[tipnode->nodenum - 1][position]){ 
			case 1:
				tipnode->A = 1;
				break;
			case 2:
				tipnode->C = 1;
				break;
			case 3:
				tipnode->G = 1;
				break;
			case 4:
				tipnode->T = 1;
				break;
			case 0:
				tipnode->A = tipnode->C = tipnode->G = tipnode->T = 1;
				break;
				
		}
		return;
	}
	
	
	switch (y[0][0][tipnode->nodenum - 1][position]){
		case 'A':
			tipnode->A = 1;
			break;
		case 'C':
			tipnode->C = 1;
			break;
		case 'G':
			tipnode->G = 1;
			break;
		case 'T':
			tipnode->T = 1;
			break;
		case 'N':
			tipnode->A = tipnode->C = tipnode->G = tipnode->T = 1;
			break;
			
	}
	return;
}


void checkbase_type(nodeinfo *tipnode, int seqtype){


	
	tipnode->A = tipnode->C = tipnode->G = tipnode->T = 0;
	
	switch (seqtype){
		case 1:
			tipnode->A = 1;
			break;
		case 2:
			tipnode->C = 1;
			break;
		case 3:
			tipnode->G = 1;
			break;
		case 4:
			tipnode->T = 1;
			break;
		case 0:
			tipnode->A = tipnode->C = tipnode->G  = tipnode->T = 1;
			break;
	}
	return;
}


void freesubtree(treeshape *target){
	freenodeinfo(target->rootnode);
	if(target->probs != NULL)
		free(target->probs);
	free(target);
	target = NULL;
}

void freenodeinfo(nodeinfo *rootnode){
	if (rootnode->rightnode != NULL) freenodeinfo(rootnode->rightnode);
	if (rootnode->leftnode !=NULL) freenodeinfo(rootnode->leftnode);
	
	free(rootnode);
}


void randompathsampling_temp(double *oldrecarray, double *newrecarray, double *weight, int *recnums, long start, long end){
	double prob_cold, prob_hot, cold_ini, hot_ini, *backward_cold, *backward_hot;
	long i;
	char lastone, endstate = 'C';
	
	memcpy(newrecarray,oldrecarray,seq_length*sizeof(double));
	
	backward_cold = (double *)calloc(seq_length,sizeof(double));
	backward_hot = (double *)calloc(seq_length,sizeof(double));
	
	if (end != seq_length - 1){
		if (oldrecarray[end + 1] == recrates[0]) endstate = 'C';
		else endstate = 'H';
	}
	
	backward_temp(backward_cold, backward_hot, weight, recnums, start, end, endstate);
	
	if (start == 0){
		cold_ini = 1 / (1 - lamda[0]) /( 1 / (1 - lamda[0]) + 1 / (1 - lamda[1]));
		hot_ini = 1 - cold_ini;
		
		prob_cold = cold_ini * exp( - recrates[0] * weight[0]) * pow(recrates[0],recnums[0]) * backward_cold[0];
		prob_hot = hot_ini * exp( - recrates[1] * weight[0]) * pow(recrates[1],recnums[0]) * backward_hot[0];
		
		if (randum() < (prob_cold/(prob_cold + prob_hot))){
			lastone = 'C';
			newrecarray[0] = recrates[0];
		}
		else{
			lastone = 'H';
			newrecarray[0] = recrates[1];
		}
	}
	else{
		if (oldrecarray[start - 1] == recrates[0]) lastone = 'C';
		else lastone = 'H';
	}
	
	
	for (i = start; i <= end ; i++){
		if (lastone == 'C'){
			prob_cold = lamda[0] * exp( - recrates[0] * weight[i]) * pow(recrates[0], recnums[i]) * backward_cold[i];//checked
			prob_hot = (1 - lamda[0]) * exp( - recrates[1] * weight[i]) * pow(recrates[1], recnums[i]) * backward_hot[i];//checked
		}
		else{
			prob_cold = (1 - lamda[1]) * exp( - recrates[0] * weight[i]) * pow(recrates[0],recnums[i]) * backward_cold[i];//checked
			prob_hot = lamda[1] * exp( - recrates[1] * weight[i]) * pow(recrates[1], recnums[i]) * backward_hot[i];//checked
		}
		if (randum() < (prob_cold/(prob_cold + prob_hot))){
			lastone = 'C';
			newrecarray[i] = recrates[0];
		}
		else{
			lastone = 'H';
			newrecarray[i] = recrates[1];
		}
	}
	free(backward_hot);
	free(backward_cold);
	return;
}


void backward_temp(double *cold, double *hot, double *weight, int *recnum, long start, long end, char lastone){
	long i;
	double cold_temp, hot_temp;
	
	if (end == seq_length - 1){
		//cold[end] = 1 / (1 - lamda[0]) /( 1 / (1 - lamda[0]) + 1 / (1 - lamda[1]));
		//hot[end] = 1 - cold[end];
		cold[end] = 1;
		hot[end] = 1;
	}
	else{
		if (lastone == 'C'){
			cold_temp = lamda[0] * exp(- recrates[0] * weight[end]) * pow(recrates[0],recnum[end]);
			hot_temp = (1 - lamda[1]) * exp(- recrates[1] * weight[end]) * pow(recrates[1], recnum[end]);
		}
		else{
			cold_temp = (1- lamda[0]) * exp(- recrates[0] * weight[end]) * pow(recrates[0],recnum[end]);
			hot_temp = lamda[1] * exp(- recrates[1] * weight[end]) * pow(recrates[1],recnum[end]);
		}
		cold[end] = cold_temp/(cold_temp + hot_temp);//checked
		hot[end] = hot_temp/(cold_temp + hot_temp);//checked
	}
	
	for(i = end - 1; i >= start; i--){
		cold_temp = lamda[0] * exp(- recrates[0] * weight[i]) * pow(recrates[0],recnum[i]) * cold[i + 1]
		+ (1 - lamda[0]) * exp(- recrates[1] * weight[i]) * pow(recrates[1],recnum[i]) * hot[i + 1]; //checked
		hot_temp = (1 - lamda[1]) * exp(- recrates[0] * weight[i]) * pow(recrates[0],recnum[i]) * hot[i + 1]
		+ lamda[1] * exp(- recrates[1] * weight[i]) * pow(recrates[1], recnum[i]) * hot[i + 1]; //checked
		cold[i] = cold_temp/(cold_temp + hot_temp);//checked
		hot[i] = hot_temp/(cold_temp + hot_temp);//checked
	}
	return;
}


// generate recombination array (numrec_array) and site array (weight_array)
void scorerecs_temp(tree *target){
	tlist *curtyme;
	long i;
	
	if (target->weight_array == NULL) target->weight_array = (double *)calloc(seq_length,sizeof(double));
	if (target->numrec_array == NULL) target->numrec_array = (double *)calloc(seq_length,sizeof(double));
	
	
	for(i = 0; i<seq_length;i++){
		target->weight_array[i] = target->numrec_array[i] = 0;
	}
	
	curtyme = target->tymelist->succ;
	
	while(1){
		if (iscoal(curtyme->eventnode)){
			addtime(findunique(curtyme->eventnode)->next->back, target->weight_array);
			addtime(findunique(curtyme->eventnode)->next->next->back, target->weight_array);
		}
		else{
			addtime(findunique(curtyme->eventnode)->back, target->weight_array);
			if (curtyme->eventnode->recstart == 0) target->numrec_array[curtyme->eventnode->recend]++;
			else target->numrec_array[curtyme->eventnode->recstart - 1]++;
		}
		if (curtyme->succ == NULL) return;
		curtyme = curtyme->succ;
	}
}

void addtime(node *target, double *weights){
	long y;
	
	for(y = target->coal[1]; y < target->coal[target->coal[0] * 2]; y++){
		weights[y] += target->length;
	}
	return;
}



void realrandompathsampling_temp(double *oldrecarray, double *newrecarray, double *weight, int *recnums, long start, long end){
	double prob_cold, prob_hot, cold_ini, hot_ini, *backward_cold, *backward_hot;
	long i,length;
	char lastone, endstate = 'C';
	
	length = seq_length;
	
	memcpy(newrecarray,oldrecarray,length*sizeof(double));
	
	backward_cold = (double *)calloc(length,sizeof(double));
	backward_hot = (double *)calloc(length,sizeof(double));
	
	if (end != seq_length - 1){
		if (oldrecarray[end + 1] == recrates[0]) endstate = 'C';
		else endstate = 'H';
	}
	
	realbackward_temp(backward_cold, backward_hot, weight, recnums, start, end, endstate);
	
	if (start == 0){
		cold_ini = 1 / (1 - lamda[0]) /( 1 / (1 - lamda[0]) + 1 / (1 - lamda[1]));
		hot_ini = 1 - cold_ini;
		
		prob_cold = cold_ini * backward_cold[0];
		prob_hot = hot_ini  * backward_hot[0];
		
		if (randum() < (prob_cold/(prob_cold + prob_hot))){
			lastone = 'C';
			newrecarray[0] = recrates[0];
		}
		else{
			lastone = 'H';
			newrecarray[0] = recrates[1];
		}
	}
	else{
		if (oldrecarray[start - 1] == recrates[0]) lastone = 'C';
		else lastone = 'H';
	}
	
	
	for (i = start; i <= end ; i++){
		if (lastone == 'C'){
			prob_cold = lamda[0] * backward_cold[i];//checked
			prob_hot = (1 - lamda[0]) * backward_hot[i];//checked
		}
		else{
			prob_cold = (1 - lamda[1])* backward_cold[i];//checked
			prob_hot = lamda[1] * backward_hot[i];//checked
		}
		if (randum() < (prob_cold/(prob_cold + prob_hot))){
			lastone = 'C';
			newrecarray[i] = recrates[0];
		}
		else{
			lastone = 'H';
			newrecarray[i] = recrates[1];
		}
	}
	return;
}


void realbackward_temp(double *cold, double *hot, double *weight, int *recnum, long start, long end, char lastone){
	long i;
	double cold_temp, hot_temp;
	
	if (end == seq_length - 1){
		cold[end] = 1 / (1 - lamda[0]) /( 1 / (1 - lamda[0]) + 1 / (1 - lamda[1]));
		hot[end] = 1 - cold[end];
	}
	else{
		if (lastone == 'C'){
			cold_temp = lamda[0];
			hot_temp = (1 - lamda[1]);
		}
		else{
			cold_temp = (1- lamda[0]);
			hot_temp = lamda[1];
		}
		cold[end] = cold_temp/(cold_temp + hot_temp);//checked
		hot[end] = hot_temp/(cold_temp + hot_temp);//checked
	}
	
	for(i = end - 1; i >= start; i--){
		cold_temp = lamda[0] + (1 - lamda[0]);
		hot_temp = (1 - lamda[1]) + lamda[1];
		cold[i] = cold_temp/(cold_temp + hot_temp);//checked
		hot[i] = hot_temp/(cold_temp + hot_temp);//checked
	}
	return;
}

void whereisrec(int *recnum){
	long i;
	for(i = 0;i<100; i++) part1 += recnum[i];
	for(i = 100;i<200; i++) part2 += recnum[i];
	for(i = 200;i<300; i++) part3 += recnum[i];
	for(i = 300;i<400; i++) part4 += recnum[i];
	for(i = 400;i<500; i++) part5 += recnum[i];
	for(i = 500;i<600; i++) part6 += recnum[i];
	for(i = 600;i<700; i++) part7 += recnum[i];
	for(i = 700;i<800; i++) part8 += recnum[i];
	for(i = 800;i<900; i++) part9 += recnum[i];
	for(i = 950; i<1000;i++) part10 += recnum[i];
	
	return;
}
double Sprob(tree *target){
	double prob = 0,cold_ini,hot_ini;
	long i;
	
	cold_ini = 1 / (1 - lamda[0]) /( 1 / (1 - lamda[0]) + 1 / (1 - lamda[1]));
	hot_ini = 1 - cold_ini;
	
	if (target->recrates[0] == recrates[0]) prob = log(cold_ini);
	else prob = log(hot_ini);
	
	for (i = 1; i < seq_length; i++){
		if (target->recrates[i - 1] == recrates[0]){
			if (target->recrates[i] == recrates[0]) prob += log(lamda[0]);
			else prob += log(1 - lamda[0]);
		}
		else{
			if (target->recrates[i] == recrates[1]) prob += log(lamda[1]);
			else prob += log(1 - lamda[1]);
		}
	}
	
	return prob;
}

double probGS(tree *target, double *rates){
	double tk, tyme, prob;
	long i;
	
	tk = tyme = prob = 0;
	
	for (i = 0; i < seq_length; i++){
		prob+= target->numrec_array[i] * log(rates[i]) - target->weight_array[i] * rates[i];
	}
	return prob;
}



void generateuniformS(){
	double cold_ini,hot_ini,rate;
	long i;
	cold_ini = 1 / (1 - lamda[0]) /( 1 / (1 - lamda[0]) + 1 / (1 - lamda[1]));
	hot_ini = 1 - cold_ini;
	
	rate = cold_ini * recrates[0] + hot_ini * 2 * recrates[1];
	
	for(i=0;i<seq_length;i++) recarray[i] = rate;
	return;
}


boolean structureheating(tree *target, double *newS, double *oldS, double temperature){
	double x, old, new;
	
	old = probGS(target,oldS) + probS(oldS);
	new = probGS(target,newS) + probS(newS);
	
	if (new < old) return TRUE;
	x = (1.0/temperature - 1) * (new - old);
	if (x > log(randum())) return TRUE;
	else return FALSE;
}


boolean simulatedtemperingswap(double prob, double C1, double C2, double temperature1, double temperature2){
	double x;
	
	x = log(C1 - C2) + prob * (1.0/temperature1 - 1.0/temperature2);
	if (log(randum() < x)) return TRUE;
	else return FALSE;
	
}

boolean ratiowithheat(tree *oldtree, tree *newtree, double *oldS, double *newS, double temperature){
	double oldtreelikelihood, newtreelikelihood, oldGoldS, newGoldS, oldGnewS, newGnewS, poldS, pnewS, x, randomvalue;
	
	oldtreelikelihood = treelikelihood(oldtree);
	newtreelikelihood = treelikelihood(newtree);
	
	oldGoldS = probGS(oldtree, oldS);
	oldGnewS = probGS(oldtree, newS);
	newGoldS = probGS(newtree, oldS);
	newGnewS = probGS(newtree, newS);
	
	poldS = probS(oldS);
	pnewS = probS(newS);
	
	
	x = (1.0 / temperature) * (newtree->dlikelihood - oldtree->dlikelihood)
    + (1.0 / temperature - 1) * (newGnewS + pnewS - oldGoldS - poldS)
    + (newGoldS - oldGnewS) + (oldtreelikelihood - newtreelikelihood);
	
	randomvalue = randum();
	
	if (x > randomvalue) return TRUE;
	else return FALSE;
}




double probS(double *recstructure){ // return log(prob(S))
	long i;
	double prob = 0, cold_ini, hot_ini;
	
	
	cold_ini = 1 / (1 - lamda[0]) /( 1 / (1 - lamda[0]) + 1 / (1 - lamda[1]));
	hot_ini = 1 - cold_ini;
	
	if (recstructure[0] == recrates[0]) prob = log(cold_ini);
	else prob = log(hot_ini);
	
	if (seq_length == 1) return prob;
	
	for (i = 1; i<(seq_length - 1); i++){
		if (recstructure[i - 1] == recrates[0]){
			if (recstructure[i] == recrates[0]) prob += log(lamda[0]);
			else prob += log(1 - lamda[0]);
		}
		else{
			if (recstructure[i] == recrates[0]) prob += log(1 - lamda[1]);
			else prob += log(lamda[1]);
		}
	}
	
	return prob;
}

double treelikelihood(tree *target){//calculate P(G|S,theta)P(S|HMM parameters)
	return probS(target->recrates) + probGS(target,target->recrates);
}



void initrecarray(double r0){
	long i;
	
	for(i = 0; i < seq_length; i++){
		recarray[i] = r0;
	}
	return;
}


void randompathgenerator(){
	long i;
	
	
	if (randum() < (1 / (1 - lamda[0]) /( 1 / (1 - lamda[0]) + 1 / (1 - lamda[1]))))   recarray[0] = recrates[0];
	else recarray[0] = recrates[1];
	
	for (i = 1; i < seq_length; i++){
		if (recarray[i - 1] == recrates[0]){
			if (randum() < lamda[0]) recarray[i] = recrates[0];
			else recarray[i] = recrates[1];
		}
		else{
			if (randum() < lamda[1]) recarray[i] = recrates[1];
			else recarray[i] = recrates[0];
		}
	}
}

double forward(tree *target){
	long i;
	double n_c = 0, n_h = 0, n1_c = 1,n1_h = 1;
	double minvalue;
	
	for (i = 0; i<seq_length; i++){
		if (n1_c < n1_h) minvalue = n1_c;
		else minvalue = n1_h;
		
		n1_c = exp(n1_c - minvalue);
		n1_h = exp(n1_h - minvalue);
		
		n_c = (n1_c * lamda[0] + n1_h * (1.0 - lamda[1])) * pow(recrates[0],target->numrec_array[i]) * exp(target->weight_array[i] * recrates[0]);
		n_h = (n1_c * (1 - lamda[0]) + n1_h * lamda[1]) * pow(recrates[1],target->numrec_array[i]) * exp(target->weight_array[i] * recrates[1]);
		
		n1_c = log(n_c) + minvalue;
		n1_h = log(n_h) + minvalue;
	}
	
	if (n_c < n_h) minvalue = n_c;
	else minvalue = n_h;
	
	n_c = exp(n_c - minvalue);
	n_h = exp(n_h - minvalue);
	
	return log(n_c + n_h) + minvalue;
	
	
}

double gethmmlikelihoodtwotree(treerec *tr1, treerec *tr2){
	double alpha0 = 0.0;
	double alpha1 = 0.0;
	double temp0,temp1;
	long i;
	double ini[2];
	
	ini[0] = (1 - lamda[1]) /  (1 - lamda[0] + 1 - lamda[1]);
	ini[1] = 1 - ini[0];
	
	
	alpha0 = ini[0] * pow(recrates[0],tr1->rec[0]) * exp(-tr1->sumsite[0] * recrates[0]) * pow(recrates[0],tr2->rec[0]) * exp(-tr2->sumsite[0] * recrates[0]);
	alpha1 = ini[1] * pow(recrates[1],tr1->rec[0]) * exp(-tr1->sumsite[0] * recrates[1]) * pow(recrates[1],tr2->rec[0]) * exp(-tr2->sumsite[0] * recrates[1]);
	
	for (i = 1; i<(seq_length - 1); i++){
		temp0 = (alpha0 * lamda[0] + alpha1 * (1 - lamda[1])) * pow(recrates[0],tr1->rec[i]) * exp(-tr1->sumsite[i] * recrates[0]) * pow(recrates[0],tr2->rec[i]) * exp(-tr2->sumsite[i] * recrates[0]);
		temp1 = (alpha0 * (1 - lamda[0]) + alpha1 * lamda[1]) * pow(recrates[1],tr1->rec[i]) * exp(-tr1->sumsite[i] * recrates[1]) * pow(recrates[1],tr2->rec[i]) * exp(-tr2->sumsite[i] * recrates[1]);
		
		alpha0 = temp0;
		alpha1 = temp1;
	}
	
	return alpha0 + alpha1;
}

double hmmlikelihood(tree *target){
	double logalpha0 = 0.0;
	double logalpha1 = 0.0;
	double temp0,temp1;
	long i;
	double ini[2];
	double mean;
	
	ini[0] = (1 - lamda[1]) / ((1 - lamda[0]) + (1 - lamda[1]));
	ini[1] = 1 - ini[0];
	
	
	logalpha0 = log(ini[0] * pow(recrates[0],target->numrec_array[0]) * exp(-target->weight_array[0] * recrates[0]));
	logalpha1 = log(ini[1] * pow(recrates[1],target->numrec_array[0]) * exp(-target->weight_array[0] * recrates[1]));
	
	for (i = 1; i<(seq_length - 1); i++){
		mean = (logalpha0 + logalpha1) / 2.0;
		temp0 = (exp(logalpha0 - mean) * lamda[0] + exp(logalpha1 - mean) * (1 - lamda[1])) * pow(recrates[0],target->numrec_array[i]) * exp(-target->weight_array[i] * recrates[0]);
		temp1 = (exp(logalpha0 - mean) * (1 - lamda[0]) + exp(logalpha1 - mean) * lamda[1]) * pow(recrates[1],target->numrec_array[i]) * exp(-target->weight_array[i] * recrates[1]);
		
		logalpha0 = log(temp0) + mean;
		logalpha1 = log(temp1) + mean;
	}
	
	mean = (logalpha0 + logalpha1) / 2.0;
	
	return log(exp(logalpha0 - mean) + exp(logalpha1 - mean)) + mean;
}



void getbetaratio(tree *target, double *betaarray0, double *betaarray1){
	long i;
	double temp0,temp1;
	
	betaarray0[seq_length - 1] = 1;
	betaarray1[seq_length - 1] = 1;
	
	for (i = seq_length - 2; i >= 0; i--){
		temp0 = lamda[0] * pow(recrates[0],target->numrec_array[i + 1]) * exp(-target->weight_array[i + 1] * recrates[0])
		+ (1 - lamda[0]) * pow(recrates[1],target->numrec_array[i + 1]) * exp(-target->weight_array[i + 1] * recrates[1]) * betaarray1[i + 1] / betaarray0[i + 1];
		
		temp1 = (1 - lamda[1]) * pow(recrates[0],target->numrec_array[i + 1]) * exp(-target->weight_array[i + 1] * recrates[0])
		+ lamda[1] * pow(recrates[1],target->numrec_array[i + 1]) * exp(-target->weight_array[i + 1] * recrates[1]) * betaarray1[i + 1] / betaarray0[i + 1];
		
		betaarray0[i] = temp0 / (temp0 + temp1);
		betaarray1[i] = temp1 / (temp0 + temp1);
	}
	return;
}


void generaterecstructure(tree *target){
	double ini[2], p0, p1, *beta0, *beta1;
	long i;
	
	beta0 = (double *)calloc(seq_length, sizeof(double));
	beta1 = (double *)calloc(seq_length, sizeof(double));
	
	// beta ratio - ratio P(i == cold| sequence i+1 .....n), P(i == hot| sequence i+1.....n)
	getbetaratio(target, beta0, beta1);
	
	//getbeta(target, 0, beta0);
	//getbeta(target, 1, beta1);
	
	ini[0] = (1 - lamda[1]) / ((1 - lamda[0]) + (1 - lamda[1]));
	ini[1] = 1 - ini[0];
	
	
	// P(1 == state1, 2 == state2,......n == state n) == P(1 == state1 | sequence 1....n) * P(2 == state2| 1 == state1) * P(2== state2 | sequence 2....n) ..... P(n == state n | n - 1 == state n-1) * P(n == state n | sequence n)
	
	// P(i == state i|i - 1 == state i-1) * P(i == state i | sequence i+1 ... n)
	// = P(state i-1 -> state i) [lamda] * P(state i| sequence i) [prob rec] * beta[state i]
	// = lamda * (recrates)^number of recombination * exp(-sum of sites * recrates) * beta
	
	
	p0 = ini[0] * pow(recrates[0],target->numrec_array[0]) * exp(-target->weight_array[0] * recrates[0]) * beta0[0];
	p1 = ini[1] * pow(recrates[1],target->numrec_array[0]) * exp(-target->weight_array[0] * recrates[1]) * beta1[0];
	
	if (randum() <  p0/(p0+p1)) recarray[0] = recrates[0];
	else recarray[0] = recrates[1];
	
	for (i = 1; i < seq_length; i++){
		// i - 1 == cold state
		if (recarray[i - 1] == recrates[0]){
			
			p0 = lamda[0] * pow(recrates[0],target->numrec_array[i]) * exp(-target->weight_array[i] * recrates[0]) * beta0[i];
			p1 = (1 - lamda[0]) * pow(recrates[1],target->numrec_array[i]) * exp(-target->weight_array[i] * recrates[1]) * beta1[i];
			
			if (randum() < p0/(p0+p1)) recarray[i] = recrates[0];
			else recarray[i] = recrates[1];
		}
		// i - 1 == hot state
		else{
			p0 = (1 - lamda[1]) * pow(recrates[0],target->numrec_array[i]) * exp(-target->weight_array[i] * recrates[0]) * beta0[i];
			p1 = lamda[1] * pow(recrates[1],target->numrec_array[i]) * exp(-target->weight_array[i] * recrates[1]) * beta1[i];
			
			if (randum() < p0/(p0+p1)) recarray[i] = recrates[0];
			else recarray[i] = recrates[1];
		}
	}
	free(beta0);
	free(beta1);
	
	return;
}

void getbeta(tree *target, int state, double *betaarray){
	long i;
	double beta0 = 1.0;
	double beta1 = 1.0;
	double temp0,temp1;
	
	betaarray[seq_length - 1] = 1;
	
	for (i = seq_length - 2; i >= 0; i--){
		temp0 = lamda[0] * pow(recrates[0],target->numrec_array[i + 1]) * exp(-target->weight_array[i + 1] * recrates[0]) * beta0 +
		(1 - lamda[0]) * pow(recrates[1],target->numrec_array[i + 1]) * exp(-target->weight_array[i + 1] * recrates[1]) * beta1;
		
		temp1 = (1 - lamda[1]) * pow(recrates[0],target->numrec_array[i + 1]) * exp(-target->weight_array[i + 1] * recrates[0]) * beta0 +
		lamda[1] * pow(recrates[1],target->numrec_array[i + 1]) * exp(-target->weight_array[i + 1] * recrates[1]) * beta1;
		
		beta0 = temp0;
		beta1 = temp1;
		
		if (state == 0)  betaarray[i] = beta0;
		else betaarray[i] = beta1;
	}
	return;
}

double treepartlikelihood(tree *newtree, tree *oldtree){
	double *oldgpartrec, *oldgpartweight;
	long i;
	double mean, temp0, temp1, logalpha0, logalpha1, ini[2], oldglikelihood, newglikelihood;
	
	
	oldgpartweight = (double *)calloc(1,sizeof(seq_length));
	oldgpartrec = (double *)calloc(1,sizeof(seq_length));
	
	for (i = 0; i<seq_length; i++){
		oldgpartweight[i] = oldtree->weight_array[i] + gpartweight[i] - newtree->weight_array[i];
		oldgpartrec[i] = oldtree->numrec_array[i] + gpartrec[i] - newtree->numrec_array[i];
	}
	
	ini[0] = (1 - lamda[1]) / ((1 - lamda[0]) + (1 - lamda[1]));
	ini[1] = 1 - ini[0];
	
	
	logalpha0 = log(ini[0] * pow(recrates[0],oldgpartrec[0]) * exp(-oldgpartweight[0] * recrates[0]));
	logalpha1 = log(ini[1] * pow(recrates[1],oldgpartrec[0]) * exp(-oldgpartweight[0] * recrates[1]));
	
	for (i = 1; i<(seq_length - 1); i++){
		mean = (logalpha0 + logalpha1) / 2.0;
		temp0 = (exp(logalpha0 - mean) * lamda[0] + exp(logalpha1 - mean) * (1 - lamda[1])) * pow(recrates[0],oldgpartrec[i]) * exp(-oldgpartweight[i] * recrates[0]);
		temp1 = (exp(logalpha0 - mean) * (1 - lamda[0]) + exp(logalpha1 - mean) * lamda[1]) * pow(recrates[1],oldgpartrec[i]) * exp(-oldgpartweight[i] * recrates[1]);
		
		logalpha0 = log(temp0) + mean;
		logalpha1 = log(temp1) + mean;
	}
	
	mean = (logalpha0 + logalpha1) / 2.0;
	
	oldglikelihood = log(exp(logalpha0 - mean) + exp(logalpha1 - mean)) + mean;
	
	
	logalpha0 = log(ini[0] * pow(recrates[0],gpartrec[0]) * exp(-gpartweight[0] * recrates[0]));
	logalpha1 = log(ini[1] * pow(recrates[1],gpartrec[0]) * exp(-gpartweight[0] * recrates[1]));
	
	for (i = 1; i<(seq_length - 1); i++){
		mean = (logalpha0 + logalpha1) / 2.0;
		temp0 = (exp(logalpha0 - mean) * lamda[0] + exp(logalpha1 - mean) * (1 - lamda[1])) * pow(recrates[0],gpartrec[i]) * exp(-gpartweight[i] * recrates[0]);
		temp1 = (exp(logalpha0 - mean) * (1 - lamda[0]) + exp(logalpha1 - mean) * lamda[1]) * pow(recrates[1],gpartrec[i]) * exp(-gpartweight[i] * recrates[1]);
		
		logalpha0 = log(temp0) + mean;
		logalpha1 = log(temp1) + mean;
	}
	
	mean = (logalpha0 + logalpha1) / 2.0;
	
	newglikelihood = log(exp(logalpha0 - mean) + exp(logalpha1 - mean)) + mean;
	
	return oldglikelihood - newglikelihood;
}


void getbetaratio_temp(tree *target, double *betaarray0, double *betaarray1,double *rec, double *tran){
	long i;
	double temp0,temp1;
	
	betaarray0[seq_length - 1] = 1;
	betaarray1[seq_length - 1] = 1;
	
	for (i = seq_length - 2; i >= 0; i--){
		temp0 = tran[0] * pow(rec[0],target->numrec_array[i + 1]) * exp(-target->weight_array[i + 1] * rec[0])
		+ (1 - tran[0]) * pow(rec[1],target->numrec_array[i + 1]) * exp(-target->weight_array[i + 1] * rec[1]) * betaarray1[i + 1] / betaarray0[i + 1];
		
		temp1 = (1 - tran[1]) * pow(rec[0],target->numrec_array[i + 1]) * exp(-target->weight_array[i + 1] * rec[0])
		+ tran[1] * pow(rec[1],target->numrec_array[i + 1]) * exp(-target->weight_array[i + 1] * rec[1]) * betaarray1[i + 1] / betaarray0[i + 1];
		
		betaarray0[i] = temp0 / (temp0 + temp1);
		betaarray1[i] = temp1 / (temp0 + temp1);
	}
	return;
}


void generaterecstructure_temp(tree *target, double *rec, double *tran){
	double ini[2], p0, p1, *beta0, *beta1;
	long i;
	
	beta0 = (double *)calloc(seq_length, sizeof(double));
	beta1 = (double *)calloc(seq_length, sizeof(double));
	
	// beta ratio - ratio P(i == cold| sequence i+1 .....n), P(i == hot| sequence i+1.....n)
	getbetaratio_temp(target, beta0, beta1, rec, tran);
	
	//getbeta(target, 0, beta0);
	//getbeta(target, 1, beta1);
	
	ini[0] = (1 - tran[1]) / ((1 - tran[0]) + (1 - tran[1]));
	ini[1] = 1 - ini[0];
	
	// P(1 == state1, 2 == state2,......n == state n) == P(1 == state1 | sequence 1....n) * P(2 == state2| 1 == state1) * P(2== state2 | sequence 2....n) ..... P(n == state n | n - 1 == state n-1) * P(n == state n | sequence n)
	
	// P(i == state i|i - 1 == state i-1) * P(i == state i | sequence i+1 ... n)
	// = P(state i-1 -> state i) [tran] * P(state i| sequence i) [prob rec] * beta[state i]
	// = tran * (rec)^number of recombination * exp(-sum of sites * rec) * beta
	
	
	p0 = ini[0] * pow(rec[0],target->numrec_array[0]) * exp(-target->weight_array[0] * rec[0]) * beta0[0];
	p1 = ini[1] * pow(rec[1],target->numrec_array[0]) * exp(-target->weight_array[0] * rec[1]) * beta1[0];
	
	if (randum() <  p0/(p0+p1)) recarray[0] = rec[0];
	else recarray[0] = rec[1];
	
	for (i = 0; i < seq_length; i++){
		// i - 1 == cold state
		if (recarray[i - 1] == rec[0]){
			
			p0 = tran[0] * pow(rec[0],target->numrec_array[i]) * exp(-target->weight_array[i] * rec[0]) * beta0[i];
			p1 = (1 - tran[0]) * pow(rec[1],target->numrec_array[i]) * exp(-target->weight_array[i] * rec[1]) * beta1[i];
			
			if (randum() < p0/(p0+p1)) recarray[i] = rec[0];
			else recarray[i] = rec[1];
		}
		// i - 1 == hot state
		else{
			p0 = (1 - tran[1]) * pow(rec[0],target->numrec_array[i]) * exp(-target->weight_array[i] * rec[0]) * beta0[i];
			p1 = tran[1] * pow(rec[1],target->numrec_array[i]) * exp(-target->weight_array[i] * rec[1]) * beta1[i];
			
			if (randum() < p0/(p0+p1)) recarray[i] = rec[0];
			else recarray[i] = rec[1];
		}
	}
	free(beta0);
	free(beta1);
	
	return;
}



void changerecrates(long step){
	double cold_ini, hot_ini, rate;
	
	cold_ini = 1 / (1 - lamda[0]) /( 1 / (1 - lamda[0]) + 1 / (1 - lamda[1]));
	hot_ini = 1 - cold_ini;
	
	rate = cold_ini * recrates[0] + hot_ini * recrates[1];
	
	if (step < 5000){
		recrates[0] = recrates[1] = rate;
		return;
	}
	else{
		recrates[1] = (recrates1 + rate) / 2.0;
		recrates[0] = (recrates[1] * hot_ini + rate * cold_ini) / cold_ini;
		return;
	}
}



double calcr2(option_struct *op, data_fmt *data, long left, long right){ // just for diallelic markers
	double AB, Ab, aB, ab;
	char site1, site2;
	long i;
	char ****dna;
	
	dna = data->dnaptr->seqs;
	AB = Ab = aB = ab = 0;
	
	site1 = dna[0][0][0][left];
	site2 = dna[0][0][0][right];
	
	AB = 1;
	
	for(i = 1; i < getdata_numseq(op,data); i++){
		if (dna[0][0][i][left] == site1){
			if (dna[0][0][i][right] == site2) AB++;
			else Ab++;
		}
		else{
			if (dna[0][0][i][right] == site2) aB++;
			else ab++;
		}
	}
	AB = AB/ (double)(i);
	Ab = Ab/ (double)(i);
	aB = aB/ (double)(i);
	ab = ab/ (double)(i);
	return (AB * ab - Ab * aB) * (AB * ab - Ab * aB)/((AB+Ab) * (aB+ab) * (AB+aB) * (Ab + ab));
	
}
double calcD(option_struct *op, data_fmt *data, long left, long right){ // just for diallelic markers
	double AB, Ab, aB, ab, delta, min, D;
	char site1, site2;
	long i;
	char ****dna;
	
	dna = data->dnaptr->seqs;
	AB = Ab = aB = ab = 0;
	
	site1 = dna[0][0][0][left];
	site2 = dna[0][0][0][right];
	
	AB = 1;
	
	for(i = 1; i < getdata_numseq(op,data); i++){
		if (dna[0][0][i][left] == site1){
			if (dna[0][0][i][right] == site2) AB++;
			else Ab++;
		}
		else{
			if (dna[0][0][i][right] == site2) aB++;
			else ab++;
		}
	}
	AB = AB/ (double)(i);
	Ab = Ab/ (double)(i);
	aB = aB/ (double)(i);
	ab = ab/ (double)(i);
	
	delta = AB * ab - Ab * aB;
	
	if (delta < 0){
		if ((AB + Ab) * (AB + aB) < (aB + ab) * (Ab + ab)) min = (AB + Ab) * (AB + aB);
		else min = (aB + ab) * (Ab + ab);
	}
	else{
		if ((AB + Ab) * (Ab * ab) < (aB * ab) * (AB + aB)) min = (AB + Ab) * (Ab + ab);
		else min = (aB + ab) * (AB + aB);
	}
	
	D = delta / min;
	
	if (D < 0) D = -D;
	
	//printf(" %lf %lf %lf %lf %lf %lf %lf\n",AB,Ab,aB,ab,(AB + Ab), (AB + aB),D);
	
	return D;
	
}

void generaterecstructure_with_temp(tree *target,double temperature){
	double ini[2], p0, p1, *beta0, *beta1;
	long i;
	
	beta0 = (double *)calloc(seq_length, sizeof(double));
	beta1 = (double *)calloc(seq_length, sizeof(double));
	
	// beta ratio - ratio P(i == cold| sequence i+1 .....n), P(i == hot| sequence i+1.....n)
	getbetaratio_with_temp(target, beta0, beta1,temperature);
	
	//getbeta(target, 0, beta0);
	//getbeta(target, 1, beta1);
	
	ini[0] = (1 - lamda[1]) / ((1 - lamda[0]) + (1 - lamda[1]));
	ini[1] = 1 - ini[0];
	
	
	// P(1 == state1, 2 == state2,......n == state n) == P(1 == state1 | sequence 1....n) * P(2 == state2| 1 == state1) * P(2== state2 | sequence 2....n) ..... P(n == state n | n - 1 == state n-1) * P(n == state n | sequence n)
	
	// P(i == state i|i - 1 == state i-1) * P(i == state i | sequence i+1 ... n)
	// = P(state i-1 -> state i) [lamda] * P(state i| sequence i) [prob rec] * beta[state i]
	// = lamda * (recrates)^number of recombination * exp(-sum of sites * recrates) * beta
	
	
	p0 = ini[0] * pow(recrates[0],target->numrec_array[0] / temperature) * exp(-target->weight_array[0] * recrates[0] / temperature) * beta0[0];
	p1 = ini[1] * pow(recrates[1],target->numrec_array[0] / temperature) * exp(-target->weight_array[0] * recrates[1] / temperature) * beta1[0];
	
	if (randum() <  p0/(p0+p1)) recarray[0] = recrates[0];
	else recarray[0] = recrates[1];
	
	for (i = 1; i < seq_length; i++){
		// i - 1 == cold state
		if (recarray[i - 1] == recrates[0]){
			
			p0 = lamda[0] * pow(recrates[0],target->numrec_array[i] / temperature) * exp(-target->weight_array[i] * recrates[0] / temperature) * beta0[i];
			p1 = (1 - lamda[0]) * pow(recrates[1],target->numrec_array[i] / temperature) * exp(-target->weight_array[i] * recrates[1] / temperature) * beta1[i];
			
			if (randum() < p0/(p0+p1)) recarray[i] = recrates[0];
			else recarray[i] = recrates[1];
		}
		// i - 1 == hot state
		else{
			p0 = (1 - lamda[1]) * pow(recrates[0],target->numrec_array[i] / temperature) * exp(-target->weight_array[i] * recrates[0] / temperature) * beta0[i];
			p1 = lamda[1] * pow(recrates[1],target->numrec_array[i] / temperature) * exp(-target->weight_array[i] * recrates[1] / temperature) * beta1[i];
			
			if (randum() < p0/(p0+p1)) recarray[i] = recrates[0];
			else recarray[i] = recrates[1];
		}
	}
	free(beta0);
	free(beta1);
	
	return;
}


void getbetaratio_with_temp(tree *target, double *betaarray0, double *betaarray1, double temperature){
	long i;
	double temp0,temp1;
	
	betaarray0[seq_length - 1] = 1;
	betaarray1[seq_length - 1] = 1;
	
	for (i = seq_length - 2; i >= 0; i--){
		temp0 = lamda[0] * pow(recrates[0],target->numrec_array[i + 1] / temperature) * exp(-target->weight_array[i + 1] * recrates[0] / temperature)
		+ (1 - lamda[0]) * pow(recrates[1],target->numrec_array[i + 1] / temperature) * exp(-target->weight_array[i + 1] * recrates[1] / temperature) * betaarray1[i + 1] / betaarray0[i + 1];
		
		temp1 = (1 - lamda[1]) * pow(recrates[0],target->numrec_array[i + 1] / temperature) * exp(-target->weight_array[i + 1] * recrates[0] / temperature)
		+ lamda[1] * pow(recrates[1],target->numrec_array[i + 1] / temperature) * exp(-target->weight_array[i + 1] * recrates[1] / temperature) * betaarray1[i + 1] / betaarray0[i + 1];
		
		betaarray0[i] = temp0 / (temp0 + temp1);
		betaarray1[i] = temp1 / (temp0 + temp1);
	}
	return;
}

double hmmlikelihood_with_temp(tree *target, double temperature){
	double logalpha0 = 0.0;
	double logalpha1 = 0.0;
	double temp0,temp1;
	long i;
	double ini[2];
	double mean;
	
	ini[0] = (1 - lamda[1]) / ((1 - lamda[0]) + (1 - lamda[1]));
	ini[1] = 1 - ini[0];
	
	
	logalpha0 = log(ini[0] * pow(recrates[0],target->numrec_array[0] / temperature) * exp(-target->weight_array[0] * recrates[0] / temperature));
	logalpha1 = log(ini[1] * pow(recrates[1],target->numrec_array[0] / temperature) * exp(-target->weight_array[0] * recrates[1] / temperature));
	
	for (i = 1; i<(seq_length - 1); i++){
		mean = (logalpha0 + logalpha1) / 2.0;
		temp0 = (exp(logalpha0 - mean) * lamda[0] + exp(logalpha1 - mean) * (1 - lamda[1])) * pow(recrates[0],target->numrec_array[i] / temperature) * exp(-target->weight_array[i] * recrates[0] / temperature);
		temp1 = (exp(logalpha0 - mean) * (1 - lamda[0]) + exp(logalpha1 - mean) * lamda[1]) * pow(recrates[1],target->numrec_array[i] / temperature) * exp(-target->weight_array[i] * recrates[1] / temperature);
		
		logalpha0 = log(temp0) + mean;
		logalpha1 = log(temp1) + mean;
	}
	
	mean = (logalpha0 + logalpha1) / 2.0;
	
	return log(exp(logalpha0 - mean) + exp(logalpha1 - mean)) + mean;
}

double hmmlikelihood_two(tree *target, int *otherrecarray, double *otherweightarray, double temperature){
	double logalpha0 = 0.0;
	double logalpha1 = 0.0;
	double temp0,temp1;
	long i;
	double ini[2];
	double mean;
	
	ini[0] = (1 - lamda[1]) / ((1 - lamda[0]) + (1 - lamda[1]));
	ini[1] = 1 - ini[0];
	
	
	logalpha0 = log(ini[0] * pow(recrates[0],target->numrec_array[0] / temperature) * exp(-target->weight_array[0] * recrates[0] / temperature) * pow(recrates[0], otherrecarray[0]) * exp( - otherweightarray[0] * recrates[0]));
	logalpha1 = log(ini[1] * pow(recrates[1],target->numrec_array[0] / temperature) * exp(-target->weight_array[0] * recrates[1] / temperature)  * pow(recrates[1], otherrecarray[0]) * exp( - otherweightarray[0] * recrates[1]));
	
	for (i = 1; i<(seq_length - 1); i++){
		mean = (logalpha0 + logalpha1) / 2.0;
		temp0 = (exp(logalpha0 - mean) * lamda[0] + exp(logalpha1 - mean) * (1 - lamda[1])) * pow(recrates[0],target->numrec_array[i] / temperature) * exp(-target->weight_array[i] * recrates[0] / temperature) * pow(recrates[0], otherrecarray[i]) * exp( - otherweightarray[i] * recrates[0]);
		temp1 = (exp(logalpha0 - mean) * (1 - lamda[0]) + exp(logalpha1 - mean) * lamda[1]) * pow(recrates[1],target->numrec_array[i] / temperature) * exp(-target->weight_array[i] * recrates[1] / temperature)  * pow(recrates[1], otherrecarray[i]) * exp( - otherweightarray[i] * recrates[1]);
		
		logalpha0 = log(temp0) + mean;
		logalpha1 = log(temp1) + mean;
	}
	
	mean = (logalpha0 + logalpha1) / 2.0;
	
	return log(exp(logalpha0 - mean) + exp(logalpha1 - mean)) + mean;
}

void getoldg(tree *oldtree, tree *newtree){
	long i;
	
	for (i = 0; i < seq_length; i++){
		oldgrec[i] = oldtree->numrec_array[i] - newtree->numrec_array[i] + newgrec[i];
		oldgweight[i] = oldtree->weight_array[i] - newtree->weight_array[i] + newgweight[i];
	}
	return;
}

void temptreeswap_rec(option_struct *op, data_fmt *data, boolean *changed,long chain)
{
	long chl, chh;
	double chance, num, denom, templ, temph;
	boolean picked;
	tree *trh, *trl;
	
	if (op->numtempchains < 2) return;
	
	/* first pick which two chains to swap */
	
	picked = FALSE;
	while(!picked) {
		chl = (long)(randum()*op->numtempchains);
		chh = (long)(randum()*op->numtempchains);
		if (chl == op->numtempchains || chh == op->numtempchains) continue;
		if (chh - chl > 0) {
			if (chh - chl == 1) picked = TRUE;
		} else {
			if (chl - chh == 1) picked = TRUE;
		}
	}
	
	trh = temptrees[chh];
	temph = op->temperature[chh];
	trl = temptrees[chl];
	templ = op->temperature[chl];
	if (temph == 99 || temph == 98 || temph == 101) temph = 1;
	if (templ == 98 || templ == 99) templ = 1;
	
	
	num = hmmlikelihood_with_temp(trh,templ) +  hmmlikelihood_with_temp(trl,temph);
	denom = hmmlikelihood_with_temp(trh,temph) +  hmmlikelihood_with_temp(trl,templ); 
	
	chance = num - denom;
	
	swap++;
	if (chance >= log(randum())) {
		temptrees[chh] = trl;
		temptrees[chl] = trh;
		changed[chh] = TRUE;
		changed[chl] = TRUE;
		//printf("swap\n");
		swacc++;
	}
	
} 

void randomhotspotputter(){
	long hotsize, hotspot,start,end,i;
	
	hotsize = (long)(1.0 / (1.0 - lamda[1]));
	hotspot = (long)(seq_length * randum());
	
	for (i = 0; i < seq_length; i++) recarray[i] = recrates[0];
	
	if (hotspot - (long)((double)hotsize / 2.0) < 0) start = 0;
	else start = hotspot - (long)((double)hotsize / 2.0);
	
	if (hotspot + (long)((double)hotsize / 2.0) >= seq_length - 1) end = seq_length - 1;
	else end = hotspot + (long)((double)hotsize / 2.0);		   
	
	//printf("%ld %ld\n",start,end);
	
	for (i = start; i < end; i++) recarray[i] = recrates[1];
	
	return;
}

double probgS(int *recg, double *weightg){
	double result  = 0;
	long i;
	
	for (i=0; i<seq_length; i++){
		result += log(recarray[i]) * recg[i] - (recarray[i] * weightg[i]);
	}
	
	return result;
}


boolean randomhottestratio(tree *oldtree, tree *newtree)
{
	double test, x, numoldbranches, numnewbranches,num = 0, denom = 0;
	long i;
	
	test = (newtree->likelihood - oldtree->likelihood);  
	
	numoldbranches = (double)oldtree->numcoals * 2.0 + (double)oldtree->numrecombs;
	numnewbranches = (double)curtree->numcoals * 2.0 + (double)curtree->numrecombs;
	test += log(numoldbranches/numnewbranches);
	
	test += hmmlikelihood(newtree) - hmmlikelihood(oldtree);
	
	for (i = 0; i<seq_length; i++){
		num +=treeS(i,oldgrec,oldgweight);
		denom += treeS(i,newgrec,newgweight);
	}
	test += log(num) - log(denom);
	
	if (test >= 0) return TRUE;
	
	x = log(randum());
	if (x <= test)     return TRUE;
	
	else       return FALSE;
}


void stronghotspotputter(tree *target){
	long hotsize, hotspot,start,end,i,pos;
	double *alpha, *beta, *hotprob, total = 0;
	double hlikelihood;
	
	hotprob = (double *)calloc(seq_length,sizeof(double));
	
	hlikelihood = hmmlikelihood(target);
	alpha = getalpha1(target);
	beta = getbeta1(target);
	
	hotsize = (long)((1.0 / (1.0 - lamda[1])) * MAG);
	
	for (i = 0; i < seq_length; i++) recarray[i] = recrates[0];
	
	for (i = 0; i < seq_length; i++){
		hotprob[i] = exp(alpha[i] + beta[i] - hlikelihood);
		total += hotprob[i];
	}
	
	pos = randum() * total;
	
	for (i = 0; i < seq_length; i++){
		pos = pos - hotprob[i];
		if (pos > 0) break;
	}
	hotspot = i;
	
	if (i > 350 && i < 550)inhot++;
	else outhot++;
	
	if (hotspot - (long)((double)hotsize / 2.0) < 0) start = 0;
	else start = hotspot - (long)((double)hotsize / 2.0);
	
	if (hotspot + (long)((double)hotsize / 2.0) >= seq_length - 1) end = seq_length - 1;
	else end = hotspot + (long)((double)hotsize / 2.0);		   
	
	//printf("%ld %ld\n",start,end);
	
	for (i = start; i < end; i++) recarray[i] = recrates[1];
	
	free(alpha);
	free(beta);
	free(hotprob);
	return;
	
}

boolean hastingswithS(tree *oldtree, tree *newtree){
	long i;
	double *alpha, *beta, *hotprob, total = 0,oldtonew = 0,newtoold = 0;
	double ratio, x, hlikelihood;
	
	
	hotprob = (double *)calloc(seq_length,sizeof(double));
	hlikelihood = hmmlikelihood(newtree);
	
	alpha = getalpha1(oldtree);
	beta = getbeta1(oldtree);
	
	for (i = 0; i < seq_length; i++){
		hotprob[i] = exp(alpha[i] + beta[i] - hlikelihood);
		total += hotprob[i];
	}
	
	for (i = 0; i < seq_length; i++){
		oldtonew += (hotprob[i]/total) * treeS(i, newgrec, newgweight);
	}
	free(alpha);
	free(beta);
	
	hlikelihood = hmmlikelihood(newtree);
	alpha = getalpha1(newtree);
	beta = getbeta1(newtree);
	
	for (i = 0; i < seq_length; i++){
		hotprob[i] = exp(alpha[i] + beta[i] - hlikelihood);
		total += hotprob[i];
	}
	
	for (i = 0; i < seq_length; i++){
		newtoold += hotprob[i]/total * treeS(i, oldgrec, oldgweight);
	}
	
	free(alpha);
	free(beta);
	free(hotprob);
	
	ratio = (newtree->likelihood - oldtree->likelihood);  
	
	ratio += log(((double)oldtree->numcoals * 2.0 + (double)oldtree->numrecombs) / ((double)curtree->numcoals * 2.0 + (double)curtree->numrecombs));
	
	ratio += log(newtoold/ oldtonew);
	
	if (ratio >= 0.0)    return TRUE;
	else {
		x = log(randum());
		if (x <= ratio)     return TRUE;
		
		else       return FALSE;
	}
}



double treeS(long pos, int *recs, double *weights){
	double *recratesarray;
	double result = 0;
	long hotsize, start,end,i, hotspot;
	
	recratesarray = (double *)calloc(seq_length,sizeof(double));
	
	hotsize = (long)(1.0 / (1.0 - lamda[1]) * MAG);
	
	for (i = 0; i < seq_length; i++) recratesarray[i] = recrates[0];
	
	hotspot = pos;
	
	if (hotspot - (long)((double)hotsize / 2.0) < 0) start = 0;
	else start = hotspot - (long)((double)hotsize / 2.0);
	
	if (hotspot + (long)((double)hotsize / 2.0) >= seq_length - 1) end = seq_length - 1;
	else end = hotspot + (long)((double)hotsize / 2.0);		   
	
	//printf("%ld %ld\n",start,end);
	
	for (i = start; i < end; i++) recratesarray[i] = recrates[1];
	
	for (i = 0; i < seq_length; i++){
		result = log(recratesarray[i]) * recs[i] - (weights[i] * recratesarray[i]);
	}
	
	free(recratesarray); 
	return exp(result);
}


double *getalpha1(tree *tr){
	double mean;
	long i;
	double ini[2];
	double *a0,*a1;
	
	a0 = (double *)calloc(seq_length,sizeof(double));
	a1 = (double *)calloc(seq_length,sizeof(double));
	
	ini[0] = (1 - lamda[1]) / (1 - lamda[0] + 1 - lamda[1]);
	ini[1] = 1 - ini[0];
	
	a0[0] = log(ini[0]) + log(recrates[0]) * tr->numrec_array[0] - (tr->weight_array[0] * recrates[0]);
	a1[0] = log(ini[1]) + log(recrates[1]) * tr->numrec_array[0] - (tr->weight_array[0] * recrates[1]);
	
	for (i = 1; i<(seq_length - 1); i++){
		
		mean = (a0[i-1] + a1[i-1]) / 2.0;
		
		a0[i] = log(exp(a0[i-1] - mean) * lamda[0]  + exp(a1[i-1] - mean) * (1 - lamda[1])) + log(recrates[0]) * tr->numrec_array[i] - (tr->weight_array[i] * recrates[0]) + mean;
		
		a1[i] = log(exp(a0[i-1] - mean) * (1 - lamda[0])  + exp(a1[i-1] - mean) * lamda[1]) + log(recrates[1]) * tr->numrec_array[i] - (tr->weight_array[i] * recrates[1]) + mean;					
		
		
	}
	free(a0);
	return a1;
}

double *getbeta1(tree *tr){
	long i;
	double mean;
	double *b0,*b1;
	
	b0 = (double *)calloc(seq_length,sizeof(double));
	b1 = (double *)calloc(seq_length,sizeof(double));
	
	b0[seq_length - 1] = log(1);
	b1[seq_length - 1] = log(1);
	
	for (i = seq_length - 2; i >= 0; i--){
		mean = (b0[i+1] + b1[i+1]) / 2.0;
		
		b0[i] = log(lamda[0] * pow(recrates[0],tr->numrec_array[i+1]) * exp(- tr->weight_array[i+1] * recrates[0]) * exp(b0[i+1] - mean) +
					(1 - lamda[0]) * pow(recrates[1],tr->numrec_array[i+1]) * exp(- tr->weight_array[i+1] * recrates[1]) * exp(b1[i+1] - mean)) + mean;
		
		b1[i] = log(lamda[1] * pow(recrates[1],tr->numrec_array[i+1]) * exp(- tr->weight_array[i+1] * recrates[1]) * exp(b1[i+1] - mean) +
					(1 - lamda[1]) * pow(recrates[0],tr->numrec_array[i+1]) * exp(- tr->weight_array[i+1] * recrates[0]) * exp(b0[i+1] - mean)) + mean;
		
	}
	free(b0);
	return b1;
}

void printprobhot(tree *tr){
	long i;
	double *a, *b, hlikelihood;
	
	a = getalpha1(tr);
	b = getbeta1(tr);
	hlikelihood = hmmlikelihood(tr);
	
	for(i=0;i<seq_length;i++){
		printf("%lf\n",exp(a[i] + b[i] - hlikelihood));
	}
}

void sumreclength(tree *target, double *rec, double *length, long whichtree, long start, long end){
	long i;
	rec[whichtree] = length[whichtree] = 0;
	
	for (i = start; i<=end; i++){
		rec[whichtree] += target->numrec_array[i];
		length[whichtree] += target->weight_array[i];
	}
	return;
}



double sumbranch(tree *target){
	tlist *t;
	double tyme, tk = 0;
	for (t = target->tymelist; t->succ != NULL; t=t->succ) {
		tyme = t->age - t->eventnode->tyme;
		tk += t->numbranch * (t->numbranch - 1.0) * tyme;
	}
	return tk;
}

double erecs(option_struct *op, data_fmt *data, double *recs, double *lengths, double recvalue, long howmany){
	double  newrec , recnum, chances, *weight, newlikelihood,oldlikelihood, trec,twei;
	long i,j;
	
	newrec = recvalue;
	weight = (double *)calloc(howmany,sizeof(double));
	
	trec = twei = 0;
	
	oldlikelihood = 0;
	for (j = 0; j<5; j++){
		recnum = chances = 0;
		newlikelihood = 0;
		for (i = 0; i < howmany; i++){
			weight[i] = recs[i] * log(newrec) - lengths[i] * newrec - (recs[i] * log(recvalue) - lengths[i] * recvalue); 
			
		}
		for (i = 0; i < howmany; i++){
			recnum += exp(weight[i] - weight[0]) * recs[i];
			chances += exp(weight[i] - weight[0]) * lengths[i];
		}
		
		trec += recs[i];
		twei += lengths[i];
		//printf("recs %lf chance %lf\n",recnum,chances);
		
		if (recnum == 0) return 0.01;
		
		newrec = recnum/chances;
		oldlikelihood = newlikelihood;
	}
	for (i = 0; i < howmany; i++){
		trec += recs[i];
		twei += lengths[i];
		
	}
	
	//printf("tree %lf\n",trec/twei);
	free(weight);
	
	return newrec;   
}


void findhotspots(option_struct *op, data_fmt *data){
	long i,steps,numsamplings,j,k,*markerpos,step, sampling, tempnum, x;  
	double *iter,temp;
	boolean accepted;
	double *averecs, *recnums,*lengths, recvalue;
	FILE *recpattern, *userpattern;
	tree *temptree, *bestlikelihoodtree = NULL, *norectree;
	
	double  highrec = 0, lowrec = 100;
	
	
	norectree = copytree(op,data,curtree);
	
	steps = 11;
	
	numsamplings = tempnum = 10000;
	
	population = 0;
	constructseqinfo(op, data, 0.7);
	if (num_markers > 50){
		temp = (double)seq_length/50.0;
		x = 0;
		for (i = 0; i < seq_length; i++){
			if (dna_sequence[i] == 'S'){
				if (i - x < temp){
					dna_sequence[i] = 'A';
					num_markers --;
				}
				else{
					x = i;
				}
			}	  
		}
	}
	else{
		for (x = 0; x < 20; x ++){
			constructseqinfo(op, data,0.7 + 0.01 * x);
			if (num_markers > 25) break;
		}
	}
	
	iter = (double *)calloc(seq_length,sizeof(double));
	averecs = (double *)calloc(seq_length,sizeof(double));
	markerpos = (long *)calloc(num_markers,sizeof(long));
	//branches = (double *)calloc(numsamplings,sizeof(double));
	recnums = (double *)calloc(numsamplings,sizeof(double));
	lengths = (double *)calloc(numsamplings,sizeof(double));
	
	recpattern = fopen("inirates","w");
	
	for (i = 0 ; i<seq_length; i++){
		averecs[i] = 0;
		iter[i] = 0;
	}
	
	for (i = 0, j = 0; i <seq_length; i++){
		if (dna_sequence[i] == 'S'){
			markerpos[j] = i;
			j++;
		}
	}
	
	constructseqinfo(op, data,1.0);
	
	if (op->userrec == 0){
		
		
		
		for (i = 0 ; i <  j - 3 ; i++){
			recvalue = 0.1;
			marker1pos = markerpos[i];
			marker2pos = markerpos[i+1];
			marker3pos = markerpos[i+2];
			marker4pos = markerpos[i+3];
			
			curtree = copytree(op,data, norectree);
			
			for (k = 0; k<500; k++) accepted = makedropH(op,data);
			
			for (step = 0; step<steps ; step++){
				israndom = 0;
				for (k = 0; k <seq_length; k ++) recarray[k] = 0;
				for (k = marker1pos; k<marker4pos;k++) recarray[k] = recvalue;
				
				for (k = 0; k<1000; k++){
					accepted = makedropH(op,data); //for burnin	 
				}
				
				for (sampling = 0; sampling < numsamplings; sampling++){
					accepted = makedropH(op,data);
					sumreclength(curtree, recnums, lengths, sampling, marker1pos, marker4pos);
				}
				recvalue = erecs(op,data,recnums,lengths,recvalue, numsamplings);
				
				//printf("rec %lf  \n",recvalue);
			}
			printf("%ld %ld %lf\n",marker1pos, marker4pos,recvalue);
			for (k = marker1pos; k <= marker4pos; k++){
				averecs[k] += recvalue;
				iter[k]++;
				//printf("%ld %lf\n",k,thetas[k]);
			}
			
		}
		for (i = 0; i<seq_length;i++){
			if (iter[i] != 0) averecs[i] = averecs[i]/iter[i];
			if (highrec < averecs[i]) highrec= averecs[i];
			if (lowrec > averecs[i]) lowrec = averecs[i];
			//printf("%lf\n",averecs[i]);
			fprintf(recpattern,"%lf\n",averecs[i]);
		}
		israndom = 0;
		//free(branches);
		
		constructseqinfo(op, data,1);
		//recrates[0] = lowrec;
		//recrates[1] = highrec;
		
		// estimate inial parameters and make new initial tree. 
	}
	else{
		
		// if usertree
		userpattern = fopen("userrates","r");
		for (i = 0; i<seq_length; i++){
			fscanf(userpattern,"%lf",&averecs[i]);
		}
		fclose(userpattern);
	}
	fclose(recpattern);
	
	marker1pos = marker2pos = -99;
	constructseqinfo(op, data,1);
	temptree = (tree *)calloc(1,sizeof(tree));
	
	temptree->numrec_array = (double *)calloc(seq_length,sizeof(double));
	temptree->weight_array = (double *)calloc(seq_length,sizeof(double));
	
	
	for (i = 0; i <seq_length; i++){
		temptree->weight_array[i] = 1.0;
		temptree->numrec_array[i] = averecs[i];
	}
	// estimate initial parameters.
	/*   getinirate(temptree->numrec_array, temptree->weight_array); */
	
	/*   printf("%lf %lf %lf %lf %lf\n",theta0, recrates[0], recrates[1], lamda[0], lamda[1]); */
	
	
	curtree = copytree(op,data,norectree);
	freetree(op,data,norectree);
	if (bestlikelihoodtree != NULL) freetree(op,data,bestlikelihoodtree);
	bestlikelihoodtree = copytree(op,data,curtree);
	
	
	op->ctemp = 1;
	
	for (i = 0; i < seq_length; i++){
		if (highrec < averecs[i]) highrec= averecs[i];
	}
	recrates[1] = highrec;
	
	generaterecstructure_with_temp(temptree, op->ctemp);
	// get initial tree.
	for (i = 0 ; i< 5000; i++){
		accepted = makedropH(op,data);
		//printf(" %lf\n ",curtree->likelihood);
		if (curtree->likelihood > bestlikelihoodtree->likelihood){
			freetree(op,data,bestlikelihoodtree);
			bestlikelihoodtree = copytree(op,data,curtree);
		}
	}
	
	op->ctemp = 1; 
	
	freetree(op,data,curtree);
	free(markerpos);
	free(averecs);
	free(iter);
	free(recnums);
	free(lengths);
	free(temptree->numrec_array);
	free(temptree->weight_array);
	free(temptree);
	curtree = copytree(op,data,bestlikelihoodtree);
	freetree(op,data,bestlikelihoodtree);
	
	printf("rec %ld\n",curtree->numrecombs);
	
	return;
}

char checkSNP(option_struct *op, data_fmt *data, long position, double freq){ // check SNP
	long i;
	char ****y;
	int nA, nC, nG, nT,nN;
	
	nA = nC = nG = nT = nN = 0;
	//num_seq = getdata_numseq(op,data);
	y = data->dnaptr->seqs;
	for (i = 0; i < seq_num; i++) {
		switch(y[0][0][i][position]){
			case 'N':
				nN ++;
				break;
			case 'A':
				nA++; 
				break; 
			case 'C':
				nC ++;
				break;
			case 'G':
				nG ++;
				break;
			case 'T':
				nT++;
				break;
		} 
	}
	if ((nN > 0) && (nN + nA == seq_num)) return '1';
	if ((nN > 0) && (nN + nC == seq_num)) return '2';
	if ((nN > 0) && (nN + nG == seq_num)) return '3';
	if ((nN > 0) && (nN + nT == seq_num)) return '4';
	
	if (nN == seq_num) return 'N';
	if (nA == seq_num) return 'A';
	if (nC == seq_num) return 'C';
	if (nG == seq_num) return 'G';
	if (nT == seq_num) return 'T';
	
	return 'S';
	
}
void getinirate(double *recs, double *weights){
	double *alpha0, *alpha1, *beta0, *beta1;
	double temprecrates[2], templamda[2];
	double recevents, sumsites, stay, transfer,temp;
	long i,j;
	
	alpha0 = (double *)calloc(seq_length,sizeof(double));
	alpha1 = (double *)calloc(seq_length,sizeof(double));
	beta0 = (double *)calloc(seq_length,sizeof(double));
	beta1 = (double *)calloc(seq_length,sizeof(double));
	
	temprecrates[0] = recrates[0];
	temprecrates[1] = recrates[1];
	templamda[0] = lamda[0];
	templamda[1] = lamda[1];
	
	
	for (j = 0; j < 10; j++){
		
		recevents = sumsites = 0;
		getalphabeta(recs, weights, temprecrates, templamda, alpha0, alpha1, beta0, beta1);
		for (i = 0; i < seq_length; i++){
			recevents += 1.0 / (1.0 + exp(alpha0[i] + beta0[i] - alpha1[i] - beta1[i])) * recs[i];
			sumsites += 1.0 / (1.0 + exp(alpha0[i] + beta0[i] - alpha1[i] - beta1[i])) * weights[i];
		}
		temprecrates[1] = recevents/sumsites;
		
		recevents = sumsites = 0;
		getalphabeta(recs, weights, temprecrates, templamda, alpha0, alpha1, beta0, beta1);
		for (i = 0; i < seq_length; i++){
			recevents += 1.0 / (1.0 + exp(alpha1[i] + beta1[i] - alpha0[i] - beta0[i])) * recs[i];
			sumsites += 1.0 / (1.0 + exp(alpha1[i] + beta1[i] - alpha0[i] - beta0[i])) * weights[i];
		}
		temprecrates[0] = recevents/sumsites;
		
		
		stay = transfer = 0;
		getalphabeta(recs, weights, temprecrates, templamda, alpha0, alpha1, beta0, beta1);
		for (i = 0; i < seq_length - 1; i++){
			stay += exp(alpha0[i] + log(temprecrates[0]) * recs[i + 1] - (temprecrates[0] * weights[i + 1]) + log(templamda[0]) + beta0[i + 1] - beta0[0]);
			
			transfer += exp(alpha0[i] + log(temprecrates[1]) * recs[i + 1] - (temprecrates[1] * weights[i + 1]) + log(1 - templamda[0]) + beta1[i + 1] - beta0[0]);
			
		}
		//printf("%lf %lf\n",stay,transfer);
		templamda[0] = stay/(stay + transfer);
		
		stay = transfer = 0;
		getalphabeta(recs, weights, temprecrates, templamda, alpha0, alpha1, beta0, beta1);
		for (i = 0; i < seq_length - 1; i++){
			stay += exp(alpha1[i] + log(temprecrates[1]) * recs[i + 1] - (temprecrates[1] * weights[i + 1]) + log(templamda[1]) + beta1[i + 1] - beta0[0]);
			
			transfer += exp(alpha1[i] + log(temprecrates[0]) * recs[i + 1] - (temprecrates[0] * weights[i + 1]) + log(1 - templamda[1]) + beta0[i + 1] - beta0[0]);
			
		}
		templamda[1] = stay/(stay + transfer);
		
		//printf("%lf %lf %lf %lf\n",temprecrates[0], temprecrates[1], templamda[0], templamda[1]);
	}
	recrates[0] = temprecrates[0];
	recrates[1] = temprecrates[1];
	lamda[0] = templamda[0];
	lamda[1] = templamda[1];
	
	if (recrates[0] > recrates[1]){
		temp = recrates[0];
		recrates[0] = recrates[1];
		recrates[1] = temp;
		temp = lamda[0];
		lamda[0] = lamda[1];
		lamda[1] = temp;
	}
	
	free(alpha0);
	free(alpha1);
	free(beta0);
	free(beta1);
	
	return;
}


void getalphabeta(double *recs, double *weights, double *recombrates, double *lamdas, double *alpha0, double *alpha1, double *beta0, double *beta1){
	double ini[2],mean;
	long i;
	
	// get alpha
	ini[0] = (1 - lamdas[1]) / (1 - lamdas[0] + 1 - lamdas[1]);
	ini[1] = 1 - ini[0];
	
	alpha0[0] = log(ini[0]) + log(recombrates[0]) * recs[0] - (weights[0] * recombrates[0]);
	alpha1[0] = log(ini[1]) + log(recombrates[1]) * recs[0] - (weights[0] * recombrates[1]);
	
	for (i = 1; i<(seq_length - 1); i++){
		
		mean = (alpha0[i-1] + alpha1[i-1]) / 2.0;
		
		alpha0[i] = log(exp(alpha0[i-1] - mean) * lamdas[0]  + exp(alpha1[i-1] - mean) * (1 - lamdas[1])) + log(recombrates[0]) * recs[i] - (weights[i] * recombrates[0]) + mean;
		
		alpha1[i] = log(exp(alpha0[i-1] - mean) * (1 - lamdas[0])  + exp(alpha1[i-1] - mean) * lamdas[1]) + log(recombrates[1]) * recs[i] - (weights[i] * recombrates[1]) + mean;					
	}
	
	// get beta
	
	beta0[seq_length - 1] = log(1);
	beta1[seq_length - 1] = log(1);
	
	for (i = seq_length - 2; i >= 0; i--){
		mean = (beta0[i+1] + beta1[i+1]) / 2.0;
		
		beta0[i] = log(lamdas[0] * pow(recombrates[0],recs[i+1]) * exp(- weights[i+1] * recombrates[0]) * exp(beta0[i+1] - mean) +
					   (1 - lamdas[0]) * pow(recombrates[1],recs[i+1]) * exp(- weights[i+1] * recombrates[1]) * exp(beta1[i+1] - mean)) + mean;
		
		beta1[i] = log(lamdas[1] * pow(recombrates[1],recs[i+1]) * exp(- weights[i+1] * recombrates[1]) * exp(beta1[i+1] - mean) +
					   (1 - lamdas[1]) * pow(recombrates[0],recs[i+1]) * exp(- weights[i+1] * recombrates[0]) * exp(beta0[i+1] - mean)) + mean;
		
	}
	
	return;
	
}



double flatrecombinationrateratio(tree *oldtree, tree *newtree){
	double test,x, sumoldrec=0, sumoldweight=0, sumnewrec = 0, sumnewweight=0;
	long i;
	
	test = (newtree->likelihood - oldtree->likelihood);  
	
	test += hmmlikelihood(newtree) - hmmlikelihood(oldtree);
	
	for (i = 0; i<seq_length; i++){
		sumoldrec += oldgrec[i];
		sumoldweight += oldgweight[i];
		sumnewrec += newgrec[i];
		sumnewweight += newgweight[i];
	}
	test += (sumoldrec * log(recarray[0]) - (sumoldweight * recarray[0])) - (sumnewrec * log(recarray[0]) - (sumnewweight * recarray[0]));
	
	if (test >= 0) return TRUE;
	
	x = log(randum());
	if (x <= test)     return TRUE;
	
	else       return FALSE;
}



// generate structure with different parameters - bigger lamdas & bigger recH
void generaterecstructure_with_big_lamda(tree *target, double temperature){
	double ini[2], p0, p1, *beta0, *beta1, templamda[2],temprecrates_h;
	long i;
	
	templamda[0] = lamda[0];
	templamda[1] = lamda[1];
	
	lamda[0] = (templamda[0] + 0.25) / 1.25;
	lamda[1] = (templamda[1] + 0.25) / 1.25;
	temprecrates_h = recrates[1] * 1.5;
	
	beta0 = (double *)calloc(seq_length, sizeof(double));
	beta1 = (double *)calloc(seq_length, sizeof(double));
	
	// beta ratio - ratio P(i == cold| sequence i+1 .....n), P(i == hot| sequence i+1.....n)
	getbetaratio_with_temp(target, beta0, beta1,temperature);
	
	//getbeta(target, 0, beta0);
	//getbeta(target, 1, beta1);
	
	ini[0] = (1 - lamda[1]) / ((1 - lamda[0]) + (1 - lamda[1]));
	ini[1] = 1 - ini[0];
	
	
	// P(1 == state1, 2 == state2,......n == state n) == P(1 == state1 | sequence 1....n) * P(2 == state2| 1 == state1) * P(2== state2 | sequence 2....n) ..... P(n == state n | n - 1 == state n-1) * P(n == state n | sequence n)
	
	// P(i == state i|i - 1 == state i-1) * P(i == state i | sequence i+1 ... n)
	// = P(state i-1 -> state i) [lamda] * P(state i| sequence i) [prob rec] * beta[state i]
	// = lamda * (recrates)^number of recombination * exp(-sum of sites * recrates) * beta
	
	
	p0 = ini[0] * pow(recrates[0],target->numrec_array[0] / temperature) * exp(-target->weight_array[0] * recrates[0] / temperature) * beta0[0];
	p1 = ini[1] * pow(recrates[1],target->numrec_array[0] / temperature) * exp(-target->weight_array[0] * recrates[1] / temperature) * beta1[0];
	
	if (randum() <  p0/(p0+p1)) recarray[0] = recrates[0];
	else recarray[0] = temprecrates_h;
	
	for (i = 1; i < seq_length; i++){
		// i - 1 == cold state
		if (recarray[i - 1] == recrates[0]){
			
			p0 = lamda[0] * pow(recrates[0],target->numrec_array[i] / temperature) * exp(-target->weight_array[i] * recrates[0] / temperature) * beta0[i];
			p1 = (1 - lamda[0]) * pow(recrates[1],target->numrec_array[i] / temperature) * exp(-target->weight_array[i] * recrates[1] / temperature) * beta1[i];
			
			if (randum() < p0/(p0+p1)) recarray[i] = recrates[0];
			else recarray[i] = recrates[1];
		}
		// i - 1 == hot state
		else{
			p0 = (1 - lamda[1]) * pow(recrates[0],target->numrec_array[i] / temperature) * exp(-target->weight_array[i] * recrates[0] / temperature) * beta0[i];
			p1 = lamda[1] * pow(recrates[1],target->numrec_array[i] / temperature) * exp(-target->weight_array[i] * recrates[1] / temperature) * beta1[i];
			
			if (randum() < p0/(p0+p1)) recarray[i] = recrates[0];
			else recarray[i] = temprecrates_h;
		}
	}
	free(beta0);
	free(beta1);
	
	lamda[0] = templamda[0];
	lamda[1] = templamda[1];
	
	
	return;
}


double hmmlikelihood_two_with_bigger_rec(tree *target, int *otherrecarray, double *otherweightarray){ // hastings ratio for bigger & hotter case.
	double logalpha0 = 0.0;
	double logalpha1 = 0.0;
	double temp0,temp1; 
	long i;
	double ini[2], templamda[2];
	double mean;
	double highrec[2];
	
	templamda[0] = lamda[0];
	templamda[1] = lamda[1];
	
	lamda[0] = (templamda[0] + 0.25) / 1.25;
	lamda[1] = (templamda[1] + 0.25) / 1.25;
	
	highrec[0] = recrates[0];
	highrec[1] = recrates[1] * 1.5;
	
	ini[0] = (1 - lamda[1]) / ((1 - lamda[0]) + (1 - lamda[1]));
	ini[1] = 1 - ini[0];
	
	
	logalpha0 = log(ini[0] * pow(recrates[0],target->numrec_array[0]) * exp(-target->weight_array[0] * recrates[0]) * pow(highrec[0], otherrecarray[0]) * exp( - otherweightarray[0] * highrec[0]));
	logalpha1 = log(ini[1] * pow(recrates[1],target->numrec_array[0]) * exp(-target->weight_array[0] * recrates[1])  * pow(highrec[1], otherrecarray[0]) * exp( - otherweightarray[0] * highrec[1]));
	
	for (i = 1; i<(seq_length - 1); i++){
		mean = (logalpha0 + logalpha1) / 2.0;
		temp0 = (exp(logalpha0 - mean) * lamda[0] + exp(logalpha1 - mean) * (1 - lamda[1])) * pow(recrates[0],target->numrec_array[i]) * exp(-target->weight_array[i] * recrates[0]) * pow(highrec[0], otherrecarray[i]) * exp( - otherweightarray[i] * highrec[0]);
		temp1 = (exp(logalpha0 - mean) * (1 - lamda[0]) + exp(logalpha1 - mean) * lamda[1]) * pow(recrates[1],target->numrec_array[i]) * exp(-target->weight_array[i] * recrates[1]) * pow(highrec[1], otherrecarray[i]) * exp( - otherweightarray[i] * highrec[1]);
		
		logalpha0 = log(temp0) + mean;
		logalpha1 = log(temp1) + mean;
	}
	
	mean = (logalpha0 + logalpha1) / 2.0;
	
	return log(exp(logalpha0 - mean) + exp(logalpha1 - mean)) + mean;
}


void simnodebase(option_struct *op, data_fmt *data, nodeinfo *target, nodeinfo *ancestor, double length, int site){
	double probs[4];
	double tempprob[4];
	double prand, tempv; 
	double total;
	int newbase;
	double lw1 = 0.0;
	
	double yy1, ww1zz1, vv1zz1, vv1zz1_sumr1, vv1zz1_sumy1, sum1, sumr1, sumy1;
	dnadata *dna;
	
	
	dna = data->dnaptr;
	
	probs[0] = 1.0;
	probs[1] = 1.0;
	probs[2] = 1.0;
	probs[3] = 1.0;
	
	if (target->istip == 1){
		if (target->nodenum <= old_hap){
			return;
		}
	}
	
	if (ancestor !=NULL){ // probability to top ?
		
		
		lw1 = -1.0 * length / data->dnaptr->fracchange;
		
		tbl[0].ww1 = exp(tbl[0].rat_xi * lw1);
		tbl[0].zz1 = exp(tbl[0].rat_xv * lw1);
		tbl[0].ww1zz1 = tbl[0].ww1 * tbl[0].zz1;
		tbl[0].vv1zz1 = (1.0 - tbl[0].ww1) * tbl[0].zz1;
		
		
		
		dna = data->dnaptr;
		
		
		ww1zz1 = tbl[0].ww1zz1;
		vv1zz1 = tbl[0].vv1zz1;
		yy1 = 1.0 - tbl[0].zz1;
		
		// 'A' P(parent|A')
		
		sum1 = yy1 * dna->freqa;
		
		sumr1 = dna->freqar;
		sumy1 = 0;
		
		vv1zz1_sumr1 = vv1zz1 * sumr1;
		vv1zz1_sumy1 = vv1zz1 * sumy1;
		
		tempprob[0] =  (sum1 + ww1zz1 + vv1zz1_sumr1);
		tempprob[1] =  (sum1 + vv1zz1_sumy1);
		tempprob[2] =  (sum1 + vv1zz1_sumr1);
		tempprob[3] =  (sum1 + vv1zz1_sumy1);
		
		total = tempprob[0] + tempprob[1] + tempprob[2] + tempprob[3];
		probs[0] = (tempprob[0] * ancestor->A + tempprob[1] * ancestor->C + tempprob[2] * ancestor->G + tempprob[3] * ancestor->T) /total; 
		//if (site == ts) printf("%lf %lf %lf %lf %lf\n",tempprob[0], tempprob[1], tempprob[2], tempprob[3],probs[0]);
		// P(A'|offspring(s))
		probs[0] = probs[0] * target->A;
		
		// 'C'
		sum1 = yy1 * dna->freqc;
		
		sumr1 = 0;
		sumy1 = dna->freqcy;
		
		vv1zz1_sumr1 = vv1zz1 * sumr1;
		vv1zz1_sumy1 = vv1zz1 * sumy1;
		
		tempprob[0] =  (sum1 + vv1zz1_sumr1);
		tempprob[1] =  (sum1 + ww1zz1 + vv1zz1_sumy1);
		tempprob[2] =  (sum1 + vv1zz1_sumr1);
		tempprob[3] =  (sum1 + vv1zz1_sumy1);
		
		total = tempprob[0] + tempprob[1] + tempprob[2] + tempprob[3];
		probs[1] = (tempprob[0] * ancestor->A + tempprob[1] * ancestor->C + tempprob[2] * ancestor->G + tempprob[3] * ancestor->T) /total; 
		//if (site == ts) printf("%lf %lf %lf %lf %lf\n",tempprob[0], tempprob[1], tempprob[2], tempprob[3],probs[1]);
		// P(G'|offspring(s))
		probs[1] = probs[1] * target->C;
		
		// 'G'
		sum1 = yy1 * dna->freqg;
		
		sumr1 = dna->freqgr;
		sumy1 = 0;
		
		vv1zz1_sumr1 = vv1zz1 * sumr1;
		vv1zz1_sumy1 = vv1zz1 * sumy1;
		
		tempprob[0] =  (sum1 + vv1zz1_sumr1);
		tempprob[1] =  (sum1 + vv1zz1_sumy1);
		tempprob[2] =  (sum1 + ww1zz1 + vv1zz1_sumr1);
		tempprob[3] =  (sum1 + vv1zz1_sumy1);
		
		total = tempprob[0] + tempprob[1] + tempprob[2] + tempprob[3];
		probs[2] = (tempprob[0] * ancestor->A + tempprob[1] * ancestor->C + tempprob[2] * ancestor->G + tempprob[3] * ancestor->T) /total; 
		//if (site == ts) printf("%lf %lf %lf %lf %lf\n",tempprob[0], tempprob[1], tempprob[2], tempprob[3],probs[2]);
		
		
		// P(G'|offspring(s))
		probs[2] = probs[2] * target->G;
		
		// 'T'
		sum1 = yy1 * dna->freqt;
		
		sumr1 = 0;
		sumy1 = dna->freqty;
		
		vv1zz1_sumr1 = vv1zz1 * sumr1;
		vv1zz1_sumy1 = vv1zz1 * sumy1;
		
		tempprob[0] =  (sum1 + vv1zz1_sumr1);
		tempprob[1] =  (sum1 + vv1zz1_sumy1);
		tempprob[2] =  (sum1 + vv1zz1_sumr1);
		tempprob[3] =  (sum1 + ww1zz1 + vv1zz1_sumy1);
		
		total = tempprob[0] + tempprob[1] + tempprob[2] + tempprob[3];
		probs[3] = (tempprob[0] * ancestor->A + tempprob[1] * ancestor->C + tempprob[2] * ancestor->G + tempprob[3] * ancestor->T) /total; 
		//if (site == ts) printf("%lf %lf %lf %lf %lf\n",tempprob[0], tempprob[1], tempprob[2], tempprob[3],probs[3]);
		
		// P(G'|offspring(s))
		probs[3] = probs[3] * target->T;
		
	    
	}
	else{
		probs[0] = target->A;
		probs[1] = target->C;
		probs[2] = target->G;
		probs[3] = target->T;
	}
	
	total = probs[0] + probs[1] + probs[2] + probs[3];
	for (newbase = 0; newbase < 4; newbase++){
		probs[newbase] = probs[newbase]/total;
	}
	prand = randum();
	tempv = 0;
	for (newbase = 0; newbase < 4; newbase++){
		tempv = tempv + probs[newbase];
		if (prand <= tempv)
			break;
	}
	//printf("%ld: %d, %lf, %lf, %lf, %lf, %lf\n",target->nodenum, i, p, probs[0]/total,probs[1]/total,probs[2]/total,probs[3]/total);
	
	
	target->A = target->C = target->G = target->T = 0;
	
	switch (newbase){
		case 0:
			target->A = 1;
			break;
		case 1:
			target->C = 1;
			break;
		case 2:
			target->G = 1;
			break;
		case 3:
			target->T = 1;
			break;
	}
	
	
	if (target->istip == 1){ // tipnode
		if (target->nodenum > old_hap){ // new haplotype
			//newsimdata[target->nodenum - old_hap][site] = i+1;
			newhaps[target->nodenum - old_hap] = newbase + 1;
		}
	}
	
	return;
}

void simnodebase_old(option_struct *op, data_fmt *data, nodeinfo *target, nodeinfo *ancestor, double length, int site){
	double probs[4];
	double tempprob[4];
	double prand, tempv; 
	double total;
	int newbase;
	double lw1 = 0.0;
	
	double yy1, ww1zz1, vv1zz1, vv1zz1_sumr1, vv1zz1_sumy1, sum1, sumr1, sumy1;
	dnadata *dna;
	
	
	dna = data->dnaptr;
	
	probs[0] = 1.0;
	probs[1] = 1.0;
	probs[2] = 1.0;
	probs[3] = 1.0;
	
	if (ancestor != NULL)
		target->A = target->C = target->G = target->T = 1;
	
	if (ancestor !=NULL){ // probability to top ?
		
		
		lw1 = -1.0 * length / data->dnaptr->fracchange;
		
		tbl[0].ww1 = exp(tbl[0].rat_xi * lw1);
		tbl[0].zz1 = exp(tbl[0].rat_xv * lw1);
		tbl[0].ww1zz1 = tbl[0].ww1 * tbl[0].zz1;
		tbl[0].vv1zz1 = (1.0 - tbl[0].ww1) * tbl[0].zz1;
		
		
		
		dna = data->dnaptr;
		
		
		ww1zz1 = tbl[0].ww1zz1;
		vv1zz1 = tbl[0].vv1zz1;
		yy1 = 1.0 - tbl[0].zz1;
		
		// 'A' P(parent|A')
		
		sum1 = yy1 * dna->freqa;
		
		sumr1 = dna->freqar;
		sumy1 = 0;
		
		vv1zz1_sumr1 = vv1zz1 * sumr1;
		vv1zz1_sumy1 = vv1zz1 * sumy1;
		
		tempprob[0] =  (sum1 + ww1zz1 + vv1zz1_sumr1);
		tempprob[1] =  (sum1 + vv1zz1_sumy1);
		tempprob[2] =  (sum1 + vv1zz1_sumr1);
		tempprob[3] =  (sum1 + vv1zz1_sumy1);
		
		total = tempprob[0] + tempprob[1] + tempprob[2] + tempprob[3];
		probs[0] = (tempprob[0] * ancestor->A + tempprob[1] * ancestor->C + tempprob[2] * ancestor->G + tempprob[3] * ancestor->T) /total; 
		//if (site == ts) printf("%lf %lf %lf %lf %lf\n",tempprob[0], tempprob[1], tempprob[2], tempprob[3],probs[0]);
		// P(A'|offspring(s))
		probs[0] = probs[0] * target->A;
		
		// 'C'
		sum1 = yy1 * dna->freqc;
		
		sumr1 = 0;
		sumy1 = dna->freqcy;
		
		vv1zz1_sumr1 = vv1zz1 * sumr1;
		vv1zz1_sumy1 = vv1zz1 * sumy1;
		
		tempprob[0] =  (sum1 + vv1zz1_sumr1);
		tempprob[1] =  (sum1 + ww1zz1 + vv1zz1_sumy1);
		tempprob[2] =  (sum1 + vv1zz1_sumr1);
		tempprob[3] =  (sum1 + vv1zz1_sumy1);
		
		total = tempprob[0] + tempprob[1] + tempprob[2] + tempprob[3];
		probs[1] = (tempprob[0] * ancestor->A + tempprob[1] * ancestor->C + tempprob[2] * ancestor->G + tempprob[3] * ancestor->T) /total; 
		//if (site == ts) printf("%lf %lf %lf %lf %lf\n",tempprob[0], tempprob[1], tempprob[2], tempprob[3],probs[1]);
		// P(G'|offspring(s))
		probs[1] = probs[1] * target->C;
		
		// 'G'
		sum1 = yy1 * dna->freqg;
		
		sumr1 = dna->freqgr;
		sumy1 = 0;
		
		vv1zz1_sumr1 = vv1zz1 * sumr1;
		vv1zz1_sumy1 = vv1zz1 * sumy1;
		
		tempprob[0] =  (sum1 + vv1zz1_sumr1);
		tempprob[1] =  (sum1 + vv1zz1_sumy1);
		tempprob[2] =  (sum1 + ww1zz1 + vv1zz1_sumr1);
		tempprob[3] =  (sum1 + vv1zz1_sumy1);
		
		total = tempprob[0] + tempprob[1] + tempprob[2] + tempprob[3];
		probs[2] = (tempprob[0] * ancestor->A + tempprob[1] * ancestor->C + tempprob[2] * ancestor->G + tempprob[3] * ancestor->T) /total; 
		//if (site == ts) printf("%lf %lf %lf %lf %lf\n",tempprob[0], tempprob[1], tempprob[2], tempprob[3],probs[2]);
		
		
		// P(G'|offspring(s))
		probs[2] = probs[2] * target->G;
		
		// 'T'
		sum1 = yy1 * dna->freqt;
		
		sumr1 = 0;
		sumy1 = dna->freqty;
		
		vv1zz1_sumr1 = vv1zz1 * sumr1;
		vv1zz1_sumy1 = vv1zz1 * sumy1;
		
		tempprob[0] =  (sum1 + vv1zz1_sumr1);
		tempprob[1] =  (sum1 + vv1zz1_sumy1);
		tempprob[2] =  (sum1 + vv1zz1_sumr1);
		tempprob[3] =  (sum1 + ww1zz1 + vv1zz1_sumy1);
		
		total = tempprob[0] + tempprob[1] + tempprob[2] + tempprob[3];
		probs[3] = (tempprob[0] * ancestor->A + tempprob[1] * ancestor->C + tempprob[2] * ancestor->G + tempprob[3] * ancestor->T) /total; 
		//if (site == ts) printf("%lf %lf %lf %lf %lf\n",tempprob[0], tempprob[1], tempprob[2], tempprob[3],probs[3]);
		
		// P(G'|offspring(s))
		probs[3] = probs[3] * target->T;
		
	    
	}
	else{
		probs[0] = target->A;
		probs[1] = target->C;
		probs[2] = target->G;
		probs[3] = target->T;
	}
	
	total = probs[0] + probs[1] + probs[2] + probs[3];
	for (newbase = 0; newbase < 4; newbase++){
		probs[newbase] = probs[newbase]/total;
	}
	prand = randum();
	tempv = 0;
	for (newbase = 0; newbase < 4; newbase++){
		tempv = tempv + probs[newbase];
		if (prand <= tempv)
			break;
	}
	
	
	//printf("%ld:, %lf, %lf, %lf, %lf\n",target->nodenum, probs[0],probs[1],probs[2],probs[3]);
	
	//printf("site %ld node %ld:, %lf, %lf, %lf, %lf, random %lf\n",site,target->nodenum, probs[0],probs[1],probs[2],probs[3],prand);
	
	
	
	target->A = target->C = target->G = target->T = 0;
	
	switch (newbase){
		case 0:
			target->A = 1;
			break;
		case 1:
			target->C = 1;
			break;
		case 2:
			target->G = 1;
			break;
		case 3:
			target->T = 1;
			break;
	}
	//if ((target->A != ancestor->A) ||(target->C != ancestor->C) ||(target->G != ancestor->G) ||(target->T != ancestor->T)){
	//  printf("new mutation site %d\n",site);
	//}
	
	
	if (target->istip == 1){ // tipnode
		olddata[target->nodenum][site] = newbase + 1;
		
	}
	
	return;
}


void simnodes_old(option_struct *op, data_fmt *data, nodeinfo *root, int site){  
	// going left
	
	simnodebase_old(op, data, root->leftnode, root, root->leftlength, site); 
	
	if (root->leftnode->istip != 1){
		simnodes_old(op, data, root->leftnode, site);
	}
	// going right
	simnodebase_old(op, data, root->rightnode, root, root->rightlength, site); 
	if (root->rightnode->istip != 1){
		simnodes_old(op, data, root->rightnode, site);
	}
	
	return;
}


void simnodes_old_fast(nodeinfo *root, int site, double *weight){  
	// going left

	int newbase;
	if (weight[0] > 0){ // going left
		weight[0] = weight[0] - root->leftlength;
		if (weight[0] < 0){ //mutation at leftnode
			
			while(1){
				newbase = (int)(4 * randum());
				if (newbase == 0){
					if (root->A == 0){
						root->leftnode->A = root->leftnode->C = root->leftnode->G = root->leftnode->T = 0;
						root->leftnode->A = 1;
						break;
					}
				}
				if (newbase == 1){
					if (root->C == 0){
						root->leftnode->A = root->leftnode->C = root->leftnode->G = root->leftnode->T = 0;
						root->leftnode->C = 1;
						break;
					}
				}
				if (newbase == 2){
					if (root->G == 0){
						root->leftnode->A = root->leftnode->C = root->leftnode->G = root->leftnode->T = 0;
						root->leftnode->G = 1;
						break;
					}
				}
				if (newbase == 3){
					if (root->T == 0){
						root->leftnode->A = root->leftnode->C = root->leftnode->G = root->leftnode->T = 0;
						root->leftnode->T = 1;
						break;
					}
				}
			}
		}
		else{ // no mutation
			root->leftnode->A = root->A;
			root->leftnode->C = root->C;
			root->leftnode->G = root->G;
			root->leftnode->T = root->T;
		}
	}
	else{
		root->leftnode->A = root->A;
		root->leftnode->C = root->C;
		root->leftnode->G = root->G;
		root->leftnode->T = root->T;
	}
    
	if (root->leftnode->istip == 1){
		if (root->leftnode->A == 1)
			olddata[root->leftnode->nodenum][site] = 1;
		if (root->leftnode->C == 1)
			olddata[root->leftnode->nodenum][site] = 2;
		if (root->leftnode->G == 1)
			olddata[root->leftnode->nodenum][site] = 3;
		if (root->leftnode->T == 1)
			olddata[root->leftnode->nodenum][site] = 4;
	}
	else{
		simnodes_old_fast(root->leftnode,site,weight);
	}
	
	// going right
	if (weight[0] > 0){ // going right
		weight[0] = weight[0] - root->rightlength;
		if (weight[0] < 0){ //mutation at rightnode
			while(1){
				newbase = (int)(4 * randum());
				if (newbase == 0){
					if (root->A == 0){
						root->rightnode->A = root->rightnode->C = root->rightnode->G = root->rightnode->T = 0;
						root->rightnode->A = 1;
						break;
					}
				}
				if (newbase == 1){
					if (root->C == 0){
						root->rightnode->A = root->rightnode->C = root->rightnode->G = root->rightnode->T = 0;
						root->rightnode->C = 1;
						break;
					}
				}
				if (newbase == 2){
					if (root->G == 0){
						root->rightnode->A = root->rightnode->C = root->rightnode->G = root->rightnode->T = 0;
						root->rightnode->G = 1;
						break;
					}
				}
				if (newbase == 3){
					if (root->T == 0){
						root->rightnode->A = root->rightnode->C = root->rightnode->G = root->rightnode->T = 0;
						root->rightnode->T = 1;
						break;
					}
				}
			}
		}
		else{ // no mutation
			root->rightnode->A = root->A;
			root->rightnode->C = root->C;
			root->rightnode->G = root->G;
			root->rightnode->T = root->T;
		}
	}
	else{
		root->rightnode->A = root->A;
		root->rightnode->C = root->C;
		root->rightnode->G = root->G;
		root->rightnode->T = root->T;
	}
    
	if (root->rightnode->istip == 1){
		if (root->rightnode->A == 1)
			olddata[root->rightnode->nodenum][site] = 1;
		if (root->rightnode->C == 1)
			olddata[root->rightnode->nodenum][site] = 2;
		if (root->rightnode->G == 1)
			olddata[root->rightnode->nodenum][site] = 3;
		if (root->rightnode->T == 1)
			olddata[root->rightnode->nodenum][site] = 4;
	}
	else{
		simnodes_old_fast(root->rightnode,site,weight);
	}	
	return;
}



void simnodes(option_struct *op, data_fmt *data, nodeinfo *root, int site){  
	// going left
	
	simnodebase(op, data, root->leftnode, root, root->leftlength, site); 
	
	if (root->leftnode->istip != 1){
		simnodes(op, data, root->leftnode, site);
	}
	// going right
	simnodebase(op, data, root->rightnode, root, root->rightlength, site); 
	if (root->rightnode->istip != 1){
		simnodes(op, data, root->rightnode, site);
	}
	
	return;
}
void printpolysite(int site){
	FILE *newhapfile;
	
	int haptype;
	newhapfile = fopen("simhap","a");
	
	fprintf(newhapfile,"%d : ",site);
	
	for (haptype = 1; haptype <=new_hap; haptype++){  
		switch(newhaps[haptype]){
			case 1:
				fprintf(newhapfile, "A");
				break;
			case 2:
				fprintf(newhapfile,"C");
				break;
			case 3:
				fprintf(newhapfile,"G");
				break;
			case 4:
				fprintf(newhapfile,"T");
				break;
		}
	}
	fprintf(newhapfile,"\n");
	fclose(newhapfile);
}

void countpoly(data_fmt *data,int site){
	int ind;
	int key = 0;
	char ****y;
	
	
	y = data->dnaptr->seqs;
	
	for (ind = 1; ind < old_hap; ind++){
		if (olddata[ind][site] != olddata[0][site]){
			key = 1;
			break;
		}
	}
	
	if (key == 1){
		for (ind = 1; ind < old_hap; ind++){
			if (y[0][0][ind][site] != y[0][0][0][site]){
				key = 2;
				break;
			}
		}
	}
	if (key == 0)
		S_N++;
	if (key == 2)
		S_O++;
	if (key == 1)
		S_H++;
	
	return;
}


void printnewsite(data_fmt *data,int site){
	int haptype, con;
	
	con = newhaps[1];
	
	for (haptype = 1; haptype <= new_hap; haptype++){
		if (newhaps[1] != newhaps[haptype]){ // polymorphic site
			printpolysite(site);
			countpoly(data, site);
			
			return;
		}
	}
	return;
} 

void simhaplotype_old(option_struct *op, data_fmt *data){ // simulate old haplotye data
	long i, site;
	char *olddna;
	int key;
	treeshape *subtree = NULL;
	double temp;
	double sumweight, rweight;
	int nums,hapi;
	double weight[1];
	
	// find mutated sites
	
	olddna = (char *)calloc(seq_length,sizeof(char));
	sumweight = 0;		
	for(i = 0; i < seq_length; i++){ 
		olddna[i] = dna_sequence[i];
		if (dna_sequence[i] != 'S'){
			sumweight = sumweight + curtree->weight_array[i];
		}
	}
	// pick #expsnpnum SNPs
	nums = 0;
	while(1){
		rweight = sumweight * randum();
		for (i = 0; i < seq_length; i++){
			if ((dna_sequence[i] != 'S') && (dna_sequence[i] != 'U')){
				rweight = rweight - curtree->weight_array[i];
				
				if (rweight <= 0){
					dna_sequence[i] = 'U';
					nums++;
					break;
				}
			}
		}
		if (nums >= expsnpnum)
			break;
	}
	
	
	
	for (site = 0; site < seq_length; site++){
		
		// construct subtree
		
		if (dna_sequence[site] == 'U'){
			
			subtree = findsubs(curtree->nodep[1],site);
			// introduce new mutations. 
			temp = dlikelihood_type(op,data, site,subtree->rootnode, 0);
			subtree->rootnode->A = subtree->rootnode->C = subtree->rootnode->G = subtree->rootnode->T = 0;
			if (olddna[site] == 'A'){
				subtree->rootnode->A = 1;
			}      
			if (olddna[site] == 'C'){
				subtree->rootnode->C = 1;
			}
			if (olddna[site] == 'G'){
				subtree->rootnode->G = 1;
			}
			if (olddna[site] == 'T'){
				subtree->rootnode->T = 1;
			}
			
			//simnodebase_old(op, data, subtree->rootnode,NULL,0,site);
			// descent to the connected nodes
			
			key = 0;
			while(key == 0){
				//simnodes_old_fast(op, data, subtree->rootnode,site, curtree->weight_array[site]); //single hit model
				//simnodes_old(op, data, subtree->rootnode,site); // likelihood model
				weight[0] = curtree->weight_array[site] * randum();
				//simnodes_old_fast(subtree->rootnode,site,weight);
				
				
				
				for (hapi = 0; hapi < old_hap; hapi++){
					if (olddata[hapi][site] != olddata[0][site]){
						key = 1;
						break;
					}
				}
				
			}
			dna_sequence[site] = 'S';
			freesubtree(subtree);
			subtree = NULL;
		}
	}	   
	
	
	
	// done !
	
	
	
}


void simhaplotype_old_2(option_struct *op, data_fmt *data){ // simulate old haplotye data over all sites
	long site;
	char *olddna;
	treeshape *subtree = NULL;

	int hapi;
	int *mafarray;
	int countp;

	double depthsum = 0;
	// find mutated sites
	
	olddna = (char *)calloc(seq_length,sizeof(char));
	mafarray = (int *)calloc(seq_length,sizeof(int));
	
	for (site = 0; site < seq_length; site++)
		olddna[site] = dna_sequence[site];
	
	// pick #expsnpnum SNPs	
	
	
	for (site = 0; site < seq_length; site++){
		
		if (dna_sequence[site] != 'S'){
			
			subtree = findsubs(curtree->nodep[1],site);
			// introduce new mutations. 
			//temp = dlikelihood_type(op,data, site,subtree->rootnode, 0);
			subtree->rootnode->A = subtree->rootnode->C = subtree->rootnode->G = subtree->rootnode->T = 0;
			if (dna_sequence[site] == 'A'){
				subtree->rootnode->A = 1;
			}      
			if (dna_sequence[site] == 'C'){
				subtree->rootnode->C = 1;
			}
			if (dna_sequence[site] == 'G'){
				subtree->rootnode->G = 1;
			}
			if (dna_sequence[site] == 'T'){
				subtree->rootnode->T = 1;
			}
			
			simnodes_old(op,data,subtree->rootnode,site);
			
			depthsum = depthsum + curtree->weight_array[site];
			
			//simnodebase_old(op, data, subtree->rootnode,NULL,0,site);
			// descent to the connected nodes
			
			freesubtree(subtree);
			subtree = NULL;
		}
	}	
	// post processing
	for (site = 0; site < seq_length; site++){
		countp = 0;
		if (dna_sequence[site] != 'S'){
			for (hapi = 0; hapi < old_hap; hapi++){
				if (olddata[hapi][site] != olddata[0][site]){
					countp++;
				}
			}
		}
		if (countp > 0){	
			if (countp > old_hap/2.0){
				mafarray[site] = old_hap - countp;
			}
			else{
				mafarray[site] = countp;
			}
		}
		else{
			mafarray[site] = 0;
		}
	}
	
	if (filtered == 1){//using fixed cutoff freq
		for (site = 0; site < seq_length; site++){
			//if (mafarray[site] >0){
			//printf("mut %d freq %d\n",site,mafarray[site]);
			//}
			if ((mafarray[site] > cutofffreq) && (dna_sequence[site] != 'S')){ // discard
				for (hapi = 0; hapi < old_hap; hapi++){
					if (olddna[site] == 'A')
						olddata[hapi][site] = 1;
					if (olddna[site] == 'C')
						olddata[hapi][site] = 2;
					if (olddna[site] == 'G')
						olddata[hapi][site] = 3;					
					if (olddna[site] == 'T')
						olddata[hapi][site] = 4;
				}
				mafarray[site] = 0;
				//printf("discarded mut site %ld %d %c freq %lf\n",site,countp,dna_sequence[site],mafarray[site]);
			}
			else{
				if (mafarray[site] > 0){
					//printf("hidden mutation : site %ld, numsel %d, marker %d \n",key, numselected, num_markers);
					//printf("hidden mututation site %ld variable %d\n",site,mafarray[site]);
					dna_sequence[site] = 'S';
				}
			}
		}
	}
	else{
		for (site = 0; site < seq_length; site++){
			if ((calling[mafarray[site]] > randum()) && (dna_sequence[site] != 'S')){ // discard
				for (hapi = 0; hapi < old_hap; hapi++){
					if (olddna[site] == 'A')
						olddata[hapi][site] = 1;
					if (olddna[site] == 'C')
						olddata[hapi][site] = 2;
					if (olddna[site] == 'G')
						olddata[hapi][site] = 3;					
					if (olddna[site] == 'T')
						olddata[hapi][site] = 4;
				}
				mafarray[site] = 0;
				//printf("discarded mut site %ld %d %c freq %lf\n",site,countp,dna_sequence[site],mafarray[site]);
			}
			else{
				if (mafarray[site] > 0){
					//printf("hidden mutation : site %ld, numsel %d, marker %d \n",key, numselected, num_markers);
					//printf("hidden mututation site %ld variable %d\n",site,mafarray[site]);
					dna_sequence[site] = 'S';
				}
			}
		}
	}
    
	/*
	 else{ // use frequency
	 // remove #num_markers new SNPs
	 numselected = 0;
	 selectarray = (int *)calloc(seq_length, sizeof(int));
	 for (site = 0; site< seq_length; site++)
	 selectarray[site] = 0;
	 
	 while(1){
	 tempv = 0;
	 key = 0;
	 for (site = 0; site < seq_length; site++)
	 tempv = tempv + mafarray[site];
	 rtempv = tempv * randum();
	 for (site = 0; site < seq_length; site++){
	 
	 rtempv = rtempv - mafarray[site];
	 if (rtempv <= 0){
	 key = site;
	 
	 }
	 if (key != 0)
	 break;
	 }
	 printf("hidden mutation : site %ld, numsel %d, marker %d \n",key, numselected, num_markers);
	 if (selectarray[key] == 0){ // discard
	 selectarray[key] = 1;
	 for (hapi = 0; hapi < old_hap; hapi++){
	 if (olddna[site] == 'A')
	 olddata[hapi][site] = 1;
	 if (olddna[site] == 'C')
	 olddata[hapi][site] = 2;
	 if (olddna[site] == 'G')
	 olddata[hapi][site] = 3;					
	 if (olddna[site] == 'T')
	 olddata[hapi][site] = 4;
	 }
	 mafarray[key] = 0;
	 numselected++;
	 }
	 if (numselected > num_markers)
	 break;
	 }
	 for (site = 0; site < seq_length; site++){
	 
	 if (mafarray[i] > 0){
	 numnewsnp++;
	 dna_sequence[site] = 'S';
	 }
	 }
	 
	 }
	 */
	
	/*
	 if (cutofffreq !=0){//using cutoff freq
	 for (site = 0; site < seq_length; site++){
	 if ((mafarray[site] > cutofffreq) && (dna_sequence[site] != 'S')){ // discard
	 for (hapi = 0; hapi < old_hap; hapi++){
	 if (olddna[site] == 'A')
	 olddata[hapi][site] = 1;
	 if (olddna[site] == 'C')
	 olddata[hapi][site] = 2;
	 if (olddna[site] == 'G')
	 olddata[hapi][site] = 3;					
	 if (olddna[site] == 'T')
	 olddata[hapi][site] = 4;
	 }
	 //printf("discarded mut site %ld %d %c freq %lf\n",site,countp,dna_sequence[site],mafarray[site]);
	 }
	 else{
	 if (mafarray[site] > 0){
	 //printf("new mut site %ld %d %c freq %lf\n",site,countp,dna_sequence[site],mafarray[site]);
	 dna_sequence[site] = 'S';
	 }
	 }
	 }
	 }
	 else{ // use frequency
	 // remove #num_markers new SNPs
	 numselected = 0;
	 selectarray = (int *)calloc(seq_length, sizeof(int));
	 for (site = 0; site< seq_length; site++)
	 selectarray[site] = 0;
	 
	 while(1){
	 tempv = 0;
	 key = 0;
	 for (site = 0; site < seq_length; site++)
	 tempv = tempv + mafarray[site];
	 rtempv = tempv * randum();
	 for (site = 0; site < seq_length; site++){
	 
	 rtempv = rtempv - mafarray[site];
	 if (rtempv <= 0){
	 key = site;
	 
	 }
	 if (key != 0)
	 break;
	 }
	 //printf("site %ld, numsel %d, marker %d \n",key, numselected, num_markers);
	 if (selectarray[key] == 0){ // discard
	 selectarray[key] = 1;
	 for (hapi = 0; hapi < old_hap; hapi++){
	 if (olddna[site] == 'A')
	 olddata[hapi][site] = 1;
	 if (olddna[site] == 'C')
	 olddata[hapi][site] = 2;
	 if (olddna[site] == 'G')
	 olddata[hapi][site] = 3;					
	 if (olddna[site] == 'T')
	 olddata[hapi][site] = 4;
	 }
	 mafarray[key] = 0;
	 numselected++;
	 }
	 if (numselected > num_markers)
	 break;
	 }
	 for (site = 0; site < seq_length; site++){
	 
	 if (mafarray[i] > 0){
	 numnewsnp++;
	 dna_sequence[site] = 'S';
	 }
	 }
	 
	 }
	 */
	
	//printf("depth \t %lf\t snp \t%d exp \t%d weight \t%lf \n",depthsum/(seq_length - num_markers - 0.0),num_markers,expsnpnum,freqw1);
	// done !
	free(olddna);
	free(mafarray);
	return;
	
}



void simhaplotype(option_struct *op, data_fmt *data){ // simulate haplotye data
	long i, site;
	
	treeshape *subtree = NULL;
	double temp,prob;
	
	
	for (site = 0; site < seq_length; site++){
		for (i = 0; i <= new_hap; i++){
			newhaps[i] = 0;
		}
		// construct subtree
		subtree = findsubs(curtree->nodep[1],site);
		
		
		if ((dna_sequence[site] != 'A') && (dna_sequence[site] != 'C') && (dna_sequence[site] !='G') && (dna_sequence[site]!='T')){
			
			temp = dlikelihood(op, data, site, subtree->rootnode);
			
			// compute likelihood of the root
		}
		else{ // introduce new mutations. 
			prob = expsnpnum / (double)(seq_length - num_markers);
			if (randum() < prob){ // polymorphic site 
				temp = dlikelihood_type(op,data, site,subtree->rootnode, 0);
				subtree->rootnode->A = subtree->rootnode->C = subtree->rootnode->G = subtree->rootnode->T = 0;
				if (dna_sequence[site] == 'A'){
					subtree->rootnode->A = 1;
				}	
				if (dna_sequence[site] == 'C'){
					subtree->rootnode->C = 1;
				}
				if (dna_sequence[site] == 'G'){
					subtree->rootnode->G = 1;
				}
				if (dna_sequence[site] == 'T'){
					subtree->rootnode->T = 1;
				}
			}
			else{
				temp = dlikelihood(op,data,site, subtree->rootnode);
			}
		}
		
		simnodebase(op, data, subtree->rootnode,NULL,0,site);
		// descent to the connected nodes
		simnodes(op, data, subtree->rootnode,site);
		
		// done !
		
		printnewsite(data,site);
		freesubtree(subtree);	
		subtree = NULL;
		
		
	}
}


void computehapprobs(nodeinfo *root, double *probarray){
	double prob;

	int numsample;
	// down left
	prob = root->leftlength;
	numsample = root->leftnum;
	probarray[numsample] = probarray[numsample] + prob;
	if (root->leftnode->istip != 1){
		computehapprobs(root->leftnode, probarray);
	}
	// down right 
	prob = root->rightlength;
	numsample = root->rightnum;
	probarray[numsample] = probarray[numsample] + prob;
	if (root->rightnode->istip != 1){
		computehapprobs(root->rightnode, probarray);
	}
	return;
	
}

int computedaughter(nodeinfo *root){
	// down left
	if (root->istip != 1){
		root->leftnum = computedaughter(root->leftnode);
		root->rightnum = computedaughter(root->rightnode);
	}
	if (root->istip == 1){
		return 1;
	}
	else{
		return root->leftnum + root->rightnum;
	}
} 

void addprobs(treeshape *result){
	int i;
	
	// for filtered option
	if (filtered != 0){
		result->probs=(double *)calloc(old_hap,sizeof(double));
		for (i = 0; i < old_hap; i++){
			result->probs[i] = 0;
		}
		computedaughter(result->rootnode);
		computehapprobs(result->rootnode, result->probs);
		for (i = 0; i < old_hap; i++){
			result->probs[0] = result->probs[0] + result->probs[i];
		}
		result->probs[0] = 1- result->probs[0];
	}
}




