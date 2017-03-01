#ifndef DROP_INCLUDE
#include "jdrop.h"
#endif
#ifdef DMEMbrDEBUG
#define MEMDEBUG
#include "memdebug.h"
#endif

#ifdef DMALLOC_FUNC_CHECK
#include "/usr/local/include/dmalloc.h"
#endif

extern double growth;

extern boolean sim_mode;
extern int old_hap, new_hap, new_tip;
extern short int **olddata;
extern double *gpartweight;
extern long **gpartrec;

extern long part1, part2,part3, part4,part5;
long HtH = 0,CtC = 0,HtC = 0,CtH = 0;
extern long sumhot[1000];

extern tree *curtree;
extern double theta0, rec0;
extern long population, locus;
extern FILE *simlog;
extern double *recarray;
extern long marker1pos,marker2pos;
extern valrec 	*tbl;
extern long seq_length;
extern boolean istip(node *p);
extern node *otherparent(node *p);
extern void fixlength(option_struct *op, data_fmt *data, node *p);
extern void ltov(option_struct *op, data_fmt *data, node *p);
extern boolean branchsub(tlist *t, node *oldbranch, node *newbranch);
extern void init_ranges_alloc(long **ranges, long numranges);
extern boolean testratio(option_struct *op, data_fmt *data, tree 
			 *oldtree, tree *newtree, char ratiotype);
extern void ranges_Malloc(node *p, boolean allokate, long numrangepairs);
extern void subtymelist(tlist *t, node *branch, boolean all);
extern double getdata_space(option_struct *op, data_fmt *data,
			    long target);
extern long countrbanches(option_struct *op, data_fmt *data,
			  tree *tr);
extern void localeval(option_struct *op, data_fmt *data, node *p,
		      boolean first);
extern void scoretree(option_struct *op, data_fmt *data, long chain);
extern double coalprob(option_struct *op, data_fmt *data, tree *tr,
		       double theta, double r);

extern long locus; /* which locus are we currently working on? */
extern long numdropped; /* how many trees are dropped due to excessive
                           amounts of recombination */
extern long indecks, apps; /* what step and chain are we on? */
/* This file contains all the functions involved in modifying the
   tree */
extern void printtreesum();
extern void scorerecs_temp(tree *target);
extern float findrate(long pos);
extern Rs *R0;
extern long startposition,endposition;
extern treeinfo **treesummary;
//test
extern struct path *rearrangestruct(recnumb *currec);
extern recnumb *scorerecs(option_struct *op, data_fmt *data, tree *target);
extern float pGS(recnumb *currec, Rs *struc);
extern struct path *copypath(path *paths);
extern void calcratio(option_struct *op, data_fmt *data, tree *oldtree, tree *newtree, long *ranges);
extern void randompathsampling_temp(double *oldrecarray, double *newrecarray, double *weight, int *recnums, long start, long end);
extern void realrandompathsampling_temp(double *oldrecarray, double *newrecarray, double *weight, int *recnums, long start, long end);
extern path *path0;

extern node *markednode;
extern double sum_prob,firstevent,firstcount;
extern long targetsite;


extern double lamda[2];
extern double recrates[2];

extern int *newgrec;
extern double *newgweight;
extern long notip;
long targettip = -99;
extern double *tiplikelist;

long countsites(option_struct *op, data_fmt *data)
{
  long nummarkers;

  nummarkers = getdata_nummarkers(op,data);

  switch(op->datatype) {
  case 'a':
    break;
  case 'b':
  case 'm':
    break;
  case 'n':
    return(data->dnaptr->sitecount[locus][nummarkers-1]);
    break;
  case 's':
    return(nummarkers);
    break;
  default:
    fprintf(ERRFILE,"countsites, found an unknown datatype %c\n",
	    op->datatype);
    break;
  }

  return(FLAGLONG);

} /* countsites */


void rename_branch(tlist *start, node *newbranch, node *oldbranch)
{
  tlist *t;
  boolean found;

  if (!newbranch) return;

  found = TRUE;
  for(t = start; found && t != NULL; t = t->succ)
    found = branchsub(t,oldbranch,newbranch);

} /* rename_branch */









void remove_branch(tlist *start, node *target)
{
  long i, j;
  node *p;
  tlist *t;
  boolean found;

  if (!target) return;

  if(!target->top) p = target->back;
  else p = target;

  for(t = start, found = TRUE; t != NULL && found; t = t->succ) {
    found = FALSE;
    for(i = 0; i < t->numbranch && !found; i++) {
      if(t->branchlist[i] == p) {
	t->numbranch--;
	for(j = i; j < t->numbranch; j++)
	  t->branchlist[j] = t->branchlist[j+1];
	found = TRUE;
      }
    }
  }

} /* remove_branch */


/* OBSOLETE */
/*********************************************************************
 * subtymenode removes the tymenode containing the node "p" from the *
 * tymelist.  The branch pointed to by the nodelet "p" will be       *
 * completely erased from the branchlist.                            */
void subtymenode(tree *tr, node *p)
{
  tlist *t;
  node *q, *r;

  t = gettymenode(tr,p->number);

  q = findunique(p);
  for(r = q->next; r == p; r = r->next) ;
  if (isrecomb(p)) rename_branch(t,q->back,r);
  else rename_branch(t,r->back,q);
  remove_branch(t,p);

} /* subtymenode */


/********************************************************************
 * fix_treenodep() needs to called everytime a node is removed from *
 * the tree.  "Number" is the number of the node to be removed.     *
 * fix_treenodep assumes that the field  containing the number of   *
 * coalescences is correct.                                         */
void fix_treenodep(tree *tr, long number)
{
  long i, numnodes;

  numnodes = 2 * tr->numcoals + 2; /* strange but true! */
  for(i = number; i < numnodes-1; i++) tr->nodep[i] = tr->nodep[i + 1];

} /* fix_treenodep */


boolean isdead(option_struct *op, node *p)
{

  if (op->fc) return(!p->coal[0]);
  else return(!p->ranges[0]);

} /* isdead */


node *findpaired_coalnode(option_struct *op, node *p)
{

  if (istip(p)) {
    fprintf(ERRFILE,"ERROR:findpaired_coalnode found a tip!\n");
    exit(-1);
  }

  if (!isrecomb(p)) return(p);

  printf("ERROR:Shouldn't be here in findpaired_coalnode!\n");
  if (isdead(op,p->next)) findpaired_coalnode(op,p->next->back);
  else if (isdead(op,p->next->next))
    findpaired_coalnode(op,p->next->next->back);
  else {fprintf(ERRFILE,"ERROR:findpaired_coalnode failed!\n"); exit(-1);}

  fprintf(ERRFILE,"ERROR:findpaired_coalnode failure.  %ld %ld\n",indecks,apps);
  return(NULL);

} /* findpaired_coalnode */


/**********************************************************************
 * remove_pairedcoalnode removes node "p" from tree "tr", and assumes *
 * that the branch to be removed is the one which has it's bottom at  *
 * the specific nodelet pointed to by "p".                            */
void remove_pairedcoalnode(option_struct *op, data_fmt *data,
			   tree *tr, node *p)
{
  tlist *t;
  node *q;

  tr->numcoals--;

  hookup(p->next->next->back,p->next->back);
  fixlength(op,data,p->next->back);
  fix_treenodep(tr,p->number);

  t = gettymenode(tr,p->number);
  q = (p->next->top) ? p->next->next->back : p->next->back;
  subtymelist(t,q,FALSE);
  freenode(p);
} /* remove_pairedcoalnode */


/****************************************************************
 * remove_pairedrecombnode removes node "p" from tree "tr", and *
 * assumes that "p" points to the unique nodelet of the         *
 * recombination node.                                          */
void remove_pairedrecombnode(option_struct *op, data_fmt *data,
			     tree *tr, node *p, node *dead)
{
  node *q;
  tlist *t;

  tr->numrecombs--;

  for(q = p->next; q == dead; q = q->next) ;
  hookup(p->back,q->back);
  fixlength(op,data,p->back);
  fix_treenodep(tr,p->number);

  t = gettymenode(tr,p->number);
  subtymelist(t,dead,FALSE);
  freenode(p);
} /* remove_pairedrecombnode */


void remove_pairednodes(option_struct *op, data_fmt *data, tree *tr,
			node *p, node *remove_thisway)
{
  node *q;

  if (!remove_thisway->top) {
    fprintf(ERRFILE,"ERROR:remove_pairednodes passed a non-top!\n");
    exit(-1);
  }
  q = findpaired_coalnode(op,remove_thisway->back);
  remove_pairedcoalnode(op,data,tr,q);
  remove_pairedrecombnode(op,data,tr,p,remove_thisway);

} /* remove_pairednodes */


long setnodenumber(boolean **nodenumber, long *numnodes)
{
  long i;
  //printf("node %ld\n",*numnodes);
  for (i = 0; i < (*numnodes); i++)
    if (!(*nodenumber)[i]) {
      (*nodenumber)[i] = TRUE;
      return(i);
    }

  (*numnodes)++;
  *nodenumber = (boolean *)realloc(*nodenumber,(*numnodes)*sizeof(boolean));
  (*nodenumber)[(*numnodes) - 1] = TRUE;

  curtree->nodep = (node **)
    realloc(curtree->nodep,((*numnodes)+1)*sizeof(node *));

  return((*numnodes) - 1);

} /*setnodenumber */

boolean foundbranch(tlist *t, node *p)
{
  long i;

  for(i = 0; i < t->numbranch; i++) 
    if (t->branchlist[i] == p) return(TRUE);

  return(FALSE);

} /* foundbranch */

void readdactive(tlist *stop, node *p)
{
  long i, j;
  tlist *t;

  t = gettymenode(curtree,p->number);
  do {
    if (!foundbranch(t,p)) {
      t->numbranch++;
      t->branchlist = 
	(node **)realloc(t->branchlist, t->numbranch*sizeof(node *));
      t->branchlist[t->numbranch-1] = p;
    }
    t = t->succ;
  } while (t != stop && t != NULL);

  /* remove the FALSE branchlist entry from the rest of the tymelist */
  while (t != NULL) {
    if (!foundbranch(t,p)) break;
    for (i = 0; i < t->numbranch; i++)
      if (t->branchlist[i] == p) break;
    for (j = i; j < t->numbranch - 1; j++)
      t->branchlist[j] = t->branchlist[j+1];
    t->numbranch--;
    t = t->succ;
  }

} /* readdactive */

void newlin(linlist **t)
{
  *t = (linlist *)(calloc(1,sizeof(linlist)));
  //*t = (linlist *)(malloc(sizeof(linlist)));
  (*t)->prev = NULL;
  (*t)->succ = NULL;
} /* newlin */

void freelin(linlist *t)
{
  free(t);
} /* freelin */


void freelinlist(linlist *t)
{
  linlist *u, *usucc;

  for(u = t; u != NULL; u = usucc) {usucc = u->succ; freelin(u);}

} /* freelinlist */


void addlinlist(linlist **t, node *target, double activelinks)
{
  linlist *r, *s, *u;

  s = *t;
  if (s == NULL) {
    newlin(&s);  
    s->branchtop = target;
    s->activesites = activelinks;
    s->succ = NULL;
    s->prev = NULL;
    (*t)=s;
  } else {
    newlin(&r);
    r->branchtop = target;
    r->activesites = activelinks;
    /* insert new entry right after s */
    u = s->succ;
    s->succ = r;
    r->prev = s;
    r->succ = u;
    if(u!=NULL) u->prev = r;
  }
} /* addlinlist */


void sublinlist(linlist **t, node *target)
{
  boolean found;
  linlist *u;

  u = *t;
  found = FALSE;
  while (!found) {
    if (u->branchtop==target) {
      found = TRUE;
      if (u->succ != NULL) u->succ->prev = u->prev;
      if (u->prev != NULL) u->prev->succ = u->succ;
      /* if this was the first entry, reset the lineage pointer */
      if (u==(*t)) *t = u->succ; 
      freelin(u);
    } else if (u->succ != NULL) u = u->succ;
    else fprintf(ERRFILE, "ERROR:SUBLINLIST failed to find target\n");
  }
} /* sublinlist */


/******************************************************
 * printlinlist() prints a linlist:  a debug function */
void printlinlist(linlist *u)
{
  if (u) {
    do {
      if(u->branchtop!=NULL) 
	fprintf(ERRFILE,"Branchtop %ld activesites %g\n",u->branchtop->number,
		u->activesites);
      else fprintf(ERRFILE,"This entry has no branchtop!\n");
      u=u->succ;
    } while (u!=NULL);
  }

}  /* printlinlist */


/***************************************************************
 * foundlin() finds whether a given nodelet is already present *
 * in the passed linlist.  TRUE = found, FALSE = not found.    */
boolean foundlin(linlist *u, node *p)
{
  for(;u != NULL; u = u->succ)
    if(u->branchtop == p) return (TRUE);

  return(FALSE);

} /* foundlin */


void findlin(linlist **u, long which)
{
  long i;

  i = 1;
  do {
    if (i==which) return;
    (*u) = (*u)->succ;
    i++;
  } while ((*u)!=NULL);
  printf("ERROR:Unable to find linlist entry!\n");
} /* findlin */


/************************************************************************
 * getnumlinks() returns the number of links between the passed "sites" */
double getnumlinks(option_struct *op, data_fmt *data, long start, long end)
{
  long i;
  double sum, *spaces;

  spaces = NULL;

  switch(op->datatype) {
  case 'a':
    break;
  case 'b':
  case 'm':
    spaces = data->msptr->mspace[population][locus];
    break;
  case 'n':
    return((double)(end - start));
    break;
  case 's':
    spaces = data->dnaptr->sspace[population][locus];
    for(i = start, sum = 0.0; i < end; i++)
      sum += spaces[i];
    return(sum);
    break;
  default:
    fprintf(ERRFILE,"unknown datatype in getnumlinks\n");
    break;
  }

  return((double)FLAGLONG);

} /* getnumlinks */


/*******************************************************************
 * count_active() counts the number of active links leaving a node */
double count_active(option_struct *op, data_fmt *data, node *p)
{
  double value;
  node *q;

  q = p;
  if (!q->top) q = q->back;

  if (istip(q)) 
    return(getnumlinks(op,data,0L,countsites(op,data)-1));

  if (isrecomb(q)) {
    value = count_rec_active(op,data,q);
  }
  else {
    value = count_coal_active(op,data,q);
  }

  return(value);

} /* count_active */


/*********************************************************************
 * count_coal_active() returns the number of active links along the  *
 * branch "p", where the node containing "p" must be coalescent; and *
 * the branch "p" may not exist yet.                                 */
double count_coal_active(option_struct *op, data_fmt *data, node *p)
{
  long *range1, *range2, newstart, newend;
      
  range1 = p->next->back->ranges;
  range2 = p->next->next->back->ranges;

  if (!range1[0] && !range2[0])
    fprintf(ERRFILE, "ERROR:Dead branch in count_coal_active! %ld %ld\n",
	    indecks,apps);
 
  if (!range1[0]) return(range2[2*range2[0]] - range2[1]);
  if (!range2[0]) return(range1[2*range1[0]] - range1[1]);

  newstart = (range1[1] < range2[1]) ? range1[1] : range2[1];
  newend = (range1[2*range1[0]] > range2[2*range2[0]]) ?
    range1[2*range1[0]] : range2[2*range2[0]];

  if (newend - newstart < 0)
    fprintf(ERRFILE,"ERROR:negative link count in count_coal_active! %ld %ld\n",
	    indecks,apps);

  return(getnumlinks(op,data,newstart,newend));

} /* count_coal_active */


/*******************************************************************
 * count_rec_active() returns the number of active links along the *
 * branch "p", where the node containing "p" must be recombinant;  *
 * and the branch "p" may not exist yet.                           */
double count_rec_active(option_struct *op, data_fmt *data, node *p)
{
  long i, newstart, newend;
  node *q;

  q = findunique(p)->back;

  if (!q->ranges[0]);
  //fprintf(ERRFILE,"ERROR:count_rec_active counted dead 3, %ld %ld\n",indecks,apps);

  for(i = 1, newstart = -1L; q->ranges[i] != FLAGLONG; i+=2) {
    if(p->recstart > q->ranges[i+1]) continue;
    newstart = (q->ranges[i] > p->recstart) ? q->ranges[i] : p->recstart;
    break;
  }

  if(newstart == -1L) return(0L);
  for(i=2*q->ranges[0], newend = -1L; i != 0; i-=2) {
    if(p->recend < q->ranges[i-1]) continue;
    newend = (q->ranges[i] < p->recend) ? q->ranges[i] : p->recend;
    break;
  }

  if(newend == -1L) return(0L);

  if (newend - newstart < 0)
    fprintf(ERRFILE,"ERROR:negative link count in count_rec_active! %ld %ld\n",
	    indecks,apps);

  return(getnumlinks(op,data,newstart,newend));

} /* count_rec_active */


/*********************************************************************
 * count_activefc() counts the number of active links leaving a node *
 * taking into account which sites have achieved final coalescence.  */
double count_activefc(option_struct *op, data_fmt *data, node *p)
{
  double value;
  node *q;

  q = p;
  if (!q->top) q = q->back;

  if (istip(q))
    return(getnumlinks(op,data,0L,countsites(op,data)-1));

  if (op->datatype == 'n')
    value = ((q->coal[0]) ? q->coal[2*q->coal[0]] - q->coal[1] : 0L);
  else 
    value = ((q->coal[0]) ? 
	     getnumlinks(op,data,q->coal[1],q->coal[2*q->coal[0]]) : 0L);

  return(value);

} /* count_activefc */


/******************************************************************
 * is_fc() returns TRUE if site has achieved fc by treesec, FALSE *
 * otherwise.                                                     */
boolean is_fc(tlist *treesec, long site)
{
  long i, count;

  for(i = 0, count = 0; i < treesec->numbranch; i++) {
    if(inrange(treesec->branchlist[i]->ranges,site)) count++;
    if (count > 1) return(FALSE);
  }

  return(TRUE);

} /* is_fc */


/***********************************************************
 * fix_coal() upodates the coal array for the passed node. *
 * "tfc" is the tymeslice used for FC calculations and all *
 * tipwards tymeslices plus their eventnode's coal arrays  *
 * are assumed to be correct.                              */
void fix_coal(option_struct *op, data_fmt *data, tree *tr, tlist *tfc,
	      node *p)
{
  long *subtrees;
  node *q, *r, *s;

  if (isrecomb(p)) {
    q = findunique(p);
    r = q->next;
    s = q->next->next;
    q = q->back;
    copycoal(q,r);
    copycoal(q,s);
    subrangefc(&(r->coal),s->recstart,s->recend);
    subrangefc(&(s->coal),r->recstart,r->recend);
  } else {
    subtrees = NULL;
    q = findunique(p);
    findsubtrees_node(op,data,tfc,tr,&subtrees);
    makecoal(op,data,tr,tfc,&(q->coal),q,subtrees);
    free(subtrees);
  }

} /* fix_coal */


/******************************************************************
 * makecoal() creates the appropiate coal array for the passed in *
 * tymeslice.                                                     *
 * makecoal() is nodelet specific, and the passed in nodelet must *
 * have correct information contained in its ranges array.        *
 * makecoal() is nodelet specific iff "p" is passed in as non-NULL*/
void makecoal(option_struct *op, data_fmt *data, tree *tr,
	      tlist *tfc, long **coal, node *p, long *subtrees)
{
  long i;

  init_coal_alloc(coal,1L);
  (*coal)[1] = 0;
  (*coal)[2] = countsites(op,data)-1;

  for(i = 0; i < subtrees[0]; i++) {
    if (p)
      if (!inrange(p->ranges,subtrees[2*i+1])) {
	subrangefc(coal,subtrees[2*i+1],subtrees[2*i+2]);
	continue;
      }

    if (is_fc(tfc,subtrees[2*i+1]))
      subrangefc(coal,subtrees[2*i+1],subtrees[2*i+2]);

  }

} /* makecoal */


/* three functions involving updating the alias array "siteptr" */

void edit_alias(option_struct *op, data_fmt *data, long *sp, long cutpoint)
{
  long i, alias_site, cutmarker, numsites;

  /* this looks awful because sp[i] does not actually contain
     the alias site, it contains the alias site + 1 (to avoid zero). */

  /* cutpoint is the first site *after* the recombination */

  numsites = getdata_nummarkers(op,data);
  cutmarker = sitetorightmarker(op,data,cutpoint);
  if (cutmarker == FLAGLONG) return;

  for(i=cutmarker;i<numsites;i++) {
    alias_site = sp[i]-1;
    if(alias_site < cutmarker && alias_site+1 > 0) sp[i]*= -1;
  }
} /* edit_alias */


void traverse_rebuild_alias(option_struct *op, data_fmt *data, 
			    tlist *tstart, long *sp)
{
  tlist *t;
  node *p;

  for(t = tstart; t != NULL; t = t->succ) {
    p = t->eventnode;
    if (!isrecomb(p)) continue;
    if (p->recstart == 0) edit_alias(op,data,sp,p->recend+1);
    else edit_alias(op,data,sp,p->recstart);
  }

} /* traverse_rebuild_alias */


void rebuild_alias(option_struct *op, data_fmt *data, long *sp)
{
  long i, numsites;
  
  numsites = getdata_nummarkers(op,data);

  for(i=0;i<numsites;i++) if (sp[i]<0) sp[i]*= -1;
  traverse_rebuild_alias(op,data,curtree->tymelist,sp);
} /* rebuild_alias */


/****************************************************************
 * contrib() sets the passed in "ranges" field, using "p".  "p" *
 * points to the nodelet at the top of the relevant branch.     */
void contrib(option_struct *op, data_fmt *data, node *p, long **newranges)
{
  node *q, *r;
  long i, numremove;


  if(istip(p)) {
    (*newranges)[0] = 1L;
    (*newranges)[1] = 0L;
    (*newranges)[2] = countsites(op,data)-1;
    (*newranges)[3] = FLAGLONG;
    return;
  }

  if(isrecomb(p)) {
    q = findunique(p)->back;
    init_ranges_alloc(newranges,q->ranges[0]);
    memcpy((*newranges),q->ranges,(2*q->ranges[0]+2)*sizeof(long));
    /* first remove the leading excess ranges */
    for(i = 1, numremove = 0; (*newranges)[i] != FLAGLONG; i+=2) {
      if(p->recstart > (*newranges)[i+1]) {numremove+=2; continue;}
      if(p->recstart > (*newranges)[i]) (*newranges)[i] = p->recstart;
      break;
    }
    memmove(&(*newranges)[1],&(*newranges)[1+numremove],
	    (2*(*newranges)[0]-numremove+1)*sizeof(long));
    (*newranges)[0] -= numremove/2;
    /* then removing the trailing excess ranges */
    for(i = 2*(*newranges)[0]-1, numremove = 0; i+1 != 0; i-=2) {
      if(p->recend < (*newranges)[i]) {numremove++; continue;}
      if(p->recend < (*newranges)[i+1]) (*newranges)[i+1] = p->recend;
      break;
    }
    (*newranges)[i+2] = FLAGLONG;
    (*newranges)[0] -= numremove;
  } else {
    q = p->next->back;
    r = p->next->next->back;
    init_ranges_alloc(newranges,q->ranges[0]);
    memcpy((*newranges),q->ranges,(2*q->ranges[0]+2)*sizeof(long));
    for(i = 1; r->ranges[i] != FLAGLONG; i+=2)
      addrange(newranges,r->ranges[i],r->ranges[i+1]);
  }

} /* contrib */


/***************************************************************
 * fix_coal_ranges() updates the "ranges" of a coalescent node */
void fix_coal_ranges(option_struct *op, data_fmt *data, node *p)
{
  p = findunique(p);
  contrib(op,data,p,&p->ranges);
} /* fix_coal_ranges */


/***************************************************************
 * fix_rec_ranges() updates the "ranges" of a recombinant node */
void fix_rec_ranges(option_struct *op, data_fmt *data, node *p)
{
  p = findunique(p);
  contrib(op,data,p->next,&p->next->ranges);
  contrib(op,data,p->next->next,&p->next->next->ranges);
} /* fix_rec_ranges */


boolean popbranch(tlist *t, long branch)
{
  long i;
  boolean found;

  found = FALSE;
  for(i=0;i<t->numbranch;i++) {
    if (t->branchlist[i]!=NULL) {
      if (t->branchlist[i]->number==branch) {
        t->branchlist[i]=NULL;
        found = TRUE;
      }
    }
  }
  return(found);
} /* popbranch */


boolean subbranch(tlist *t, long oldbranch, node *newbranch)
{
  long i;
  boolean found;

  found = FALSE;
  for(i=0;i<t->numbranch;i++) {
    if (t->branchlist[i]!=NULL) {
      if (t->branchlist[i]->number==oldbranch) {
        t->branchlist[i]=newbranch;
        found = TRUE;
      }
    }
  }
  return(found);
} /* subbranch */


void poptymelist (node *p)
{  /* remove futile tymelist entry */
  boolean done;
  node *q;
  tlist *s, *t, *u;

  t = gettymenode(curtree, p->number); 
  u = t;
  if(isrecomb(p)) {
    s = u->succ;
    do {
      if(s==NULL) break;
      done = !popbranch(s,p->number);   
      s = s->succ;
    } while (!done);
    u->prev->age=t->age;
    u->prev->succ=u->succ;
    if (t->succ != NULL) u->succ->prev=u->prev;
    freetymenode(u);
  } else {
    s = u->succ;
    if(p->next->top) q=p->next->next->back;
    else q=p->next->back;
    do {
      if(s==NULL) break; 
      done = !subbranch(s,p->number,q);   
      s = s->succ;
    } while (!done);
    u->prev->age=t->age;
    u->prev->succ=u->succ;
    if (u->succ != NULL) u->succ->prev=u->prev;
    freetymenode(u);
  }
} /* poptymelist */
   

void fixbranch(tlist *t)
{
  long i, j, nullcnt, newnum;
  node **newarray;

  nullcnt = 0;
  for(i=0;i<t->numbranch;i++) {
    if(t->branchlist[i]==NULL) nullcnt++;
  }
  if(nullcnt==0) return;
  /* workaround to avoid callocing 0 (or less) elements */
  newnum = (t->numbranch-nullcnt < 1) ? 1 : t->numbranch-nullcnt;
  newarray = (node **)calloc(1,newnum*sizeof(node *));
  j = 0;
  for(i=0;i<t->numbranch;i++) {
    if(t->branchlist[i]!=NULL) {
      newarray[j]=t->branchlist[i];
      j++;
    }
  }
  free(t->branchlist);
  t->numbranch = newnum;
  t->branchlist=newarray;
} /* fixbranch */


/***********************************************************
 * is_futile() checks to see if a node has been tagged for *
 * futile removal.                                         */
boolean is_futile(node *p)
{

  return((p->futileflag == FLAGLONG));

} /* is_futile */


/**********************************************************
 * find_nonfutile() finds the first non-futile node along *
 * the path of nodes going in the direction indicated by  *
 * upwards (TRUE = up, FALSE = down).                     */
node *find_nonfutile(node *p, node *cutnode, boolean upwards)
{
  node *q;

  if (p == cutnode) return(NULL);
  if (istip(p) || !is_futile(p)) return(p);

  p = findunique(p);
  if (isrecomb(p)) {
    if (upwards) q = find_nonfutile(p->back,cutnode,upwards);
    else {
      q = find_nonfutile(p->next->back,cutnode,upwards);
      if (!q) q = find_nonfutile(p->next->next->back,cutnode,upwards);
    }
  } else {
    if (upwards) {
      q = find_nonfutile(p->next->back,cutnode,upwards);
      if (!q) q = find_nonfutile(p->next->next->back,cutnode,upwards);
    } else q = find_nonfutile(p->back,cutnode,upwards);
  }

  return(q);

} /* find_nonfutile */


/**********************************************************
 * tag_futile() sets a flag to indicate that it should be *
 * removed.  The flag set is for recstart & recend set to *
 * FLAGLONG;                                              */
void tag_futile(node *p)
{

  if (isrecomb(p)) {
    p->futileflag = p->next->futileflag = FLAGLONG;
    p->next->next->futileflag = FLAGLONG;
    tag_futile(p->next->back);
    if (p->next->next->back != curtree->root)
      tag_futile(p->next->next->back);
    return;
  } else {
    if (!is_futile(p)) {
      p->futileflag = p->next->futileflag = FLAGLONG;
      p->next->next->futileflag = FLAGLONG;
      return;
    } else {
      if (p->next->top) tag_futile(p->next->back);
      else tag_futile(p->next->next->back);
    }
  }

} /* tag_futile */


/**********************************************************************
 * renamebrlistb() goes through a brlist replacing old branchtops with *
 * the new branchtop.                                                 */
void renamebrlist(brlist **bigbr, brlist *start, node *oldtop, node *newtop)
{
  brlist *br, *brsucc;

  for(br = start; br != NULL; br = brsucc) {
    brsucc = br->succ;
    if (br->branchtop == oldtop) {
      br->branchtop = newtop;
      if (newtop->tyme >= br->endtyme) {
	br_remove(bigbr,br);
	continue;
      }
      if (br->starttyme < newtop->tyme) {
	br->starttyme = newtop->tyme;
	if (br->prev) br->prev->succ = br->succ;
	else (*bigbr) = br->succ;
	if (br->succ) br->succ->prev = br->prev;
	hookup_brlist(bigbr,br);
      }
    }
  }

} /* renamebrlist */


/**********************************************************************
 * drop_renamebrlist() is a coalesce specific driver for renamebrlist */
void drop_renamebrlist(brlist **bigbr, node *oldtop, node *newtop)
{
  brlist *br;

  if (!bigbr) return;

  br = (*bigbr);

  renamebrlist(bigbr,br,oldtop,newtop);

} /* drop_renamebrlist */


/********************************************************************
 * drop_brfix() updates the ranges & coal information, and changes  *
 * brlist info, for all branches which start at or after the passed *
 * in tymelist's eventnode.                                         */
void drop_brfix(option_struct *op, data_fmt *data, tree *newtree,
		tlist *start, brlist **oldbr)
{
  node *p, *q, *r;
  tlist *t;
  long *oldcoal1, *oldcoal2;
  dnadata *dna;

  oldcoal1 = oldcoal2 = NULL;
  dna = data->dnaptr;
 
  for (t = start; t != NULL; t = t->succ) {
    p = findunique(t->eventnode);
    if (istip(p)) continue;
    if (isrecomb(p)) {
      q = p->next;
      r = p->next-> next;
      if (op->fc) {
	fix_rec_ranges(op,data,p);
	init_coal_alloc(&oldcoal1,q->coal[0]);
	init_coal_alloc(&oldcoal2,r->coal[0]);
	memcpy(oldcoal1,q->coal,(2*q->coal[0]+2)*sizeof(long));
	memcpy(oldcoal2,r->coal,(2*r->coal[0]+2)*sizeof(long));
	fix_coal(op,data,newtree,t,p);
	if (!sameranges(oldcoal1,q->coal))
	  addbrlistH(op,data,oldbr,q,oldcoal1,q->coal,NULL,q->tyme,q->back->tyme,FALSE);
	if (!sameranges(oldcoal2,r->coal))
	  addbrlistH(op,data,oldbr,r,oldcoal2,r->coal,NULL,r->tyme,r->back->tyme,FALSE);
      } else {
	init_ranges_alloc(&oldcoal1,q->ranges[0]);
	init_ranges_alloc(&oldcoal2,r->ranges[0]);
	memcpy(oldcoal1,q->ranges,(2*q->ranges[0]+2)*sizeof(long));
	memcpy(oldcoal2,r->ranges,(2*r->ranges[0]+2)*sizeof(long));
	fix_rec_ranges(op,data,p);
	if (!sameranges(oldcoal1,q->ranges))
	  addbrlistH(op,data,oldbr,q,oldcoal1,q->ranges,NULL,q->tyme,
		     q->back->tyme,FALSE);
	if (!sameranges(oldcoal2,r->ranges))
	  addbrlistH(op,data,oldbr,r,oldcoal2,r->ranges,NULL,r->tyme,
		     r->back->tyme,FALSE);
      }
    } else {
      if (op->fc) {
	fix_coal_ranges(op,data,p);
	init_coal_alloc(&oldcoal1,p->coal[0]);
	memcpy(oldcoal1,p->coal,(2*p->coal[0]+2)*sizeof(long));
	fix_coal(op,data,newtree,t,p);
	if (!sameranges(oldcoal1,p->coal))
	  addbrlistH(op,data,oldbr,p,oldcoal1,p->coal,NULL,
		     p->tyme,p->back->tyme,FALSE);
      } else {
	init_ranges_alloc(&oldcoal1,p->ranges[0]);
	memcpy(oldcoal1,p->ranges,(2*p->ranges[0]+2)*sizeof(long));
	fix_coal_ranges(op,data,p);
	if (!sameranges(oldcoal1,p->ranges))
	  addbrlistH(op,data,oldbr,p,oldcoal1,p->ranges,NULL,
		     p->tyme,p->back->tyme,FALSE);
      }
    }
  }

  if (oldcoal1) free(oldcoal1);
  if (oldcoal2) free(oldcoal2);

} /* drop_brfix */


/************************************************************
 * samesites() checks to see if the segranges array and the *
 * ranges array could be describing the same set of sites.  */
boolean samesites(int *segranges, long *ranges, long numsites)
{
  long i;

  for(i = 0; i < numsites; i++) {
    if(segranges[i] == 1 || segranges[i] == 2) {
      if(!inrange(ranges,i)) return(FALSE);
    } else {
      if(inrange(ranges,i)) return(FALSE);
    }
  }

  return(TRUE);

} /* samesites */


/****************************************************************
 * treelegalbrlist() checks to make sure that a brlist entry is *
 * legal given the tree.                                        */
void treelegalbrlist(option_struct *op, data_fmt *data, tree *tr,
		     brlist **bigbr, brlist *target)
{
  long *comp;
  node *top, *bottom, *treenode;

  top = target->branchtop;
  bottom = top->back;
  treenode = tr->nodep[top->number];

  if (isrecomb(top)) {
    if(treenode != top && treenode != otherparent(top)) {
      br_remove(bigbr,target);
      return;
    }
  } else {
    if(treenode != top) {
      br_remove(bigbr,target);
      return;
    }
  }

  if (target->endtyme > bottom->tyme) {
    if (isrecomb(bottom)) {
      bottom = findunique(bottom);
      if (!bottom->next->back) top = bottom->next->next;
      else if (!bottom->next->next->back) top = bottom->next;
      else {
	comp = ((op->fc) ? bottom->next->coal : bottom->next->ranges);
	top = 
	  (samesites(target->segranges,comp,countsites(op,data))
	   ? bottom->next : bottom->next->next);
      }
    } else top = findunique(bottom);

    if (target->endtyme < curtree->root->back->tyme)
      addbrlistH(op,data,bigbr,top,NULL,NULL,target->segranges,bottom->tyme,
		 target->endtyme,FALSE);
    target->endtyme = bottom->tyme;
  }

} /* treelegalbrlist */


/******************************************************************
 * samebranch() checks to see if the two segments are on the same *
 * branch (i.e. whether they have the same tops).                 */
boolean samebranch(brlist *br1, brlist *br2)
{

  return(br1->branchtop == br2->branchtop);

} /* samebranch */


/********************************************************************
 * checktreebranch() checks to see if a segment is properly part of *
 * the branch whose top is p.                                       *
 *                                                                  *
 * It returns: 0 if no match is found                               *
 *             1 if starts match                                    *
 *            -1 if ends match                                      *           
 *             2 if both start and end match                        */
long checktreebranch(brlist *br, node *p, double starttyme,
		     double endtyme)
{

  if (br->branchtop != p) return(0L);

  if (br->starttyme == starttyme && br->endtyme == endtyme) return(2L);
  if (br->starttyme == starttyme) return(1L);
  if (br->endtyme == endtyme) return(-1L);

  return(0L);

} /* checktreebranch */


/********************************************************************
 * sametreebranch() checks to see if a segment was ever on the same *
 * branch as the branch possibly bisected by the node p.            *
 *                                                                  *
 * It returns: 0 if no match is found                               *
 *             1 if starts match                                    *
 *            -1 if ends match                                      *           
 *             2 if both start and end match                        */
long sametreebranch(brlist *br, node *p, double starttyme,
		    double endtyme)
{
  node *target;
  long match;

  match = checktreebranch(br,p,starttyme,endtyme);

  if (match /*|| istip(p)*/) return(match);

  target = findunique(p);
  if (isrecomb(target)) {
    match = checktreebranch(br,target->back,starttyme,endtyme);
  } else {
    match = checktreebranch(br,target->next->back,starttyme,endtyme);
    if (match) return(match);
    match = checktreebranch(br,target->next->next->back,starttyme,endtyme);
  }

  return(match);

} /* sametreebranch */


/****************************************************************
 * isoverlap() checks to see if two brlist entries overlap each *
 * other on the tree.                                           */
boolean isoverlap(brlist *br1, brlist *br2)
{
  double br1start, br2start;

  if (!samebranch(br1,br2)) return(FALSE);

  br1start = br1->starttyme;
  br2start = br2->starttyme;

  if (((br1start > br2start) && (br1start < br2->endtyme)) ||
      ((br2start > br1start) && (br2start < br1->endtyme))) return(TRUE);

  return(FALSE);


} /* isoverlap */


/******************************************************************
 * iscontainedin() checks to see if one brlist entry is entirely  *
 * within another.                                                */
boolean iscontainedin(brlist *br1, brlist *br2)
{

  if (!samebranch(br1,br2)) return(FALSE);

  if(((br1->starttyme > br2->starttyme) && (br1->endtyme < br2->endtyme))
     || ((br2->starttyme > br1->starttyme) && (br2->endtyme < br1->endtyme)))
    return(TRUE);

  return(FALSE);

} /* iscontainedin */


/****************************************************************
 * mash_brlist() searches a brlist for a match to its "target". *
 * If it finds one, it mashes the two entries together.         */
void mash_brlist(tree *tr, brlist *start, brlist *target)
{
  brlist *br, *brupper, *brlower;
  double temp;

  for(br = start; br != NULL; br = br->succ) {
    if(!isoverlap(target,br) && !iscontainedin(target,br)) continue;
    if(target->starttyme < br->starttyme) {
      brupper = target; 
      brlower = br;
    } else {
      brupper = br; 
      brlower = target;
    }
    if (iscontainedin(target,br)) {
      brupper->endtyme = brlower->starttyme;
      continue;
    }
    temp = brupper->endtyme;
    brupper->endtyme = brlower->starttyme;
    brlower->starttyme = temp;
  }

} /* mash_brlist */


/*****************************************************************
 * consolidate_brlist() attempts to remove and replace redundate *
 * entries from a brlist.  The tree must be correct.             */
void consolidate_brlist(option_struct *op, data_fmt *data, tree *tr,
			brlist **bigbr)
{
  brlist *br, *brsucc;

  if (!bigbr) return;

  for(br = (*bigbr); br != NULL; br = brsucc) {
    brsucc = br->succ;
    treelegalbrlist(op,data,tr,bigbr,br);
  }

  for(br = (*bigbr); br != NULL; br = br->succ) {
    mash_brlist(tr,br->succ,br);
  }

} /* consolidate_brlist */


/******************************************************************
 * br_remove() performs a very basic removal of the passed brlist *
 * element from whatever brlist it is attached to, then frees it. */
void br_remove(brlist **bigbr, brlist *target)
{
  if (*bigbr == NULL){
    free(target);
    return;
  }

  if (target->prev) target->prev->succ = target->succ;
  else (*bigbr) = target->succ;
  if (target->succ) target->succ->prev = target->prev;

  freebrlist(target);

} /* br_remove */


/*************************************
 * printbrlist() prints out a brlist */
void printbrlist(brlist *start)
{
  brlist *br;

  fprintf(ERRFILE,"Brlist follows\n");
  for(br = start; br != NULL; br = br->succ) {
    fprintf(ERRFILE,"   top %3ld at tyme %12.8f [%12.8f/%12.8f]",
	    br->branchtop->number,br->branchtop->tyme,br->starttyme, br->endtyme);
    fprintf(ERRFILE," with newlinks %f\n",br->numnewactives);
  }

} /* printbrlist */


/*********************************************************************
 * find_futilecoaldtr() returns the "traversed" daughter node from a *
 * the passed node "p".                                              */
node *find_futilecoaldtr(node *cutnode, node *p)
{
  node *q;

  q = findtop(p);

  if (is_futile(q->next->back) || q->next->back == cutnode)
    return(q->next->next->back);
  else return(q->next->back);

} /* find_futilecoaldtr */


/********************************************************************
 * remove_futile() removes the nodes from the tree that have been   *
 * tagged by tag_futile, while maintaing brlist entries for them.   *
 * remove_futile() assumes that the first entry in the tymelist     *
 * contains the tips and so will never be removed.  It also assumes *
 * that coalescent nodes require tree clean-up, whereas recombinant *
 * nodes may be simply removed.                                     *
 * Traverse the tymelist in reverse order so that path of removal   *
 * info is preserved.                                               */
void remove_futile(option_struct *op, data_fmt *data, tree *newtree,
		   tree *oldtree, node *cutnode, boolean *nodenumber, brlist **br)
{
  long i, numnodes, **oranges, *oldranges, topcount, *actives;
  double btyme;
  tlist *t, *tsucc, *firstfutile;
  node *p, *q, *dtr, *par, **obnode, **otnode, *oldnode, *newtopnode,
    **tops;
  dnadata *dna;

  dna = data->dnaptr;

  numnodes = 2 * oldtree->numcoals + 2;
  obnode = (node **)calloc(numnodes,sizeof(node *));
  otnode = (node **)calloc(numnodes,sizeof(node *));
  tops = (node **)calloc(numnodes,sizeof(node *));
  oranges = (long **)calloc(numnodes,sizeof(long *));

  topcount = 0;
  firstfutile = NULL;

  /* scanning tymelist from the bottom up */
  for(t = newtree->tymelist; t->succ != NULL; t = t->succ) ;

  /* modify the basic underlying tree structure */
  for(; t->prev != NULL; t = t->prev) {
    obnode[t->eventnode->number] = otnode[t->eventnode->number] = NULL;
    if (!is_futile(t->eventnode)) continue;
    firstfutile = t;
    p = t->eventnode;
    nodenumber[p->number] = FALSE;
    if (isrecomb(p)) {
      newtree->numrecombs--;
      q = findunique(p);
      otnode[p->number] = (is_futile(q->next->back)) ? 
	q->next : q->next->next;
    } else {
      obnode[p->number] = p->back;
      dtr = find_nonfutile(p,cutnode,TRUE);
      otnode[p->number] = dtr;
      if (dtr != NULL) {
        tops[topcount] = dtr;
        topcount++;
      }
      actives = (op->fc) ? p->coal : p->ranges;
      init_coal_alloc(&oranges[p->number],actives[0]);
      memcpy(oranges[p->number],actives,(2*actives[0]+2)*sizeof(long));
      if (!is_futile(p->back)) {
	par = p->back;
	hookup(dtr,par);
	fixlength(op,data,dtr);
      }
      newtree->numcoals--;
    }
  }

  /* remove futile branches from the branchlists */
  for(t = firstfutile; t != NULL; t = tsucc) {
    tsucc = t->succ;
    if (!is_futile(t->eventnode)) continue;
    p = t->eventnode;
    if (isrecomb(p)) {
      if (p == cutnode || otherparent(p) == cutnode) {
	if (p == cutnode) {
	  remove_branch(t,p);
	  rename_branch(t,findunique(p)->back,otherparent(p));
	} else {
	  remove_branch(t,otherparent(p));
	  rename_branch(t,findunique(p)->back,p);
	}
      } else {
	remove_branch(t,p);
	remove_branch(t,otherparent(p));
      }
    } else {
      if (otnode[p->number]) rename_branch(t,otnode[p->number],p);
      else remove_branch(t,p);
    }
  }

  extbranch(curtree,gettymenode(newtree,cutnode->number),cutnode);

  /* scanning tymelist from the top down, updating coal and ranges arrays */
  for (t = firstfutile; t != NULL; t = t->succ) {
    p = findunique(t->eventnode);
    if (!is_futile(p)) {
      if (isrecomb(p)) {
	fix_rec_ranges(op,data,p);
      } else {
	fix_coal_ranges(op,data,p);
      }
      if (op->fc) fix_coal(op,data,newtree,t,p);
    }
  }


  /* scanning tymelist from the top down again, building the brlist,
     and finally removing the remnant tymelist entries */
  for (t = firstfutile; t != NULL; t = tsucc) {
    p = findunique(t->eventnode);
    tsucc = t->succ;
    if (!is_futile(p)) {
      oldnode = findunique(oldtree->nodep[p->number]);
      if (isrecomb(p)) {
	if (op->fc) {
	  oldranges = oldnode->next->coal;
	  newtopnode = p->next;
	  btyme = oldnode->next->back->tyme;
	  addbrlistH(op,data,br,newtopnode,oldranges,newtopnode->coal,
		     NULL,newtopnode->tyme,btyme,TRUE);

	  oldranges = oldnode->next->next->coal;
	  newtopnode = p->next->next;
	  btyme = oldnode->next->next->back->tyme;
	  addbrlistH(op,data,br,newtopnode,oldranges,newtopnode->coal,
		     NULL,newtopnode->tyme,btyme,TRUE);
	} else {
	  oldranges = oldnode->next->ranges;
	  newtopnode = p->next;
	  btyme = oldnode->next->back->tyme;
	  addbrlistH(op,data,br,newtopnode,oldranges,newtopnode->ranges,
		     NULL,newtopnode->tyme,btyme,TRUE);

	  oldranges = oldnode->next->next->ranges;
	  newtopnode = p->next->next;
	  btyme = oldnode->next->next->back->tyme;
	  addbrlistH(op,data,br,newtopnode,oldranges,newtopnode->ranges,
		     NULL,newtopnode->tyme,btyme,TRUE);
	}

      } else {
	if (op->fc) {
	  oldranges = oldnode->coal;
	  newtopnode = p;
	  btyme = oldnode->back->tyme;
	  addbrlistH(op,data,br,newtopnode,oldranges,newtopnode->coal,
		     NULL,newtopnode->tyme,btyme,TRUE);
	} else {
	  oldranges = oldnode->ranges;
	  newtopnode = p;
	  btyme = oldnode->back->tyme;
	  addbrlistH(op,data,br,newtopnode,oldranges,newtopnode->ranges,
		     NULL,newtopnode->tyme,btyme,TRUE);
	}
      }
    } else {
      if (!isrecomb(p)) {
	if (op->fc) {
	  addbrlistH(op,data,br,otnode[p->number],oranges[p->number],
		     obnode[p->number]->back->coal,NULL,p->tyme,
		     obnode[p->number]->tyme,TRUE);
	} else {
	  addbrlistH(op,data,br,otnode[p->number],oranges[p->number],
		     obnode[p->number]->back->ranges,NULL,p->tyme,
		     obnode[p->number]->tyme,TRUE);
	}
      }
      if (isrecomb(p)) subtymelist(t,p->next,TRUE);
      else subtymelist(t,p,TRUE);
      freenode(p);
    }
  }

  /* The following loop deals with the case (and probably others as well)
     where the cutnode is the top of a recombinant loop and the non-cut
     branch needs to be extended down to replace the bottom nodes' branch
  */
  for(i=0; i<topcount; i++) {
    t = gettymenode(newtree,tops[i]->number);
    readdbranch(newtree,t,tops[i]);
  }

  cutnode->back = NULL;

  free(obnode);
  free(otnode);
  free(tops);
  for(i = 0; i < numnodes; i++) if (oranges[i]) free(oranges[i]);
  free(oranges);

} /* remove_futile */


/**************************************************************
 *  Re-introduces a branch into the tymelist after it has      *
 *  been stripped out during remove_futile.  This code is      *
 *  quite specific to remove_futile; don't use it generically. */
void readdbranch(tree *newtree, tlist *u, node *addbranch)
{
  long i;
  tlist *t;
  boolean found;
  
  for (t = u; t!=NULL; t=t->succ) {
    if(t->eventnode->number==addbranch->back->number) return;
    found = FALSE;
    for (i = 0; i < t->numbranch; i++) {
      if (t->branchlist[i]==addbranch) found = TRUE;
    }
    if (!found) {
      t->numbranch++;
      t->branchlist = (node **)realloc(t->branchlist,
				       t->numbranch*sizeof(node *));
      t->branchlist[t->numbranch - 1] = addbranch;
    }
  }
} /* addbranch */


/***********************************************************************
 * rembranch() removes the passed branch from all tymelist branchlists *
 * from "start" to the bottom/end of the tree.                         */
void rembranch(tree *newtree, tlist *start, node *sbranch)
{
  long i;
  tlist *t;
  boolean found;
  node **newlist;

  for (t = start; t!=NULL; t=t->succ) {
    newlist = (node **)calloc(t->numbranch,sizeof(node *));
    found = FALSE;
    for(i = 0; i < t->numbranch; i++) {
      if (!found) newlist[i] = t->branchlist[i];
      else newlist[i-1] = t->branchlist[i];
      if (t->branchlist[i] == sbranch) found = TRUE;
    }
    free(t->branchlist);
    t->branchlist = newlist;
    if (found) t->numbranch--;
  }

} /* rembranch */



/******************************************************************
 * extbranch() adds the passed branch to all tymelist branchlists *
 * from "start" to the bottom/end of the tree.                    */
void extbranch(tree *newtree, tlist *start, node *addbranch)
{
  long i;
  tlist *t;
  boolean found;

  for (t = start; t!=NULL; t=t->succ) {
    found = FALSE;
    for (i = 0; i < t->numbranch; i++) {
      if (t->branchlist[i]==addbranch) found = TRUE;
    }
    if (!found) {
      t->numbranch++;
      t->branchlist = (node **)realloc(t->branchlist,
				       t->numbranch*sizeof(node *));
      t->branchlist[t->numbranch - 1] = addbranch;
    }
  }

} /* extbranch */


boolean isactive(node *p, linlist *u)
{
  linlist *t;

  for(t = u; t != NULL; t = t->succ)
    if (p == t->branchtop) return(TRUE);

  return(FALSE);

} /* isactive */


long countinactive(option_struct *op, tlist *t, linlist *u)
{
  long i, numinactive, *actives;
  tlist *activetime = NULL;
  
  activetime = t;
  while (1){
    if (activetime->update != 1) break;
    activetime = activetime->prev;
  }
  numinactive = activetime->numbranch;
  for (i = 0; i < activetime->numbranch; i++) {
    if (isactive(activetime->branchlist[i],u)) {numinactive--; continue;}
    if (activetime->branchlist[i] != NULL){
      actives = activetime->branchlist[i]->coal;
      if (!actives[0]) {numinactive--; continue;}
    }
    else numinactive --;
  }

  return(numinactive);

} /* countinactive */


/*********************************************************************
 * pickbrtarget() returns a randum brlist segment from the passed in *
 * brlist using the passed in number of active sites as a guide.     */
brlist *pickbrtarget(brlist *start, double numactive)
{
  double pick;
  brlist *br;

  pick = randum() * numactive - start->numnewactives;
  for(br = start; pick > 0; br = br->succ) pick -= br->succ->numnewactives;

  return(br);

} /* pickbrtarget */


/*********************************************************************
 * pickbrcross() picks a randum point of crossover within the passed *
 * in brsegment.                                                     */
long pickbrcross(option_struct *op, data_fmt *data, brlist *source)
{
  long cross;
  double numactive;
  boolean found;

  numactive = randum() * source->numnewactives;
  found = FALSE;

  if (source->nfsite < source->ofsite && source->nlsite > source->olsite) {
    for(cross = source->nfsite; cross < source->ofsite && !found; cross++) {
      if (op->datatype == 'n') numactive--;
      else numactive -= getdata_space(op,data,cross);
      if (numactive < 0) {found = TRUE; break;}
    }
    if (!found)
      for(cross = source->olsite; cross < source->nlsite && !found; cross++) {
	if (op->datatype == 'n') numactive--;
	else numactive -= getdata_space(op,data,cross);
	if (numactive < 0) {found = TRUE; break;}
      }
  } else {
    if (source->nlsite > source->olsite && source->nfsite < source->olsite)
      for(cross = source->olsite; cross < source->nlsite && !found; cross++) {
	if (op->datatype == 'n') numactive--;
	else numactive -= getdata_space(op,data,cross);
	if (numactive < 0) {found = TRUE; break;}
      }
    else
      for(cross = source->nfsite; cross < source->nlsite && !found; cross++) {
	if (op->datatype == 'n') numactive--;
	else numactive -= getdata_space(op,data,cross);
	if (numactive < 0) {found = TRUE; break;}
      }
  }

  return(cross);


} /* pickbrcross */


/*********************************************************************
 * eventprobbr() is an auxiliary function to eventprob that helps it *
 * deal with brlist segments.                                        *
 *                                                                   *
 * Due to brsegment endpoint problems, (what endpoints to include?), *
 * eventprobbr() now also picks both the brlist segment that will    *
 * have an event, if one is desired; and what the crossover point    *
 * will be.  Segment and point info are passed in the bseg struct.   */
double eventprobbr(option_struct *op, data_fmt *data, brlist **bigbr,
		   tlist *tint, double pc, double pr, double *pb, bseg *brs)
{
  double ttyme, btyme, tyme, offset, brlength, stoptyme, bsitesum;
  brlist *br, *brsucc, *finish;

  bsitesum = 0.0;
  offset = 0.0;
  stoptyme = tint->age;
  tyme = BOGUSTREETYME;
  brlength = BOGUSTREETYME - 1.0;
  btyme = tint->eventnode->tyme;
  finish = NULL;

  while (tyme > brlength) {
    if (!(*bigbr)) {
      ttyme = btyme;
      btyme = stoptyme;
      (*pb) = 0.0;
    } else {
      if ((*bigbr)->starttyme > btyme) {
	ttyme = btyme;
	btyme = (((*bigbr)->starttyme < stoptyme) ? 
		 (*bigbr)->starttyme : stoptyme);
	(*pb) = 0.0;
	finish = (*bigbr);
      } else {
	ttyme = (*bigbr)->starttyme;
	btyme = DBL_MAX;
	bsitesum = 0;
	for(br = (*bigbr); br->starttyme == ttyme; br = br->succ) {
	  if (br->endtyme < btyme) btyme = br->endtyme;
	  bsitesum += br->numnewactives;
	  if (br->succ == NULL) {br = NULL; break;}
	}
	if (btyme > stoptyme) btyme = stoptyme;
	/* now also find anything that starts before the shortest entry ends;
	   finish will end up pointing to the first entry we don't want to do
	   anything with.  Assumption of contingous brlist entries contingent
	   upon the list being sorted by starttymes, least to greatest! */
	finish = br;
	if (br) if (br->starttyme < btyme) btyme = br->starttyme;
	*pb = bsitesum*rec0;
      }
    }

    if (*pb) {
      brs->target = pickbrtargetH((*bigbr),bsitesum);
      brs->cross = pickbrcrossH(op,data,brs->target);
    }

    /* when? */
    if (pc || pr || *pb) tyme = -log(randum())/(pc + pr + (*pb));
    else tyme = btyme + 0.1;

    brlength = btyme - ttyme;
    if (tyme > brlength) {
      if (*bigbr)
	for(br = (*bigbr); br != finish; br = brsucc) {
	  brsucc = br->succ;
	  br->starttyme = btyme;
	  if (br->starttyme == br->endtyme) br_remove(bigbr,br);
	}
      offset += btyme - ttyme;
    }

    if (btyme == stoptyme) break;

  } 

  return(tyme + offset);

} /* eventprobbr */


/**********************************************************************
 * eventprob() calculates both what kind of event is the "next" event *
 * and what the time is to that event.                                */
double eventprob(option_struct *op, data_fmt *data, linlist *alines,
		 tlist *t, brlist **bigbr, char *event, bseg *brs)
{
  long numilines, numalines;
  double pc, pr, pb, tester, tyme, asitesum;
  linlist *u;

  if (!alines && (!*bigbr)) {(*event) = '0'; return(0.0);}

  /* number of inactive lineages */
  numilines = countinactive(op,t,alines);

  /* count the active lineages and sum active sites */
  numalines = 0;
  asitesum = 0.0;
  for (u = alines; u != NULL; u = u->succ) {
    numalines++;
    asitesum += u->activesites;
  }

  pc = (numalines*(numalines-1) + 2.0*numalines*numilines) / theta0;
  pr = asitesum*rec0;

  if (!(*bigbr)) {pb = 0.0; tyme = -log(randum())/(pc + pr);}
  else tyme = eventprobbr(op,data,bigbr,t,pc,pr,&pb,brs);

  /* what kind of event? */
  if (!pc && !pr && !pb) {(*event) = '0'; return(tyme);}
  tester = randum();
  if (tester <= pr/(pc+pr+pb)) (*event) = 'r';
  else if (tester <= (pc+pr)/(pc+pr+pb)) (*event) = 'c';
  else (*event) = 'b';

  return(tyme);

} /* eventprob */


/******************************************************************
 * traverse_flagbelow sets the "updated" field of the passed node *
 * and all nodes rootwards of it to FALSE.                        */
void traverse_flagbelow(tree *tr, node *p)
{
  tlist *t;

  for(t = gettymenode(tr,p->number); t != NULL; t = t->succ) {
    t->eventnode->updated = FALSE;
    if (!istip(t->eventnode)) {
      t->eventnode->next->updated = FALSE;
      t->eventnode->next->next->updated = FALSE;
    }
  }

} /* traverse_flagbelow */


/*********************************************************************
 * flag_below marks the node "p" iff p is not a tip, and in any case *
 *    marks all nodes "rootwards" of p.                              */
void flag_below(node *p)
{

  if (istip(p) && p != curtree->root) flag_below(p->back);

  if(!istip(p)) {
    while (p->top) p = p->next;
    p->updated = FALSE;
    p->next->updated = FALSE;
    if(p->next->top) flag_below(p->next->back);
    p->next->next->updated = FALSE;
    if(p->next->next->top) flag_below(p->next->next->back);
  }
} /* flag_below */


/*********************************************************************
 * traverse_unflag() sets the "updated" field of the passed node and *
 * all nodes rootwards of it to TRUE.                                */
void traverse_unflag(tree *tr, node *p)
{
  tlist *t;

  for(t = gettymenode(tr,p->number); t != NULL; t = t->succ) {
    t->eventnode->updated = TRUE;
    if (!istip(t->eventnode)) {
      t->eventnode->next->updated = TRUE;
      t->eventnode->next->next->updated = TRUE;
    }
  }

} /* traverse_unflag */


void unflag(node *p)
     /* remove flags.  call with curtree->root->back. */
{
  if(!istip(p)) {
    if(!p->next->top) unflag(p->next->back);
    if(!p->next->next->top) unflag(p->next->next->back);
    p->updated = TRUE;
    p->next->updated = TRUE;
    p->next->next->updated = TRUE;
  }
} /* unflag */

node *pickactive(long numalin, linlist **lineages, node *prohibited)
{
  long daughter;
  linlist *u;

  do {
    daughter = (long)(randum() * numalin) + 1;
    u = *lineages;
    findlin(&u,daughter);
  } while (u->branchtop == prohibited);

  return(u->branchtop);

} /* pickactive */

node *pickinactive(option_struct *op, tlist *t, node *prohibited,
		   linlist *alines)
{
  node *daughter;

  do {
    daughter = t->branchlist[(long)(randum() * t->numbranch)];
  } while (daughter == prohibited || isactive(daughter,alines) || (isrecomb(daughter) && isdead(op,daughter)) || daughter->coal[0] == 0);

  return(daughter);

} /* pickinactive */

void coalesce(option_struct *op, data_fmt *data, tlist *t,
	      linlist **lineages, long *lines, double tyme, boolean **nodenumber,
	      long *numnodes, boolean rootdropping, brlist **br)
{
  long numalin, numilin, *active1, *active2;
  double Pbactive, activelinks;
  node *p, *daughter1, *daughter2;
  boolean bothactive, rootpick;
  dnadata *dna;

  dna = data->dnaptr;

  /* find the number of active (numalin) and inactive (numilin) lineages,
     when rootdropping, then the last inactive lineage (the root)
     becomes active; leaving no inactive lineages. */
  numalin = *lines;
  if (!rootdropping) numilin = countinactive(op,t,*lineages);
  else numilin = 0;
  

  /* pick the 2 lineages to coalesce, daughter1 and daughter2
     Pbactive = the chance that they are both active */
  Pbactive = (numalin)*(numalin-1) /
    ((numalin)*(numalin-1) + 2.0*numalin*numilin);

  if (Pbactive >= randum()) {
    daughter1 = pickactive(numalin,lineages,NULL);
    daughter2 = pickactive(numalin,lineages,daughter1);
    bothactive = TRUE;
  } else {
    daughter1 = pickactive(numalin,lineages,NULL);
    daughter2 = pickinactive(op,t,daughter1,*lineages);
    bothactive = FALSE;
  }

  rootpick = FALSE;
  if (daughter1 == curtree->root->back ||
      daughter2 == curtree->root->back)
    rootpick = TRUE;

  /* coalesce the two partners */
  newnode(&p);
  p->number = setnodenumber(nodenumber,numnodes);
  p->next->number = p->number;
  p->next->next->number = p->number;
  p->top = FALSE;
  free_x(op,p);
  p->next->top = FALSE;
  free_x(op,p->next);
  p->next->next->top = TRUE;
  //allocate_x(op,data,p->next->next);
  if (op->map) allocate_z(op,data,p->next->next);
  ranges_Malloc(p->next->next,TRUE,0L);
  p->tyme = tyme;
  p->next->tyme = tyme;
  p->next->next->tyme = tyme;
  p->updated = FALSE;
  p->next->updated = FALSE;
  p->next->next->updated = FALSE;
  p->type = p->next->type = p->next->next->type = 'c';

  curtree->numcoals += 1;
  curtree->nodep[p->number] = p->next->next;

  if(!bothactive) { 
    hookup(p->next->next,daughter2->back);
    fixlength(op,data,p->next->next); /* but only if not bothactive! */
  }
  hookup(p,daughter1);
  fixlength(op,data,p);
  hookup(p->next,daughter2);
  fixlength(op,data,p->next);
  if (rootpick) {
    hookup(p->next->next,curtree->root);
    fixlength(op,data,p->next->next);
  }

  /* now set the ranges for the new node */
  contrib(op,data,p->next->next,&p->next->next->ranges);

  /* update the tymelist */
  /*readdactive(t->succ,daughter1);*/
  /*if (bothactive) readdactive(t->succ,daughter2);*/
  rembranch(curtree,t->succ,daughter1);
  if (bothactive) {
    rembranch(curtree,t->succ,daughter2);
  }
  insertaftertymelist(t,p);
  if (bothactive) {
    extbranch(curtree,t->succ,p->next->next);
  }

  /* now set the coal array for the new node */
  if (op->fc) fix_coal(op,data,curtree,t->succ,p);

  sublinlist(lineages,daughter1);
  if (bothactive) {
    if (op->fc) activelinks = count_activefc(op,data,p->next->next);
    else activelinks = count_active(op,data,p->next->next);
    addlinlist(lineages,p->next->next, activelinks);
    sublinlist(lineages,daughter2);
  }
  (*lines)--;  

  if(!bothactive) {
    /* rename the changed branchtop segments */
    drop_renamebrlist(br,daughter2,p->next->next);
    /* then fix the rest of the ranges and all brlist segments */
    if (op->fc) {
      active1 = p->next->next->coal;
      active2 = daughter2->coal;
    } else {
      active1 = p->ranges;
      active2 = daughter2->ranges;
      //active1 = p->next->next->ranges;
      //active2 = daughter2->ranges;
    }
    if(!sameranges(active1,active2))
      addbrlistH(op,data,br,p->next->next,active2,
		 active1,NULL,p->next->next->tyme,
		 p->next->next->back->tyme,FALSE);
    drop_brfix(op,data,curtree,t->succ,br);
    consolidate_brlist(op,data,curtree,br);
  }

} /* coalesce */


/***********************************************************************
 * init_numnewactive() sets the numnewactives field for an entire      *
 * brlist, removing any brlist segments whose tyme has passed.         */
void init_numnewactive(option_struct *op, data_fmt *data, tree *tr,
		       brlist **bigbr, node *branchtop)
{
  brlist *br, *brsucc;

  for(br = (*bigbr); br != NULL; br = brsucc) {
    brsucc = br->succ;
    if (branchtop->tyme >= br->endtyme) {br_remove(bigbr,br); continue;}
    if (branchtop->tyme > br->starttyme) {
      br->starttyme = branchtop->tyme;
      if (br->starttyme >= br->endtyme) {br_remove(bigbr,br); continue;}
    }
    br->numnewactives = count_activebr(op,data,br);
  }

} /* init_numnewactive */


/**********************************************************************
 * pickbr_recomb() picks the proper segment to use for re-droping the *
 * new lineage, given the tyme at which the drop must occur.  It then *
 * removes/truncates appropiately all non-chosen segments in the list.*
 * It also picks which site will be the point of crossover.           */
brlist *pickbr_recomb(dnadata *dna, brlist **bigbr, double droptyme, long *cross)
{
  brlist *br, *brsucc, *start, *picked;
  long numbr;
  double pick, numactive;
  boolean found;

  /* initializations to make lint happy */
  start = NULL;
  numbr = 0;
  numactive = 0.0;

  for(br = (*bigbr); br != NULL; br = brsucc) {
    brsucc = br->succ;
    if(!isintyme(br,droptyme)) {br_remove(bigbr,br); continue;}
    numbr = 0;
    numactive = 0.0;
    for(start = br; isintyme(br,droptyme); br = br->succ) {
      numbr++;
      numactive += br->numnewactives;
      if (br->succ == NULL) break;
    }
    break;
  }
  if (!numbr) 
    fprintf(ERRFILE,"ERROR:pickbr_recomb failed to do any picking!\n");

  if (!numactive)
    fprintf(ERRFILE,"ERROR:pickbr_recomb has no active links!\n");


  pick = randum() * numactive - start->numnewactives;
  for(br = start; pick > 0; br = br->succ) pick -= br->succ->numnewactives;
  picked = br;

  numactive = pick+br->numnewactives;
  found = FALSE;

  if (br->nfsite < br->ofsite && br->nlsite > br->olsite) {
    for(*cross = br->nfsite; *cross < br->ofsite && !found; (*cross)++) {
      numactive -= dna->sspace[population][locus][*cross];
      if (numactive < 0) {found = TRUE; break;}
    }
    if (!found)
      for(*cross = br->olsite; *cross < br->nlsite && !found;
          (*cross)++) {
	numactive -= dna->sspace[population][locus][*cross];
	if (numactive < 0) {found = TRUE; break;}
      }
  } else {
    if (br->nlsite > br->olsite && br->nfsite < br->olsite) 
      for(*cross = br->olsite; *cross < br->nlsite && !found;
	  (*cross)++) {
	numactive -= dna->sspace[population][locus][*cross];
	if (numactive < 0) {found = TRUE; break;}
      }
    else
      for(*cross = br->nfsite; *cross < br->nlsite && !found;
	  (*cross)++) {
	numactive -= dna->sspace[population][locus][*cross];
	if (numactive < 0) {found = TRUE; break;}
      }
  }

  for(pick = 0, br = start; pick < numbr; pick++, br = brsucc) {
    brsucc = br->succ;
    if (br != picked) {
      if(droptyme < br->endtyme) br->starttyme = droptyme;
      else br_remove(bigbr,br);
    }
  }

  return(picked);

} /* pickbr_recomb */


/*********************************************************************
 * boolean brlistrecomb() inserts a new recombination into a partial *
 * tree, using a brlist entry as the basis for the recombination.    */
boolean brlistrecomb(option_struct *op, data_fmt *data, tlist *t,
		     linlist **lineages, long *lines, brlist **bigbr, double tyme,
		     long *siteptr, boolean **nodenumber, long *numnodes, bseg *brs, tree *tr, seglist *cursegment)
{
  long target, numsites;
  brlist *br;
  node *d, *p;
  dnadata *dna;
  brlist *newbr;

  dna = data->dnaptr;
  numsites = countsites(op,data);

  br = brs->target;
  target = brs->cross;
  d = br->branchtop;

  /* create recombination node */
  newnode(&p);
  p->number = setnodenumber(nodenumber,numnodes);
  p->next->number = p->number;
  p->next->next->number = p->number;
  p->top = FALSE;
  free_x(op,p);
  p->next->top = TRUE;
  //allocate_x(op,data,p->next);
  if(op->map)allocate_z(op,data,p->next);
  ranges_Malloc(p->next,TRUE,0L);
  p->next->next->top = TRUE;
  //allocate_x(op,data,p->next->next);
  if(op->map)allocate_z(op,data,p->next->next);
  ranges_Malloc(p->next->next,TRUE,0L);
  p->tyme = tyme;
  p->next->tyme = tyme;
  p->next->next->tyme = tyme;
  p->updated = FALSE;
  p->next->updated = FALSE;
  p->next->next->updated = FALSE;
  p->type = p->next->type = p->next->next->type = 'r';

  curtree->numrecombs += 1;
  curtree->nodep[p->number] = p->next;

  newgrec[target]++;

  hookup(p->next,d->back);
  fixlength(op,data,p->next);
  hookup(p,d);
  fixlength(op,data,p);
  fix_rec_ranges(op,data,p);


  /* partition the sites */

  if (p->next->back->oldbackranges[1] > target){
    p->next->recstart = target + 1;
    p->next->recend = numsites - 1;
    p->next->next->recstart = 0;
    p->next->next->recend = target;
  }
  else{
    if(p->next->back->oldbackranges[p->next->back->oldbackranges[0] * 2] < (target + 1) ){
      p->next->recstart = 0;
      p->next->recend = target;
      p->next->next->recstart = target + 1;
      p->next->next->recend = numsites - 1;
    }
    else printf("error at brlistrecomb\n");
  }

  /* now update the ranges for the new node */
  contrib(op,data,p->next,&p->next->ranges);
  contrib(op,data,p->next->next,&p->next->next->ranges);

  /* update the tymelist */
  insertaftertymelist(t,p);
  //extbranch(curtree,t->succ,p->next->next);


  // new fc code

  init_ranges_alloc(&p->ranges,p->back->ranges[0]);
  memcpy((p->ranges),p->back->ranges,(p->back->ranges[0]*2+2)*sizeof(long));

  //contrib_temp(p,&p->ranges);

  contrib_temp(p,&p->coal);
  contrib_temp(p->next,&p->next->coal);
  contrib_temp(p->next->next, &p->next->next->coal);
  

  modifynode(p,t,cursegment);

  newbr = initbrlist(p->next, p->next->tyme, p->next->back->tyme,countsites(op,data));
  updatesegranges(op,data,newbr,NULL,p->next->back->oldbackranges,p->next->coal);
  //updatesegranges(op,data,newbr,NULL,p->back->oldcoal,p->next->coal);
  
  if (newbr->weight > 0){
    hookup_brlist(bigbr,newbr);
    p->next->branch = newbr;
  }
  else freebrlist(newbr);
  //hookup_brlist(bigbr,newbr);
  //p->next->branch = newbr;

  if (p->back->branch != NULL){
    br_remove(bigbr,p->back->branch);
    p->back->branch = NULL;
  }


  /* fix alias entries */
  if (op->datatype == 'n' || op->datatype == 's')
    edit_alias(op,data,siteptr,target);

  /* add new entry to linlist */
  addlinlistH(op,data,lineages,p->next->next,0);
  (*lines)++;

  if (curtree->numrecombs > RECOMB_MAX) return(FALSE);
  return(TRUE);

} /* brlistrecomb */


/********************************************************************
 * choosepsite() chooses the sequence site that a new recombination *
 * will start at.                                                   */
long choosepsite(option_struct *op, data_fmt *data, node *p,
		 double targetlink)
{
  long target;


  for(target = p->ranges[1]; target < p->ranges[2*p->ranges[0]]; target++) {
    if (op->datatype == 'n') targetlink--;
    else targetlink -= getdata_space(op,data,target);
    if (targetlink <= 0) break;
  }

  if (target == p->ranges[2*p->ranges[0]]) {
    fprintf(ERRFILE,"ERROR:choosepsite: failure to choose %ld %ld\n",indecks,
	    apps);
  }

  return(target);

} /* choosepsite */


/**********************************************************************
 * choosepsitefc() chooses the sequence site that a new recombination *
 * will start at.                                                     */
long choosepsitefc(option_struct *op, data_fmt *data, node *p,
		   double targetlink)
{
  long target;

  for(target = p->coal[1]; target < p->coal[2*p->coal[0]]; target++) {
    if (op->datatype == 'n') targetlink--;
    else targetlink -= getdata_space(op,data,target);
    if (targetlink <= 0) break;
  }

  if (target == p->coal[2*p->coal[0]]) {
    fprintf(ERRFILE,"ERROR:choosepsitefc: failure to choose %ld %ld\n",indecks,
	    apps);
  }

  return(target);

} /* choosepsite */


boolean recomb(option_struct *op, data_fmt *data, tlist *t,
	       linlist **lineages, long *lines, double tyme, long *siteptr,
	       boolean **nodenumber, long *numnodes)
{
  long target, numsites;
  double ntarget, numactive;
  linlist *u;
  node *d, *p;
  boolean rootpick;
  dnadata *dna;

  dna = data->dnaptr;
  numsites = countsites(op,data);

  u = *lineages;
  /* find splitting lineage, which must be active */
  /* we first count up all active sites on active lineages -- */
  numactive = 0;
  do {
    numactive += u->activesites;
    u = u->succ;
  } while (u != NULL);
  /* -- then we choose a lineage proportional to the active sites */
  ntarget = randum()*numactive;
  numactive = 0;
  u = *lineages;
  do {
    numactive += u->activesites;
    if (u->activesites != 0 && numactive >= ntarget) break;
    u = u->succ;
  } while(1);
  d = u->branchtop;

  rootpick = FALSE;
  if (d == curtree->root->back) rootpick = TRUE;

  /* create recombination node */
  newnode(&p);
  p->number = setnodenumber(nodenumber,numnodes);
  p->next->number = p->number;
  p->next->next->number = p->number;
  p->top = FALSE;
  free_x(op,p);
  p->next->top = TRUE;
  //allocate_x(op,data,p->next);
  if(op->map)allocate_z(op,data,p->next);
  ranges_Malloc(p->next,TRUE,0L);
  p->next->next->top = TRUE;
  //allocate_x(op,data,p->next->next);
  if(op->map)allocate_z(op,data,p->next->next);
  ranges_Malloc(p->next->next,TRUE,0L);
  p->tyme = tyme;
  p->next->tyme = tyme;
  p->next->next->tyme = tyme;
  p->updated = FALSE;
  p->next->updated = FALSE;
  p->next->next->updated = FALSE;
  p->type = p->next->type = p->next->next->type = 'r';

  curtree->numrecombs += 1;
  curtree->nodep[p->number] = p->next;

  hookup(p,d);
  fixlength(op,data,p);
  if (rootpick) {
    hookup(p->next,curtree->root);
    fixlength(op,data,p->next);
  }

  /* partition the sites */
  numactive = randum()*u->activesites;
  if (op->fc) target = choosepsitefc(op,data,d,numactive);
  else target = choosepsite(op,data,d,numactive);
  if(randum()>0.5) {
    p->next->recstart = 0;
    p->next->recend = target;
    p->next->next->recstart = target+1;
    p->next->next->recend = numsites-1;
  } else {
    p->next->recstart = target+1;
    p->next->recend = numsites-1;
    p->next->next->recstart = 0;
    p->next->next->recend = target;
  }

  /* set the ranges fields for the 2 rootwards branches */
  contrib(op,data,p->next,&p->next->ranges);
  contrib(op,data,p->next->next,&p->next->next->ranges);

  /* update the tymelist */
  /*readdactive(t->succ,d);*/
  rembranch(curtree,t->succ,d);
  insertaftertymelist(t,p);
  extbranch(curtree,t->succ,p->next->next);
  extbranch(curtree,t->succ,p->next);

  /* set the coal fields for the 2 rootwards branches */
  if (op->fc) fix_coal(op,data,curtree,t->succ,p);

  /* fix alias entries */
  if (op->datatype == 'n' || op->datatype == 's')
    edit_alias(op,data,siteptr,target+1);

  /* remove old entry from linlist */
  sublinlist(lineages,d);

  /* add two new entries to linlist */
  if (op->fc) {
    addlinlist(lineages,p->next,count_activefc(op,data,p->next));
    addlinlist(lineages,p->next->next,count_activefc(op,data,p->next->next));
  } else {
    addlinlist(lineages,p->next,count_active(op,data,p->next));
    addlinlist(lineages,p->next->next,count_active(op,data,p->next->next));
  }
  (*lines)++;

  if (curtree->numrecombs > RECOMB_MAX) return(FALSE);
  return(TRUE);
} /* recomb */


boolean find_rootloop(tlist *start, tlist **found)
{
  tlist *t;

  for (t = start; t != NULL; t = t->succ)
    if (t->numbranch == 1)
      if (t->eventnode != curtree->root->back) {
	(*found) = t;
	return(TRUE); 
      }

  (*found) = NULL;
  return(FALSE);

} /* find_rootloop */


/******************************************************************
 * free_from_tymelist truncates the tymelist starting at argument *
 * "start", then cleans up the revenant time slices.              */ 
void free_from_tymelist(tlist *start)
{
  tlist *t, *temp;

  start->prev->succ = NULL;

  for (t = start; t != NULL; t = temp) {
    temp = t->succ;
    if (isrecomb(t->eventnode)) curtree->numrecombs--;
    else curtree->numcoals--;
    fix_treenodep(curtree,t->eventnode->number);
    freenode(t->eventnode);
    freetymenode(t);
  }

} /* free_from_tymelist */


/***********************************************************************
 * update_rangecoal() updates all ranges and coal arrays in the passed *
 * in tree.  It assumes that both the tree and tymelist are correct.   */
void update_rangecoal(option_struct *op, data_fmt *data, tree *tr)
{
  tlist *t;
  node *p;

  for(t = tr->tymelist->succ; t != NULL; t = t->succ) {
    p = findunique(t->eventnode);
    if (isrecomb(p)) {
      fix_rec_ranges(op,data,p);
    } else {
      fix_coal_ranges(op,data,p);
    }
    if (op->fc) fix_coal(op,data,tr,t,p);
  }

} /* update_rangecoal */


void remove_rootloop(option_struct *op, data_fmt *data, long *sp)
{
  tlist *t;
  dnadata *dna;

  dna = data->dnaptr;

  if (curtree->numrecombs == 0) return;

  if (find_rootloop(curtree->tymelist,&t)) {
    hookup(t->eventnode,curtree->root);
    t->age = ROOTLENGTH + t->eventnode->tyme;
    curtree->root->tyme = t->age;
    curtree->root->back->length = ROOTLENGTH;
    ltov(op,data,curtree->root->back);
    free_from_tymelist(t->succ);
    if (op->datatype == 'n' || op->datatype == 's')
      rebuild_alias(op,data,sp);
  }
     
} /* remove_rootloop */


void remove_deadrecombs(option_struct *op, data_fmt *data, tree *tr)
{
  long deadcount;
  node *p;
  tlist *t;

  if (tr->numrecombs < 2) return;

  deadcount = 0;
  for(t = tr->tymelist; t->succ != NULL; t = t->succ) ;
  for(; t != NULL; t = t->prev) {
    p = findunique(t->eventnode);
    if (isrecomb(p)) {
      if (isdead(op,p->next)) {
	remove_pairednodes(op,data,tr,p,p->next);
	deadcount++;
	for(t = tr->tymelist; t->succ != NULL; t = t->succ) ;
      } else {
	if (isdead(op,p->next->next)) {
	  remove_pairednodes(op,data,tr,p,p->next->next);
	  deadcount++;
	  for(t = tr->tymelist; t->succ != NULL; t = t->succ) ;
	}
      }
    }
  }

  if (deadcount) update_rangecoal(op,data,tr);

} /* remove_deadrecombs */


/*********************************************************************
 * remove_excess_tree() removes "dead" branches from the tree:       *
 *    all rootloops are removed, all nodes and branches after the    *
 *        total number of branches in an interval has reached 1,     *
 *    all branches that no longer contain any "live" sites, those    *
 *        sites that eventually reach a tip.                         *
 *                                                                   *
 * Note: these two procedures destroy the exact correspondence       *
 * the "nodep" array and the tree, renumber_nodes() should be called *
 * aftercalling this proceedure.                                     */
void remove_excess_tree(option_struct *op, data_fmt *data, tree *tr,
			long *siteptr)
{
  remove_deadrecombs(op,data,tr);
  remove_rootloop(op,data,siteptr);
} /* remove_excess_tree */


long counttymelist(tlist *t)
     /* count number of entries in tymelist */
{
  long cntr;

  cntr = 0;
  do {
    cntr++;
    t = t->succ;
  } while (t!=NULL);
  return(cntr);
} /* counttymelist */


boolean rootdrop(option_struct *op, data_fmt *data, tlist *t, double starttyme, linlist **lineages, 
		 long *lines, long *siteptr, boolean **nodenumber, long *numnodes,
		 brlist **br)
{
  double basetyme, newtyme;
  char eventtype;
  boolean not_aborted;
  node *p;

  /* if necessary, add the root to the active lineages */
  if (!foundlin(*lineages,curtree->root->back)) {
    if (op->fc) 
      addlinlist(lineages,curtree->root->back,
		 count_activefc(op,data,curtree->root->back));
    else 
      addlinlist(lineages,curtree->root->back,
		 count_active(op,data,curtree->root->back));
    (*lines)++;
  }

  basetyme = starttyme;
  p = curtree->root->back;
  not_aborted = TRUE;
  delete_brlist(br);

  while (*lines > 1) {
    newtyme = basetyme + eventprob(op,data,*lineages,t,br,&eventtype,NULL);
    switch ((int)eventtype) {
    case 'c' :
      coalesce(op,data,t,lineages,lines,newtyme,nodenumber,
	       numnodes,TRUE,br);
      break;
    case 'r' :
      not_aborted = recomb(op,data,t,lineages,lines,
			   newtyme,siteptr,nodenumber,numnodes);
      break;
    case 'b' :
      fprintf(ERRFILE,"ERROR:rootdrop case b: can't get here\n");
					break;
    case '0' :
      fprintf(ERRFILE,"ERROR:rootdrop case 0: can't get here\n");
					break;
    default :
      fprintf(ERRFILE,"ERROR:rootdrop case: can't get here\n");
      break;
    }
    if (!not_aborted) break;
    t = t->succ;
    basetyme = newtyme;

  }


  if (not_aborted) traverse_flagbelow(curtree,p);
  return(not_aborted);

} /* rootdrop */


node *pickbranch(option_struct *op, data_fmt *data, tree *source, 
		 tlist **tymeslice, double *offset)
{
  long branch, numbranches;
  node *p;
  tlist *t;

  (*offset) = 0.0;

  if (recarray[seq_length] == 1)numbranches = getdata_numtips(op,data);
  else;
  numbranches = countbranches(op,data,source);

  branch = (long)(randum() * numbranches) + 1;

  if (branch <= getdata_numtips(op,data)) {
    if (tymeslice) (*tymeslice) = source->tymelist;
    return(source->nodep[branch]);
  }
  branch -= getdata_numtips(op,data);

  for (t = source->tymelist->succ; t != NULL; t = t->succ) {
    if (isrecomb(t->eventnode)) branch = branch - 2;
    else if (findunique(t->eventnode)->coal[0] != 0) branch = branch - 1;
    if (branch < 1) {
      (*tymeslice) = t;
      p = findunique(t->eventnode);
      if (isrecomb(p)) 
	return((randum() < 0.5) ? p->next : p->next->next);
      else
	return(findtop(p));
    }
  }

  fprintf(ERRFILE,"ERROR:pickbranch--failure to find branch\n");
  return(NULL);

} /* pickbranch */


void renumber_nodes(option_struct *op, data_fmt *data, tree *trii)
{
  tlist *t;
  long nodenumber;
	int i;

  /* since the root is "0", the first non-tip node will be "numseq+1" */
  nodenumber = getdata_numtips(op,data) + 1;
	for (i = 0; i < nodenumber; i++){
		trii->nodep[i]->id = i;
		trii->nodep[i]->number = i;
	}

  for (t = trii->tymelist->succ; t != NULL; t = t->succ) {
    t->eventnode->number = nodenumber;
    t->eventnode->next->number = nodenumber;
	  t->eventnode->next->next->number = nodenumber;
	  t->eventnode->id = nodenumber;
	  t->eventnode->next->id = nodenumber;
	  t->eventnode->next->next->id = nodenumber;
    trii->nodep[nodenumber] = t->eventnode;
    nodenumber++;
  }

} /* renumber_nodes */


void renumber_nodes_S(option_struct *op, data_fmt *data, tree *trii)
{
	tlist *t;
	long nodenumber;
	int i;
	node **templist;
	
	//printtreesum_list();

	/* since the root is "0", the first non-tip node will be "numseq+1" */
	nodenumber = old_hap + new_hap + 1;
	
	for (t = trii->tymelist->succ; t != NULL; t = t->succ) {
		nodenumber++;
	}
	//printf("node check %d %d\n",nodenumber, curtree->numcoals);
	templist = (node **)calloc(nodenumber,sizeof(node *));
	
	nodenumber = old_hap + new_hap + 1;
	for (i = 0; i < nodenumber; i++){ // copy tips and root
	  templist[i] = trii->nodep[i];
  
		templist[i]->id = i;
		templist[i]->number = i;
	}
       
	
	for (t = trii->tymelist->succ; t != NULL; t = t->succ) { // update internal nodes
	
	  templist[nodenumber] = t->eventnode;
		templist[nodenumber]->number = nodenumber;
		templist[nodenumber]->next->number = nodenumber;
		templist[nodenumber]->next->next->number = nodenumber;
		templist[nodenumber]->id = nodenumber;
		templist[nodenumber]->next->id = nodenumber;
		templist[nodenumber]->next->next->id = nodenumber;
		  nodenumber++;
	}
	free(curtree->nodep);
	curtree->nodep = templist;
	//makeprinttreesum_list();
} /* renumber_nodes */

void renumber_nodes_S_newtip(option_struct *op, data_fmt *data)
{
	tlist *t;
	long nodenumber;
	int i;
	node **templist;
	
	//printtreesum_list();

	/* since the root is "0", the first non-tip node will be "numseq+1" */
	nodenumber = getdata_numtips(op,data) + 1;
	if (sim_mode)
		nodenumber = nodenumber + new_hap;
	
	for (t = curtree->tymelist->succ; t != NULL; t = t->succ) {
		nodenumber++;
	}
	nodenumber = nodenumber + new_hap; // 
	//printf("node check %d %d\n",nodenumber, curtree->numcoals);
	templist = (node **)calloc(nodenumber,sizeof(node *));
	
	nodenumber = getdata_numtips(op,data) + 1;
	

	for (i = 0; i < nodenumber; i++){ // copy tips and root
	  templist[i] = curtree->nodep[i];
  
		templist[i]->id = i;
		templist[i]->number = i;
	}
	
	nodenumber = nodenumber + new_hap;
	for (t = curtree->tymelist->succ; t != NULL; t = t->succ) { // update internal nodes
	  
	  templist[nodenumber] = t->eventnode;
		templist[nodenumber]->number = nodenumber;
		templist[nodenumber]->next->number = nodenumber;
		templist[nodenumber]->next->next->number = nodenumber;
		templist[nodenumber]->id = nodenumber;
		templist[nodenumber]->next->id = nodenumber;
		templist[nodenumber]->next->next->id = nodenumber;
		nodenumber++;
		
	}
	free(curtree->nodep);
	curtree->nodep = templist;
	//makeprinttreesum_list();
} /* renumber_nodes */


void finishbadtree(option_struct *op, data_fmt *data, tree *tr, linlist **lineages, long *numlins, 
		   boolean **nodenumber, long *numnodes, boolean rootdropped, 
		   brlist **brlines)
{
  double ntyme;
  tlist *t;

  /* go to the end of the tymelist */
  for(t = tr->tymelist; t->succ != NULL; t = t->succ) ;

  while (*numlins > 1) {
    ntyme = t->eventnode->tyme + 0.01;
    coalesce(op,data,t,lineages,numlins,ntyme,nodenumber,numnodes,rootdropped,
	     brlines);
    t = t->succ;
  }
  delete_brlist(brlines);

} /* finishbadtree */


/***********************************************************
 * initbrlist() allocates and initializes a brlist element */
brlist *initbrlist(node *branchtop, double starttyme, double endtyme,
		   long numsites)
{
  brlist *newbr;

  newbr = (brlist *)calloc(1,sizeof(brlist));
  newbr->branchtop = branchtop;
  newbr->segranges = (int *)calloc(numsites,sizeof(int));
  newbr->starttyme = starttyme;
  newbr->endtyme = endtyme;
  newbr->updated = FALSE;
  return(newbr);

} /* initbrlist */


/**************************************************
 * freebrlist() just frees an element of a brlist */
void freebrlist(brlist *br)
{

  free(br->segranges);
  free(br);

} /* freebrlist */


/*********************************************************
 * delete_brlist() deletes an entire brlist, setting the *
 * base pointer to NULL.                                 */
void delete_brlist(brlist **bigbr)
{
  brlist *br, *brsucc;

  if (!bigbr) return;

  for(br = (*bigbr); br != NULL; br = brsucc) {
    brsucc = br->succ;
    br_remove(bigbr,br);
  }

  (*bigbr) = NULL;

} /* delete_brlist */


/******************************************************************
 * hookup_brlist() mechanically adds a new element into a brlist. *
 * New elements are inserted so that the list stays time-sorted   *
 * (least to greatest).                                           */
void hookup_brlist(brlist **bigbr, brlist *newbr)
{
  brlist *br;

  if (!*bigbr) { /* new list being created */
    (*bigbr) = newbr;
    newbr->succ = newbr->prev = NULL;
  } else {
    br = getbrlistbytyme((*bigbr),newbr->starttyme);
    if (!br) { /* adding to end of list */
      for (br = (*bigbr); br->succ != NULL; br = br->succ)
	;
      br->succ = newbr;
      newbr->prev = br;
      newbr->succ = NULL;
    } else {
      newbr->succ = br;
      newbr->prev = br->prev;
      if (br->prev) br->prev->succ = newbr;
      br->prev = newbr;
      if (*bigbr == br) (*bigbr) = newbr;
    }
  }

} /* hookup_brlist */


/*******************************************************************
 * set_sitesbr() sets the "sites" information for a brlist segment */
void set_sitesbr(option_struct *op, data_fmt *data, brlist *br)
{
  long i, numsites;
  boolean done1, done2;

  numsites = countsites(op,data);

  done1 = FALSE;
  done2 = FALSE;

  br->nfsite = numsites;
  br->ofsite = numsites;

  for (i = 0; i < numsites; i++) {
    if (!done1) 
      if (br->segranges[i] >= 1) {
	br->nfsite = i;
	done1 = TRUE;
      } 
    if (!done2)
      if (br->segranges[i] == -1 || 
	  br->segranges[i] == 2) {
	br->ofsite = i;
	done2 = TRUE;
      }
    if (done1 && done2) break;
  }
      
  done1 = FALSE;
  done2 = FALSE;

  br->nlsite = br->nfsite - 1;
  br->olsite = br->ofsite - 1;

  for (i = numsites-1; i >= 0; i--) {
    if (!done1)
      if (br->segranges[i] >= 1) {
	br->nlsite = i;
	done1 = TRUE;
      }
    if (!done2)
      if (br->segranges[i] == -1 || 
	  br->segranges[i] == 2) {
	br->olsite = i;
	done2 = TRUE;
      }
    if (done1 && done2) break;
  }

} /* set_sitesbr */


/***************************************************************
 * count_activebr() returns the number of active links present *
 * in a brlist segment.                                        */
double count_activebr(option_struct *op, data_fmt *data, brlist *br)
{

  if (!br->updated) {
    set_sitesbr(op,data,br);
    br->updated = TRUE;
  } 

  if (br->nfsite == countsites(op,data)) return(0L);

  if (br->olsite <= br->nfsite || br->ofsite >= br->nlsite) {
    return(getnumlinks(op,data,br->nfsite,br->nlsite));
  }

  if (br->nfsite < br->ofsite && br->nlsite > br->olsite) {
    return(getnumlinks(op,data,br->nfsite,br->ofsite) +
	   getnumlinks(op,data,br->olsite,br->nlsite));
  }

  if (br->nfsite < br->ofsite) {
    return(getnumlinks(op,data,br->nfsite,br->ofsite));
  }

  if (br->nlsite > br->olsite) {
    return(getnumlinks(op,data,br->olsite,br->nlsite));
  }

  return(0.0);


} /* count_activebr */


/***************************************************************
 * isintyme() checks to see if a given tyme is within the range *
 * covered by a particular brlist entry.                       */
boolean isintyme(brlist *br, double target)
{

  return((target >= br->starttyme) && (target <= br->endtyme) ? 
	 TRUE : FALSE);

} /* isintyme */


/******************************************************************
 * getbrlistbystart() returns a pointer to the brlist entry which *
 * contains the nodelet; returning NULL if nothing is found.      */
brlist *getbrlistbystart(brlist *start, node *target, double starttyme)
{
  brlist *br;
  long match;

  for(br = start; br != NULL; br = br->succ) {
    match = sametreebranch(br,target,starttyme,BOGUSTREETYME);
    if (match == 1 || match == 2) return(br);
  }

  return(NULL);

} /* getbrlistbystart */


/****************************************************************
 * getbrlistbyend() returns a pointer to the brlist entry which *
 * contains the nodelet; returning NULL if nothing is found.    */
brlist *getbrlistbyend(brlist *start, node *target, double endtyme)
{
  brlist *br;
  long match;

  for(br = start; br != NULL; br = br->succ) {
    match = sametreebranch(br,target,BOGUSTREETYME,endtyme);
    if (match == -1 || match == 2) return(br);
  }

  return(NULL);

} /* getbrlistbyend */


/********************************************************************
 * getbrlistbytyme() finds the first element in a brlist that has a *
 * starttyme >= the searchtyme, returning NULL if nothing is found. */
brlist *getbrlistbytyme(brlist *start, double searchtyme)
{
  brlist *br;

  for(br = start; br != NULL; br = br->succ)
    if (br->starttyme >= searchtyme) return(br);

  return(NULL);

} /* getbrlistbytyme */


/*******************************************************
 * updatesegranges() updates the segranges of a brlist *
 * given the passed arrays.                            *
 *                                                     *
 * segranges is 0 if site has remained dead            *
 *           is +1 for becoming live                   *
 *           is -1 for becoming dead                   *
 *           is +2 for remaining alive                 */
void updatesegranges(option_struct *op, data_fmt *data, brlist *br,
		     int *oldsegranges, long *oldranges, long *newranges)
{
  long i, newi, newstart, newend, numsites;

  long oldind, newind, cstart, cend, csize, oval, nval;
  boolean oldlive, newlive;
  dnadata *dna;
  int *bstart;

  br->updated = FALSE;

  numsites = countsites(op,data);

  if (!oldranges && !newranges && oldsegranges) {
    memcpy(br->segranges,oldsegranges,numsites*sizeof(int));
    return;
  }

  dna = data->dnaptr;
  newi = 1;
  newstart = newranges[newi];
  if (newstart == FLAGLONG) {
    newstart = numsites;
    newend = numsites;
  } else newend = newranges[newi+1];
  if (oldranges) {
    if (oldranges[0]) {
      oldlive = FALSE;
      newlive = FALSE;
      oldind = newind = 1;
      oval = oldranges[oldind];
      nval = ((newstart != numsites) ? newranges[newind] : numsites);
      cstart = 0;
      cend = ((nval < oval) ? nval : oval);
      while (1) {
	csize = cend - cstart;
	bstart = &(br->segranges[cstart]);
	if (csize) {
	  if (oldlive && newlive)
	    memcpy(bstart,dna->segranges2,csize*sizeof(int));
	  else if (!oldlive && newlive)
	    memcpy(bstart,dna->segranges1,csize*sizeof(int));
	  else if (!oldlive && !newlive)
	    memcpy(bstart,dna->segranges0,csize*sizeof(int));
	  else if (oldlive && !newlive)
	    memcpy(bstart,dna->segranges3,csize*sizeof(int));
	}
	if (cend == numsites) break;
	if (oval == cend) {
	  oldind++;
	  oldlive = !oldlive;
	  if (oldlive) oval = oldranges[oldind] + 1;
	  else oval = oldranges[oldind];
	  if (oval == FLAGLONG)
	    oval = numsites;
	}
	if (nval == cend) {
	  newind++;
	  newlive = !newlive;
	  if (newlive) nval = newranges[newind] + 1;
	  else nval = newranges[newind];
	  if (nval == FLAGLONG)
	    nval = numsites;
	}
	cstart = cend;
	cend = (nval < oval) ? nval : oval;
      }
    }
  } else {
    for(i = 0; i < numsites; i++) {
      if (i > newend) {
	newi+=2;
	if (newranges[newi] == FLAGLONG) {
	  newstart = numsites;
	  newend = numsites;
	} else {
	  newstart = newranges[newi];
	  newend = newranges[newi+1];
	}
      }
      if (i >= newstart && i <= newend) {
	if (oldsegranges[i] == 0) br->segranges[i] = 1;
	else if (oldsegranges[i] == -1) br->segranges[i] = 2;
	else br->segranges[i] = oldsegranges[i];
      } else {
	if (oldsegranges[i] == 2) br->segranges[i] = -1;
	else if (oldsegranges[i] == 1) br->segranges[i] = 0;
	else br->segranges[i] = oldsegranges[i];
      }
         
    }
  }

} /* updatesegranges */


/************************************************************
 * segranges_to_ranges() converts from segranges to ranges. *
 *                                                          *
 * This code may be wrong--it converts segranges to current *
 * range info, not the original?                            */
void segranges_to_ranges(int *segranges, long **ranges, long numsites)
{
  long i, newstart, newend;
  boolean live;

  /* initialization for lint */
  newend = 0;

  for(i = 0, newstart = FLAGLONG, live = FALSE; i < numsites; i++) {
    if((!live && (segranges[i] == -1 || segranges == 0)) ||
       (live && (segranges[i] == 1 || segranges[i] == 2)))
      continue;
    live = !live;
    if(live) {newstart = i; newend = FLAGLONG;}
    if(!live) newend = i-1;
    if(newstart != FLAGLONG && newend != FLAGLONG) {
      addrange(ranges,newstart,newend);
      newstart = FLAGLONG;
    }
  }

} /* segranges_to_ranges */


/*********************************************************
 * addbrlist() adds an element to an existing brlist and *
 * depends on both the tree and tymelist being correct   *
 * from the tips to the tyme of the nodelet newtop.      */
void addbrlist(option_struct *op, data_fmt *data, brlist **bigbr,
	       node *newtop, long *oldranges, long *newranges, int *segranges,
	       double starttyme, double endtyme, boolean building)
{
  double oldendtyme;
  int *tsegranges, *oldsegs;
  brlist *tnewbr, *bnewbr, *newbr;
  long numsites;

  numsites = countsites(op,data);
  if (starttyme == endtyme || !newtop) return;

  tnewbr = getbrlistbystart(*bigbr,newtop,starttyme);
  bnewbr = getbrlistbyend(*bigbr,newtop,endtyme);

  if (!tnewbr && !bnewbr) {
    newbr = initbrlist(newtop,starttyme,endtyme,countsites(op,data));
    if (oldranges && newranges) {
      updatesegranges(op,data,newbr,NULL,oldranges,newranges);
    } else {
      updatesegranges(op,data,newbr,segranges,NULL,NULL);
    }
    hookup_brlist(bigbr,newbr);
    return;
  }

  /* we have the same start & end */
  if (tnewbr && bnewbr) {
    if (tnewbr->endtyme == endtyme) { /* we have the same segment */
      if (newranges) {
	if (building) {
	  updatesegranges(op,data,tnewbr,NULL,oldranges,newranges);
	} else {
	  updatesegranges(op,data,tnewbr,tnewbr->segranges,NULL,newranges);
	}
      } else {
	updatesegranges(op,data,tnewbr,segranges,NULL,NULL);
      }
    } else {
      if (newranges) {
	for(; ; tnewbr = tnewbr->succ) {
	  if (tnewbr->branchtop != bnewbr->branchtop) continue;
	  if (building) {
	    updatesegranges(op,data,tnewbr,NULL,oldranges,newranges);
	  } else {
	    updatesegranges(op,data,tnewbr,tnewbr->segranges,NULL,
			    newranges);
	  }
	  if (tnewbr == bnewbr) break;
	  if (tnewbr->endtyme != bnewbr->starttyme)
	    addbrlistH(op,data,bigbr,newtop,oldranges,newranges,segranges,
		       tnewbr->endtyme,bnewbr->starttyme,building);
	}
      }
    }
    return;
  }

  /* we only have the same start */
  if (tnewbr) {
    oldendtyme = tnewbr->endtyme;
    if (oldendtyme > endtyme) {
      tnewbr->starttyme = endtyme;
      addbrlistH(op,data,bigbr,newtop,oldranges,newranges,segranges,
		 starttyme,endtyme,building);
    } else {
      updatesegranges(op,data,tnewbr,tnewbr->segranges,NULL,newranges);
      addbrlistH(op,data,bigbr,newtop,oldranges,newranges,segranges,
		 oldendtyme,endtyme,building);
    }
    return;
  } else { /* we only have the same end */
    if (starttyme > bnewbr->starttyme) {
      bnewbr->endtyme = starttyme;
      oldsegs = (int *)calloc(numsites,sizeof(int));
      memcpy(oldsegs,bnewbr->segranges,numsites*sizeof(int));
      updatesegranges(op,data,bnewbr,bnewbr->segranges,NULL,newranges);
      tsegranges = (int *)calloc(numsites,sizeof(int));
      memcpy(tsegranges,bnewbr->segranges, numsites*sizeof(int));
      memcpy(bnewbr->segranges,oldsegs,numsites*sizeof(int));
      addbrlistH(op,data,bigbr,newtop,NULL,NULL,tsegranges,starttyme,
		 endtyme,building);
      free(oldsegs);
      free(tsegranges);
    } else {
      addbrlistH(op,data,bigbr,newtop,oldranges,newranges,segranges,
		 starttyme,bnewbr->starttyme,building);
      updatesegranges(op,data,bnewbr,bnewbr->segranges,NULL,newranges);
    }
    return;
  }

} /* addbrlist */


boolean makeevent(option_struct *op, data_fmt *data, tree *growtree,
		  tlist *tymeslice, brlist **brlines, linlist **alines, long *numalines,
		  double newtyme, boolean **nodenumber, long *numnodes, long *siteptr,
		  char eventtype, boolean rootdropped, bseg *brs)
{
  boolean succeeded;

  succeeded = TRUE;
  switch((int)eventtype) {
  case 'c' :
    coalesce(op,data,tymeslice,alines,numalines,newtyme,nodenumber,
	     numnodes,rootdropped,brlines);
    break;
  case 'r' :
    succeeded = recomb(op,data,tymeslice,alines,numalines,newtyme,
		       siteptr,nodenumber,numnodes);
    break;
  case 'b' :
    succeeded = brlistrecomb(op,data,tymeslice,alines,numalines,
			     brlines,newtyme,siteptr,nodenumber,numnodes,brs,growtree, NULL);
    break;
  case '0' :
    break;
  default :
    fprintf(ERRFILE,"ERROR:makeevent impossible case: can't get here\n");
    break;
  }

  return(succeeded);
} /* makeevent */


boolean growlineages(option_struct *op, data_fmt *data, tree *growtree,
		     tree *oldtree, tlist *tstart, brlist **brlines, linlist **alines,
		     long *numalines, double offsetstart, boolean **nodenumber,
		     long *numnodes, long *siteptr)
{
  tlist *tymeslice;
  double btyme, obtyme, newtyme, offset;
  boolean succeeded, rootdropped;
  char eventtype;
  bseg brs;

  tymeslice = tstart;
  offset = offsetstart;
  rootdropped = FALSE;
  succeeded = TRUE;

  while(tymeslice != NULL) {

    /* just drop the root if you're in the last interval but don't make the
       root active til you're below the oldtree->root->back */
    if (tymeslice->succ == NULL) {
      btyme = curtree->root->back->tyme;
      obtyme = oldtree->root->back->tyme;
      while (obtyme > btyme) {
	init_numnewactive(op,data,growtree,brlines,tymeslice->eventnode);
	newtyme = tymeslice->eventnode->tyme + offset +
	  eventprob(op,data,(*alines),tymeslice,brlines,&eventtype,&brs);
	if (newtyme <= obtyme) {
	  succeeded = makeevent(op,data,growtree,tymeslice,brlines,alines,
				numalines,newtyme,nodenumber,numnodes,siteptr,eventtype,
				rootdropped,&brs);
	  if (!succeeded) break;
	  if (succeeded && eventtype != '0') tymeslice = tymeslice->succ;
	}
	btyme = newtyme;
	offset = 0.0;
	if ((*numalines) == 0 && !(*brlines)) break;
      }
      if ((*numalines) == 0 && !(*brlines)) break;
      if (!succeeded) break;
      succeeded = rootdrop(op,data,tymeslice,obtyme,alines,numalines,
			   siteptr,nodenumber,numnodes,brlines);
      rootdropped = TRUE;
      break;
    }

    init_numnewactive(op,data,growtree,brlines,tymeslice->eventnode);
    newtyme = tymeslice->eventnode->tyme + offset +
      eventprob(op,data,(*alines),tymeslice,brlines,&eventtype,&brs);
    if (newtyme <= tymeslice->age) {
      succeeded = makeevent(op,data,growtree,tymeslice,brlines,alines,
			    numalines,newtyme,nodenumber,numnodes,siteptr,eventtype,
			    rootdropped,&brs);
      if ((*numalines) == 0 && !(*brlines)) break;
    }
    if ((*numalines) == 0 && !(*brlines)) break;
    if (!succeeded) break;
    if (tymeslice->succ) tymeslice = tymeslice->succ;
    offset = 0.0;
  }

  if (!succeeded) {
    numdropped++;
    finishbadtree(op,data,growtree,alines,numalines,nodenumber,numnodes,rootdropped,brlines);
  }

  return(succeeded);

} /* growlineages */

boolean makedrop(option_struct *op, data_fmt *data)
{
  long i, numnodes, numalines, nummarkers, *oldsiteptr, *siteptr;
  double offset;
  linlist *alines;
  brlist *brlines;
  tlist *tymeslice;
  boolean accept_change, *nodenumber;
  node *p;
  tree *oldtree;

  oldsiteptr = NULL;
  siteptr = data->siteptr;

  /* copy the original tree configuration */
  oldtree = copytree(op,data,curtree);
  numnodes = 2 * oldtree->numcoals + 2;
  nodenumber = (boolean *)calloc(1,numnodes*sizeof(boolean));
  nummarkers = getdata_nummarkers(op,data);
  for(i = 0; i < numnodes; i++) nodenumber[i] = TRUE;
  if(op->datatype == 'n' || op->datatype == 's') {
    oldsiteptr = (long *)calloc(nummarkers,sizeof(long));
    memcpy(oldsiteptr,siteptr,nummarkers*sizeof(long));
  }

  /* pick a spot to cut */
  alines = NULL;
  p = pickbranch(op,data,curtree,&tymeslice,&offset);
  if (op->fc) addlinlist(&alines,p,count_activefc(op,data,p));
  else addlinlist(&alines,p,count_active(op,data,p));
  numalines = 1;

  brlines = NULL;
  tag_futile(p->back);
  traverse_flagbelow(curtree,p->back);
  remove_futile(op,data,curtree,oldtree,p,nodenumber,&brlines);
  accept_change = TRUE;

  accept_change = growlineages(op,data,curtree,oldtree,tymeslice,&brlines,
			       &alines,&numalines,offset,&nodenumber,&numnodes,siteptr);

  if (accept_change) {
    renumber_nodes(op,data,curtree);
    remove_excess_tree(op,data,curtree,siteptr);
    renumber_nodes(op,data,curtree);
#if ALWAYS_ACCEPT
    curtree->likelihood = oldtree->likelihood;
    accept_change = testratio(op,data,oldtree,curtree,'d');
#endif
#if !ALWAYS_ACCEPT
    traverse_flagbelow(curtree,p->back);
    if (op->datatype == 'n' || op->datatype == 's')
      rebuild_alias(op,data,siteptr);
    localeval(op,data,curtree->root->back,FALSE);
    if (marker1pos == - 99){
      curtree->coalprob = coalprob(op,data,curtree,theta0,rec0);
      accept_change = testratio(op,data,oldtree,curtree,'d');
    }
    //else //testratio(op, data, oldtree, curtree);
#endif

  }

#if ALWAYS_REJECT
  op->numout[1]++;
  scoretree(op,data,0L);
  accept_change = FALSE;
#endif

  free(nodenumber);
  freelinlist(alines);

  if (accept_change) {
    freetree(op,data,oldtree);
    if (op->datatype == 'n' || op->datatype == 's') free(oldsiteptr);
    return(TRUE);
  } else {
    freetree(op,data,curtree);
    curtree = oldtree;
    if (op->datatype == 'n' || op->datatype == 's') {
      memcpy(siteptr,oldsiteptr,nummarkers*sizeof(long));
      free(oldsiteptr);
    }
    return(FALSE);
  }

} /* makedrop */


node *addtip_S(option_struct *op, data_fmt *data, tree *source, 
			   tlist **tymeslice, double *offset)
{
	node *p;
	int nodenumber, curtip, i;
	
	nodenumber = old_hap + new_tip + 1;
	(*offset) = 0.0;
	newnode(&p);
	p->number = nodenumber;
	p->next->number = nodenumber;
	p->next->next->number = nodenumber;
	curtree->nodep[nodenumber] = p;
	copynode(op,data,curtree->nodep[1], p);
	p->back = NULL;
	p->next->back = NULL;
	p->next->next->back = NULL;
	
	source->tymelist->numbranch = old_hap + new_tip;
	if (tymeslice) (*tymeslice) = source->tymelist;	
	
	
	curtip = old_hap + new_tip;
	
	source->nodep[curtip] = p;
	source->tymelist->numbranch = curtip;
	free(source->tymelist->branchlist);
     source->tymelist->branchlist = (node **)calloc(source->tymelist->numbranch, sizeof(node *));	
	for (i = 1; i <= curtip; i++){
		source->tymelist->branchlist[i - 1] = source->nodep[i];
	}

	curtree->tymelist->segments->numsam = old_hap + new_tip;
	
	return(p);
	
	
}


void addnewtips(option_struct *op, data_fmt *data)
{
	node *p;
	int nodenumber, curtip, i;
	// update nodep array;
	renumber_nodes_S_newtip(op,data);
	
	// make newtip nodes
	
	for (new_tip = 1; new_tip <= new_hap; new_tip++){
	  nodenumber = old_hap + new_tip;
	  newnode(&p);

	  copynode(op,data,curtree->nodep[1], p);
	  p->number = nodenumber;
	  p->next->number = nodenumber;
	  p->next->next->number = nodenumber;
	  p->id = nodenumber;
	  p->next->id = nodenumber;
	  p->next->next->id = nodenumber;
	  p->back = NULL;
	  p->next->back = NULL;
	  p->next->next->back = NULL;

	  curtree->nodep[new_tip + old_hap] = p;
	}
	
	
	curtree->tymelist->numbranch = old_hap + new_hap;
	curtip = old_hap + new_hap;
	
	
	curtree->tymelist->numbranch = curtip;
	free(curtree->tymelist->branchlist);
        
	curtree->tymelist->branchlist = (node **)calloc(curtree->tymelist->numbranch, sizeof(node *));	
	for (i = 1; i <= curtip; i++){
		curtree->tymelist->branchlist[i - 1] = curtree->nodep[i];
	}

	curtree->tymelist->segments->numsam = old_hap + new_hap;
	
	return;
}


void updatefirsttyme_S(){
	int curtip;
	
	curtip = old_hap + new_tip;
	
	curtree->tymelist->numbranch = curtip;
	free(curtree->tymelist->branchlist);
    curtree->tymelist->branchlist = (node **)calloc(curtree->tymelist->numbranch, sizeof(node *));	
	return;
}
  

boolean makedrop_S_new(option_struct *op, data_fmt *data){ // simulate whole new tips once.
	long i, numnodes, numalines, nummarkers, *oldsiteptr, *siteptr; 
	double offset = 0;
	linlist *alines;
	brlist *brlines;
	tlist *tymeslice;
	boolean accept_change, *nodenumber;

	tree *oldtree = NULL;
	
	
	oldsiteptr = NULL;
	siteptr = data->siteptr;
	// test
	curtree->root->back->coal[0] = 0;
	
	constructoldranges(curtree);
		if (oldsiteptr != NULL) free(oldsiteptr);
      	numnodes = curtree->numcoals + curtree->numrecombs + *data->dnaptr->numseq + 1;
	
	curtree->oldnodep = NULL;
	nodenumber = (boolean *)calloc(1,numnodes*sizeof(boolean));
	nummarkers = getdata_nummarkers(op,data);
	for(i = 0; i < numnodes; i++) nodenumber[i] = TRUE;
	if(op->datatype == 'n' || op->datatype == 's') {
		oldsiteptr = (long *)calloc(nummarkers,sizeof(long));
		memcpy(oldsiteptr,siteptr,nummarkers*sizeof(long));
	}
	/*add new tips */
	alines = NULL;
	numalines = 0;
       
	addnewtips(op,data);

	// add new tips to alines
	for (i = old_hap + 1; i <= old_hap + new_hap; i++){
	  addlinlistH(op,data,&alines,curtree->nodep[i],0);
	  numalines++;
	}

	tymeslice = curtree->tymelist;
	brlines = NULL;
	
	//marknode(p);
	
	accept_change = TRUE;
	accept_change = growlineagesH(op,data,curtree,oldtree,tymeslice,&brlines,&alines,&numalines,offset,&nodenumber,&numnodes,siteptr);
	
	//	printf("3\n");
	//checktymelist();
	
	if (accept_change) {
		renumber_nodes_S(op,data,curtree);
		//traverse_flagbelow(curtree,p->back);
		//if (op->datatype == 'n' || op->datatype == 's')
		//	rebuild_alias(op,data,siteptr);
		//summerize tree
		//scorerecs_temp(curtree);
		//calculate P(G|parameters)
		//if (marker1pos == -99 && op->ctemp != 1) curtree->coalprob = treeprob(curtree);
		//calculate P(D|G)
		
		//calcratio(op,data,oldtree,curtree,p->coal);     
	} 


	//checktymelist();
	
	free(nodenumber);
	freelinlist(alines);
	
	//accept_change = testratio(op,data,oldtree,curtree,'d'); 
	
	//for testing
	//accept_change = TRUE;
	
	
	curtree->recsummary = NULL;
	if (op->datatype == 'n' || op->datatype == 's')  free(oldsiteptr);
	return TRUE;
	
}



boolean makedrop_S(option_struct *op, data_fmt *data){
	long i, numnodes, numalines, nummarkers, *oldsiteptr, *siteptr; 
	double offset;
	linlist *alines;
	brlist *brlines;
	tlist *tymeslice;
	boolean accept_change, *nodenumber;
	node *p;
	tree *oldtree = NULL;
	
	
	oldsiteptr = NULL;
	siteptr = data->siteptr;
	// test
	curtree->root->back->coal[0] = 0;
	
	//printf("1\n");
	//checktymelist();

	constructoldranges(curtree);
	
	//checkstruct();
	//	printf("2\n");
	//	checktymelist();
	// test
	if (oldsiteptr != NULL) free(oldsiteptr);
      
	//constructoldranges(curtree);
	//numnodes = 2 * oldtree->numcoals + 2;
	numnodes = curtree->numcoals + curtree->numrecombs + *data->dnaptr->numseq + 1;
	
	curtree->oldnodep = NULL;
	nodenumber = (boolean *)calloc(1,numnodes*sizeof(boolean));
	nummarkers = getdata_nummarkers(op,data);
	for(i = 0; i < numnodes; i++) nodenumber[i] = TRUE;
	if(op->datatype == 'n' || op->datatype == 's') {
		oldsiteptr = (long *)calloc(nummarkers,sizeof(long));
		memcpy(oldsiteptr,siteptr,nummarkers*sizeof(long));
	}
	/* pick a spot to cut */
	alines = NULL;
	
	p = addtip_S(op,data,curtree, &tymeslice,&offset);
	numnodes ++;
	
	tymeslice = findtymeslice(p);
		
	
	if (op->fc) addlinlistH(op,data,&alines,p,0);
	else addlinlistH(op,data,&alines,p,0);
	numalines = 1;
	
	brlines = NULL;
	
	//marknode(p);
	
	accept_change = TRUE;
	accept_change = growlineagesH(op,data,curtree,oldtree,tymeslice,&brlines,&alines,&numalines,offset,&nodenumber,&numnodes,siteptr);
	
	//	printf("3\n");
	//checktymelist();
	
	if (accept_change) {
		renumber_nodes_S(op,data,curtree);
		traverse_flagbelow(curtree,p->back);
		if (op->datatype == 'n' || op->datatype == 's')
			rebuild_alias(op,data,siteptr);
		//summerize tree
		//scorerecs_temp(curtree);
		//calculate P(G|parameters)
		//curtree->coalprob = treeprobS(curtree);
		//calculate P(D|G)
		
		//calcratio(op,data,oldtree,curtree,p->coal);     
	} 


	//checktymelist();
	
	free(nodenumber);
	freelinlist(alines);
	
	//accept_change = testratio(op,data,oldtree,curtree,'d'); 
	
	//for testing
	//accept_change = TRUE;
	
	
	curtree->recsummary = NULL;
	if (op->datatype == 'n' || op->datatype == 's')  free(oldsiteptr);
	return TRUE;
	
}





long countrec(tlist *start)
{
  long accum;
  tlist *t;

  accum = 0;
  for (t = start; t != NULL; t = t->succ)
    if (isrecomb(t->eventnode)) accum++;

  return(accum);

} /* countrec */

node *findrec(tlist *start, long target)
{
  long accum;
  tlist *t;

  accum = 0;
  for (t = start; t != NULL; t = t->succ) {
    if (isrecomb(t->eventnode)) accum++;
    if (accum == target) return(findunique(t->eventnode));
  }

  fprintf(ERRFILE,"ERROR:findrec--unable to find recombination #%ld\n",target);
  exit(-1);

} /* findrec */

boolean twiddle(option_struct *op, data_fmt *data)
{
  long i, numrecs, target, *oldsiteptr, numnodes, numalines, *siteptr;
  double activelinks;
  node *p;
  boolean accept_change, *nodenumber;
  tlist *tstart;
  brlist *brlines;
  linlist *alines;
  tree *oldtree;

  numrecs = countrec(curtree->tymelist);
  if (!numrecs) return(FALSE);

  oldtree = copytree(op,data,curtree);

  target = (long)(randum()*numrecs) + 1;
  p = findrec(curtree->tymelist,target);
  tstart = gettymenode(curtree,p->number);


  siteptr = data->siteptr;
  oldsiteptr = (long *)calloc(getdata_nummarkers(op,data),sizeof(long));
  memcpy(oldsiteptr,siteptr,getdata_nummarkers(op,data)*sizeof(long));

  /* now choose the new recombination point, adjusting for the fact that
     the active sites don't necessairly start at the beginning of the
     sequence */
  if (op->fc) {
    activelinks = count_activefc(op,data,p);
    activelinks = randum()*activelinks;
    target = choosepsitefc(op,data,p->back,activelinks);
  } else {
    activelinks = count_active(op,data,p);
    activelinks = randum()*activelinks;
    target = choosepsite(op,data,p->back,activelinks);
  }

 

  if (!p->back->ranges[0])
    fprintf(ERRFILE,"ERROR:twiddle chose a dead branch!\n");

  /* 50% chance which branch represents the "first" part of the sites */
  if (randum() > 0.5) {
    p->next->recstart = 0;
    p->next->recend = target;
    p->next->next->recstart = target + 1;
    p->next->next->recend = countsites(op,data) - 1;
  } else {
    p->next->recstart = target + 1;
    p->next->recend = countsites(op,data) - 1;
    p->next->next->recstart = 0;
    p->next->next->recend = target;
  }

  /* fix up the ranges and siteptr info, construct the brlist;
     startting drop_brifx_ranges at tstart->prev because a recombination
     should never be the eventnode of the first entry of the tymelist */
  brlines = NULL;
  drop_brfix(op,data,curtree,tstart->prev,&brlines);
  if (op->datatype == 'n' || op->datatype == 's')
    rebuild_alias(op,data,siteptr);

  /* now reevaluate the tree for added recombinations */
  alines = NULL;
  numalines = 0;
  numnodes = getdata_numtips(op,data) + counttymelist(curtree->tymelist);
  nodenumber = (boolean *)calloc(numnodes,sizeof(boolean));
  for(i = 0; i < numnodes; i++) nodenumber[i] = TRUE;

  accept_change = growlineages(op,data,curtree,oldtree,tstart,&brlines,
			       &alines,&numalines,0.0,&nodenumber,&numnodes,siteptr);

  if (accept_change) {
    renumber_nodes_S(op,data,curtree);
    remove_excess_tree(op,data,curtree,siteptr);
    renumber_nodes_S(op,data,curtree);
    /* root->back's range info may not be set; and must be hand set here. */
    if(curtree->root->back->ranges[0] != 1L)
      addrange(&curtree->root->back->ranges,0L,countsites(op,data)-1);
#if ALWAYS_ACCEPT
    curtree->likelihood = oldtree->likelihood;
    accept_change = testratio(op,data,oldtree,curtree,'t');
#endif
#if !ALWAYS_ACCEPT
    traverse_flagbelow(curtree,p->back);
    if (op->datatype == 'n' || op->datatype == 's')
      rebuild_alias(op,data,siteptr);
    localeval(op,data,curtree->root->back,FALSE);
    curtree->coalprob = coalprob(op,data,curtree,theta0,rec0);
    accept_change = testratio(op,data,oldtree,curtree,'t');
#endif
  } 

  /* if the change was accepted then nothing further needs to be done
     except for cleanup */
  if (!accept_change) {
    freetree(op,data,curtree);
    curtree = oldtree;
    memcpy(siteptr,oldsiteptr,getdata_nummarkers(op,data)*sizeof(long));
  } else freetree(op,data,oldtree);

  free(nodenumber);
  free(oldsiteptr);

  return(accept_change);

} /* twiddle */


/******************************************************************
 * isprohibited() returns TRUE if the "value" is contained in the *
 * passed "badstuff" and FALSE otherwise.                         */
boolean isprohibited(long value, long numbad, long *badstuff)
{
  long bad;

  for(bad = 0; bad < numbad; bad++)
    if (value == badstuff[bad]) return(TRUE);

  return(FALSE);

} /* isprohibited */


/*******************************************************************
 * isinvariant() returns TRUE if the site is invariant for a given *
 * creature, otherwise FALSE                                       */
boolean isinvariant(option_struct *op, data_fmt *data, long site,
		    creature *cr)
{
  char **dna;

  switch(op->datatype) {
  case 'a':
    break;
  case 'b':
    break;
  case 'm':
    break;
  case 'n':
  case 's':
    dna = data->dnaptr->seqs[population][locus];
    if (dna[cr->haplotypes[0]->number-1][site] ==
	dna[cr->haplotypes[1]->number-1][site])
      return(TRUE);
    break;
  }

  return(FALSE);

} /* isinvariant */


/**********************************************
 * choose_flipsite() returns a valid flipsite */
long choose_flipsite(option_struct *op, data_fmt *data, creature *cr,
		     long numprohibited, long *prohibited)
{
  long choice;

  do {
    choice = (long)(randum()*cr->numflipsites);
  } while (isprohibited(choice,numprohibited,prohibited));

  return(cr->flipsites[choice]);

} /* choose_flipsite */


/********************************************************************
 * datachangesite() flips a passed site between two elements of the *
 * actual data matrix.                                              */
void datachangesite(option_struct *op, data_fmt *data, node *p1, node *p2,
		    long site)
{
  char tempbase;

  tempbase = data->dnaptr->seqs[population][locus][p1->number-1][site];
  data->dnaptr->seqs[population][locus][p1->number-1][site] = 
    data->dnaptr->seqs[population][locus][p2->number-1][site];
  data->dnaptr->seqs[population][locus][p2->number-1][site] = tempbase;

} /* datachangesite */


/*****************************************************************
 * fliphap() flips NUMHAPFLIP sites among the haplotypes of a    *
 * randomly chosen creature of the passed in tree.  The sites to *
 * be flipped are chosen randomly amongst the eligible sites of  *
 * that creature.                                                */
boolean fliphap(option_struct *op, data_fmt *data, tree *tr)
{
  long i, j, flipsite, numcreatures, numslice, numcategs, datanumsites,
    *oldsiteptr, *flippedsites, *siteptr;
  double ***temp;
  creature *cr;
  tree *oldtree;
  boolean accept_change;

  datanumsites = getdata_nummarkers(op,data);
  siteptr = data->siteptr;

  oldtree = copytree(op,data,tr);
  oldsiteptr = (long *)calloc(datanumsites,sizeof(long));
  memcpy(oldsiteptr,siteptr,datanumsites*sizeof(long));

  numslice = (op->panel) ? NUMSLICE : 1L;
  numcategs = op->categs;

  temp = (double ***)calloc(op->categs,sizeof(double **));
  temp[0] = (double **)calloc(op->categs*numslice,sizeof(double *));
  for(i = 1; i < op->categs; i++) temp[i] = temp[0] + i*numslice;
  temp[0][0] = (double *)calloc(op->categs*numslice*4,sizeof(double));
  for(i = 0; i < op->categs; i++)
    for(j = 0; j < numslice; j++)
      temp[i][j] = temp[0][0] + i*numslice*4 + j*4;


  flippedsites = (long *)calloc(NUMHAPFLIP,sizeof(long));

  numcreatures = getdata_numtips(op,data)/NUMHAPLOTYPES;
  do {
    cr = &(tr->creatures[(long)(randum()*numcreatures)]);
  } while(cr->numflipsites == 0);

  /* debug DEBUG warning WARNING--this code is two haplotype specific */
  /* and it is also dna/snp specific!!!!                              */
  for(i = 0; i < NUMHAPFLIP; i++) {
    flippedsites[i] = flipsite = choose_flipsite(op,data,cr,i,flippedsites);
    /* first flip the x-arrays */
    memcpy(temp[0][0],cr->haplotypes[0]->x->s[flipsite][0][0],
	   numcategs*numslice*4*sizeof(double));
    memcpy(cr->haplotypes[0]->x->s[flipsite][0][0],
	   cr->haplotypes[1]->x->s[flipsite][0][0],
	   numcategs*numslice*4*sizeof(double));
    memcpy(cr->haplotypes[1]->x->s[flipsite][0][0],temp[0][0],
	   numcategs*numslice*4*sizeof(double));
    /* then flip the actual data */
    datachangesite(op,data,cr->haplotypes[0],cr->haplotypes[1],flipsite);
  }

  /* can just rebuild_alias() because flippable sites can no
     longer be aliased */
  if (op->datatype == 'n' || op->datatype == 's')
    rebuild_alias(op,data,siteptr);

  for(i = 0; i < NUMHAPLOTYPES; i++)
    traverse_flagbelow(tr,cr->haplotypes[i]->back);
  localeval(op,data,tr->root->back,FALSE);
  tr->coalprob = coalprob(op,data,tr,theta0,rec0);
  accept_change = testratio(op,data,oldtree,tr,'f');
#if ALWAYS_ACCEPT
  tr->likelihood = oldtree->likelihood;
  accept_change = testratio(op,data,oldtree,tr,'f');
#endif

  if (!accept_change) {
    for(i = 0; i < NUMHAPFLIP; i++)
      datachangesite(op,data,cr->haplotypes[1],cr->haplotypes[0],
		     flippedsites[i]);
    freetree(op,data,curtree);
    curtree = oldtree;
    memcpy(siteptr,oldsiteptr,datanumsites*sizeof(long));
  } else freetree(op,data,oldtree);

  free(oldsiteptr);
  free(flippedsites);
  free(temp[0][0]);
  free(temp[0]);
  free(temp);

  return(accept_change);

} /* fliphap */


/***********************************************************
 * isintymelist() returns TRUE if the nodelet "p" is among *
 * the nodelets in the tymelist; FALSE otherwise.          */
boolean isintymelist(option_struct *op, data_fmt *data,
		     tlist *start, node *p)
{
  tlist *t;

  for(t = start; t != NULL; t = t->succ) {
    if (p == t->eventnode) return(TRUE);
    if (p->next) {
      if (p->next == t->eventnode) return(TRUE);
      if (p->next->next)
	if (p->next->next == t->eventnode) return(TRUE);
    }
  }

  return(FALSE);

} /* isintymelist */


/**************************************************************
 * flipdrop_prunebr() removes from the brlist the bad entries *
 * created by serial-calling remove_futile().                 */
void flipdrop_prunebr(option_struct *op, data_fmt *data, tree *tr,
		      brlist **brlines)
{
  brlist *br, *brsucc;

  for(br = (*brlines); br != NULL; br = brsucc) {
    brsucc = br->succ;
    if (isintymelist(op,data,tr->tymelist,br->branchtop)) continue;
    br_remove(brlines,br);
  }

} /* flipdrop_prunebr */


/*****************************************************************
 * flipdrop() flips NUMHAPFLIP sites among the haplotypes of a   *
 * randomly chosen creature of the passed in tree, and then      *
 * does a drop on one or both affected tips.  The sites to       *
 * be flipped are chosen randomly amongst the eligible sites of  *
 * that creature.                                                */

/* WARNING:  the logic for dropping multiple lines here works    *
 * only for two tips, not in the general case!                   */

boolean flipdrop(option_struct *op, data_fmt *data, tree *tr, 
		 long numdrop)
{
  long i, j, flipsite, numcreatures, numslice, numcategs, datanumsites,
    *oldsiteptr, *flippedsites, *siteptr;
  double ***temp;
  creature *cr;
  tree *oldtree;
  linlist *alines;
  brlist *brlines;
  boolean accept_change, *nodenumber;
  node *p;
  long numnodes, numalines;

  datanumsites = getdata_nummarkers(op,data);
  siteptr = data->siteptr;

  oldtree = copytree(op,data,tr);
  oldsiteptr = (long *)calloc(datanumsites,sizeof(long));
  memcpy(oldsiteptr,siteptr,datanumsites*sizeof(long));
  numnodes = 2 * oldtree->numcoals + 2;
  nodenumber = (boolean *)calloc(1,numnodes*sizeof(boolean));
  for(i = 0; i < numnodes; i++) nodenumber[i] = TRUE;

  numslice = (op->panel) ? NUMSLICE : 1L;
  numcategs = op->categs;

  temp = (double ***)calloc(op->categs,sizeof(double **));
  temp[0] = (double **)calloc(op->categs*numslice,sizeof(double *));
  for(i = 1; i < op->categs; i++) temp[i] = temp[0] + i*numslice;
  temp[0][0] = (double *)calloc(op->categs*numslice*4,sizeof(double));
  for(i = 0; i < op->categs; i++)
    for(j = 0; j < numslice; j++)
      temp[i][j] = temp[0][0] + i*numslice*4 + j*4;


  flippedsites = (long *)calloc(NUMHAPFLIP,sizeof(long));

  numcreatures = getdata_numtips(op,data)/NUMHAPLOTYPES;
  do {
    cr = &(tr->creatures[(long)(randum()*numcreatures)]);
  } while(cr->numflipsites == 0);

  /* debug DEBUG warning WARNING--this code is two haplotype specific */
  /* and it is also dna/snp specific!!!!                              */
  for(i = 0; i < NUMHAPFLIP; i++) {
    flippedsites[i] = flipsite = choose_flipsite(op,data,cr,i,flippedsites);
    /* first flip the x-arrays */
    memcpy(temp[0][0],cr->haplotypes[0]->x->s[flipsite][0][0],
	   numcategs*numslice*4*sizeof(double));
    memcpy(cr->haplotypes[0]->x->s[flipsite][0][0],
	   cr->haplotypes[1]->x->s[flipsite][0][0],
	   numcategs*numslice*4*sizeof(double));
    memcpy(cr->haplotypes[1]->x->s[flipsite][0][0],temp[0][0],
	   numcategs*numslice*4*sizeof(double));
    /* then flip the actual data */
    datachangesite(op,data,cr->haplotypes[0],cr->haplotypes[1],flipsite);
  }

  /* now cut off the two tips */
  alines = NULL;
  brlines = NULL;
  numalines = 0;
  if (numdrop == 1) { /* drop one lineage at random */
    if (randum()>= 0.5) i=0;
    else i=1;
    p = cr->haplotypes[i];
    if (op->fc) addlinlist(&alines,p,count_activefc(op,data,p));
    else addlinlist(&alines,p,count_active(op,data,p));
    numalines++;
    tag_futile(p->back);
    traverse_flagbelow(curtree,p->back);
    remove_futile(op,data,curtree,oldtree,p,nodenumber,&brlines);
  } else { /* drop both lineages */
    for (i = 0; i < NUMHAPLOTYPES ; i++) { 
      if (randum()>= 0.5) i=0;
      else i=1;
      p = cr->haplotypes[i];
      if (op->fc) addlinlist(&alines,p,count_activefc(op,data,p));
      else addlinlist(&alines,p,count_active(op,data,p));
      numalines++;
      tag_futile(p->back);
      traverse_flagbelow(curtree,p->back);
      remove_futile(op,data,curtree,oldtree,p,nodenumber,&brlines);
    }
    flipdrop_prunebr(op,data,curtree,&brlines); 
  }

  accept_change = growlineages(op,data,curtree,oldtree,
			       curtree->tymelist,&brlines,&alines,&numalines,0L,&nodenumber,
			       &numnodes,siteptr);
  if (accept_change) {
    renumber_nodes(op,data,curtree);
    remove_excess_tree(op,data,curtree,siteptr);
    renumber_nodes(op,data,curtree);
    rebuild_alias(op,data,siteptr);
    for(i = 0; i < NUMHAPLOTYPES; i++) 
      traverse_flagbelow(curtree,cr->haplotypes[i]);
#if !ALWAYS_ACCEPT
    localeval(op,data,tr->root->back,FALSE);
    tr->coalprob = coalprob(op,data,tr,theta0,rec0);
#endif
#if ALWAYS_ACCEPT
    tr->likelihood = oldtree->likelihood;
#endif
    accept_change = testratio(op,data,oldtree,tr,'f');
  }

  free(nodenumber);
  freelinlist(alines);

  if (!accept_change) {
    for(i = 0; i < NUMHAPFLIP; i++)
      datachangesite(op,data,cr->haplotypes[1],cr->haplotypes[0],
		     flippedsites[i]);
    freetree(op,data,curtree);
    curtree = oldtree;
    memcpy(siteptr,oldsiteptr,datanumsites*sizeof(long));
  } else freetree(op,data,oldtree);

  free(oldsiteptr);
  free(flippedsites);
  free(temp[0][0]);
  free(temp[0]);
  free(temp);

  return(accept_change);

} /* flipdrop */


/*******************************************************************
 * addfractrecomb() adds the constant FRACTRECOMB to the number of *
 * recombinations scored for the last tree.                        *
 * This function should only be called at the end of a chain, and  *
 * then only if the trees scored during that chain have no         *
 * recombinations in any of them!                                  */
void addfractrecomb(option_struct *op, data_fmt *data, long chain,
		    treerec ***treessum)
{
  long refchain, chaintype;
  treerec *lasttree;

  refchain = REF_CHAIN(chain);
  chaintype = TYPE_CHAIN(chain);

  lasttree = &treessum[locus][refchain][op->numout[chaintype]-1];

  lasttree->numrecombs += FRACTRECOMB;

} /* addfractrecomb */


/***********************************************************************
 * add_one_recomb() runs a special rearrangement that always begins    *
 * with a recombination event.  The resulting tree is always accepted, *
 * and then is forcibly sampled (replacing the last tree of the chain).*/
void add_one_recomb(option_struct *op, data_fmt *data, tree *tr,
		    long *siteptr, long chain)
{
  long i, numalines, numnodes, *oldsiteptr = NULL, chaintype;
  double offset, newtyme, pc, pr;
  linlist *alines;
  brlist *brlines;
  tlist *tymeslice;
  node *cutpt;
  boolean *nodenumber, succeeded;
  tree *oldtree;

  oldtree = copytree(op,data,curtree);
  numnodes = 2 * curtree->numcoals + 2;
  nodenumber = (boolean *)calloc(numnodes,sizeof(boolean));
  for(i = 0; i < numnodes; i++) nodenumber[i] = TRUE;
  if(op->datatype == 'n' || op->datatype == 's') {
    oldsiteptr = (long *)calloc(getdata_nummarkers(op,data),sizeof(long));
    memcpy(oldsiteptr,siteptr,getdata_nummarkers(op,data)*sizeof(long));
  }
  alines = NULL;
  brlines = NULL;

  cutpt = tr->nodep[1];
  tymeslice = tr->tymelist;
  offset = 0.0;

  if (op->fc) addlinlist(&alines,cutpt,count_activefc(op,data,cutpt));
  else addlinlist(&alines,cutpt,count_active(op,data,cutpt));
  numalines = 1;

  tag_futile(cutpt->back);
  traverse_flagbelow(curtree,cutpt->back);
  remove_futile(op,data,curtree,oldtree,cutpt,nodenumber,&brlines);

  /* startoff with 1 recombination */
  /* generate the time the recombination will occur at */
  pc = 2.0 * countinactive(op,tymeslice,alines)/theta0;
  pr = rec0 * alines->activesites;
  newtyme = tymeslice->eventnode->tyme + offset -
    log(randum())/(pc + pr);

  /* find the tymeslice that corresponds with that time */
  for(; tymeslice->age < newtyme; tymeslice = tymeslice->succ)
    ;

  /* put in the recombination */
  succeeded = recomb(op,data,tymeslice,&alines,&numalines,newtyme,
		     siteptr,&nodenumber,&numnodes);
  tymeslice = tymeslice->succ;

  succeeded = growlineages(op,data,curtree,oldtree,tymeslice,&brlines,
			   &alines,&numalines,offset,&nodenumber,&numnodes,siteptr);


  free(nodenumber);
  freelinlist(alines);

  if (succeeded) {
    renumber_nodes(op,data,curtree);
    remove_excess_tree(op,data,curtree,siteptr);
    renumber_nodes(op,data,curtree);
    traverse_flagbelow(curtree,cutpt->back);
    if (op->datatype == 'n' || op->datatype == 's')
      rebuild_alias(op,data,siteptr);
    localeval(op,data,curtree->root->back,FALSE);
    curtree->coalprob = coalprob(op,data,curtree,theta0,rec0);

    /* now actually sample the generated tree */
    chaintype = TYPE_CHAIN(chain);
    op->numout[chaintype]--;

    scoretree(op,data,chain);
    op->numout[chaintype]++;
  }

  freetree(op,data,curtree);
  curtree = oldtree;
  if (op->datatype == 'n' || op->datatype == 's') {
    memcpy(siteptr,oldsiteptr,getdata_nummarkers(op,data)*sizeof(long));
    free(oldsiteptr);
  }


} /* add_one_recomb */













/*-----------------------------------------------------------------------*/






/*added for heterogenous recombination rate model */


boolean twiddleH(option_struct *op, data_fmt *data)
{
  long i, numrecs, target, *oldsiteptr, numnodes, numalines, *siteptr;
  //double activelinks;
  double weight;
  node *p;
  boolean accept_change, *nodenumber;
  tlist *tstart;
  brlist *brlines;
  linlist *alines;
  tree *oldtree;

  numrecs = countrec(curtree->tymelist);
  if (!numrecs) return(FALSE);

  oldtree = copytree(op,data,curtree);

  target = (long)(randum()*numrecs) + 1;
  p = findrec(curtree->tymelist,target);
  tstart = gettymenode(curtree,p->number);


  siteptr = data->siteptr;
  oldsiteptr = (long *)calloc(getdata_nummarkers(op,data),sizeof(long));
  memcpy(oldsiteptr,siteptr,getdata_nummarkers(op,data)*sizeof(long));


  // pick new recombination site

  if (op->fc) weight = randum() * count_activefc(op,data,p);
  else weight = randum() * count_activeH(op,data,p);

  for(target = p->ranges[1]; target < p->ranges[2*p->ranges[0]]; target++) {
    weight = weight - recarray[target];
    if (weight <= 0) break;
  }

  
  if (!p->back->ranges[0])
    fprintf(ERRFILE,"ERROR:twiddle chose a dead branch!\n");

  /* 50% chance which branch represents the "first" part of the sites */
  if (randum() > 0.5) {
    p->next->recstart = 0;
    p->next->recend = target;
    p->next->next->recstart = target + 1;
    p->next->next->recend = countsites(op,data) - 1;
  } else {
    p->next->recstart = target + 1;
    p->next->recend = countsites(op,data) - 1;
    p->next->next->recstart = 0;
    p->next->next->recend = target;
  }

  /* fix up the ranges and siteptr info, construct the brlist;
     startting drop_brifx_ranges at tstart->prev because a recombination
     should never be the eventnode of the first entry of the tymelist */
  brlines = NULL;
  drop_brfix(op,data,curtree,tstart->prev,&brlines);
  if (op->datatype == 'n' || op->datatype == 's')
    rebuild_alias(op,data,siteptr);

  /* now reevaluate the tree for added recombinations */
  alines = NULL;
  numalines = 0;
  numnodes = getdata_numtips(op,data) + counttymelist(curtree->tymelist);
  nodenumber = (boolean *)calloc(numnodes,sizeof(boolean));
  for(i = 0; i < numnodes; i++) nodenumber[i] = TRUE;

  accept_change = growlineages(op,data,curtree,oldtree,tstart,&brlines,
			       &alines,&numalines,0.0,&nodenumber,&numnodes,siteptr);

  if (accept_change) {
    renumber_nodes(op,data,curtree);
    remove_excess_tree(op,data,curtree,siteptr);
    renumber_nodes(op,data,curtree);
    /* root->back's range info may not be set; and must be hand set here. */
    if(curtree->root->back->ranges[0] != 1L)
      addrange(&curtree->root->back->ranges,0L,countsites(op,data)-1);

    curtree->likelihood = oldtree->likelihood; 
    accept_change = testratio(op,data,oldtree,curtree,'t');     
#if ALWAYS_ACCEPT
    curtree->likelihood = oldtree->likelihood;
    accept_change = testratio(op,data,oldtree,curtree,'t');
#endif
#if !ALWAYS_ACCEPT
    traverse_flagbelow(curtree,p->back);
    if (op->datatype == 'n' || op->datatype == 's')
      rebuild_alias(op,data,siteptr);
    //localeval(op,data,curtree->root->back,FALSE);
    //curtree->coalprob = coalprob(op,data,curtree,theta0,rec0);
    accept_change = testratio(op,data,oldtree,curtree,'t');
#endif
  } 

  /* if the change was accepted then nothing further needs to be done
     except for cleanup */
  if (!accept_change) {
    freetree(op,data,curtree);
    curtree = oldtree;
    memcpy(siteptr,oldsiteptr,getdata_nummarkers(op,data)*sizeof(long));
  } else freetree(op,data,oldtree);

  free(nodenumber);
  free(oldsiteptr);

  return(accept_change);

} /* twiddle */

double eventprobbrH_growth(option_struct *op, data_fmt *data, brlist **bigbr,
			   tlist *tint, double pc, double pr, double *pb, bseg *brs) // return only time of brlist recom
{
  double btyme, tyme, offset, brlength, stoptyme, bsitesum,bweight=0;
  brlist *finish, *curbr;

  bsitesum = 0.0;
  offset = 0.0;
  stoptyme = tint->age;
  tyme = BOGUSTREETYME;
  brlength = BOGUSTREETYME - 1.0;
  btyme = tint->eventnode->tyme;
  finish = NULL;
  (*pb) = 0.0;
  bweight = 0;
  curbr = *bigbr;
  while(curbr != NULL){
    bweight += curbr->weight;
    curbr = curbr->succ;
  }
  
  *pb = bweight;

  if (*pb) {
    brs->target = pickbrtargetH((*bigbr),bweight);
    brs->cross = pickbrcrossH(op,data,brs->target);
  }
  if (pc || pr || *pb) tyme = -log(randum())/((*pb));
  else tyme = btyme + 0.1;

  return(tyme + offset);

} 


double eventprobH_growth(option_struct *op, data_fmt *data, linlist *alines,
		  tlist *t, brlist **bigbr, char *event, bseg *brs)
{
  long numilines, numalines;
  double pc = 0, pr = 0, pb = 0, tyme_coal, tyme_rec, tyme_br, lasttyme, tyme;
  double g_time;
  linlist *u;

  if (!alines && (!*bigbr)) {(*event) = '0'; return(0.0);}

  /* number of inactive lineages */
  numilines = countinactive(op,t,alines);

  /* count the active lineages and sum active sites */
  numalines = 0;
  for (u = alines; u != NULL; u = u->succ) {
    numalines++;
    pr += u->weight;
  }
  lasttyme = t->eventnode->tyme;
  // rec time;
  if (pr != 0)
    tyme_rec = -log(randum())/(pr);
  else
    tyme_rec = 10000;

  // coal time
  pc = (numalines*(numalines-1) + 2.0*numalines*numilines);
  if (growth == 0){
    tyme_coal = -log(randum()) * theta0 /pc;
  }
  else{ // growth
    // 1. - alphag[pop]*size[pop]*exp(-alphag[pop]*(t - tlast[pop] ) )* log(rdum) / coal_prob
    // alphag - growth rate, size[pop] = theta0, t - current time (real time), tlast[pop] = 0
    g_time = 1.0 - growth * theta0 * exp(-growth * lasttyme) * log(randum()) / pc  ;
    tyme_coal = log(g_time) / growth;
  }
    

  if (!(*bigbr) || brs == NULL) {
    pb = -1; 
  }
  else{
      tyme_br = eventprobbrH_growth(op,data,bigbr,t,pc,pr,&pb,brs);
  }

  /* what kind of event? */
  if (pb < 0){
    if (tyme_coal < tyme_rec){
      tyme = tyme_coal;
      (*event) = 'c';
    }
    else{
      tyme = tyme_rec;
      (*event) = 'r';
    }
  }
  else{
    if (tyme_br < tyme_coal && tyme_br < tyme_rec){
      tyme = tyme_br;
      (*event) = 'b';
    }
    else{
      if (tyme_coal < tyme_rec){
	tyme = tyme_coal;
	(*event) = 'c';
      }
      else{
	tyme = tyme_rec;
	(*event) = 'r';
      }
    }
  }  
  
  if ((*event) != 'c' && (*event) != 'b' &&(*event) != 'r' ){
    printf("error");
  }
  return(tyme);

} 



double eventprobH(option_struct *op, data_fmt *data, linlist *alines,
		  tlist *t, brlist **bigbr, char *event, bseg *brs)
{
  long numilines, numalines;
  double pc = 0, pr = 0, pb = 0, tester, tyme;
  linlist *u;

  if (!alines && (!*bigbr)) {(*event) = '0'; return(0.0);}

  /* number of inactive lineages */
  numilines = countinactive(op,t,alines);

  /* count the active lineages and sum active sites */
  numalines = 0;
  for (u = alines; u != NULL; u = u->succ) {
    numalines++;
    pr += u->weight;
  }
  

  pc = (numalines*(numalines-1) + 2.0*numalines*numilines) / theta0;


  if (!(*bigbr) || brs == NULL) {
    pb = 0.0; 
    tyme = -log(randum())/(pc + pr);
  }
  else{
      tyme = eventprobbrH(op,data,bigbr,t,pc,pr,&pb,brs);
  }

  /* what kind of event? */
  if (!pc && !pr && !pb) {(*event) = '0'; return(tyme);}
  tester = randum();
  if (tester <= pr/(pc+pr+pb)) (*event) = 'r';
  else if (tester <= (pc+pr)/(pc+pr+pb)) (*event) = 'c';
  else (*event) = 'b';
  
  return(tyme);

} 

 
boolean makeeventH(option_struct *op, data_fmt *data, tree *growtree,tlist *tymeslice, brlist **brlines, linlist **alines, long *numalines,double newtyme, boolean **nodenumber, long *numnodes, long *siteptr, char eventtype, boolean rootdropped, bseg *brs, seglist *curseg)
{ 
  boolean succeeded;

  succeeded = TRUE;
  switch((int)eventtype) {
  case 'c' : 
    coalesceH(op,data,tymeslice,alines,numalines,newtyme,nodenumber,
	      numnodes,rootdropped,brlines,curseg);
    break;
  case 'r' :  
    succeeded = recombH(op,data,tymeslice,alines,numalines,newtyme,
			siteptr,nodenumber,numnodes,curtree,curseg);
    break;
  case 'b' :
    succeeded = brlistrecomb(op,data,tymeslice,alines,numalines,
			     brlines,newtyme,siteptr,nodenumber,numnodes,brs, growtree,curseg);
    growtree->numrecombs ++;
    break;
  case '0' :
    break;
  default :
    fprintf(ERRFILE,"ERROR:makeevent impossible case: can't get here\n");
    break;
  }

  return(succeeded);
} /* makeevent */


boolean recombH(option_struct *op, data_fmt *data, tlist *t,linlist **lineages, long *lines, double tyme, long *siteptr, boolean **nodenumber, long *numnodes, tree *tr, seglist *cursegs)
{
  long target = 0, numsites,i;
  double ntarget, numactive;
  double weighted_active = 0.0;
  linlist *u;
  node *d, *p;
  boolean rootpick;
  dnadata *dna;

  dna = data->dnaptr;
  numsites = countsites(op,data);

    

  u = *lineages;
  /* find splitting lineage, which must be active */
  /* we first count up all active sites on active lineages -- */
  numactive = 0;
  do {
    //numactive += u->activesites;
    //weighted_active += getweight(u->start,u->end);
    weighted_active += u->weight;
    u = u->succ;
  } while (u != NULL);


  ntarget = randum()*weighted_active;
  weighted_active = 0;
  u = *lineages;
  do {
    //weighted_active += getweight(u->start,u->end);
    weighted_active += u->weight;
    if (/*u->activesites != 0 && */weighted_active >= ntarget) break;
    u = u->succ;
  } while(1);
  d = u->branchtop;


  rootpick = FALSE;
  if (d == curtree->root->back) rootpick = TRUE;

  /* create recombination node */
  newnode(&p);
  p->number = setnodenumber(nodenumber,numnodes);
  p->next->number = p->number;
  p->next->next->number = p->number;
  p->top = FALSE;
  free_x(op,p);
  p->next->top = TRUE;
  //allocate_x(op,data,p->next);
  if(op->map)allocate_z(op,data,p->next);
  ranges_Malloc(p->next,TRUE,0L);
  p->next->next->top = TRUE;
  //allocate_x(op,data,p->next->next);
  if(op->map)allocate_z(op,data,p->next->next);
  ranges_Malloc(p->next->next,TRUE,0L);
  p->tyme = tyme;
  p->next->tyme = tyme;
  p->next->next->tyme = tyme;
  p->updated = FALSE;
  p->next->updated = FALSE;
  p->next->next->updated = FALSE;
  p->type = p->next->type = p->next->next->type = 'r';

  ranges_Malloc(p,TRUE,0L);

  curtree->numrecombs += 1;
  curtree->nodep[p->number] = p->next;

  hookup(p,d);
  fixlength(op,data,p);


  if (rootpick) {
    hookup(p->next,curtree->root);
    fixlength(op,data,p->next);
  }

  /* partition the sites */
  /*   numactive = randum()*u->activesites; */
  /*   numactive = getsites(u->start,u->end) - u->start; */

  weighted_active = randum() * u->weight;
  for (i=u->start;i<u->end;i++){
    weighted_active = weighted_active - recarray[i];
    if (weighted_active<0){
      target = i;
      break;
    }
  }

  newgrec[i]++;

  //printf("rec target %ld\n",target);

  if(randum()>0.5) {
    p->next->recstart = 0;
    p->next->recend = target;
    p->next->next->recstart = target+1;
    p->next->next->recend = numsites-1;
  } else {
    p->next->recstart = target+1;
    p->next->recend = numsites-1;
    p->next->next->recstart = 0;
    p->next->next->recend = target;
  }

  /* set the ranges fields for the 2 rootwards branches */
  contrib(op,data,p->next,&p->next->ranges);
  contrib(op,data,p->next->next,&p->next->next->ranges);

  //contrib_temp(p,&p->ranges);

  contrib_temp(p,&p->coal);
  contrib_temp(p->next,&p->next->coal);
  contrib_temp(p->next->next, &p->next->next->coal);

  insertaftertymelist(t,p);
  //extbranch(curtree,t->succ,p->next->next);
  //extbranch(curtree,t->succ,p->next);

  /* set the coal fields for the 2 rootwards branches */
  if (op->fc) fix_coal(op,data,curtree,t->succ,p);

  /* fix alias entries */
  if (op->datatype == 'n' || op->datatype == 's')
    edit_alias(op,data,siteptr,target+1);

  /* remove old entry from linlist */
  sublinlist(lineages,d);

  //test
  modifynode(p,t,cursegs);


  /* add two new entries to linlist */
  if (op->fc) {
    addlinlistH(op,data,lineages,p->next,0);
    addlinlistH(op,data,lineages,p->next->next,0);
  } else {
    addlinlistH(op,data,lineages,p->next,0);
    addlinlistH(op,data,lineages,p->next->next,0);
  }
  (*lines)++;

  //if (curtree->numrecombs > RECOMB_MAX) return(FALSE);
  return(TRUE);
} /* recomb */


/* getsites - pick recombination site between start and end */

long getsites(long start, long end){
  long i;
  float sumweight = 0.0;
  float weight = 0.0;
  
  sumweight = getweight(start,end);
  sumweight = randum() * sumweight;

  for (i = start; i<=end ; i++){
    weight += recarray[i];
    if(weight>sumweight) return i;
  }
  return 0;
}

void addlinlistH(option_struct *op, data_fmt *data,linlist **t, node *target, double activelinks){
  linlist *r, *s, *u;

  s = *t;
  if (s == NULL) {
    newlin(&s);   
    s->branchtop = target;
    s->activesites = activelinks;
    s->succ = NULL;
    s->prev = NULL;
    
    s->start = -1;
    s->end = -1;
/*     if (target->ranges[0]!=0){ */
/*       s->start = target->ranges[1]; */
/*       s->end = target->ranges[2 * target->ranges[0]]; */
/*     } */
    
    if (target->coal[0]!=0){
      s->start = target->coal[1];
      s->end = target->coal[2 * target->coal[0]];
    }
    if ((s->end - s->start) != count_activefc(op,data,target)) printf("DIFFERENT NUMBER1\n");
    s->weight = getweight(s->start,s->end);
    (*t)=s;
  } 
  else {
    newlin(&r);
    r->branchtop = target;
    r->activesites = activelinks;

    r->start = -1;
    r->end = -1;
/*     if (target->ranges[0]!=0){ */
/*       r->start = target->ranges[1]; */
/*       r->end = target->ranges[2 * target->ranges[0]]; */
/*     } */

    if (target->coal[0]!=0){
      r->start = target->coal[1];
      r->end = target->coal[2 * target->coal[0]];
    }
    r->weight = getweight(r->start,r->end);

    /* insert new entry right after s */
    u = s->succ;
    s->succ = r;
    r->prev = s;
    r->succ = u;
    if ((r->end - r->start) != count_activefc(op,data,target)) printf("DIFFERENT NUMBER2\n");
   
    if(u!=NULL) u->prev = r;


  }
} /* addlinlist */

void coalesceH(option_struct *op, data_fmt *data, tlist *t,linlist **lineages, long *lines, double tyme, boolean **nodenumber,long *numnodes, boolean rootdropping, brlist **br, seglist *cursegs)
{
  long numalin, numilin, *active1, *active2;
  double Pbactive;
  node *p, *daughter1, *daughter2;
  boolean bothactive, rootpick;
  dnadata *dna;
  brlist *newbr;

  tlist *tcopy;

  dna = data->dnaptr;

  /* find the number of active (numalin) and inactive (numilin) lineages,
     when rootdropping, then the last inactive lineage (the root)
     becomes active; leaving no inactive lineages. */

  tcopy = t;
  while (1){
    if (t->update != 1) break;
    else t = t->prev;
  }

  numalin = *lines;
  if (!rootdropping){
    numilin = countinactive(op,t,*lineages);
  }
  else numilin = 0;
  

  /* pick the 2 lineages to coalesce, daughter1 and daughter2
     Pbactive = the chance that they are both active */
  Pbactive = (numalin)*(numalin-1) / ((numalin)*(numalin-1) + 2.0*numalin*numilin);

  if (Pbactive >= randum()) {
    daughter1 = pickactive(numalin,lineages,NULL);
    daughter2 = pickactive(numalin,lineages,daughter1);
    bothactive = TRUE;
  } else {
    daughter1 = pickactive(numalin,lineages,NULL);
    daughter2 = pickinactive(op,t,daughter1,*lineages);
    bothactive = FALSE;
  }
  t = tcopy;


  rootpick = FALSE;
  if (daughter1 == curtree->root->back || daughter2 == curtree->root->back)
    rootpick = TRUE;

  /* coalesce the two partners */
  newnode(&p);
  p->number = setnodenumber(nodenumber,numnodes);
  p->next->number = p->number;
  p->next->next->number = p->number;
  p->top = FALSE;
  free_x(op,p);
  p->next->top = FALSE;
  free_x(op,p->next);
  p->next->next->top = TRUE;
  //allocate_x(op,data,p->next->next);
  if (op->map) allocate_z(op,data,p->next->next);
  ranges_Malloc(p->next->next,TRUE,0L);
  p->tyme = tyme;
  p->next->tyme = tyme;
  p->next->next->tyme = tyme;
  p->updated = FALSE;
  p->next->updated = FALSE;
  p->next->next->updated = FALSE;
  p->type = p->next->type = p->next->next->type = 'c';

  curtree->numcoals += 1;
  curtree->nodep[p->number] = p->next->next;

  if(!bothactive) { 
    hookup(p->next->next,daughter2->back);
    fixlength(op,data,p->next->next); /* but only if not bothactive! */
  }
  hookup(p,daughter1);
  fixlength(op,data,p);
  hookup(p->next,daughter2);
  fixlength(op,data,p->next);
  if (rootpick) {
    hookup(p->next->next,curtree->root);
    fixlength(op,data,p->next->next);
  } 
  

  /* now set the ranges for the new node */
  contrib(op,data,p->next->next,&p->next->next->ranges);

  insertaftertymelist(t,p);

  sublinlist(lineages,daughter1);

  //new fc code	
    
  updatecoalpart(p);

  // handle coal stuff.
  modifynode(p, t, cursegs);

  if (p->next->next->coal[0] == 0){
    t->succ->numbranch --;
    t->succ->branchlist[0] = t->succ->branchlist[t->succ->numbranch];
    t->succ->branchlist[t->succ->numbranch] = findunique(p);
  }
    
  if (bothactive) {    
    if (p->next->next->coal[0] != 0) addlinlistH(op,data,lineages,p->next->next,0);
    else (*lines)--;
	  
    sublinlist(lineages,daughter2);
  }
  (*lines)--;  
    
  if(!bothactive) {
    /* rename the changed branchtop segments */
    //drop_renamebrlist(br,daughter2,p->next->next);
    /* then fix the rest of the ranges and all brlist segments */
    if (daughter2->branch != NULL) {
      br_remove(br,daughter2->branch);
      daughter2->branch = NULL;
    }
    active1 = p->next->next->coal;
    active2 = daughter2->coal;
    
    if (p->next->next->coal[0] !=0){
          
      if (p->next->next->back->oldbackranges[1] > p->next->next->coal[p->next->next->coal[0] * 2] || p->next->next->back->oldbackranges[p->next->next->back->oldbackranges[0] * 2] < p->next->next->coal[1]){
	addlinlistH(op,data,lineages,p->next->next,0);
	p->next->next->back->update = 2;
	(*lines)++;
	return;
      }
      else{    
	if (!sameranges(p->next->next->back->oldbackranges,p->next->next->coal)){
	  newbr = initbrlist(p->next->next, p->next->next->tyme,p->next->next->back->tyme,countsites(op,data));
	  updatebrs(newbr, p->next->next->back->oldbackranges, p->next->next->coal);
	  //updatesegranges(op,data,newbr,NULL,p->next->next->back->oldbackranges,p->next->next->coal);
	  newbr->weight = count_weightbr(op,data,newbr);
	  if (newbr->weight > 0){
	    hookup_brlist(br,newbr);
	    p->next->next->branch = newbr;
	  }
	  else freebrlist(newbr);
	  
	}
      }
    }
    else p->next->next->back = curtree->root;
  }
    
} /* coalesce */


long findrecstart(option_struct *op, data_fmt *data, node *p){
  long i, newstart;
  long *range1, *range2;
  node *q;

  if (istip(p)) return 0;
  
  if (iscoal(p)) {
    range1 = p->next->back->ranges;
    range2 = p->next->next->back->ranges;

    if (!range1[0] && !range2[0])
      fprintf(ERRFILE, "ERROR:Dead branch in count_coal_active! - findrecstart %ld %ld\n",indecks,apps);

    if (!range1[0]){
      newstart = range2[1];
      return newstart;
    }
    if (!range2[0]){
      newstart = range1[1];
      return newstart;
    }
    newstart = (range1[1] < range2[1]) ? range1[1] : range2[1];
    return newstart;
  }
  q = findunique(p)->back;

  if (!q->ranges[0])
    fprintf(ERRFILE,"ERROR:count_rec_active counted dead 1, %ld %ld\n",indecks,apps);
  i = count_rec_active(op,data,p);
  for(i = 1, newstart = -1L; q->ranges[i] != FLAGLONG; i+=2) {
    if(p->recstart > q->ranges[i+1]) continue;
    newstart = (q->ranges[i] > p->recstart) ? q->ranges[i] : p->recstart;
    break;
  }

  return newstart;

} 


long findrecend(option_struct *op, data_fmt *data, node *p)
{
  long i, newend;
  long *range1, *range2;
  node *q;

  if (istip(p)) return *(data->dnaptr->sites) - 1;

  if (iscoal(p)){
    range1 = p->next->back->ranges;
    range2 = p->next->next->back->ranges;
    if (!range1[0] && !range2[0])
      fprintf(ERRFILE, "ERROR:Dead branch in count_coal_active! findrecend %ld %ld\n",indecks,apps);
      
    if (!range1[0]){
      newend = range2[2*range2[0]];
      return newend;
    }

    if (!range2[0]){
      newend = range1[2*range1[0]];
      return newend;
    }
    newend = (range1[2*range1[0]] > range2[2*range2[0]]) ? range1[2*range1[0]] : range2[2*range2[0]];
    return newend;
  }



  q = findunique(p)->back;
  for(i=2*q->ranges[0], newend = -1L; i != 0; i-=2) {
    if(p->recend < q->ranges[i-1]) continue;
    newend = (q->ranges[i] < p->recend) ? q->ranges[i] : p->recend;
    break;
  }

  if(newend == -1L) return(0L);

  return newend;

} 

double eventprobbrH(option_struct *op, data_fmt *data, brlist **bigbr,
		    tlist *tint, double pc, double pr, double *pb, bseg *brs)
{
  double btyme, tyme, offset, brlength, stoptyme, bsitesum,bweight=0;
  brlist *finish, *curbr;

  bsitesum = 0.0;
  offset = 0.0;
  stoptyme = tint->age;
  tyme = BOGUSTREETYME;
  brlength = BOGUSTREETYME - 1.0;
  btyme = tint->eventnode->tyme;
  finish = NULL;
  (*pb) = 0.0;
  bweight = 0;
  curbr = *bigbr;
  while(curbr != NULL){
    bweight += curbr->weight;
    curbr = curbr->succ;
  }
  
  *pb = bweight;

  if (*pb) {
    brs->target = pickbrtargetH((*bigbr),bweight);
    brs->cross = pickbrcrossH(op,data,brs->target);
  }
  if (pc || pr || *pb) tyme = -log(randum())/(pc + pr + (*pb));
  else tyme = btyme + 0.1;

  return(tyme + offset);

} 

brlist *pickbrtargetH(brlist *start, long weight)
{
  double pick;
  brlist *br;
  
  br = start;
  pick = randum() * weight;

  while (1){
    pick = pick - br->weight;
    if (pick < 0) return br;
    br = br->succ;
  }
} 

                                                 
long pickbrcrossH(option_struct *op, data_fmt *data, brlist *source)
{
  long cross;
  double weight;
  boolean found;

  weight = randum() * source->weight;
  found = FALSE;
  
  
  for(cross = source->nfsite; cross < source->ofsite;  cross++) {
    weight = weight - recarray[cross];
    if (weight < 0) return cross;
  }
  for(cross = source->olsite; cross < source->nlsite; cross++) {
      weight = weight - recarray[cross];
      if (weight < 0) return cross;
  }
  printf("error at picking site\n");
  return -1;

} 

void init_numnewactiveH(option_struct *op, data_fmt *data, tree *tr,
			brlist **bigbr, node *branchtop)
{
  brlist *br, *brsucc;

  for(br = (*bigbr); br != NULL; br = brsucc) {
    brsucc = br->succ;
    if (branchtop->tyme >= br->endtyme) {br_remove(bigbr,br); continue;}
    if (branchtop->tyme > br->starttyme) {
      br->starttyme = branchtop->tyme;
      if (br->starttyme >= br->endtyme) {br_remove(bigbr,br); continue;}
    }
    br->numnewactives = count_activebr(op,data,br);
    br->weight = count_weightbr(op,data,br);
  }

} 
                                   
double count_weightbr(option_struct *op, data_fmt *data, brlist *br)
{

    return(getweight(br->nfsite,br->ofsite) + getweight(br->olsite,br->nlsite));
    

} 



float count_weighted_tlist(option_struct *op, data_fmt *data, tlist *t)
{
  long i;  
  float weighted; // = sum of rec_rates * sites

  weighted = 0.0;
  for (i = 0; i < t->numbranch; i++) 
    if (op->fc) weighted += count_activefc(op,data,t->branchlist[i]->back);
    else weighted += count_activeH(op,data,t->branchlist[i]->back);

  return (weighted);

}    
  

float count_activeH(option_struct *op, data_fmt *data, node *p)
{
  double value;
  node *q;

  q = p;
  if (!q->top) q = q->back;

  if (istip(q)) 
    return(getweightedlinks(op,data,0L,countsites(op,data)-1));
  
  if (isrecomb(q)) value = count_rec_activeH(op,data,q);
  else value = count_coal_activeH(op,data,q);

  return(value);

} /* count_active */

double getweight(long start, long end){
  long i;
  double weight = 0.0;
  for (i = start; i<end ; i++){
    weight += recarray[i];
  }
  return weight;
}
  

float getweightedlinks(option_struct *op, data_fmt *data, long start, long end) // just for tips? - no ,just for nucleotide -for now
{
  long i;
  double sum, *spaces;
  float weighted;

  spaces = NULL;

  switch(op->datatype) {
  case 'a':
    break;
  case 'b':
  case 'm':
    spaces = data->msptr->mspace[population][locus];
    break;
  case 'n':
    weighted = getweight(start,end);
    return(weighted);
    break;
  case 's': // SNP data?
    spaces = data->dnaptr->sspace[population][locus];
    for(i = start, sum = 0.0; i < end; i++)
      sum += spaces[i];
    return(sum);
    break;
  default:
    fprintf(ERRFILE,"unknown datatype in getnumlinks\n");
    break;
  }

  return((double)FLAGLONG);

}


float count_coal_activeH(option_struct *op, data_fmt *data, node *p)
{
  long *range1, *range2, newstart, newend;

  range1 = p->next->back->ranges;
  range2 = p->next->next->back->ranges;

  if (!range1[0] && !range2[0])
    fprintf(ERRFILE, "ERROR:Dead branch in count_coal_active! - count_coal_activeH %ld %ld\n",
	    indecks,apps);

  if (!range1[0]) return(getweightedlinks(op,data,range2[1],range2[2*range2[0]]));
  if (!range2[0]) return(getweightedlinks(op,data,range1[1],range1[2*range1[0]]));

  newstart = (range1[1] < range2[1]) ? range1[1] : range2[1];
  newend = (range1[2*range1[0]] > range2[2*range2[0]]) ?
    range1[2*range1[0]] : range2[2*range2[0]];

  if (newend - newstart < 0)
    fprintf(ERRFILE,"ERROR:negative link count in count_coal_active! %ld %ld\n",
	    indecks,apps);

  return(getweightedlinks(op,data,newstart,newend));

} 


float count_rec_activeH(option_struct *op, data_fmt *data, node *p)
{
  long i, newstart, newend;
  node *q;

  q = findunique(p)->back;

  if (!q->coal[0])
    fprintf(ERRFILE,"ERROR:count_rec_active counted dead 2, %ld %ld\n",
	    indecks,apps);

  for(i = 1, newstart = -1L; q->coal[i] != FLAGLONG; i+=2) {
    if(p->recstart > q->coal[i+1]) continue;
    newstart = (q->coal[i] > p->recstart) ? q->coal[i] : p->recstart;
    break;
  }

  if(newstart == -1L) return(0L);

  for(i=2*q->coal[0], newend = -1L; i != 0; i-=2) {
    if(p->recend < q->coal[i-1]) continue;
    newend = (q->coal[i] < p->recend) ? q->coal[i] : p->recend;
    break;
  }

  if(newend == -1L) return(0L);

  if (newend - newstart < 0)
    fprintf(ERRFILE,"ERROR:negative link count in count_rec_active! %ld %ld\n",
	    indecks,apps);

  return(getweightedlinks(op,data,newstart,newend));

} /* count_rec_active */
  






long findcoalpos(option_struct *op, data_fmt *data, node *p,int key){
  long *range1, *range2, newstart, newend;

  range1 = p->next->back->coal;
  range2 = p->next->next->back->coal;

  if (!range1[0] && !range2[0])
    fprintf(ERRFILE, "ERROR:Dead branch in count_coal_active! findcoalpos %ld %ld\n",
	    indecks,apps);

  if (!range1[0]){
    newstart = range2[1];
    newend = range2[2*range2[0]];
    if (key == 0) return newstart;
    return newend;
  }

  if (!range2[0]){
    newstart = range1[1];
    newend = range1[2*range2[0]];
    if (key == 0) return newstart;
    return newend;
  }
  newstart = (range1[1] < range2[1]) ? range1[1] : range2[1];
  newend = (range1[2*range1[0]] > range2[2*range2[0]]) ?
    range1[2*range1[0]] : range2[2*range2[0]];

  if (key == 0) return newstart;
  return newend;

}

boolean rootdropH(option_struct *op, data_fmt *data, tlist *t, linlist **lineages, long *lines, long *siteptr, boolean **nodenumber, long *numnodes, brlist **br, seglist *curseg, tlist *activet, double basetyme){
  tlist *temptyme;
  double /*basetyme,*/ newtyme;
  char eventtype;
  boolean not_aborted;
bseg *nobr = NULL;	
  node *p;
  long i;

  temptyme = activet;
  
  if (activet->eventnode->next->update == 2 || activet->eventnode->next->next->update == 2) activet->update = 1;

  while (1){
    if (activet->update != 1) break;
    activet = activet->prev;
  }


  /* if necessary, add the other nodes to the active lineages */
  for (i = 0; i < activet->numbranch; i++){
    if (!isactive(activet->branchlist[i], *lineages) && activet->branchlist[i]->coal[0] != 0){
      addlinlistH(op,data,lineages,activet->branchlist[i],0);
      (*lines)++;
    }
  }
  activet = temptyme;

  //basetyme = temptyme->eventnode->tyme;
  p = curtree->root->back;
  not_aborted = TRUE;
  delete_brlist(br);
	

  while (*lines > 1) {
    newtyme = basetyme + eventprobH_growth(op,data,*lineages,activet,br,&eventtype,nobr);

    updateweights(lineages,br,newtyme - basetyme);

    switch ((int)eventtype) {
    case 'c' :
      coalesceH(op,data,t,lineages,lines,newtyme,nodenumber,numnodes,TRUE,br, curseg);
      break;
    case 'r' :
      not_aborted = recombH(op,data,t,lineages,lines, newtyme,siteptr,nodenumber,numnodes,curtree, curseg);
      break;
    case 'b' :
      fprintf(ERRFILE,"ERROR:rootdrop case b: can't get here\n");
					break;
    case '0' :
      fprintf(ERRFILE,"ERROR:rootdrop case 0: can't get here\n");
					break;
    default :
      fprintf(ERRFILE,"ERROR:rootdrop case: can't get here\n");
      break;
    }
    //printf("%ld %c %d\n",*lines,eventtype,(int)not_aborted);
    //printf("eventype %c node %ld\n", eventtype, t->succ->eventnode->id);
    if (!not_aborted) break;
    t = t->succ;
    activet = t;
    curseg = t->segments;
    basetyme = newtyme;
  }
  //if (not_aborted) traverse_flagbelow(curtree,p);
  return(not_aborted);
}

//addbrlist - adding "weight" part
void addbrlistH(option_struct *op, data_fmt *data, brlist **bigbr,
		node *newtop, long *oldranges, long *newranges, int *segranges,
		double starttyme, double endtyme, boolean building)
{
  double oldendtyme;
  int *tsegranges, *oldsegs;
  brlist *tnewbr, *bnewbr, *newbr;
  long numsites;

  numsites = countsites(op,data);
  if (starttyme == endtyme || !newtop) return;

  tnewbr = getbrlistbystart(*bigbr,newtop,starttyme);
  bnewbr = getbrlistbyend(*bigbr,newtop,endtyme);
  if (tnewbr) tnewbr->weight = count_weightbr(op,data,tnewbr);
  if (bnewbr) bnewbr->weight = count_weightbr(op,data,bnewbr);

  if (!tnewbr && !bnewbr) {
    newbr = initbrlist(newtop,starttyme,endtyme,countsites(op,data));
    if (oldranges && newranges) {
      updatesegranges(op,data,newbr,NULL,oldranges,newranges);
    } else {
      updatesegranges(op,data,newbr,segranges,NULL,NULL);
    }
    newbr->weight = count_weightbr(op,data,newbr);
    if (newbr->weight > 0){
      hookup_brlist(bigbr,newbr);
    }
    else freebrlist(newbr);
    return;
  }

  /* we have the same start & end */
  if (tnewbr && bnewbr) {
    if (tnewbr->endtyme == endtyme) { /* we have the same segment */
      if (newranges) {
	if (building) {
	  updatesegranges(op,data,tnewbr,NULL,oldranges,newranges);
	} 
	else {
	  updatesegranges(op,data,tnewbr,tnewbr->segranges,NULL,newranges);
	}
      } 
      else {
	updatesegranges(op,data,tnewbr,segranges,NULL,NULL);
      }
    }
    else {
      if (newranges) {
	for(; ; tnewbr = tnewbr->succ) {
	  if (tnewbr->branchtop != bnewbr->branchtop) continue;
	  if (building) {
	    updatesegranges(op,data,tnewbr,NULL,oldranges,newranges);
	  }
	  else {
	    updatesegranges(op,data,tnewbr,tnewbr->segranges,NULL,newranges);
	  }
	  if (tnewbr == bnewbr) break;
	  if (tnewbr->endtyme != bnewbr->starttyme)
	    addbrlistH(op,data,bigbr,newtop,oldranges,newranges,segranges,tnewbr->endtyme,bnewbr->starttyme,building);
	}
      }
    }
    return;
  }

  /* we only have the same start */
  if (tnewbr) {
    oldendtyme = tnewbr->endtyme;
    if (oldendtyme > endtyme) {
      tnewbr->starttyme = endtyme;
      addbrlistH(op,data,bigbr,newtop,oldranges,newranges,segranges,
		 starttyme,endtyme,building);
    }
    else {
      updatesegranges(op,data,tnewbr,tnewbr->segranges,NULL,newranges);
      addbrlistH(op,data,bigbr,newtop,oldranges,newranges,segranges,
		 oldendtyme,endtyme,building);
    }
    return;
  }
  else { /* we only have the same end */
    if (starttyme > bnewbr->starttyme) {
      bnewbr->endtyme = starttyme;
      oldsegs = (int *)calloc(numsites,sizeof(int));
      memcpy(oldsegs,bnewbr->segranges,numsites*sizeof(int));
      updatesegranges(op,data,bnewbr,bnewbr->segranges,NULL,newranges);
      tsegranges = (int *)calloc(numsites,sizeof(int));
      memcpy(tsegranges,bnewbr->segranges, numsites*sizeof(int));
      memcpy(bnewbr->segranges,oldsegs,numsites*sizeof(int));
      addbrlistH(op,data,bigbr,newtop,NULL,NULL,tsegranges,starttyme,
		 endtyme,building);
      free(oldsegs);
      free(tsegranges);
    } 
    else {
      addbrlistH(op,data,bigbr,newtop,oldranges,newranges,segranges,
		 starttyme,bnewbr->starttyme,building);
      updatesegranges(op,data,bnewbr,bnewbr->segranges,NULL,newranges);
    }
    return;
  }

} /* addbrlistH */

boolean makedropH(option_struct *op, data_fmt *data){
  long i, numnodes, numalines, nummarkers, *oldsiteptr, *siteptr; 
  double offset;
  linlist *alines;
  brlist *brlines;
  tlist *tymeslice;
  boolean accept_change, *nodenumber;
  node *p;
  tree *oldtree = NULL;
  

  oldsiteptr = NULL;
  siteptr = data->siteptr;
  // test
  curtree->root->back->coal[0] = 0;

  constructoldranges(curtree);

  //checkstruct();

  /* copy the original tree configuration */
  if (oldtree != NULL) freetree(op,data,oldtree);
  oldtree = copytree(op,data,curtree);
    
  // test
  if (oldsiteptr != NULL) free(oldsiteptr);
    
  //constructoldranges(curtree);
  //numnodes = 2 * oldtree->numcoals + 2;
  numnodes = curtree->numcoals + curtree->numrecombs + *data->dnaptr->numseq + 1;

  curtree->oldnodep = NULL;
  nodenumber = (boolean *)calloc(1,numnodes*sizeof(boolean));
  nummarkers = getdata_nummarkers(op,data);
  for(i = 0; i < numnodes; i++) nodenumber[i] = TRUE;
  if(op->datatype == 'n' || op->datatype == 's') {
    oldsiteptr = (long *)calloc(nummarkers,sizeof(long));
    memcpy(oldsiteptr,siteptr,nummarkers*sizeof(long));
  }
  /* pick a spot to cut */
  alines = NULL;
  
  while (1){
    p = pickbranch(op,data,curtree,&tymeslice,&offset);
    if (p->coal[0] != 0) break;
  }
  
  tymeslice = findtymeslice(p);
  
  if (op->fc) addlinlistH(op,data,&alines,p,0);
  else addlinlistH(op,data,&alines,p,0);
  numalines = 1;
  
  brlines = NULL;
  
  marknode(p);
  
  accept_change = TRUE;
  accept_change = growlineagesH(op,data,curtree,oldtree,tymeslice,&brlines,&alines,&numalines,offset,&nodenumber,&numnodes,siteptr);
  //accept_change = TRUE;
  
  if (accept_change) {
    renumber_nodes(op,data,curtree);
    traverse_flagbelow(curtree,p->back);
    if (op->datatype == 'n' || op->datatype == 's')
      rebuild_alias(op,data,siteptr);
    //summerize tree
    //scorerecs_temp(curtree);
    //calculate P(G|parameters)
    //if (marker1pos == -99 && op->ctemp != 1) curtree->coalprob = treeprob(curtree);
      //calculate P(D|G)
    
    calcratio(op,data,oldtree,curtree,p->coal);     
  } 

  free(nodenumber);
  freelinlist(alines);

  accept_change = testratio(op,data,oldtree,curtree,'d'); 

  //for testing
  //accept_change = TRUE;

  if (accept_change){
    curtree->recsummary = NULL;
    freerecinfo(oldtree->recsummary);
    freepath(oldtree->recpath);
    freetree(op,data,oldtree);
    if (op->datatype == 'n' || op->datatype == 's')  free(oldsiteptr);
    return TRUE;
  }
  else{
    freerecinfo(curtree->recsummary);
    freetree(op,data,curtree);
    curtree = oldtree;
    if (op->datatype == 'n' || op->datatype == 's') {
      memcpy(siteptr,oldsiteptr,nummarkers*sizeof(long));
      free(oldsiteptr);
    }
    return FALSE;
  }
}

    
struct seglist *copyseglist(seglist *original){ // copy segments
  seglist *curseg = NULL,*oldseg, *newlist,*tempseg;
    
  newlist = (struct seglist *)malloc(sizeof (struct seglist));
  newlist->start = original->start;
  newlist->end = original->end;
  newlist->numsam = original->numsam;
  newlist->prev = NULL;
  newlist->next = NULL;
    
  if (original->next != NULL){
    oldseg = newlist;
    for (tempseg = original->next; tempseg != NULL; tempseg = tempseg->next){
      curseg = (struct seglist *)malloc(sizeof (struct seglist));
      curseg->start = tempseg->start;
      curseg->end = tempseg->end;
      curseg->numsam = tempseg->numsam;
      curseg->prev = oldseg;
	    
      oldseg->next = curseg;
      oldseg = curseg;
    }
    curseg->next = NULL;
  }
  return newlist;
} 
    
void updateranges(tree *tr, node *eventnode){ // ONLY mark connected nodes     

  //tyme = gettymenode(tr, eventnode->back->number);

  if (!isrecomb(eventnode)){
	
    if (eventnode->coal[0] == 0){ // all segments are finally coalesced
      eventnode->back->update = 2; // mark the connected node as "dead"
      eventnode->back = tr->root;
      return;
    }
  }
}

void removenode(node *eventnode, tlist *tyme){
  tlist *tymelist;
  node *target;
  //update tymelist - remove eventnode from tymelist
  tyme->update = 1; // mark tymelist as "dead"
  
  target = findunique(eventnode);
  if (isrecomb(target)){

    // mark connected node as "dead"
    target->next->back->update = 2;
    findunique(target->next->back)->update = 2;
    target->next->next->back->update = 2;
    findunique(target->next->next->back)->update = 2;

    target->next->next->back->back = NULL;
    target->next->back->back = NULL;
	
    return;
  }
  else{ // coalescent

    if (target == curtree->root->back){ // is it root? - go back to find deepest live node.
      for (tymelist = tyme->prev; tymelist != NULL; tymelist = tymelist->prev){
	if (tymelist->update != 1){
	  curtree->root->back = tymelist->eventnode;
	  tymelist->eventnode->back = curtree->root;
	  return;
	}
      }
    }
    
    if (target->back == curtree->root){
      return;
    }
    
    if (target->next->update == 2 && target->next->next->update ==  2){  // both parent(daughter?) branches were removed - is it possible?
      target->back->update = 2;

      return;
    }

    else{  
      if (target->next->update == 2){
	target->back->back = target->next->next->back;
	target->next->next->back->back = target->back;
	return;
      }

      else{
	target->back->back = target->next->back;
	target->next->back->back = target->back;
		
	return;
      }
    }
  }		    
}
 
 
void updatenode(option_struct *op, data_fmt *data, brlist **br, node *eventnode, linlist **lineages, tlist *tyme, seglist *curseg, long *numlines){ // no new event case - update existing recombinant or coalescent node.
  int update = 0;
  long recsite;
  brlist *newbr;


  eventnode = findunique(eventnode);
  if (isrecomb(eventnode)){
    update = eventnode->update; 
    if (eventnode->back != NULL){
      if (eventnode->back->branch != NULL){
	br_remove(br,eventnode->back->branch);
	eventnode->back->branch = NULL;
      }
    }
  }
  else{
    if (eventnode->next->back != NULL){
      if (eventnode->next->back->branch != NULL){
	br_remove(br,eventnode->next->back->branch);
	eventnode->next->back->branch = NULL;
      }
    }
    if (eventnode->next->next->back != NULL){
      if (eventnode->next->next->back->branch != NULL){
	br_remove(br,eventnode->next->next->back->branch);
	eventnode->next->next->back->branch = NULL;
      }
    }
    if (eventnode->next->update > eventnode->next->next->update)
      update = eventnode->next->update;
    else 
      update = eventnode->next->next->update;
  }

  if (update != 2){

    freesegment(tyme->succ->segments);
    tyme->succ->segments = makenewseglist(curseg,eventnode); // update segment list and "coal" stuff.    
    if(iscoal(eventnode)){
      if (eventnode->back != curtree->root){
	if (!sameranges(eventnode->back->oldbackranges,eventnode->coal)) eventnode->update = 1;
      }
      else {
	if (eventnode->coal[0] != 0) update = 1;
      }
    }
    else{
      if (!sameranges(eventnode->next->back->oldbackranges,eventnode->next->coal)) eventnode->update = 1;
      if (!sameranges(eventnode->next->next->back->oldbackranges,eventnode->next->next->coal)) eventnode->update = 1;
    }
    update = eventnode->update;
  }
  else{
    eventnode->update = 2;
  }
  
  if (iscoal(eventnode) && eventnode->coal[0] == 0 && update != 2){
      eventnode->back->update = 2;
      eventnode->back = curtree->root;
      return;
  }


  if (update == 0) return;
 
  if (update == 1){ //just update eventnode
    if (iscoal(eventnode)){
      if (eventnode->back == curtree->root && eventnode->coal[0] != 0){ // FC -> active
	addlinlistH(op,data,lineages,eventnode,0);
	(*numlines)++;
	return;
      }
      else{
	if (eventnode->coal[1] > eventnode->back->oldbackranges[eventnode->back->oldbackranges[0] * 2] || eventnode->coal[eventnode->coal[0] * 2] < eventnode->back->oldbackranges[1]){
	  addlinlistH(op,data,lineages,eventnode,0);
	  (*numlines)++;
	  eventnode->back->update = 2;
	  return;
	}
	
	eventnode->back->update = 1;
	newbr = initbrlist(eventnode, eventnode->tyme,eventnode->back->tyme,countsites(op,data));
	updatebrs(newbr, eventnode->back->oldbackranges, eventnode->coal);
	//updatesegranges(op,data,newbr,NULL,eventnode->back->oldbackranges,eventnode->coal); 
	newbr->weight = count_weightbr(op,data,newbr);
	if (newbr->weight > 0){
	  hookup_brlist(br,newbr);
	  eventnode->branch = newbr;
	}
	else freebrlist(newbr);
      }
    }
    else{ // recombination node

      if (eventnode->next->recstart > eventnode->next->next->recstart){
	recsite = eventnode->next->recstart - 1;
      }
      else{
	recsite = eventnode->next->recend;
      }
      if (eventnode->back->coal[1] <= recsite && eventnode->back->coal[eventnode->back->coal[0] * 2] > recsite){ 
	eventnode->next->back->update = 1;
	eventnode->next->next->back->update = 1;

	if (eventnode->next->back->oldbackranges[1] > eventnode->next->coal[eventnode->next->coal[0] * 2] || eventnode->next->back->oldbackranges[eventnode->next->back->oldbackranges[0] * 2] < eventnode->next->coal[1]){ // totally new branch
	  addlinlistH(op,data,lineages,eventnode->next,0);
	  (*numlines)++;
	  eventnode->next->back->update = 2;
	}
	else{
	  newbr = initbrlist(eventnode->next, eventnode->next->tyme, eventnode->next->back->tyme,countsites(op,data));
	  updatebrs(newbr, eventnode->next->back->oldbackranges, eventnode->next->coal);
	  //updatesegranges(op,data,newbr,NULL,eventnode->next->back->oldbackranges,eventnode->next->coal);
	  //updatesegranges(op,data,newbr,NULL,eventnode->next->oldcoal,eventnode->next->coal); 
	  
	  newbr->weight = count_weightbr(op,data,newbr);
	  if (newbr->weight > 0){
	    hookup_brlist(br,newbr);
	    eventnode->next->branch = newbr;
	  }
	  else freebrlist(newbr);
	}

	if (eventnode->next->next->back->oldbackranges[1] > eventnode->next->next->coal[eventnode->next->next->coal[0] * 2] || eventnode->next->next->back->oldbackranges[eventnode->next->next->back->oldbackranges[0] * 2] < eventnode->next->next->coal[1]){
	  addlinlistH(op,data,lineages,eventnode->next->next,0);
	  (*numlines)++;
	  eventnode->next->next->back->update = 2;
	}	
	else{
	  newbr = initbrlist(eventnode->next->next, eventnode->next->next->tyme, eventnode->next->next->back->tyme,countsites(op,data));
	  updatebrs(newbr, eventnode->next->next->back->oldbackranges, eventnode->next->next->coal);
	  //updatesegranges(op,data,newbr,NULL,eventnode->next->next->back->oldbackranges,eventnode->next->next->coal);  
	  newbr->weight = count_weightbr(op,data,newbr);
	  if (newbr->weight > 0){
	    hookup_brlist(br,newbr);
	    eventnode->next->next->branch = newbr;
	  }
	  else freebrlist(newbr);
	}
	return;
      }
      else{ // remove recombinant node
	eventnode->update = 2;
	tyme->succ->update = 1;

	if (eventnode->back->coal[1] > recsite){ // kill "left" branch
	  if (eventnode->next->next->recend == recsite){ // which one is left?
	    eventnode->next->next->back->update = 2;

	    if (eventnode->next->back->oldbackranges[1] > eventnode->next->coal[eventnode->next->coal[0] * 2] || eventnode->next->back->oldbackranges[eventnode->next->back->oldbackranges[0] * 2] < eventnode->next->coal[1]){
	      eventnode->next->back->update = 2;
	      (*numlines)++;
	      addlinlistH(op,data,lineages,eventnode->back,0);
	    }
	    else{	    
	      if (!sameranges(eventnode->back->coal, findunique(eventnode->next->back)->coal)){
		eventnode->next->back->update = 1;
	      }
	      eventnode->next->back->back = eventnode->back;
	      eventnode->back->back = eventnode->next->back;
	      
	      newbr = initbrlist(eventnode->back, eventnode->back->tyme, eventnode->next->back->tyme,countsites(op,data));
	      updatebrs(newbr, eventnode->next->back->oldbackranges, eventnode->next->coal);
	      //updatesegranges(op,data,newbr,NULL,eventnode->next->back->oldbackranges,eventnode->back->coal);
	      //updatesegranges(op,data,newbr,NULL,eventnode->next->oldcoal,eventnode->back->coal); 
	      
	      newbr->weight = count_weightbr(op,data,newbr);
	      if (newbr->weight > 0){
		hookup_brlist(br,newbr);
		eventnode->back->branch = newbr;
	      }
	      else freebrlist(newbr);
	    }
	  }
	  else{
	    eventnode->next->back->update = 2;	 
	    
	    if (eventnode->next->next->back->oldbackranges[1] > eventnode->next->next->coal[eventnode->next->next->coal[0] * 2] || eventnode->next->next->back->oldbackranges[eventnode->next->next->back->oldbackranges[0] * 2] < eventnode->next->next->coal[1]){
	      eventnode->next->next->back->update = 2;
	      (*numlines)++;
	      addlinlistH(op,data,lineages,eventnode->back,0);
	    }
	    else{
	      if (!sameranges(eventnode->back->coal, findunique(eventnode->next->next->back)->coal)){
		eventnode->next->next->back->update = 1;
	      }
	      
	      eventnode->next->next->back->back = eventnode->back;
	      eventnode->back->back = eventnode->next->next->back;
	      
	      newbr = initbrlist(eventnode->back, eventnode->back->tyme, eventnode->next->next->back->tyme,countsites(op,data));
	      updatebrs(newbr,eventnode->next->next->back->oldbackranges,eventnode->back->coal);
	      //updatesegranges(op,data,newbr,NULL,eventnode->next->next->back->oldbackranges,eventnode->back->coal);	
	      //updatesegranges(op,data,newbr,NULL,eventnode->next->next->oldcoal,eventnode->coal); 
	      newbr->weight = count_weightbr(op,data,newbr);
	      if (newbr->weight > 0){
		hookup_brlist(br,newbr);
		eventnode->back->branch = newbr;
	      }
	      else freebrlist(newbr);
	    }
	  }
	  return;
	}
	else{ // kill "right" side

	  if (eventnode->next->recend == recsite){
	    eventnode->next->next->back->update = 2;

	    if (eventnode->next->back->oldbackranges[1] > eventnode->next->coal[eventnode->next->coal[0] * 2] || eventnode->next->back->oldbackranges[eventnode->next->back->oldbackranges[0] * 2] < eventnode->next->coal[1]){
	      eventnode->next->back->update = 2;
	      (*numlines)++;
	      addlinlistH(op,data,lineages,eventnode->back,0);
	    }
	    else{
	      if (!sameranges(eventnode->back->coal, findunique(eventnode->next->back)->coal)){
		eventnode->next->back->update = 1;
	      }
	      
	      eventnode->next->back->back = eventnode->back;
	      eventnode->back->back = eventnode->next->back;

	      newbr = initbrlist(eventnode->back, eventnode->back->tyme, eventnode->next->back->tyme,countsites(op,data));
	      updatebrs(newbr, eventnode->next->back->oldbackranges, eventnode->back->coal);
	      // updatesegranges(op,data,newbr,NULL,eventnode->next->back->oldbackranges,eventnode->back->coal);
	      //updatesegranges(op,data,newbr,NULL,eventnode->next->oldcoal,eventnode->back->coal); 
	      
	      newbr->weight = count_weightbr(op,data,newbr);
	      if (newbr->weight > 0){
		hookup_brlist(br,newbr);
		eventnode->back->branch = newbr;
	      }
	      else freebrlist(newbr);
	    }
	  }
	  else{
	    eventnode->next->back->update = 2;

	    if (eventnode->next->next->back->oldbackranges[1] > eventnode->next->next->coal[eventnode->next->next->coal[0] * 2] || eventnode->next->next->back->oldbackranges[eventnode->next->next->back->oldbackranges[0] * 2] < eventnode->next->next->coal[1]){
	      eventnode->next->next->back->update = 2;
	      (*numlines)++;
	      addlinlistH(op,data,lineages,eventnode->back,0);
	    }
	    else{
	      if (!sameranges(eventnode->back->coal, findunique(eventnode->next->next->back)->coal)){
		eventnode->next->next->back->update = 1;
	      }
	      
	      eventnode->next->next->back->back = eventnode->back;
	      eventnode->back->back = eventnode->next->next->back;
	      
	      newbr = initbrlist(eventnode->back, eventnode->next->next->tyme, eventnode->next->next->back->tyme,countsites(op,data));
	      updatebrs(newbr, eventnode->next->next->back->oldbackranges, eventnode->back->coal);
	      //updatesegranges(op,data,newbr,NULL,eventnode->next->next->back->oldbackranges,eventnode->back->coal);
	      //updatesegranges(op,data,newbr,NULL,eventnode->next->next->oldcoal,eventnode->back->coal); 
	      
	      newbr->weight = count_weightbr(op,data,newbr);
	      if (newbr->weight > 0){
		hookup_brlist(br,newbr);
		eventnode->back->branch = newbr;
	      }
	      else freebrlist(newbr);
	    }
	  }
	  return;
	}    
      }
    }
  }
  if (update == 2){ // remove eventnode and mark connected nodes as "futile"
    if (iscoal(eventnode)){
      if (isdeadnode(eventnode)){
	eventnode->back->update = 2;	  
	removenode(eventnode,tyme->succ);
	return;
      }
      if (eventnode->next->update != 2){ // update brlist or update linlist (fc case).  
	if (!isactive(eventnode->next->back,*lineages)){
	  if (eventnode->coal[0] == 0){ //was fc node?
	    addlinlistH(op,data,lineages,eventnode->next->back,0);
	    (*numlines)++;
	    removenode(eventnode,tyme->succ);
	    return;
	  }

	  if (eventnode->coal[1] > eventnode->next->back->coal[eventnode->next->back->coal[0] * 2] || eventnode->coal[eventnode->coal[0] * 2] <  eventnode->next->back->coal[1]){// totally new seg info?
	    eventnode->back->update = 2;
	    addlinlistH(op,data,lineages,eventnode->next->back,0);
	    (*numlines)++;
	    removenode(eventnode,tyme->succ);
	    return;
	  }
	    
	  eventnode->next->back->back = eventnode->back;
	  eventnode->back->back = eventnode->next->back;
	    
	  if (!sameranges(eventnode->next->back->coal, eventnode->coal)){ 
	    eventnode->back->update = 1;
	    newbr = initbrlist(eventnode->next->back, eventnode->next->back->tyme,eventnode->back->tyme,countsites(op,data));
	    updatebrs(newbr, eventnode->back->oldbackranges, eventnode->next->back->coal);
	    //updatesegranges(op,data,newbr,NULL,eventnode->coal,eventnode->next->back->coal);
	      
	    newbr->weight = count_weightbr(op,data,newbr);
	    if (newbr->weight > 0){
	      hookup_brlist(br,newbr);
	      eventnode->next->back->branch = newbr;
	    }
	    else freebrlist(newbr);
	  }
	}
	else{
	  eventnode->back->update = 2;
	}
      }
      if (eventnode->next->next->update != 2){ // update brlist.
	if (!isactive(eventnode->next->next->back,*lineages)){
	    
	  if (eventnode->coal[0] == 0){ //fc case?
	    addlinlistH(op,data,lineages,eventnode->next->next->back,0);
	    (*numlines)++;
	    removenode(eventnode,tyme->succ);
	    return;
	  }

	  eventnode->next->next->back->back = eventnode->back;
	  eventnode->back->back = eventnode->next->next->back;
	  

	  if (eventnode->coal[1] > eventnode->next->next->back->coal[eventnode->next->next->back->coal[0] * 2] || eventnode->coal[eventnode->coal[0] * 2] <  eventnode->next->next->back->coal[1]){// totally new seg info?
	    eventnode->back->update = 2;
	    addlinlistH(op,data,lineages,eventnode->next->next->back,0);
	    (*numlines)++;
	    removenode(eventnode,tyme->succ);
	    return;
	  }


	  if (!sameranges(eventnode->next->next->back->coal, eventnode->coal)){
	    eventnode->back->update = 1;
	    newbr = initbrlist(eventnode->next->next->back, eventnode->next->next->back->tyme,eventnode->back->tyme,countsites(op,data));
	    updatebrs(newbr,eventnode->back->oldbackranges, eventnode->next->next->back->coal);
	    //updatesegranges(op,data,newbr,NULL,eventnode->coal,eventnode->next->next->back->coal);
	      
	    newbr->weight = count_weightbr(op,data,newbr);
	    if (newbr->weight > 0){
	      hookup_brlist(br,newbr);
	      eventnode->next->next->back->branch = newbr;
	    }
	    else freebrlist(newbr);
	  }	    
	}	      
	else{
	  eventnode->back->update = 2;
	}
      }
    }
    
    if (isrecomb(eventnode)){
      eventnode->next->back->update = 2;
      eventnode->next->next->back->update = 2;
    }
	
    removenode(eventnode,tyme->succ);
    return;
  }
  return;
  
}
 
void modifynode(node *eventnode, tlist *tyme, seglist *cursegment){ // new coalescent or recombination event 
  int i;
  node *next_node = NULL;
  eventnode = findunique(eventnode);

  //printf("modify node %ld\n",eventnode->id);
  if (tyme->succ->segments) freesegment(tyme->succ->segments);
  if (isrecomb(eventnode)){ // recombination event
    tyme->succ->segments = makenewseglist(cursegment,eventnode);

    if (eventnode->next->back != NULL || eventnode->next->next->back != NULL){ // brlist recomb
      if (eventnode->next->back != NULL) eventnode->next->back->update = 1;
      else eventnode->next->next->back->update = 1;
    }
    return;
  }
  else{ // coalescent -> mark connected node as "update"
    tyme->succ->segments = makenewseglist(cursegment,eventnode);
    if (eventnode->coal[0] == 0){ //final coalesce of all
      if (eventnode->back != NULL)  eventnode->back->update = 2;
      eventnode->back = curtree->root;
      return;
    }
    else{
      if (eventnode->back != NULL) next_node = findunique(eventnode->back);
      else return;
      if (next_node != curtree->root && eventnode->back != NULL && next_node->update != 2){
	next_node->update = 0;
	if (eventnode->coal[0] == next_node->coal[0]){
	  for (i = 1; i <= eventnode->coal[0]; i++){
	    if (eventnode->coal[2 * i - 1] != next_node->coal[2 * i - 1] || eventnode->coal[2 * i] != next_node->coal[2 * i]){
	      next_node->update = 1;
	      break;
	    }
	  }
	  return;
	}
	next_node->update = 1;
	return;
      }
      return;
    }
  }
}
    
void marknode(node *p){
  p->back->update = 2;
  return;
}

void removetymelist(option_struct *op,tree *target){
  tlist *tyme = NULL,*rmtyme = NULL;
   
  target->numcoals = target->numrecombs = 0;
  tyme = target->tymelist;
  while (1){
    if (tyme->update == 1){
      if (tyme->prev != NULL) {
	tyme->prev->age = tyme->age;
	tyme->prev->succ = tyme->succ;
      }
      else curtree->tymelist = tyme->succ;
      if (tyme->succ != NULL) tyme->succ->prev = tyme->prev;
      else tyme->age = 10000;
      rmtyme = tyme;
      tyme = tyme->succ;
      freeeventnode(op, rmtyme->eventnode);
      freesegment(rmtyme->segments);
      free(rmtyme->branchlist);				     
      free(rmtyme);
    }
    else{
      if (iscoal(tyme->eventnode)) {
	target->numcoals ++;
      }
      if (isrecomb(tyme->eventnode)) {
	target->numrecombs ++;
      }
      tyme->eventnode->update = tyme->eventnode->next->update = tyme->eventnode->next->next->update = 0;
      if (tyme->succ != NULL) tyme = tyme->succ;
      else {
	curtree->root->back = tyme->eventnode;
	return;
      }
    }
    if (tyme == NULL) return;
  }
}

boolean isnewrange(long *old, long *new){
  if (old[1] > new[1]) return (FALSE);
  if (old[old[0] * 2] < new[new[0] * 2]) return (FALSE);
  return (TRUE);
}

void removebranches(tlist *tyme){
  curtree->root->back = tyme->eventnode;
  while (tyme == NULL){
    tyme = tyme->succ;
    tyme->update = 1;
    free(tyme->eventnode->next->next);
    free(tyme->eventnode->next);
    free(tyme->eventnode);
  }
  return;
}
  
    
void updatetymelist(tlist *tymelist){ // update branchlist part of tymelist
  int i,j;
  tlist *prevtime;

  if (tymelist->eventnode->next->update == 2 || tymelist->eventnode->next->next->update == 2) tymelist->update = 1;
  
  if (tymelist->prev == NULL || tymelist->update == 1) return;

  prevtime = tymelist->prev;
  while (prevtime->update == 1){
    prevtime = prevtime->prev;
  }
  //prevtime = activet;
  
  tymelist->numbranch = prevtime->numbranch;
  i = 0;

  if (iscoal(tymelist->eventnode) && tymelist->eventnode->coal[0] == 0) i = 1;
  
  free(tymelist->branchlist);
  if (iscoal(tymelist->eventnode)){
    tymelist->numbranch --;
    if (tymelist->numbranch == 0) tymelist->numbranch++; // final root

    tymelist->branchlist = (node **)calloc(tymelist->numbranch, sizeof(node *));
    if (i != 1){
      tymelist->branchlist[0] = tymelist->eventnode;
      for (j = 1 ; i < prevtime->numbranch; i++){
	if (prevtime->branchlist[i] != tymelist->eventnode->next->back && prevtime->branchlist[i] != tymelist->eventnode->next->next->back){
	  tymelist->branchlist[j] = prevtime->branchlist[i];
	  j ++;
	}
      }
    }
    else{
      tymelist->numbranch --;
      if (tymelist->numbranch == 0){
	tymelist->numbranch++; // final root  
	tymelist->branchlist[0] = tymelist->eventnode;
      }
      else{	
	for (j = 0, i = 0 ; i < prevtime->numbranch; i++){
	  if (prevtime->branchlist[i] != tymelist->eventnode->next->back && prevtime->branchlist[i] != tymelist->eventnode->next->next->back){
	    tymelist->branchlist[j] = prevtime->branchlist[i];
	    j ++;
	  }
	}
	tymelist->branchlist[tymelist->numbranch] = tymelist->eventnode;
      }
    }
    return;
  }	
  else{
    tymelist->numbranch ++;
    //tymelist->branchlist = (node **)realloc(tymelist->branchlist,tymelist->numbranch * sizeof(node *));
    tymelist->branchlist = (node **)calloc(tymelist->numbranch,sizeof(node *));
    tymelist->branchlist[0] = findunique(tymelist->eventnode)->next;
    tymelist->branchlist[1] = findunique(tymelist->eventnode)->next->next;
    for (j = 2; i < prevtime->numbranch; i++){
      if (prevtime->branchlist[i] != findunique(tymelist->eventnode)->back){
	tymelist->branchlist[j] = prevtime->branchlist[i];
	j ++;
      }
    }
    return;
  }
}

void updaterootnode(tree *tr, node *rootnode, tlist *tyme, linlist **lineages, seglist *curseg){

  int update = 0;
  
  if (rootnode->next->update >= rootnode->next->next->update)
    update = rootnode->next->update;
  else
    update = rootnode->next->next->update;

  rootnode->update = update;

  if (update != 2){
    freesegment(tyme->segments);
    tyme->segments = makenewseglist(curseg,rootnode);
    return;
  }
  
  if (update == 1){
    tr->numcoals++;
    updateranges(tr, rootnode); //update ranges and mark connected nodes as "upate"
    return;
  }
}


void contrib_temp(node *eventnode, long **newcoal)
{
  node *q, *r;
  long i, numremove;
 
  if (iscoal(eventnode)) {
    q = eventnode->next->back;
    r = eventnode->next->next->back;
    init_ranges_alloc(newcoal,q->coal[0]);
    memcpy((*newcoal),q->coal,(2*q->coal[0]+2)*sizeof(long));
    for(i = 1; r->coal[i] != FLAGLONG; i+=2) addrange(newcoal,r->coal[i],r->coal[i+1]);
    return;
  }

  if(eventnode != findunique(eventnode)) {
    q = findunique(eventnode)->back;
    init_ranges_alloc(newcoal,q->coal[0]);
    memcpy((*newcoal),q->coal,(2*q->coal[0]+2)*sizeof(long));
    /* first remove the leading excess ranges */
    for(i = 1, numremove = 0; (*newcoal)[i] != FLAGLONG; i+=2) {
      if(eventnode->recstart > (*newcoal)[i+1]) {numremove+=2; continue;}
      if(eventnode->recstart > (*newcoal)[i]) (*newcoal)[i] = eventnode->recstart;
      break;
    }
    memmove(&(*newcoal)[1],&(*newcoal)[1+numremove],
	    (2*(*newcoal)[0]-numremove+1)*sizeof(long));
    (*newcoal)[0] -= numremove/2;
    /* then removing the trailing excess ranges */
    for(i = 2*(*newcoal)[0]-1, numremove = 0; i+1 != 0; i-=2) {
      if(eventnode->recend < (*newcoal)[i]) {numremove++; continue;}
      if(eventnode->recend < (*newcoal)[i+1]) (*newcoal)[i+1] = eventnode->recend;
      break;
    }
    (*newcoal)[i+2] = FLAGLONG;
    (*newcoal)[0] -= numremove;
  }
  else {
    q = eventnode->back;
    init_ranges_alloc(newcoal,q->coal[0]);

    memcpy((*newcoal),q->coal,(2*q->coal[0]+2)*sizeof(long));
  }
  return;
}
    
      

void copyoldranges(node *source, long **oldranges){
  init_ranges_alloc(oldranges,source->coal[0]+1);
  memcpy((*oldranges),source->ranges,(source->coal[0]*2+2)*sizeof(long));
  return;
}

struct seglist *makenewseglist(seglist *origin, node *eventnode){// update seglist and coal stuff.

  node *node1, *node2;
  seglist *newseglist;
  seglist *tempsegment, *cursegment;
  int index1,index2,key,key1,key2;
  long recsite , i;
  long *tempcoal, *oldcoal;
    
  eventnode = findunique(eventnode);
  newseglist = copyseglist(origin);

  if (iscoal(eventnode)){
    node1 = eventnode->next->back;
    node2 = eventnode->next->next->back;

    tempcoal = (long *)calloc((node1->coal[0] + node2->coal[0]) * 100 + 2,sizeof(long));
      
    index1 = index2 = 1;
    tempcoal[0] = 0;
    for (cursegment = newseglist; cursegment != NULL; cursegment = cursegment->next){
      key1 = key2 = 0;
      for (i = index1; i <= node1->coal[0]; i++){
	if (node1->coal[2 * i - 1] <= cursegment->start && node1->coal[2 * i] >= cursegment->start){
	  key1 = 1;
	  if (node1->coal[2 * i] == cursegment->end) index1 ++;
	  break;
	}
      }
      for (i = index2; i <= node2->coal[0]; i++){
	if (node2->coal[2 * i - 1] <= cursegment->start && node2->coal[2 * i] >= cursegment->start){
	  key2 = 1;
	  if (node2->coal[2 * i] == cursegment->end) index2++;
	  break;
	}
      }
	
      if (key1 == 1 || key2 == 1){
	  
	if (key1 == 1 && key2 == 1){
	  if (cursegment->numsam == 2){//final coalescent
	    cursegment->numsam = 0;
	    eventnode->update = 1;
	  }
	  else{
	    cursegment->numsam --;
	  }
	}
	  
	if (cursegment->numsam != 0){
	  if (tempcoal[0] == 0){
	    tempcoal[0] = 1;
	    tempcoal[1] = cursegment->start;
	    tempcoal[2] = cursegment->end;
	  }
	  else{
	    if (tempcoal[2 * tempcoal[0]] + 1 == cursegment->start){//continuous segment
	      tempcoal[2 * tempcoal[0]] = cursegment->end;
	    }
	    else{
	      tempcoal[0]++;
	      tempcoal[2 * tempcoal[0] - 1] = cursegment->start;
	      tempcoal[2 * tempcoal[0]] = cursegment->end;
	    }
	  }
	}
      }
    }    
    oldcoal = eventnode->coal;
    eventnode->coal = NULL;
    
    init_ranges_alloc(&eventnode->coal,tempcoal[0]);
    memcpy(eventnode->coal,tempcoal,(tempcoal[0]*2+2)*sizeof(long));
    free(tempcoal);
    //if (eventnode->coal != NULL) free(eventnode->coal);
    //eventnode->coal = NULL;
    //init_ranges_alloc(&eventnode->coal,tempcoal[0]);
    //coalpart = (long *)realloc(tempcoal,(tempcoal[0] * 2 + 2)*sizeof(long));
    //free(eventnode->coal);
    //eventnode->coal = tempcoal;
    eventnode->coal[eventnode->coal[0] * 2 + 1] =  FLAGLONG;
    if (!sameranges(eventnode->coal,oldcoal)) eventnode->update = 1;
    free(oldcoal);
    return newseglist;
  }
  else{ // recombination node 
    if (eventnode->update == 1){	
      free(eventnode->coal);
      free(eventnode->next->coal);
      free(eventnode->next->next->coal);
      eventnode->coal = eventnode->next->coal = eventnode->next->next->coal = NULL;
      contrib_temp(eventnode, &eventnode->coal);
      contrib_temp(eventnode->next,&eventnode->next->coal);
      contrib_temp(eventnode->next->next,&eventnode->next->next->coal);
	
    }	


    if (eventnode->next->recstart == 0) recsite = eventnode->next->recend;
    else recsite = eventnode->next->next->recend;
      
    key = 0;
    i = 1;
    while (1){
      if (eventnode->coal[2 * i - 1] <= recsite && recsite < eventnode->coal[2 * i]){
	key = 1;
	break;
      }
      if (i >= eventnode->coal[0]) break;
      if (recsite < eventnode->coal[2 * i]) break;
      i ++;
    }
      
    if (key == 1){
      for (cursegment = newseglist; cursegment != NULL ; cursegment = cursegment->next){
	if (cursegment->start <= recsite && cursegment->end > recsite){
	  tempsegment = (struct seglist *)malloc(sizeof (struct seglist));
	  tempsegment->start = recsite + 1;
	  tempsegment->end = cursegment->end;
	  tempsegment->numsam = cursegment->numsam;
	  tempsegment->next = cursegment->next;
	  tempsegment->prev = cursegment;
	    
	  cursegment->end = recsite;
	  cursegment->next = tempsegment;
	  return newseglist;
	}
      }
    }

    return newseglist;
  }
}

boolean isdeadnode(node *target){
  if (target->next->back == NULL && target->next->next->back == NULL)
    return (TRUE);
  else{
    if (target->next->back == NULL || target->next->next->back == NULL){
      if (target->next->back == NULL){
	if (target->next->next->back->back != target->next->next)
	  return (TRUE);
	else
	  return (FALSE);
      }
      else{
	if (target->next->back->back != target->next)
	  return (TRUE);
	else
	  return (FALSE);
      }
    }
    if (target->next->back->back != target->next && target->next->next->back->back != target->next->next)
      return (TRUE);
    else
      return (FALSE);
  }
}
	
void freesegment(seglist *target){
  if (target->next != NULL) freesegment(target->next);
  
  free(target);
  target = NULL;
}


void freeeventnode(option_struct *op,node *source){
  int i;
  node *target;
  
  target = source;
  for (i = 0; i < 3; i++){
    free_x(op,target);
    if(op->map)free_z(op,target);
    if (target->nayme) free(target->nayme);
    
    if (target->ranges != NULL) free(target->ranges);
    if (target->oldbackranges != NULL) free(target->oldbackranges);
    if (target->coal != NULL) free(target->coal);
    target->coal = target->ranges = target->oldbackranges = NULL;

    if (target->branch != NULL) free(target->branch);
    target = target->next;
  }
  free(target->next->next);
  free(target->next);
  free(target);
  return;
}


//update "coal"
void updatecoalpart(node *eventnode){
  if (iscoal(eventnode)){
    eventnode = findunique(eventnode);
    init_ranges_alloc(&eventnode->next->coal,eventnode->next->back->coal[0]);
    memcpy((eventnode->next->coal),eventnode->next->back->coal,(eventnode->next->back->coal[0]*2+2)*sizeof(long));

    init_ranges_alloc(&eventnode->next->next->coal,eventnode->next->next->back->coal[0]);
    memcpy((eventnode->next->next->coal),eventnode->next->next->back->coal,(eventnode->next->next->back->coal[0]*2+2)*sizeof(long));
    
    contrib_temp_coal(eventnode, &eventnode->coal);
    return;
  }
  if (isrecomb(eventnode)){
    eventnode = findunique(eventnode);
    init_ranges_alloc(&eventnode->coal,eventnode->back->coal[0]);
    memcpy((eventnode->coal),eventnode->back->coal,(eventnode->back->coal[0]*2+2)*sizeof(long));
    contrib_temp_coal(eventnode->next, &eventnode->next->coal);
    contrib_temp_coal(eventnode->next->next, &eventnode->next->next->coal);
    return;
  }
}


void contrib_temp_coal(node *eventnode, long **newranges)
{
  node *q, *r;
  long i, numremove;
 
  //eventnode = findunique(eventnode);

  if (iscoal(eventnode)) {
    q = eventnode->next->back;
    r = eventnode->next->next->back;
    init_ranges_alloc(newranges,q->coal[0]);
    memcpy((*newranges),q->coal,(2*q->coal[0]+2)*sizeof(long));
    for(i = 1; r->coal[i] != FLAGLONG; i+=2)
      addrange(newranges,r->coal[i],r->coal[i+1]);
    return;
  }

  if(isrecomb(eventnode)) {
    q = findunique(eventnode)->back;
    init_ranges_alloc(newranges,q->coal[0]);
    memcpy((*newranges),q->coal,(2*q->coal[0]+2)*sizeof(long));
    /* first remove the leading excess ranges */
    for(i = 1, numremove = 0; (*newranges)[i] != FLAGLONG; i+=2) {
      if(eventnode->recstart > (*newranges)[i+1]) {numremove+=2; continue;}
      if(eventnode->recstart > (*newranges)[i]) (*newranges)[i] = eventnode->recstart;
      break;
    }
    memmove(&(*newranges)[1],&(*newranges)[1+numremove], (2*(*newranges)[0]-numremove+1)*sizeof(long));
    (*newranges)[0] -= numremove/2;
    /* then removing the trailing excess ranges */
    for(i = 2*(*newranges)[0]-1, numremove = 0; i+1 != 0; i-=2) {
      if(eventnode->recend < (*newranges)[i]) {numremove++; continue;}
      if(eventnode->recend < (*newranges)[i+1]) (*newranges)[i+1] = eventnode->recend;
      break;
    }
    (*newranges)[i+2] = FLAGLONG;
    (*newranges)[0] -= numremove;
  }
  else {
    q = eventnode->back;
    init_ranges_alloc(newranges,q->coal[0]);

    memcpy((*newranges),q->coal,(2*q->coal[0]+2)*sizeof(long));
  }
  return;
}

boolean growlineagesH(option_struct *op, data_fmt *data, tree *growtree,tree *oldtreeup, tlist *tstart, brlist **brlines, linlist **alines,long *numalines, double offsetstart, boolean **nodenumber,long *numnodes, long *siteptr) {
  tlist *tymeslice, *activetyme;
  double newtyme, offset;
  double curtyme = 0;
  boolean succeeded, rootdropped;
  char eventtype;
  bseg brs;
  seglist *cursegment;


  activetyme = tymeslice = tstart;
  offset = offsetstart;
  rootdropped = FALSE;
  succeeded = TRUE;
  growtree->numrecombs = growtree->numcoals = 0;
  cursegment = tymeslice->segments;
  curtyme = tymeslice->eventnode->tyme;

  while(tymeslice != NULL) {
    if (tymeslice->succ == NULL){ // root node
      if ((*numalines) == 0 && !(*brlines)) break;
      succeeded = rootdropH(op,data,tymeslice,alines,numalines,siteptr,nodenumber,numnodes,brlines,cursegment,activetyme,tymeslice->eventnode->tyme);
      rootdropped = TRUE;
      break;
    }
    if ((*numalines) != 0 || (*brlines)){
      newtyme = tymeslice->eventnode->tyme + eventprobH_growth(op,data,(*alines),activetyme,brlines,&eventtype,&brs);

      if (newtyme <= tymeslice->age) { // new node between two tymeslices
	updateweights(alines,brlines,newtyme - curtyme);
	succeeded = makeeventH(op,data,growtree,tymeslice,brlines,alines,numalines,newtyme,nodenumber,numnodes,siteptr,eventtype,rootdropped,&brs, cursegment);
	curtyme = newtyme;
      }
      else{ // no new node - update existing node at tymeslice->succ
	updateweights(alines,brlines,tymeslice->age - curtyme);

	if (tymeslice->succ->succ == NULL && tymeslice->succ->eventnode->update != 0) updaterootnode(growtree,tymeslice->succ->eventnode,tymeslice->succ,alines,cursegment);
	else updatenode(op,data,brlines,tymeslice->succ->eventnode,alines,tymeslice,cursegment,numalines);

	if (findunique(tymeslice->succ->eventnode)->update != 2) {
	  updatetymelist(tymeslice->succ);//update branchlist part
	  tymeslice->succ->update = 0;

	}
	else  tymeslice->succ->update = 1;
      }

      if (tymeslice->succ->update == 0){
	cursegment = tymeslice->succ->segments;
	activetyme = tymeslice->succ;
	curtyme = tymeslice->age;
      }
    }
    else{ // no active branches - update remaining tree.
      if (tymeslice->succ->succ == NULL && tymeslice->succ->eventnode->update != 0) updaterootnode(growtree,tymeslice->succ->eventnode,tymeslice->succ,alines,cursegment);
      else updatenode(op,data,brlines,tymeslice->succ->eventnode,alines,tymeslice,cursegment,numalines);
      if (findunique(tymeslice->succ->eventnode)->update != 2) {
	updatetymelist(tymeslice->succ);//update branchlist part
	tymeslice->succ->update = 0;
      }
      else  tymeslice->succ->update = 1;
      
      if (tymeslice->succ->update != 1){
	cursegment = tymeslice->succ->segments;
	activetyme = tymeslice->succ;
	curtyme = tymeslice->age;
      }
    }
    if (!succeeded) break;
    if (tymeslice->succ != NULL) tymeslice = tymeslice->succ;
  }
  removetymelist(op, growtree);
  checkbranch(growtree);
  if (!succeeded) {
    numdropped++;
    finishbadtree(op,data,growtree,alines,numalines,nodenumber,numnodes,rootdropped,brlines);
  }
  //succeeded = checktree(growtree);
  return(succeeded);
}

void checkbranch(tree *target){
  tlist *current;
  long i;
  node *targetnode;
  
  for (i = 0; i < target->tymelist->numbranch; i++){
      target->tymelist->branchlist[i]->branch = NULL;
  } 
  current = target->tymelist->succ;
  while (1){
    if (current->eventnode->oldbackranges != NULL)
      free(current->eventnode->oldbackranges);    
    if (current->eventnode->next->oldbackranges != NULL)
      free(current->eventnode->next->oldbackranges);
    if (current->eventnode->next->next->oldbackranges)
      free(current->eventnode->next->next->oldbackranges);
    current->eventnode->oldbackranges = current->eventnode->next->oldbackranges = current->eventnode->next->next->oldbackranges = NULL;
    if (iscoal(current->eventnode)){
      targetnode = findunique(current->eventnode);
      targetnode->next->back->length = targetnode->next->length = targetnode->tyme - targetnode->next->back->tyme;
      targetnode->next->next->back->length = targetnode->next->next->length = targetnode->tyme - targetnode->next->next->back->tyme;
      if (targetnode->coal[0] == 0) targetnode->back = curtree->root;
    }
    else{
      targetnode = findunique(current->eventnode);
      targetnode->back->length = targetnode->length = targetnode->tyme - targetnode->back->tyme;
    }      

    if (current->eventnode->branch != NULL) current->eventnode->branch  = NULL;
    if (current->eventnode->next->branch != NULL) current->eventnode->next->branch  = NULL;
    if (current->eventnode->next->next->branch != NULL) current->eventnode->next->next->branch  = NULL;
    if (current->succ != NULL) current = current->succ;
    else{
      targetnode->length = 10000 - targetnode->tyme;
      curtree->root->back = targetnode;
      curtree->root->back->length = targetnode->length;
      curtree->root->back = targetnode;
      return;
    }    
  }
}



boolean checktree(tree *target){
  tlist *cur;
  node *curnode;
  
  cur = target->tymelist->succ;

  while (1){
    curnode = findunique(cur->eventnode);
    if (isrecomb(curnode)){
      if (curnode->back->back != curnode){
	printf("error 1\n");
	return FALSE;
      }
      if (curnode->next->back->back != curnode->next){	
	printf("error 1.2\n");
	return FALSE;
      }
      if (curnode->next->next->back->back != curnode->next->next){
	printf("error 1.3\n");
	return FALSE;
      }      
      if (curnode->next->coal[1] < curnode->next->next->coal[1]){
	if (curnode->next->coal[1] != curnode->back->coal[1]) {
	  printf("erroe 1.4\n");
	  return FALSE;
	}
	if (curnode->next->next->coal[curnode->next->next->coal[0] * 2] != curnode->back->coal[curnode->back->coal[0] * 2]){
	  printf("error 1.5\n");
	  return FALSE;
	}
      }
      else{
	if (curnode->next->next->coal[1] != curnode->back->coal[1]) {
	  printf("erroe 1.6\n");
	  return FALSE;
	}
	if (curnode->next->coal[curnode->next->coal[0] * 2] != curnode->back->coal[curnode->back->coal[0] * 2]){
	  printf("error 1.7\n");
	  return FALSE;
	}
      }
	



    }
    if (iscoal(curnode)){
      if (curnode->back == NULL){
	printf("error 5\n");
	return FALSE;
      }
      if (curnode->next->back->back != curnode->next || curnode->next->next->back->back != curnode->next->next) {
	printf("error 4\n");   
	return FALSE;   
      }
    }
    if (cur->succ == NULL){
      if (target->root->back != curnode){
	printf("error 5\n");
	return FALSE;
      }
      break;
    }
    cur = cur->succ;
  }
  return TRUE;
}
  
  

void freepath(path *target){
  if (target == NULL) return;
  if (target->next != NULL){
    freepath(target->next);
  }
  free(target);
  return;
}
    
void freerecinfo(recnumb *target){
  if (target == NULL) return;
  if (target->next != NULL){
    freerecinfo(target->next);
  }
  free(target);
  return;
}
  
treeinfo *treesum(tree *targettree, recnumb *recinfo, path *pathinfo){
  recnumb *curseg;
  path *curpath;
  treeinfo *result;
  tlist *temptlist;
  long i;

  result = (treeinfo *)calloc(1,sizeof(treeinfo));

  for(curseg = recinfo; curseg != NULL; curseg = curseg->next){
    for (i = curseg->beg; i < curseg->end; i++){
      if (findstate(pathinfo,i) == 0) result->weight_cold += curseg->swt;
      else result->weight_hot += curseg->swt;
    }
    if (i == curseg->end){
      if (findstate(pathinfo,i) == 0){
	result->weight_cold += curseg->last_weight;
	result->numrec_cold += curseg->numrec;
      }
      else{
	result->weight_hot += curseg->last_weight;
	result->numrec_hot += curseg->numrec;
      }
    }
  }

  for(temptlist = targettree->tymelist; temptlist->succ != NULL; temptlist = temptlist->succ){
    result->branch += temptlist->numbranch * (temptlist->numbranch - 1.0) * (temptlist->age - temptlist->eventnode->tyme);
  }

  result->numcoal = targettree->numcoals;

  for(curpath = pathinfo; curpath != NULL; curpath = curpath->next){
    if (curpath->next != NULL){
      if (curpath->state == 0){
	result->CC += curpath->end - curpath->beg;
	result->CH ++;
      }
      else{
	result->HH += curpath->end - curpath->beg;
	result->HC ++;
      }
    }
    else{   
      if (curpath->state == 0){
	result->CC += curpath->end - curpath->beg;
      }
      else{
	result->HH += curpath->end - curpath->beg;
      }
    }
  }

  return result;
}

//test
void *EM(treeinfo **tree_sum, double th0, double recH0, double recC0, double CC0, double HH0, long numtree){
  //parmset *parmn;
  double th,recH,recC,CC,HH;
  double sumweight_hot,sumweight_cold;
  double sumrec_hot,sumrec_cold;
  double sumbranch, sumcoal;
  double sumCC,sumCH,sumHH,sumHC;
  double weight;
  long i, numC, numH;
  double com_L, com_L_old;

  th = th0;
  recH = recH0;
  recC = recC0;
  CC = CC0;
  HH = HH0;
  numC = numH =0;

  while (1){
    sumweight_hot = sumweight_cold = sumrec_hot = sumrec_cold = sumbranch = sumcoal = sumCC = sumCH = sumHH = sumHC = 0;
    com_L = com_L_old = 0;
    //test

    for (i = 100; i<numtree;i++){
      numC += tree_sum[i]->numrec_cold;
    
      numH += tree_sum[i]->numrec_hot;
    }


    for (i = 100; i<numtree;i++){
      weight = gettreeweight(tree_sum[i], th, recH, recC, HH, CC, th0, recH0, recC0, HH0, CC0);
      com_L_old += weight;
      sumbranch += weight * tree_sum[i]->branch;
      sumcoal += weight * tree_sum[i]->numcoal;
    }
    th = sumbranch/sumcoal;

    //bigweight = find_bigweight(treesum,numtree,recC0,recC,recH0,recH);

    for (i = 100; i<numtree;i++){
      weight = gettreeweight(tree_sum[i], th, recH, recC, HH, CC, th0, recH0, recC0, HH0, CC0);
      sumweight_hot += weight * tree_sum[i]->weight_hot;
      sumrec_hot += weight * tree_sum[i]->numrec_hot;
    }
    recH = sumrec_hot/sumweight_hot;



    for (i = 100; i<numtree;i++){
      weight = gettreeweight(tree_sum[i], th, recH, recC, HH, CC, th0, recH0, recC0, HH0, CC0);
      sumweight_cold += weight * tree_sum[i]->weight_cold;
      sumrec_cold += weight * tree_sum[i]->numrec_cold;
    }
    recC = sumrec_cold/sumweight_cold;
    for (i = 100; i<numtree;i++){
      weight = gettreeweight(tree_sum[i], th, recH, recC, HH, CC, th0, recH0, recC0, HH0, CC0);
      sumHH += weight * tree_sum[i]->HH;
      sumHC += weight * tree_sum[i]->HC;
    }
    HH = sumHH / (sumHH + sumHC);

    for (i = 100; i<numtree;i++){
      weight = gettreeweight(tree_sum[i], th, recH, recC, HH, CC, th0, recH0, recC0, HH0, CC0);
      sumCC += weight * tree_sum[i]->CC;
      sumCH += weight * tree_sum[i]->CH;
    }
    CC = sumCC / (sumCC + sumCH);


    for (i = 100; i<numtree;i++){
      
      com_L += gettreeweight(tree_sum[i], th, recH, recC, HH, CC, th0, recH0, recC0, HH0, CC0);
    }
    //printf("L = %lf\n",com_L);
    printf("theta = %lf, recH = %lf\n",th, recH);
    if (com_L < com_L_old) printf ("error in EM\n");
    //if (com_L == com_L_old) return;
    
  }
}
double gettreeweight(treeinfo *tree_rec, double th, double recH, double recC, double HH, double CC, double th0, double recH0, double recC0, double HH0, double CC0){
  double prob_r, prob_theta, prob_lamda;

  prob_r = pow((recH / recH0), tree_rec->numrec_hot) * exp( - tree_rec->weight_hot * (recH - recH0))
    * pow((recC / recC0), tree_rec->numrec_cold) * exp( - tree_rec->weight_cold * (recC - recC0)) ;
  
  prob_theta = pow((th0 / th), tree_rec->numcoal) * exp( - tree_rec->branch * (1.0/th - 1.0/th0));
  
  prob_lamda = pow((HH/HH0), tree_rec->HH) * pow((1 - HH) / (1 - HH0), tree_rec->HC) * pow((CC/CC0), tree_rec->CC) * pow((1 - CC) / (1 - CC0), tree_rec->CH);

  return prob_r * prob_theta * prob_lamda;
}
  
	
    
    




double find_bigweight(treeinfo **treessum, long numtree,double recC0, double recC, double recH0, double recH){
  double biggest = 0;
  long i;

  for (i = 0; i<numtree; i++){
    if (biggest < - treessum[i]->weight_cold * (recC0 - recC) - treessum[i]->weight_hot * (recH0 - recH))
      biggest =- treessum[i]->weight_cold * (recC0 - recC) - treessum[i]->weight_hot * (recH0 - recH);
    
  }
  return biggest;
}

void freeR(Rs *target){
  if (target == NULL) return;
  if (target->next != NULL) freeR(target->next);
  free(target);
  return;
}


void checkpath(path *target){
  long i;
  if (target->next != NULL){
    if (target->state == 1){
      HtH += target->end - target->beg;
      for (i = target->beg; i <= target->end; i++) sumhot[i]++;
      HtC ++;
    }
    else{
      CtC+= target->end - target->beg;
      CtH ++;
    }
    checkpath(target->next);
  }
  else{
    if (target->state == 1){
      HtH += target->end - target->beg;
      for (i = target->beg; i <= target->end; i++) sumhot[i]++;
    }
    else{
      CtC+= target->end - target->beg;
    }
    return;
  }
  return;
}

void checkrecnum(recnumb *target){
  if (target->next != NULL){
    checkrecnum(target->next);
  }
/*   if (target->end == 0) part1 = part1 + target->numrec; */
/*   if (target->end == 1) part2 = part2 + target->numrec; */
/*   if (target->end == 2) part3 = part3 + target->numrec; */
/*   if (target->end > 600 && target->end <= 800) part4 = part4 + target->numrec; */
/*   if (target->end > 800) part5  = part5 + target->numrec; */
  if (target->end < 101) part1 = part1 + target->numrec;
  if (target->end > 100 && target->end <= 200) part2 = part2 + target->numrec;
  if (target->end > 200 && target->end <= 300) part3 = part3 + target->numrec;
  if (target->end > 300 && target->end <= 400) part4 = part4 + target->numrec;
  if (target->end > 400) part5  = part5 + target->numrec;
  
  return;
    
}

void constructoldranges(tree *target){
  tlist *cur;
  long *oldranges1, *oldranges2;
  int tip;
  cur = target->tymelist->succ;

  for (tip = 1; tip <= old_hap; tip++){
    curtree->nodep[tip]->coal[0] = 1;
  }


  while(1){
    if (iscoal(cur->eventnode)){      
      if (findunique(cur->eventnode)->next->oldbackranges != NULL) free(findunique(cur->eventnode)->next->oldbackranges);    
      if (findunique(cur->eventnode)->next->next->oldbackranges != NULL) free(findunique(cur->eventnode)->next->next->oldbackranges);
      
      oldranges1 = calloc(cur->eventnode->next->back->coal[0] * 2 + 2, sizeof(long));
      oldranges2 = calloc(cur->eventnode->next->next->back->coal[0] * 2 + 2, sizeof(long));
      
      memcpy(oldranges1, cur->eventnode->next->back->coal,(cur->eventnode->next->back->coal[0] * 2 + 2)*sizeof(long));
      memcpy(oldranges2, cur->eventnode->next->next->back->coal,(cur->eventnode->next->next->back->coal[0] * 2 + 2)*sizeof(long));
      
      cur->eventnode->next->oldbackranges = oldranges1;
      cur->eventnode->next->next->oldbackranges = oldranges2;
    }
    else{
      if (findunique(cur->eventnode)->oldbackranges != NULL) free(findunique(cur->eventnode)->oldbackranges);

      oldranges1 = calloc(findunique(cur->eventnode)->back->coal[0] * 2 + 2, sizeof(long));
      memcpy(oldranges1, findunique(cur->eventnode)->back->coal,(findunique(cur->eventnode)->back->coal[0] * 2 + 2)*sizeof(long));
      
      findunique(cur->eventnode)->oldbackranges = oldranges1;
    }
    if (cur->succ == NULL) return; 
    cur = cur->succ;
  }
  return;
}
      
void updatebrs(brlist *target, long *old, long *new){
  target->ofsite = old[1];
  target->olsite = old[old[0] * 2];
  
  target->nfsite = new[1];
  target->nlsite = new[new[0] * 2];
  
  return;
}

     
      
node *pickbranch_temp(option_struct *op, data_fmt *data, tree *source, tlist **tymeslice, double *offset) {
  long i;
  double depth, sum = 0, target;
  tlist *t;
  node *p;

  depth = curtree->root->back->tyme;

  sum = depth * getdata_numtips(op,data);
  
  for (t = source->tymelist->succ; t != NULL; t = t->succ) {
    if (isrecomb(t->eventnode))  sum += (depth - t->eventnode->tyme) * 2; 
    else if (findunique(t->eventnode)->coal[0] != 0) sum += depth - t->eventnode->tyme;
  }

  target = sum * randum();

  for (i = 0; i < getdata_numtips(op,data); i++){
    target = target -  depth;
    if (target < 0) {
      (*tymeslice) = source->tymelist;
      return source->nodep[i+1];
    }
  }


  for (t = source->tymelist->succ; t != NULL; t = t->succ) {
    if (isrecomb(t->eventnode)) {
      target = target - ((depth - t->eventnode->tyme) * 2); 
      if (target < 0){
	(*tymeslice) = t;
	p = findunique(t->eventnode);
	return ((randum() < 0.5) ? p->next : p->next->next);
      }
    }
    else{
      if (findunique(t->eventnode)->coal[0] != 0){
	target = target - depth - t->eventnode->tyme;       
	if (target < 0){
	  (*tymeslice) = t;
	  return findunique(t->eventnode);
	}
      }
    }
  }
  printf("problem at pick branch\n");
  return NULL;
}



//void setupinitialtree(){
  
boolean isin(long *coal, long site){
  long y;
  
  if (coal == NULL) return FALSE;
    
  for (y = 1; y <= coal[0]; y++){
    if (coal[2 * y - 1] <= site && coal[2 * y] >= site){
      return TRUE;
    }
    if (coal[2 * y - 1] > site) return FALSE;
  }
  return FALSE;
}

long howmanynodes(tree *targettree, long site){
  long y,numnodes;
  tlist *t;
  node *targetnode;

  numnodes = targettree->tymelist->numbranch;

  for (t = targettree->tymelist->succ; t->succ != NULL; t = t->succ){
    targetnode = findunique(t->eventnode);
    if (iscoal(targetnode) && targetnode->coal[0] != 0){
      if (targetnode->coal[1] <= site && targetnode->coal[targetnode->coal[0] * 2] >= site){
	for (y = 1; y <= targetnode->coal[0]; y++){
	  if (targetnode->coal[2 * y - 1] <= site && targetnode->coal[2 * y] >= site){
	    numnodes ++;
	    break;
	  }
	}
      }
    }
    if (isrecomb(targetnode)){
      targetnode = targetnode->back;    
      if (targetnode->coal[1] <= site && targetnode->coal[targetnode->coal[0] * 2] >= site){
	for (y = 1; y <= targetnode->coal[0]; y++){
	  if (targetnode->coal[2 * y - 1] <= site && targetnode->coal[2 * y] >= site){
	    numnodes ++;
	    break;
	  }
	}
      }
    }
  }
  return numnodes;
}

node *selectwithD(tree *target){
  double x = 0;
  long i, site = -1, targetn,y;
  tlist *t;
  node *targetnode;
        
  x = -target->likelihood * randum();
  
  for (i = 0; i < seq_length; i++){
    x = x + target->dlikelihood[i];
    if (x < 0){
      site = i;
      break;
    }
  }
  
  targetsite = site;

  targetn = howmanynodes(target,site) * randum();

  if (targetn < target->tymelist->numbranch){
    return target->tymelist->branchlist[targetn];
  }  
  else targetn = targetn - target->tymelist->numbranch + 1;
  
  for (t = target->tymelist->succ; t->succ != NULL; t = t->succ){
    targetnode = findunique(t->eventnode);
    if (iscoal(targetnode) && targetnode->coal[0] != 0){
      if (targetnode->coal[1] <= targetsite && targetnode->coal[targetnode->coal[0] * 2] >= targetsite){
	for (y = 1; y <= targetnode->coal[0]; y++){
	  if (targetnode->coal[2 * y - 1] <= targetsite && targetnode->coal[2 * y] >= targetsite){
	    targetn --;
	    if (targetn <= 0) return targetnode;
	    break;
	  }
	}
      }
    }
    if (isrecomb(targetnode)){
      targetnode = targetnode->back;    
      if (targetnode->coal[1] <= targetsite && targetnode->coal[targetnode->coal[0] * 2] >= targetsite){
	for (y = 1; y <= targetnode->coal[0]; y++){
	  if (targetnode->coal[2 * y - 1] <= targetsite && targetnode->coal[2 * y] >= targetsite){
	    targetn --;
	    if (targetn <= 0){
	      targetnode = targetnode->back;
	      if (targetnode->next->coal[0] > targetsite) return targetnode->next->next;
	      else return targetnode->next;
	    break;
	    }
	  }
	}
      }
    }  
  }
  printf("no nodes selected\n");
  return NULL;
}





double Hastingsratio(tree *oldtree, tree *newtree, long site){
  double ratio;
  
  
  ratio = (newtree->dlikelihood[site]  - newtree->likelihood) 
    - (oldtree->dlikelihood[site] - oldtree->likelihood);

  //printf("%ld,%ld\n",howmanynodes(newtree,site),howmanynodes(oldtree,site));
  ratio = ratio + log(((1.0 / (double)howmanynodes(newtree,site)) / (1.0/ (double)howmanynodes(oldtree,site))));

  return ratio;
}
  
tlist *findtymeslice(node *target){
  tlist *temp;
  
  if (istip(target)) return curtree->tymelist;

  for (temp = curtree->tymelist; temp->succ != NULL; temp = temp->succ){
    if (temp->eventnode->tyme == target->tyme) return temp;
  }
  printf("found no tymelist\n");
  return NULL;
}



void checkstruct(){
  long i, hots=0;
  for (i=0;i<seq_length;i++){
    if (recarray[i] > 1 && i>7000 && i<8000) hots++;
  }
  
  printf(" hot on hot %ld\n", hots);
}
      
void checktarget(node *target){
  long i,j, hot = 0, hotonhot = 0, hotoncold = 0;
  
  for(i = 1; i <= target->coal[0]; i++){
    for (j = target->coal[2 * i - 1]; j < target->coal[2 * i]; j++){
      if (j > 400 && j <600){
	hot++;
	if (recarray[j] > 0.05) hotonhot++;
      }
      else{
	if (recarray[j] > 0.05) hotoncold++;
      }
      
    }
  }
  printf("  time %lf, target: hot %ld, hots on hot %ld hots on cold %ld \n", target->tyme, hot, hotonhot,hotoncold);
  
}

double treeprob(tree *target){
  tlist *t;
  double tk=0;
  double tyme, prob;

  for (t = target->tymelist; t->succ != NULL; t=t->succ) {
    tyme = t->age - t->eventnode->tyme;
    tk += t->numbranch * (t->numbranch - 1.0) * tyme;
  }

  prob = -(tk/theta0) + (log(2.0/theta0))*target->numcoals;

  prob += hmmlikelihood(target);
   
  return prob;
}

double treesprobS(tree *oldtree, tree *newtree, double *rates){
  double ratio = 0;
  long i;

  for (i = 0; i<seq_length; i++){
    ratio += log(recarray[i]) * (oldtree->numrec_array[i] - newtree->numrec_array[i]) 
      - rates[i] * (oldtree->weight_array[i] - newtree->weight_array[i]);
  }

  return ratio;
}
  
void updateweights(linlist **actives, brlist **brs, double tymeinterval){
  linlist *active;
  brlist *branch;
  for (active = *actives; active != NULL; active = active->succ){
    addlinweights(active, tymeinterval);
  }
  for (branch = *brs; branch != NULL; branch = branch->succ){
    addbrweights(branch, tymeinterval);
  }

  return;
}

void addlinweights(linlist *target, double tyme){
  long i;

  for (i = target->start; i < target->end; i++) newgweight[i] += tyme;
 
  return;
}

void addbrweights(brlist *target, double tyme){
  long i;

  for(i = target->nfsite; i < target->ofsite; i++) newgweight[i] += tyme;

  for(i = target->olsite; i < target->nlsite; i++) newgweight[i] += tyme;

  return;
}

node *picktipbranch(option_struct *op, data_fmt *data, tree *target, tlist **tymeslice, double *offset){
  long i, numtips;
  double *dlike, totaldlike = 0;
  tree *temptree;

  (*tymeslice) = target->tymelist;
  *offset = 0.0;
  numtips = getdata_numtips(op,data);
  temptree = copytree(op,data,target);
  dlike = (double *)calloc(numtips,sizeof(double));

  for (i = 0; i < numtips; i++){
    notip = i;
    calcratio(op, data, temptree, temptree, temptree->nodep[1]->coal);
    tiplikelist[notip] = temptree->likelihood;
    
    totaldlike += tiplikelist[notip];
  }
  
  freetree(op,data,temptree);

  notip = -99;

  totaldlike = totaldlike * randum();
  for (i = 0; i < numtips; i++){
    totaldlike = totaldlike - tiplikelist[i];
    if (totaldlike > 0){
      targettip = i;
      return target->nodep[i + 1];
    }
  }
  printf("problem at pick branch\n");
  return NULL;
}

boolean tipstestratio(option_struct *op, data_fmt *data, tree *oldtree, tree *newtree){
  long i, numtips;
  double totaloldlike = 0, totalnewlike = 0, oldtip, newtip, *newlist;
  double test = 0, x;

  numtips = getdata_numtips(op,data);
  newlist = (double *)calloc(numtips,sizeof(double));

  for (i = 0; i < numtips; i++){
    totaloldlike += tiplikelist[i];
    
    notip = i;
    calcratio(op, data, newtree, newtree, newtree->nodep[1]->coal);
    newlist[i] = newtree->likelihood;
    totalnewlike += newlist[i];
  }
  notip = -99;

  oldtip = tiplikelist[targettip];
  newtip = newlist[targettip];

  test += newtree->likelihood - oldtree->likelihood;
  
  test += log(oldtip/newtip) + log(totalnewlike/totaloldlike);

  if (test >= 0){
    for (i = 0; i<numtips; i++) tiplikelist[i] = newlist[i];
    free(newlist);
    return TRUE;
  }

  x = log(randum());

  if (x <= test){
    for (i = 0; i<numtips; i++) tiplikelist[i] = newlist[i];
    free(newlist);
    return TRUE;
  }
  else{
    free(newlist);
    return FALSE;
  }
}
