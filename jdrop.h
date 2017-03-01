#define DROP_INCLUDE
#ifndef RECOMBINE_INCLUDE
#include "recombine.h"
#endif

typedef struct brlist {
  struct brlist *prev, *succ;
  node *branchtop;
  double starttyme, endtyme;
  int *segranges; /* int to save space */
  double numnewactives;
  long nfsite, nlsite, ofsite, olsite;
  double weight;
  boolean updated;
} brlist;

typedef struct linlist {
  struct linlist *prev, *succ;
  node *branchtop;
  double activesites;
  long start;
  long end;
  double weight;
} linlist;

typedef struct bseg {
   brlist *target;
   long cross;
} bseg;


long countsites(option_struct *op, data_fmt *data);

brlist *initbrlist(node *branchtop, double starttyme, double endtyme,
   long numsites);
void freebrlist(brlist *br);
void delete_brlist(brlist **bigbr);
void hookup_brlist(brlist **bigbr, brlist *newbr);
boolean isintyme(brlist *br, double target);
void set_sitesbr(option_struct *op, data_fmt *data, brlist *br);
double count_activebr(option_struct *op, data_fmt *data, brlist *br);
brlist *getbrlistbystart(brlist *start, node *target, double starttyme);
brlist *getbrlistbyend(brlist *start, node *target, double endtyme);
brlist *getbrlistbytyme(brlist *start, double searchtyme);
void updatesegranges(option_struct *op, data_fmt *data, brlist *br,
   int *oldsegranges, long *oldranges, long *newranges);
void segranges_to_ranges(int *segranges, long **ranges, long numsites);
void addbrlist(option_struct *op, data_fmt *data, brlist **bigbr,
   node *newtop, long *oldranges, long *newranges, int *segmembers,
   double starttyme, double endtyme, boolean building);

void renamebrlist(brlist **bigbr, brlist *start, node *oldtop, node *newtop);
void drop_renamebrlist(brlist **bigbr, node *oldtop, node *newtop);
void drop_brfix(option_struct *op, data_fmt *data, tree *newtree,
   tlist *start, brlist **oldbr);
void br_remove(brlist **bigbr, brlist *target);
void printbrlist(brlist *start);
void treelegalbrlist(option_struct *op, data_fmt *data, tree *tr,
   brlist **bigbr, brlist *target);
boolean samebranch(brlist *br1, brlist *br2);
long checktreebranch(brlist *br, node *p, double starttyme, double endtyme);
long sametreebranch(brlist *br, node *p, double starttyme, double endtyme);
boolean isoverlap(brlist *br1, brlist *br2);
boolean iscontainedin(brlist *br1, brlist *br2);
void mash_brlist(tree *tr, brlist *start, brlist *target);
void consolidate_brlist(option_struct *op, data_fmt *data, tree *tr,
   brlist **bigbr);
boolean samesites(int *segranges, long *ranges, long numsites);

void init_numnewactive(option_struct *op, data_fmt *data, tree *tr,
   brlist **bigbr, node *branchtop);
brlist *pickbr_recomb(dnadata *dna, brlist **bigbr, double droptyme, long *cross);
boolean brlistrecomb(option_struct *op, data_fmt *data, tlist *t,
   linlist **lineages, long *lines, brlist **bigbr, double tyme,
   long *siteptr, boolean **nodenumber, long *numnodes, bseg *brs, tree *tr, seglist *cursegment);

boolean is_futile(node *p);
node *find_nonfutile(node *p, node *cutnode, boolean upwards);
void tag_futile(node *p);
node *find_futilecoaldtr(node *cutnode, node *p);
void remove_futile(option_struct *op, data_fmt *data, tree *newtree,
   tree*oldtree, node *start, boolean *nodenumber, brlist **br);
void readdbranch(tree *newtree, tlist *u, node *addbranch);
void extbranch(tree *newtree, tlist *u, node *addbranch);
void rembranch(tree *newtree, tlist *u, node *subbranch);

void remove_coalnode(tree *tr, node *p);
void remove_recombnode(tree *tr, node *p, node *remove_this);

void rename_branch(tlist *start, node *source, node *target);
void remove_branch(tlist *start, node *target);
void subtymenode(tree *tr, node *p);
void fix_treenodep(tree *tr, long number);

node *findpaired_coalnode(option_struct *op, node *p);
void remove_pairedcoalnode(option_struct *op, data_fmt *data, tree *tr,
  node *p);
void remove_pairedrecombnode(option_struct *op, data_fmt *data, tree *tr,
  node *p, node *dead);
void remove_pairednodes(option_struct *op, data_fmt *data, tree *tr,
   node *p, node *remove_thisway);

boolean foundbranch(tlist *t, node *p);
boolean isactive(node *p, linlist *u);
long countinactive(option_struct *op, tlist *t, linlist *u);
boolean rootdrop(option_struct *op, data_fmt *data, tlist *t, double starttyme, linlist **lineages, 
   long *lines, long *siteptr, boolean **nodenumber, long *numnodes,
   brlist **br);
node *pickbranch (option_struct *op, data_fmt *data, tree *source, 
   tlist **tymeslice, double *offset);
long setnodenumber(boolean **nodenumber, long *numnodes);
void readdactive(tlist *stop, node *p);
node *pickactive(long numalin, linlist **lineages, node *prohibited);
node *pickinactive(option_struct *op, tlist *t, node *prohibited,
   linlist *alines);
void edit_alias(option_struct *op, data_fmt *data, long *sp, long cutpoint);
void traverse_rebuild_alias(option_struct *op, data_fmt *data, 
   tlist *tstart, long *sp);
void rebuild_alias(option_struct *op, data_fmt *data, long *sp);
boolean popbranch(tlist *t, long branch);
boolean subbranch(tlist *t, long oldbranch, node *newbranch);
void fixbranch(tlist *t);
void poptymelist (node *p);
boolean find_rootloop(tlist *start, tlist **found);
boolean isdead(option_struct *op, node *p);
void free_from_tymelist(tlist *t);
void update_rangecoal(option_struct *op, data_fmt *data, tree *tr);
void remove_rootloop(option_struct *op, data_fmt *data, long *sp);
void remove_deadrecombs(option_struct *op, data_fmt *data, tree *tr);
void remove_excess_tree(option_struct *op, data_fmt *data, tree *tr,
   long *siteptr);
brlist *pickbrtarget(brlist *start, double numactive);
long pickbrcross(option_struct *op, data_fmt *data, brlist *source);
double eventprobbr(option_struct *op, data_fmt *data, brlist **bigbr,
   tlist *tint, double pc, double pr, double *pb, bseg *brs);
double eventprob(option_struct *op, data_fmt *data, linlist *alines,
   tlist *t, brlist **bigbr, char *event, bseg *brs);
long counttymelist(tlist *t);
void renumber_nodes(option_struct *op, data_fmt *data, tree *tr);
void finishbadtree(option_struct *op, data_fmt *data, tree *tr, linlist **lineages, long *numlins, 
   boolean **nodenumber, long *numnodes, boolean rootdropped, 
   brlist **brlines);
boolean makeevent(option_struct *op, data_fmt *data, tree *growtree,
   tlist *tymeslice, brlist **brlines, linlist **alines, long *numalines,
   double newtyme, boolean **nodenumber, long *numnodes, long *siteptr,
   char eventtype, boolean rootdropped, bseg *brs);
boolean growlineages(option_struct *op, data_fmt *data, tree *growtree,
   tree *oldtree, tlist *tstart, brlist **brlines, linlist **alines,
   long *numalines, double offsetstart, boolean **nodenumber,
   long *numnodes, long *siteptr);
boolean makedrop(option_struct *op, data_fmt *data);
void coalesce(option_struct *op, data_fmt *data, tlist *t, linlist **lineages, long *lines, double tyme,
   boolean **nodenumber, long *numnodes, boolean rootdropping,
   brlist **br);
long choosepsite(option_struct *op, data_fmt *data, node *p,
   double targetlink);
long choosepsitefc(option_struct *op, data_fmt *data, node *p,
   double targetlink);
boolean recomb(option_struct *op, data_fmt *data, tlist *t, linlist **lineages, long *lines, double tyme,
  long *siteptr, boolean **nodenumber, long *numnodes);
void newlin(linlist **t);
void freelin(linlist *t);
void freelinlist(linlist *t);
void addlinlist(linlist **t, node *target, double activelinks);
void sublinlist(linlist **t, node *target);
void printlinlist(linlist *u);
void findlin(linlist **u, long which);
boolean foundlin(linlist *u, node *p);
double getnumlinks(option_struct *op, data_fmt *data, long start, long end);
double count_active(option_struct *op, data_fmt *data, node *p);
double count_coal_active(option_struct *op, data_fmt *data, node *p);
double count_rec_active(option_struct *op, data_fmt *data, node *p);
double count_activefc(option_struct *op, data_fmt *data, node *p);
double count_coal_activefc(option_struct *op, data_fmt *data, node *p);
double count_rec_activefc(option_struct *op, data_fmt *data, node *p);
void makecoal(option_struct *op, data_fmt *data, tree *tr,
   tlist *tfc, long **coal, node *p, long *subtrees);
boolean is_fc(tlist *treesec, long site);
void contrib(option_struct *op, data_fmt *data, node *p, long **newranges);
void traverse_flagbelow(tree *tr, node *p);
void flag_below(node *p);
void traverse_unflag(tree *tr, node *p);
void unflag(node *p);
long countrec(tlist *start);
node *findrec(tlist *start, long target);
boolean twiddle(option_struct *op, data_fmt *data);
void datachangesite(option_struct *op, data_fmt *data, node *p1, node *p2,
   long site);
boolean isprohibited(long value, long numbad, long *badstuff);
boolean isinvariant(option_struct *op, data_fmt *data, long site,
   creature *cr);
long choose_flipsite(option_struct *op, data_fmt *data, creature *cr,
   long numprohibited, long *prohibited);
boolean fliphap(option_struct *op, data_fmt *data, tree *tr);
boolean flipdrop(option_struct *op, data_fmt *data, tree *tr, 
   long numdrop);
boolean isintymelist(option_struct *op, data_fmt *data, tlist *tylist,
   node *p);
void flipdrop_prunebr(option_struct *op, data_fmt *data, tree *tr,
   brlist **br);
void addfractrecomb(option_struct *op, data_fmt *data, long chain,
   treerec ***treessum);
void add_one_recomb(option_struct *op, data_fmt *data, tree *tr,
   long *siteptr, long chain);

void fix_coal_ranges(option_struct *op, data_fmt *data, node *p);
void fix_rec_ranges(option_struct *op, data_fmt *data, node *p);
void fix_ranges(option_struct *op, data_fmt *data, tree *tr, node *p);

void fix_coal(option_struct *op, data_fmt *data, tree *tr, tlist *tfc,
   node *p);




boolean twiddleH(option_struct *op, data_fmt *data);
boolean growlineagesH(option_struct *op, data_fmt *data, tree *growtree,tree *oldtree, tlist *tstart, brlist **brlines, linlist **alines,long *numalines, double offsetstart, boolean **nodenumber,long *numnodes, long *siteptr);
double eventprobH(option_struct *op, data_fmt *data, linlist *alines, tlist *t, brlist **bigbr, char *event, bseg *brs);
boolean makeeventH(option_struct *op, data_fmt *data, tree *growtree,tlist *tymeslice, brlist **brlines, linlist **alines, long *numalines,double newtyme, boolean **nodenumber, long *numnodes, long *siteptr, char eventtype, boolean rootdropped, bseg *brs, seglist *curseg);
boolean recombH(option_struct *op, data_fmt *data, tlist *t,linlist **lineages, long *lines, double tyme, long *siteptr, boolean **nodenumber, long *numnodes, tree *tr, seglist *cursegs);
long getsites(long start, long end);
void addlinlistH(option_struct *op, data_fmt *data,linlist **t, node *target, double activelinks);
void coalesceH(option_struct *op, data_fmt *data, tlist *t,linlist **lineages, long *lines, double tyme, boolean **nodenumber,long *numnodes, boolean rootdropping, brlist **br, seglist *cursegs);
long findrecstart(option_struct *op, data_fmt *data, node *p);
long findrecend(option_struct *op, data_fmt *data, node *p);
boolean makedropH(option_struct *op, data_fmt *data);
double eventprobbrH(option_struct *op, data_fmt *data, brlist **bigbr,
		    tlist *tint, double pc, double pr, double *pb, bseg *brs);
brlist *pickbrtargetH(brlist *start, long weight);
long pickbrcrossH(option_struct *op, data_fmt *data, brlist *source);
void init_numnewactiveH(option_struct *op, data_fmt *data, tree *tr,
			brlist **bigbr, node *branchtop);
double count_weightbr(option_struct *op, data_fmt *data, brlist *br);
float count_weighted_tlist(option_struct *op, data_fmt *data, tlist *t);
float count_activeH(option_struct *op, data_fmt *data, node *p);
double getweight(long start, long end);
float getweightedlinks(option_struct *op, data_fmt *data, long start, long end);
float count_coal_activeH(option_struct *op, data_fmt *data, node *p);
double count_activefcH(option_struct *op, data_fmt *data, node *p);
float count_rec_activeH(option_struct *op, data_fmt *data, node *p);
double count_activefcH(option_struct *op, data_fmt *data, node *p);
double eventprobbrH(option_struct *op, data_fmt *data, brlist **bigbr,
		    tlist *tint, double pc, double pr, double *pb, bseg *brs);
brlist *pickbrtargetH(brlist *start, long weight);
long pickbrcrossH(option_struct *op, data_fmt *data, brlist *source);
void init_numnewactiveH(option_struct *op, data_fmt *data, tree *tr,
			brlist **bigbr, node *branchtop);
double count_weightbr(option_struct *op, data_fmt *data, brlist *br);
long findcoalpos(option_struct *op, data_fmt *data, node *p,int key);
boolean rootdropH(option_struct *op, data_fmt *data, tlist *t, linlist **lineages, long *lines, long *siteptr, boolean **nodenumber, long *numnodes, brlist **br, seglist *curseg, tlist *activetyme, double basetyme);

void addbrlistH(option_struct *op, data_fmt *data, brlist **bigbr,
   node *newtop, long *oldranges, long *newranges, int *segmembers,
   double starttyme, double endtyme, boolean building);
void updaterange(node *target);
void modifynode(node *eventnode,tlist *tyme, seglist *cursegment);
void updatenode(option_struct *op, data_fmt *data, brlist **br, node *eventnode, linlist **lineages, tlist *tyme, seglist *curseg, long *numlines);

void removetymelist(option_struct *op, tree *target);
void updaterootnode(tree *tr, node *rootnode, tlist *tyme, linlist **lineages, seglist *curseg);
void updatetymelist(tlist *tymelist);
void contrib_temp(node *eventnode, long **newranges);
void copyoldranges(node *source, long **oldranges);
boolean isnewrange(long *old, long *new);
seglist *makenewseglist(seglist *origin, node *eventnode);// update segl
void marknode(node *p);
boolean isdeadnode(node *target);
void freesegment(seglist *target);
void freeeventnode(option_struct *op, node *source);
void updatecoalpart(node *eventnode);
void contrib_temp_coal(node *eventnode,long **newranges);
void checkbranch(tree *target);
boolean checktree(tree *target);

seglist *copyseglist(seglist *original);
void freepath(path *target);
void freerecinfo(recnumb *target);
treeinfo *treesum(tree *targettree, recnumb *recinfo, path *pathinfo);
double gettreeweight(treeinfo *tree_rec, double th, double recH, double recC, double HH, double CC, double th0, double recH0, double recC0, double HH0, double CC0);
void *EM(treeinfo **tree_sum, double th0, double recH0, double recC0, double CC0, double HH0, long numtree);

void removecreatedtymelist(option_struct *op, tree *target);

void removefutiletymelist(option_struct *op, tree *target);
void cleanfutiletlist(option_struct *op, tree *target);
void removeoldtlist(option_struct *op, tlist *target);
void copynodep(option_struct *op, data_fmt *data, tree *target);
void initoldstuff(option_struct *op, data_fmt *data, tree *target);
double find_bigweight(treeinfo **treessum, long numtree,double recC0, double recC, double recH0, double recH);
void freeR(Rs *target);
void checkpath(path *target);
void constructoldranges(tree *target);
void updatebrs(brlist *target, long *old, long *new);
      
node *pickbranch_temp(option_struct *op, data_fmt *data, tree *source, tlist **tymeslice, double *offset);
node *selectwithD(tree *target);
long howmanynodes(tree *targettree, long site);
boolean isin(long *coal, long site);
double Hastingsratio(tree *oldtree, tree *newtree, long site);
tlist *findtymeslice(node *target);

void add(node *targetnode, double *recs, double *weights);

double calcbranchprob(node *targetnode);
node *pickbranchwithS(option_struct *op, data_fmt *data, tree *source, tlist **tymeslice, double *offset);
double tran_prob(tree *targettree, node *targetnode);
double treeprob(tree *target);
void updateg(node *targetnode, node *mothernode, brlist *targetbr1, brlist *targetbr2, int eventtype);  
void updateweights(linlist **actives, brlist **brs, double tymeinterval);
void addlinweights(linlist *target, double tyme);
void addbrweights(brlist *target, double tyme);
node *picktipbranch(option_struct *op, data_fmt *data, tree *target, tlist **tymeslice, double *offset);
boolean tipstestratio(option_struct *op, data_fmt *data, tree *oldtree, tree *newtree);


