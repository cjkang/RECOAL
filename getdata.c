#ifndef GETDATA_INCLUDE
#include "getdata.h"
#endif

#ifdef DMEMDEBUG
#define MEMDEBUG
#include "memdebug.h"
#endif

unsigned int baseA=0, baseC=1, baseG=2, baseT=3;

extern FILE *simlog, *weightfile;
extern long locus, population;

extern void read_line(FILE *source, char *line);
extern void setupmsatdata(msatdata **ms, long numpop, long numloci,
   long numind);
extern void freemsatdata(msatdata *ms);
extern void calculate_steps(msatdata *ms);
extern boolean isinvariant(option_struct *op, data_fmt *data, long site,
  creature *cr);


/**********************************************************
 * setupdata() initializes the appropiate data structures */
void setupdata(data_fmt *data, char datatype, long numpop, long numloci,
  long numind, long *numsites, char *title)
{

switch(datatype) {
   case 'a':
      data->dnaptr = NULL;
      data->msptr = NULL;
      break;
   case 'b':
   case 'm':
      data->dnaptr = NULL;
      setupmsatdata(&(data->msptr),numpop,numloci,numind);
      strcpy(data->msptr->title,title);
      calculate_steps(data->msptr);
      break;
   case 'n':
   case 's':
      setupdnadata(&(data->dnaptr),numsites,numind,numloci,numpop);
      strcpy(data->dnaptr->title,title);
      data->msptr = NULL;
      break;
   default:
      fprintf(ERRFILE,"setupdata, found an unknown datatype %c\n",
         datatype);
      exit(-1);
}
} /* setupdata */


/***********************************************
 * freedata() frees the general data structure */
void freedata(data_fmt *data)
{

if (data->dnaptr) {freednadata(data->dnaptr); data->dnaptr = NULL;}
if (data->msptr) {freemsatdata(data->msptr); data->msptr = NULL;}
} /* freedata */


/********************************************************************
 * getdata_numpop() returns the number of populations stored in the *
 * general data structure or FLAGLONG if it fails.                  */
long getdata_numpop(option_struct *op, data_fmt *data)
{

switch(op->datatype)
{
   case 'a':
      break;
   case 'b':
   case 'm':
      return(data->msptr->numpop);
      break;
   case 'n':
   case 's':
      return(data->dnaptr->numpop);
      break;
   default:
      return(FLAGLONG);
      break;
}

return(FLAGLONG);

} /* getdata_numpop */


/**************************************************************
 * getdata_numloci() returns the number of loci stored in the *
 * general data structure or FLAGLONG if it fails.            */
long getdata_numloci(option_struct *op, data_fmt *data)
{

switch(op->datatype)
{
   case 'a':
      break;
   case 'b':
   case 'm':
      return(data->msptr->numloci);
      break;
   case 'n':
   case 's':
      return(data->dnaptr->numloci);
      break;
   default:
      return(FLAGLONG);
      break;
}

return(FLAGLONG);

} /* getdata_numloci */


/********************************************************************
 * getdata_nummarkers() returns the number of markers that the user *
 * provides in the input data.  FLAGLONG is returned on failure.    */
long getdata_nummarkers(option_struct *op, data_fmt *data)
{

switch(op->datatype)
{
   case 'a':
      break;
   case 'b':
   case 'm':
      return(data->msptr->numloci);
      break;
   case 'n':
   case 's':
      return(data->dnaptr->sites[locus]);
      break;
   default:
      return(FLAGLONG);
      break;
}

return(FLAGLONG);

} /* getdata_nummarkers */


/**********************************************************************
 * getdata_numseq() returns the number of individual sequences in the *
 * data stored in the general data structure or FLAGLONG if it fails. */
long getdata_numseq(option_struct *op, data_fmt *data)
{

switch(op->datatype)
{
   case 'a':
      break;
   case 'b':
   case 'm':
      return(data->msptr->numind[population]*NUMCHROM);
      break;
   case 'n':
   case 's':
      return(data->dnaptr->numseq[0]);
      break;
   default:
      return(FLAGLONG);
      break;
}

return(FLAGLONG);

} /* getdata_numseq */


/********************************************************************
 * getdata_numtips() returns the number of tips present in the tree *
 * stored in the general data structure or FLAGLONG if it fails.    */
long getdata_numtips(option_struct *op, data_fmt *data)
{

switch(op->datatype)
{
   case 'a':
      break;
   case 'b':
   case 'm':
      return(data->msptr->numind[population]*NUMCHROM);
      break;
   case 'n':
      if (op->panel)
         return(data->dnaptr->numseq[population] +
                op->numpanel[population]);
   case 's':
      return(data->dnaptr->numseq[population]);
      break;
   default:
      return(FLAGLONG);
      break;
}

return(FLAGLONG);

} /* getdata_numtips */


/*******************************************************************
 * getdata_sitecount() returns the contents of the sitecount array *
 * at the passed marker or FLAGLONG if it fails.                   */
long getdata_sitecount(option_struct *op, data_fmt *data, long marker)
{

switch(op->datatype)
{
   case 'a':
      break;
   case 'b':
   case 'm':
      break;
   case 'n':
   case 's':
      return(data->dnaptr->sitecount[locus][marker]);
      break;
   default:
      return(FLAGLONG);
      break;
}

return(FLAGLONG);

} /* getdata_sitecount */


/*********************************************************************
 * getdata_markersite() returns the contents of the markersite array *
 * at the passed marker or FLAGLONG if it fails.                     */
long getdata_markersite(option_struct *op, data_fmt *data, long marker)
{

switch(op->datatype)
{
   case 'a':
      break;
   case 'b':
   case 'm':
      break;
   case 'n':
   case 's':
      return(data->dnaptr->markersite[locus][marker]);
      break;
   default:
      return(FLAGLONG);
      break;
}

return(FLAGLONG);

} /* getdata_markersite */


/*******************************************************************
 * getdata_space() returns the length of the passed link stored in *
 * the general data structure or FLAGLONG if it fails.             */
double getdata_space(option_struct *op, data_fmt *data, long target)
{

switch(op->datatype)
{
   case 'a':
      break;
   case 'b':
   case 'm':
      return(data->msptr->mspace[population][locus][target]);
      break;
   case 'n':
      if (target == 0) {
         return(data->dnaptr->sspace[population][locus][target] +
 data->dnaptr->sspace[population][locus][getdata_nummarkers(op,data)]);
         break;
      }
   case 's':
      return(data->dnaptr->sspace[population][locus][target]);
      break;
   default:
      return(FLAGLONG);
      break;
}

return(FLAGLONG);

} /* getdata_space */


/***************************************************************
 * read_spacefile() reads the spacing info from the designated *
 * spacefile.                                                  *
 *                                                             *
 * entry "0" from a spacefile is read as the number of spaces  *
 * before the first marker.                                    */
void read_spacefile(FILE *file, option_struct *op, data_fmt *data)
{
  long numgaps, gapcnt, markernum;
  char *input,*tok;
  dnadata *dna;

  input = (char *)calloc(LINESIZE,sizeof(char));
  *input = '\0';

  dna = data->dnaptr;

  read_line(file,input);
  sscanf(input,"%ld%c[^\n]",&numgaps,&(dna->sdlm));

  for(gapcnt = 0; gapcnt < numgaps; gapcnt++) {
    read_line(file,input);
    tok = strtok(input,&dna->sdlm);
    markernum = atol(tok) - 1;
    /* if the entry is supposed to be the spacing before the first
       marker, then put it at its proper place */
    if (markernum == -1) markernum = getdata_nummarkers(op,data);
    dna->sspace[population][locus][markernum] = 
      (double)atof(strtok(NULL," "));
  }

} /* read_spacefile */


/*****************************************************************
 * longcmp() is a qsort helper for longs, thanks to Peter Beerli */
int longcmp (const void *v1, const void *v2)
{
  if (*(long *) v1 < *(long *) v2)
    {
      return -1;
    }
  else
    {
      if (*(long *) v1 > *(long *) v2)
        {
          return 1;
        }
      else
        return 0;
    }
} /* longcmp */


/*********************************************
 * doublecmp() is a qsort helper for doubles */
int doublecmp (const void *v1, const void *v2)
{
  if (*(double *) v1 < *(double *) v2)
    {
      return -1;
    }
  else
    {
      if (*(double *) v1 > *(double *) v2)
        {
          return 1;
        }
      else
        return 0;
    }
} /* doublecmp */


/********************************************************************
 * read_flipfile() reads the questionable site information from     *
 * the passed file.                                                 *
 *                                                                  *
 * read_flipfile() can only be called after treesetup().            *
 *                                                                  *
 * Expected file format of a flipfile:                              *
 *    1st line the number of individuals covered by the flipfile    *
 *                                                                  *
 *    then, for each individual in the file:                        *
 *       a line containing one of the following                     *
 *          "all" -> use if all the site are flippable.             *
 *                                                                  *
 *          numbers seperated by spaces -> use if you're specifying *
 *             which sites are flippable.                           *
 *             (e.g. "5 4 67 89 32 100" would mean that there are 5 *
 *             flippable sites with position numbers 4, 32, 67, 89  *
 *             and 100)                                             *
 *                                                                  *
 *          an exclamation point immediately followed by a list of  *
 *          numbers seperated by spaces -> use if you're specifying *
 *             which sites are NOT flippable.                       *
 *             (e.g. "!5 4 67 89 32 100" means that there are 5     *
 *             non-flippable sites with position numbers 4, 32, 67, *
 *             89 and 100)                                          *
 *                                                                  *
 *    the last line should just contain the word "end".             */
void read_flipfile(option_struct *op, data_fmt *data, tree *tr, FILE *file)
{
long i, j, whichcreature, numcreatures, site, numsites, temp, *tempa;
char *fileline, *fline, *num, char1;
creature *cr;

fscanf(file,"%ld",&numcreatures);
if (numcreatures != getdata_numtips(op,data)/NUMHAPLOTYPES) {
   fprintf(ERRFILE,"\nERROR--number of sequences in infile and number");
   fprintf(ERRFILE," of individuals in flipfile don't match!\n\n");
   exit(-1);
}

numsites = getdata_nummarkers(op,data);

for(whichcreature = 0; whichcreature < numcreatures; whichcreature++) {
   fline = fileline = (char *)calloc(LINESIZE,sizeof(char));
   read_line(file,fileline);
   char1 = fileline[0];
   if (char1 != 'a' && char1 != '!' && !isdigit(char1)) {
      fprintf(ERRFILE,"\nERROR--not a valid flipfile entry for line");
      fprintf(ERRFILE," containing:\n %s\n",fileline);
      exit(-1);
   }
   cr = &(tr->creatures[whichcreature]);
   switch(char1) {
      case 'a':
         cr->numflipsites = numsites;
         cr->flipsites = (long *)calloc(numsites,sizeof(long));
         for(site = 0; site < numsites; site++) cr->flipsites[site] = site;
         break;
      case '!':
         fileline++;
         temp = atol(strtok(fileline," "));
         cr->numflipsites = numsites-temp;
         tempa = (long *)calloc(temp,sizeof(long));
         for(i = 0; i < temp; i++) tempa[i] = atol(strtok(NULL," "));
         qsort((void *)tempa,temp,sizeof(long),longcmp);
         cr->flipsites = (long *)calloc(numsites-temp,sizeof(long));
         for(site = 0, i = 0, j = 0; site < numsites; site++, j++) {
            if (i < temp) {
               if (site == tempa[i]) {
                  j--;
                  i++;
                  continue;
               }
            }
            cr->flipsites[j] = site;
         }
         free(tempa);
         break;
      default: /* assume we're dealing with simply digits */
         cr->numflipsites = atol(strtok(fileline," "));
         cr->flipsites = (long *)calloc(cr->numflipsites,sizeof(long));
         for(site = 0; site < cr->numflipsites; site++) {
            num = strtok(NULL," ");
            if (!num) {
               fprintf(ERRFILE,"ERROR--in flipfile, creature %ld",
                  whichcreature);
               fprintf(ERRFILE," %ldth entry, failure to read.\n",
                  site+1);
               exit(-1);
            }
            cr->flipsites[site] = atol(num);
            /* cr->flipsites[site] = atol(strtok(NULL," ")); */
         }
         qsort((void *)cr->flipsites,cr->numflipsites,sizeof(long),longcmp);
         break;
   }
   free(fline);
}


} /* read_flipfile */


/*******************************************************************
 * pruneflips() removes "uninteresting" entries from the flippable *
 * list.  It also "verifies" the marker aliasing for likelihood.   *
 *                                                                 *
 * "uninteresting" = invariant                                     */
void pruneflips(option_struct *op, data_fmt *data, tree *tr)
{
long i, whichcreature, numcreatures, site, nummarkers, *fmarkers, temp,
  *siteptr;
creature *cr;
boolean noflips;

numcreatures = getdata_numtips(op,data)/NUMHAPLOTYPES;
nummarkers = getdata_nummarkers(op,data);
siteptr = data->siteptr;
fmarkers = (long *)calloc(nummarkers,sizeof(long));

for(whichcreature = 0; whichcreature < numcreatures; whichcreature++) {
   cr = &(tr->creatures[whichcreature]);
   for(site = cr->numflipsites-1; site > -1; site--) {
      if (isinvariant(op,data,cr->flipsites[site],cr)) {
         cr->numflipsites--;
         for(i = site; i < cr->numflipsites; i++) {
            cr->flipsites[i] = cr->flipsites[i+1];
         }
      } else fmarkers[cr->flipsites[site]] = 1;
   }
}

/* now to disable datalikelihood aliasing on all remaining flippable
markers-if removed then fliphap must do a complete redo of siteptr
array, not just rebuild_alias() */
for(site = 0; site < nummarkers; site++) {
    if (!fmarkers[site]) continue;

    temp = siteptr[site];
    siteptr[site] = 0;   
    for(i = site+1; i < nummarkers; i++) {
       if (siteptr[i] == site+1) {
          siteptr[i] = temp;
          break;
       }
    }
}

free(fmarkers);

noflips = TRUE;
for(whichcreature = 0; whichcreature < numcreatures; whichcreature++) {
   cr = &(tr->creatures[whichcreature]);
   if (cr->numflipsites != 0) noflips = FALSE;
}
if (noflips) {
   fprintf(ERRFILE,"\n\nAll sites specified as phase-unknown are ");
   fprintf(ERRFILE,"invariant, probable error in specifying flipfile.\n\n");
   exit(-1);
}

} /* pruneflips */


/*********************************************************
 * read_indname() reads the individual names of the data */
void read_indname(FILE *file, data_fmt *data, long pop, long lowcus,
   long ind, long nmlength, char datatype)
{
long i=0;

switch(datatype) {
   case 'a':
      break;
   case 'b':
   case 'm':
      while(i<nmlength) 
         data->msptr->indnames[pop][lowcus][ind][i++]=getc(file);
      data->msptr->indnames[pop][lowcus][ind][nmlength]='\0';
      break;
   case 'n':
   case 's':
      while(i<nmlength)
         data->dnaptr->indnames[pop][lowcus][ind][i++]=getc(file);
      data->dnaptr->indnames[pop][lowcus][ind][nmlength]='\0';
      break;
   default:
      fprintf(ERRFILE,"read_indname, can't get here!\n");
      break;
}

} /* read_indname */


/******************************************************************
 * read_microalleles() actually reads in the microsattellite data */
void read_microalleles(FILE *infile, msatdata *data, long pop, long ind)
{
char *input, dlm[2],ddlm[2], *a, *a1, *a2;
long lowcus,i;

input = (char *) calloc(1,sizeof(char)*(LINESIZE+1));
a = (char *) calloc(1,sizeof(char)*LINESIZE);
a1 = (char *) calloc(1,sizeof(char)*LINESIZE);
a2 = (char *) calloc(1,sizeof(char)*LINESIZE);
dlm[0]=data->dlm, dlm[1]='\0';
ddlm[0]=' ', ddlm[1]='\0';

read_line(infile,input);

for(lowcus=0; lowcus<data->numloci; lowcus++) {
   if(input[0]=='\0') read_line(infile,input);
   while(isspace((int) *input)) input++;
   i=0;
   while(input[i]!=' ' && input[i]!=dlm[0]) {
      a1[i]=input[i];
      i++;
   }
   a1[i]='\0';
   input += i;
   i=0;
   if(input[i]==dlm[0]){
      input++;
      while(input[i]!=' ' && input[i]!='\0') {
        a2[i]=input[i];
        i++;
      }
      a2[i]='\0';
      if(a2[0]=='\0') strcpy(a2,a1);
   input += i;
   } else strcpy(a2,a1);
   data->msats[pop][ind][lowcus][0] = atol(a1);
   data->msats[pop][ind][lowcus][1] = atol(a2);
}

free(a);
free(a1);
free(a2);
free(input);

} /* read_microalleles */


/*************************************************************
 * read_ind_seq() reads in the first line of a dna sequence, *
 * subsequent lines will be read by finish_read_seq()        */
long read_ind_seq(FILE *infile, dnadata *data, option_struct *op,
   long lowcus, long pop, long ind, long baseread)
{
long j;
char charstate;

j = (op->interleaved) ? baseread : 0;
charstate = getc(infile);
ungetc((int) charstate,infile);
while (j < data->sites[lowcus] && 
       !(op->interleaved && charstate=='\n')) {
   charstate = getc(infile);
   if (charstate == '\n') {
      if(op->interleaved) return j;
      else charstate=' ';
   }
   if (isspace(charstate) || isdigit(charstate))
        continue;
   charstate = uppercase(charstate);
   if ((strchr("ABCDGHKMNRSTUVWXY?O-",(int) charstate)) == NULL){
        fprintf(ERRFILE,"ERROR: BAD BASE: %c AT POSITION %5ld OF",
           charstate,j);
        fprintf(ERRFILE," INDIVIDUAL %3li in POPULATION %ld\n",ind,pop);
        exit(-1);
   }
   data->seqs[pop][lowcus][ind][j++] = charstate;
}

charstate = getc(infile); /* swallow the \n */
return j;

} /* read_ind_seq */


/***************************************************************
 * finish_read_seq() finishes reading in the dna data begun by *
 * read_ind_seq().                                             */
void finish_read_seq(FILE *infile, data_fmt *data, option_struct *op,
   long pop, long baseread)
{

long ind, baseread2=0, lowcus=0;

if(op->interleaved){
   while(baseread<data->dnaptr->sites[0]) {
      for(ind=0; ind < data->dnaptr->numseq[pop]; ind++) {
         baseread2 =
            read_ind_seq(infile,data->dnaptr,op,lowcus,pop,ind,baseread);
      }
      baseread = baseread2;
   }
}

for(lowcus = 1; lowcus < getdata_numloci(op,data); lowcus++){
   baseread=0;
/* "population" doesn't exist at this point so... */
   for(ind=0; ind < data->dnaptr->numseq[pop]; ind++) {
      read_indname(infile,data,pop,lowcus,ind,NMLNGTH,op->datatype);
      baseread =
         read_ind_seq(infile,data->dnaptr,op,lowcus,pop,ind,0);
   }
   if(op->interleaved){
      while(baseread<data->dnaptr->sites[lowcus]){
         for (ind=0; ind < data->dnaptr->numseq[pop]; ind++) {
             baseread2 =
                read_ind_seq(infile,data->dnaptr,op,lowcus,pop,ind,baseread);
         }
         baseread=baseread2;
      }
   }
}

} /* finish_read_seq */


/*****************************************************
 * read_popdata() reads in the data from file INFILE */
void read_popdata(FILE *infile, data_fmt *data, long pop,
   option_struct *op)
{
long ind, baseread=0, lowcus=0, numentries;

/* initialization for lint */
numentries = 0;

switch(op->datatype) {
   case 'a':
      break;
   case 'b':
   case 'm':
      numentries = data->msptr->numind[pop];
      break;
   case 'n':
   case 's':
      numentries = data->dnaptr->numseq[pop];
      break;
   default:
      fprintf(ERRFILE,"read_popdata: can't be here!\n");
      exit(-1);
}

for(ind=0; ind < numentries; ind++) {
   //read_indname(infile,data,pop,lowcus,ind,NMLNGTH,op->datatype);
   switch(op->datatype) {
   case 'a':
     /* read_alleles(infile, data,pop, ind); */
      break;
   case 'b':
   case 'm':
     /* if (data->dlm == '\0') read_alleles(infile, data,pop, ind);
      else */ read_microalleles(infile,data->msptr,pop,ind);
      break;
   case 'n':
   case 's':
      baseread=read_ind_seq(infile,data->dnaptr,op,lowcus,pop,ind,0L);
      break;
   default: 
      fprintf(stderr,"Wrong datatype, only the types a, m, s");
      fprintf(stderr," (electrophoretic alleles, \nmicrosatellite data,");
      fprintf(stderr,"sequence data) are allowed.\n");
      exit(-1);
      break;
   }
}
if(op->datatype=='s' || op->datatype=='n')
   finish_read_seq(infile,data,op,pop,baseread);

} /* read_popdata */


/***********************************************************************
 * setpopstuff() sets the data structure's population specific fields: *
 *    # of individuals per population, population title.               */
void setpopstuff(data_fmt *data, long pop, long numind, char *poptitle,
   char datatype)
{

switch(datatype) {
   case 'a':
      break;
   case 'b':
   case 'm':
      data->msptr->numind[pop] = numind;
      strcpy(data->msptr->popnames[pop],poptitle);
      break;
   case 'n':
   case 's':
      data->dnaptr->numseq[pop] = numind;
      strcpy(data->dnaptr->popnames[pop],poptitle);
      break;
   default:
      fprintf(ERRFILE,"setpopstuff, found an unknown datatype %c\n",
         datatype);
      exit(-1);
}

} /* setpopstuff */


/*****************************************************
 * setupdnadata() initializes the dna data structure */
void setupdnadata(dnadata **dna, long *numsites, long numseq,
   long numloci, long numpop)
{
long i, j, k, maxsites=0;

(*dna) = (dnadata *)calloc(1,sizeof(dnadata));

(*dna)->title = (char *)calloc(LINESIZE,sizeof(char));
(*dna)->popnames = (char **)calloc(numpop,sizeof(char *));
(*dna)->popnames[0] = (char *)calloc(numpop*LINESIZE,sizeof(char));
for(i = 1; i < numpop; i++)
   (*dna)->popnames[i] = (*dna)->popnames[0] + i*LINESIZE;

(*dna)->numpop = numpop;
(*dna)->numloci = numloci;

(*dna)->numseq = (long *)calloc(numpop,sizeof(long));
/* This is reset for all populations, except the first,
   by setpopstuff() */
for(i = 0; i < numpop; i++) (*dna)->numseq[i] = numseq;

(*dna)->sites = (long *)calloc(numloci,sizeof(long));
(*dna)->sitecount = (long **)calloc(numloci,sizeof(long *));
(*dna)->markersite = (long **)calloc(numloci,sizeof(long *));
for(i = 0; i < numloci; i++) {
   (*dna)->sites[i] = numsites[i];
   if (numsites[i] > maxsites) maxsites = numsites[i];
}
(*dna)->sitecount[0] = (long *)calloc(numloci*maxsites,sizeof(long));
(*dna)->markersite[0] = (long *)calloc(numloci*maxsites,sizeof(long));
for(i = 1; i < numloci; i++) {
   (*dna)->sitecount[i] = (*dna)->sitecount[0] + i*maxsites;
   (*dna)->markersite[i] = (*dna)->markersite[0] + i*maxsites;
}

(*dna)->seqs = (char ****)calloc(numpop,sizeof(char ***));
(*dna)->indnames = (char ****)calloc(numpop,sizeof(char ***));
(*dna)->sspace = (double ***)calloc(numpop,sizeof(double **));

(*dna)->seqs[0] = (char ***)calloc(numpop*numloci,sizeof(char **));
(*dna)->indnames[0] = (char ***)calloc(numpop*numloci,sizeof(char **));
(*dna)->sspace[0] = (double **)calloc(numpop*numloci,sizeof(double *));
for(i = 1; i < numpop; i++) {
   (*dna)->seqs[i] = (*dna)->seqs[0] + i*numloci;
   (*dna)->indnames[i] = (*dna)->indnames[0] + i*numloci;
   (*dna)->sspace[i] = (*dna)->sspace[0] + i*numloci;
}

(*dna)->seqs[0][0] = 
   (char **)calloc(numpop*numloci*numseq,sizeof(char *));
(*dna)->indnames[0][0] = 
   (char **)calloc(numpop*numloci*numseq,sizeof(char *));
/* use maxsites+1 in allocating sspace to include spacing before
   the first marker */
(*dna)->sspace[0][0] = 
   (double *)calloc(numpop*numloci*(maxsites+1),sizeof(double));
for(i = 0; i < numpop; i++)
   for(j = 0; j < numloci; j++) {
      (*dna)->seqs[i][j] = (*dna)->seqs[0][0] + i*numloci*numseq
         + j*numseq;
      (*dna)->indnames[i][j] = (*dna)->indnames[0][0] + i*numloci*numseq
         + j*numseq;
      (*dna)->sspace[i][j] = (*dna)->sspace[0][0] + i*numloci*(maxsites+1)
         + j*(maxsites+1);
   }

for(i = 0; i < numpop; i++)
   for(j = 0; j < numloci; j++) {
      for(k = 0; k < maxsites; k++)
         (*dna)->sspace[i][j][k] = 1.0;
/* initialize the most leftward space entry to zero, no additional sites
   before the first marker. */
      (*dna)->sspace[i][j][(*dna)->sites[i]] = 0.0;
   }

(*dna)->seqs[0][0][0] = 
   (char *)calloc(numpop*numloci*numseq*maxsites,sizeof(char));
(*dna)->indnames[0][0][0] = 
   (char *)calloc(numpop*numloci*numseq*(NMLNGTH+1),sizeof(char));
for(i = 0; i < numpop; i++)
   for(j = 0; j < numloci; j++)
      for(k = 0; k < numseq; k++) {
         (*dna)->seqs[i][j][k] = (*dna)->seqs[0][0][0]
             + i*numloci*numseq*maxsites + j*numseq*maxsites + k*maxsites;
         (*dna)->indnames[i][j][k] = (*dna)->indnames[0][0][0]
             + i*numloci*numseq*(NMLNGTH+1) + j*numseq*(NMLNGTH+1)
             + k*(NMLNGTH+1);
      }

(*dna)->dnaweight = (long *)calloc(maxsites,sizeof(long));
for(i = 0; i < maxsites; i++) (*dna)->dnaweight[i] = 1;

(*dna)->segranges0 = NULL;
(*dna)->segranges1 = NULL;
(*dna)->segranges2 = NULL;
(*dna)->segranges3 = NULL;

} /* setupdnadata */


/**********************************************
 * freednadata() frees the dna data structure */
void freednadata(dnadata *dna)
{
free(dna->popnames[0]);
free(dna->popnames);
free(dna->title);
free(dna->dnaweight);
free(dna->numseq);
free(dna->sites);
free(dna->sspace[0]);
free(dna->sspace);
free(dna->sitecount[0]);
free(dna->sitecount);
free(dna->markersite[0]);
free(dna->markersite);
free(dna->indnames[0][0][0]);
free(dna->indnames[0][0]);
free(dna->indnames[0]);
free(dna->indnames);
free(dna->seqs[0][0][0]);
free(dna->seqs[0][0]);
free(dna->seqs[0]);
free(dna->seqs);
free(dna->segranges0);
free(dna->segranges1);
free(dna->segranges2);
free(dna->segranges3);
free(dna);
} /* freednadata */


/*******************************************************************
 * printdnadata() prints the dna sequences contained in the passed *
 * data structure to the passed file in sequential form.           */
void printdnadata(option_struct *op, data_fmt *data, FILE *out)
{
long pop, numpop, loc, numloci, ind, numind, site, numsite, siteout;
dnadata *dna;

/* can't use standard functions to set all num* due to use of global
   "population" */
dna = data->dnaptr;
numpop = getdata_numpop(op,data);
for(pop = 0; pop < numpop; pop++) {
   fprintf(out,"----------------------------\n");
   fprintf(out,"Sequences for Population %ld\n",pop+1);
   fprintf(out,"----------------------------\n");
   numloci = getdata_numloci(op,data);
   for(loc = 0; loc < numloci; loc++) {
      fprintf(out,"Locus %ld -------------------\n",loc+1);
      numind = dna->numseq[pop];
      for(ind = 0; ind < numind; ind++) {
         fprintf(out,"%s",dna->indnames[pop][loc][ind]);
         numsite = dna->sites[loc];
         for(site = 0, siteout = NMLNGTH; site < numsite; site++) {
            fprintf(out,"%c",dna->seqs[pop][loc][ind][site]);
            siteout++;
            if (siteout == OUTLINESIZE) {fprintf(out,"\n"); siteout = 0;}
         }
         fprintf(out,"\n");
      }
   }
   fprintf(out,"\n");
}
fprintf(out,"\n\n");

} /* printdnadata */


/********************************************************************
 * initinvartips() sets up the special fractional tips for SNP data */
void initinvartips(option_struct *op, data_fmt *data,  long categs,
   tree *curtree)
{
long i, j, k, m, pos, n, numseq, numtips;

pos = getdata_nummarkers(op,data);
numseq = getdata_numseq(op,data);
numtips = getdata_numtips(op,data);

if (op->panel) {
   /* real tips */
   for(i = 1; i <= numseq; i++) {
      for (m = pos; m < pos+NUMINVAR; m++) {
         for(j=0; j < categs; j++) {
            for (n = 0; n < NUMSLICE; n++) {
               for (k = baseA; k <= baseT; k++) 
                  curtree->nodep[i]->x->s[m][j][n][k] = 1.0;
            }
         }
      }
   }
   /* fake tips */
   for (i=numseq+1; i <= numtips; i++) {
      for (m = 0; m < pos+NUMINVAR; m++) {
         for(j = 0; j < categs; j++) {
            for(k = baseA; k <= baseT; k++) {
               curtree->nodep[i]->x->s[m][j][0][k] = 1.0;
               for(n = 1; n < NUMSLICE; n++)
                  curtree->nodep[i]->x->s[m][j][n][k] = 0.0;
            }
            curtree->nodep[i]->x->s[m][j][1][baseA] = 1.0;
            curtree->nodep[i]->x->s[m][j][2][baseC] = 1.0;
            curtree->nodep[i]->x->s[m][j][3][baseG] = 1.0;
            curtree->nodep[i]->x->s[m][j][4][baseT] = 1.0;
         }
      }
   }
} else {
   for(i = 1; i <= numtips; i++) {
      for(j = 0; j < categs; j++) {
         for (m = pos; m < pos+NUMINVAR; m++) {
            for(k = baseA; k <= baseT; k++) {
               curtree->nodep[i]->x->s[m][j][0][k] = 0.0;
            }      
         }
         curtree->nodep[i]->x->s[pos][j][0][baseA] = 1.0;
         curtree->nodep[i]->x->s[pos+1][j][0][baseC] = 1.0;
         curtree->nodep[i]->x->s[pos+2][j][0][baseG] = 1.0;
         curtree->nodep[i]->x->s[pos+3][j][0][baseT] = 1.0;
      }
   }
}

} /* initinvartips */


/*************************************************************
 * makednavalues() set up fractional likelihoods at the tips */
void makednavalues(option_struct *op, data_fmt *data, long categs,
   tree *curtree)
{
  long i, k, l, m, numslice;
  long b;
  char ****y;

  y = data->dnaptr->seqs;
  numslice = (op->panel) ? NUMSLICE : 1L;

  for (k = 0; k < getdata_nummarkers(op,data); k++) {
    /* real tips are 1 to numseq */
    for (i = 1; i <= getdata_numseq(op,data); i++) {
      strcpy(curtree->nodep[i]->nayme,
          data->dnaptr->indnames[population][locus][i-1]);
      for (l = 0; l < categs; l++) {
        for (m = 0; m < numslice; m++) {
	  for (b = baseA; b <= baseT; b = b + 1)
	    curtree->nodep[i]->x->s[k][l][m][b] = 0.0;

  	switch (y[population][locus][i-1][k]) {
  
  	case 'A':
  	  curtree->nodep[i]->x->s[k][l][m][baseA] = 1.0;
  	  break;
  
  	case 'C':
  	  curtree->nodep[i]->x->s[k][l][m][baseC] = 1.0;
  	  break;
  
  	case 'G':
  	  curtree->nodep[i]->x->s[k][l][m][baseG] = 1.0;
  	  break;
  
  	case 'T':
  	  curtree->nodep[i]->x->s[k][l][m][baseT] = 1.0;
  	  break;
  
  	case 'U':
  	  curtree->nodep[i]->x->s[k][l][m][baseT] = 1.0;
  	  break;
  
  	case 'M':
  	  curtree->nodep[i]->x->s[k][l][m][baseA] = 1.0;
  	  curtree->nodep[i]->x->s[k][l][m][baseC] = 1.0;
  	  break;
  
  	case 'R':
  	  curtree->nodep[i]->x->s[k][l][m][baseA] = 1.0;
  	  curtree->nodep[i]->x->s[k][l][m][baseG] = 1.0;
  	  break;
  
  	case 'W':
  	  curtree->nodep[i]->x->s[k][l][m][baseA] = 1.0;
  	  curtree->nodep[i]->x->s[k][l][m][baseT] = 1.0;
  	  break;
  
  	case 'S':
  	  curtree->nodep[i]->x->s[k][l][m][baseC] = 1.0;
  	  curtree->nodep[i]->x->s[k][l][m][baseG] = 1.0;
  	  break;
  
  	case 'Y':
  	  curtree->nodep[i]->x->s[k][l][m][baseC] = 1.0;
  	  curtree->nodep[i]->x->s[k][l][m][baseT] = 1.0;
  	  break;
  
  	case 'K':
  	  curtree->nodep[i]->x->s[k][l][m][baseG] = 1.0;
	  curtree->nodep[i]->x->s[k][l][m][baseT] = 1.0;
	  break;

	case 'B':
	  curtree->nodep[i]->x->s[k][l][m][baseC] = 1.0;
	  curtree->nodep[i]->x->s[k][l][m][baseG] = 1.0;
	  curtree->nodep[i]->x->s[k][l][m][baseT] = 1.0;
	  break;

	case 'D':
	  curtree->nodep[i]->x->s[k][l][m][baseA] = 1.0;
	  curtree->nodep[i]->x->s[k][l][m][baseG] = 1.0;
	  curtree->nodep[i]->x->s[k][l][m][baseT] = 1.0;
	  break;

	case 'H':
	  curtree->nodep[i]->x->s[k][l][m][baseA] = 1.0;
	  curtree->nodep[i]->x->s[k][l][m][baseC] = 1.0;
	  curtree->nodep[i]->x->s[k][l][m][baseT] = 1.0;
	  break;

	case 'V':
	  curtree->nodep[i]->x->s[k][l][m][baseA] = 1.0;
	  curtree->nodep[i]->x->s[k][l][m][baseC] = 1.0;
	  curtree->nodep[i]->x->s[k][l][m][baseG] = 1.0;
	  break;

	case 'N':
	  for (b = baseA; b <= baseT; b++)
	    curtree->nodep[i]->x->s[k][l][m][b] = 1.0;
	  break;

	case 'X':
	  for (b = baseA; b <= baseT; b++)
	    curtree->nodep[i]->x->s[k][l][m][b] = 1.0;
	  break;

	case '?':
	  for (b = baseA; b <= baseT; b++)
	    curtree->nodep[i]->x->s[k][l][m][b] = 1.0;
	  break;

	case 'O':
	  for (b = baseA; b <= baseT; b++)
	    curtree->nodep[i]->x->s[k][l][m][b] = 1.0;
	  break;

	case '-':
	  for (b = baseA; b <= baseT; b++)
	    curtree->nodep[i]->x->s[k][l][m][b] = 1.0;
	  break;

        default:
          fprintf(ERRFILE,"**ERROR in setting up tips, probable error");
          fprintf(ERRFILE," in input data file.\nInterleaved/sequential");
          fprintf(ERRFILE," status may be incorrectly set.\n\n");
          exit(-1);
          break;
        } /* sorry about the indenting */
	}
      }
    }
  }
}  /* makednavalues */


/*************************************************************************
 * empiricaldnafreqs calculates empirical base frequencies from the data */
void empiricaldnafreqs(option_struct *op, data_fmt *data, tree *curtree)
{
  long i, j;
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
  for (i = 1; i <= getdata_numseq(op,data); i++) {
     for (j = 0; j < getdata_nummarkers(op,data); j++) {
	w = dna->dnaweight[j];
	temp = dna->freqa * curtree->nodep[i]->x->s[j][0][0][baseA];
	temp += dna->freqc * curtree->nodep[i]->x->s[j][0][0][baseC];
	temp += dna->freqg * curtree->nodep[i]->x->s[j][0][0][baseG];
	temp += dna->freqt * curtree->nodep[i]->x->s[j][0][0][baseT];
	suma += w * dna->freqa * curtree->nodep[i]->x->s[j][0][0][baseA] / temp;
	sumc += w * dna->freqc * curtree->nodep[i]->x->s[j][0][0][baseC] / temp;
	sumg += w * dna->freqg * curtree->nodep[i]->x->s[j][0][0][baseG] / temp;
	sumt += w * dna->freqt * curtree->nodep[i]->x->s[j][0][0][baseT] / temp;
     }
  }
  temp = suma + sumc + sumg + sumt;
  dna->freqa = suma / temp;
  dna->freqc = sumc / temp;
  dna->freqg = sumg / temp;
  dna->freqt = sumt / temp;

}  /* empiricaldnafreqs */


void getbasednafreqs(dnadata *dna, option_struct *op, double locus_ttratio,
  FILE *outfile)
{
  double aa, bb;

  putc('\n', outfile);
  if (op->freqsfrom)
    fprintf(outfile, "Empirical ");
  fprintf(outfile, "Base Frequencies:\n\n");
  fprintf(outfile, "   A    %10.5f\n", dna->freqa);
  fprintf(outfile, "   C    %10.5f\n", dna->freqc);
  fprintf(outfile, "   G    %10.5f\n", dna->freqg);
  fprintf(outfile, "  T(U)  %10.5f\n", dna->freqt);
  dna->freqr = dna->freqa + dna->freqg;
  dna->freqy = dna->freqc + dna->freqt;
  dna->freqar = dna->freqa / dna->freqr;
  dna->freqcy = dna->freqc / dna->freqy;
  dna->freqgr = dna->freqg / dna->freqr;
  dna->freqty = dna->freqt / dna->freqy;
  fprintf(outfile, "Transition/transversion ratio = %10.6f\n", locus_ttratio);
  aa = locus_ttratio * dna->freqr * dna->freqy - 
       dna->freqa * dna->freqg - dna->freqc * dna->freqt;
  bb = dna->freqa * dna->freqgr + dna->freqc * dna->freqty;
  dna->xi = aa / (aa + bb);
  dna->xv = 1.0 - dna->xi;
  dna->ttratio = dna->xi / dna->xv;
  if (dna->xi <= 0.0) {
    printf("WARNING: This transition/transversion ratio\n");
    printf("is impossible with these base frequencies!\n");
    dna->xi = 3.0 / 5;
    dna->xv = 2.0 / 5;
    fprintf(outfile, " Transition/transversion parameter reset\n\n");
  }
  fprintf(outfile, "(Transition/transversion parameter = %10.6f)\n",
	  dna->xi / dna->xv);
  dna->fracchange = dna->xi * (2 * dna->freqa * dna->freqgr + 
      2 * dna->freqc * dna->freqty) +
      dna->xv * (1.0 - dna->freqa * dna->freqa - 
      dna->freqc * dna->freqc - dna->freqg * dna->freqg - 
      dna->freqt * dna->freqt);
}  /* getbasednafreqs */


char uppercase(char ch)
{
return((islower((int)ch) ? (char)toupper((int)ch) : (ch)));

} /* uppercase */


/******************************************************************
 * inputdnaweights() reads the sitewise weights for the dna data. *
 * 0-9 and A-Z for weights 0-35.                                  */
void inputdnaweights(long numchars, dnadata *dna, option_struct *op)
{
char ch;
long i;

for (i = 0; i < numchars; i++) {
   do {
      if (eoln(weightfile)) {
         fscanf(weightfile, "%*[^\n]");
         getc(weightfile);
      }
      ch = getc(weightfile);
      if (ch == '\n') ch = ' ';
   } while (ch == ' ');

   dna->dnaweight[i] = 1;
   if (isdigit((int)ch))
      dna->dnaweight[i] = (long)(ch - '0');
   else if (isalpha((int)ch)) {
      ch = (char)uppercase(ch);
      dna->dnaweight[i] = (long)(ch - 'A' + 10);
   } else {
      printf("BAD WEIGHT CHARACTER: %c\n", ch);
      exit(-1);
   }
}

fscanf(weightfile, "%*[^\n]");
getc(weightfile);
op->weights = TRUE;
fclose(weightfile);

}  /* inputweights */


/************************************************
 * hapdist computes a distance between two sets *
 * of haplotype resolutions.                    */
long hapdist(option_struct *op,data_fmt *hap1, data_fmt *hap2, creature *critter)
{
long mkr, cr, numcreatures, markers;
long score[2];
long fullscore;

numcreatures = getdata_numtips(op,hap1)/NUMHAPLOTYPES;
markers = getdata_nummarkers(op,hap1);

fullscore = 0;
for (cr=0; cr<numcreatures; cr++) {
   score[0] = 0;
   score[1] = 0;
   for (mkr=0; mkr < markers; mkr++) {
/* matching the two 0 and the two 1 haplotypes */
     if (hap1->dnaptr->seqs[population][locus]
       [critter[cr].haplotypes[0]->number-1][mkr] !=
         hap2->dnaptr->seqs[population][locus]
       [critter[cr].haplotypes[0]->number-1][mkr])
       score[0]++;
     if (hap1->dnaptr->seqs[population][locus]
       [critter[cr].haplotypes[1]->number-1][mkr] !=
         hap2->dnaptr->seqs[population][locus]
       [critter[cr].haplotypes[1]->number-1][mkr])
       score[0]++;
/* matching 0-1 and 1-0 haplotypes */
     if (hap1->dnaptr->seqs[population][locus]
       [critter[cr].haplotypes[0]->number-1][mkr] !=
         hap2->dnaptr->seqs[population][locus]
       [critter[cr].haplotypes[1]->number-1][mkr])
       score[1]++;
     if (hap1->dnaptr->seqs[population][locus]
       [critter[cr].haplotypes[1]->number-1][mkr] !=
         hap2->dnaptr->seqs[population][locus]
       [critter[cr].haplotypes[0]->number-1][mkr])
       score[1]++;
  }
  if (score[0]<score[1]) fullscore += score[0];
  else fullscore += score[1];
}
return(fullscore);

} /* hapdist */

/********************************************************
 * This reads in the "true" haplotypes for the purposes *
 * of comparing them with derived ones.  Input file is  *
 * truehaps.                                            */
void readtruehaps(option_struct *op, data_fmt *hap, long numpop, long
  numloci, long numind, long *numsites)
{
char *s;
FILE *truehaps;

s = (char *)calloc(LINESIZE,sizeof(char));

truehaps = fopen("truehaps","r");
/* read off first 3 lines of standard input file as unneeded */
read_line(truehaps,s); read_line(truehaps,s); read_line(truehaps,s);

setupdata(hap,op->datatype,numpop,numloci,numind,numsites,"True Haplotypes");
setpopstuff(hap,0L,numind,"poptitle",op->datatype);
read_popdata(truehaps,hap,0L,op);

free(s);

} /* readtruehaps */

/********************************************************
 *  This copies a set of haplotypes into a new data     *
 *  structure which it allocates.                       */
void copyhaps(option_struct *op, data_fmt *oldhap, data_fmt *newhap,
  long numpop, long numloci, long numind, long *numsites)
{
long pop, seq, site, loc;

setupdata(newhap,op->datatype,numpop, numloci, numind,
  numsites, "Copy Haplotypes");
setpopstuff(newhap,0L,numind,"poptitle",op->datatype);


for (pop = 0; pop < numpop; pop++) {
  for (loc = 0; loc < numloci; loc++) {
    for (seq = 0; seq < numind; seq++) {
      for (site = 0; site < numsites[locus]; site++) { 
        newhap->dnaptr->seqs[pop][locus][seq][site] = 
        oldhap->dnaptr->seqs[pop][locus][seq][site];
      }
    }
  }
}

} /* copyhaps */
