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
  recnumb *recs
  double llike;
} treerec;

/***********************************************************
 * scoretree saves the values necessary for evaluating the *
 * model_likelihood for a single tree in the 'sum' array   */
void scoretree(option_struct *op, data_fmt *data, long chain)
{
  tlist *t;
  treerec *trii;
  recnumb *recnums,*cur,*temp,*curtemp;
  long i, j, refchain, chaintype, entries;
  double temp,tyme;

  temp = malloc(sizeof(recnumb));
  recnums = temp;
  cur = recnum0;
  temp->beg = cur->beg;
  temp->end = cur->end;
  temp->next = NULL;
  do{
    if (cur->next == NULL) break;
    cur = cur->next;
    curtemp = malloc(sizeof(recnumb));
    temp->next = curtemp;
    temp = temp->next;
    temp->beg = cur->beg;
    temp->end = cur->end;
    temp->next = NULL;
  } 
    
 


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
  if(trii->eventtype) free(trii->eventtype);
  trii->eventtype = (long *)calloc(1,entries*sizeof(long));
  if (trii->sitescore) free(trii->sitescore);
  /* the +1 is to guarantee a minimum allocation of size 1 */
  trii->sitescore = (long *)calloc(curtree->numrecombs+1,sizeof(long));

  trii->tk = 0.0;
  trii->ts = 0.0;

  for(i=0,j=0,temp=0.0;i<entries;i++) {
    if (t->numbranch == 1) break;
    tyme = t->age - t->eventnode->tyme;
    trii->eventtype[i] = isrecomb(t->eventnode);
    temp += t->numbranch * (t->numbranch - 1) * tyme;
    trii->tk = temp;
    trii->ts += count_active_tlist(op,data,t) * tyme;
    tsr += count_weighted_active_tlist(op,data,t) * tyme;
    addrec(op,data,t,tyme, &recnums);

    if (trii->eventtype[i]) {
      trii->sitescore[j] = findlink(t->eventnode);
      j++;
    }
    t = t->succ;
  }

/*   for (t = tr->tymelist; t->succ != NULL; t=t->succ) { */
/*     tyme = t->age - t->eventnode->tyme; */
/*     tk += t->numbranch * (t->numbranch - 1.0) * tyme; */
/*     tsr += count_weighted_active_tlist(op,data,t) * tyme; */
/*     addrec(op,data,t,tyme, &recnums); */
/*   } */

/*   prob = -tk/th - tsr + (log(2.0/th))*tr->numcoals; */
/*   for (cur = recnums; cur->next !=NULL ; cur = cur->next){ */
/*     prob += log(findrate(cur->site))*cur->numrecs; */
/*   } */
/*   return(prob); */

  trii->numcoals = curtree->numcoals;
  trii->numrecombs = curtree->numrecombs;
  trii->llike = curtree->likelihood;

  trii->recs = recums;
#if !ALWAYS_REJECT
  if (temp == 0.0) fprintf(ERRFILE,"WARNING:  Tree has become length zero\n");
#endif
  if (trii->numcoals - trii->numrecombs != getdata_numtips(op,data) - 1) {
    fprintf(ERRFILE,"ERROR:scoretree says: bad nodes in the tree!\n");
    exit(-1);
  }

}  /* scoretree */
