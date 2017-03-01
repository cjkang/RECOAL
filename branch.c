
void upateg(node *eventnode, node *sourcenode, brlilst *targetbr1, brlist *targetbr2, int case){
  long i, start, end;
  double tyme;

  eventnode = findunique(eventnode);
  sourcenode = findunique(sourcenode);
 
  if (case == 0){ // new rec case at active branch
    tyme = eventnode->tyme - eventnode->back->tyme;
    for (i = eventnode->coal[1]; i < eventnode->coal[eventnode->coal[0] * 2]; i++)  gpartweight[i] += tyme;

    if (eventnode->next->recstart == 0) recpos = eventnode->next->recend;
    else recpos = eventnode->next->recstart - 1;
    
    gpartrec[recpos] ++;
    return;
  }
  if (case == 1){ // new coal case - coalesce of two active branches
    tyme = (eventnode->tyme - eventnode->next->back->tyme) + (eventnode->tyme - eventnode->next->next->back->tyme);
    
    targetnode = findunique(eventnode->next->back);
    for (i = targetnode->coal[1]; i < targetnode->coal[targetnode->coal[0] * 2]; i++)  gpartweight[i] += tyme;
    targetnode = findunique(eventnode->next->next->back);
    for (i = targetnode->coal[1]; i < targetnode->coal[targetnode->coal[0] * 2]; i++)  gpartweight[i] += tyme;

    return;
  }
  if (case == 2){ // new coal case - coalesce of one active with inactive.
    tyme = eventnode->tyme - sourcenode->tyme;
    
    for (i = sourcenode->coal[1]; i < sourcenode->coal[sourcenode->coal[0] * 2]; i++)  gpartweight[i] += tyme;
    
    if (targetbr1 != NULL){

      tyme = targetbr1->endtyme - targetbr1->starttyme;
    
      for (i = targetbr1->nfsite; i < targetbr1->ofsitel; i ++)  gpartweight[i] += tyme;
      for (i = targetbr1->olsite; i < targetbr1->nlsitel; i ++)  gpartweight[i] += tyme;
    }
    
    return;
  }
  if (case == 3){ // new rec at brlist
    tyme = targetbr1->endtyme - targetbr1->starttyme;
    
    for (i = targetbr1->nfsite; i < targetbr1->ofsitel; i ++)  gpartweight[i] += tyme;
    for (i = targetbr1->olsite; i < targetbr1->nlsitel; i ++)  gpartweight[i] += tyme;

    gpartrec[brs->cross] ++;

    return;
  }
  if (case == 4){ // no new event - upate changes
    
    if (targetbr1 != NULL){
      tyme = targetbr1->endtyme - targetbr1->starttyme;
      for (i = targetbr1->nfsite; i < targetbr1->ofsitel; i ++)  gpartweight[i] += tyme;
      for (i = targetbr1->olsite; i < targetbr1->nlsitel; i ++)  gpartweight[i] += tyme;
    }
    if (targetbr2 != NULL){
      tyme = targetbr2->endtyme - targetbr2->starttyme;
      for (i = targetbr2->nfsite; i < targetbr2->ofsitel; i ++)  gpartweight[i] += tyme;
      for (i = targetbr2->olsite; i < targetbr2->nlsitel; i ++)  gpartweight[i] += tyme;
    }

    return;
  }
}

    
  
