struct treeinfo{
  long numrec_hot, numrec_cold;
  double weight_hot, weight_cold;
  long numcoal;
  double branch;
  long HH,HC,CC,CH;
}
struct parmset{
}
double gettreeweight(treeinfo *treerec, parmset *parm_n, parmset *parm_not){

  prob_r = pow((parm_n->recH / parm_not->recH), treerec->numrec_hot) * exp(weight_hot * (parm_n->recH - parm_not->recH)) 
    * pow((parm_n->recC / parm_not->recC), treerec->numrec_cold) * exp(weight_cold * (parm_n->recC - parm_not->recC)) ;
  
  prob_theta = pow((parm_not->theta / parm_n->theta), treeinfo->numcoal) * exp(treeinfo->branch * (1.0/parm_n->theta - 1.0/parm_not->theta));
  
  prob_lamda = pow((parm_n->lamda_hot/parm_not->lamda_hot), treeinfo->HH) * pow((1 - parm_n->lamda_hot) / (1 - parm_not->lamda_hot), treeinfo->HC)
    * pow((parm_n->lamda_cold/parm_not->lamda_cold), treeinfo->CC) * pow((1 - parm_n->lamda_cold) / (1 - parm_not->lamda_cold), treeinfo->CH);

  return prob_r * prob_theta * prob_lamda;
}

parmset *EM(treeinfo **treesum, parmset *parm0, long numtree){
  parmset *parmn;
  double sumweight_hot,sumweight_cold;
  double sumrec_hot,sumrec_cold;
  double sumbranch, sumcoal;
  double sumCC,sumCH,sumHH,sumHC;
  long i;
  double com_L, com_L_old = 0;

  tree = *curtree;
  recnumb = **sumtree;
  recnum = calloc();

  parmn->theta = parm0->theta;
  parmn->recH = parm0->recH;
  parmn->recC = parm0->recC;
  parmn->lamda_hot = parm0->lamda_hot;
  parmn->lamda_cold = parm0->lamda_cold;
  
  while (1){
    sumweight_hot = sumweight_cold = sumrec_hot = sumrec_cold = sumbranch = sumcoal = sumCC = sumCH = sumHH = sumHC = 0;
    com_L = com_L_old = 0;
    for (i = 1; i<numtree;i++){
      weight = getttreeweight(tree_sum, parmn, parm0);
      com_L_old += weight;
      sumbranch += weight * treesum[i]->branch;
      sumcoal += weight * treesum[i]->numcoal;
    }
    pamrn->theta = sumbranch/sumcoal;

    for (i = 1; i<numtree;i++){
      weight = getttreeweight(tree_sum, parmn, parm0);
      sumweight_hot += weight * treesum[i]->weight_hot;
      sumrec_hot += weight * treesum[i]->numrec_hot;
    }
    parmn->recH = sumrec_hot/sumweight_hot;

    for (i = 1; i<numtree;i++){
      weight = getttreeweight(tree_sum, parmn, parm0);
      sumweight_cold += weight * treesum[i]->weight_cold;
      sumrec_cold += weight * treesum[i]->numrec_cold;
    }
    parmn->recC = sumrec_cold/sumweight_cold;
    for (i = 1; i<numtree;i++){
      weight = getttreeweight(tree_sum, parmn, parm0);
      sumHH += weight * treesum[i]->HH;
      sumHC += weight * treesum[i]->HC;
    }
    parmn->lamda_hot = sumHH / (sumHH + sumHC);

    for (i = 1; i<numtree;i++){
      weight = getttreeweight(tree_sum, parmn, parm0);
      sumCC += weight * treesum[i]->CC;
      sumCH += weight * treesum[i]->CH;
    }
    parmn->lamda_hot = sumCC / (sumCC + sumCH);


    for (i = 1; i<numtree;i++){
      weight = getttreeweight(tree_sum, parmn, parm0);
      com_L += weight;
    }
    
    if (com_L < com_L_old) printf ("error in EM\n");
    if (com_L == com_L_old) return parmn;
    
  }
}

