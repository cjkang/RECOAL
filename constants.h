/* constants.h -- constants for RECOMBINE.  Changes to this file will
   take effect only when the program is recompiled.

   Some of these constants (listed first) can usefully be changed to
   adapt the program to conditions.  A few may have to be changed if 
   you use a computer with less than 32-bit words.  The remainder
   should not be changed.
*/

/***************************************************/
/* The following constants can usefully be changed */
/***************************************************/

/**** Input format constants ****/
#define NMLNGTH     0L       /* required # of characters in sequence name */
#define MENU       TRUE      /* put up an interactive menu at start? */
#define TRUTHKNOWN FALSE     /* for a mapping run, is the truth known,
                                and an appropiate "truetrait" file
                                present? */

/**** Runtime constants ****/
#define ITERS      100       /* how many times to call the NR curve
                                maximiser at the end of each chain */
#define REC_THRESHOLD 1.0    /* minimum allowable # of recombinations in
                                a proposed tree */
#define RECOMB_MAX  100000      /* maximum number of recombinations allowed
                                in a tree */
#define TWIDDLE_PROB 0.00    /* chance we will re-evaluate an old
                                recombination rather than doing a normal
                                rearrangement */
#define FLIP_PROB    0.00    /* chance we will flip an individual's
                                haplotype assignment rather than doing a
                                normal rearrangement.
                                NOTE: FLIP_PROB will be treated as 0.0
                                if haplotyping has not been selected.  */
#define FRACTRECOMB  1.00    /* the number of recombinations that are
                                added to a chain that records no
                                recombinations in its sampled trees. */
#define NUMBURNIN    2000 /* the number of initial trees to discard
                                at the beginning of each chain prior to
                                sampling */

#define INIBURNIN    10000

/**** Output format constants ****/
#define ONEBESTREE TRUE      /* "TRUE" = only single best tree in file "bestree"
                                "FALSE" = running tally of best trees in file
                                "bestree" */
#define LINESIZE    10240    /* this is actually (maximum linesize)+1
                                because of trailing '\0' */
#define STEPMAX    1000      /* the maximum number of steps that Newton
                                Raphson will attempt to do */
#define STEPPRINT    8       /* this number should be greater than the
                                order, base 10, of the greatest number
                                of steps that will be attempted by the
                                Newton-Raphson maximizer (STEPMAX) */
#define DISTRIBUTION FALSE   /* should a summary table of the number
                                of recombinations present in the trees
                                of the last run chain be printed into
                                the outfile? */

/**** Microsattellite parameters ****/
#define MICRO_ALLELEMAX 200  /* the maximum number of differing alleles
                                that we distinguish in microsats */
#define MICRO_MAXCHANGE 10   /* the maximum number of differences that
                                will be examined in the likelihood
                                model */
#define NUMCHROM 2L          /* the number of entries per microsat
                                marker--"ploidy" */

/***** disease trait likelihood parameters ******/
#define NUMTRAIT 2L          /* the number of distinct, scored, types
                                present in the data */

/***** haplotyping parameters ******/
#define NUMHAPLOTYPES 2L     /* the number of haplotypes present in the
                                data for each individual, MUST be the
                                same for every individual in the data! */
#define NUMHAPFLIP  1L       /* the number of sites that are flipped
                                during each fliphap() attempt */
#define TRUEHAP  FALSE       /* are we reading true haplotypes from file? */     

/* histogram output control constants */
#define HIST_MAXHT 60 /* maximum number of characters high */
#define HIST_MAXWD 60 /* maximum number of characters wide */
#define LEFT_MARGIN 5 /* width of left margin in characters */
#define YOUT 5 /* interval spacing for y-axis labelling */
#define XOUT 5 /* interval spacing for x-axis labelling, must be
                  bigger than the number of digits in the "number
                  of sites in the data." */


/**********************************************************************/
/* The following constants may have to be changed for small computers */
/**********************************************************************/

/**** Arithmetic precision constants ****/
#ifndef DBL_MAX
   #define DBL_MAX 99e15
#endif
#define epsilon    0.0000001     /* a small number, less than epsilon1,
                                    also used as minimum parameter
                                    values in BFGS */
#define epsilon1   0.00001       /* a small enough number? */
#define EXPMIN    -40.0          /* minimum value of "x" in "exp(x)" */
#define EXPMAX     200.0         /* maximum value of "x" in "exp(x)" */
#define EPSILON    99e-15        /* a really small number, also used
                                    as minimum parameter values in NR */
#define DBL_EPSILON 2.2204460492503131e-16 /* another small number */
#define POSMAX       DBL_MAX     /* largest positive number */
#define NEGMAX  -(DBL_MAX - 1.0) /* smallest negative number */

/*************************************************/
/* The following constants should not be changed */
/*************************************************/

#define VERSION_NUM   0.99       /* current version number */
#define GROWTHUSED    FALSE      /* are growth-enabled summary trees
                                    in use? */
#define DF2  2.0                 /* critical value for 2 DF */
#define REPEAT_TOLERANCE 5       /* how many times (sequentially) should the
                                    exact same prior-likelihood value be
                                    encountered in point estimation before
                                    giving up. */
#define ALIASING   TRUE          /* are we using sitewise aliasing? */
#define ROOTNUM        0         /* the number of the root in the
                                    "nodep" array of nodes */
#define ROOTLENGTH   10000       /* length of root branch must be much
                                    larger than the total height of the
                                    tree */
#define NUMBOOL       22         /* number of boolean tokens in parameter
                                    file */
#define NUMNUMBER     18         /* number of numeric (double) tokens in
                                    parameter file */
#define NUM_TYPE_CHAINS 2        /* different lengths of chains */
#define NUMPARAMETERS   2        /* number of parameters of interest */
#define XARRAYSIZE 100           /* number of bins for recycling x arrays */
#define ALWAYS_ACCEPT   0        /* should we always accept all proposed
                                    changes to the tree? */
#define ALWAYS_REJECT   0        /* should we always reject all proposed
                                    changes to the tree?
                                    NOTE: REJECT will supersede ACCEPT
                                    when it is set TRUE.
                                    NOTE2: only use with 1 long chain,
                                    sampling increment 1, and NO short
                                    chains */
#define ARBITRARY_THETA 0.1      /* arbitrary (unused) theta value used
                                    in ALWAYS_REJECT workaround */
#define FLAGLONG       -99L      /* a arbitrary special value */
#define BOGUSTREETYME  -99.0     /* an arbitrary impossible tyme for a
                                    tree node */
#define NUMINVAR        4        /* number of invariant sites used
                                    in SNP data-likelihood calculation */
#define NUMSLICE        5        /* number of different ways of
                                    calculating P(D|G) in panel-SNP case
                                  */
#define LOG2 0.693147180559945309417232121458 /*N[Log[2] ,30]*/
#define THETAMIN   EPSILON       /* minimum value that theta may attain
                                    before final estimation */
