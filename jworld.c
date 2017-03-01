/*------------------------------------------------------
 Maximum likelihood estimation 
 of migration rate  and effectice population size
 using a Metropolis-Hastings Monte Carlo algorithm                            
 -------------------------------------------------------                        
 W O R L D   R O U T I N E S 

 creates tree structures,
 reads tree [has to be done
 
 prints results,
 and finally helps to destroy itself.
                                                                                                               
 Peter Beerli 1996, Seattle
 beerli@genetics.washington.edu
 $Id$

 modifed for use with recombine by Jon Yamato (1/8/97)
-------------------------------------------------------*/

#include "world.h"

#ifndef REC_MODELLIKE_INCLUDE
#include "rec_modellike.h"
#endif

#ifdef DMALLOC_FUNC_CHECK
#include "dmalloc.h"
#endif

#ifdef DMEMDEBUG
#define MEMDEBUG
#include "memdebug.h"
#endif

#define PLANESIZE 36
#define PLANEBIGTICKS 6
#define PLANEBIGTICKVALUES {-3, -2, -1, 0, 1, 2}
#define PLANETICKS   {0.001, 0.0013895, 0.0019307, 0.0026827, 0.00372759, \
                     0.00517947, 0.00719686, 0.01, 0.013895, 0.019307, \
                     0.026827, 0.0372759, 0.0517947, 0.0719686, 0.1, \
                     0.13895, 0.19307, 0.26827, 0.372759, 0.517947, \
                     0.719686, 1., 1.3895, 1.9307, 2.6827, 3.72759, \
                     5.17947, 7.19686, 10., 13.895, 19.307, 26.827, 37.2759, \
                     51.7947, 71.9686, 100.}
#define CONTOURLEVELS 8
#define CONTOURS_LOCUS {0.0,-1.38629/2., -5.99147/2., -9.21034/2.,\
                        0.0,-1.38629/2., -5.99147/2., -9.21034/2.}
#if 0
#define CONTOURS_LOCUS {0.0,-3.35670/2., -9.48773/2., -13.2767/2.,\
                        0.0,-3.35670/2., -9.48773/2., -13.2767/2.}
#define CONTOURS_LOCI  {0.0,-4.35146/2., -11.0705/2., -15.0863/2.,\
                        0.0,-4.35146/2., -11.0705/2., -15.0863/2.}
#endif

#define PAGEFEED fprintf(outfile,"\n\n\n\n");

extern long totchains;
extern FILE *outfile;

extern double model_likelihood(option_struct *op, data_fmt *data, double theta,
   double rec, long lowcus, long chain);
extern double fmodel_likelihood(option_struct *op, data_fmt *data, double theta,
   double rec, double **denom, long lowcus, long chain);


void plot_surface(option_struct *op, data_fmt *data, long lowcus, char ***plane, long x)
{
   long i;
   long loci;

   loci = getdata_numloci(op,data);

   fprintf(outfile, "\nLog-Likelihood surface\n");
   fprintf(outfile, "---------------------------------------\n\n");
   fprintf(outfile, "Legend:\n");
   fprintf(outfile, "   X = Maximum likelihood\n");
   fprintf(outfile, "   * = in approximative 50%% confidence limit\n");
   fprintf(outfile, "   + = in approximative 95%% confidence limit\n");
   fprintf(outfile, "   - = in approximative 99%% confidence limit\n");
   if (lowcus != -1) { /* single locus curve */
	fprintf(outfile, "\n\nLocus %li\n(x-axis= Theta\n", lowcus + 1);
	fprintf(outfile, " y-axis = Recombination-rate,\nunits = log10)\n");
	for (i = x+1; i >= 0; i--)
		 fprintf(outfile, "%42.42s\n", plane[lowcus][i]);
	fprintf(outfile, "%42.42s\n", plane[lowcus][x+2]);
	PAGEFEED;
   } else {
	fprintf(outfile, "\nOver all loci\n\n");
	lowcus = loci;
	fprintf(outfile, "(x-axis= recombination-rate,");
	fprintf(outfile, " y-axis = Theta, units = log10)\n");
	for (i = x; i >= 0; i--)
		fprintf(outfile, "%42.42s\n", plane[lowcus][i]);
	PAGEFEED;
   }
} /* plot_surface */

#if 0
void create_loci_plot(world_fmt * world, char **plane, timearchive_fmt * atl, long loci)
{
   long intervals = PLANESIZE;
   long i, g = 1;
   nr_fmt *nr;
   double **pl;
   double contours[CONTOURLEVELS] = CONTOURS_LOCI;

   nr = (nr_fmt *) calloc(1, sizeof(nr_fmt));
   for (i = 1; i < loci + 1; i++) {
	  if (g < atl[i].T)
		 g = atl[i].T;
   }
   create_nr(nr, atl[loci].numpop, g);
   pl = (double **) calloc(1, sizeof(double *) * intervals);
   for (i = 0; i < intervals; i++) {
	  pl[i] = (double *) calloc(1, sizeof(double) * 2 * intervals);
   }
   calc_loci_plane(world, nr, atl, pl, loci, contours);
   fill_plotplane(plane, pl, contours);
   destroy_nr(nr);
   print_mathematica(world, pl, intervals, intervals);
   for (i = 0; i < intervals; i++) {
	  free(pl[i]);
   }
   free(pl);
}
#endif


void create_locus_plot(option_struct *op, data_fmt *data, long lowcus, char **plane)
{
   long intervals = PLANESIZE;
   long i;
   double **pl;
   double contours[CONTOURLEVELS] = CONTOURS_LOCUS;

   pl = (double **) calloc(1, sizeof(double *) * intervals);
   for (i = 0; i < intervals; i++) {
	  pl[i] = (double *) calloc(1, sizeof(double) * 2 * intervals);
   }
   calc_locus_plane(op,data,lowcus,pl,contours);
   fill_plotplane(plane, pl, contours);
   print_mathematica(pl, intervals, intervals);
   for (i = 0; i < intervals; i++) {
	  free(pl[i]);
   }
   free(pl);
}


/*private functions========================================== */
/* ---------------------------------------------------------
   creates memory for archive of timelists for each locus
   the "locus 0" is reserved for for the current locus, this
   fake locus is use to speed up the NR-estimation for all chains
   whereas the loci 1 .. n are used to save the results of the last
   chain for the combined estimate at the end of the program */

void create_plotplane(option_struct *op, data_fmt *data, char**** plane)
{
   short locus, i;
   long numloci;

   numloci = getdata_numloci(op,data);
   (*plane) = (char ***) calloc(1, sizeof(char **) * (numloci + 1));
   for (locus = 0; locus < numloci + 1; locus++) {
	  (*plane)[locus] = (char **) calloc(1, sizeof(char *)
						 * (PLANESIZE + 3));
	  for (i = 0; i < PLANESIZE + 3; i++) {
		 (*plane)[locus][i] = (char *) calloc(1, sizeof(char)
											 * (PLANESIZE + PLANESIZE + 20));
	  }
   }
}

#if 0
void calc_loci_plane(world_fmt * world, nr_fmt * nr, timearchive_fmt * atl, double **pl, long loci, double contours[])
{
   long i, j;
   double max1 = -DBL_MAX;
   double max2 = -DBL_MAX;
   double values[PLANESIZE] = PLANETICKS;
   long intervals = PLANESIZE;
   if (world->options->gamma) {
	  nr->param[4] = world->atl[loci + 1].param[4];
   }
   for (i = 0; i < intervals; i++) {
	  nr->param[0] = values[i];
	  if (world->options->gamma)
		 calc_gamma(nr);
	  for (j = 0; j < intervals; j++) {
		 nr->param[0] = values[i];
		 nr->param[1] = world->atl[loci + 1].param[1];
		 nr->param[2] = values[j] / values[i];
		 nr->param[3] = world->atl[loci + 1].param[3];
		 calc_loci_like(nr, atl, loci, world->options->gamma);
		 pl[i][j] = nr->llike;
		 if (max1 < nr->llike)
			max1 = nr->llike;
	  }
   }
   nr->param[0] = world->atl[loci + 1].param[0];
   if (world->options->gamma)
	  calc_gamma(nr);
   for (i = 0; i < intervals; i++) {
	  for (j = 0; j < intervals; j++) {
		 nr->param[0] = world->atl[loci + 1].param[0];
		 nr->param[1] = values[i];
		 nr->param[2] = world->atl[loci + 1].param[2];
		 nr->param[3] = values[j] / values[i];
		 calc_loci_like(nr, atl, loci, world->options->gamma);
		 pl[i][j + intervals] = nr->llike;
		 if (max2 < nr->llike)
			max2 = nr->llike;
	  }
   }
   contours[0] += max1;
   contours[1] += max1;
   contours[2] += max1;
   contours[3] += max1;
   contours[4] += max2;
   contours[5] += max2;
   contours[6] += max2;
   contours[7] += max2;
}
#endif

void calc_locus_plane(option_struct *op, data_fmt *data, long lowcus, double **pl,
   double contours[])
{
   long intervals = PLANESIZE;
   long lastchain = totchains - 1;
   long i, j, chaintype, numloci;
   double max1 = -DBL_MAX;
   double max2 = -DBL_MAX;
   double values[PLANESIZE] = PLANETICKS;
   double **llike0;

   chaintype = TYPE_CHAIN(lastchain);
   numloci = getdata_numloci(op,data);

   llike0 = (double **)calloc(1,numloci*sizeof(double *));
   llike0[0] = (double *)
      calloc(1,numloci*op->numout[chaintype]*sizeof(double));
   for(i = 1; i < numloci; i++)
      llike0[i] = llike0[0] + i*op->numout[chaintype];

   precalc_llike0(op,data,lowcus,lastchain,llike0);

   for (i = 0; i < intervals; i++) {
	  for (j = 0; j < intervals; j++) {
                 pl[i][j] =
                    fmodel_likelihood(op,data,values[j],values[i],llike0,
                       lowcus,lastchain);
		 if (max1 < pl[i][j])
			max1 = pl[i][j];
	  }
   }
#if 0
/* this stuff is for a second population */
   for (i = 0; i < intervals; i++) {
	  for (j = 0; j < intervals; j++) {
		 nr->param[0] = param[0];
		 nr->param[1] = values[i];
		 nr->param[2] = param[2];
		 nr->param[3] = values[j] / values[i];
		 calc_like(nr, tl, G);
		 pl[i][j + intervals] = nr->llike;
		 if (max2 < nr->llike)
			max2 = nr->llike;
	  }
   }
#endif
   contours[0] += max1;
   contours[1] += max1;
   contours[2] += max1;
   contours[3] += max1;
   contours[4] += max2;
   contours[5] += max2;
   contours[6] += max2;
   contours[7] += max2;

   free(llike0[0]);
   free(llike0);

} /* calc_locus_plane */

void fill_plotplane(char **plane, double **pl, double contours[])
{
   long i, j, zz = 0;

   char line[100];
   long val[PLANEBIGTICKS] = PLANEBIGTICKVALUES;
   long intervals = PLANESIZE;
   for (i = 0; i < intervals; i++) {
	  if (i % 7)
		 line[i] = '-';
	  else
		 line[i] = '+';
   }
   line[i] = '\0';
   sprintf(plane[0], "     +%s+    +%s+   ", line, line);
   for (i = 0; i < intervals; i++) {
	  memset(plane[i + 1], ' ', sizeof(char) * (intervals + intervals + 20));
	  plane[i + 1][intervals + intervals + 19] = '\0';
	  if (!((i) % 7)) {
		 sprintf(plane[i + 1], "  %2.0li +", val[zz++]);
		 if(val[zz-1]==0)
		   plane[i+1][3]='0';
	  }
	  else
		 plane[i + 1][5] = '|';
	  for (j = 0; j < intervals; j++) {
		 if (pl[i][j] < contours[1]) {
			if (pl[i][j] < contours[2]) {
			   if (pl[i][j] < contours[3]) {
				  plane[i + 1][j + 6] = ' ';
			   }
			   else {
				  plane[i + 1][j + 6] = '-';
			   }
			}
			else {
			   plane[i + 1][j + 6] = '+';
			}
		 }
		 else {
			if (pl[i][j] < contours[0] - EPSILON)
			   plane[i + 1][j + 6] = '*';
			else {
			   plane[i + 1][j + 6] = 'X';
			}
		 }
	  }
	  if ((i) % 7) {
		 plane[i + 1][j + 6] = '|';
		 plane[i + 1][j + 11] = '|';
	  }
	  else {
		 plane[i + 1][j + 6] = '+';
		 plane[i + 1][j + 11] = '+';
	  }
#if 0 /* this is stuff for the second population */
	  for (j = intervals + 7 + 5 /*- 1*/; j < intervals + 7 + 5 + intervals /*- 2*/; j++) {
		 z = j - 7 - 5 /*+ 1*/;
		 if (pl[i][z] < contours[5]) {
			if (pl[i][z] < contours[6]) {
			   if (pl[i][z] < contours[7]) {
				  plane[i + 1][j] = ' ';
			   }
			   else {
				  plane[i + 1][j] = '-';
			   }
			}
			else {
			   plane[i + 1][j] = '+';
			}
		 }
		 else {
			if (pl[i][z] < contours[4] - EPSILON)
			   plane[i + 1][j] = '*';
			else {
			   plane[i + 1][j] = 'X';
			}
		 }
	  }
	  if ((i) % 7) {
		 plane[i + 1][j] = '|';
		 plane[i + 1][j + 1] = '\0';
	  }
	  else {
		 plane[i + 1][j] = '+';
		 plane[i + 1][j + 1] = '\0';
	  }
#endif
   }
   sprintf(plane[intervals+1], "     +%s+    +%s+", line, line);
   sprintf(plane[intervals+2], "     -3     -2     -1      0      1      2     -3     -2     -1      0      1      2");
}

void print_mathematica(double **plane, long x, long y)
{
   long i, j;
   static long number = 1;
   FILE *mathfile;

   mathfile = fopen("mathfile","w+");

   fprintf(mathfile, "\n\nlocus%-li={{\n", number++);
   for (i = 0; i < x-1; i++) {
      fprintf(mathfile, "{");
      for (j = 0; j < y-1; j++) {
         fprintf(mathfile, "%20.20g,", plane[i][j]);
      }
      fprintf(mathfile, "%20.20g},\n", plane[i][j]);
   }
   fprintf(mathfile, "{");
   for (j = 0; j < y-1; j++) {
      fprintf(mathfile, "%20.20g,", plane[i][j]);
   }
   fprintf(mathfile, "%20.20g}},{\n", plane[i][j]);
   for (i = 0; i < x-1; i++) {
	 fprintf(mathfile, "{");
	 for (j = y; j < 2 * y-1; j++) {
		fprintf(mathfile, "%20.20g,", plane[i][j]);
	 }
	 fprintf(mathfile, "%20.20g},\n", plane[i][j]);
  }
  fprintf(mathfile, "{");
  for (j = y; j < 2 * y-1; j++) {
	 fprintf(mathfile, "%20.20g,", plane[i][j]);
  }
  fprintf(mathfile, "%20.20g}}}\n", plane[i][j]);

   fclose(mathfile);

} /* print_mathematica */


void free_plotplane(option_struct *op, data_fmt *data, char ***plane)
{
long locus, i, numloci;

numloci = getdata_numloci(op,data);

for (locus = 0; locus < numloci + 1; locus++) {
   for (i = 0; i < PLANESIZE + 3; i++) free(plane[locus][i]);
   free(plane[locus]);
}

free(plane);

} /* free_plotplane */


void print_locusplot (option_struct *op, data_fmt *data, long lowcus)
{
char ***plane;
long numloci;

numloci = getdata_numloci(op,data);

create_plotplane(op,data,&plane);
if (lowcus!= -1)
   create_locus_plot(op,data,lowcus,plane[lowcus]);
else create_locus_plot(op,data,lowcus,plane[numloci]);
plot_surface(op,data,lowcus,plane,(long)PLANESIZE);
free_plotplane(op,data,plane);

} /* print_locusplot */















