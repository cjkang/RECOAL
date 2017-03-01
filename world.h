#ifndef WORLD_INCLUDE
#define WORLD_INCLUDE
/*------------------------------------------------------
 Maximum likelihood estimation 
 of migration rate  and effectice population size
 using a Metropolis-Hastings Monte Carlo algorithm                            
-------------------------------------------------------                        
 W O R L D   R O U T I N E S 

 creates tree structures,
 calculates smple parameter estimates (FST,...)
 reads tree [has to be done]
 
 prints results,
 and finally helps to destroy itself.
                                                                                                               
 Peter Beerli 1996, Seattle
 beerli@genetics.washington.edu
 $Id$
-------------------------------------------------------*/

#ifndef RECOMBINE_INCLUDE
#include "recombine.h"
#endif

void plot_surface(option_struct *op, data_fmt *data, long lowcus, char ***plane, long x);
void create_locus_plot(option_struct *op, data_fmt *data, long lowcus, char **plane);
void create_plotplane(option_struct *op, data_fmt *data, char**** plane);
void calc_locus_plane(option_struct *op, data_fmt *data, long lowcus, double **pl,
   double contours[]);
void fill_plotplane(char **plane, double **pl, double contours[]);
void print_mathematica(double **plane, long x, long y);
void print_locusplot(option_struct *op, data_fmt *data, long lowcus);
void free_plotplane(option_struct *op, data_fmt *data, char ***plane);

#endif /*WORLD_INCLUDE*/
