RECOAL - coalescent simulator using reference haplotype data

1. Install.

After extraction, recoal is compiled by typing "make"

2. The input files

recoal requires 3 files (+2 additional files)

 - reference haplotype file, option file and random seed file

 - fine scale recombination rate file (optional)

 - SNP calling probability file (optional)

The parameters and simulation options are read from "simoption" file.

At the "simoption" file.

// number of reference haploytes
// size of the region
// number of new simulated haplotypes
// theta (4Neu per site - NOT PER REGION)
// growth rate (Nea) - Don't give huge minus growth, samples might not coalesce in finite time
// rho (4Ner per site) - If you want to use fine scale variation of recombination rate put "-1" here
// length of burn-in for MCMC - if you give larger size, high recombination rate or big reference haplotypes, put larger number here.
// number of simulation replicates 
// intervals between replicates
// SNP calling status
// (optional) SNP calling threshold

If you want to use recombination hotspot structure, prepare "recrates" file.

For each line of recrates represents the recombination rate at a single site. (number of rec rates should be same as sequence length)

If the reference haplotype only has common SNPs, put "1" or "-1" for SNP calling status.

If SNPs were called with a certain threshold (number of variant), put "1" and give the value of threshold at the next line.

If there's information for SNP calling probability (calling probability in the function of the number of variants), put "-1" and use "probs" file.

Remember, recoal would be much slower with this hidden SNP model.  


At the "refhap" file

Each line represents the dna sequence (not SNP) of a single haplotype. 


3. Running

recoal spends most of the running time to reach the stationary distribution. As a typical MCMC searching, there is no exact answer for when to stop.

It provides you the posterior probabilities - you can guess whether it converges to the stationary.

The convergence time increases with sample size (number of samples) and number of recombination events.
   

4. Output

After searching trees for given haplotypes, a genealogy for new haplotypes are simulated based on the proposed tree. Then the new haplotypes are simulated.

The output will be written at "simhap" file.


5. Problems

Currently, simulation for larger region ( > 1Mb) or from the larger reference sample size (>100) with high recombination rates is not recommended.

- For those kinds of data, the acceptance ratio of the MCMC is very low.

- The differences between new ARG and old ARG become smaller - you need longer MCMC iteration and burn in time.

- With above problems, MCMC will never reach the stationary distribution.



6. Further development

 - Simulation for the larger regions 

 - Simulation with larger reference size.

 - Simulation with other evolutionary forces - migration, gene conversion and population structures.

 - Simulation with reference genotypes rather than haplotype data



7. Note

RECOAL is based on RECOMBINE-H which is MCMC sampler for recombination hotspot.



8. Thanks to

 - Paul Marjoram
 
 - Joe Felsenstein

 - Mary Kuhner

 - John Yamato

