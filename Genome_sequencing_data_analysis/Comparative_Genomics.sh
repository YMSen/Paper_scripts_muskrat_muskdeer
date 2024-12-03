
#1. Gene family cluster analysis
orthofinder -f <the protein sequence files from 20 focal mammals> -o <OUTDIR>

#We minimized the impact of multiple sequence alignment errors and divergent regions by applying the Gblocks (v0.91b) package

Gblocks alignment.fasta -e=.fa -b5=n -t=c -b0=10

#The high-quality alignments were obtained from single-copy orthologs with PRANK (v.170427)

prank -d=alignment.fasta.fa -o=$i -f=paml -convert

#We discarded alignments shorter than 150 nt and used the remaining 6,658 single-copy orthologs for subsequent analysis.


#2. Phylogeny construction and divergence time estimation.

#We performed maximum likelihood (ML) phylogenetic reconstruction based on our set of 6,658 loci for 20 focal mammals in RaxML (v8.2.12)
/raxmlHPC-PTHREADS -s sp20.phy -n nwk -m GTRGAMMA -p 12345 -x 12345 -# 1000 -f ad -T 30 -o human,mouse 


#To generate a time-calibrated tree, we estimated divergence times in MCMCTree in PAML 4.9
perl ./script/run_mcmctree_estimate.pl ./input.phy.4d ./tree.nwk ./in.desc --output ./mcmc --rootage 120 --clock 3

#3. Expansion and contraction analysis
We performed expansion and contraction analysis follows the methodology outlined in section 3.1.3 of the Café5 tutorial (https://github.com/hahnlab/CAFE5/blob/master/docs/tutorial/tutorial.md)

#4. Identification of genes with undergoing both positive selection and rapid evolution

# Analysis of rapid evolution using branch model using CodeML
codeml alte.ctl
codeml null.ctl


#we compared a two-ratio branch model in which the foreground branch can evolve at a faster rate than the background to model M0 in which ω is fixed across the tree. 

#The alte.ctl and null.ctl were configured as follows:

seqfile = seqfile * sequence data filename
treefile = treefile      * tree structure file name
outfile = alte.mlc         * main result file name

noisy = 9  * 0,1,2,3,9: how much rubbish on the screen
verbose = 1  * 0: concise; 1: detailed, 2: too much
runmode = 0  * 0: user tree;  1: semi-automatic;  2: automatic
* 3: StepwiseAddition; (4,5):PerturbationNNI; -2: pairwise

seqtype = 1  * 1:codons; 2:AAs; 3:codons-->AAs
CodonFreq = 2  * 0:1/61 each, 1:F1X4, 2:F3X4, 3:codon table

*ndata = 10
clock = 0  * 0:no clock, 1:clock; 2:local clock; 3:CombinedAnalysis
aaDist = 0  * 0:equal, +:geometric; -:linear, 1-6:G1974,Miyata,c,p,v,a
aaRatefile = dat/jones.dat  * only used for aa seqs with model=empirical(_F)
* dayhoff.dat, jones.dat, wag.dat, mtmam.dat, or your own

model = 2
* models for codons:
* 0:one, 1:b, 2:2 or more dN/dS ratios for branches
* models for AAs or codon-translated AAs:
* 0:poisson, 1:proportional, 2:Empirical, 3:Empirical+F
* 6:FromCodon, 7:AAClasses, 8:REVaa_0, 9:REVaa(nr=189)

NSsites = 0  * 0:one w;1:neutral;2:selection; 3:discrete;4:freqs;
* 5:gamma;6:2gamma;7:beta;8:beta&w;9:beta&gamma;
* 10:beta&gamma+1; 11:beta&normal>1; 12:0&2normal>1;
* 13:3normal>0

icode = 0  * 0:universal code; 1:mammalian mt; 2-10:see below
Mgene = 0
* codon: 0:rates, 1:separate; 2:diff pi, 3:diff kapa, 4:all diff
* AA: 0:rates, 1:separate

fix_kappa = 0  * 1: kappa fixed, 0: kappa to be estimated
kappa = 2  * initial or fixed kappa
fix_omega = 0  * 1: omega or omega_1 fixed, 0: estimate
omega = .4 * initial or fixed omega, for codons or codon-based AAs

fix_alpha = 1  * 0: estimate gamma shape parameter; 1: fix it at alpha
alpha = 0. * initial or fixed alpha, 0:infinity (constant rate)
Malpha = 0  * different alphas for genes
ncatG = 8  * # of categories in dG of NSsites models

getSE = 0  * 0: don't want them, 1: want S.E.s of estimates
RateAncestor = 1  * (0,1,2): rates (alpha>0) or ancestral states (1 or 2)

Small_Diff = .5e-6
cleandata = 1  * remove sites with ambiguity data (1:yes, 0:no)?
*  fix_blength = -1  * 0: ignore, -1: random, 1: initial, 2: fixed
method = 0  * Optimization method 0: simultaneous; 1: one branch a time




#To configure a template file (null.ctl) for running " null model"

seqfile = seqfile * sequence data filename
treefile = treefile      * tree structure file name
outfile = null.mlc         * main result file name

noisy = 9  * 0,1,2,3,9: how much rubbish on the screen
verbose = 1  * 0: concise; 1: detailed, 2: too much
runmode = 0  * 0: user tree;  1: semi-automatic;  2: automatic
                        * 3: StepwiseAddition; (4,5):PerturbationNNI; -2: pairwise

seqtype = 1  * 1:codons; 2:AAs; 3:codons-->AAs
CodonFreq = 2  * 0:1/61 each, 1:F1X4, 2:F3X4, 3:codon table

*ndata = 10
clock = 0  * 0:no clock, 1:clock; 2:local clock; 3:CombinedAnalysis
aaDist = 0  * 0:equal, +:geometric; -:linear, 1-6:G1974,Miyata,c,p,v,a
aaRatefile = dat/jones.dat  * only used for aa seqs with model=empirical(_F)
                        * dayhoff.dat, jones.dat, wag.dat, mtmam.dat, or your own

model = 0 * models for codons:
                 * 0:one, 1:b, 2:2 or more dN/dS ratios for branches
                 * models for AAs or codon-translated AAs:
                * 0:poisson, 1:proportional, 2:Empirical, 3:Empirical+F
                * 6:FromCodon, 7:AAClasses, 8:REVaa_0, 9:REVaa(nr=189)

NSsites = 0  * 0:one w;1:neutral;2:selection; 3:discrete;4:freqs;
                * 5:gamma;6:2gamma;7:beta;8:beta&w;9:beta&gamma;
                * 10:beta&gamma+1; 11:beta&normal>1; 12:0&2normal>1;
                * 13:3normal>0

icode = 0  * 0:universal code; 1:mammalian mt; 2-10:see below
Mgene = 0 * codon: 0:rates, 1:separate; 2:diff pi, 3:diff kapa, 4:all diff
                 * AA: 0:rates, 1:separate

fix_kappa = 0  * 1: kappa fixed, 0: kappa to be estimated
kappa = 2  * initial or fixed kappa
fix_omega = 0  * 1: omega or omega_1 fixed, 0: estimate
omega = .4 * initial or fixed omega, for codons or codon-based AAs

fix_alpha = 1  * 0: estimate gamma shape parameter; 1: fix it at alpha
alpha = 0. * initial or fixed alpha, 0:infinity (constant rate)
Malpha = 0  * different alphas for genes
ncatG = 8  * # of categories in dG of NSsites models

getSE = 0  * 0: don't want them, 1: want S.E.s of estimates
RateAncestor = 1  * (0,1,2): rates (alpha>0) or ancestral states (1 or 2)
method = 0  * Optimization method 0: simultaneous; 1: one branch a time
Small_Diff = .5e-6
cleandata = 1  * remove sites with ambiguity data (1:yes, 0:no)?
*  fix_blength = -1  * 0: ignore, -1: random, 1: initial, 2: fixed



# Analysis of positive selection using branch-site model using CodeML
codeml alte.ctl
codeml null.ctl

#To configure a template file (alte.ctl) for running "branch-site model A"
eqfile = seqfile
treefile = treefile
outfile = alte.mlc

noisy = 3
 verbose = 1
 runmode = 0

 seqtype = 1
CodonFreq = 2
 clock = 0
model = 2

NSsites = 2
icode = 0

fix_kappa = 0
kappa = 2.5
fix_omega = 0
 omega = 2

 fix_alpha = 1
 alpha = .0
 Malpha = 0
 ncatG = 4

 getSE = 0
 RateAncestor =0

 fix_blength = -1
 method = 0



#To configure a template file (null.ctl) for running " null model"

seqfile = seqfile
treefile = treefile
outfile = null.mlc

noisy = 3
 verbose = 1
 runmode = 0

 seqtype = 1
CodonFreq = 2
 clock = 0
model = 2

NSsites = 2
icode = 0

fix_kappa = 0
kappa = 2.5
fix_omega = 1
 omega = 1

 fix_alpha = 1
 alpha = .0
 Malpha = 0
 ncatG = 4

 getSE = 0
 RateAncestor = 0

 fix_blength = -1
 method = 0
