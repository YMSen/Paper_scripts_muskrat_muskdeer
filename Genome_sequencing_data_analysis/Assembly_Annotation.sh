##1. Genome assembly

#Step 1 the initial contigs were assembled using NextDenovo v2.5.0 (https://github.com/Nextomics/NextDenovo) and corrected with NextPolish v1.4.0

#(1) Estimating the recommended minimum seed length.
nextdenovo-2.5.0/bin/seq_stat -g 3000M -d 45 input.fofn
#(2) To configure a template file (run.cfg) for running NextDenovo

[General]
job_type = sge
job_prefix = nextDenovo
task = all  # 'all', 'correct', 'assemble'
rewrite = yes # yes/no
deltmp = yes
rerun = 3
parallel_jobs = 30
input_type = raw
read_type = ont
input_fofn = input.fofn
workdir = ./
cluster_options = -cwd -V -l vf=20g -q all.q -w n#for sge

[correct_option]
read_cutoff = 1k
seed_cutoff = <Estimated seed length>
blocksize = 1g
pa_correction = 3
seed_cutfiles = 3
sort_options = -m 20g -t 10 -k 40
minimap2_options_raw = -x ava-ont -t 20
correction_options = -p 20
genome_size = 3g
seed_depth = 45

[assemble_option]
random_round = 20
minimap2_options_cns = -x ava-ont -t 20 -k17 -w17
nextgraph_options = -a 1 -q 10

#(3) Running NextDenovo

nextdenovo-2.5.0/bin/nextDenovo run.cfg

#(4) Running Nextpolish

#To configure a template file (nextpolish.cfg)

[General]
job_type = sge
job_prefix = nextPolish
task = default
rewrite = yes
rerun = 3
parallel_jobs = 100
multithread_jobs = 8
genome = ./final.fa
genome_size = auto
workdir = ./01_rundir
cluster_options = -l vf=15G,p=8  -w n

[lgs_option]
lgs_fofn = ./lgs.fofn
lgs_options = -min_read_len 1k -max_read_len 100k -max_depth 100
lgs_minimap2_options = -x map-ont

[sgs_option]
sgs_fofn = ./sgs.fofn
sgs_options = -max_depth 200


#Running Nextpolish

nextPolish nextpolish.cfg


#Step 2: Cluster initial contigs by Hi-C interaction pair. Linked contigs were then clustered by calculating Hi-C interaction frequencies using an agglomerative hierarchical clustering algorithm implemented in ALLHiC software (v0.9.13)
#(1) To configure a template file (hic.cfg) for mapping

genome=genome.nextpolish.fasta
readsfile=read.lst  ##left.fq	right.fq
enz=DpnII
workdir=./hicup
fill=yes
format_=Sanger
longest_=800
shortest_=100
filter=yes


#(2) Running HiC mapping and generating HiC bam file(group.clean.bam)

./script/HiC_main.sh hic.cfg

#(3) Calculating Hi-C interaction frequencies
./ALLHiC-v0.9.13/allhic extract group.clean.bam genome.nextpolish.fasta --RE GATC

#(4) Cluster initial contigs 
perl ./script/partition.pl --pairsfile group.clean.pairs.txt --contigfile group.clean.counts_GATC.txt -K <chromosome number> --minREs 20 --maxlinkdensity 3 --NonInformativeRabio 0 --lencut 1.00 

#(5) Used Juicebox v1.22 (https://github.com/aidenlab/Juicebox) to visualise and adjust the placement and orientation of contigs. It will generate a cluster.fasta file.

fa=genome.nextpolish.fasta
bam=group.clean.bam
tourdir=./
touch merged_abnorm.sam merged_unmapped.sam merged_norm.txt
./samtools-1.5/samtools view $bam | awk -v "fname1"=merged_norm.txt -v "fname2"=merged_abnorm.sam -v "fname3"=merged_unmapped.sam -f ${scripts}//chimeric_blacklist.awk - 
awk '{printf("%s %s %s %d %s %s %s %d", $1, $2, $3, 0, $4, $5, $6, 1); for (i=7; i<=NF; i++) {printf(" %s",$i);}printf("\n");}' merged_norm.txt > merged.frag.txt
sort -T ./ -k2,2d -k6,6d -k4,4n -k8,8n -k1,1n -k5,5n -k3,3n merged.frag.txt > merged.sort.txt
touch dups.txt optdups.txt merged_nodups.txt
awk -f ${scripts}/dups.awk -v name=./ merged.sort.txt
python ./script/makecprops.py $fa".fai" > $fa".cprops"
perl ./scripts/tour2asm.pl $fa".cprops" $tourdir
bash ./3d-dna-master/visualize/run-asm-visualizer.sh -p true $fa".cprops" genome.asm merged_nodups.txt


#Step 3: Classify the Nanopore reads to perform local assembly
#(1) Mapping the ONT reads to cluster.fasta
minimap2-2.17-r941/bin/minimap2 -x map-pb -a -t 40 cluster.fasta ont.fa -o ont.sam
./samtools-1.5/samtools view -bS ont.sam -@ 20 -o ont.bam
./samtools-1.5/samtools sort -@ 40 -m 2G ont.bam -o ont.sort.bam 

#(2) Select Mapped ONT reads for each cluster contigs.

#For each cluster:
samtools view -@ 2 -h -bF 4 ont.sort.bam -L cluster[number].bed >cluster_[number].bam
samtools fasta cluster_[number].bam >cluster_[number].fasta

#(3) Performing Local assemebly

Repeat the assembly of Step 1 for each cluster.

#step 4: Anchor contigs onto chromosomes

Repeat the assembly of Step 2.




##2. Repetitive sequences identification

#(1) Homology-based prediction
# RepeatMasker (v4.0.5)
RepeatMasker -a -nolow -no_is -norna -parallel 40 -lib ./data/RepeatMasker.lib -s genome.fasta >genome.fa.1.log 2>genome.fa.2.log

#RepeatProteinMask (v4.0.5)
RepeatProteinMask -noLowSimple -pvalue 0.0001 -engine ncbi genome.fasta

#(2) Ab initio prediction

#Using RepeatModeler (v1.0.8) to generate a de novo repeat library (repeatmodeler.library)
BuildDatabase -engine ncbi -name genome genome.fasta
RepeatModeler -database genome -engine ncbi -threads 15 >&run.out

#Using LTR FINDER (v1.0.7) to generate a de novo repeat library (ltr.library)
ltr_finder -C -w 2 -s ./data/Sc-tRNAs.fa genome.fasta >./ltr.library 2>./ltr.log

#Using RepeatScout (v1.0.5) to generate a de novo repeat library (RepeatScout.library)
build_lmer_table -sequence genome.fasta -freq ./genome.freq
RepeatScout -sequence ./genome.fasta -freq ./genome.freq -output ./genome.out
perl ./scripts/filter-stage-1.prl ./genome.out >./RepeatScout.library

#Using PILER (v3.3.0) to generate a de novo repeat library (RepeatScout.library)
sh ./script/Piler.sh genome.fa 20 ./Piler.out

#Using RepeatProteinMask (v4.0.5) generate Repetitive sequences using 4 de novo repeat libraries

cat repeatmodeler.library ltr.library RepeatScout.library RepeatScout.library >denovo.lib
perl ./script/changeAa.pl denovo.lib >denovo.lib.convert
perl ./script/uclust.pl denovo.lib.convert >denovo.lib.convert.unredundance.fa
RepeatMasker -a -nolow -no_is -norna -parallel 40 -lib denovo.lib.convert.unredundance.fa -s genome.fasta >genome.fa.1.log 2>genome.fa.2.log

#(3) Using Identifying tandem repeats sequences using TRF
trf genome.fasta 2 7 7 80 10 50 2000 -d -h



##3. Gene structure

#(1) Homology-based prediction. We aligned these two non-redundant uniform protein sets (from the suborder Yangochiroptera and Yinpterochiroptera) plus proteins from horse, human and mouse, to the respective bat genome assemblies

blastall -p tblastn -e 1e-05 -F T -m 8 -d genome.fasta -i query.pep.fasta -o query.pep.blast
perl ./script/solar.pl -a prot2genome2 -z -f m8 query.pep.blast >query.pep.blast.solar
perl ./script/struct_solar_filter.pl query.pep.blast.solar 0.25
perl ./script/struct_solar_to_table.pl query.pep.blast.solar.filter 
perl ./script/struct_genomic_cluster.pl -overlap_percent 0.5 query.pep.blast.solar.filter.table > query.pep.blast.solar.filter.table.nonredundance
perl ./script/fishInWinter.pl -bf table -ff table  query.pep.blast.solar.filter.table.nonredundance query.pep.blast.solar.filter >query.pep.blast.solar.filter.nr
./wise2.4.1/genewise -trev -genesf -gff -sum each.query.gene.fa each.genome.region.fa >each.genome.region.genewise


#(2) Predict genes using RNA-seq data

#predicted the isoform structure (gene model)
./bowtie2-2.3.5/bin/bowtie2-build -q -f genome.fasta genome
./tophat-2.0.13/bin/tophat -p 6 --max-intron-length 500000 -m 2 --library-type fr-unstranded  -o ./outdir genome tissue_1.fq.gz tissue_2.fq.gz
./Cufflinks-2.1.1/cufflinks  -I 500000 -p 1 --library-type fr-unstranded -L CUFF -o ./outdir accepted_hits.bam
perl ./script/Cufflink_gtf2gff3.pl transcripts.gtf >transcripts.gff3

# Assembled transcripts from RNA-seq data using Trinity (v2.1.1)

./trinity-2.1.1/bin/Trinity --seqType fq --left tissue_1.fq.gz --right tissue_2.fq.gz --CPU 20 --max_memory 200G --normalize_reads --full_cleanup --min_glue 2 --min_kmer_cov 2 --KMER_SIZE 25 --output trinity.out


#Assembled transcripts were aligned to the assembled genome using the software PASA to generate the training set
perl ./script/pipeline_pasa.pl --genome_seq genome.fasta  --round 1 --cpu 10 --sql_name genome --trinity_seq trinity.out.fasta --outdir pasa1

perl ./script/training_gff3Togff.pl pasa1.gff3 >pasa1.gff
perl ./script/training_grep_complete.pl pasa1.pep pasa1.complete.pep
perl ./script/training_blast_database.pl --swissprot --evalue 1e-05 --cutf 20 --cpu 20 pasa1.complete.pep
./blast-2.2.26/bin/blastall -b 5  -F F  -p blastp -e 1e-05  -d uniprot_sprot.fasta -a 10 -i pasa1.complete.pep -o pasa1.swissprot.blast
perl script/blast_parser.pl ./pasa1.swissprot.blast > ./pasa1.swissprot.blast.tab;
mv pasa1.swissprot.blast.tab pasa1.complete.pep.nr.blast.tab
awk '$13>=95' pasa1.complete.pep.nr.blast.tab >pasa1.complete.pep.nr.blast.tab.score95
perl ./script/training_blast_stat.pl pasa1.complete.pep.nr.blast.tab.score95 >pasa1.complete.pep.nr.blast.tab.ratio
awk '$4>=0.8' pasa1.complete.pep.nr.blast.tab.ratio >pasa1.complete.pep.nr.blast.tab.ratio.cvg0.8
perl ./script/fastaDeal.pl -attr id:len pasa1.complete.pep >pasa1.complete.pep.len
perl ./script/training_find_nosame_with_head.pl pasa1.complete.pep.nr.blast.tab.ratio pasa1.complete.pep.len >pasa1.complete.pep.nr.blast.tab.ratio.noalign
awk '$2>=1000' pasa1.complete.pep.nr.blast.tab.ratio.noalign >pasa1.complete.pep.nr.blast.tab.ratio.noalign.longORF
cat pasa1.complete.pep.nr.blast.tab.ratio.cvg0.8 pasa1.complete.pep.nr.blast.tab.ratio.noalign.longORF >pasa1.train
perl ./script/training_get.gff.pl pasa1.train pasa1.gff >pasa1.train.gff
perl ./script/training_perfect_gene4pasa.pl --start 10 --stop 10 genome.fasta pasa1.train.gff
perl ./script/training_clustergff.pl pasa1.train.gff
perl ./script/training_trainset_uniq.pl pasa1.train.gff.nr.gff
perl ./script/training_trainset_random.pl pasa1.train.gff.nr.gff.uniq.gff 1000


#(3) Ab initio-based prediction predict using the soft-masked genome using Augustus (v2.5.5), GlimmerHMM (v3.0.1) and SNAP (Semi-HMM-based Nucleic Acid Parser)

#  Augustus (v2.5.5)
perl ./augustus/scripts/autoAugTrain.pl --genome=genome.masked.fasta --trainingset=PASA4training/pasa1.train.gff.nr.gff.uniq.gff.random.gff --optrounds=0 --species=pasa1
./augustus/bin/augustus --species=pasa1 --AUGUSTUS_CONFIG_PATH=training/Augustus/conf --extrinsicCfgFile=training/Augustus/conf/extrinsic/extrinsic.cfg --uniqueGeneId=true --noInFrameStop=true --gff3=on --genemodel=complete --strand=both  genome.masked.fasta >genome.augustus


# GlimmerHMM (v3.0.1)
perl ./script/training_auto_train_glimmerhmmm.pl pasa1.train.gff.nr.gff.uniq.gff.random.gff genome.masked.fasta
./script/glimmerhmm genome.masked.fasta -d training/GlimmerHMM/pasa1 -f -g >genome.gff

# SNAP (Semi-HMM-based Nucleic Acid Parser)
perl ./script/training_train_snap_new.pl --name genome.masked.fasta pasa1 pasa1.train.gff.nr.gff.uniq.gff.random.gff
./script/snap-2013-11-29/snap -gff training/SNAP/pasa1.hmm genome.masked.fasta >genome.snap

#GeneID (v1.4)
geneid -P ./data/homo_sapiens.param -v -G -p geneid genome.masked.fasta >genome.Geneid

#GeneScan (v1.0)
./script/genscan ./data/HumanIso.smat genome.masked.fasta >genome.genscan

#(4) EVidenceModeler (v1.1.1)  was used to integrate all predictions

#To configure a weight file (weights.txt) for running NextDenovo
PROTEIN GeneWise        100
TRANSCRIPT      CUFF    100
ABINITIO_PREDICTION     Augustus        10
ABINITIO_PREDICTION     genscan 1
ABINITIO_PREDICTION     GeneID  1
ABINITIO_PREDICTION     GlimmerHMM      1
ABINITIO_PREDICTION     SNAP    1
OTHER_PREDICTION        transdecoder    100

# Running EVidenceModeler (v1.1.1) 
We used the parameters recommended by the software to run the EVidenceModeler.
Please see https://github.com/EVidenceModeler/EVidenceModeler/wiki#running-evm
We generated the file "evm.gff3".

#(5) PASA2 was used to predict untranslated regions and alternative splicing variations

perl ./script/pipeline_pasa.pl --round 2 --trinity_seq pasa1/genome.assemblies.fasta.transdecoder.cds --genome_seq genome.fasta  --sql_name genome --mysqllib pasa1/mysql_bin/data/  --evm evm.gff3 --outdir 08.pasa 

