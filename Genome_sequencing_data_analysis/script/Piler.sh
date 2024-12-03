#!/bin/sh
#a small pipeline of PILER-DF

date
if [ $# -lt 1 ]
then 
	echo " Usage:$0 <genome.fa><cut><outdir> "
	exit 1
fi


script=/Data/User/tianshilin/Pipeline/pipeline_v2.0/scripts
genome=`basename $1`
cut=$2
Outdir=$3
hit=final.all.hit.gff
trs=$genome.hit.trs.gff
fams=$Outdir/$genome\_fams
cons=$Outdir/$genome\_cons
aligned_fams=$Outdir/$genome\_aligned_fams
library=$genome.piler_library.fa
log=$genome.piler.log
while [ -e $log ]
do
	rm -r $log
done
#local align
perl $script/fastaDeal_piler.pl -cutf $cut $1
perl $script/rpt_piler_pipeline.pl $Outdir/$genome.cut/
perl $script/qsub-sge.pl -maxjob $cut -interval 0 -reqsub -convert no $Outdir/find_pals.sh
find find_pals.sh.*.qsub -name '*.hit.gff' -exec cat {} \; > final.all.hit.gff

echo local-align-finished;
#trs
./piler2 -trs $Outdir/$hit -out $Outdir/$trs 1>>$Outdir/$log 2>&1
echo trs-finished;
while [ -d $fams ] 
do
 	rm -r $fams
done
until [ -d $fams ]
do
	mkdir $fams
done

#get seq
./piler2 -trs2fasta $Outdir/$trs -seq $1 -path $fams 1>>$Outdir/$log 2>&1
echo piler-finished;
while [ -d $aligned_fams ]
do
	rm -r $aligned_fams
done
until [ -d $aligned_fams ]
do
	mkdir $aligned_fams
done

cd $fams
#muscle
for i in *
do ./muscle -in $fams/$i -maxiters 1 -diags -out $aligned_fams/$i  1>>$Outdir/$log 2>&1
done 
echo muscle-finished;
cd ..

while [ -d $cons ]
do
	rm -r $cons
done
until [ -d $cons ]
do
	mkdir $cons
done
cd $aligned_fams

#get cons seq 
for i in *
do ./piler2   -cons $aligned_fams/$i -out $cons/$i -label $i 1>>$Outdir/$log 2>&1
done

cd $cons
cat * > $Outdir/$library
echo cons-finished;
cd ..

date
