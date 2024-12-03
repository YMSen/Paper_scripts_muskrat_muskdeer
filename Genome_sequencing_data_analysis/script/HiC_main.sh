#!/bin/bash
#########################################################################
# File Name: run.sh
#########################################################################
set -vex
set -o pipefail

[ -n "$1" ] && source $1

export PAHT=./script/:$PATH
cutfile=split_reads
genomefile=genome
hicupfile=hicup
resultfile=result
allbam=bam
rlreads='left.fq right.fq'
enzstr=`grep -i ^$enz'\s' ./script/bionet.902 | cut -f 2 | awk '{print $NF","$1}'`
enzm=`grep -i ^$enz'\s' ./script/bionet.902 | cut -f 2 | awk '{print $1}'`


if [ "$filter" == "yes" ];then
    hicupBin=./HiCUP-0.8.0/hicup_ray_raw
else
    hicupBin=./HiCUP-0.8.0/hicup_ray
fi

if [ "$fill" == "yes" ];then
    nofill=""
else
    nofill="--nofill"
fi

mkdir -p $workdir

cd $workdir
mkdir -p $genomefile $cutfile $hicupfile $resultfile $allbam

{
# bowtie-build
if [ ! -f $workdir/$genomefile/build.done ];then
	cd $workdir/$genomefile
	./bowtie2-2.3.5/bin/bowtie2-build $genome ref && touch build.done
fi
} & 

{
# digester
if [ ! -f $workdir/$genomefile/digester.done ];then
	cd $workdir/$genomefile
	./HiCUP-0.8.0/hicup_digester  --re1 $enzstr $genome && touch digester.done
fi
} & 


# fqSplit
if [ ! -f $workdir/$cutfile/split.done ];then
	(
	cd $workdir/$cutfile
    cat <<EOF | awk -f /dev/stdin $readsfile > split.cmd
{
    if(\$1~/gz$/){
        cat="zcat";
    }else{
        cat="cat";
    }
    print cat, \$1" | LANG=en_US.UTF-8 ./script/split -l 13000000 -d --verbose --additional-suffix _1 - $workdir/$cutfile/r_"NR"_ > $workdir/$cutfile/log_"NR"_1";
    if(\$2~/gz$/){
        cat="zcat";
    }else{
        cat="cat";
    }
    print cat, \$2" | LANG=en_US.UTF-8 ./scipt/split -l 13000000 -d --verbose --additional-suffix _2 - $workdir/$cutfile/r_"NR"_ > $workdir/$cutfile/log_"NR"_2";
}
EOF
    cat split.cmd | ./script/parallel -j 10 
    sleep 30

    cat split.cmd | awk '{print $NF}' | xargs cat | awk '{print $3}' | perl -pe "s/[‘’\']//g" | sort | paste - - | perl -lane 'chomp; $F[0]=~s#.*/##; $F[0]=~s/_1$//; print "$_\t$F[0]"' > split.log  && [ -s split.log ]
    ) && touch $workdir/$cutfile/split.done
fi

wait

list='perl -e '\''for(<'$workdir/$genomefile'/*>){if(/Digest_unspecified_genome_'${enzm}'_None_\d\d\-\d\d\-\d\d_\d\d\-\d\d\-\d\d\d\d\.txt$/){print "$_\n"}}'\'
enzlist=`eval $list`

if [ ! -f $workdir/$hicupfile/hicup.done ];then
    (
	cd $workdir/$hicupfile
    mkdir -p log

    cat  <<EOF | awk -f /dev/stdin $workdir/$cutfile/split.log | tee run_hicup.sh | ./script/parallel -j 100 --pipe -N 5 qsub -V -cwd -l vf=10g,p=6 -sync  y -o log/ -e log/
{
    print "set -vex"; 
    print "mkdir -p $workdir/$hicupfile/"\$3" && cd $workdir/$hicupfile/"\$3; 
    print "ln -sf "\$1" left.fq && ln -sf "\$2" right.fq";
    print "export PATH=./HiCUP-0.8.0/:./script/:\$PATH";
	print "$hicupBin --zip --bowtie2 ./bowtie2-2.3.5/bin/bowtie2 --format Sanger $nofill --longest $longest_ --shortest $shortest_ --thread 10 --keep --index $workdir/$genomefile/ref --digest $enzlist left.fq right.fq"
}
EOF
    [ -s run_hicup.sh ]
    sleep 30
    ) && touch $workdir/$hicupfile/hicup.done
fi

if [ ! -f $workdir/$allbam/catbam.done ];then
	(
	cd $workdir/$allbam
    awk '{print "ls '$workdir/$hicupfile/'"$3"/*.hicup.bam"}' $workdir/$cutfile/split.log | bash > bam.list && [ -s bam.list ]
	./samtools-1.5/bin/samtools view -H `head -n1 bam.list` | head -n -4 > bam.header
	echo "./samtools-1.5/bin/samtools cat -h bam.header -o all.bam -b bam.list " | tee _catBam.sh | qsub -V -cwd -l vf=8g,p=0 -sync  y 
    sleep 30
    ) && touch $workdir/$allbam/catbam.done
fi

if [ ! -f $workdir/$resultfile/result.done ];then		
	(
	cd $workdir/$resultfile
    awk '{print "'$workdir/$hicupfile/'"$3}' $workdir/$cutfile/split.log  > alldir.txt 
	perl ./script/mergeResult.pl alldir.txt
    cd $workdir/$resultfile/merged_report
	./HiCUP-0.8.0/hicup_report --longest $longest_ --shortest $shortest_
    ) && touch $workdir/$resultfile/result.done
fi	

touch run.sh.done

