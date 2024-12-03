#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use lib './script/perllib';
use pply;
use File::Basename qw(basename dirname);
use Pod::Usage;

=head1 SYNOPSIS

	--round <int>       the times to run PASA,we usually run round1 after Trinity,run round2 after EVM.
	--trinity_seq       input the sequence after trinity. Please use the absolutly directory. 
	--genome_seq        input the genome sequence please use the absolutly directory.
	--sql_name          set the SQL name,like species name please use the absolutly directory.
	--evm               input the EVM result,use for the round-2,please use the absolutly directory.
	--strand-specific   use strand-specific para to run PASA, if the trinity seq is strand-specific you'd better use this param.
	--mysqllib          set the mysql libdir for pasa1 (eg:/path/pasa_round1/mysql_bin/data/).(required)
	--outdir            set the output directory.default is ./
	--cpu               set the cpu for pasa round 1. default is 10
	--aligner           blat(default) / blat_ly(for > 4G(2^32) Genome) / gmap(not tested!), default blat
	                    blat_ly: cpu( < 5), resorce(vf=1g,p=1). (cut Genome into 20)
                        

    version(2015-03-04)

=cut

my($round,$tseq,$gseq,$name,$outdir,$evm,$ss,$tname,$mysqllib,$cpu,$aligner);

$aligner = 'blat';
$cpu = 10;
GetOptions(
       "round:i"=>\$round,
       "cpu:i"=>\$cpu,
       "trinity_seq:s"=>\$tseq,
       "genome_seq:s"=>\$gseq,
       "sql_name:s"=>\$name,
       "evm:s"=>\$evm,
       "mysqllib:s"=>\$mysqllib,
       "outdir:s"=>\$outdir,
       "aligner:s"=>\$aligner,
       "strand-specific" =>\$ss,
);

if(!defined $round || ! defined $name || !defined $tseq){
	pod2usage(VERBOSE=>1);
}

if(( ! defined $mysqllib) && $round==2){
	print "Mysql Lib needed for round2\n";
	pod2usage(VERBOSE=>1);
}

if(defined $tseq){
	$tname=basename($tseq);
	$tseq =abspath $tseq;
}

if(defined $gseq){
	$gseq =abspath $gseq;
}

if($round == 2){
    $mysqllib = abspath $mysqllib;
    $evm = abspath $evm;
}

$outdir ||=".";
$outdir = abspath $outdir;
`mkdir -m 755 -p $outdir` unless(-d "$outdir");
 
my $param_ss = '';
$param_ss = "--transcribed_is_aligned_orient" if(defined $ss);
my $cmd;
if($round==1){
	my $seqclean_cpu = $cpu > 10 ? 10 : $cpu;
	$cmd = "
echo start at time `date +%F'  '%H:%M`
export PATH=$outdir/mysql_bin/bin/:./blast-2.2.26/bin/:\$PATH
cd $outdir
ln -s $tseq $tname
echo -e 'MYSQLDB=$name\nvalidate_alignments_in_db.dbi:--MIN_PERCENT_ALIGNED=75\nvalidate_alignments_in_db.dbi:--NUM_BP_PERFECT_SPLICE_BOUNDARY=0' > $outdir/alignAssembly.config
./script/seqclean $tname -v ./data/UniVec -c $seqclean_cpu
./script/prepare_pasa.sh $outdir;
$outdir/pasa_bin/scripts/drop_mysql_db_if_exists.dbi -c alignAssembly.config
$outdir/pasa_bin/scripts/Launch_PASA_pipeline_ly.pl -c alignAssembly.config -C -R -g $gseq -t $tname.clean -T -u $tname --ALIGNERS $aligner --CPU $cpu -N 2 $param_ss
$outdir/pasa_bin/scripts/pasa_asmbls_to_training_set.dbi --pasa_transcripts_fasta $name.assemblies.fasta --pasa_transcripts_gff3 $name.pasa_assemblies.gff3
echo finish at time `date +%F'  '%H:%M`
";

}elsif($round==2){
	$cmd = "
echo start at time `date +%F'  '%H:%M`
export PATH=$outdir/mysql_bin/bin/:\$PATH
cd $outdir
echo -e 'MYSQLDB=$name\nvalidate_alignments_in_db.dbi:--MIN_PERCENT_ALIGNED=75\nvalidate_alignments_in_db.dbi:--NUM_BP_PERFECT_SPLICE_BOUNDARY=0' > $outdir/alignAssembly.config
echo 'MYSQLDB=$name' > $outdir/annotCompare.config
./script/pasa_cds.pl $evm >$outdir/evm.cds
./script/pasa_gff3_validator.pl $outdir/evm.cds
./script/prepare_pasa.sh $outdir
cp -R $mysqllib/$name $outdir/mysql_bin/data/$name;
$outdir/pasa_bin/scripts/build_comprehensive_transcriptome.dbi -c alignAssembly.config -t $tseq
$outdir/pasa_bin/scripts/Launch_PASA_pipeline_ly.pl -c annotCompare.config -g $gseq -t $tseq -A -L --annots_gff3 $outdir/evm.cds $param_ss
./script/pasa2_longest.pl *.gene_structures_post_PASA_updates.*.gff3;
echo finish at time `date +%F'  '%H:%M`
";

}

open SH,">$outdir/PASA_round\_$round.sh";
print SH $cmd;
close SH;

system "sh $outdir/PASA_round\_$round.sh";

# check
if($round == 1){
	&checkfile("$outdir/genome.assemblies.fasta.transdecoder.genome.gff3");
}elsif($round == 2){
	&checkfile("$outdir/pasa2.longest.gff");
}


sub checkfile{
	my ($nfile) = @_;
	die "$nfile is not exist or is null\n" unless (-e $nfile && -s $nfile > 0) ;
}

