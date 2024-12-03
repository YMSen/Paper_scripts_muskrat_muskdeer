#!/usr/bin/perl

=head1 Name

 run_mcmctree_estimate.pl - perform divergence time estimation with PAML mcmctree.

=head1 Version

 Author: wangzhuo@genomics.org.cn;
         xiaofei@genomics.org.cn;
	 jiangwenkai@novogene.cn;

 Version 2.0, last change: 2012-06-28

=head1 Usage

 perl run_mcmctree.pl <in.phy> <in.tree> <in.desc> [options]
  <in.phy>           sequential phylip format nucleotide sequence file (*.phy)
  <in.tree>          input file contains newick format rooted tree with calibration time
  <in.desc>          input file describes scientific/common names for each leaf node in tree
  --output <dir>     output results to the directory, default ./
  --clock <num>      which molecular clock model to be used, default 3
                     1: global clock; 2: independent rates; 3: correlated rates
  --model <num>      nucleotide substitution model, default 0
                     0:JC69, 1:K80, 2:F81, 3:F84, 4:HKY85, 5:T92, 6:TN93, 7:REV(GTR), 8:UNREST 
  --rootage <num>    the age of input tree root node is less than this value,
                     default 100 (mya);
  --alpha <num>      alpha value for gamma rates at sites, default 0
  --ndata <num>      number of loci (or site partitions) in a combined analysis
  --finetune="<str>" finetune for MCMC process, default "1: .1 .1 .1 .1 .1 .1"
  --burnin <num>     the first <num> of iterations will be dicarded, default 10000
  --sampfreq <num>   the sample frequency parameter (default 2)
  --nsample <num>    the number of samples for MCMC process, default 100000\

  --sampfreq_est <num> Estimate Parameter: the sample frequency parameter (default 2)
  --nsample_est <num>  Estimate Parameter: the number of samples for MCMC process, default 100000
  --clean            delete temporary files, default not;
  --help             show this help imformation

=head1 Example

 perl run_mcmctree.pl in.phy in.tree in.desc --rootage 200

=cut

use strict;
use Getopt::Long;
use Cwd qw(abs_path);
use FindBin qw($Bin);
use lib "$Bin";

my $mcmctree = "$Bin/mcmctree";
my $mcmctree_ft = "$Bin/mcmctree_ft";
my $drawtree = "perl $Bin/draw_tree.pl";
my ($output, $clock, $model, $burnin, $nsample, $rootage, $alpha, $ndata,$sampfreq,$nsample_est,$sampfreq_est) 
	= (".", 3, 0, 10000, 100000, 100, 0, 1, 2,100000,2);
my ($help, $clean);
my $finetune = "1: .1 .1 .1 .1 .1 .1";
GetOptions(
    "output:s"  => \$output,
    "clock:i"   => \$clock,
    "model:i"   => \$model,
    "burnin:i"  => \$burnin,
    "nsample:i" => \$nsample,
    "sampfreq:i" => \$sampfreq,
    "sampfreq_est:i"=>\$sampfreq_est,
    "nsample_est:i"=>\$nsample_est,
    "rootage:f" => \$rootage,
    "alpha:f" => \$alpha,
    "ndata:i"=>\$ndata,
    "clean"    => \$clean,
    "help"    => \$help,
);

die `pod2text $0` if (@ARGV != 3 || $help);

$rootage /= 100;
$rootage = sprintf "%.1f", $rootage;
$output =~ s/\/+$//;
mkdir $output if (not -e $output);
$output = abs_path($output);

my $input     = shift;
my $time_tree = shift;
my $indesc    = shift;
$input = abs_path($input);
$time_tree = abs_path($time_tree);

# Read in tree/phy file;
open IN, "$time_tree" or die "Can't open file $time_tree: $!";
my $tmpstr = <IN>;
$/ = ";";
my $treestr = <IN>;
close IN;
$/ = "\n";
$treestr = "$tmpstr $treestr" if ($tmpstr !~ /^[\d\s]+$/);

my $newtree = $treestr;
$newtree =~ s/\'[LB\(\d\.\s\>\<\)e\-,]+\'//g;
print $newtree,"\n";
my @spnames;
push @spnames, $1 while ($newtree =~ m/[(),\s]*([\w]+)[(),\s]+/g);
my $spnum = @spnames;
system("cp $input $output/tmp.phy");

open OUT, ">$output/tmp.tree" or die "Can't write to file $output/tmp.tree: $!\n";
print OUT "   $spnum   1\n$treestr\n";
close OUT;

# Decide finetune parameters;
my $cwd = Cwd::getcwd();
chdir $output;

# Perform mcmc;
write_ctl_estimate($output, $clock, $rootage, $model, $alpha, $finetune, $burnin, $nsample_est,$sampfreq_est);
open(TO,">tmpenv.sh");
print TO "export PATH=$Bin:\$PATH";
close(TO);
system "sh tmpenv.sh";
system "$mcmctree mcmctree_estimate.ctl";
system "cp out.BV in.BV";
write_ctl($output, $clock, $rootage, $model, $alpha, $finetune, $burnin, $nsample,$sampfreq);
system "$mcmctree mcmctree.ctl";
chdir $cwd;

# Extract tree information & output to file;
open(IN, "$output/mcmctree.out") or die "Can't open file $output/mcmctree.out: $!\n";
while (<IN>) {
	if (/^Species/) {
		open (OUT, ">$output/divtree.newick") 
			or die "Can't write to file $output/divtree.newick: $!\n";
		my $line1 = <IN>;
		my $line2 = <IN>;
		$line2 = <IN>;
		my ($a, $b);
		$line2 =~ s/:\s+([\d\.]+)/$a=$1*100; ":$a";/ge;
		print OUT "$line2";
		my $line3 = <IN>;
		$line3 = <IN>;
		$line3 =~ s/:\s+([\d\.]+)/$a=$1*100; ":$a";/ge;
		$line3 =~ s/([\d\.]+)-([\d\.]+)/$a=$1*100; $b=$2*100; "\'$a-$b\'";/ge;
		print OUT "$line3";
		close OUT;
		last;
	}
}
close IN;

my $num = `wc -l $output/mcmctree.out`;
print "Original result:\n $output/mcmctree.out\n" if ($num > 0);
my $num = `wc -l $output/divtree.newick`;
print "Divergence-time tree(Mya):\n $output/divtree.newick\n" if ($num > 0);

`$drawtree $output/divtree.newick $indesc -cali $time_tree > $output/divtree.svg`;

`rm SeedUsed mcmc.out $output/tmp.tree` if ($clean);
`rm $output/mcmctree.ctl` if ($clean);


sub write_ctl_estimate {
	my ($output, $clock, $rootage, $model, $alpha, $finetune, $burnin, $nsample, $sampfreq) = @_;

	#my $sampfreq = ;
	my $ctl = <<END;
         seed =  -1
      seqfile =  $output/tmp.phy
     treefile =  $output/tmp.tree
      outfile =  $output/mcmctree.out
        ndata =  $ndata
      usedata =  3    * 0: no data; 1:seq like; 2:use in.BV; 3: out.BV
        clock =  $clock    * 1: global clock; 2: independent rates; 3: correlated rates
      RootAge =  <$rootage  * safe constraint on root age, used if no fossil for root.
        model =  $model    * 0:JC69, 1:K80, 2:F81, 3:F84, 4:HKY85
        alpha =  $alpha    * alpha for gamma rates at sites
        ncatG =  4    * No. categories in discrete gamma
    cleandata =  0    * remove sites with ambiguity data (1:yes, 0:no)?
      BDparas =  1 1 0    * birth, death, sampling
  kappa_gamma =  6 2      * gamma prior for kappa
  alpha_gamma =  1 1      * gamma prior for alpha
  rgene_gamma =  2 2   * gamma prior for overall rates for genes
 sigma2_gamma =  1 10    * gamma prior for sigma^2     (for clock=2 or 3)
     finetune =  1: .1 .1 .1 .1 .1 .1  * times, rates, mixing, paras, RateParas
        print =  1
       burnin =  $burnin
     sampfreq =  $sampfreq 
      nsample =  $nsample

END

	open (CTL, ">$output/mcmctree_estimate.ctl") or die "Can't write to file $output/mcmctreei_estimate.ctl: $!\n";
	print CTL $ctl;
	close CTL;
}

sub write_ctl {
        my ($output, $clock, $rootage, $model, $alpha, $finetune, $burnin, $nsample,$sampfreq) = @_;

        #my $sampfreq = 2;
        my $ctl = <<END;
         seed =  -1
      seqfile =  $output/tmp.phy
     treefile =  $output/tmp.tree
      outfile =  $output/mcmctree.out
        ndata =  $ndata
      usedata =  2    * 0: no data; 1:seq like; 2:use in.BV; 3: out.BV
        clock =  $clock    * 1: global clock; 2: independent rates; 3: correlated rates
      RootAge =  <$rootage  * safe constraint on root age, used if no fossil for root.
        model =  $model    * 0:JC69, 1:K80, 2:F81, 3:F84, 4:HKY85
        alpha =  $alpha    * alpha for gamma rates at sites
        ncatG =  4    * No. categories in discrete gamma
    cleandata =  0    * remove sites with ambiguity data (1:yes, 0:no)?
      BDparas =  1 1 0    * birth, death, sampling
  kappa_gamma =  6 2      * gamma prior for kappa
  alpha_gamma =  1 1      * gamma prior for alpha
  rgene_gamma =  2 2   * gamma prior for overall rates for genes
 sigma2_gamma =  1 10    * gamma prior for sigma^2     (for clock=2 or 3)
     finetune =  1: .1 .1 .1 .1 .1 .1  * times, rates, mixing, paras, RateParas
        print =  1
       burnin =  $burnin
     sampfreq =  $sampfreq
      nsample =  $nsample

END

        open (CTL, ">$output/mcmctree.ctl") or die "Can't write to file $output/mcmctree.ctl: $!\n";
        print CTL $ctl;
        close CTL;
}

