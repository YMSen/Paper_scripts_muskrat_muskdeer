#! usr/bin/perl -w
use strict;
use Getopt::Long;

my $usage="\nDescription:\n".
	"=========================\n".
	"$0 [options]\n\n".
	"--inputDir <str> the dir of lastz output\n".
	"--seqLen <str> the secquence length file format:id len\n".
	"--output <str>\n".
	"--sdLen <num> the minimal SD length to accept. Default: --sdLen 1000\n".
	"--identityType [1,2] 1 fot blast identity 2 for lastz identity. Default: blast identity --identity 1\n".
	"--identityCut <min-max> set the  SD identity range to accept. Default: --identityCut 0.8-1\n".
	"--innerOverlap <num> the maxnum percentage of overlap between one pair of SD to accept(for example 0.2), Defult: --innerOverlap 0\n".
	"--help this information\n".
	"\n\nExample:\n".
	"Use default parameter:\nperl count_sd_update.pl  --inputDir /share/project002/huangzhy/svn_test/repeatmasker/lastz_out --seqLen /share/project002/huangzhy/svn_test/repeatmasker/hsal_v3.scaffoldLen  --output OUT.SD\n".
	"Set identity range:\nperl count_sd_update.pl  --inputDir ../lastz_out --seqLen ../hsal_v3.scaffoldLen --identityCut 0.9-0.98 --output OUT.SD\n";
	my($inputDir,$seqLen,$output,$sdLen,$identityType,$identityCut,$innerOverlap,$help);
	GetOptions("inputDir:s"=>\$inputDir,"seqLen:s"=>\$seqLen,"output:s"=>\$output,"sdLen:s"=>\$sdLen,"identityType:i"=>\$identityType,"identityCut:s"=>\$identityCut,"innerOverlap:f"=>\$innerOverlap,"help:s"=>\$help);

$sdLen=1000 unless defined $sdLen;
$identityType ||=1;
$identityCut ||="0.8-1";
$innerOverlap=0 unless defined $innerOverlap;
die $usage if (! defined $inputDir || ! defined $seqLen ||! defined $output || $help);
	
#get identity range;
my @identityCut=split(/-/,$identityCut);
my $minIdentity=$identityCut[0];
my $maxIdentity=$identityCut[1];

#count sequence length

my %seqLen;
open(IN,$seqLen)||die;
while(<IN>)
{
	chomp;
	my @cut=split /\s+/;
	$seqLen{$cut[0]}=$cut[1];
}
close IN;

print "scaffoldLen loaded\n";
## read lastz output and count Sd


sub get_identity;
open(OUT,">$output")||die;
my $title=join("\t",("#Seq1","Start1","End1","Seq2","Start2","End2","Strand","Seq1Len","Seq2Len","BlockLen","IndelNumber","IndelSpece","AlignBase","MatchBase","MismatchBase","Score","Identity"));
print OUT "$title\n";
opendir (DIR,$inputDir)||die;
while(my $dir=readdir DIR)
{
	next if $dir !~ /\w/;
	next if $dir !~ /\.axt$/;
	my $file = $inputDir."/".$dir;
	open (FILE,$file)||die;
	while(<FILE>)
	{
		next if (/#/);
		next if ($_ eq "");
		chomp;
		if(/\+/ || /\-/)
		{
			my $message=$_;
			my $seq1=<FILE>;
			chomp $seq1;
			my $seq2=<FILE>;
			chomp $seq2;
			#count identity
			get_identity($message,$seq1,$seq2);
		}

	}
	close FILE;
}


print "original SD data complete\n";


##sort SD data

my %sort;
open(IN,$output)||die;
open(OUT,">$output.sort")||die;
my $head=<IN>;
print OUT $head;
while(<IN>)
{
	
	chomp;
	my @cut=split /\s+/;
	$sort{$cut[0]}{$cut[3]}{$cut[1]}{$cut[4]}{$cut[2]}{$cut[5]}{$cut[6]}=$_;
}
close IN;
foreach my $key1(sort keys %sort)
{
	foreach my $key2(sort keys %{$sort{$key1}})
	{
		foreach my $key3(sort {$a<=>$b} keys%{$sort{$key1}{$key2}})
		{
			foreach my $key4(sort{$a<=>$b} keys %{$sort{$key1}{$key2}{$key3}})
			{
				foreach my $key5(sort{$a<=>$b} keys %{$sort{$key1}{$key2}{$key3}{$key4}})
				{
					foreach my $key6(sort{$a<=>$b} keys %{$sort{$key1}{$key2}{$key3}{$key4}{$key5}})
					{
						foreach my $key7(keys %{$sort{$key1}{$key2}{$key3}{$key4}{$key5}{$key6}})
						{
							print OUT "$sort{$key1}{$key2}{$key3}{$key4}{$key5}{$key6}{$key7}\n";
						}
					}
				}
			}
		}
	}
}
close OUT;
print "SD data sorted\n";



##grep SD by extenalOverlap
open(IN,"$output.sort")||die;
open(OUT,">$output.sort.merge")||die;
my %SD;
print "#Seq1\tStart1\tEnd1\tSeq2\tStart2\tEnd2\tStrand\tSeq1Len\tSeq2Len\tidentity\n";
<IN>;
$head=<IN>;
my @head=split(/\s+/,$head);
my $seq1=$head[0];
my $start1=$head[1];
my $end1=$head[2];
my $seq2=$head[3];
my $start2=$head[4];
my $end2=$head[5];
my $stand=$head[6];
my $identity=$head[-1];
my $num=1;
while(<IN>)
{
	 chomp;
	 my @cut=split /\s+/;
	 my $cov1_start=$start1>$cut[1]?$start1:$cut[1];
	 my $cov1_end=$end1<$cut[2]?$end1:$cut[2];
	 my $cov2_start=$start2>$cut[4]?$start2:$cut[4];
	 my $cov2_end=$end2<$cut[5]?$end2:$cut[5];
	 if($cut[0] ne $seq1 || $cut[3] ne $seq2 ||$cut[6] ne $stand || $cov1_start>$cov1_end || $cov2_start>$cov2_end)
	 {
		 my $seqlen1=$end1-$start1+1;
		 my $seqlen2=$end2-$start2+1;
		 $identity=$identity/$num;
		 print OUT "$seq1\t$start1\t$end1\t$seq2\t$start2\t$end2\t$stand\t$seqlen1\t$seqlen2\t$identity\n";
		 $seq1=$cut[0];
		 $start1=$cut[1];
		 $end1=$cut[2];
		 $seq2=$cut[3];
		 $start2=$cut[4];
		 $end2=$cut[5];
		 $stand=$cut[6];
		 $identity=$cut[-1];
		 $num=1;
	 }
	 else
	 {
		 $start1=$start1<$cut[1]?$start1:$cut[1];
		 $end1=$end1>$cut[2]?$end1:$cut[2];
		 $start2=$start2<$cut[4]?$start2:$cut[4];
		 $end2=$end2>$cut[5]?$end2:$cut[5];
		 $identity+=$cut[-1];
		 $num++;
	 }
}
close IN;
my $seqlen1=$end1-$start1+1;
my $seqlen2=$end2-$start2+1;
$identity=$identity/$num;
print OUT "$seq1\t$start1\t$end1\t$seq2\t$start2\t$end2\t$stand\t$seqlen1\t$seqlen2\t$identity\n";
close IN;
close OUT;



sub get_identity()
{
	my ($info,$sequence1,$sequence2)=@_;
	#my @cut=split(/ /,$info);
	my @cut=split(/\s+/,$info);
	my $scaffold1=$cut[1];
	my $start1=$cut[2];
	my $end1=$cut[3];
	my $scaffold2=$cut[4];
	my $start2=$cut[5];
	my $end2=$cut[6];
	my $strand=$cut[7];
	my $score=$cut[8];
	my $seqLen1=$end1-$start1+1;
	my $seqLen2=$end2-$start2+1;
	my $indelN=0;                 #indel number
	my $indelS=0;                 #indel spece
	my $blockLen=0;				  #total sequence length
	my $alignB=0;                 #bases align
	my $matchB=0;                 #aligned bases that match
	my $mismatchB=0;              #aligned bases that do not match 
	my $identity=0;

	
	## first filter Sd length
	if($seqLen1<$sdLen || $seqLen2 < $sdLen)
	{
		##do notthing 
	}
	
	else
	{
			
			for(my $i=0;$i<length($sequence1);$i++)
			{
				

				if(substr($sequence1,$i,1) ne "-" && substr($sequence2,$i,1) ne "-" &&  uc(substr($sequence1,$i,1)) ne "N" && uc(substr($sequence2,$i,1)) ne "N"  && uc(substr($sequence1,$i,1)) eq uc(substr($sequence2,$i,1)))
				{
					$matchB++;
					$blockLen++;
					$alignB++;
				}
				elsif(substr($sequence1,$i,1) ne "-" && substr($sequence2,$i,1) ne "-" &&  uc(substr($sequence1,$i,1)) ne "N" && uc(substr($sequence2,$i,1)) ne "N"   && uc(substr($sequence1,$i,1)) ne uc(substr($sequence2,$i,1)))
				{
					$mismatchB++;
					$blockLen++;
					$alignB++;
				}
				elsif(substr($sequence1,$i,1) eq "-" || substr($sequence2,$i,1) eq "-" || uc(substr($sequence1,$i,1)) eq "N" || uc(substr($sequence2,$i,1)) eq "N" )
				{
					$indelS++;
					$blockLen++;
				}
			}
			my $indelN1=$sequence1=~s/-{1,}//g;
			my $indelN2=$sequence2=~s/-{1,}//g;
			$indelN=$indelN1+$indelN2;
			if($strand eq "-")
			{
				my $b=$start2;
				$start2=$seqLen{$scaffold2}-$end2+1;
				$end2=$seqLen{$scaffold2}-$b+1;
			}
			if($identityType==1)
			{
				$identity=$matchB/$blockLen;
			}
			if($identityType==2)
			{
				$identity=$matchB/$alignB;
			}			

			#sort the SD result
			if(($scaffold1 eq $scaffold2 && $start1>$start2)||($scaffold1 gt $scaffold2))
			{	
				#change scaffold1 and scaffold2
				my $scaffold=$scaffold1;
				my $start=$start1;
				my $end=$end1;
				my $seqLen=$seqLen1;
				$scaffold1=$scaffold2;
				$start1=$start2;
				$end1=$end2;
				$seqLen1=$seqLen2;
				$scaffold2=$scaffold;
				$start2=$start;
				$end2=$end;
				$seqLen2=$seqLen;
			}
			
			#cut innr overlap
		if($scaffold1 eq $scaffold2 && ($end1-$start2+1> $innerOverlap*$seqLen1 || $end1-$start2+1> $innerOverlap*$seqLen2))
		{
			
		}
			
		elsif($identity<$minIdentity||$identity>$maxIdentity)
		{
		
		}
		else
		{	
			my $join=join("\t",($scaffold1,$start1,$end1,$scaffold2,$start2,$end2,$strand,$seqLen1,$seqLen2,$blockLen,$indelN,$indelS,$alignB,$matchB,$mismatchB,$score,$identity));
			print OUT "$join\n";
		}

	}
}









			#my $title=$join("\t",("#Seq1","Start1","End1","Seq2","Start2","End2","Strand","Seq1Len","Seq2Len","BlockLen","IndelNumber","IndelSpece","AlignBase","MatchBase","MismatchBase","Score","Identity"));
		










