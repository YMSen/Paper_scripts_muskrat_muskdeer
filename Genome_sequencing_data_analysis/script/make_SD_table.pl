#! usr/bin/perl -w
use strict;

use Getopt::Long;

my $usage="\nDescription:\n".
    "=========================\n".
	"$0 [options]\n\n".
	"--input <str> input SD file\n".
	"--seqLen <str> the secquence length file format:id len\n".
	"--output <str>\n".
	"--help this information\n".
	"\n\nExample:\n".
	"Use default parameter:\nperl make_SD_table.pl --input /share/project002/huangzhy/svn_test/repeatmasker/count_SD_single/OUT.SD.sort.cutoverlap.0.8 --output SD.table\n".
	"Set seqLen:\nperl make_SD_table.pl --input SD.out --seqLen seqence.Len --output SD.table\n";
	my($input,$seqLen,$output,$help);
	GetOptions("input:s"=>\$input,"seqLen:s"=>\$seqLen,"output:s"=>\$output,"help:s"=>\$help);
	$seqLen||="NULL";
	die $usage if (! defined $input ||! defined $output || $help);
	sub sortGff;
	
	

my %seqLen;
my $total_seqLen;


if($seqLen ne "NULL")
{
	##read scaffold len
	open(IN,$seqLen)||die;
	while(<IN>)
	{
		chomp;
		my @cut=split /\s+/;
		$seqLen{$cut[0]}=$cut[1];
		$total_seqLen+=$cut[1];
	}
	close IN;
}
##count_SD_coverage
my %SD_len;
my %SD_len_cut_overlap;
my @len;
my @identity;
my %Len_in_scaffold;
my $total_len;
open(IN,$input)||die;
while(<IN>)
{
next if(/#/);
chomp;
my @cut=split /\s+/;
push @len,$cut[7];
push @len,$cut[8];
push @identity,$cut[-1];
push @{$SD_len{$cut[0]}},[$cut[1],$cut[2]];
push @{$SD_len{$cut[3]}},[$cut[4],$cut[5]];
}
close IN;

=top
#draw photo
open(OUT,">$output.R.identity.source")||die;
for(my $i=0;$i<@identity;$i++)
{
	print OUT "$identity[$i]\n";
}
close OUT;
@identity=();
open(OUT,">$output.R.SDlen.source")||die;
for(my $i=0;$i<@len;$i++)
{
	 print OUT "$len[$i]\n";
}
close OUT;
@len=();

#make R.shell
open(OUT,">$output.R")||die;
print OUT "identity<-read.table("$output.R.identity.source",head=FALSE)\n"
print OUT "SDlen<-read.table("$output.R.SDlen.source",head=FALSE)\n"

=cut











##sort SD_len data
foreach my $key(keys %SD_len)
{
	sortGff(\@{$SD_len{$key}});
	#cut overlap
	my $start=${$SD_len{$key}}[0][0];
	my $end=${$SD_len{$key}}[0][1];
	if(scalar(@{$SD_len{$key}})>1)
	{
		for(my $i=1;$i<@{$SD_len{$key}};$i++)
		{
			if(${$SD_len{$key}}[$i][0]<=$end+1)
			{
				$end=$SD_len{$key}[$i][1]>=$end?$SD_len{$key}[$i][1]:$end;
			}
			else
			{
				push @{$SD_len_cut_overlap{$key}},[$start,$end];
				$start=${$SD_len{$key}}[$i][0];
				$end=${$SD_len{$key}}[$i][1];
			}
		}
		push @{$SD_len_cut_overlap{$key}},[$start,$end];
	}
	else
	{
		push @{$SD_len_cut_overlap{$key}},[$start,$end];
	}

	for(my $i=0;$i<@{$SD_len_cut_overlap{$key}};$i++)
	{
		$Len_in_scaffold{$key}+=(${$SD_len_cut_overlap{$key}}[$i][1]-${$SD_len_cut_overlap{$key}}[$i][0]+1);
		$total_len+=(${$SD_len_cut_overlap{$key}}[$i][1]-${$SD_len_cut_overlap{$key}}[$i][0]+1);
	}
}

#draw table1
if($seqLen ne "NULL")		
{
	open(OUT,">$output.1")||die;
	print OUT "#sequence_id\tsequence_len\tSD_len\t\%\n";

	#change $total_len in the reign of scaffold
	$total_len=0;
	foreach my $key(sort keys %seqLen)
	{
		if(exists $Len_in_scaffold{$key})
		{
			$total_len+=$Len_in_scaffold{$key};
		}
	}

	my $per=100*$total_len/$total_seqLen;
	print OUT "total\t$total_seqLen\t$total_len\t$per\n";
	foreach my $key(sort keys %seqLen)
	{
		if(exists $Len_in_scaffold{$key})
		{
			my $percent=100*$Len_in_scaffold{$key}/$seqLen{$key};
			print OUT "$key\t$seqLen{$key}\t$Len_in_scaffold{$key}\t$percent\n";
		}
		else 
		{
			my $percent=0;
			print OUT "$key\t$seqLen{$key}\t0\t$percent\n";
		}
	}

}
if($seqLen eq "NULL")
{
	open(OUT,">$output.1")||die;
	print OUT "#sequence_id\tSD_len\n";
	print OUT "total\t$total_len\n";
	foreach my $key(sort keys %Len_in_scaffold)
	{
		print OUT "$key\t$Len_in_scaffold{$key}\n";
	}
}




	## function
sub sortGff()
{
	my ($gff_hash)=@_;
	@{$gff_hash}=sort{$a->[0] <=> $b->[0]} @{$gff_hash};
}
					
