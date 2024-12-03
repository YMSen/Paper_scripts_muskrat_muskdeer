#!/usr/bin/perl

use strict;

die"Usage:\n\tperl $0 <gff> > <out.file>\n" if(@ARGV != 1);

my $file=shift;
my %gene;

read_gff($file,\%gene);
my $out;
for my $name(keys %gene){
	for my $id(keys %{$gene{$name}}){
		@{$gene{$name}{$id}{cds}}=sort {$a->[1] <=> $b->[1]} @{$gene{$name}{$id}{cds}};
		if($gene{$name}{$id}{sign} eq '+'){
			for my $cds(@{$gene{$name}{$id}{cds}}){
				$out.="$cds->[0]\t$cds->[1]\t$cds->[2]\n";
			}
		}else{
			for my $cds(reverse @{$gene{$name}{$id}{cds}}){
				$out.="$cds->[0]\t$cds->[2]\t$cds->[1]\n";
			}
		}
		$out.="\n";
	}
}

print $out;


##############################################
sub read_gff{
	my ($file,$gene)=@_;
	open(IN,$file)||die"$!\n";
	while(<IN>){
		next if(/^#/);
		my @a=split /\t+/;
		my ($name,$id,$sign);
		next unless($a[2] =~ /CDS/i);
		$name=$a[0];
		@a[3,4]=@a[4,3] if ($a[3] > $a[4]);
                #$id=$1 if($a[8]=~/Parent=([^;\s]+)/);
		#$gene{$id}{name}=$name;
		$id=$1 if($a[8]=~/Parent=([^;]*)/);
                $id=~s/\s$//;
                $sign=$a[6];
		$gene{$id}{$name}{sign}=$sign;
		push @{$gene{$id}{$name}{cds}},[$name,$a[3],$a[4]];
		$gene{$id}{$name}{num}++;
	}
	close IN;
}
