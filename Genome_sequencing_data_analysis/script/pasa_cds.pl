#!/usr/bin/env perl
while(<>){
	if(/CDS/){
		my @aa = split;
		$_ =~ s/ID=(\S+)/ID=$aa[3]_$aa[4]$1/;
	}
	print $_;
}
