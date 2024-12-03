#!/usr/bin/perl
use strict;
my$dir=$ARGV[1];
my%cprops;
my%groupctg;
open IN1,"$ARGV[0]"||die $!;
open OUT,">genome.asm";
while(<IN1>){
	chomp;
	my($id,$num)=(split)[0,1];
	$cprops{$id}=$num;
}
my$toura=`ls $dir | grep '.tour\$'|wc -l`;
foreach my$t(1..$toura){
    my$tour="group".$t.".tour";
#while(my$tour=glob "$dir/*.tour"){
    my $tour1=$dir."/".$tour;
	my $last_line=`tail -n 1 $tour1`;
	my @ctgdb=split (/\s+/,$last_line);
	my @ctg_num;
	foreach my$i(0..$#ctgdb){
		my $ctg;my$fdir;my$rdir;my$ctgnew;
		if($ctgdb[$i]=~/(.*)\+$/){
			$ctg=$1;
			$ctgnew=$cprops{$ctg};
		}
		if($ctgdb[$i]=~/(.*)\-$/){
			$ctg=$1;
                        $ctgnew="-".$cprops{$ctg};
		}
		push(@ctg_num,$ctgnew);
		$groupctg{$ctg}=();
	}
	my$group=join(" ",@ctg_num);
	print OUT "$group\n";
}
foreach my$ctgid (keys %cprops){
	if(!exists $groupctg{$ctgid}){
		my$nogroup=$cprops{$ctgid};
		print OUT "$nogroup\n";
	}
}
close IN1;
close OUT;
