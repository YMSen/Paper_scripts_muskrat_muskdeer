#!/usr/bin/perl
use strict;
my $dir=shift;
opendir THIS,$dir or die "$!";
my @fa=grep {/\d+$/}readdir THIS;
close THIS;
open OUT,">find_pals.sh";

my %seq_num;
for my $file(@fa){
    my $count;
    open FL,$dir.$file;
    while(<FL>){
        $count ++ if(/^>/);
    }
    close FL;
    $seq_num{$file} = $count;
}

@fa = sort {$seq_num{$b}<=>$seq_num{$a}} @fa;

foreach my $i(0..$#fa)
{
    my $filename=$dir.$fa[$i];
    print OUT "./pals -self $filename -out $fa[$i].hit.gff\n";
    for (my $j=$i+1;$j<@fa;$j++)
    {
        my $file=$dir.$fa[$j];
        print OUT "./pals -target $filename -query $file -out $fa[$i].$fa[$j].hit.gff\n";
    }		
}
close OUT;

