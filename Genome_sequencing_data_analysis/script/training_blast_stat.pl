#!/usr/bin/perl -w
use strict;
die "use <tab>" unless @ARGV==1;

my %geneAlign;

my $tabfile=shift;

open(IN,"$tabfile")||die;
<IN>;
while(<IN>){
      chomp;
      my @cut=split /\t+/;
      $geneAlign{$cut[0]}{l} = $cut[1];
      push @{$geneAlign{$cut[0]}{a}},[$cut[2],$cut[3]];
}
close IN;

foreach my $key(keys %geneAlign){
        @{$geneAlign{$key}{a}}=sort{$a->[0]<=>$b->[0]}@{$geneAlign{$key}{a}};
        my $st=${$geneAlign{$key}{a}}[0][0];
        my $en=${$geneAlign{$key}{a}}[0][1];
        my $sum=0;
        my $site="";
        my $mostlen=$en-$st+1;
        for(my $i=1;$i<@{$geneAlign{$key}{a}};$i++){
                    my $start=${$geneAlign{$key}{a}}[$i][0];
                    my $end=${$geneAlign{$key}{a}}[$i][1];
                    if(($en-$st+1)>$mostlen){$mostlen=$en-$st+1;}
                    if($start>$en){
                             #print "$key\t$st\t$en\n";
                             $sum+=($en-$st+1);
                             $site.="$st,$en;";
                             $st=$start;
                             $en=$end;
                    }elsif($end>$en){
                             $en=$end;
                    }
        }
        #print "$key\t$st\t$en\n";
        $sum+=($en-$st+1);
        $site.="$st,$en;";
        my $ratio=$sum/$geneAlign{$key}{l};
        my $mostratio=$mostlen/$geneAlign{$key}{l};
        print "$key\t$geneAlign{$key}{l}\t$sum\t$ratio\t$site\t$mostlen\t$mostratio\n";
}
