#!/usr/bin/env perl
my %trans;
my $fname = shift;

open FL,$fname or die("Cannot open pasa gff3 file: '$fname'");
open PEP,">pasa2.longest.pep";
open GFF,">pasa2.longest.gff";
open CDS,">pasa2.longest.cds";

my @geneline;
while(<FL>){
    chomp;
    my @tem = split /\t/;
    if(@tem == 9){
        if($tem[2] eq 'gene'){
            # print GFF "$_\n";
            @geneline = @tem;
        }elsif($tem[2] eq 'mRNA'){
            if($tem[8] =~ /ID=([^;]+)/){
                $trans{$1}{'block'} .= "$_\n";
                @geneline_m = @geneline;
                $geneline_m[3] = $tem[3];
                $geneline_m[4] = $tem[4];
                $trans{$1}{'geneline'} = [@geneline_m];
            }else{die;}
        }else{
            if($tem[8] =~ /Parent=([^;]+)/){
                $trans{$1}{'block'} .= "$_\n";
                if($tem[2] eq 'CDS'){
                    if($tem[8] =~ /Parent=([^;]+)/){
                        $trans{$1}{'len'} += abs($tem[4] - $tem[3]) + 1;
                    }else{die;}
                }
            }else{die;}
        }
    }elsif(/^#CDS/){
        my @tem = split /\s+/;
        $trans{$tem[1]}{'cds'} = $tem[3];
    }elsif(/^#PROT/){
        my @tem = split /\s+/;
        $trans{$tem[1]}{'prot'} = $tem[3];
    }elsif(/^# /){
        for(sort {$trans{$b}{'len'}<=>$trans{$a}{'len'} || $a cmp $b} keys %trans){
            my $gline = join "\t",@{$trans{$_}{'geneline'}};
            print GFF "$gline\n$trans{$_}{'block'}";
            print PEP ">$_\n$trans{$_}{'prot'}\n";
            print CDS ">$_\n$trans{$_}{'cds'}\n";
            last;
        }
        undef %trans;
    }
}
close FL;

for(sort {$trans{$b}{'len'}<=>$trans{$a}{'len'} || $a cmp $b} keys %trans){
    my $gline = join "\t",@{$trans{$_}{'geneline'}};
    print GFF "$gline\n$trans{$_}{'block'}";
    print PEP ">$_\n$trans{$_}{'prot'}\n";
    print CDS ">$_\n$trans{$_}{'cds'}\n";
    last; 
}   
close PEP;
close GFF;
close CDS;

