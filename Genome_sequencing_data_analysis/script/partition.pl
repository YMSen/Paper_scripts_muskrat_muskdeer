#!/usr/bin/perl

# partition.pl  10/15/19 14:59:40 

use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;

my ($help,$contigfile,$pairsfile,$K,$minREs,$MaxLinkDensity,$MinAvgLinkage,$NonInformativeRatio,$LengthCutOff);
GetOptions(
	'help'                  =>\$help,
    'contigfile=s'          =>\$contigfile,
    'pairsfile=s'           =>\$pairsfile,
    'K=i'                   =>\$K,
    'minREs:i'              =>\$minREs,
    'maxlinkdensity:i'      =>\$MaxLinkDensity,
    'minavglinkage:i'       =>\$MinAvgLinkage,
    'NonInformativeRabio:i' =>\$NonInformativeRatio,
    'lencut:f'              =>\$LengthCutOff,
);

pod2usage 1 if($help);

$minREs ||= 25;
$MinAvgLinkage ||= 0;
$MaxLinkDensity ||= 2;
$NonInformativeRatio ||= 2;
$LengthCutOff ||= 0.98;

my $largestRE = 0;
my @matrix;
my %contigstoidx;
my %contigs;
my %pairs;
my %clusters;


&getRE();
&makeMatrix();
&skipRepeats();
&cluster();
&printClusters();

sub getRE{
    open CTF,$contigfile || die "no contigs file\n";
    my $No = 0;
    my $genomeTotalLen = 0;
    my $chooseLen = 0;
    while (<CTF>){
        chomp;
        next if ($_ =~ /^#/);
        my @line = split;
        $contigstoidx{$line[0]} = $No;
        if ($line[1] > $largestRE){
            $largestRE = $line[1];
        }
        $contigs{$No}->{"name"} = $line[0];
        $contigs{$No}->{"REcounts"} = $line[1];
        $contigs{$No}->{"length"} = $line[2];
        $genomeTotalLen += $line[2];
        if ($line[1] < $minREs){
            print "contig $line[0] has $line[1] RE sites ------ marked short\n";
            $contigs{$No}->{"skip"} = "YES";
        }else{
            $contigs{$No}->{"skip"} = "NO";
        }
        $No++;
    }
    close CTF;
    $chooseLen = $genomeTotalLen * $LengthCutOff;
    my $addLen = 0;
    foreach my $i (sort{$contigs{$b}{"length"} <=> $contigs{$a}{"length"}} keys %contigs){
        $addLen += $contigs{$i}{"length"};
        if ($addLen > $chooseLen){
            $contigs{$i}->{"skip"} = "YES";
        }
    }
}

sub makeMatrix{
    my $allidx = scalar(keys %contigs) - 1;
    foreach my $i(0..$allidx){
        foreach my $l(0..$allidx){
            $matrix[$i][$l] = 0;
        }
    }
    my $largestSquared = $largestRE * $largestRE;
    open PAIR,$pairsfile || die "no pairs file \n";
    my $No = 0;
    while (<PAIR>){
        chomp;
        next if ($_ =~ /^#/);
        my @line = split;
        my $ctgAidx = $contigstoidx{$line[2]};
        my $ctgBidx = $contigstoidx{$line[3]};
        $pairs{$No}->{"contigA"} = $line[2];
        $pairs{$No}->{"contigB"} = $line[3];
        $pairs{$No}->{"RE1"} = $line[4];
        $pairs{$No}->{"RE2"} = $line[5];
        $pairs{$No}->{"observedlinks"} = $line[6];
        if ($line[4] ne $line[5]){
            my $w = $line[6] * $largestSquared / ($line[4] * $line[5]);
            $matrix[$ctgAidx][$ctgBidx] = $w;
            $matrix[$ctgBidx][$ctgAidx] = $w;
        }
    }
    close PAIR;
}
        
sub skipRepeats{
    my $allidx = scalar(keys %contigs) - 1;
    my %ctgLinks;
    my $totalLinks = 0;
    foreach my $i(0..$allidx){
        foreach my $j($i+1..$allidx){
            $totalLinks += $matrix[$i][$j];
            $ctgLinks{$i} += $matrix[$i][$j];
            $ctgLinks{$j} += $matrix[$i][$j];
        }
    }
    my $ctgLinksAvg = 2.0 * $totalLinks / ($allidx + 1);
    foreach my $l(0..$allidx){
        my $factor = $ctgLinks{$l} / $ctgLinksAvg;
        foreach my $g(0..$allidx){
            if ($matrix[$l][$g] != 0){
                $matrix[$l][$g] = $matrix[$l][$g] / $factor;
            }
        }
        if ($factor >= $MaxLinkDensity){
            my $repeatcontig = $contigs{$l}->{"name"};
            print "contig $repeatcontig has $factor x the average number of Hi-C links ------ marked repeat\n";
            $contigs{$l}->{"skip"} = "YES";
        }
    }
}

sub cluster{
    my %clusterID;
    my %clusterSize;
    my %clusterExists;
    my $nonSingletonClusters = 0;
    my $skipedContigs = 0;
    my $allidx = scalar(keys %contigs) - 1;
    foreach my $i(0..$allidx){
        if ($contigs{$i}->{"skip"} eq "YES"){
            $clusterID{$i} = -1;
            $skipedContigs++;
        }else{
            $clusterID{$i} = $i;
            $clusterSize{$i} = 1;
            $clusterExists{$i} = "YES";
        }
    }
    my $nonskipedContigs = scalar(keys %contigs) - $skipedContigs;
    my %merges;
    my $nMerges = 0;
    foreach my $i(0..$allidx){
        next if ($contigs{$i}->{"skip"} eq "YES");
        foreach my $l($i+1..$allidx){
            if ($contigs{$l}->{"skip"} eq "NO" && $matrix[$i][$l] > $MinAvgLinkage){
                $merges{$i}->{$l} = $matrix[$i][$l];
            }
        }
    }
    my $tenu = 1;
    while(1){
        if (scalar(keys %merges) == 0){
            print "No more merges to do since the queue is empty\n";
            last;
        }

        #choose the best score clusters;
        my %bestMerge;
        my $topa = (keys %merges)[0];
        my $topb = (keys %{$merges{$topa}})[0];
        $bestMerge{"clusterA"} = $topa;
        $bestMerge{"clusterB"} = $topb;
        $bestMerge{"score"} = $merges{$topa}->{$topb};
        foreach my $i(keys %merges){
            foreach my $j(keys %{$merges{$i}}){
                if ($merges{$i}->{$j} > $bestMerge{"score"}){
                    $bestMerge{"clusterA"} = $i;
                    $bestMerge{"clusterB"} = $j;
                    $bestMerge{"score"} = $merges{$i}->{$j};
                }
            }
        }
        
        #foreach my $i(0..$allidx){
        #    if ($clusterID{$i} == $bestMerge{"clusterA"}){
        #        my $absc = $contigs{$i}{"name"};
        #        print "$absc\t";
        #    }
        #}
        #print "\n";
        #foreach my $i(0..$allidx){
        #    if ($clusterID{$i} == $bestMerge{"clusterB"}){
        #        my $absc = $contigs{$i}{"name"};
        #        print "$absc\t";
        #    }
        #}
        #print "\n";
        #my $asds = $bestMerge{"score"};
        #print "$tenu-------score:$asds\n";
        #$tenu++;

        #generate the new cluster;
        my $newClusterID = scalar(keys %contigs) + $nMerges;
        $clusterExists{$bestMerge{"clusterA"}} = "NO";
        $clusterExists{$bestMerge{"clusterB"}} = "NO";
        $clusterExists{$newClusterID} = "YES";
        $clusterSize{$newClusterID} = $clusterSize{$bestMerge{"clusterA"}} + $clusterSize{$bestMerge{"clusterB"}};
        if ($bestMerge{"clusterA"} < scalar(keys %contigs)){
            $nonSingletonClusters++;
        }
        if ($bestMerge{"clusterB"} < scalar(keys %contigs)){
            $nonSingletonClusters++;
        }
        $nonSingletonClusters--;
        my @newCluster;
        foreach my $i(0..$allidx){
            if ($clusterID{$i} == $bestMerge{"clusterA"} || $clusterID{$i} == $bestMerge{"clusterB"}){
                $clusterID{$i} = $newClusterID;
                push @newCluster,$i;
            }
        }
        $nMerges++;
        
        #Calculate new score entries for the new cluster and remove all used clusters;
        my %newMerges;
        foreach my $i(keys %merges){
            foreach my $j(keys %{$merges{$i}}){
                if ($clusterExists{$i} eq "YES" && $clusterExists{$j} eq "YES"){
                    $newMerges{$i}->{$j} = $merges{$i}->{$j};
                }
            }
        }
        my %LinktoNewCluster;
        my %clustersContig;
        my %AvgLinkAll;
        foreach my $i(0..$allidx){
            next if ($clusterID{$i} == $newClusterID || $clusterID{$i} == -1);
            push @{$clustersContig{$clusterID{$i}}},$i;
        }
        my $clusterLenB = 0;
        foreach my $cB(@newCluster){
            $clusterLenB += $contigs{$cB}->{"length"};
        }
        foreach my $i(keys %clustersContig){
            my $clusterANu = scalar(@{$clustersContig{$i}});
            my $clusterBNu = scalar(@newCluster);
            my $clusterLenA = 0;
            foreach my $cA(@{$clustersContig{$i}}){
                $clusterLenA += $contigs{$cA}->{"length"};
            }
            my %linkToCluster;
            if ($clusterLenA > $clusterLenB){
                foreach my $j(@{$clustersContig{$i}}){
                    foreach my $l(@newCluster){
                        $linkToCluster{$j} += $matrix[$j][$l];
                    }
                }
                my $ActgNu = 0;
                my $ActgLenAll = 0;
                foreach my $j(sort{$linkToCluster{$b} <=> $linkToCluster{$a} || $b <=> $a} keys %linkToCluster){
                    if ($ActgLenAll <= $clusterLenB){
                        foreach my $m(@newCluster){
                            if($j<$m){
                                $LinktoNewCluster{$i} += $matrix[$j][$m];
                            }else{
                                $LinktoNewCluster{$i} += $matrix[$m][$j];
                            }
                            #$LinktoNewCluster{$i} += $matrix[$j][$m];
                            $ActgLenAll += $contigs{$j}->{"length"};
                            $ActgNu ++;
                        }
                    }else{
                        last;
                    }
                }
                $AvgLinkAll{$i} = $LinktoNewCluster{$i} / $clusterBNu / $ActgNu;
            }else{
                foreach my $j(@newCluster){
                    foreach my $l(@{$clustersContig{$i}}){
                        $linkToCluster{$j} += $matrix[$l][$j];
                    }
                }
                my $BctgNu = 0;
                my $BctgLenAll = 0;
                foreach my $j(sort{$linkToCluster{$b} <=> $linkToCluster{$a} || $b <=> $a} keys %linkToCluster){
                    if ($BctgLenAll <= $clusterLenA){
                        foreach my $m(@{$clustersContig{$i}}){
                            if($j<$m){
                                $LinktoNewCluster{$i} += $matrix[$j][$m];
                            }else{
                                $LinktoNewCluster{$i} += $matrix[$m][$j];
                            }
                            #$LinktoNewCluster{$i} += $matrix[$m][$j];
                            $BctgLenAll += $contigs{$j}->{"length"};
                            $BctgNu ++;
                        }
                    }else{
                        last;
                    }
                }
                $AvgLinkAll{$i} = $LinktoNewCluster{$i} / $clusterANu / $BctgNu;
            }
        }
        my $dallidx = 2 * $allidx;
        foreach my $i(0..$dallidx){
            next if (!exists $LinktoNewCluster{$i} || $LinktoNewCluster{$i} <= 0);
            my $AvgLink = $AvgLinkAll{$i};
            next if ($AvgLink < $MinAvgLinkage);
            my $clusta;
            my $clustb;
            if ($i <= $newClusterID ){
                $clusta = $i;
                $clustb = $newClusterID;
            }else{
                $clusta = $newClusterID;
                $clustb = $i;
            }
            $newMerges{$clusta}->{$clustb} = $AvgLink;
        }

        #Analyze the current clusters if enough merges occurred;
        if ($nMerges > $nonskipedContigs/2 && $nonSingletonClusters <= $K){
            if ($nonSingletonClusters == $K){
                last;
            }
        }
        %merges = %newMerges;
    }
    &skipContigClusters(\%clusterID);
}

sub skipContigClusters{
    my %clusterID = %{$_[0]};
    my $allidx = scalar(keys %contigs) - 1;
    my %skipContigClusters;
    foreach my $i(sort{$a <=> $b} keys %clusterID){
        next if ($clusterID{$i} == $i || $clusterID{$i} == -1);
        push @{$clusters{$clusterID{$i}}},$i;
    }
    if (!($NonInformativeRatio == 0 || $NonInformativeRatio > 1)){
        print "NonInformativeRatio needs to either 0 or > 1\n";
        exit(1);
    }
    foreach my $i(0..$allidx){
        next if ($clusterID{$i} != -1);
        my %linkages = &ctgToClusterLinkages($i);
        next if (scalar(keys %linkages) == 0);
        my $largestLinkageCluster = (sort{$linkages{$b} <=> $linkages{$a}} keys %linkages)[0];
        my $largestLinkage = $linkages{$largestLinkageCluster};
        my $top2LinkageCluster;
        my $top2Linkage = 0;
        if (scalar(keys %linkages) > 1){
            $top2LinkageCluster = (sort{$linkages{$b} <=> $linkages{$a}} keys %linkages)[1];
            $top2Linkage = $linkages{$top2LinkageCluster};
        }
        if ($largestLinkage >= $NonInformativeRatio && (scalar(keys %linkages) ==1 || $top2Linkage == 0 || $largestLinkage/$top2Linkage >= $NonInformativeRatio)){
            $skipContigClusters{$i} = $largestLinkageCluster;
        }
    }
    foreach my $i(keys %skipContigClusters){
        push @{$clusters{$skipContigClusters{$i}}},$i;
    }
}

sub ctgToClusterLinkages{
    my %linkages;
    my $skipContig = $_[0];
    foreach my $i(keys %clusters){
        my $clusterSize = scalar(@{$clusters{$i}});
        my $totalLinkage = 0;
        foreach my $j(@{$clusters{$i}}){
            if ($skipContig == $j){
                $clusterSize--;
            }else{
                $totalLinkage += $matrix[$skipContig][$j];
            }
        }
        if ($totalLinkage > 0){
            $linkages{$i} = $totalLinkage / $clusterSize;
        }
    }
    return %linkages;
}

sub printClusters{
    my %clusterLens;
    foreach my $i(keys %clusters){
        my @contigList = @{$clusters{$i}};
        foreach my $j(@contigList){
            $clusterLens{$i} += $contigs{$j}->{"length"};
        }
    }
    open CLU,">group.cluster.txt";
    print CLU "#Group\tnContigs\tContigs\n";
    my $cluN = 1;
    foreach my $i(sort{$clusterLens{$b} <=> $clusterLens{$a}} keys %clusterLens){
        open GROUP,">group$cluN.txt";
        my @contigList = @{$clusters{$i}};
        my $contigN = scalar(@contigList);
        print GROUP "#Contig\tRECounts\tLength\n";
        print CLU "$cluN\t$contigN\t";
        foreach my $j(@contigList){
            my $contigName = $contigs{$j}->{"name"};
            print CLU "$contigName ";
            my $contigREs = $contigs{$j}->{"REcounts"};
            my $contigLen = $contigs{$j}->{"length"};
            print GROUP "$contigName\t$contigREs\t$contigLen\n";
        }
        print CLU "\n";
        $cluN += 1;
        close GROUP;
    }
    close CLU;
    print "cluster is done!\n";
}

#sub printClusters{
#    my %clusterLens;
#    foreach my $i(keys %clusters){
#        my @contigList = @{$clusters{$i}};
#        foreach my $j(@contigList){
#            $clusterLens{$i} += $contigs{$j}->{"length"};
#        }
#    }
#    open CLU,">group.cluster.txt";
#    print CLU "#Group\tnContigs\tContigs\n";
#    my $cluN = 1;
#    foreach my $i(sort{$clusterLens{$b} <=> $clusterLens{$a}} keys %clusterLens){
#        my @contigList = @{$clusters{$i}};
#        my $contigN = scalar(@contigList);
#        print CLU "$cluN\t$contigN\t";
#        foreach my $j(@contigList){
#            my $contigName = $contigs{$j}->{"name"};
#            print CLU "$contigName ";
#        }
#        print CLU "\n";
#        $cluN += 1;
#    }
#    close CLU;
#    print "cluster is done!\n";
#}


        






=head1 SYNOPSIS

partition.pl --help

=head1 OPTIONS

 --help        help
=cut

