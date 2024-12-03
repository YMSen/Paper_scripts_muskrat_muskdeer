#!/usr/bin/perl

# mergeResult.pl  03/04/19 09:14:29  

use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;

my ($help);
GetOptions(
	'help'=>\$help,
);

pod2usage 1 if($help);

my $dirlist = shift;

open FL, $dirlist or die "Cannot read dirlist";
my @flist = <FL>;
close FL;
chomp for @flist;


&mergeTrunc(@flist);
&mergeMapping(@flist);
&mergeFilter(@flist);
&mergeDedup(@flist);
&mergeDitagSize(@flist);


=head1 SYNOPSIS

mergeResult.pl --help

=head1 OPTIONS

 --help        help
=cut


sub mergeTrunc{
    my @flist = @_;
    my %data;
    for my $d(@flist){
        open RES,"$d/hicup_truncater_summary_work.txt";
        $data{'title'} = <RES>;
        while(<RES>){
            chomp;
            my @tem = split /\t/,$_;
            my $file;
            if($tem[0] =~ /^left/){
                $file = 'left';
            }elsif($tem[0] =~ /^right/){
                $file = 'right';
            }else{
                die "filename err\n";
            }
            $data{$file}{'Total_Reads_Processed'} += $tem[1];
            $data{$file}{'Truncated'} += $tem[2];
            $data{$file}{'Not_truncated'} += $tem[4];
            $data{$file}{'length_truncated'} += $tem[2] * $tem[6];
        }
        close RES;
    }
    system("mkdir merged_report") unless(-d "merged_report");

    my $file;
    open FLS,">merged_report/hicup_truncater_summary_work.txt";
    print FLS $data{'title'};
    $file = 'left';
    print FLS +(join "\t",(
        $file,
        $data{$file}{'Total_Reads_Processed'},
        $data{$file}{'Truncated'}, 
        sprintf("%.2f", $data{$file}{'Truncated'} / $data{$file}{'Total_Reads_Processed'} * 100),
        $data{$file}{'Not_truncated'},
        sprintf("%.2f", $data{$file}{'Not_truncated'} / $data{$file}{'Total_Reads_Processed'} * 100),
        sprintf("%.2f", $data{$file}{'length_truncated'} / $data{$file}{'Truncated'}),
    )),"\n";

    $file = 'right';
    print FLS +(join "\t",(
        $file,
        $data{$file}{'Total_Reads_Processed'},
        $data{$file}{'Truncated'}, 
        sprintf("%.2f", $data{$file}{'Truncated'} / $data{$file}{'Total_Reads_Processed'} * 100),
        $data{$file}{'Not_truncated'},
        sprintf("%.2f", $data{$file}{'Not_truncated'} / $data{$file}{'Total_Reads_Processed'} * 100),
        sprintf("%.2f", $data{$file}{'length_truncated'} / $data{$file}{'Truncated'}),
    )),"\n";
    close FLS;
}


sub mergeMapping{
    my @flist = @_;
    my %data;
    for my $d(@flist){
        open RES,"$d/hicup_mapper_summary_work.txt";
        $data{'title'} = <RES>;
        while(<RES>){
            chomp;
            my @tem = split /\t/,$_;
            my $file;
            if($tem[0] =~ /^left/){
                $file = 'left';
            }elsif($tem[0] =~ /^right/){
                $file = 'right';
            }else{
                die "filename err\n";
            }
            $data{$file}{'Total_Reads_Processed'} += $tem[1];
            $data{$file}{'Reads_too_short_to_map'} += $tem[2];
            $data{$file}{'Unique_alignments'} += $tem[4];
            $data{$file}{'Multiple_alignments'} += $tem[6];
            $data{$file}{'Failed_to_align'} += $tem[8];
            $data{$file}{'Paired'} += $tem[10];
        }
        close RES;
    }
    system("mkdir merged_report") unless(-d "merged_report");

    my $file;
    open FLS,">merged_report/hicup_mapper_summary_work.txt";
    print FLS $data{'title'};
    $file = 'left';
    print FLS +(join "\t",(
        $file,
        $data{$file}{'Total_Reads_Processed'},
        $data{$file}{'Reads_too_short_to_map'},
        sprintf("%.1f", $data{$file}{'Reads_too_short_to_map'}/$data{$file}{'Total_Reads_Processed'}*100),
        $data{$file}{'Unique_alignments'},
        sprintf("%.1f", $data{$file}{'Unique_alignments'}/$data{$file}{'Total_Reads_Processed'}*100),
        $data{$file}{'Multiple_alignments'},
        sprintf("%.1f", $data{$file}{'Multiple_alignments'}/$data{$file}{'Total_Reads_Processed'}*100),
        $data{$file}{'Failed_to_align'},
        sprintf("%.1f", $data{$file}{'Failed_to_align'}/$data{$file}{'Total_Reads_Processed'}*100),
        $data{$file}{'Paired'},
        sprintf("%.1f", $data{$file}{'Paired'}/$data{$file}{'Total_Reads_Processed'}*100),
    )),"\n";

    $file = 'right';
    print FLS +(join "\t",(
        $file,
        $data{$file}{'Total_Reads_Processed'},
        $data{$file}{'Reads_too_short_to_map'},
        sprintf("%.1f", $data{$file}{'Reads_too_short_to_map'}/$data{$file}{'Total_Reads_Processed'}*100),
        $data{$file}{'Unique_alignments'},
        sprintf("%.1f", $data{$file}{'Unique_alignments'}/$data{$file}{'Total_Reads_Processed'}*100),
        $data{$file}{'Multiple_alignments'},
        sprintf("%.1f", $data{$file}{'Multiple_alignments'}/$data{$file}{'Total_Reads_Processed'}*100),
        $data{$file}{'Failed_to_align'},
        sprintf("%.1f", $data{$file}{'Failed_to_align'}/$data{$file}{'Total_Reads_Processed'}*100),
        $data{$file}{'Paired'},
        sprintf("%.1f", $data{$file}{'Paired'}/$data{$file}{'Total_Reads_Processed'}*100),
    )),"\n";
    close FLS;
}

sub mergeFilter{
    my @flist = @_;
    my %data;
    my @data;
    $data[0] = 'pair';
    for my $d(@flist){
        open RES,"$d/hicup_filter_summary_work.txt";
        $data{'title'} = <RES>;
        while(<RES>){
            chomp;
            my @tem = split /\t/,$_;
            for(1 .. $#tem){
                $data[$_] += $tem[$_]; 
            }
        }
        close RES;
    }
    system("mkdir merged_report") unless(-d "merged_report");
    open FLS,">merged_report/hicup_filter_summary_work.txt";
    print FLS $data{'title'};
    print FLS +(join "\t", @data),"\n";
    close FLS;
}

sub mergeDedup{
    my @flist = @_;
    my %data;
    my @data;
    $data[0] = 'pairs';
    for my $d(@flist){
        open RES,"$d/hicup_deduplicator_summary_work.txt";
        $data{'title'} = <RES>;
        while(<RES>){
            chomp;
            my @tem = split /\t/,$_;
            for(1 .. $#tem){
                $data[$_] += $tem[$_]; 
            }
        }
        close RES;
    }
    system("mkdir merged_report") unless(-d "merged_report");
    open FLS,">merged_report/hicup_deduplicator_summary_work.txt";
    print FLS $data{'title'};
    print FLS +(join "\t", @data),"\n";
    close FLS;
}

sub mergeDitagSize{
    my @flist = @_;
    my %data;
    my $title;
    for my $d(@flist){
        open RES,"$d/left.right.ditag_size_distribution_HTML_report.temp";
        $title = <RES>;
        while(<RES>){
            chomp;
            my @tem = split /\t/,$_;
            $data{$tem[0]} += $tem[1];
        }
        close RES;
    }
    system("mkdir merged_report") unless(-d "merged_report");
    open FLS,">merged_report/left.right.ditag_size_distribution_HTML_report.temp";
    print FLS $title;
    for(sort {$a<=>$b} keys %data){
        print FLS "$_\t$data{$_}\n";
    }
    close FLS;
}

