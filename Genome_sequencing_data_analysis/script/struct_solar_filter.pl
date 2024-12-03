#!/usr/bin/perl
my $infile = shift;
my $align_rate = shift;
my %solardata;
my $output;
open IN, "$infile" || die "fail $infile";
while (<IN>) {
  chomp;
  s/^\s+//;
  my @t = split /\s+/;
  my $query = $t[0];
  my $score = $t[10];
  next if($score < 25);
  my $query_size = $t[1];
  my $align_size;
  while ($t[11]=~/(\d+),(\d+);/g) {
     $align_size += abs($2 - $1) + 1;
  }
  next if($align_size / $query_size < $align_rate);
  push @{$solardata{$query}},[$score,$_]; ## hits that better than cutoff
}
  open OUT, ">$infile.filter" || die "fail $infile.filter";
     foreach my $query (sort keys %solardata){
             my $pp = $solardata{$query};
                @$pp = sort {$b->[0] <=> $a->[0]} @$pp;
                for (my $i=0; $i<@$pp; $i++) {
                   last if(defined $Tophit && $i>=$Tophit);
                     my $query_Dup = "$query-D".($i+1);
                     $pp->[$i][1] =~ s/$query/$query_Dup/ if ($i>0);
                     print OUT $pp->[$i][1],"\n";
             }
}
 close OUT;
