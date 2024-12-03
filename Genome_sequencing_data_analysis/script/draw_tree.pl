#!/usr/bin/env perl
use strict;
use FindBin qw($Bin);
use lib "$Bin/lib";
use Tree::nhx_svg;
use Getopt::Long;

my $usage = <<END;

 perl $0 <in.nhx> <in.desc> [option]
  <in.nhx>           nhx or newick format tree file;
  <in.desc>          name, scientific name, common name;
  --width <num>      tree figure width;
  --cali <file>      input tree for mcmctree with calibrations;

END

my ($width, $cali) = (600, "");
GetOptions (
	'width:f' => \$width,
	'cali:s' => \$cali,
);

die $usage if (@ARGV != 2);
my $innhx = shift;
my $indesc = shift;

# Read in tree;
$/ = ";";
open (IN,"$innhx") or die "Can't open file $innhx: $!\n";
my $str = <IN>;
my $str2 = <IN>;
if (length($str2) > length($str)) 
{
	$str = $str2;
	$str =~ s/^\s+//g;
	print STDERR "$str\n";
}
close IN;

my $nhx = Tree::nhx_svg->new('show_B',0,'show_W',0,'show_ruler',1,'show_interval',1,
	'c_line',"#000000",'line_width',1.5,'c_W',"#000000",
	'c_B',"#000000",'c_node',"#880000",'fsize',16, "width", $width, "node_fsize", 12);
$nhx->parse($str);

# Read in cali tree;
if (defined $cali && $cali ne "") {
	$/ = ";";
	open (IN,"$cali") or die "Can't open file $cali: $!\n";
	$str = <IN>;
	close IN;
	my $nhx_tmp = Tree::nhx_svg->new('show_B',0,'show_W',0,'show_ruler',1,'show_interval',1,
		'c_line',"#000000",'line_width',1.5,'c_W',"#000000",
		'c_B',"#000000",'c_node',"#FF0000",'fsize',16, "width", $width);
	$nhx_tmp->parse($str);
	$nhx->add_cali($nhx_tmp);
}

# Read in scientific names;
$/ = "\n";
open (IN,"$indesc") or die "Can't open file $indesc: $!\n";
my @arr1 = <IN>;
close IN;

my %scinames;
my %comnames;
foreach (@arr1) {
	my @arr2 = split /\t+/;
	my @arr3 = split /\s+/, $arr2[1];
	my $sci = substr($arr3[0], 0, 1).".".substr($arr3[-1], 0, 3);
	$scinames{$arr2[0]} = "\u$sci";
	$comnames{$arr2[0]} = "\u$arr2[2]" if (defined $arr2[2] && $arr2[2] ne "");
}
$nhx->add_character("scientific", \%scinames);
$nhx->add_character("common", \%comnames);

#print $nhx->info();
print $nhx->plot();

