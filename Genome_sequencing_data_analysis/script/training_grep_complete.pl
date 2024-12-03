#! use/bin/perl -w
use strict;
die "open<fa><OUT>\n" unless @ARGV==2;

if($ARGV[0]=~/\.gz$/){
        open(DATA,"gunzip -dc $ARGV[0]|")||die;
}else{
        open(DATA,$ARGV[0])||die;
}

open(OUT,">$ARGV[1]")||die;
$/=">";
<DATA>;
while(<DATA>)
{
        /(.*)\n/;
        my $info=$1;
        my ($type) = $info =~ /type:([^\s]+)/;
        if($type eq 'complete')
        {
                s/.+\n//;
                s/>//g;
	print OUT ">$info\n";
        print OUT $_;
        }
}
$/ = "\n";
close DATA;
close OUT;

