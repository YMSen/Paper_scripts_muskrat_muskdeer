package pply;

use strict;
use Exporter;
use IPC::Open2;
use File::Spec;

use vars qw($VERSION @ISA @EXPORT @EXPORT_OK %EXPORT_TAGS);

$VERSION     = 1.00;
@ISA         = qw(Exporter);
@EXPORT      = qw(cmdexe waitjob abspath);
@EXPORT_OK   = qw();
%EXPORT_TAGS = (
    DEFAULT => [qw(&func1)],
    Both    => [qw(&func1 &func2)]
);

sub cmdexe {
    my ($cmd,$count) = @_;
    my $prefix = (split /\s+/,$cmd)[0];
    $prefix =~ s/.*\///;
    
    open2(*Reader, *Writer, "qsub -N $prefix -cwd -V");
    print Writer $cmd;
    close Writer;
    my $jobid = <Reader>;
    close Reader;
    
    ($jobid) = $jobid =~ /(\d+)/;
    my $maxmem = waitjob($jobid);
    print "MAXMEM\t$maxmem\n";
}

sub waitjob{
    my ($jobid) = @_;
    my $maxmem = 0;
    while(my $status = `qstat -j $jobid 2> /dev/null`){
        if($status =~ /(usage\s+1:.*?maxvmem=(.*))/){
            if($2 =~ /^([\d\.]+[BKMGT])/){
                $maxmem = $1;
            }
        }
        sleep 1;
    }
    return $maxmem;
}

sub abspath {  
    my ($path) = @_;  
    $path =~ s/^~/$ENV{HOME}/;
    $path = File::Spec->rel2abs($path);
    my $newpath;  
    do {  
        $newpath = $path;  
        $path =~ s|[^/][^/]*\/\.\.\/||;  
    } while $newpath ne $path;  
    return $newpath;  
}

1;
