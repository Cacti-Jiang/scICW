#!/usr/bin/perl
use strict;
use warnings;

die "\nUsage:\n\tperl  $0  Sample1.log  Sample2.log ...  >Reads_stat.txt\n\n" if @ARGV < 1;

print "SampleID\tRawReads/Pairs\tR1adapter\tR1adp(%)\tR2adapter\tR2adp(%)\tManyN\tManyNp(%)\tCleanReads/Pairs\tCp(%)\n";

foreach my $file (@ARGV) {
    my ( $samp, $total, $R1ad, $R1adp, $R2ad, $R2adp, $N, $Np, $pass, $pp ) = ("-") x 10;
    ($samp) = $file =~ /([^\/]+).log$/;
    open IN, $file or die $!;
    while (<IN>) {
        $total = $1 if /Total read pairs processed:\s+(\S+)/                or /Total reads processed:\s+(\S+)/;
        $R1ad  = $1 if /Read 1 with adapter:\s+(\S+)/                       or /Reads with adapters:\s+(\S+)/;
        $R1adp = $1 if /Read 1 with adapter:\s+\S+.*\((\S+)%/               or /Reads with adapters:\s+\S+.*\((\S+)%/;
        $R2ad  = $1 if /Read 2 with adapter:\s+(\S+)/;
        $R2adp = $1 if /Read 2 with adapter:\s+\S+.*\((\S+)%/;
        $N     = $1 if /Pairs with too many N:\s+(\S+)/                     or /Reads with too many N:\s+(\S+)/;
        $Np    = $1 if /Pairs with too many N:\s+\S+.*\((\S+)%/             or /Reads with too many N:\s+\S+.*\((\S+)%/;
        $pass  = $1 if /Pairs written \(passing filters\):\s+(\S+)/         or /Reads written \(passing filters\):\s+(\S+)/;
        $pp    = $1 if /Pairs written \(passing filters\):\s+\S+.*\((\S+)%/ or /Reads written \(passing filters\):\s+\S+.*\((\S+)%/;
        last if ( $_ =~ /===/ && $_ !~ /Summary/ );
    }
    close IN;
	$total =~ s/,//g;
	$R1ad =~ s/,//g;
	$R2ad =~ s/,//g;
	$N =~ s/,//g;
	$pass =~ s/,//g;
    print "$samp\t$total\t$R1ad\t$R1adp\t$R2ad\t$R2adp\t$N\t$Np\t$pass\t$pp\n";
}
