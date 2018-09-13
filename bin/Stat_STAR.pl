#!/usr/bin/perl
use strict;
use warnings;

die "\nUsage:\n\tperl  $0  Sample1.Log.final.out  Sample2.Log.final.out ...  >Align_stat.txt\n\n" if @ARGV < 1;

print "SampleID\tInputReads/Pairs\ttotalMapped\ttotalMapped(%)\tuniqMapped(%)\tunMapped(%)\n";

foreach my $file (@ARGV) {
    my ( $samp, $input, $totalM, $totalMp, $uniqMp, $unMp ) = ("0") x 6;
    ($samp) = $file =~ /([^\/]+)Log.final.out$/;
    open IN, $file or die $!;
    my @lines = <IN>;
    foreach (@lines) {
        $input  = $1 if /Number of input reads \|\s+(.*)\n/;
        $uniqMp = $1 if /Uniquely mapped reads % \|\s+(.*)%/;
        $totalMp += $1 if /Uniquely mapped reads % \|\s+(.*)%/ or /of reads mapped to multiple loci \|\s+(.*)%/;
        $totalM += $1 if /Uniquely mapped reads number \|\s+(\d+)\n/ or /Number of reads mapped to multiple loci \|\s+(\d+)\n/;
    }
    close IN;
    $unMp = sprintf( "%.2f", 100 - $totalMp );
    print "$samp\t$input\t$totalM\t$totalMp\t$uniqMp\t$unMp\n";
}
