#!/usr/bin/perl
use strict;
use warnings;

die "\nUsage:\n\tperl  $0  rnaseq_qc_results.txt1  rnaseq_qc_results.txt2  ...  >Distribution_stat.txt\n\n" if @ARGV < 1;

print "SampleID\tMappedReads/Pairs\tExonic(%)\tIntronic(%)\tIntergenic(%)\n";

foreach my $file (@ARGV) {
    my ( $samp, $totalM, $exon, $intron, $intergene ) = ("-") x 5;
    open IN, $file or die $!;
    while (<IN>) {
        if (/bam file\s+=\s+(\S+)/){
                my $cap = $1;
                $cap =~ s/\.bam$//i;
                $cap =~ s/\.sorted$//;    # pipeline: hisat2
                $cap =~ s/Aligned\.sortedByCoord\.out$//;    # pipeline: STAR
                ($samp) = $cap =~ /([^\s^\/]+)$/;
                next;
        }
        if (/reads aligned\s+\(left\/right\)\s+=\s+([\d,]+)\s+\/\s+([\d,]+)/){
                #PE
                my ($num1,$num2) = ($1,$2);
                $num1 =~ s/,//g;
                $num2 =~ s/,//g;
                $totalM = sprintf ("%.0f", ($num1+$num2)/2 );
                next;
        }elsif(/reads aligned\s+=\s+([\d,]+)/){
                #SE
                $totalM = $1;
                $totalM =~ s/,//g;
                next;
        }
        $exon = $1 if /exonic\s+=\s+[\d,]+\s+\(([\d\.]+)%\)/;
        $intron = $1 if /intronic\s+=\s+[\d,]+\s+\(([\d\.]+)%\)/;
    }
    $intergene = sprintf ("%.2f", 100-$exon-$intron);
    close IN;
    print "$samp\t$totalM\t$exon\t$intron\t$intergene\n";
}
