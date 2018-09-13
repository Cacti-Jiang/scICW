#!/usr/bin/perl
use warnings;
use strict;
die "perl $0 [merged.vcf]" if (@ARGV!=1);
my ($vcf)=@ARGV;

if ($vcf =~/\.gz$/){
	open IN,"zcat $vcf |" or die $!;
}else{
	open IN,"$vcf" or die $!;
}

##DP>=20, ALT>=3
while(<IN>){
	next if /^\#\#/;
	chomp;
	my @F=split;
	if (/^\#CH/){
		for (9..$#F){
			my $samp = $1 if($F[$_] =~ /(\S+)\.variant/);
			print "$samp\t";
		}
		print "\n";
		next;
	}
	my @G=split /\:/,$F[8];
	next if(!($#G==4));
	my $sign = $1 if ($F[0] =~ /chr(\d+)/);
	my $pos=$sign.".$F[1]";
#	my $gt=0;
	my $count_AL = 0;
	print "$pos\t";
	for my $n(9..$#F){
		my $gt=0;
		@G = split /\:/,$F[$n];
		my @AL=split(",",$G[1]);
		if ($F[$n] =~ /\.\/\./){
			$gt="NA";
		}elsif($G[2] >= 20){
			$gt=$AL[0];
			if ($G[0] =~ /0\// && $AL[1] >=3){
				$gt = 1;
			}else{
				for(1..$#AL){
					$count_AL = 1 if($AL[$_]>=3);
				}		
				$gt=2 if($count_AL == 1);
			}
		}
		print "$gt\t";
	}
	print "\n";
}
close IN;
