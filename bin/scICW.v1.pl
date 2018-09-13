#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use FindBin qw /$RealBin/;
use Cwd qw/abs_path/;

my ( $Config, $samplelist, $template, $help);
my ( $outdir, $steps) = ( "./", "1,2,3,4,5,6,7,8,9,10");
GetOptions(
		"conf:s"     => \$Config,
		"list:s"     => \$samplelist,
		"outdir:s"   => \$outdir,
		"steps:s"    => \$steps,
		"template"   => \$template,
		"help|?"     => \$help,
);

if ($template){
	Generate_template();
	print "\nA template configure file was generated.\n";
	print "Please modify this file and rerun the script with configure file.\n\n";
	exit;
}

if (!$Config || !$samplelist || !$steps || $help) {
  die <<USAGE;
================================================================================================================================
Description:    RNA-seq Pipeline
Usage:          perl $0 [options]
Version:        1.0
Options:
	Basic:
			* --conf      Configure file
			* --list      Sample list, format 'SampleID  fqPath'
			                (If PE reads, use comma to separat two fq files)
			  --outdir    Output directory, default [./]
	Step:
			  --steps     Step choosen with comma separating, default [1,2,3,4,5,6,7,8,9,10]
			                1. Inspect raw data
			                2. Filter adapter contained, high Ns, and low-quality sequences
			                3. Align clean reads to reference
			                4. Quantify expression level based on ReadCount/FPKM/TPM
			                5. Merge RC/FPKM/TPM of samples to matrixes
					6. Clustering cells by Seurat
					7. Detection of alternative splicing events
					8. SNP calling
					9. Detection of fusion genes
					10. Detection of RNA-editing events

	Func:
			  --template  Generate default configure file
			  --help      Print this help information

E.g.:
			perl $0 --template
			perl $0 --conf config.txt --list sample.list --outdir ./
================================================================================================================================
USAGE
}

### Read Parameter Configure
$outdir = abs_path $outdir;
my %steps = map {$_, 1} split (",",$steps);
my %hash_Config;
open CONF, $Config or die "File not exists: $Config\n";
while (<CONF>){
        next if (/^\s*$/ || /^\s*#/ || !/=/);
        s/^([^#]+)#.*$/$1/;
        my ($key, $value) = split /=/,$_,2;
        $key =~ s/^\s*(.*?)\s*$/$1/;
        $value =~ s/^\s*(.*?)\s*$/$1/;
        next unless $value;
        $hash_Config{$key} = $value;
}
close CONF;

### Read Sample List
my ($SampleNumber,@SampleName,@inFastq,@seqtype) = (0);
open SAMP, $samplelist  or die $!;
while(<SAMP>){
	$SampleNumber++;
	my @F = split /\s+/;
	die "The input lst must be checked\n" if(scalar @F != 2);
	push @SampleName,$F[0];
	push @inFastq,$F[1];
	$F[1]=~/,/ ? (push @seqtype,"pe") : (push @seqtype,"se") ;
	for(split /,/,$F[1]){
		die "The input files not exist: $_\n" unless(-e $_);
	}
}
close SAMP;


######################################Main Program######################################
print "=========================================Program Start~!!=========================================\n";
system("date");
system("mkdir -p $outdir");
system("mkdir -p $outdir/shell");
system("mkdir -p $outdir/Stat");

### Step1. Inspect raw data
if( exists $steps{1} ){
	system("mkdir -p $outdir/01.RawData");
	if ($hash_Config{Step1_FastQC} eq "True"){
		open TMP1,">$outdir/shell/01-1.RawFastQC.sh" or die $!;
	
		for(0..$#SampleName){
			my $dir = "$outdir/01.RawData/$SampleName[$_]";
			system("mkdir -p $dir");
			if ($seqtype[$_] eq "se"){
				my $fq = abs_path $inFastq[$_];
				system("ln -sf $fq $dir/$SampleName[$_].fq.gz");
				$fq = "$dir/$SampleName[$_].fq.gz";
				$inFastq[$_] = $fq;
				if ($hash_Config{Step1_FastQC} eq "True"){
					print TMP1 "$hash_Config{FastQC} --outdir $dir --threads $hash_Config{Step1_threads} --quiet $fq \n";
				}
			}else{
				my @fqs = split ",",$inFastq[$_];
				$fqs[0] = abs_path $fqs[0];
				$fqs[1] = abs_path $fqs[1];
				system("ln -sf $fqs[0] $dir/$SampleName[$_]_1.fq.gz");
				system("ln -sf $fqs[1] $dir/$SampleName[$_]_2.fq.gz");
				$fqs[0] = "$dir/$SampleName[$_]_1.fq.gz";
				$fqs[1] = "$dir/$SampleName[$_]_2.fq.gz";
				$inFastq[$_] = "$fqs[0],$fqs[1]";
				if ($hash_Config{Step1_FastQC} eq "True"){
					print TMP1 "$hash_Config{FastQC} --outdir $dir --threads $hash_Config{Step1_threads} --quiet $fqs[0] $fqs[1] \n";
				}
			}
		}
		close TMP1;
	}
	
	# MultiQC
	if ($hash_Config{Step1_MultiQC} eq "True"){
		system("mkdir -p $outdir/Stat/MultiQC/01.RawData");
		open TMP2,">$outdir/shell/01-2.RawMultiQC.sh" or die $!;
		print TMP2 "export PYTHONPATH=\"$hash_Config{python_env}\" && $hash_Config{MultiQC} -o $outdir/Stat/MultiQC/01.RawData --quiet $outdir/01.RawData \n";
		close TMP2;
	}
}

### Step2. Filter reads
if( exists $steps{2} ){
	system("mkdir -p $outdir/02.CleanData/");
	open TMP1,">$outdir/shell/02-1.FilterReads.sh" or die $!;
	if ($hash_Config{Step2_FastQC} eq "True"){
		open TMP2,">$outdir/shell/02-2.CleanFastQC.sh" or die $!;
	}
	for(0..$#SampleName){
		my $dir = "$outdir/02.CleanData/$SampleName[$_]";
		system("mkdir -p $dir");
		my $pars = "";
		($hash_Config{phred} eq "phred33") ? ($pars = "--quality-base=33") : ($pars = "--quality-base=64");
		if ($hash_Config{Clean_Model} eq "Mask" ){
			$pars = "--mask-adapter $pars --max-n=$hash_Config{MaxN}";
		}elsif($hash_Config{Clean_Model} eq "Discard"){
			$pars = "$pars --max-n=$hash_Config{MaxN} --discard";
		}else{
			$pars = "--quality-cutoff $hash_Config{Quality_cutoff},$hash_Config{Quality_cutoff} $pars -m $hash_Config{MinLen} --max-n=$hash_Config{MaxN}";
		}
		if($seqtype[$_] eq "se"){
			# cutadapt with python2 version only used 1 thread
			$pars = "-j 1 -a $hash_Config{R1_3ad} $pars";
			print TMP1 "export PYTHONPATH=\"$hash_Config{python_env}\" && $hash_Config{Cutadapt} $pars -o $dir/$SampleName[$_].clean.fq.gz $inFastq[$_] >$dir/$SampleName[$_].log \n";
            $inFastq[$_] = "$dir/$SampleName[$_].clean.fq.gz";
			if ($hash_Config{Step2_FastQC} eq "True"){
				print TMP2 "$hash_Config{FastQC} --outdir $dir --threads $hash_Config{Step2_threads} --quiet $dir/$SampleName[$_].clean.fq.gz \n";
			}
		}
		else{
			my @fqs = split ",",$inFastq[$_];
			# cutadapt with python2 version only used 1 thread
			$pars = "-j 1 -a $hash_Config{R1_3ad} -A $hash_Config{R2_3ad} $pars";
			print TMP1 "export PYTHONPATH=\"$hash_Config{python_env}\" && $hash_Config{Cutadapt} $pars -o $dir/$SampleName[$_].clean_1.fq.gz -p $dir/$SampleName[$_].clean_2.fq.gz $fqs[0] $fqs[1] >$dir/$SampleName[$_].log \n";
            $inFastq[$_] = "$dir/$SampleName[$_].clean_1.fq.gz,$dir/$SampleName[$_].clean_2.fq.gz";
			if ($hash_Config{Step2_FastQC} eq "True"){
				print TMP2 "$hash_Config{FastQC} --outdir $dir --threads $hash_Config{Step2_threads} --quiet $dir/$SampleName[$_].clean_1.fq.gz $dir/$SampleName[$_].clean_2.fq.gz \n";
			}
		}
	}
	close TMP1;
	# FastQC
	if ($hash_Config{Step2_FastQC} eq "True"){
		close TMP2;
	}
	# MultiQC
	if ($hash_Config{Step2_MultiQC} eq "True"){
		system("mkdir -p $outdir/Stat/MultiQC/02.CleanData");
		open TMP3,">$outdir/shell/02-3.CleanMultiQC.sh" or die $!;
		print TMP3 "export PYTHONPATH=\"$hash_Config{python_env}\" && $hash_Config{MultiQC} -o $outdir/Stat/MultiQC/02.CleanData --quiet $outdir/02.CleanData \n";
		close TMP3;
	}
}

### Step3. Align clean reads
if( exists $steps{3} ){
	system("mkdir -p $outdir/03.Align/");
	open TMP1, ">$outdir/shell/03-1.Align2genome.sh" or die $!;
	if ($hash_Config{RunQualiMap} eq "True"){
		open TMP2, ">$outdir/shell/03-2.AlignQualimap.sh";
	}
	for(0..$#SampleName){
		my $dir = "$outdir/03.Align/$SampleName[$_]";
		system("mkdir -p $dir");
		## STAR
		my $par1 = "--runMode alignReads --runThreadN $hash_Config{Step3_threads} --genomeDir $hash_Config{STAR_index}";
		my $par2 = "--readFilesCommand zcat --outFileNamePrefix $dir/$SampleName[$_] --outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif";
		$par2 = "$par2 --outSAMunmapped Within" if ($hash_Config{outSAMunmapped} eq "True");
		$par2 = "$par2 --outFilterType BySJout --outFilterMultimapNmax $hash_Config{MultimapNmax} --quantMode TranscriptomeSAM";
		my $par3 = "--twopassMode Basic --outReadsUnmapped None --chimSegmentMin 12 --chimJunctionOverhangMin 12 --alignSJDBoverhangMin 10 --alignMatesGapMax 100000 --alignIntronMax 100000 --chimSegmentReadGapMax 3 --alignSJstitchMismatchNmax 5 -1 5 5";
		if($seqtype[$_] eq "se"){
			print TMP1 "$hash_Config{STAR} $par1 --readFilesIn $inFastq[$_] $par2 $par3 && $hash_Config{Samtools} index $dir/$SampleName[$_]Aligned.sortedByCoord.out.bam \n";
		}else{
			my @fqs = split ",",$inFastq[$_];
			print TMP1 "$hash_Config{STAR} $par1 --readFilesIn $fqs[0] $fqs[1] $par2 $par3 && $hash_Config{Samtools} index $dir/$SampleName[$_]Aligned.sortedByCoord.out.bam \n";
		}
		my $bam_prefix = "$dir/$SampleName[$_]Aligned.sortedByCoord.out";
		# QualiMap
		if ($hash_Config{RunQualiMap} eq "True"){
			my $command = "$hash_Config{QualiMap} rnaseq --java-mem-size=$hash_Config{QM_vf} -a proportional -bam $bam_prefix.bam -gtf $hash_Config{QM_GTF} -outdir $dir/QualiMap";
			$command = "$command --paired" if $seqtype[$_] eq "pe";
			print TMP2 "$command \n";
		}
	# MultiQC
		if ($hash_Config{Step3_MultiQC} eq "True"){
			system("mkdir -p $outdir/Stat/MultiQC/03.Align");
			open TMP3,">$outdir/shell/03-3.AlignMultiQC.sh" or die $!;
			print TMP3 "export PYTHONPATH=\"$hash_Config{python_env}\" && $hash_Config{MultiQC} -o $outdir/Stat/MultiQC/03.Align --quiet $outdir/03.Align \n";
			close TMP3;
		}
	}
	close TMP1;
	close TMP2;		
}

### Step4. Calculate expression values
if( exists $steps{4} ) {
	system ("mkdir -p $outdir/04.QuantExpr/");	
	open TMP,">$outdir/shell/04.CalExpr.sh" or die $!;
	for (0..$#SampleName){
		my $dir1 = "$outdir/04.QuantExpr/$SampleName[$_]";
		system("mkdir -p $dir1");
		my $bam1 = "$outdir/03.Align/$SampleName[$_]/$SampleName[$_]Aligned.toTranscriptome.out.bam";
		my $pars = "--alignments --no-bam-output -p $hash_Config{Step4_threads}";
		$pars = "$pars --paired-end" if ($seqtype[$_] eq "pe");
		$pars = "$pars --single-cell-prior" if ($hash_Config{Single_cell_prior} eq "True");
		print TMP "$hash_Config{RSEM} $pars $bam1 $hash_Config{RSEM_index} $dir1/$SampleName[$_] \n";
	}
	close TMP;
}


### Step5. Generate expression matrix
if( exists $steps{5} ) {
	my $dir1 = "$outdir/05.ExprMat";
	system ("mkdir -p $dir1");
	my $Lines1 = "";
	# RSEM results
	open L1, ">$dir1/gene_result.list" or die $!;
	for (0..$#SampleName){
		print L1 "$SampleName[$_]\t$outdir/04.QuantExpr/$SampleName[$_]/$SampleName[$_].genes.results\n";
	}
	close L1;
	open TMP1, ">$outdir/shell/05.ExprMat.sh" or die $!;
	print TMP1 "$hash_Config{Rscript} $hash_Config{RSEM_Merge2Matrix} $dir1/gene_result.list $dir1/ReadCount.gene.matrix.txt $dir1/TPM.gene.matrix.txt $dir1/FPKM.gene.matrix.txt \n";
	close TMP1;
}

### Step6. Clustering by Seurat
if( exists $steps{6} ) {
	system ("mkdir -p $outdir/06.Clustering/");
	open TMP,">$outdir/shell/06.Cluster.sh" or die $!;
	print TMP "$hash_Config{Rscript} $hash_Config{ClusterBySeurat} $outdir/05.ExprMat/TPM.gene.matrix.txt -o $outdir/06.Clustering/AllCells\n";
	close TMP;
}

### Step7. AS detection
if( exists $steps{7} ) {
	system ("mkdir -p $outdir/07.Outrigger/");
	open TMP,">$outdir/shell/07.Outrigger.sh" or die $!;
	print TMP "source activate outrigger-env\n";
	print TMP "outrigger index --low-memory --n-jobs $hash_Config{Step7_threads} -resume --gtf $hash_Config{QM_GTF} --sj-out-tab $outdir/03.Align/*/*SJ.out.tab -o $outdir/07.Outrigger\n";
	print TMP "outrigger validate --low-memory -f $hash_Config{Ref_fasta} -g $hash_Config{Ref_len} -o $outdir/07.Outrigger\n";
	print TMP "outrigger psi --low-memory --n-jobs $hash_Config{Step7_threads} -o $outdir/07.Outrigger\n";
	close TMP;
}

### Step8. SNP calling
if( exists $steps{8} ) {
	system ("mkdir -p $outdir/08.SNP/");
	system("mkdir -p $outdir/08.SNP/java_tmp");
 	system("mkdir -p $outdir/08.SNP/08-1.preBAM");
 	system("mkdir -p $outdir/08.SNP/08-2.SplitTrim");
 	system("mkdir -p $outdir/08.SNP/08-3.BQSR");
 	system("mkdir -p $outdir/08.SNP/08-4.VarCall");
 	system("mkdir -p $outdir/08.SNP/08-5.VarFilter");
 	system("mkdir -p $outdir/08.SNP/08-6.MergeVcfs");
 	system("mkdir -p $outdir/08.SNP/08-7.Sifit");

 	open TMP1, ">$outdir/shell/08-1.preBAM.sh";
	open TMP2, ">$outdir/shell/08-2.SplitTrim.sh";
	open TMP3, ">$outdir/shell/08-3.BQSR.sh";
	open TMP4, ">$outdir/shell/08-4.VarCall.sh";
	open TMP5, ">$outdir/shell/08-5.VarFilter.sh";
	open L1, ">$outdir/08.SNP/vcf.list";
	for(0..$#SampleName){
		#01.preBAM
 	   	system("mkdir -p $outdir/08.SNP/08-1.preBAM/$SampleName[$_]");
    		#add @RG tag
    		print TMP1 "$hash_Config{java} -Djava:q.io.tmpdir=$outdir/08.SNP/java_tmp -Xmx5g -jar $hash_Config{Picard} AddOrReplaceReadGroups I=$outdir/03.Align/$SampleName[$_]/$SampleName[$_]Aligned.sortedByCoord.out.bam O=$outdir/08.SNP/08-1.preBAM/$SampleName[$_]/$SampleName[$_].rg.bam SO=coordinate RGID=$SampleName[$_] RGLB=$SampleName[$_] RGPL=$hash_Config{Platform} RGPU=$hash_Config{Platform} RGSM=$SampleName[$_]\n";
    		#deDup
    		print TMP1 "$hash_Config{java} -Djava:q.io.tmpdir=$outdir/08.SNP/java_tmp -Xmx5g -jar $hash_Config{Picard} MarkDuplicates I=$outdir/08.SNP/08-1.preBAM/$SampleName[$_]/$SampleName[$_].rg.bam O=$outdir/08.SNP/08-1.preBAM/$SampleName[$_]/$SampleName[$_].rg.dedup.bam CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT M=$outdir/08.SNP/08-1.preBAM/$SampleName[$_]/$SampleName[$_].rg.dedup.metrics\n";

    		#02.SplitTrim
    		system("mkdir -p $outdir/08.SNP/08-2.SplitTrim/$SampleName[$_]");
    		#20180523, BUG in quality score error, please add "--fix_misencoded_quality_scores" option
    		print TMP2 "$hash_Config{java} -Djava:q.io.tmpdir=$outdir/08.SNP/java_tmp -Xmx16g -jar $hash_Config{GATK} -T SplitNCigarReads -R $hash_Config{Ref_fasta} -I $outdir/08.SNP/08-1.preBAM/$SampleName[$_]/$SampleName[$_].rg.dedup.bam -o $outdir/08.SNP/08-2.SplitTrim/$SampleName[$_]/$SampleName[$_].split.bam -rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 -U ALLOW_N_CIGAR_READS --fix_misencoded_quality_scores\n";

    		#03.BQSR
	    	#generate bqsr.table
    		system("mkdir -p $outdir/08.SNP/08-3.BQSR/$SampleName[$_]");
    		print TMP3 "$hash_Config{java} -Djava:q.io.tmpdir=$outdir/08.SNP/java_tmp -Xmx16g -jar $hash_Config{GATK} -T BaseRecalibrator -R $hash_Config{Ref_fasta} -I $outdir/08.SNP/08-2.SplitTrim/$SampleName[$_]/$SampleName[$_].split.bam -knownSites $hash_Config{dbSNP} -o $outdir/08.SNP/08-3.BQSR/$SampleName[$_]/$SampleName[$_].bqsr.table\n";
    		#PrintReads
    		print TMP3 "$hash_Config{java} -Djava:q.io.tmpdir=$outdir/08.SNP/java_tmp -Xmx16g -jar $hash_Config{GATK} -T PrintReads -R $hash_Config{Ref_fasta} -I $outdir/08.SNP/08-2.SplitTrim/$SampleName[$_]/$SampleName[$_].split.bam -BQSR $outdir/08.SNP/08-3.BQSR/$SampleName[$_]/$SampleName[$_].bqsr.table -o $outdir/03.BQSR/$SampleName[$_]/$SampleName[$_].bqsr.bam\n";
    		#new ApplyBQSR (gatk4.0 will replace 'PrintReads')
    		#print TMP3 "$hash_Config{java} -Djava:q.io.tmpdir=$outdir/08.SNP/java_tmp -Xmx16g -jar $hash_Config{GATK} -T ApplyBQSR -R $hash_Config{Ref_fasta} -I $outdir/08.SNP/08-2.SplitTrim/$SampleName[$_]/$SampleName[$_].split.bam --bqsr-recal-file $outdir/03.BQSR/$SampleName[$_]/$SampleName[$_].bqsr.table -O $outdir/08.SNP/08-3.BQSR/$SampleName[$_]/$SampleName[$_].bqsr.bam\n";

    		#04.VarCall
	    	system("mkdir -p $outdir/08.SNP/08-4.VarCall/$SampleName[$_]");
    		print TMP4 "$hash_Config{java} -Djava:q.io.tmpdir=$outdir/08.SNP/java_tmp -Xmx16g -jar $hash_Config{GATK} -T HaplotypeCaller -R $hash_Config{Ref_fasta} -I $outdir/08.SNP/08-3.BQSR/$SampleName[$_]/$SampleName[$_].bqsr.bam -dontUseSoftClippedBases -stand_call_conf 20.0 -o $outdir/08.SNP/08-4.VarCall/$SampleName[$_]/$SampleName[$_].vcf\n";

    		#05.VarFilter
    		system("mkdir -p $outdir/08.SNP/08-5.VarFilter/$SampleName[$_]");
    		print TMP5 "$hash_Config{java} -Djava:q.io.tmpdir=$outdir/08.SNP/java_tmp -Xmx16g -jar $hash_Config{GATK} -T VariantFiltration -R $hash_Config{Ref_fasta} -V $outdir/08.SNP/08-4.VarCall/$SampleName[$_]/$SampleName[$_].vcf -window 35 -cluster 3 -filterName FS -filter \"FS > 30.0\" -filterName QD -filter \"QD < 2.0\" -o $outdir/05.VarFilter/$SampleName[$_]/$SampleName[$_].filter.vcf\n";
		
		##vcf list
		print L1 "$outdir/05.VarFilter/$SampleName[$_]/$SampleName[$_].filter.vcf\n";
	}
	close TMP1;
	close TMP2;
	close TMP3;
	close TMP4;
	close TMP5;
	close L1;
	
	#Merge vcfs(SNPs only)
	open TMP6, ">$outdir/shell/08-6.MergeVcfs.sh";
	print TMP6 "$hash_Config{java} -jar $hash_Config{GATK} -T CombineVariants -R $hash_Config{Ref_fasta} --variant $outdir/08.SNP/vcf.list -o $outdir/08.SNP/08-6.MergeVcfs/Merged.vcf -genotypeMergeOptions UNIQUIFY\n";
	print TMP6 "$hash_Config{Vcftools} --vcf $outdir/08.SNP/08-6.MergeVcfs/Merged.vcf --remove-indels --out $outdir/08.SNP/08-6.MergeVcfs/Merged.snps.vcf --recode --recode-INFO-all\n";
	close TMP6;

	open TMP7, ">$outdir/shell/08-7.Sifit.sh";
	print TMP7 "$hash_Config{Rscript} $hash_Config{Vcf2Mat} $outdir/08.SNP/08-6.MergeVcfs/Merged.snps.vcf > $outdir/08.SNP/08-7.Sifit/Merged.snps.raw.mat\n";
	print TMP7 "head -n 1 $outdir/08.SNP/08-7.Sifit/Merged.snps.raw.mat | perl -wlane 'chomp; my \$samp=join(\" \", \@F); print \$samp;' > $outdir/08.SNP/08-7.Sifit/cellname.txt\n";
	print TMP7 "$hash_Config{Rscript} $hash_Config{BiMat} $outdir/08.SNP/08-7.Sifit/Merged.snps.raw.mat 0.3 $outdir/08.SNP/08-7.Sifit/Merged.snps.0308.bi\n";
	print TMP7 "$hash_Config{java} -jar $hash_Config{Sifit} -m 49 -n 610 -fp 0.002 -fn 0.2 -iter 10000 -df 0 -ipMat $outdir/08.SNP/08-7.Sifit/Merged.snps.0308.bi.mat -cellNames $outdir/08.SNP/08-7.Sifit/cellname.txt > $outdir/08.SNP/08-7.Sifit/sifit.log\n";
	close TMP7;

}

### Step9. Detection of fusion genes
if(exists $steps{9}) {
	system ("mkdir -p $outdir/09.StarFusion");
	open TMP, ">$outdir/shell/09.StarFusion.sh";
	for(0..$#SampleName){
		system("mkdir -p $outdir/09.StarFusion/$SampleName[$_]");
		print TMP "$hash_Config{StarFusion} --genome_lib_dir $hash_Config{CTAT_lib} -J $outdir/03.Align/$SampleName[$_]/$SampleName[$_]Chimeric.out.junction --output_dir $outdir/09.StarFusion/$SampleName[$_] --CPU $hash_Config{Step9_threads} --FusionInspector inspect\n";
	}
	close TMP;
}

### Step10. Detection of RNA-editing events
if(exists $steps{10}) {
	system ("mkdir -p $outdir/10.Editing");
	open TMP, ">$outdir/shell/10.Editing.sh";
	for(0..$#SampleName){
		system("mkdir -p $outdir/10.Editing/$SampleName[$_]");
		print TMP "$hash_Config{perl} $hash_Config{red_ML} --rnabam $outdir/03.Align/$SampleName[$_]/$SampleName[$_]Aligned.sortedByCoord.out.bam --reference $hash_Config{Ref_fasta} --dbsnp $hash_Config{dbSNP} --simpleRepeat $hash_Config{SimpleRepeat} --alu $hash_Config{ALU} --outdir $outdir/10.Editing/$SampleName[$_]\n";
	}
	close TMP;
}

### Final. Statistics
open STAT, ">$outdir/shell/Statistics.sh" or die $!;
my $StatDir = "$outdir/Stat";
my @plot_density_cols = ();
my @plot_bar_cols = ();
if(exists $steps{2}){
	print STAT "$hash_Config{perl} $hash_Config{Stat_Cutadapt} $outdir/02.CleanData/*/*.log >$StatDir/01.Reads_stat.txt \n";
	push @plot_density_cols, "RawReads/Pairs,CleanReads/Pairs,Cp(%)";
	push @plot_bar_cols, "CLEAN,Cp(%)";
	# generate CleanFq list
	open LFQ, ">$outdir/Stat/CleanFq.list" or die $!;
	print LFQ "SampleID\tCleanFq\n";
	for (0..$#SampleName){
		print LFQ "$SampleName[$_]\t$inFastq[$_]\n";
	}
	close LFQ;
}
if(exists $steps{3}){
	push @plot_density_cols, "totalMapped,totalMapped(%),uniqMapped(%)";
	push @plot_bar_cols, "totalMapped,MAPPING";
	print STAT "$hash_Config{perl} $hash_Config{Stat_STAR} $outdir/03.Align/*/*Log.final.out >$StatDir/02-1.Align_stat.txt \n";
	# generate Bam list
	open LBAM, ">$outdir/Stat/Bam.list" or die $!;
	print LBAM "SampleID\tBam\n";
	for (0..$#SampleName){
		print LBAM "$SampleName[$_]\t$outdir/03.Align/$SampleName[$_]/$SampleName[$_]Aligned.sortedByCoord.out.bam\n";
	}
	close LBAM;
	if ($hash_Config{RunQualiMap} eq "True"){
		# Reads Distribution
		print STAT "$hash_Config{perl} $hash_Config{Stat_QualiMap} $outdir/03.Align/*/QualiMap/rnaseq_qc_results.txt >$StatDir/02-2.Distribution_stat.txt \n";
		push @plot_density_cols, "Exonic(%),Intronic(%),Intergenic(%)";
		push @plot_bar_cols, "REGION";
	}
}
if (exists $steps{5}){
	print STAT "$hash_Config{Rscript} $hash_Config{Stat_Expr} $outdir/05.ExprMat/noSPIKEIN/TPM.gene.matrix.txt $StatDir/03.Expr_stat.txt \n";
	push @plot_density_cols, "TPM0,TPM1";
	push @plot_bar_cols, "TPM";
}
if (exists $steps{7}){
	print STAT "awk -F ',' '{print \$(NF-2)\"\\t\"\$(NF-1)}' psi/outrigger_summary.csv | sed 1d | awk '\$1!=\"NA\"' |cut -f 2  |sort |uniq -c | awk '{print \$2\"\\t\"\$1}' >$StatDir/07.AS_stat.txt\n";
}
my $plot_density_cols = join(",", @plot_density_cols);
my $plot_bar_cols = join(",", @plot_bar_cols);
print STAT "$hash_Config{Rscript} $hash_Config{Stat_Merge} $outdir/Stat/*_stat.txt $outdir/Stat/*.list >$outdir/Stat/Merge_info.xls \n";
close STAT;

system("date");
print "=========================================Program finished=========================================\n";
###-----------END-----------


sub Generate_template{
	open OUTCONF, ">./default_config.txt" or die $!;
	print OUTCONF "#### Configure file
### Common tools
java = /path/to/bin/java
perl = /path/to/ActivePerl-5.18.4/bin/perl
python = /path/to/Anaconda2-4.3.0/bin/python
python_env = /path/to/Anaconda2-4.3.0/lib/python2.7     # seperated by colon
Rscript = /path/to/R-3.4.3/bin/Rscript
Samtools = /path/to/samtools-1.5/bin/samtools
Vcftools = /path/to/vcftools-0.1.15/bin/vcftools
FastQC = /path/to/FastQC/fastqc
MultiQC = /path/to/multiqc
QualiMap = /path/to/qualimap_v2.2.1/qualimap


### Sequencing information
phred = phred33  # phred33 or phred64
read_length = 100

### Step1: Inspect raw data
Step1_FastQC = True   # Run FastQC for raw reads
Step1_MultiQC = True  # Run MultiQC to summary results of FastQC
Step1_threads = 1

### Step2: Filter reads
Cutadapt = /path/to/Anaconda2-4.3.0/bin/cutadapt
Clean_Model = Mask  # Mask, Discard or Trim
R1_3ad = CTGTCTCTTATACACATCTCCGAGCCCACGAGAC   # 3' end contamination for Reads1, default is Nextera adapter
R2_3ad = CTGTCTCTTATACACATCTGACGCTGCCGACGA    # 3' end contamination for Reads2, default is Nextera adapter; SE reads will be ignored
MaxN = 0.2  # Maximum percentage of N bases
Quality_cutoff = 5  # Cutoff for trimming low-quality ends; Only used for TRIM clean
MinLen = 80 # Minimum length for cleaned reads; Only used for TRIM clean
Step2_FastQC = True   # Run FastQC for clean reads
Step2_MultiQC = True  # Run MultiQC to summary reports
Step2_threads = 1

### Step3: Align clean reads
RunQualiMap = True  # Run QualiMap to calculate reads mapping distribution
QM_vf = 20G  # Specific memory setting for QualiMap
QM_GTF = /path/to/Ref/gencode.v27.primary_assembly.annotation.gtf #Gencode_human.release27
Step3_MultiQC = True  # Run MultiQC to summary reports
Step3_threads = 4
## STAR MODEL
STAR = /path/to/STAR-2.6.0c/bin/Linux_x86_64/STAR
STAR_index = /path/to/Ref/STAR_sj100/ 
MultimapNmax = 10  # Maximum multi-map number
outSAMunmapped = False  # If 'True', also output unmapped reads in BAM file

### Step4: Quantify expression level
## RSEM MODEL
RSEM = /path/to/RSEM-1.3.0/rsem-calculate-expression
RSEM_index = /path/to/Ref/RSEM_STAR_sj100/gencode.v27
Single_cell_prior = True  # Used SingleCell parameter in RSEM
Step4_threads = 2

### Step5: Generate expression matrix
RSEM_Merge2Matrix = /path/to/bin/Merge2Matrix_from_RSEM.R   # For STAR->RSEM model

### Step6: Clustering
ClusterBySeurat = /path/to/bin/ClusterBySeurat.R

### Step7: Outrigger
Ref_fasta = /path/to/Ref/GRCh38.primary_assembly.genome.fa
Ref_len = /path/to/Ref/GRCh38.primary_assembly.genome.fa.len
Step7_threads = 8

### Step8: SNP calling
Picard = /path/to/picard_2.10.10_data/picard.jar
GATK = /path/to/GenomeAnalysisTK.jar
Platform = COMPLETE # choose 'ILLUMINA' or 'COMPLETE'
dbSNP = /path/to/Ref/dbsnp_150.hg38.vcf.gz
Vcf2Mat = /path/to/bin/Mvcf2RawMat.nochr.pl
BiMat = /path/to/bin/RawMat2Mat0.R
Sifit = /path/to/SiFit.jar

### Step9: Detection of fusion genes
StarFusion = /path/to/STAR-Fusion-v1.4.0/STAR-Fusion
CTAT_lib = /path/to/Ref/ctat_genome_lib_build_dir
Step9_threads = 4

### Step10: Detection of RNA-editing events
red_ML = /path/to/red_ML.v2.pl
SimpleRepeat = /path/to/Ref/simpleRepeat.bed
ALU = /path/to/Ref/hg38.alu.bed

### Statistics
Stat_Cutadapt = /path/to/bin/Stat_Cutadapt.pl
Stat_STAR = /path/to/bin/Stat_STAR.pl  # For STAR->RSEM model
Stat_QualiMap = /path/to/bin/Stat_QualiMap_RNAseq.pl
Stat_Expr = /path/to/bin/Stat_ExprMat.R
Stat_Merge = /path/to/bin/Stat_Merge.R
";
}

######################################Program End######################################
