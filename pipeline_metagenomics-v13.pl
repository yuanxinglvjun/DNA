#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use Config::IniFiles;
use Cwd qw/getcwd abs_path/;
use FindBin qw($Bin $Script $RealBin);
use File::Basename qw(basename dirname);
use PBS::Queue;

my $BEGIN_TIME=time();
my $script   = "$Bin/script";
my $software = "$Bin/software";
my $database = "$Bin/database";

my $workdir  = getcwd;

my ($help, $fq_list, $file_meta, $file_comp, $outdir, $species, $show, $refer, $insertsize);
my ($single_sample);
my ($max_cpus, $max_jobs, $queue_name, $job_prefix, $job_local, $queue_rerun);
my %samples;
my %group_info;
my @sample_names;

sub usage{
	my $program = basename($0);
    die "
Usage:	$program [arguments]
Arguments:
  Input Options:
    -input          <string>    clean data fq file list.[require]
    -outdir         <string>    output dir, default is 'out'.
    -meta           <string>    meta information for grouping comparison. [ optional ]
    -comp           <string>    comparison configure file. [ optional ]
    -refer          <string>    bwa index of host reference genome. [ optional ]

  Step Options:
    -single_sample  <NA>        only perform analysis for single sample, do not analyze intergrated analysis.

  Anvanced Options:
    -insertsize     <number>    insert size of library, default is 400
    
    -max_cpus       <number>    max cpu number limitation, default is 200 [ optional ]
    -max_jobs       <number>    max job number limitation, default is 20  [ optional ]
    -queue_name     <string>    queue name, default is zh
    -local          <NA>        local model
    -rerun          <NA>        rerun model
###################
#PS1:<fq.lst> contain 3 colums, sample_id, read1.fq(gz), read2.fq(gz)
#		sample_1	/path/sample_1_R1.fq	/path/sample_1_R2.fq
#		sample_2	/path/sample_2_R1.fq	/path/sample_2_R2.fq
#		sample_3	/path/sample_3_R1.fq	/path/sample_3_R2.fq
#		......
#		For different lines, samples name must be unique !!!
#PS2:[meta] contain at leat 2 colums, 1st is [Sample], 2nd, 3rd, 4th ... are [Group Standard]
#		SampleID	Group_STD_1	Group_STD_2
#		sample_1	g1	w
#		sample_2	g1	d
#		sample_3	g2	d
#		sample_4	g2	w
#		sample_5	g3	d
#		sample_6	g3	w
#		For different lines, samples name must be unique !!!
#PS3:[comp] 1st is [Sample] or [Group Standard], 2nd, 3rd, 4th ... are group name
#		Sample	sample_1	sample_3
#		Group_STD_1	g1	g3
#		Group_STD_1	g1	g2	g3
#		Group_STD_2	w	d
#		For different lines, samples name must be unique !!!
###################
\n";
}


&GetOptions(
	"help!"        => \$help,
	"input:s"      => \$fq_list,
	"meta:s"       => \$file_meta,
	"comp:s"       => \$file_comp,
	"outdir:s"     => \$outdir,
	"refer:s"      => \$refer,
	"insertsize:i" => \$insertsize,

	
	"single_sample!"=> \$single_sample,
	
	"max_cpus:i"   => \$max_cpus,
	"max_jobs:i"   => \$max_jobs,
	"queue_name:s" => \$queue_name,
	"job_prefix:s" => \$job_prefix,
	"local!"       => \$job_local,
	"rerun!"       => \$queue_rerun,
);

&parsing_parameters();


my $queue = PBS::Queue->new(
	{
		'queue_name'  => $queue_name,
		'job_prefix'  => $job_prefix,
		'max_cpus'    => $max_cpus,
		'max_jobs'    => $max_jobs,
		'job_local'   => $job_local,
		'queue_rerun' => $queue_rerun,
	}
);

&parsing_sample_info;
&parsing_meta_and_comp_info;
&basic_analysis_per_sample;
#&rawdata_qc();
#&assembly();
#&gene_predict();
&unigene_set();
&gene_profile();
&gene_annotation();
&profile_collection();
&composition();
&comparison();
&difference();
&result_collection();
&html_report();

$queue->jointhreads();

print " all jobs done !!\n";

################################################################################

sub parsing_parameters{

	&usage if ($help);
	
	unless(defined( $fq_list )){
		print STDERR "Error:-input must be specified!\n";
		&usage;
	}
	
	$insertsize  ||= 400;
	$queue_name  ||= "yxsw";
	$job_prefix  ||= "metagenome";
	$queue_rerun ||= 0;
	$job_local   ||= 1;
	$max_cpus    ||= 200;
	$max_jobs    ||= 20;
	$outdir      ||= "out";

	$outdir = abs_path( $outdir );
	system("mkdir -p $outdir");
}

sub load_config{
	my ($configfile,$hash) = @_;
	
	unless($configfile){
		print STDERR "species cfg must be specified!\n";
		&usage;
	}
	unless(-e $configfile){
		print STDERR "species cfg file $configfile does not exist!";
		&usage;
	}
	my $cfg = Config::IniFiles->new( -file => $configfile );
	my @sec = $cfg->Sections();
	foreach my $s (@sec) {
		my @pa = $cfg->Parameters($s);
		foreach my $p (@pa) {
			$$hash{$s}{$p} = $cfg->val( $s, $p );
		}
	}
}

sub parsing_sample_info{
	open  FIN,$fq_list or die "can not open raw fq list file!";
	while(<FIN>){
		chomp;
		next if /^#/;
		my ($sample,$fq1,$fq2) = split /\t/,$_;
		if(exists $samples{$sample}){
			die "sample_name [$sample] is duplicate! please check input file[$fq_list]!";
		}
		$samples{$sample}{R1} = $fq1;
		$samples{$sample}{R2} = $fq2;
		push @sample_names,$sample;
	}
	close FIN;
}

sub parsing_meta_and_comp_info{
	system("mkdir -p $outdir/group_info");
	my %meta;
	my %mat_info;
	my %all_std;
	if(defined $file_meta && defined $file_comp){
		open FIN , "<$file_meta" or die "can not open meta information file [$file_meta]!";
		my $head = <FIN>; chomp($head);
		unless($head =~ "^#SampleID\t"){
			die "Meta file format error[header], please check!";
		}
		my ($a,@gs) = split /\t/,$head;
		while(<FIN>){
			chomp;
			my ($sample,@G) = split /\t/,$_;
			unless(exists $samples{$sample}){
				die "sample_name[$sample] in meta file does not existe in fasta file!";
			}
			if(scalar @G != scalar @gs){
				print $sample,"\t",scalar(@G),"\t",scalar(@gs),"\n";
				die "meta file format error[content], please check!";
			}
			for(my $i = 0; $i < scalar(@gs); $i++){
				if($G[$i] ne "null"){
					$meta{$gs[$i]}{$G[$i]}{$sample}++;
				}
			}
			$meta{"Sample"}{$sample}{$sample}++;
		}
		close FIN;
		
		open FIN, "<$file_comp" or die "can not open comp information file [$file_comp]!";
		while(<FIN>){
			chomp;
			next if /^#/;
			my ($name,$class,$tmp) = split /\t/,$_;
			unless(exists $meta{$class}){
				die "group std name [$class] does not exists in meta info file [$file_meta]!";
			}
			my @G =split /,/,$tmp;
			if( scalar(@G) == 1 ){
				die "the number of elements in comparison should not be 1!";
			}elsif( scalar(@G) == 0 ){
				@G = sort keys %{$meta{$class}};
				$name = $class."__all" unless $name;
			}else{
				$name = $class."__".join("_vs_",@G) unless $name;
			}
			if(exists $all_std{$name}){
				die "std[$name] was exists!!!";
			}
			$all_std{$name}++;
			foreach my $ele(@G){
				die "group name [$ele] does not exists in group std [$class] !!" unless exists $meta{$class}{$ele};
				foreach my $sam(sort keys %{$meta{$class}{$ele}}){
					$group_info{$name}{$ele}{$sam}++;
					$mat_info{$sam}{$name} = $ele;
				}
			}
		}
		close FIN;
		
		open  FOUT, ">$outdir/group_info/group_info.tsv";
		print FOUT  "#SampleID\t".join("\t",sort keys %all_std)."\n";
		foreach my $sam(sort keys %mat_info){
			print FOUT $sam;
			foreach my $name(sort keys %all_std){
				if(exists $mat_info{$sam}{$name}){
					print FOUT "\t",$mat_info{$sam}{$name};
				}else{
					print FOUT "\t";
				}
			}
			print FOUT "\n";
		}
		close FOUT;
		
		foreach my $class(sort keys %group_info){
			open  F1,">$outdir/group_info/$class.info";
			open  F2,">$outdir/group_info/$class.txt";
			print F1 "#id\t$class\n";
			foreach my $type(sort keys %{$group_info{$class}}){
				foreach my $sample(sort keys %{$group_info{$class}{$type}}){
					print F1 $sample,"\t",$type,"\n";
					print F2 $sample,"\t",$type,"\n";
				}
			}
			close F1;
			close F2;
		}
	}elsif(defined $file_meta && !(defined $file_comp)){
		open FIN , "<$file_meta" or die "can not open meta information file [$file_meta]!";
		# #SampleID\tBarcodeSequence\tLinkerPrimerSequence\tDescription
		my $head = <FIN>; chomp($head);
		unless($head =~ "^#SampleID\t"){
			die "Meta file format error, please check!";
		}
		my ($a,@gs) = split /\t/,$head;
		while(<FIN>){
			chomp;
			my ($sample,@G) = split /\t/,$_;
			unless(exists $samples{$sample}){
				die "sample_name[$sample] in meta file does not existe in fasta file!";
			}
			if(scalar @G != scalar @gs){
				die "meta file format error, please check!";
			}
			for(my $i = 0; $i < scalar(@gs); $i++){
				if($G[$i] ne ""){
					$group_info{$gs[$i]}{$G[$i]}{$sample}++;
					$mat_info{$sample}{$gs[$i]} = $G[$i];
					$all_std{$gs[$i]}++;
				}
			}
			$group_info{"Sample"}{$sample}{$sample}++;
			$mat_info{$sample}{"Sample"} = $sample;
			$all_std{"Sample"}++;
		}
		close FIN;
		close FOUT;
		
		open  FOUT, ">$outdir/group_info/group_info.tsv";
		print FOUT  "#SampleID\t".join("\t",sort keys %all_std)."\n";
		foreach my $sam(sort keys %mat_info){
			print FOUT $sam;
			foreach my $name(sort keys %all_std){
				if(exists $mat_info{$sam}{$name}){
					print FOUT "\t",$mat_info{$sam}{$name};
				}else{
					print FOUT "\t";
				}
			}
			print FOUT "\n";
		}
		close FOUT;
		
		foreach my $class(sort keys %group_info){
			open  F1,">$outdir/group_info/$class.info";
			open  F2,">$outdir/group_info/$class.txt";
			print F1 "#id\t$class\n";
			foreach my $type(sort keys %{$group_info{$class}}){
				foreach my $sample(sort keys %{$group_info{$class}{$type}}){
					print F1 $sample,"\t",$type,"\n";
					print F2 $sample,"\t",$type,"\n";
				}
			}
			close F1;
			close F2;
		}
	}elsif(defined $file_comp && !(defined $file_meta)){
		open FIN,"$fq_list";
		while(<FIN>){
			chomp;
			my ($sample) = split /\t/,$_;
			$meta{"Sample"}{$sample}{$sample}++;
			#$mat_info{$sample}{"Sample"} = $sample;
		}
		close FIN;
		
		open FIN, "<$file_comp" or die "can not open comp information file [$file_comp]!";
		while(<FIN>){
			chomp;
			my ($name,$class,$tmp) = split /\t/,$_;
			if($class ne "Sample"){
				die "class must be \"Sample\"!";
			}
			my @G =split /,/,$tmp;
			if( scalar(@G) == 1 ){
				die "the number of elements in comparison should not be 1!";
			}elsif( scalar(@G) == 0 ){
				@G = sort keys %{$meta{$class}};
				$name = $class."__all" unless $name;
			}else{
				$name = $class."__".join("_vs_",@G) unless $name;
			}
			$all_std{$name}++;
			foreach my $ele(@G){
				die "group name [$ele] does not exists in group std [$class] !!" unless exists $meta{$class}{$ele};
				foreach my $sam(sort keys %{$meta{$class}{$ele}}){
					$group_info{$name}{$ele}{$sam}++;
					$mat_info{$sam}{$name} = $ele;
				}
			}
		}
		close FIN;
		
		open  FOUT, ">$outdir/group_info/group_info.tsv";
		print FOUT  "#SampleID\t".join("\t",sort keys %all_std)."\n";
		foreach my $sam(sort keys %mat_info){
			print FOUT $sam;
			foreach my $name(sort keys %all_std){
				if(exists $mat_info{$sam}{$name}){
					print FOUT "\t",$mat_info{$sam}{$name};
				}else{
					print FOUT "\t";
				}
			}
			print FOUT "\n";
		}
		close FOUT;
		
		foreach my $class(sort keys %group_info){
			open  F1,">$outdir/group_info/$class.info";
			open  F2,">$outdir/group_info/$class.txt";
			print F1 "#id\t$class\n";
			foreach my $type(sort keys %{$group_info{$class}}){
				foreach my $sample(sort keys %{$group_info{$class}{$type}}){
					print F1 $sample,"\t",$type,"\n";
					print F2 $sample,"\t",$type,"\n";
				}
			}
			close F1;
			close F2;
		}
		
	}else{
		open FIN,"$fq_list";
		open F1,">$outdir/group_info/Sample_all.info";
		open F2,">$outdir/group_info/Sample_all.txt";
		print F1 "#id\tSample\n";
		while(<FIN>){
			chomp;
			my ($sample) = split /\t/,$_;
			$group_info{"Sample_all"}{$sample}{$sample}++;
			$mat_info{$sample}{"Sample_all"} = $sample;
			print F1 $sample,"\t",$sample,"\n";
			print F2 $sample,"\t",$sample,"\n";
		}
		close F1;
		close F2;
		close FIN;
		
		open FOUT,">$outdir/group_info/group_info.tsv";
		print FOUT "#SampleID\tSample_all\n";
		foreach my $sample(sort keys %{$group_info{"Sample_all"}}){
			print FOUT "$sample\t$sample\n";
		}
		close FOUT;
	}
	$file_meta = "$outdir/group_info/group_info.tsv";
}

sub basic_analysis_per_sample{
	return unless $single_sample;
	foreach my $sn(@sample_names){
		$queue->set_job_cpu(10);
		$queue->set_job_mem(10);
		$queue->set_job_name("样本 $sn 基本分析(指控、拼接、基因预测)");
		$queue->set_work_dir($workdir);
		
		$queue->addcommond("mkdir -p $outdir/fastq_QC/");
		foreach my $x(qw/R1 R2/){
			my @fqs = split /,/,$samples{$sn}{$x};
			$queue->addcommond("rm -f $outdir/fastq_QC/$sn.$x.fq");
			foreach my $fq(@fqs){
				if($fq =~ /.gz$/){
					$queue->addcommond("zcat $fq >> $outdir/fastq_QC/$sn.$x.fq");
				}else{
					$queue->addcommond("cat  $fq >> $outdir/fastq_QC/$sn.$x.fq");
				}
			}
		}
		$queue->addcommond("echo -e '$sn\\t$outdir/fastq_QC/$sn.R1.fq\\t$outdir/fastq_QC/$sn.R2.fq\\n' > $outdir/fastq_QC/$sn.rawfq.list");
		#$queue->addcommond("java -jar $script/FastqStat.jar -i $outdir/fastq_QC/$sn.rawfq.list > $outdir/fastq_QC/$sn.rawdata.stat");
		#$queue->addcommond("cut -f 1-3,11-14 $outdir/fastq_QC/$sn.rawdata.stat > $outdir/fastq_QC/$sn.rawdata.stat.xls");
		#$queue->addcommond("rm $outdir/fastq_QC/$sn.rawdata.stat");
		
		$queue->addcommond("cutadapt -q 20 -e 0.1 -n 1 -m 20 -a AGATCGGAAGAGC -g GCTCTTCCGATCT -o $outdir/fastq_QC/$sn.trim.R1.fq -p $outdir/fastq_QC/$sn.trim.R2.fq $outdir/fastq_QC/$sn.R1.fq $outdir/fastq_QC/$sn.R2.fq &> $outdir/fastq_QC/$sn.cutadapt.log");
		
		$queue->addcommond("echo -e '$sn\\t$outdir/fastq_QC/$sn.trim.R1.fq\\t$outdir/fastq_QC/$sn.trim.R2.fq\\n' > $outdir/fastq_QC/$sn.trimfq.list");
		#$queue->addcommond("java -jar $script/FastqStat.jar -i $outdir/fastq_QC/$sn.trimfq.list > $outdir/fastq_QC/$sn.trimdata.stat");
		#$queue->addcommond("cut -f 1-3,11-14 $outdir/fastq_QC/$sn.trimdata.stat > $outdir/fastq_QC/$sn.trimdata.stat.xls");
		#$queue->addcommond("rm $outdir/fastq_QC/$sn.trimdata.stat");
		
		if($refer){
			$queue->addcommond("bwa mem -t 20 $refer $outdir/fastq_QC/$sn.trim.R1.fq $outdir/fastq_QC/$sn.trim.R2.fq > $outdir/fastq_QC/$sn.sam");
			$queue->addcommond("samtools view -f 12 -S -h -b -\@ 20 $outdir/fastq_QC/$sn.sam | samtools sort -n -\@ 20 -o $outdir/fastq_QC/$sn.sorted.bam");
			$queue->addcommond("bedtools bamtofastq -i $outdir/fastq_QC/$sn.sorted.bam -fq $outdir/fastq_QC/$sn.clean.R1.fq  -fq2 $outdir/fastq_QC/$sn.clean.R2.fq");
			$queue->addcommond("rm $outdir/fastq_QC/$sn.sam;rm $outdir/fastq_QC/$sn.sorted.bam");
			
			$queue->addcommond("echo -e '$sn\\t$outdir/fastq_QC/$sn.clean.R1.fq\\t$outdir/fastq_QC/$sn.clean.R2.fq\\n' > $outdir/fastq_QC/$sn.cleanfq.list");
			#$queue->addcommond("java -jar $script/FastqStat.jar -i $outdir/fastq_QC/$sn.cleanfq.list > $outdir/fastq_QC/$sn.cleandata.stat");
			#$queue->addcommond("cut -f 1-3,11-14 $outdir/fastq_QC/$sn.cleandata.stat > $outdir/fastq_QC/$sn.cleandata.stat.xls");
			$queue->addcommond("rm $outdir/fastq_QC/$sn.cleandata.stat");
			$queue->addcommond("rm $outdir/fastq_QC/$sn.trim.*.fq");
			$queue->addcommond("rm $outdir/fastq_QC/$sn.R1.fq;rm $outdir/fastq_QC/$sn.R2.fq");
		}else{
			$queue->addcommond("ln -sf $outdir/fastq_QC/$sn.trim.R1.fq $outdir/fastq_QC/$sn.clean.R1.fq");
			$queue->addcommond("ln -sf $outdir/fastq_QC/$sn.trim.R2.fq $outdir/fastq_QC/$sn.clean.R2.fq");
			$queue->addcommond("ln -sf $outdir/fastq_QC/$sn.trimdata.stat.xls $outdir/fastq_QC/$sn.cleandata.stat.xls");
		}
		
		#$queue->addcommond("$software/fastqc/0.11.4/fastqc -t 10 --nogroup --extract $outdir/fastq_QC/$sn.R1.fq $outdir/fastq_QC/$sn.R2.fq &> $outdir/fastq_QC/$sn.fastqc.log");
		
		
		############################
		
		$queue->addcommond("mkdir -p $outdir/assembly");
		$queue->addcommond("$software/megahit/1.1.3/megahit -1 $outdir/fastq_QC/$sn.clean.R1.fq -2 $outdir/fastq_QC/$sn.clean.R2.fq -o $outdir/assembly/$sn -t 20");
		$queue->addcommond("$software/biocode/fasta/filter_fasta_by_size.pl -f $outdir/assembly/$sn/final.contigs.fa -s 500 -u $outdir/assembly/$sn/filtered_contig.fa");
		$queue->addcommond("$software/seqkit/v0.9.1/seqkit stat -a -T $outdir/assembly/$sn/filtered_contig.fa > $outdir/assembly/$sn/contig_stat.txt");
                $queue->addcommond("perl $script/Statistics_for_contigs.pl $outdir/assembly/");

		$queue->addcommond("gzip -c $outdir/fastq_QC/$sn.clean.R1.fq > $outdir/fastq_QC/$sn.clean.R1.fq.gz;gzip -c $outdir/fastq_QC/$sn.clean.R2.fq > $outdir/fastq_QC/$sn.clean.R2.fq.gz;rm $outdir/fastq_QC/$sn.clean.*fq;");
		#############################
		
		$queue->addcommond("mkdir -p $outdir/gene_predict/$sn");
		$queue->addcommond("$software/MetaGeneMark/mgm/gmhmmp  -f G  -m $software/MetaGeneMark/mgm/MetaGeneMark_v1.mod  $outdir/assembly/$sn/final.contigs.fa  -o $outdir/gene_predict/$sn/seq.gff  -A $outdir/gene_predict/$sn/seq.faa.tmp  -D $outdir/gene_predict/$sn/seq.fna.tmp");
		$queue->addcommond("$script/parse_MetaGeneMark.pl -faa $outdir/gene_predict/$sn/seq.faa.tmp -fna $outdir/gene_predict/$sn/seq.fna.tmp -name $sn -outdir $outdir/gene_predict/$sn");
		$queue->addcommond("$software/biocode/fasta/filter_fasta_by_size.pl -f $outdir/gene_predict/$sn/seq.fna  -s 100 -u $outdir/gene_predict/$sn/seq_more100.fna");
                $queue->addcommond("perl $script/Statistics_for_orf.pl $outdir/gene_predict/");
	}
}

sub rawdata_qc{
	return if $single_sample;
	foreach my $sn(sort keys %samples){
		$queue->set_job_cpu(20);
		$queue->set_job_mem(20);
		$queue->set_job_name("样本 $sn 数据质量评估与剪切");
		$queue->set_work_dir($workdir);
		
		$queue->addcommond("mkdir -p $outdir/fastq_QC/");
		foreach my $x(qw/R1 R2/){
			my @fqs = split / ,/,$samples{$sn}{$x};
			$queue->addcommond("rm -f $outdir/fastq_QC/$sn.$x.fq");
			foreach my $fq(@fqs){
				if($fq =~ /.gz$/){
					$queue->addcommond("zcat $fq >> $outdir/fastq_QC/$sn.$x.fq");
				}else{
					$queue->addcommond("cat  $fq >> $outdir/fastq_QC/$sn.$x.fq");
				}
			}
		}
		
		$queue->addcommond("echo -e '$sn\\t$outdir/fastq_QC/$sn.R1.fq\\t$outdir/fastq_QC/$sn.R2.fq\\n' > $outdir/fastq_QC/$sn.rawfq.list");
		$queue->addcommond("java -jar $script/FastqStat.jar -i $outdir/fastq_QC/$sn.rawfq.list > $outdir/fastq_QC/$sn.rawdata.stat");
		$queue->addcommond("cut -f 1-3,11-14 $outdir/fastq_QC/$sn.rawdata.stat > $outdir/fastq_QC/$sn.rawdata.stat.xls");
		$queue->addcommond("rm $outdir/fastq_QC/$sn.rawdata.stat");
		
		$queue->addcommond("cutadapt -q 20 -e 0.1 -n 1 -m 20 -a AGATCGGAAGAGC -g GCTCTTCCGATCT -o $outdir/fastq_QC/$sn.trim.R1.fq -p $outdir/fastq_QC/$sn.trim.R2.fq $outdir/fastq_QC/$sn.R1.fq $outdir/fastq_QC/$sn.R2.fq &> $outdir/fastq_QC/$sn.cutadapt.log");
		$queue->addcommond("rm $outdir/fastq_QC/$sn.R1.fq; rm $outdir/fastq_QC/$sn.R2.fq");
		$queue->addcommond("echo -e '$sn\\t$outdir/fastq_QC/$sn.trim.R1.fq\\t$outdir/fastq_QC/$sn.trim.R2.fq\\n' > $outdir/fastq_QC/$sn.trimfq.list");
		$queue->addcommond("java -jar $script/FastqStat.jar -i $outdir/fastq_QC/$sn.trimfq.list > $outdir/fastq_QC/$sn.trimdata.stat");
		$queue->addcommond("cut -f 1-3,11-14 $outdir/fastq_QC/$sn.trimdata.stat > $outdir/fastq_QC/$sn.trimdata.stat.xls");
		$queue->addcommond("rm $outdir/fastq_QC/$sn.trimdata.stat;");
		
		if($refer){
			$queue->addcommond ("bwa mem -t 8 $refer $outdir/fastq_QC/$sn.trim.R1.fq $outdir/fastq_QC/$sn.trim.R2.fq > $outdir/fastq_QC/$sn.sam");
			$queue->addcommond("rm $outdir/fastq_QC/$sn.trim.R1.fq $outdir/fastq_QC/$sn.trim.R2.fq");
			$queue->addcommond("samtools view -f 12 -S -h -b -\@ 20 $outdir/fastq_QC/$sn.sam | samtools sort -n -\@ 20 -o $outdir/fastq_QC/$sn.sorted.bam");
			$queue->addcommond("rm $outdir/fastq_QC/$sn.sam");
			$queue->addcommond("bedtools bamtofastq -i $outdir/fastq_QC/$sn.sorted.bam -fq $outdir/fastq_QC/$sn.clean.R1.fq  -fq2 $outdir/fastq_QC/$sn.clean.R2.fq");
			$queue->addcommond("rm $outdir/fastq_QC/$sn.sorted.bam");
			$queue->addcommond("echo -e '$sn\\t$outdir/fastq_QC/$sn.clean.R1.fq\\t$outdir/fastq_QC/$sn.clean.R2.fq\\n' > $outdir/fastq_QC/$sn.cleanfq.list");
			$queue->addcommond("java -jar $script/FastqStat.jar -i $outdir/fastq_QC/$sn.cleanfq.list > $outdir/fastq_QC/$sn.cleandata.stat");
			$queue->addcommond("cut -f 1-3,11-14 $outdir/fastq_QC/$sn.cleandata.stat > $outdir/fastq_QC/$sn.cleandata.stat.xls");
			$queue->addcommond("rm $outdir/fastq_QC/$sn.cleandata.stat");
		}else{
			$queue->addcommond("ln -sf $outdir/fastq_QC/$sn.trim.R1.fq $outdir/fastq_QC/$sn.clean.R1.fq");
			$queue->addcommond("ln -sf $outdir/fastq_QC/$sn.trim.R2.fq $outdir/fastq_QC/$sn.clean.R2.fq");
			
			$queue->addcommond("ln -sf $outdir/fastq_QC/$sn.trimdata.stat.xls $outdir/fastq_QC/$sn.cleandata.stat.xls");
		}
		
		$queue->addcommond("$software/fastqc/0.11.4/fastqc -t 10 --nogroup --extract $outdir/fastq_QC/$sn.R1.fq $outdir/fastq_QC/$sn.R2.fq &> $outdir/fastq_QC/$sn.fastqc.log");
		
		$queue->run();
	}
	$queue->wait();
	
	$queue->set_job_cpu(2);
	$queue->set_job_mem(2);
	$queue->set_job_name("所有数据质量评估整合");
	$queue->set_work_dir($workdir);

	my $cmd;
	$cmd = "$script/cat_files.pl -head 1 -outfile $outdir/fastq_QC/rawdata.stat.xls";
	foreach my $sn(sort keys %samples){
		$cmd .= " $outdir/fastq_QC/$sn.rawdata.stat.xls";
	}
	$queue->addcommond("$cmd");
	
	$cmd = "$script/cat_files.pl -head 1 -outfile $outdir/fastq_QC/trimdata.stat.xls";
	foreach my $sn(sort keys %samples){
		$cmd .= " $outdir/fastq_QC/$sn.trimdata.stat.xls";
	}
	$queue->addcommond("$cmd");
	
	$cmd = "$script/cat_files.pl -head 1 -outfile $outdir/fastq_QC/cleandata.stat.xls";
	foreach my $sn(sort keys %samples){
		$cmd .= " $outdir/fastq_QC/$sn.cleandata.stat.xls";
	}
	$queue->addcommond("$cmd");
	
	$queue->addcommond("# rm $outdir/fastq_QC/*.rawdata.stat.xls $outdir/fastq_QC/*.trimdata.stat.xls $outdir/fastq_QC/*.cleandata.stat.xls");
	
	$queue->run();
	$queue->wait();
}

sub assembly{
	return if $single_sample;
	foreach my $sn(sort keys %samples){
		$queue->set_job_cpu(20);
		$queue->set_job_mem(20);
		$queue->set_job_name("样本 $sn 宏基因组拼接");
		$queue->set_work_dir($workdir);
		
		$queue->addcommond("mkdir -p $outdir/assembly");
		$queue->addcommond("rm -rf $outdir/assembly/$sn");
		$queue->addcommond("$software/megahit/1.1.3/megahit -1 $outdir/fastq_QC/$sn.clean.R1.fq -2 $outdir/fastq_QC/$sn.clean.R2.fq -o $outdir/assembly/$sn -t 20");
		$queue->addcommond("rm -rf $outdir/assembly/$sn/intermediate_contigs");
		$queue->addcommond("$software/biocode/fasta/filter_fasta_by_size.pl -f $outdir/assembly/$sn/final.contigs.fa -s 500 -u $outdir/assembly/$sn/filtered_contig.fa");
		$queue->addcommond("$software/seqkit/v0.9.1/seqkit stat -a -T $outdir/assembly/$sn/filtered_contig.fa > $outdir/assembly/$sn/contig_stat.txt");
		$queue->run();
	}
	$queue->wait();
        `perl $script/Statistics_for_contigs.pl	$outdir/assembly`;
}

sub gene_predict{
	return if $single_sample;
	foreach my $sn(sort keys %samples){
		$queue->set_job_cpu(10);
		$queue->set_job_mem(10);
		$queue->set_job_name("样本 $sn 基因ORF预测");
		$queue->set_work_dir($workdir);
		$queue->addcommond("mkdir -p $outdir/gene_predict/$sn");
		$queue->addcommond("$software/MetaGeneMark/mgm/gmhmmp  -f G  -m $software/MetaGeneMark/mgm/MetaGeneMark_v1.mod  $outdir/assembly/$sn/filtered_contig.fa  -o $outdir/gene_predict/$sn/seq.gff  -A $outdir/gene_predict/$sn/seq.faa.tmp  -D $outdir/gene_predict/$sn/seq.fna.tmp");
		$queue->addcommond("$script/parse_MetaGeneMark.pl -faa $outdir/gene_predict/$sn/seq.faa.tmp -fna $outdir/gene_predict/$sn/seq.fna.tmp -name $sn -outdir $outdir/gene_predict/$sn");
		$queue->addcommond("$software/biocode/fasta/filter_fasta_by_size.pl -f $outdir/gene_predict/$sn/seq.fna  -s 100 -u $outdir/gene_predict/$sn/seq_more100.fna");
		#$queue->addcommond("$script/Statistics.pl $outdir/gene_predict/$sn/seq_more100.fna > $outdir/gene_predict/$sn/seq_more100.fna.stats.xls");
		$queue->run();
	}
	$queue->wait();
        `perl $script/Statistics_for_orf.pl $outdir/gene_predict")`;
}

sub unigene_set{
	my $cmd = "cat";
	foreach my $sn(sort keys %samples){
		$cmd .= " $outdir/gene_predict/$sn/seq_more100.fna";
	}
	$cmd .= " >> $outdir/unigene_set/merged_samples_orfs.fna";
	$queue->set_job_cpu(20);
	$queue->set_job_mem(20);
	$queue->set_job_name("合并所有样本基因并去冗余");
	$queue->set_work_dir($workdir);
	
	$queue->addcommond("mkdir -p $outdir/unigene_set/ ");
	$queue->addcommond("rm -f $outdir/unigene_set/merged_samples_orfs.fna");
	$queue->addcommond("$cmd");
	$queue->addcommond("$software/cd-hit/V4.6.8/cd-hit-est  -i $outdir/unigene_set/merged_samples_orfs.fna  -o $outdir/unigene_set/uniGeneSet.fna  -c 0.95  -aS 0.9  -T 20  -M 0");
	$queue->addcommond("$software/EMBOSS/6.6.0/transeq -sequence $outdir/unigene_set/uniGeneSet.fna -table 11 -trim -outseq $outdir/unigene_set/uniGeneSet.faa.tmp");
	$queue->addcommond("sed 's/_1\$//g' $outdir/unigene_set/uniGeneSet.faa.tmp > $outdir/unigene_set/uniGeneSet.faa");
	$queue->addcommond("$software/soap/2.21/2bwt-builder $outdir/unigene_set/uniGeneSet.fna");
        #$queue->addcommond("$script/Statistics.pl $outdir/unigene_set/uniGeneSet.fna > $outdir/unigene_set/uniGeneSet.fna.stats.xls");
        $queue->addcommond("$script/Statistics_for_non_orf.pl $outdir/unigene_set");
	$queue->run();

	$queue->wait();
}

sub gene_profile{
	my $min = $insertsize - 100;
	my $max = $insertsize + 100;
	foreach my $sn(sort keys %samples){
		$queue->set_job_cpu(20);
		$queue->set_job_mem(20);
		$queue->set_job_name("样本 $sn reads比对非冗余基因序列");
		$queue->set_work_dir($workdir);
		
		$queue->addcommond("mkdir -p $outdir/gene_profile/$sn");
		$queue->addcommond("$software/soap/2.21/soap  -a $outdir/fastq_QC/$sn.clean.R1.fq  -b $outdir/fastq_QC/$sn.clean.R2.fq  -D $outdir/unigene_set/uniGeneSet.fna.index -o $outdir/gene_profile/$sn/soap.pair.pe  -2 $outdir/gene_profile/$sn/soap.pair.se  -r 1  -l 35  -M 4  -p 6  -v 20  -c 0.95  -m $min  -x $max ");
		$queue->run();
	}
	$queue->wait();
	
	$queue->set_job_cpu(10);
	$queue->set_job_mem(10);
	$queue->set_job_name("Gene Profile矩阵构建");
	$queue->set_work_dir($workdir);
	
	$queue->addcommond("rm -f $outdir/gene_profile/soap.info");
	foreach my $sn(sort keys %samples){
		$queue->addcommond("echo -e '$sn\\t$insertsize\\t$outdir/gene_profile/$sn/soap.pair.pe,$outdir/gene_profile/$sn/soap.pair.se'>> $outdir/gene_profile/soap.info");
	}
	
	$queue->addcommond("$script/gene_profile.pl -refer $outdir/unigene_set/uniGeneSet.fna -infor $outdir/gene_profile/soap.info -outdir $outdir/gene_profile");
	
	$queue->run();
	$queue->wait();
}

sub gene_annotation{
	foreach my $db_name(qw/nr kegg eggnog cazy ardb card phi mvirdb tcdb qs vfdb signalp pfam t3ss/){
		$queue->set_job_cpu(30);
		$queue->set_job_mem(50);
		$queue->set_job_local(1);
		$queue->set_job_name("基因注释：$db_name");
		$queue->set_work_dir($workdir);
		
		$queue->addcommond("mkdir -p $outdir/gene_annotation/$db_name");
		
		#$queue->addcommond("cd $outdir/gene_annotation/$db_name; /home/pub/app/annotation_meta/annot-$db_name-meta-v1.0.pl -fasta $outdir/unigene_set/uniGeneSet.faa -ftype pro -profile $outdir/gene_profile/gene_profile.reads_number.txt; cd -");
		#if()
		$queue->addcommond("cd $outdir/gene_annotation/$db_name; /home/pub/app/annotation_meta/annot-$db_name-meta-v1.0.pl -fasta $outdir/unigene_set/uniGeneSet.faa -ftype pro -profile $outdir/gene_profile/gene_profile.TPM.txt; cd -");
		
		$queue->run();
		$queue->wait();
	}
}

sub profile_collection{
	$queue->set_job_cpu(10);
	$queue->set_job_mem(10);
	$queue->set_job_local(1);
	$queue->set_job_name("各种表达矩阵整理和格式转化");
	$queue->set_work_dir($workdir);
	
	$queue->addcommond("mkdir -p $outdir/profile_collection");
	
	$queue->addcommond("cp $outdir/gene_profile/gene_profile.TPM.txt   $outdir/profile_collection/gene_table.tab");
	$queue->addcommond("cp $outdir/gene_annotation/nr/outdir/gene_taxonomy_profile.xls $outdir/profile_collection/gene_taxonomy_table.tab");
	$queue->addcommond("cp $outdir/gene_annotation/nr/outdir/nr.taxonomy_profile.xls $outdir/profile_collection/taxonomy_table.tab");
	foreach my $mm(qw/kingdom phylum class order family genus species/){
		$queue->addcommond("cp $outdir/gene_annotation/nr/outdir/summary/$mm.xls $outdir/profile_collection/nr_$mm\_table.tab");
	}
	$queue->addcommond("cp $outdir/gene_annotation/kegg/outdir/pathway_table.xls $outdir/profile_collection/kegg_table.tab");
	$queue->addcommond("cp $outdir/gene_annotation/eggnog/outdir/eggnog.cog_table.xls $outdir/profile_collection/eggnog_table.tab");
	$queue->addcommond("cp $outdir/gene_annotation/cazy/outdir/cazy.family.profile.xls $outdir/profile_collection/cazy_family_table.tab");
	$queue->addcommond("cp $outdir/gene_annotation/cazy/outdir/cazy.class.profile.xls $outdir/profile_collection/cazy_class_table.tab");
	$queue->addcommond("cp $outdir/gene_annotation/ardb/outdir/ardb.class.profile.xls $outdir/profile_collection/ardb_class_table.tab");
	$queue->addcommond("cp $outdir/gene_annotation/ardb/outdir/ardb.type.profile.xls $outdir/profile_collection/ardb_type_table.tab");
	$queue->addcommond("cp $outdir/gene_annotation/card/outdir/gene.card.profile.xls $outdir/profile_collection/card_table.tab");
	
	foreach my $nn(qw/taxonomy nr_kingdom nr_phylum nr_class nr_order nr_family nr_genus nr_species eggnog kegg cazy_family cazy_class ardb_class ardb_type card gene/){
		$queue->addcommond("$script/matrix_sort.pl -in $outdir/profile_collection/$nn\_table.tab -out $outdir/profile_collection/$nn\_table.txt -type down -remove");
		if($nn eq "gene"){
			$queue->addcommond("biom convert -i $outdir/profile_collection/$nn\_table.txt -o $outdir/profile_collection/$nn\_table.biom --table-type 'Gene table' --to-hdf5");
		}elsif($nn =~ /^nr_/){
			$queue->addcommond("biom convert -i $outdir/profile_collection/$nn\_table.txt -o $outdir/profile_collection/$nn\_table.biom --table-type 'Taxon table' --to-hdf5");
		}elsif($nn =~ /^taxonomy/){
			$queue->addcommond("biom convert -i $outdir/profile_collection/$nn\_table.txt -o $outdir/profile_collection/$nn\_table.biom --table-type 'Taxon table' --to-hdf5");
		}else{
			$queue->addcommond("biom convert -i $outdir/profile_collection/$nn\_table.txt -o $outdir/profile_collection/$nn\_table.biom --table-type 'Function table' --to-hdf5");
		}
	}
	$queue->addcommond("$script/matrix_sort.pl -in $outdir/profile_collection/gene_taxonomy_table.tab -out $outdir/profile_collection/gene_taxonomy_table.txt -type down -cols taxonomy -remove");
	$queue->addcommond("biom convert -i $outdir/profile_collection/gene_taxonomy_table.txt -o $outdir/profile_collection/gene_taxonomy_table.biom --table-type 'Gene table' --to-hdf5 --process-obs-metadata taxonomy");
	
	$queue->addcommond("###############################");
	
	foreach my $nn(qw/taxonomy nr_kingdom nr_phylum nr_class nr_order nr_family nr_genus nr_species eggnog kegg cazy_family cazy_class ardb_class ardb_type card gene/){
		foreach my $std(sort keys %group_info){
			next if (scalar keys %{$group_info{$std}} < 2);
			$queue->addcommond("mkdir -p $outdir/profile_collection/group/$nn/$std");
			$queue->addcommond("$script/matrix_sort_and_select.pl --input $outdir/profile_collection/$nn\_table.txt --output $outdir/profile_collection/group/$nn/$std/data.txt --group $outdir/group_info/$std.txt -type down -remove");
			
			if($nn eq "gene"){
				$queue->addcommond("biom convert -i $outdir/profile_collection/group/$nn/$std/data.txt -o $outdir/profile_collection/group/$nn/$std/data.biom --table-type 'Gene table' --to-hdf5");
			}elsif($nn =~ /^nr_/){
				$queue->addcommond("biom convert -i $outdir/profile_collection/group/$nn/$std/data.txt -o $outdir/profile_collection/group/$nn/$std/data.biom --table-type 'Taxon table' --to-hdf5");
			}elsif($nn =~ /^taxonomy/){
				$queue->addcommond("biom convert -i $outdir/profile_collection/group/$nn/$std/data.txt -o $outdir/profile_collection/group/$nn/$std/data.biom --table-type 'Taxon table' --to-hdf5");
			}else{
				$queue->addcommond("biom convert -i $outdir/profile_collection/group/$nn/$std/data.txt -o $outdir/profile_collection/group/$nn/$std/data.biom --table-type 'Function table' --to-hdf5");
			}
		}
	}
	
	$queue->run();
	$queue->wait();
}

sub composition{
	$queue->set_job_cpu(2);
	$queue->set_job_mem(2);
	$queue->set_job_name("物种与功能组成分析:venn");
	$queue->set_work_dir($workdir);
	$queue->addcommond("mkdir -p $outdir/composition/venn/");
	#foreach my $nn(qw/nr_kingdom nr_phylum nr_class nr_order nr_family nr_genus nr_species eggnog kegg cazy_family cazy_class ardb_class ardb_type card gene/){
	foreach my $nn(qw/nr_kingdom nr_phylum nr_class nr_order nr_family nr_genus nr_species eggnog kegg cazy_family cazy_class ardb_class ardb_type card/){
		foreach my $std(sort keys %group_info){
			$queue->addcommond("mkdir -p $outdir/composition/venn/$std");
			$queue->addcommond("$script/plot-venn.pl -i $outdir/profile_collection/$nn\_table.txt -g $outdir/group_info/$std.txt -o $outdir/composition/venn/$std/venn_$nn");
		}
	}
	$queue->run();

	$queue->set_job_cpu(2);
	$queue->set_job_mem(2);
	$queue->set_job_name("物种与功能组成分析:heatmap");
	$queue->set_work_dir($workdir);
	$queue->addcommond("mkdir -p $outdir/composition/heatmap/");
	foreach my $nn(qw/nr_kingdom nr_phylum nr_class nr_order nr_family nr_genus nr_species eggnog kegg cazy_family cazy_class ardb_class ardb_type card gene/){
		foreach my $std(sort keys %group_info){
			$queue->addcommond("mkdir -p $outdir/composition/heatmap/$nn/$std");
			$queue->addcommond("$script/plot-pheatmap.pl -i $outdir/profile_collection/group/$nn/$std/data.txt -o $outdir/composition/heatmap/$nn/$std/heatmap -top 50 -sorted 1 -map $outdir/group_info/$std.info -g $std");
		}
	}
	$queue->run();
	
	$queue->set_job_cpu(2);
	$queue->set_job_mem(2);
	$queue->set_job_name("物种与功能组成分析:bubble");
	$queue->set_work_dir($workdir);
	$queue->addcommond("mkdir -p $outdir/composition/bubble/");
	foreach my $nn(qw/nr_kingdom nr_phylum nr_class nr_order nr_family nr_genus nr_species eggnog kegg cazy_family cazy_class ardb_class ardb_type card gene/){
		foreach my $std(sort keys %group_info){
			$queue->addcommond("mkdir -p $outdir/composition/bubble/$nn/$std");
			$queue->addcommond("$script/plot-bubble.pl -i $outdir/profile_collection/group/$nn/$std/data.txt -o $outdir/composition/bubble/$nn/$std/bubble -top 50 -sorted 1 -g $outdir/group_info/$std.txt -ylab $nn");
		}
	}
	$queue->run();

	$queue->set_job_cpu(2);
	$queue->set_job_mem(10);
	$queue->set_job_name("物种与功能组成分析:core_pan");
	$queue->set_work_dir($workdir);
	$queue->addcommond("mkdir -p $outdir/composition/core_pan");
	$queue->addcommond("cd $outdir/composition/core_pan");
	$queue->addcommond("compute_core_microbiome.py -i $outdir/profile_collection/gene_taxonomy_table.biom -o . --min_fraction_for_core 0.5 --max_fraction_for_core 1");
	$queue->addcommond("cd -");
	$queue->run();
	
	$queue->set_job_cpu(2);
	$queue->set_job_mem(2);
	$queue->set_job_name("物种与功能组成分析:megan");
	$queue->set_work_dir($workdir);
	$queue->addcommond("mkdir -p $outdir/composition/megan");
	
	$queue->addcommond("mkdir -p $outdir/composition/megan/eachSample");
	$queue->addcommond("$script/level_tree_sample_input.2.pl $outdir/profile_collection/taxonomy_table.txt $outdir/composition/megan/eachSample/tax");
	$queue->addcommond("cd $outdir/composition/megan/eachSample");
	$queue->addcommond("$script/level_tree_sample.float.2.py -f tax.level_tree_sample.profile.xls -t 20 -r 10 --config tax.level_tree_sample.config");
	$queue->addcommond("cd -");
	
	$queue->addcommond("mkdir -p $outdir/composition/megan/multiSamples");
	$queue->addcommond("$script/level_tree_input.2.pl $outdir/profile_collection/taxonomy_table.txt $outdir/composition/megan/multiSamples/tax");
	foreach my $std(sort keys %group_info){
		$queue->addcommond("mkdir -p $outdir/composition/megan/multiSamples/$std");
		$queue->addcommond("cd $outdir/composition/megan/multiSamples/$std");
		$queue->addcommond("$script/level_tree.float.py -f ../tax.level_tree.profile.xls -t 10 -g $outdir/group_info/$std.txt --hi 2000");
	}
	$queue->addcommond("cd -");
	
	$queue->run();
	
	## 设计分组
	$queue->set_job_cpu(2);
	$queue->set_job_mem(2);
	$queue->set_job_name("物种与功能组成分析:rarefaction");
	$queue->set_work_dir($workdir);
	$queue->addcommond("mkdir -p $outdir/composition/rarefaction");
	foreach my $std(sort keys %group_info){
		$queue->addcommond("cd $outdir/composition/rarefaction");
		
		$queue->addcommond("# multiple_rarefactions.py -i $outdir/profile_collection/gene_taxonomy_table.biom -m 10 -x 16544 -s1653 -o data_rarefaction");
		$queue->addcommond("# alpha_diversity.py -i data_rarefaction -o data_alpha_div --metrics ACE,shannon,chao1,observed_species,goods_coverage,simpson");
		$queue->addcommond("# collate_alpha.py -i data_alpha_div -o data_alpha_div_collated");
		$queue->addcommond("# make_rarefaction_plots.py -i data_alpha_div_collated -m map.txt -o data_alpha_rarefaction_plots");
		# $queue->addcommond("mkdir -p $outdir/composition/core_pan/$nn");
		# $queue->addcommond("cd $outdir/community_composition/core_pan/$nn");
		# #alpha_rarefaction.py -i otu_table/otu_table.biom -m map.txt -o div_alpha/ -p alpha_params.txt -t rep_phylo.tre
		
		# $queue->addcommond("cp $outdir/otu/original/otus.shared $outdir/alpha_diversity/");
		# $queue->addcommond("cd $outdir/alpha_diversity/");
		# $queue->addcommond("$software/mothur/1.41.0/mothur \"#summary.single(shared=otus.shared,calc=ace-chao-shannon-simpson-coverage,groupmode=f); rarefaction.single(shared=otus.shared,calc=sobs-chao-shannon,groupmode=f,freq=100)\" &> log.mothur.summary_single");
		# $queue->addcommond("$script/plot-rarefaction.pl  -d . -a t -l 0.97 -o ");
		# $queue->addcommond("$script/plot-rarefaction.pl  -d . -a t -l 0.97 -m shannon ");
		# $queue->addcommond("$script/shannon-ace-table.pl -d . -o estimators.html");
		
		# $queue->addcommond("compute_core_microbiome.py -i $outdir/profile_collection/$nn\_table.biom -o . --min_fraction_for_core 0.5 --max_fraction_for_core 1");
		# $queue->addcommond("cd -");
		$queue->addcommond("cd -");
	}
	$queue->run();
	
	$queue->set_job_cpu(2);
	$queue->set_job_mem(2);
	$queue->set_job_name("物种与功能组成分析:graphlan物种组成树");
	$queue->set_work_dir($workdir);
	$queue->addcommond("mkdir -p $outdir/composition/graphlan");
	foreach my $std(sort keys %group_info){
		$queue->addcommond("cd $outdir/composition/graphlan");
		$queue->addcommond("biom convert -i $outdir/profile_collection/gene_taxonomy_table.biom -o $outdir/composition/graphlan/gene_taxonomy_table.txt --to-tsv --header-key taxonomy --output-metadata-id 'Consensus Lineage'"); 
		$queue->addcommond("qiime2lefse.py --in $outdir/composition/graphlan/gene_taxonomy_table.txt --out $outdir/composition/graphlan/gene_taxonomy_table.txt");
		
		$queue->addcommond("cd -");
	}
	$queue->run();
	
	$queue->set_job_cpu(2);
	$queue->set_job_mem(2);
	$queue->set_job_name("物种与功能组成分析:优势物种分析");
	$queue->set_work_dir($workdir);
	$queue->addcommond("mkdir -p $outdir/composition/top");
	$queue->addcommond("cd $outdir/composition/top");
	$queue->addcommond("$script/drawTopGenus.py  -i $outdir/profile_collection/nr_genus_table.txt -taxa $outdir/otu/normalize/otu_taxa_table.xls -o bar -t 30");
	$queue->addcommond("cd -");
	$queue->run();
	
	$queue->set_job_cpu(2);
	$queue->set_job_mem(2);
	$queue->set_job_name("物种与功能组成分析:单样本多级物种组成");
	$queue->set_work_dir($workdir);
	$queue->addcommond("mkdir -p $outdir/composition/pie");
	foreach my $sam(sort keys %samples){
		$queue->addcommond("mkdir -p $outdir/composition/pie/$sam");
		$queue->addcommond("cd $outdir/composition/pie");
		$queue->addcommond(".............");
		$queue->addcommond("cd -");
	}
	
	$queue->run();
		
	#### accumulation ####
	#$queue->addcommond("mkdir -p $outdir/composition/accumulation/gene");
	#$queue->addcommond("cd $outdir/composition/accumulation/gene");
	#$queue->addcommond("$script/shared_otu_curves.pl  -in $outdir/gene_profile/gene_profile.base_percent.txt -out gene.accumulation -ylp Genes \n");
	#$queue->addcommond("cd -");
	#$queue->run();
	#$queue->wait();
}

sub comparison{
	$queue->set_job_cpu(10);
	$queue->set_job_mem(10);

	$queue->set_job_name("物种与功能比较分析:PCA主成分分析");
	$queue->set_work_dir($workdir);
	$queue->addcommond("mkdir -p $outdir/comparison/pca");
	foreach my $nn(qw/nr_kingdom nr_phylum nr_class nr_order nr_family nr_genus nr_species eggnog kegg cazy_family cazy_class ardb_class ardb_type card gene/){
		foreach my $std(sort keys %group_info){
			$queue->addcommond("mkdir -p $outdir/comparison/pca/$nn/$std");
			$queue->addcommond("$script/plot-pca-3d.pl -i $outdir/profile_collection/group/$nn/$std/data.txt -o $outdir/comparison/pca/$nn/$std -m $outdir/group_info/$std.info -g $std  -pc 1-2-3 -w 8 -h 8 -scale T -cen T -f T");
			$queue->addcommond("rename pca_  $nn\_$std\_pca_ $outdir/comparison/pca/$nn/$std/pca_*");
		}
	}
	$queue->run();

	$queue->set_job_cpu(2);
	$queue->set_job_mem(2);
	$queue->set_job_name("物种与功能比较分析:样本间距离计算");
	$queue->set_work_dir($workdir);
	$queue->addcommond("mkdir -p $outdir/comparison/distance");
	foreach my $nn(qw/nr_kingdom nr_phylum nr_class nr_order nr_family nr_genus nr_species eggnog kegg cazy_family cazy_class ardb_class ardb_type card gene/){
	    $queue->addcommond("mkdir -p $outdir/comparison/distance/$nn");
	    $queue->addcommond("beta_diversity.py -i $outdir/profile_collection/$nn\_table.biom -o $outdir/comparison/distance/$nn -m bray_curtis,euclidean,hellinger");
	    $queue->addcommond("beta_diversity.py -i $outdir/profile_collection/$nn\_table.biom -o $outdir/comparison/distance/$nn -m abund_jaccard,bray_curtis,euclidean,hellinger");
	    foreach my $dist(qw/bray_curtis euclidean hellinger/){
		$queue->addcommond("$script/plot-heatmap-dis.pl -i $outdir/comparison/distance/$nn/$dist\_$nn\_table.txt -o $outdir/comparison/distance/$nn/$dist\_$nn\_heatmap -lw 0.1:1:0.1:7.5 -lh 1:0.1:7.5:1.8 -h 6 -slas 2 -col gray-lightyellow-darkorange  -rlc 0.4  -clc 0.4");
	        $queue->addcommond("cd $outdir/comparison/distance/$nn/$dist/");
	        $queue->addcommond("rename heatmap  $nn\_$dist\_heatmap  heatmap*");
            }
	}
	$queue->run();
	$queue->wait();

	$queue->set_job_cpu(10);
	$queue->set_job_mem(10);
	$queue->set_job_name("物种与功能比较分析:PCoA主坐标分析");
	$queue->set_work_dir($workdir);
	$queue->addcommond("mkdir -p $outdir/comparison/pcoa");
	foreach my $nn(qw/nr_kingdom nr_phylum nr_class nr_order nr_family nr_genus nr_species eggnog kegg cazy_family cazy_class ardb_class ardb_type card gene/){
		foreach my $dist(qw/bray_curtis euclidean hellinger/){
			foreach my $std(sort keys %group_info){
				$queue->addcommond("mkdir -p $outdir/comparison/pcoa/$nn/$dist/$std");
				$queue->addcommond("$script/plot-pcoa.pl -i $outdir/comparison/distance/$nn/$dist\_$nn\_table.txt -o $outdir/comparison/pcoa/$nn/$dist/$std -m $outdir/group_info/$std.info -g $std -pc  1-2 -w 8 -h 8");
				$queue->addcommond("$script/plot-pcoa.pl -i $outdir/comparison/distance/$nn/$dist\_$nn\_table.txt -o $outdir/comparison/pcoa/$nn/$dist/$std -m $outdir/group_info/$std.info -g $std -pc  1-3 -w 8 -h 8");
				$queue->addcommond("$script/plot-pcoa.pl -i $outdir/comparison/distance/$nn/$dist\_$nn\_table.txt -o $outdir/comparison/pcoa/$nn/$dist/$std -m $outdir/group_info/$std.info -g $std -pc  2-3 -w 8 -h 8");
	                    $queue->addcommond("cd $outdir/comparison/pcoa/$nn/$dist/$std/");
			    $queue->addcommond("rename  pcoa_  $nn\_$dist\_$std\_pcoa_   $outdir/comparison/pcoa/$nn/$dist/$std/pcoa_*");
			}
		}
	}
	$queue->run();
	
	$queue->set_job_cpu(10);
	$queue->set_job_mem(10);
	$queue->set_job_name("物种与功能比较分析:NMDS非度量多维尺度分析");
	$queue->set_work_dir($workdir);
	$queue->addcommond("mkdir -p $outdir/comparison/nmds");
	foreach my $nn(qw/nr_kingdom nr_phylum nr_class nr_order nr_family nr_genus nr_species eggnog kegg cazy_family cazy_class ardb_class ardb_type card gene/){
		foreach my $dist(qw/bray_curtis euclidean hellinger/){
			foreach my $std(sort keys %group_info){
		            $queue->addcommond("mkdir -p $outdir/comparison/nmds/$nn/$dist/$std");
			    $queue->addcommond("$script/plot-nmds.pl -i $outdir/comparison/distance/$nn/$dist\_$nn\_table.txt -o $outdir/comparison/nmds/$nn/$dist/$std -m $outdir/group_info/$std.info -g $std");
	                    $queue->addcommond("cd $outdir/comparison/nmds/$nn/$dist/$std/");
			    $queue->addcommond("rename  nmds  $nn\_$dist\_$std\_nmds  nmds*");
			}
		}
	}
	$queue->run();

	$queue->set_job_cpu(10);
	$queue->set_job_mem(10);
	$queue->set_job_name("物种与功能比较分析:ANOSIM相似性分析");
	$queue->set_work_dir($workdir);
	$queue->addcommond("mkdir -p $outdir/comparison/anosim");
	foreach my $nn(qw/nr_kingdom nr_phylum nr_class nr_order nr_family nr_genus nr_species eggnog kegg cazy_family cazy_class ardb_class ardb_type card gene/){
		foreach my $dist(qw/bray_curtis euclidean hellinger/){
			foreach my $std(sort keys %group_info){
				$queue->addcommond("mkdir -p $outdir/comparison/anosim/$nn/$dist/$std");
				$queue->addcommond("cd $outdir/comparison/anosim/$nn/$dist/$std");
				$queue->addcommond("$script/dist_select.pl -i $outdir/comparison/distance/$nn/$dist\_$nn\_table.txt -o $outdir/comparison/anosim/$nn/$dist/$std/distance.txt -m $outdir/group_info/$std.info -g $std");
				$queue->addcommond("compare_categories.py --method anosim -i $outdir/comparison/anosim/$nn/$dist/$std/distance.txt -m $outdir/group_info/$std.info -c $std -o $outdir/comparison/anosim/$nn/$dist/$std");
			        $queue->addcommond("$script/anosim.py -i $outdir/comparison/anosim/$nn/$dist/$std/distance.txt -g $outdir/group_info/$std.txt -o $outdir/comparison/anosim/$nn/$dist/$std/anosim.pdf");
				$queue->addcommond("cd -");  
			}
		}
	}
	$queue->run();
	
	$queue->set_job_cpu(2);
	$queue->set_job_mem(2);
	$queue->set_job_name("物种与功能比较分析:Adonis多因素方差分析");
	$queue->set_work_dir($workdir);
	$queue->addcommond("mkdir -p $outdir/comparison/adonis");
	foreach my $nn(qw/nr_kingdom nr_phylum nr_class nr_order nr_family nr_genus nr_species eggnog kegg cazy_family cazy_class ardb_class ardb_type card gene/){
		foreach my $dist(qw/bray_curtis euclidean hellinger/){
			foreach my $std(sort keys %group_info){
				$queue->addcommond("mkdir -p $outdir/comparison/adonis/$nn/$dist/$std");
				$queue->addcommond("$script/dist_select.pl -i $outdir/comparison/distance/$nn/$dist\_$nn\_table.txt -o $outdir/comparison/adonis/$nn/$dist/$std/distance.txt -m $outdir/group_info/$std.info -g $std");
				$queue->addcommond("cd $outdir/comparison/adonis/$nn/$dist/$std");
				$queue->addcommond("compare_categories.py --method adonis -i $outdir/comparison/adonis/$nn/$dist/$std/distance.txt -m $outdir/group_info/$std.info -c $std -o $outdir/comparison/adonis/$nn/$dist/$std");
				$queue->addcommond("cd -");
			}
		}
	}
	$queue->run();
	
	$queue->set_job_cpu(10);
	$queue->set_job_mem(10);
	$queue->set_job_name("物种与功能比较分析:Hcluster层次聚类树");
	$queue->set_work_dir($workdir);
	$queue->addcommond("mkdir -p $outdir/comparison/hcluster_tree");
	foreach my $nn(qw/nr_kingdom nr_phylum nr_class nr_order nr_family nr_genus nr_species eggnog kegg cazy_family cazy_class ardb_class ardb_type card gene/){
		#foreach my $dist(qw/abund_jaccard bray_curtis euclidean hellinger/){
		foreach my $dist(qw/bray_curtis euclidean hellinger/){
			$queue->addcommond("mkdir -p $outdir/comparison/hcluster_tree/$nn/$dist");
			$queue->addcommond("cd $outdir/comparison/hcluster_tree/$nn/$dist");
			$queue->addcommond("$script/plot-hcluster_tree.pl -i $outdir/comparison/distance/$nn/$dist\_$nn\_table.txt -o $outdir/comparison/hcluster_tree/$nn/$dist");
			$queue->addcommond("rename hcluster_tree $nn\_$dist\_hcluster_tree hcluster_tree*");
			$queue->addcommond("cd -");
		}
	}
	$queue->run();

=head	
	$queue->set_job_cpu(10);
	$queue->set_job_mem(10);
	$queue->set_job_name("物种与功能比较分析:phylogeny_tree");
	$queue->set_work_dir($workdir);
	$queue->addcommond("mkdir -p $outdir/comparison/phylogeny_tree");
	$queue->addcommond("cd $outdir/comparison/phylogeny_tree");
	foreach my $level(qw/phylum class order family genus/){
	    $queue->addcommond("$script/spe.duo.pl $outdir/profile_collection/gene_taxonomy_table.txt $outdir/otu/normalize/otu_reps.fa $level");
            $queue->addcommond("align_seqs.py  -i $level.rep.fasta -m muscle  -o .");
	    $queue->addcommond("FastTree -nt $level.rep_aligned.fasta > $level.rep_aligned.tre");
            $queue->addcommond("$script/plot-tree.pl  -i genus.rep_aligned.fasttre -o 1.tree.pdf -h 20 -w 8 -cex 0.75 -eg 1.2");
	    $queue->addcommond("$script/plot-tree.pl  -i genus.rep_aligned.fasttre -o 2.tree.pdf -h 20 -w 20 -cex 0.8 -eg 1 -tretype fan  -ma T");
	}
	$queue->run();
	$queue->wait();
=cut

}

sub difference{
	$queue->set_job_cpu(10);
	$queue->set_job_mem(10);
	$queue->set_job_name("物种与功能差异分析:组间常规差异分析(Wilcoxon和kruskal_wallis)");
	$queue->set_work_dir($workdir);
	foreach my $nn(qw/nr_kingdom nr_phylum nr_class nr_order nr_family nr_genus nr_species eggnog kegg cazy_family cazy_class ardb_class ardb_type card gene/){
		foreach my $std(sort keys %group_info){
			next if (scalar keys %{$group_info{$std}} < 2);
			$queue->addcommond("mkdir -p $outdir/difference/classic/$nn/$std");
			if(scalar keys %{$group_info{$std}} == 2){
				$queue->addcommond("group_significance.py -i $outdir/profile_collection/group/$nn/$std/data.biom -m $outdir/group_info/$std.info -c $std -s kruskal_wallis -o $outdir/difference/classic/$nn/$std/difference.xls &> $outdir/difference/classic/$nn/$std/log.txt");
			}else{
				$queue->addcommond("group_significance.py -i $outdir/profile_collection/group/$nn/$std/data.biom -m $outdir/group_info/$std.info -c $std -s kruskal_wallis -o $outdir/difference/classic/$nn/$std/difference.xls &> $outdir/difference/classic/$nn/$std/log.txt");
			}
		}
	}
	$queue->run();
	
	$queue->set_job_cpu(10);
	$queue->set_job_mem(10);
	$queue->set_job_name("物种与功能差异分析:组建LEfSe差异判别分析-taxonomy");
	$queue->set_work_dir($workdir);
	$queue->addcommond("mkdir -p $outdir/difference/lefse");
	foreach my $std(sort keys %group_info){
		$queue->addcommond("mkdir -p $outdir/difference/lefse/taxonomy/$std");
		$queue->addcommond("cd $outdir/difference/lefse/taxonomy/$std");
		$queue->addcommond("$script/lefse.profile_for_biom.pl $outdir/profile_collection/group/taxonomy/$std/data.txt data");
		$queue->addcommond("biom convert -i data.lefse.profile_for_biom.xls -o data.lefse.profile.biom --table-type \"OTU table\" --process-obs-metadata taxonomy --to-hdf5");
		$queue->addcommond("$script/lefse.config.pl data.lefse.profile.biom $outdir/group_info/$std.txt 6 lefse.config");
		$queue->addcommond("$script/analysis_lefse.py -i data.lefse.profile.biom -c lefse.config -m $outdir/group_info/$std.txt -o tax_lefse");
		$queue->addcommond("cd -");
	}
	$queue->run();
	
	$queue->set_job_cpu(10);
	$queue->set_job_mem(10);
	$queue->set_job_name("物种与功能差异分析:组建LEfSe差异判别分析-kegg");
	$queue->set_work_dir($workdir);
	$queue->addcommond("mkdir -p $outdir/difference/lefse");
	foreach my $std(sort keys %group_info){
		$queue->addcommond("mkdir -p $outdir/difference/lefse/kegg/$std");
		$queue->addcommond("cd $outdir/difference/lefse/kegg/$std");
		$queue->addcommond("$script/lefse.profile_for_biom.pathway.pl $outdir/gene_annotation/kegg/outdir/kegg.Pathway.profile.xls data.lefse.profile_for_biom.xls");
		$queue->addcommond("biom convert -i data.lefse.profile_for_biom.xls -o data.lefse.profile.biom --table-type \"OTU table\" --process-obs-metadata taxonomy --to-hdf5");
		$queue->addcommond("$script/lefse.config.pl data.lefse.profile.biom $outdir/group_info/$std.txt 3 lefse.config");
		$queue->addcommond("$script/analysis_lefse.py -i data.lefse.profile.biom -c lefse.config -m $outdir/group_info/$std.txt -o pathway_lefse");
		$queue->addcommond("cd -");
	}
	$queue->run();
	
	$queue->set_job_cpu(10);
	$queue->set_job_mem(10);
	$queue->set_job_name("物种与功能差异分析:组建LEfSe差异判别分析-eggnog");
	$queue->set_work_dir($workdir);
	$queue->addcommond("mkdir -p $outdir/difference/lefse");
	foreach my $std(sort keys %group_info){
		$queue->addcommond("mkdir -p $outdir/difference/lefse/eggnog/$std");
		$queue->addcommond("cd $outdir/difference/lefse/eggnog/$std");
		$queue->addcommond("$script/lefse.profile_for_biom.cogFunction.pl $outdir/gene_annotation/eggnog/outdir/eggnog.function.xls data.lefse.profile_for_biom.xls");
		$queue->addcommond("biom convert -i data.lefse.profile_for_biom.xls -o data.lefse.profile.biom --table-type \"OTU table\" --process-obs-metadata taxonomy --to-hdf5");
		$queue->addcommond("$script/lefse.config.pl data.lefse.profile.biom $outdir/group_info/$std.txt 2 lefse.config");
		$queue->addcommond("$script/analysis_lefse.py -i data.lefse.profile.biom -c lefse.config -m $outdir/group_info/$std.txt -o cog_lefse");
		$queue->addcommond("cd -");
	}
	$queue->run();
	
	$queue->set_job_cpu(10);
	$queue->set_job_mem(10);
	$queue->set_job_name("物种与功能差异分析:组间STAMP差异分析");
	$queue->set_work_dir($workdir);
	$queue->addcommond("mkdir -p $outdir/difference/stamp");
	foreach my $nn(qw/nr_phylum nr_class nr_order nr_family nr_genus nr_species eggnog kegg cazy_family cazy_class ardb_class ardb_type card gene/){
		foreach my $std(sort keys %group_info){
			$queue->addcommond("mkdir -p $outdir/difference/stamp/$nn/$std");
			$queue->addcommond("$script/file_for_stamp.pl $outdir/group_info/$std.info $outdir/profile_collection/group/$nn/$std/data.txt  500 $outdir/difference/stamp/$nn/$std/data");
		}
	}
	$queue->run();

	# #### wilcoxon ####
	# $queue->addcommond("mkdir -p $outdir/Differential_analysis/wilcoxon/gene");
	# $queue->addcommond("cd $outdir/Differential_analysis/wilcoxon/gene");
	# $queue->addcommond("$script/wilcox.py -i $outdir/gene_profile/gene_profile.base_percent.txt -o gene.wilcoxon.result.xls -g NULL1 -c NULL2 ");
	# $queue->addcommond("cd -");
	# $queue->addcommond("sed -i 's/\\\"//g' $outdir/Differential_analysis/wilcoxon/gene/gene.wilcoxon.result.xls\n");
	# $queue->addcommond("sed -i 's/\\\"//g' $outdir/Differential_analysis/wilcoxon/gene/gene.wilcoxon.result.xls.filtered.xls\n");

	# #### MetagenomeSeq ####
	# $queue->addcommond("mkdir -p $outdir/Differential_analysis/MetagenomeSeq_diff/gene");
	# $queue->addcommond("biom convert -i $outdir/gene_profile/gene_profile.base_percent.txt -o $outdir/Differential_analysis/MetagenomeSeq_diff/gene/table_gene.biom --table-type \"OTU table\" --to-hdf5\n");
	# $queue->addcommond("$script/metagenomeSeq.pl -i $outdir/Differential_analysis/MetagenomeSeq_diff/gene/table_gene.biom -group NULL1 -group_name NULL2 -group_name_A NULL3 -group_name_B NULL4 -o $outdir/MetagenomeSeq_diff/MetagenomeSeq_diff/gene/gene\n");

    #biom convert -i otu.biom -o otu.txt --to-tsv --header-key Taxonomy --output-metadata-id 'Consensus Lineage'

	#### lefse  ####
	#### ternary ###
	#### circos ###
	#### Enterotyping ####
	### Gene_Taxon ####
	#
	$queue->wait();
}

sub result_collection{
	$queue->set_job_cpu(10);
	$queue->set_job_mem(10);
	$queue->set_job_name("结果整理");
	$queue->set_work_dir($workdir);

	$queue->addcommond("mkdir -p $outdir/Result/results/");
	$queue->addcommond("mkdir -p $outdir/Result/results/00.Data");
	$queue->addcommond("cp $outdir/fastq_QC/rawdata.stat.xls  $outdir/Result/results/00.Data/");
	$queue->addcommond("cp $outdir/fastq_QC/cleandata.stat.xls $outdir/Result/results/00.Data/cleandata.nohost.stat.xls");
	$queue->addcommond("cp $outdir/fastq_QC/*fastqc.zip  $outdir/Result/results/00.Data/");

	$queue->addcommond("mkdir -p $outdir/Result/results/01.Assembly");
	$queue->addcommond("mkdir -p $outdir/Result/results/02.GenePredict");
	$queue->addcommond("cp $outdir/assembly/assembly_stat.xls     $outdir/Result/results/01.Assembly/");
	$queue->addcommond("cp $outdir/gene_predict/genePredict_stat.xls  $outdir/Result/results/02.GenePredict/");

	foreach my $sn(sort keys %samples){
	     $queue->addcommond("mkdir  -p $outdir/Result/results/01.Assembly/$sn");
	     $queue->addcommond("gzip -c $outdir/assembly/$sn/filtered_contig.fa > $outdir/Result/results/01.Assembly/$sn/filtered_contig.fa.gz");
	     $queue->addcommond("cp $outdir/assembly/$sn/contig_stat.txt  $outdir/Result/results/01.Assembly/$sn/");
	     $queue->addcommond("mkdir  -p $outdir/Result/results/02.GenePredict/$sn");
             $queue->addcommond("gzip -c $outdir/gene_predict/seq.fna > $outdir/Result/results/02.GenePredict/$sn/seq.fna.gz");
             $queue->addcommond("gzip -c $outdir/gene_predict/seq.faa > $outdir/Result/results/02.GenePredict/$sn/seq.faa.gz");
	     $queue->addcommond("gzip -c $outdir/gene_predict/$sn/seq_more100.fna > $outdir/Result/results/02.GenePredict/$sn/seq_more100.fna.gz");
	     $queue->addcommond("cp $outdir/gene_predict/$sn/seq.gff > $outdir/Result/results/02.GenePredict/$sn/");
	}

        $queue->addcommond("mkdir  -p $outdir/Result/results/03.GeneSet");
        $queue->addcommond("cp $outdir/unigene_set/geneset_stat.xls  $outdir/Result/results/03.GeneSet/");
        $queue->addcommond("gzip -c $outdir/unigene_set/uniGeneSet.faa > $outdir/Result/results/03.GeneSet/uniGeneSet.faa.gz");
        $queue->addcommond("gzip -c $outdir/unigene_set/uniGeneSet.fna > $outdir/Result/results/03.GeneSet/uniGeneSet.fna.gz");

        $queue->addcommond("mkdir  -p $outdir/Result/results/04.GeneProfile");
        $queue->addcommond("gzip -c $outdir/gene_profile/gene_profile.TPM.txt > $outdir/Result/results/04.GeneProfile/gene_profile.TPM.txt.gz");
	$queue->addcommond("gzip -c $outdir/gene_profile/gene_profile.reads_number.txt > $outdir/Result/results/04.GeneProfile/gene_profile.reads_number.txt.gz");

	$queue->addcommond("mkdir -p $outdir/Result/results/05.Annotation");
	$queue->addcommond("mkdir -p $outdir/Result/results/05.Annotation/ardb");
	$queue->addcommond("mkdir -p $outdir/Result/results/05.Annotation/card");
	$queue->addcommond("mkdir -p $outdir/Result/results/05.Annotation/eggnog");
	$queue->addcommond("mkdir -p $outdir/Result/results/05.Annotation/cazy");
	$queue->addcommond("mkdir -p $outdir/Result/results/05.Annotation/kegg");
	$queue->addcommond("mkdir -p $outdir/Result/results/05.Annotation/taxon");
	#$queue->addcommond("mkdir -p $outdir/Result/results/05.Annotation/mvirdb");
	#$queue->addcommond("mkdir -p $outdir/Result/results/05.Annotation/phi");
	#$queue->addcommond("mkdir -p $outdir/Result/results/05.Annotation/qs");
	#$queue->addcommond("mkdir -p $outdir/Result/results/05.Annotation/tcdb");
	$queue->addcommond("mkdir -p $outdir/Result/results/05.Annotation/vfdb");

	$queue->addcommond("cp -r $outdir/gene_annotation/ardb/outdir   $outdir/Result/results/05.Annotation/ardb/");
	$queue->addcommond("cp -r $outdir/gene_annotation/card/outdir   $outdir/Result/results/05.Annotation/card/");
	$queue->addcommond("cp -r $outdir/gene_annotation/eggnog/outdir   $outdir/Result/results/05.Annotation/eggnog/");
	$queue->addcommond("cp -r $outdir/gene_annotation/cazy/outdir   $outdir/Result/results/05.Annotation/cazy/");
	$queue->addcommond("cp -r $outdir/gene_annotation/kegg/outdir   $outdir/Result/results/05.Annotation/kegg/");
	$queue->addcommond("cp -r $outdir/gene_annotation/nr/outdir   $outdir/Result/results/05.Annotation/taxon/");

	#$queue->addcommond("cp -r $outdir/gene_annotation/mvirdb/outdir   $outdir/Result/results/05.Annotation/mvirdb/");
	#$queue->addcommond("cp -r $outdir/gene_annotation/phi/outdir   $outdir/Result/results/05.Annotation/phi/");
	#$queue->addcommond("cp -r $outdir/gene_annotation/qs/outdir   $outdir/Result/results/05.Annotation/qs/");
	#$queue->addcommond("cp -r $outdir/gene_annotation/tcdb/outdir   $outdir/Result/results/05.Annotation/tcdb/");
	
	$queue->addcommond("cp -r $outdir/gene_annotation/vfdb/outdir   $outdir/Result/results/05.Annotation/vfdb/");
	$queue->addcommond("rm -rf $outdir/Result/results/05.Annotation/*/outdir/*tmp");
	$queue->addcommond("rm -rf $outdir/Result/results/05.Annotation/*/outdir/*r");
	$queue->addcommond("rm -rf $outdir/Result/results/05.Annotation/*/outdir/*r");
	$queue->addcommond("rm -rf $outdir/Result/results/05.Annotation/*/outdir/*blast.txt");

	$queue->addcommond("mkdir -p $outdir/Result/results/06.Composition");
	$queue->addcommond("mkdir -p $outdir/Result/results/06.Composition/Bubble");
	$queue->addcommond("mkdir -p $outdir/Result/results/06.Composition/Bubble/ardb/");
	$queue->addcommond("mkdir -p $outdir/Result/results/06.Composition/Bubble/cazy/");
	$queue->addcommond("mkdir -p $outdir/Result/results/06.Composition/Bubble/taxon/");
	$queue->addcommond("mkdir -p $outdir/Result/results/06.Composition/Heatmap");
	$queue->addcommond("mkdir -p $outdir/Result/results/06.Composition/Heatmap/ardb/");
	$queue->addcommond("mkdir -p $outdir/Result/results/06.Composition/Heatmap/cazy/");
	$queue->addcommond("mkdir -p $outdir/Result/results/06.Composition/Heatmap/taxon/");
	$queue->addcommond("mkdir -p $outdir/Result/results/06.Composition/Megan");
	$queue->addcommond("cp -r $outdir/composition/bubble/ardb_*  $outdir/Result/results/06.Composition/Bubble/ardb/");
	$queue->addcommond("cp -r $outdir/composition/bubble/nr_*   $outdir/Result/results/06.Composition/Bubble/taxon/");
	$queue->addcommond("cp -r $outdir/composition/bubble/cazy_* $outdir/Result/results/06.Composition/Bubble/cazy/");
	$queue->addcommond("cp -r $outdir/composition/bubble/card  $outdir/Result/results/06.Composition/Bubble/");
	$queue->addcommond("cp -r $outdir/composition/bubble/kegg  $outdir/Result/results/06.Composition/Bubble/");
	$queue->addcommond("cp -r $outdir/composition/bubble/eggnog  $outdir/Result/results/06.Composition/Bubble/");
	$queue->addcommond("cp -r $outdir/composition/bubble/gene  $outdir/Result/results/06.Composition/Bubble/");
	$queue->addcommond("cp -r $outdir/composition/heatmap/ardb_*  $outdir/Result/results/06.Composition/Heatmap/ardb/");
	$queue->addcommond("cp -r $outdir/composition/heatmap/nr_*   $outdir/Result/results/06.Composition/Bubble/taxon/");
	$queue->addcommond("cp -r $outdir/composition/heatmap/cazy_* $outdir/Result/results/06.Composition/Heatmap/cazy/");
	$queue->addcommond("cp -r $outdir/composition/heatmap/card  $outdir/Result/results/06.Composition/Heatmap/");
	$queue->addcommond("cp -r $outdir/composition/heatmap/kegg  $outdir/Result/results/06.Composition/Heatmap/");
	$queue->addcommond("cp -r $outdir/composition/heatmap/eggnog  $outdir/Result/results/06.Composition/Heatmap/");
	$queue->addcommond("cp -r $outdir/composition/heatmap/gene  $outdir/Result/results/06.Composition/Heatmap/");
	$queue->addcommond("cp -r $outdir/composition/megan/* $outdir/Result/results/06.Composition/Megan/");
	$queue->addcommond("rm $outdir/Result/results/06.Composition/*/*/*/*r");

	$queue->addcommond("mkdir -p $outdir/Result/results/07.Comparison");
	$queue->addcommond("mkdir -p $outdir/Result/results/07.Comparison/Adonis");
	$queue->addcommond("mkdir -p $outdir/Result/results/07.Comparison/Adonis/ardb");
	$queue->addcommond("mkdir -p $outdir/Result/results/07.Comparison/Adonis/taxon");
	$queue->addcommond("mkdir -p $outdir/Result/results/07.Comparison/Adonis/cazy");
	$queue->addcommond("mkdir -p $outdir/Result/results/07.Comparison/Anosim");
	$queue->addcommond("mkdir -p $outdir/Result/results/07.Comparison/Anosim/ardb");
	$queue->addcommond("mkdir -p $outdir/Result/results/07.Comparison/Anosim/taxon");
	$queue->addcommond("mkdir -p $outdir/Result/results/07.Comparison/Anosim/cazy");
	$queue->addcommond("mkdir -p $outdir/Result/results/07.Comparison/Distance");
	$queue->addcommond("mkdir -p $outdir/Result/results/07.Comparison/Distance/ardb");
	$queue->addcommond("mkdir -p $outdir/Result/results/07.Comparison/Distance/taxon");
	$queue->addcommond("mkdir -p $outdir/Result/results/07.Comparison/Distance/cazy");
	$queue->addcommond("mkdir -p $outdir/Result/results/07.Comparison/Nmds");
	$queue->addcommond("mkdir -p $outdir/Result/results/07.Comparison/Nmds/ardb");
	$queue->addcommond("mkdir -p $outdir/Result/results/07.Comparison/Nmds/taxon");
	$queue->addcommond("mkdir -p $outdir/Result/results/07.Comparison/Nmds/cazy");
	$queue->addcommond("mkdir -p $outdir/Result/results/07.Comparison/Hcluster_tree");
	$queue->addcommond("mkdir -p $outdir/Result/results/07.Comparison/Hcluster_tree/ardb");
	$queue->addcommond("mkdir -p $outdir/Result/results/07.Comparison/Hcluster_tree/taxon");
	$queue->addcommond("mkdir -p $outdir/Result/results/07.Comparison/Hcluster_tree/cazy");
	$queue->addcommond("mkdir -p $outdir/Result/results/07.Comparison/Pca");
	$queue->addcommond("mkdir -p $outdir/Result/results/07.Comparison/Pca/ardb");
	$queue->addcommond("mkdir -p $outdir/Result/results/07.Comparison/Pca/taxon");
	$queue->addcommond("mkdir -p $outdir/Result/results/07.Comparison/Pca/cazy");
	$queue->addcommond("mkdir -p $outdir/Result/results/07.Comparison/Pcoa");
	$queue->addcommond("mkdir -p $outdir/Result/results/07.Comparison/Pcoa/ardb");
	$queue->addcommond("mkdir -p $outdir/Result/results/07.Comparison/Pcoa/taxon");
	$queue->addcommond("mkdir -p $outdir/Result/results/07.Comparison/Pcoa/cazy");
   

	$queue->addcommond("cp -r $outdir/comparison/adonis/ardb_*  $outdir/Result/results/07.Comparison/Adonis/ardb/");
	$queue->addcommond("cp -r $outdir/comparison/adonis/nr_*   $outdir/Result/results/07.Comparison/Adonis/taxon/");
	$queue->addcommond("cp -r $outdir/comparison/adonis/cazy_* $outdir/Result/results/07.Comparison/Adonis/cazy/");
	$queue->addcommond("cp -r $outdir/comparison/adonis/card  $outdir/Result/results/07.Comparison/Adonis/");
	$queue->addcommond("cp -r $outdir/comparison/adonis/kegg  $outdir/Result/results/07.Comparison/Adonis/");
	$queue->addcommond("cp -r $outdir/comparison/adonis/eggnog  $outdir/Result/results/07.Comparison/Adonis/");
	$queue->addcommond("cp -r $outdir/comparison/adonis/gene  $outdir/Result/results/07.Comparison/Adonis/");

	$queue->addcommond("cp -r $outdir/comparison/anosim/ardb_*  $outdir/Result/results/07.Comparison/Anosim/ardb/");
	$queue->addcommond("cp -r $outdir/comparison/anosim/nr_*   $outdir/Result/results/07.Comparison/Anosim/taxon/");
	$queue->addcommond("cp -r $outdir/comparison/anosim/cazy_* $outdir/Result/results/07.Comparison/Anosim/cazy/");
	$queue->addcommond("cp -r $outdir/comparison/anosim/card  $outdir/Result/results/07.Comparison/Anosim/");
	$queue->addcommond("cp -r $outdir/comparison/anosim/kegg  $outdir/Result/results/07.Comparison/Anosim/");
	$queue->addcommond("cp -r $outdir/comparison/anosim/eggnog  $outdir/Result/results/07.Comparison/Anosim/");
	$queue->addcommond("cp -r $outdir/comparison/anosim/gene  $outdir/Result/results/07.Comparison/Anosim/");

	$queue->addcommond("cp -r $outdir/comparison/distance/ardb_*  $outdir/Result/results/07.Comparison/Distance/ardb/");
	$queue->addcommond("cp -r $outdir/comparison/distance/nr_*   $outdir/Result/results/07.Comparison/Distance/taxon/");
	$queue->addcommond("cp -r $outdir/comparison/distance/cazy_* $outdir/Result/results/07.Comparison/Distance/cazy/");
	$queue->addcommond("cp -r $outdir/comparison/distance/card  $outdir/Result/results/07.Comparison/Distance/");
	$queue->addcommond("cp -r $outdir/comparison/distance/kegg  $outdir/Result/results/07.Comparison/Distance/");
	$queue->addcommond("cp -r $outdir/comparison/distance/eggnog  $outdir/Result/results/07.Comparison/Distance/");
	$queue->addcommond("cp -r $outdir/comparison/distance/gene  $outdir/Result/results/07.Comparison/Distance/");

	$queue->addcommond("cp -r $outdir/comparison/nmds/ardb_*  $outdir/Result/results/07.Comparison/Nmds/ardb/");
	$queue->addcommond("cp -r $outdir/comparison/nmds/nr_*   $outdir/Result/results/07.Comparison/Nmds/taxon/");
	$queue->addcommond("cp -r $outdir/comparison/nmds/cazy_* $outdir/Result/results/07.Comparison/Nmds/cazy/");
	$queue->addcommond("cp -r $outdir/comparison/nmds/card  $outdir/Result/results/07.Comparison/Nmds/");
	$queue->addcommond("cp -r $outdir/comparison/nmds/kegg  $outdir/Result/results/07.Comparison/Nmds/");
	$queue->addcommond("cp -r $outdir/comparison/nmds/eggnog  $outdir/Result/results/07.Comparison/Nmds/");
	$queue->addcommond("cp -r $outdir/comparison/nmds/gene  $outdir/Result/results/07.Comparison/Nmds/");

	$queue->addcommond("cp -r $outdir/comparison/hcluster_tree/ardb_*  $outdir/Result/results/07.Comparison/Hcluster_tree/ardb/");
	$queue->addcommond("cp -r $outdir/comparison/hcluster_tree/nr_*   $outdir/Result/results/07.Comparison/Hcluster_tree/taxon/");
	$queue->addcommond("cp -r $outdir/comparison/hcluster_tree/cazy_* $outdir/Result/results/07.Comparison/Hcluster_tree/cazy/");
	$queue->addcommond("cp -r $outdir/comparison/hcluster_tree/card  $outdir/Result/results/07.Comparison/Hcluster_tree/");
	$queue->addcommond("cp -r $outdir/comparison/hcluster_tree/kegg  $outdir/Result/results/07.Comparison/Hcluster_tree/");
	$queue->addcommond("cp -r $outdir/comparison/hcluster_tree/eggnog  $outdir/Result/results/07.Comparison/Hcluster_tree/");
	$queue->addcommond("cp -r $outdir/comparison/hcluster_tree/gene  $outdir/Result/results/07.Comparison/Hcluster_tree/");


	$queue->addcommond("cp -r $outdir/comparison/pca/ardb_*  $outdir/Result/results/07.Comparison/Pca/ardb/");
	$queue->addcommond("cp -r $outdir/comparison/pca/nr_*   $outdir/Result/results/07.Comparison/Pca/taxon/");
	$queue->addcommond("cp -r $outdir/comparison/pca/cazy_* $outdir/Result/results/07.Comparison/Pca/cazy/");
	$queue->addcommond("cp -r $outdir/comparison/pca/card  $outdir/Result/results/07.Comparison/Pca/");
	$queue->addcommond("cp -r $outdir/comparison/pca/kegg  $outdir/Result/results/07.Comparison/Pca/");
	$queue->addcommond("cp -r $outdir/comparison/pca/eggnog  $outdir/Result/results/07.Comparison/Pca/");
	$queue->addcommond("cp -r $outdir/comparison/pca/gene  $outdir/Result/results/07.Comparison/Pca/");


	$queue->addcommond("cp -r $outdir/comparison/pcoa/ardb_*  $outdir/Result/results/07.Comparison/Pcoa/ardb/");
	$queue->addcommond("cp -r $outdir/comparison/pcoa/nr_*   $outdir/Result/results/07.Comparison/Pcoa/taxon/");
	$queue->addcommond("cp -r $outdir/comparison/pcoa/cazy_* $outdir/Result/results/07.Comparison/Pcoa/cazy/");
	$queue->addcommond("cp -r $outdir/comparison/pcoa/card  $outdir/Result/results/07.Comparison/Pcoa/");
	$queue->addcommond("cp -r $outdir/comparison/pcoa/kegg  $outdir/Result/results/07.Comparison/Pcoa/");
	$queue->addcommond("cp -r $outdir/comparison/pcoa/eggnog  $outdir/Result/results/07.Comparison/Pcoa/");
	$queue->addcommond("cp -r $outdir/comparison/pcoa/gene  $outdir/Result/results/07.Comparison/Pcoa/");

	$queue->addcommond("rm  $outdir/Result/results/07.Comparison/*/*/*/*r");
	$queue->addcommond("rm  $outdir/Result/results/07.Comparison/*/*/*/*/*r");

	$queue->addcommond("mkdir -p $outdir/Result/results/08.Difference");
	$queue->addcommond("mkdir -p $outdir/Result/results/08.Difference/LEfSe");
	$queue->addcommond("mkdir -p $outdir/Result/results/08.Difference/LEfSe/taxon/");
	$queue->addcommond("mkdir -p $outdir/Result/results/08.Difference/LEfSe/eggnog/");
	$queue->addcommond("mkdir -p $outdir/Result/results/08.Difference/LEfSe/kegg/");
	$queue->addcommond("mkdir -p $outdir/Result/results/08.Difference/Wilcoxon");
	$queue->addcommond("mkdir -p $outdir/Result/results/08.Difference/Wilcoxon/ardb/");
	$queue->addcommond("mkdir -p $outdir/Result/results/08.Difference/Wilcoxon/taxon/");
	$queue->addcommond("mkdir -p $outdir/Result/results/08.Difference/Wilcoxon/cazy/");

	$queue->addcommond("cp -r $outdir/difference/lefse/eggnog/group  $outdir/Result/results/08.Difference/LEfSe/eggnog/");
	$queue->addcommond("cp -r $outdir/difference/lefse/kegg/group  $outdir/Result/results/08.Difference/LEfSe/kegg/");
	$queue->addcommond("cp -r $outdir/difference/lefse/taxonomy/group  $outdir/Result/results/08.Difference/LEfSe/taxon/");

	$queue->addcommond("cp -r $outdir/difference/classic/ardb_*  $outdir/Result/results/08.Difference/Wilcoxon/ardb/");
	$queue->addcommond("cp -r $outdir/difference/classic/nr_*   $outdir/Result/results/08.Difference/Wilcoxon/taxon/");
	$queue->addcommond("cp -r $outdir/difference/classic/cazy_* $outdir/Result/results/08.Difference/Wilcoxon/cazy/");
	$queue->addcommond("cp -r $outdir/difference/classic/card  $outdir/Result/results/08.Difference/Wilcoxon/");
	$queue->addcommond("cp -r $outdir/difference/classic/kegg  $outdir/Result/results/08.Difference/Wilcoxon/");
	$queue->addcommond("cp -r $outdir/difference/classic/eggnog  $outdir/Result/results/08.Difference/Wilcoxon/");
	$queue->addcommond("cp -r $outdir/difference/classic/gene  $outdir/Result/results/08.Difference/Wilcoxon/");
        $queue->addcommond("perl $script/pdf2png.pl  $outdir/Result/results/");
	$queue->run();
	$queue->wait();
}

sub html_report{
	$queue->set_job_cpu(10);
	$queue->set_job_mem(10);
	$queue->set_job_name("网页报告整理");
	$queue->set_work_dir($workdir);

	$queue->addcommond("cp -r /home/pub/pipeline/metagenomics/html/report.html  $outdir/Result/");
	$queue->addcommond("mkdir -p $outdir/Result/web/res");
        $queue->addcommond("cp -r /home/pub/pipeline/metagenomics/html/web/css $outdir/Result/web/");
        $queue->addcommond("cp -r /home/pub/pipeline/metagenomics/html/web/js $outdir/Result/web/");
        $queue->addcommond("cp /home/pub/pipeline/metagenomics/html/web/res/*html $outdir/Result/web/res/");
	$queue->addcommond("mkdir -p $outdir/Result/web/res/00.Data");
	$queue->addcommond("perl $script/data_html.pl $outdir/Result/results/00.Data/  $outdir/Result/web/res/");

	$queue->addcommond("mkdir -p $outdir/Result/web/res/01.Assembly");
	$queue->addcommond("perl $script/assembly_html.pl $outdir/Result/results/01.Assembly/  $outdir/Result/web/res/");

	$queue->addcommond("mkdir -p $outdir/Result/web/res/02.GenePredict");
	$queue->addcommond("perl $script/genepredict_html.pl $outdir/Result/results/02.GenePredict/  $outdir/Result/web/res/");

	$queue->addcommond("mkdir -p $outdir/Result/web/res/03.GeneSet");
	$queue->addcommond("perl $script/geneset_html.pl $outdir/Result/results/03.GeneSet/  $outdir/Result/web/res/");
	
	$queue->addcommond("mkdir -p $outdir/Result/web/res/04.GeneProfile");
	$queue->addcommond("perl $script/geneprofile_html.pl $outdir/gene_profile/  $outdir/Result/web/res/");

	$queue->addcommond("mkdir -p $outdir/Result/web/res/05.Annotation");
	$queue->addcommond("perl $script/annotation_html.pl $outdir/Result/results/05.Annotation/  $outdir/Result/web/res/");

	$queue->addcommond("mkdir -p $outdir/Result/web/res/06.Composition");
	$queue->addcommond("perl $script/composition_html.pl $outdir/Result/results/06.Composition/ $outdir/Result/web/res/");

	$queue->addcommond("mkdir -p $outdir/Result/web/res/07.Comparison");
	$queue->addcommond("perl $script/comparison_html.pl $outdir/Result/results/07.Comparison/  $outdir/Result/web/res/");

	$queue->addcommond("mkdir -p $outdir/Result/web/res/08.Difference");
	$queue->addcommond("perl $script/difference_html.pl $outdir/Result/results/08.Difference/Wilcoxon  $outdir/Result/web/res/");
	$queue->addcommond("perl $script/lefse_html.pl $outdir/Result/results/08.Difference/LEfSe  $outdir/Result/web/res/");
	$queue->run();
	$queue->wait();
}

__END__
