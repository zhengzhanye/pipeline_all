#!/usr/bin/perl -w

#Elzo de Wit
#Netherlands Cancer Institute
#e.d.wit@nki.nl

#4C mapping consists of three steps:
#1. Splitting and trimming of the original fastq files based on a list of primer sequences
#2. Mapping the splitted sequences to the genome
#3. Filtering mapping file

#Results
#1. A directory splitted_seq/ containing a seperate fastq file for every
#   primer in the index file
#2. output/ directory containing a sam file for every splitted sequence
#3. cis_hits/, count_file/ and filtered_file/ various stages of filtering
#   steps with the filtered_file being the most mature form of 4C files.
#   These files can be directly used in a downstream analysis (such as
#   peakC).

#Command line options:
#index_file: file containing the minimal descriptive information for a
#            4C mapping analysis, see README for details
#run_dir: directory for the storage of the results, will be created
#fastq_file: file containing the sequencing of the 4C experiment (and
#            additional experiments, which will be removed in the splitting
#            step
#threads:    how many threads should bwa use?
#frag_map:   a fragment map directory, see generate_fragment_map.pl
#repeat_dir: repeat directory, see getRepeats.pl


#use Inline::C to speed up the calculation
use Inline C;
use strict;

my $index_file = shift @ARGV or usage();
my $run_dir = shift @ARGV or usage();
my $fastq_file = shift @ARGV or usage();
my $threads = shift @ARGV or usage();
my $frag_map = shift @ARGV;
my $repeat_dir = shift @ARGV;

#create the directory to store the results
mkdir $run_dir;
mkdir "$run_dir/splitted_seq";
mkdir "$run_dir/output";
mkdir "$run_dir/cis_hits";
mkdir "$run_dir/count_file";
mkdir "$run_dir/filtered_file" if defined $frag_map;

#split the reads in the fastq file on the fastq file
splitSeqFiles($index_file, $run_dir,$fastq_file);

#start mapping of the individual files
run_bwa( $index_file, $run_dir, $threads );

#parse SAM file so that only the correct sequences are saved
parseSAM( $index_file, $run_dir );

#create files that contain the raw 4C counts
createCountFile( $index_file, $run_dir );

#intersect the count file with the fragment map
#if this variable is defined by the user
#add zeroes for non-covered/non-ligated fragments
#remove repetitive (non-ligated) fragments if repeat file
#is given (note that this is optional, but should be done)
selectFragments( $index_file, $run_dir, $frag_map, $repeat_dir ) if defined $frag_map;



sub splitSeqFiles{
	my ($file,$run_dir,$fastq_file) = @_;

	#create a structure that holds the file handles and the primer sequences

	my %seq_len;
	open PRIMER, $file or die "Cannot open file: $!";
	while(<PRIMER>){
		chomp;
		my ($exp,$primer,$re_site) = (split /\t/)[0,1,4];
		my $out_file = "$run_dir/splitted_seq/$exp.fastq";
		my $fh;
		open $fh, ">$out_file" or die "Cannot create $out_file: $!";
		#store the name of the experiment, the sequence of the primer, the file handle and the length of the restriction site (-1), which is later used
		#in the selectReads function
		push @{$seq_len{$fastq_file}->{length($primer)}}, [$exp,$primer, $fh, length($re_site)];
	}
	#outputs (individual) fastq files
	selectReads($fastq_file, %{$seq_len{$fastq_file}});
	`gzip $run_dir/splitted_seq/*`;
}


sub selectReads{
	my ($in_file, %seq_len) = @_;
	
	if($in_file =~ /gz$/){
		open SEQ, "gunzip -c $in_file |" or die "Cannot open $in_file: $!";
	}else{
		open SEQ, $in_file or die "Cannot open $in_file: $!";
	}	
	my @fastq_line;
	while(! eof SEQ){   ##用 eof（没有圆括弧）在 while (<>) 循环里检查每个文件的文件结束
		my $header = <SEQ>;
		my $sequence = <SEQ>;
		my $plus = <SEQ>;
		my $quality = <SEQ>;
		for my $len ( keys %seq_len ){
			for my $entry( @{$seq_len{$len}} ){
				#note that the maximum difference is set to 1, which means
				#that one mismatch is allowed in the primer, note that this
				#can be increased if necessary                     
				if(return_diff($sequence,$entry->[1], $len, 1)){       ## 1 代表只允许引物中有一个mismatch ,需要时这个值可以改变
					$sequence = substr($sequence, $len-$entry->[3]);
					$quality = substr($quality, $len-$entry->[3]);
					my $fh = $entry->[2];
					print $fh ($header,$sequence,"+\n",$quality);
				}	
			}	
		}
	}	
}


#unused function for creating an input file for CutAdapt
#but cutadapt was to slow for our purposes, so we
#wrote a dedicated splitter
sub create_adapter_file{
	my ($file,$run_dir) = @_;
	my $adapter_file = "$run_dir/AdapterSeq.fa";
	open INDEX, $file or die "Cannot open Index file: $!";
	open FASTA, ">$adapter_file" or die "Cannot create adapter file: $!";
	while(<INDEX>){
		chomp;
		my ($exp,$primer,$re_site) = (split /\t/)[0,1,3];
		$primer =~ s/$re_site$//;
		print FASTA ">$exp\n$primer\n";
	}
	return $adapter_file;
}	


#run BWASW for every splitted fastq file in the output directory
#for a designated number of threads
sub run_bwa{
	my ($file, $run_dir, $threads) = @_;
	open INDEX, $file or die "Cannot open Index file: $!";
	while(<INDEX>){
		chomp;
		my ($exp,$reference) = (split /\t/)[0,2];
		my $fastq = "$run_dir/splitted_seq/$exp" . ".fastq.gz";
		my $sam_out = "$run_dir/output/$exp" . ".sam";
		my $command = "bwa bwasw -t $threads $reference $fastq > $sam_out";
		print $command, "\n";
		system($command);
	}
}	

#parse the resulting SAM file
#select mapped reads that start or end with the restriction site used
#in the initial 3C restriction and ligation step, only select reads
#that unique in the genome
sub parseSAM{
	my ($file, $run_dir) = @_;
	open INDEX, $file or die "Cannot open Index file: $!";
	while(<INDEX>){
		chomp;
		my ($exp,$chrom,$re_site) = (split /\t/)[0,5,3];
		my $grep_stream = "grep -P \"(\\t$re_site|$re_site\\t)\" $run_dir/output/$exp.sam | grep -P \"\\t$chrom\\t\"";
		open CISHITS, "$grep_stream |" or die "Cannot open stream: $!";
		open OUT, ">$run_dir/cis_hits/$exp.sam" or die "Cannot create file: $!";
		while(<CISHITS>){
			chomp;
			my ($ori,$pos,$map_qual,$cigar,$seq) = (split /\t/)[1,3,4,5,9];
			if($ori == 0 and $map_qual > 0 and $cigar =~ /^\d+M/ and $seq =~ /^$re_site/){
				print OUT $_, "\n";
			}elsif(	$ori == 16 and $map_qual > 0 and $cigar =~ /\d+M$/ and $seq =~ /$re_site$/){  ## flags = 16 ，该read其反向互补序列能够比对到参考序列 
				print OUT $_, "\n";
			}
		}
	}
}	

sub createCountFile{
	my ($file, $run_dir) = @_;
	open INDEX, $file or die "Cannot open Index file: $!";
	while(<INDEX>){
		chomp;
		my ($exp,$chrom) = (split /\t/)[0,5];
		open SAM, "$run_dir/cis_hits/$exp.sam" or die "Cannot create file: $!";
		my %counts;
		while(<SAM>){
			chomp;
			my ($ori,$pos,$map_qual,$cigar,$seq) = (split /\t/)[1,3,4,5,9];
			if($ori == 0){
				$counts{$pos}{5}++;
			}elsif($ori == 16	){
				my ($map_length) = $cigar =~ /(\d+)M$/;
				$pos += ($1-1) ; #minus one so the alignment to the restriction site is correct
				$counts{$pos}{3}++;
			}
		}
		open OUT, ">$run_dir/count_file/$exp.sam" or die "Cannot create file: $!";
		for my $pos ( sort { $a <=> $b } keys %counts ){
			for my $flank ( keys %{$counts{$pos}} ){
				print OUT join("\t", ($pos,$flank,$counts{$pos}{$flank})), "\n";
			}
		}	

	}

}	


sub selectFragments{
	my ( $index_file, $run_dir, $frag_map, $repeats ) = @_;
	
	#define a hash containing the positions of the repeat sequences
	open INDEX, $index_file or die "Cannot open Index file: $!";
	while(<INDEX>){
		chomp;
		my ($exp,$chrom) = (split /\t/)[0,5];
		open COUNT, "$run_dir/count_file/$exp.sam" or die "Cannot open count file: $!";
		open OUT, ">$run_dir/filtered_file/$exp.txt" or die "Cannot create filtered file: $!";
		#read the positions of the fragments
		my %frag_map = readFragMap( $frag_map, $chrom );
		#read in the repeats for the chromosome
		my %repeats;
		%repeats = readRepeats( $repeats, $chrom ) if defined $repeats;

		#create a structure with the counts from the 4C experiment
		my %counts;
		while(<COUNT>){
			chomp;
			my ($pos, $count ) = (split /\t/)[0,2];
			$counts{$pos} = $count;
		}
		for my $pos ( sort { $frag_map{$a} <=> $frag_map{$b} } keys %frag_map ){
			#skip the fragment if it is a repeated fragment
			next if defined $repeats{$pos};
			if(defined $counts{$pos}){
				print OUT join("\t", ($pos, $counts{$pos})), "\n";
			}else{
				print OUT join("\t", ($pos, 0)), "\n";
			}
		}
	}
}	



#read in the positions of the repeat fragments
sub readRepeats{
	my ($repeat_dir, $chrom) = @_;
	my %repeats;
	print "$repeat_dir/$chrom.txt\n";
	open REPEAT, "$repeat_dir/$chrom.txt" or die "Cannot open repeat file: $!";
	while(<REPEAT>){
		chomp;
		$repeats{$_} = 1;
	}
	return %repeats;
}	

#read in the positions of the fragments ends
sub readFragMap{
	my ( $dir, $chrom ) = @_;
	my %frag_map;
	open FRAG, "$dir/$chrom\.txt" or die "Cannot open fragment map for $chrom: $!";
	while(<FRAG>){
		chomp;
		my ( $frag_start, $frag_end, $ori ) = split /\t/;
		if($ori == 5){
			$frag_map{$frag_start} = int(($frag_start+$frag_end)/2);
		}else{
			$frag_map{$frag_end} = int(($frag_start+$frag_end)/2);
		}
	}
	return %frag_map;
}	




				

				

sub usage{
	print "mapping_pipeline.pl index_file out_dir fastq number_of_threads [fragment_map repeat_file]\n";
	exit;
}	


__END__
__C__
#include <stdio.h>

int return_diff( char* str1, char* str2, int len, int max_mm  ) {
	int i;
	int cnt_mm = 0;
	for(i = 0; i < len; i++){
		if(str1[i]!=str2[i]){
			cnt_mm++;
		}
		if(cnt_mm > max_mm){
			return 0;
		}
	}
	return 1;
}	


