#!/usr/bin/perl -w

use strict;


#Elzo de Wit
#Netherlands Cancer Institute
#e.d.wit@nki.nl


#generate a map of the fragments that are used in the 4C analysis
#In 4C two restriction enzymes are used: RE1 and RE2. Using an inverse
#PCR DNA fragments are amplified that are ligated to a region of interest: 
#the bait or viewpoint. The ligated fragments and the viewpoint are 
#fragments that are between the first and second restriction enzyme.
#In the mapping pipeline the position of these fragments is used to
#create a mapping file.

#Note: in the 4C analysis the following fragments can be ligated to
#viewpoint:
#1: RE1-RE2 (or RE2-RE1)
#2: RE1-RE1 (blind fragments)

#The blind fragments have a lower performance compared to the
#non-blind fragments which is why we leave them out of the
#fragment map (and therefore the downstream analysis).
#Future versions may include an option to include blind fragments,
#although this impacts both the repeat finding and the mapping
#script, because you would need to discriminate blind from
#non-blind
#This script only generates fragment that are formed
#between the first and the second restriction enzyme.

my $fasta = shift @ARGV or usage(); #genome fasta file  ##@ARGV 传递脚本的命令行参数列表
my $re1 = shift @ARGV or usage(); #first restriction enzyme
my $re2 = shift @ARGV or usage(); #second restriction enzyme
my $output_dir = shift @ARGV or usage(); #where to store the fragment map

$|=1;  ##不进行缓存，立即写

#create output directory with the fragment map
mkdir $output_dir;
my $chrom;
my $seq;
open FASTA, $fasta or die "Cannot open file: $!";
my $pos = 1;
my $lines = 0;
my $add_end;
my %seen;
while(<FASTA>){
	chomp;
	#fasta header
	if(/^>(.*)/){
		$chrom = $1;   #捕获匹配的第一个括号里的内容 即chr1 
		#if there are whitespace characters in the header 
		#remove everything after the whitespace
		$chrom =~ s/\s(.*)//;
		if( defined $seen{$chrom} ){
			print "$chrom has been seen before in this fasta file. Exiting\n";
			exit;
		}
		$seen{$chrom} = 1;
		#move into a new sequence and file for creating fragments for
		$seq = "";
		$pos = 1;
		#create a new file to hold the positions of the fragments
		open OUT, ">$output_dir/$chrom.txt" or die "Cannot create file: $!";
		$add_end = 0;
	}else{
		if(/^N+$/){
			$add_end += length($_);  ##length() 计算字符串的长度
		}else{
			$seq .= uc $_;   ##uc() 将所有字符都转换成大写
		}	
		while($seq =~ /($re1.*?$re1)/g){     ## $re1 第一个酶切序列
			my ($sub_seq) = $seq =~ /.*?($re1.*?$re1)/;
			my $end_pos = $pos + length($&) - 1 + $add_end;  ## $& 前一次成功模式匹配的字符串
			if(length( $seq ) < 5e3){		
				#5' fragment
				if($sub_seq =~ /^($re1.*?$re2)/ ){
					my $sub_pos = $pos + length($1) - 1;
					print OUT join("\t", ($pos,$sub_pos,5)), "\n";
				}
				#3' fragment
				if($sub_seq =~ /.*($re2.*?$re1$)/){
					my $sub_pos = $end_pos - length($1) + 1;
					print OUT join("\t", ($sub_pos,$end_pos,3)), "\n";
				}
			}	
			#substitute the entire RE1 fragment to the sequence of RE1
			$seq =~ s/.*?$re1.*?$re1/$re1/;
			$pos = $end_pos - length($re1) + 1; #add one to get the GATC start position
			$add_end = 0;
		}
	}
}	
		

sub usage{
	print "generate_fragment_map.pl genome_fasta re1 re2 output_directory\n";
	exit;
}	
