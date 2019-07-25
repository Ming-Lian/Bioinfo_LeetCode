#!/usr/bin/perl
use strict;
use warnings;

# 使用方法：
# perl extractSeqFromFasta.pl <in.fasta> <gene list> <out.fasta>

my geneList = $ARGV[0];
my fasta = $ARGV[1];
my out = $ARGV[2];

my %Hash_fasta;
my $seqId;

# 读入fasta文件，存为哈希
open FA,"<$fasta" or die "$!\n";
while(<FA>){
	chomp;
	next if(/^\s?$/);
	if(/^>(.+)$/){
		$seqId = $1;
	}else{
		$Hash_fasta{$seqId} .= $_;
	}
}
close(FA);

# 逐行读入geneList文件，并从上一步的哈希中将该序列提取出来
open LIST,"<$geneList" or die "$!\n";
open OUT,">$out" or die "$!\n";
while(<LIST>){
	chomp;
	next if(/^\s?$/);
	my $gene = $_;
	if($Hash_fasta{$gene}){
		print OUT ">$gene\n$Hash_fasta{$gene}\n";
	}
}
close(LIST);
close(OUT);
