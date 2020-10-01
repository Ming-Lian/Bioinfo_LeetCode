#!/usr/bin/perl
use strict;
use warnings;

my $infile = $ARGV[0];

open IN,"<$infile" or die "$!\n";

my %Contig;
my $seqId;

# 读入fasta文件中的每条contig序列
while(<IN>){
	chomp;
	next if(/^\s?$/);
	if(/^>/){
		$seqId = $_;
	}else{
		$Contig{$seqId} .= $_;
	}
}
close IN;

# 计算每条序列的长度，同时累加每条序列的长度得到总长度N
my %ContigLen;
my $N = 0;
foreach my $seq (keys %Contig){
	my $n = length($Contig{$seq});
	$ContigLen{$seq} = $n;
	$N += $n;
}

# 按照序列长度从大到小遍历每条序列n_i，并逐一累加序列长度L+=n_i
# 一旦L>= 0.5N，则结束遍历，则最后那条被访问的序列的长度即为N50
my $L = 0;
foreach my $seq (sort {$ContigLen{$a} <=> $ContigLen{$b}} keys %ContigLen){
	$L += $ContigLen{$seq};
	if($L >= 0.5*$N){
		print "N50=$ContigLen{$seq}\n";
		last;
	}
}
