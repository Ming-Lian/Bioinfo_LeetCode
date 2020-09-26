#!/usr/bin/perl
use strict;
use warnings;

sub BWT_index(){
	my ($infile) = @_;
	open IN, "<$infile" or die "$!\n";
	# 解析fasta文件
	my $fa_seq;
	while(<IN>){
		chomp;
		next if(/^>/);
		$fa_seq .= $_;
	}
	# 将以字符串形式组织的序列，打散成数组，并在末尾追加终止标识符
	my @seq = split //,$fa_seq;
	push @seq, "\$";
	# 循环出队入队，得到BWT matrix
	my @BWT;
	my $tail = pop @seq;
	while($seq[$#seq] ne "\$"){
		unshift @seq, $tail;
		push @BWT, join("",@seq);
		$tail = pop @seq;
	}
	# 对@BWT中的各个字符串元素按照字母顺序逐一进行遍历，截取最后一个字母
	my @LastColumn;
	foreach my $string (sort {$a cmp $b} @BWT){
		push @LastColumn, substr($string, length($string)-1);
	}
	return @LastColumn;
}

my $infile = $ARGV[0];
my @index = &BWT_index($infile);
print "@index\n";
