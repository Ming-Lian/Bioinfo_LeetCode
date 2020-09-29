#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;

# 脚本帮助文档
=head1 Description

	Thise script is used to extract a number of fastq recorders from original input fastq file

=head1 Usage

	$0 -n <totalReads> -e <base-pairs to extract> -l <readLength> -1 <input1.fastq> -2 <input2.fastq> [-o <outdir>]

=head1 Parameters

	-n	[int]	Total reads number of input.fastq
	-e	[int]	Base-pairs per end you want to extract
	-l	[int]	The reads length
	-1	[str]	Input 1st-end fastq file
	-2	[str]	Input 2nd-end fastq file
	-o	[str]	Output directory [default: current folder]

=cut

my ($totalReads,$bpNum,$length,$Input1,$Input2,$Outdir);
GetOptions(
	"n:i"=>\$totalReads,
	"e:i"=>\$bpNum,
	"l:i"=>\$length,
	"1:s"=>\$Input1,
	"2:s"=>\$Input2,
	"o:s"=>\$Outdir
	);

$Outdir=`pwd` unless (defined($Outdir));
die `pod2text $0` if ((!$totalReads) or (!$bpNum)) or (!$length) or (!$Input1) or (!$Input2);

open FQ1,"<$Input1" or die "$!\n";
open FQ2,"<$Input2" or die "$!\n";

my $Input1_basename=`basename $Input1`;
chomp $Input1_basename;
my $Input2_basename=`basename $Input2`;
chomp $Input2_basename;
open OUT1,">$Outdir/${Input1_basename}.extract" or die "$!\n";
open OUT2,">$Outdir/${Input2_basename}.extract" or die "$!\n";

my $readsRemain=$bpNum/$length;

my $remainCount=0;

while(! eof($FQ1)){
	# 读入1st-end fastq 文件的四行
	my $fq1_1=<FQ1>;
	my $fq1_2=<FQ1>;
	my $fq1_3=<FQ1>;
	my $fq1_4=<FQ1>;
	chomp($fq1_1,$fq1_2,$fq1_3,$fq1_4);
	# 读入2nd-end fastq 文件的四行
	my $fq2_1=<FQ2>;
	my $fq2_2=<FQ2>;
	my $fq2_3=<FQ2>;
	my $fq2_4=<FQ2>;
	chomp($fq2_1,$fq2_2,$fq2_3,$fq2_4);

	# 随机抽取
	if (rand()<$readsRemain/$totalReads){
		$remainCount++;
		print OUT1 "$fq1_1\n$fq1_2\n$fq1_3\n$fq1_4\n";
		print OUT2 "$fq2_1\n$fq2_2\n$fq2_3\n$fq2_4\n";
	}
}

print "Total reads: $totalReads\n";
print "Theorical remained reads: $readsRemain\n";
print "Practical remained reads: $remainCount\n";

close FQ1;
close FQ2;
close OUT1;
close OUT2;
