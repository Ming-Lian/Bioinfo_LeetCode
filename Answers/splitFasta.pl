#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use POSIX;

# 帮助文档
=head1 Description

	This script is used to split fasta file, which is too large with thosands of sequence

=head1 Usage

	$0 -i <input> -o <output_dir> [-n <seq_num_per_file>] [-m <output_file_num>]
	
=head1 Parameters

	-i	[str]	Input raw fasta file
	-o	[str]	Output file to which directory
	-n	[int]	Sequence number per file, alternate chose paramerter "-n" or "-m", if set "-n" and "-m" at the same time, only take "-n" parameter
	-m	[int]	Output file number (default:100)
=cut

my ($input,$output_dir,$seq_num,$file_num);
GetOptions(
	"i:s"=>\$input,
	"o:s"=>\$output_dir,
	"n:i"=>\$seq_num,
	"m:i"=>\$file_num
	);

die `pod2text $0` if ((!$input) or (!$output_dir));

# 设置每个文件的序列条数
if(!defined($seq_num)){
	if(!defined($file_num)){
		$file_num=100;
		my $total_seq_num=`awk 'BEGIN{n=0} /^>/{n++} END{print n}' $input`;
		chomp $total_seq_num;
		$seq_num=ceil($total_seq_num/$file_num);
	}else{
		my $total_seq_num=`awk 'BEGIN{n=0} /^>/{n++} END{print n}' $input`;
		chomp $total_seq_num;
		$seq_num=ceil($total_seq_num/$file_num);
	}
}

open IN,"<$input" or die "Cann't open $input\n";

my $n_seq=0;	# 该变量用于记录当前扫描到的序列数
my $n_file=1;	# 该变量用于记录当前真正写入的文件的计数
my $input_base=`basename $input`;
chomp $input_base;

open OUT,">$output_dir/${input_base}_${n_file}" or die "Cann't create $output_dir/${input_base}_${n_file}\n";

while(<IN>){
	next if (/^\s+$/);	# 跳过空行
	chomp;
	if (/^>/){
		$n_seq++;
		# 判断目前已经扫描到的序列数，若大于设定的split的序列数，则创建新文件
		if ($n_seq>$seq_num){
			$n_seq=1;
			$n_file++;
			close OUT;
			open OUT,">$output_dir/${input_base}_${n_file}" or die "Cann't create $output_dir/${input_base}_${n_file}\n";
			print OUT "$_\n";
		}else{
			print OUT "$_\n";
		}
	}else{
		print OUT "$_\n";
	}
}

close IN;
