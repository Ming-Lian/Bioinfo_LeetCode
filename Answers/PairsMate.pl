#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;

# 脚本帮助文档
=head1 Description

Thise script is used to match the pair-end reads from one sample

=head1 Usage

$0 -1 <input1.fastq> -2 <input2.fastq> [-o <outdir>]

=head1 Parameters

-1	[str]	Input 1st-end fastq file
-2	[str]	Input 2nd-end fastq file
-o	[str]	Output directory [default: current folder]

=cut

my ($Input1,$Input2,$Outdir);
GetOptions(
"1:s"=>\$Input1,
"2:s"=>\$Input2,
"o:s"=>\$Outdir
);

$Outdir=`pwd` unless (defined($Outdir));
die `pod2text $0` if (!$Input1) or (!$Input2);

open FQ1,"<$Input1" or die "Cann't open $Input1\n";
open FQ2,"<$Input2" or die "Cann't open $Input2\n";

my $Input1_basename=`basename $Input1`;
chomp $Input1_basename;
my $Input2_basename=`basename $Input2`;
chomp $Input2_basename;
open OUT1,">$Outdir/${Input1_basename}.match" or die "Cann't create $Outdir/${Input1_basename}.match\n";
open OUT2,">$Outdir/${Input2_basename}.match" or die "Cann't create $Outdir/${Input2_basename}.match\n";

# 载入两个fastq文件，保存成哈希
my (%fq1_seq,%fq1_qua,%fq2_seq,%fq2_qua);
while (!eof(FQ1)){
        my ($id,)=split(/\s/,<FQ1>);
        $fq1_seq{$id}=<FQ1>;
        <FQ1>;
        $fq1_qua{$id}=<FQ1>;
        chomp ($fq1_seq{$id},$fq1_qua{$id});
}
while (!eof(FQ2)){
        my ($id,)=split(/\s/,<FQ2>);
        $fq2_seq{$id}=<FQ2>;
        <FQ2>;
        $fq2_qua{$id}=<FQ2>;
        chomp ($fq2_seq{$id},$fq2_qua{$id});
}
close FQ1;
close FQ2;

# 双端配对
foreach my $key (sort keys %fq1_seq){
        if(defined($fq2_seq{$key})){
                print OUT1 "$key\n$fq1_seq{$key}\n+\n$fq1_qua{$key}\n";
                print OUT2 "$key\n$fq2_seq{$key}\n+\n$fq2_qua{$key}\n";
        }
}
close OUT1;
close OUT2;
