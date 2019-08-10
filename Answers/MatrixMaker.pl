#!/usr/bin/perl -w
use Getopt::Long;
use POSIX;

## ********************
## ** 获取并解析参数 **
## ********************

my %opts=();
GetOptions(\%opts,"i:s","o:s", "suffix:s");

if (!$opts{i} or !$opts{o} ){
	print "
    Description:
    
    This script is used to merge serveral quant files into a a matrix
    
	USAGE: perl $0
		-i input dir
		-o output file
        -suffix the suffix of input file in input dir\n\n";
	exit;
}

## **************
## ** 载入数据 **
## **************
print STDERR "\tread dir and load quantitative data ... \n";

# 获取指定目录下所有文件的文件名
$indir = $opts{i};
opendir DIR,$indir or die "can't opendir $indir: $!";
@files = readdir(DIR);
closedir DIR;

# 只对指定后缀的文件进行处理
$file_count = 0;
%Hash_Feat2Sample2Quant = ();
foreach $file (@files){
    if($file =~ /$opts{suffix}$/){
        $file_count++;
        $sampleId = `basename $file $opts{suffix}`; # 从文件名中获取样本Id
        chomp $sampleId;
        push @SampleList, $sampleId;
        
        print STDERR "\t\t$file_count: $file\n";
        open IN,"<$indir/$file" or die;
        while(<IN>){
            chomp;
            @row = split /\t/;
            $Hash_Feat2Sample2Quant{$row[0]}{$sampleId} = $row[1];
        }
        close IN;
    }
}

print STDERR "\tin total $file_count files loaded; \n\n";

## **************
## ** 写出数据 **
## **************

print STDERR "\tnow ready to generate output ... \n";

## --- open OUT file for output --
open OUT, ">$opts{o}" or die;
print OUT join("\t", "GENE", @SampleList), "\n"; ## 写入文件的表头

my $feat_count = 0;
while( my ( $feat, $hashref ) = each %Hash_Feat2Sample2Quant ){
	$feat_total ++;
    
    my @current_row;
    foreach $sample (@SampleList){
        push(@current_row, exists $$hashref{$sample}?$$hashref{ $sample } : 0);
    }
    
    print OUT join("\t", $feat, @current_row ), "\n";
    
    if($feat_total%10000==0){
        print STDERR "\t\t$feat_total features have been written into outfile\n";
    }
}
close OUT;

print STDERR "\tall jobs done; \n\n";
