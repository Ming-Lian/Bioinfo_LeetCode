#!/usr/bin/perl

use strict;
use warnings;
use POSIX;

sub TandemRepeat(){
	my ($seqId, $seq, $minUnit, $maxUnit, $minRepeat, $maxRepeat) = @_;
	my @Seq = split //,$$seq;
	my %Queue;
	my %Start;
	my %formerScore;
	my $Head = 0; # 该变量用于限制重复序列起始点，以避免同一个重复序列以倍数单元长度被多次报告
	# 初始化检索起始位点和%Quene
	foreach my $i ($minUnit..$maxUnit){
		$Start{$i} = 0;
		$formerScore{$i} = 0;
		my @a;
		$Queue{$i} = \@a;
	}
	foreach my $i (0..$#Seq){
		# 对当前单位长度的重复序列进行检索
		foreach my $currentUnit ($minUnit..$maxUnit){
			# 若队列未满，则直接入队尾
			if(scalar(@{$Queue{$currentUnit}}) < $currentUnit){
				push @{$Queue{$currentUnit}}, $Seq[$i];
			# 否则，取出队首，与当前碱基比较
			}else{
				my $outBase = shift @{$Queue{$currentUnit}};
				push @{$Queue{$currentUnit}}, $Seq[$i];
				# 若碱基组成相同，则根据上一次的score值来决定是否更新重复序列起始位置
				#   若上一次score==0，则更新重复序列起始位置
				if($outBase eq $Seq[$i]){
					if($formerScore{$currentUnit} == 0){
						$Start{$currentUnit} = $i;
					}
					$formerScore{$currentUnit} = 1;
				# 若碱基组成不同，则当前重复序列的延伸终止
				#   判断当前重复序列的长度是否满足输出要求
				}else{
					if(($i - $Start{$currentUnit} >= ($minRepeat-1)*$currentUnit) and ($i - $Start{$currentUnit} < $maxRepeat*$currentUnit)){
						# 输出时，根据队列中的元素还原出重复序列单元
						my $mode = ($i - $Start{$currentUnit}) % $currentUnit;
						my $start = $Start{$currentUnit} - $currentUnit;
						my $Unit = join("",@{$Queue{$currentUnit}}[($currentUnit-$mode-1)..($currentUnit-2)]).$outBase.join("",@{$Queue{$currentUnit}}[0..($currentUnit-$mode-2)]);
						my $repeat = floor(($i - $Start{$currentUnit})/$currentUnit) + 1;
						
						if($start > $Head){
							print "$seqId\t$start\t$Unit\t$repeat\n";
							$Head = $start;
						}
					}
					$formerScore{$currentUnit} = 0;
					$Start{$currentUnit} = $i;	
				}
			}
		}
	}
}


my $infile = $ARGV[0];
my $minUnit = $ARGV[1];
my $maxUnit = $ARGV[2];
my $minRepeat = $ARGV[3];
my $maxRepeat = $ARGV[4];

open IN, "<$infile" or die "$!\n";

my $seqId;
my %Sequence;
while(<IN>){
	chomp;
	next if(/^\s+$/);
	if(/^>(.*)$/){
		$seqId = $1;
	}else{
		$Sequence{$seqId} .= $_;
	}
}
close IN;

foreach my $currentSeqId (keys %Sequence){
	&TandemRepeat($currentSeqId, \$Sequence{$currentSeqId}, $minUnit, $maxUnit, $minRepeat, $maxRepeat);
}
