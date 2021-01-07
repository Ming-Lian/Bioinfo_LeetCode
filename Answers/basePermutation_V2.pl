#!/usr/bin/perl

use strict;
use warnings;

sub basePermutation(){
	my ($Base, $n) = @_;
	my @Queue;
	my $i = 0;
	# 重复执行以下操作，直到当前碱基序列长度达到n
	while($i < $n){
		# 当前碱基位置为0时，逐一让@Base中的碱基入队
		if($i==0){
			foreach my $b (@{$Base}){
				push @Queue,$b;
			}
		}else{
			my $l = $#Queue + 1; #获取当前的队列长度
			# 反复执行l次【1次出队-4次入队】的组合操作
			while($l--){
				# 1次出队
				my $currentBase = shift @Queue;
				# 4次入队
				foreach my $b (@{$Base}){
					push @Queue,$currentBase.$b;
				}
			}
		}
		$i++;
	}
	return @Queue;
}

my $n = $ARGV[0];
my @Base = ('A', 'T', 'C', 'G');
my @Queue = &basePermutation(\@Base, $n);

# 输出预览
print join("\n", @Queue)."\n";
