#/usr/bin/perl
use strict;
use warnings;

sub FM_index(){
	my ($infile, $SA, $OCC, $LastColumn, $Pre) = @_;
	open IN, "<$infile" or die "$!\n";
	# 解析fasta文件
	my $fa_seq;
	while(<IN>){
		chomp;
		next if(/^>/);
		$fa_seq .= $_;
	}
	close IN;
	# 将以字符串形式组织的序列，打散成数组，并在末尾追加终止标识符
	my @seq = split //,$fa_seq;
	push @seq, "\$";
	print STDERR "Input:\n";
	print STDERR "\t".join("",@seq)."\n\n";
	# 循环出队入队，得到BWT matrix，并在每行末尾加上末尾字符在原始序列的位置索引（1-base）
	my @BWT;
	my @index;
	my $N = $#seq + 1;
	push @BWT, join("", @seq);
	push @index, $N;
	my $n=0;
	print STDERR "Cyclic rotations:\n";
	print STDERR "\t$n:$BWT[$n]\t$index[$n]\n";
	my $tail = pop @seq;
	while($seq[$#seq] ne "\$"){
		$n++;
		$N--;
		unshift @seq, $tail;
		push @BWT, join("",@seq);
		push @index, $N;
		print STDERR "\t$n:$BWT[$n]\t$index[$n]\n";
		$tail = pop @seq;
	}
	print STDERR "\n";
	# 对@BWT中的各个字符串元素按照字母顺序逐一进行遍历，得到FM-index的3个数组和1个散列
	print STDERR "Sort and construct FM-index:\n";
	print STDERR "\tBWTmatrix\tSA\tOCC\n";
	$n=0;
	my %Count;
	## 获得3个数组：SA, OCC, BWT
	foreach my $i (sort {$BWT[$a] cmp $BWT[$b]} 0..$#BWT){
		my $LC = substr($BWT[$i], length($BWT[$i])-1);
		push @$LastColumn, $LC;
		push @$SA, $index[$i];
		$Count{$LC} = defined $Count{$LC} ? $Count{$LC} + 1 : 0;
		push @$OCC, $Count{$LC};
		print STDERR "\t$n:$BWT[$i]\t$$SA[$n]\t$$OCC[$n]\n";
		$n++;
	}
	print STDERR "\n";
	print STDERR "FM-index:\n";
	## 获得1个散列：Pre
	print STDERR "\tPre:";
	$n = 0;
	foreach my $base (sort {$a cmp $b} keys %Count){
		$$Pre{$base} = $n;
		$n += $Count{$base}+1;
		print STDERR "$base=$$Pre{$base};";
	}
	print STDERR "\n";
	print STDERR "\tSA :".join("", @$SA)."\n";
	print STDERR "\tOCC:".join("", @$OCC)."\n";
	print STDERR "\tBWT:".join("", @$LastColumn)."\n\n";
}

sub BWT_reverse_transformation(){
	my ($SA, $OCC, $BWT, $Pre) = @_;
	my @seq;
	my $currentBase = "\$";
	my $currentIndex = 0;
	my $index;
	do{
		unshift @seq, $currentBase;
		# Back_trace，得到当前碱基的前一位碱基的组成
		$currentBase = &Back_trace($currentIndex,$BWT);
		# LC2FC，得到LC中的当前碱基在FC对应的碱基位置
		$currentIndex = &LC2FC($currentBase, $currentIndex, $Pre, $OCC);
	}while($currentBase ne "\$");
	print STDERR "Reverse transformation:".join("", @seq)."\n\n";
}

sub Back_trace(){
	my ($index, $BWT) = @_;
	return $$BWT[$index];
}

sub LC2FC(){
	my ($base, $index, $Pre, $OCC) = @_;
	$index = $$Pre{$base} + $$OCC[$index];
	return $index;
}

my $infile = $ARGV[0];

# 构建FM-index
my (@SA, @OCC, @BWT, %Pre);
&FM_index($infile, \@SA, \@OCC, \@BWT, \%Pre);
print "FM-index:\n";
print "\tPre:";
foreach my $base (sort {$a cmp $b} keys %Pre){
	print STDERR "$base=$Pre{$base};";
}
print "\n";
print "\tSA :".join("", @SA)."\n";
print "\tOCC:".join("", @OCC)."\n";
print "\tBWT:".join("", @BWT)."\n\n";

# 从FM-index还原原始参考序列
&BWT_reverse_transformation(\@SA, \@OCC, \@BWT, \%Pre);
