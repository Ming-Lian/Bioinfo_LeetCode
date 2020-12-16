#/usr/bin/perl
use strict;
use warnings;

sub BWT_Trans(){
	my ($infile) = @_;
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
	# 循环出队入队，得到BWT matrix
	my @BWT;
	push @BWT, join("", @seq);
	my $n=0;
	print STDERR "Cyclic rotations:\n";
	print STDERR "\t$n:$BWT[$n]\n";
	my $tail = pop @seq;
	while($seq[$#seq] ne "\$"){
		$n++;
		unshift @seq, $tail;
		push @BWT, join("",@seq);
		print STDERR "\t$n:$BWT[$n]\n";
		$tail = pop @seq;
	}
	print STDERR "\n";
	# 对@BWT中的各个字符串元素按照字母顺序逐一进行遍历，截取最后一个字母
	print STDERR "BWT matrix:\n";
	my @LastColumn;
	$n=0;
	foreach my $string (sort {$a cmp $b} @BWT){
		print STDERR "\t$n:$string\n";
		push @LastColumn, substr($string, length($string)-1);
		$n++;
	}
	print STDERR "\n";
	print STDERR "BWT output:\n";
	print STDERR "\t".join("", @LastColumn)."\n\n";
	return @LastColumn;
}

sub BWT_revTrans(){
	my ($BWT) = @_;
	my @Array;
	foreach my $i (0..$#$BWT){
		# 在@Array之前追加新的一列，为@BWT
		foreach my $j (0..$#$BWT){
			$Array[$j] = defined $Array[$j] ? $$BWT[$j].$Array[$j] : $$BWT[$j];
		}
		print STDERR "Roud-$i\n\t";
		print STDERR join("\n\t", @Array)."\n\n\t";
		# 对@Array进行排序
		my @sArray;
		foreach my $string (sort {$a cmp $b} @Array){
			push @sArray, $string;
		}
		@Array = @sArray;
		print STDERR join("\n\t", @Array)."\n";
	}
	return $Array[0];
}

my $infile = $ARGV[0];
# 构建BWT的LastColumn
my @BWT = &BWT_Trans($infile);
my $seq = &BWT_revTrans(\@BWT);

print STDERR "The reconstructed sequence: $seq\n";
