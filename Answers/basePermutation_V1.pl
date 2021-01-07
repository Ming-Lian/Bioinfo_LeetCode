use POSIX;

@base = ('A', 'T', 'C', 'G');
foreach $Col (0..7){
	$maxRow = 4 ** 8;
	$subRow = $maxRow / (4**($Col+1));
	foreach $Row (0..$maxRow){
		$index = ceil(($Row+1) / $subRow) % 4; # 用1-base来算
		$index = ($index + 3) % 4; # 换算成0-base
		$str[$Row] .= $base[$index];
	}
}

foreach $substr (@str){
	print "$substr\n";
}
