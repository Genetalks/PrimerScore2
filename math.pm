##二分法查找目标值$v在数组@{$aarr}的子数组（index从left到right的子数组）中的index；
##数组aarr按从小到大排序，可以存在多个相同的值
##type为">="表示若遇到相同的值，则取最小的index，而"<="取最大的index。
sub binarySearch{
    my ($v, $aarr, $type, $left, $right)=@_;
    my $mid;
    my $maxi = scalar @{$aarr} -1 ;
    while($left<=$right){
        $mid = int(($left+$right)/2);
        if(($type eq ">=" && $aarr->[$mid]>=$v && ($mid==0 || ($mid>0 && $aarr->[$mid-1]<$v))) || ($type eq "<=" && $aarr->[$mid]<=$v && ($mid==$maxi || ($mid<$maxi && $aarr->[$mid+1]>$v)))){
            return $mid;
        }
        if($aarr->[$mid]>$v || ($aarr->[$mid]==$v && $type eq ">=")){
            $right = $mid-1;
        }else{
            $left = $mid+1;
        }
    }
    return -1;
}



sub max{#&max(lists or arry);
	#求列表中的最大值
	my $max=shift;
	my $temp;
	while (@_) {
		$temp=shift;
		$max=$max>$temp?$max:$temp;
	}
	return $max;
}


sub min{#&min(lists or arry);
	#求列表中的最小值
	my $min=shift;
	my $temp;
	while (@_) {
		$temp=shift;
		$min=$min<$temp?$min:$temp;
	}
	return $min;
}

################################################################################################################
sub N1{
	#求公共集合的大小，N(1,10)=10;N(1,10,5,10)=6;N(1,10,5,10,8,10)=3;参数是“开始位置，结束位置” 的重复。
	my $sum=0;
	my $L=shift @_;
	my $R=shift @_;
	while (@_) {

		my $s=shift @_;
		my $e=shift @_;

		$L=$L>$s?$L:$s;
		$R=$R<$e?$R:$e;
	}

	$sum=$R>=$L?($R-$L+1):0;
	return $sum;
}


sub N2{
	#求公共集合，N(1,10)=(1,10);N(1,10,5,10)=(5,10);N(1,10,5,10,8,10)=(8,10);参数是“开始位置，结束位置” 的重复。
	my $sum=0;
	my $L=shift @_;
	my $R=shift @_;
	while (@_) {

		my $s=shift @_;
		my $e=shift @_;

		$L=$L>$s?$L:$s;
		$R=$R<$e?$R:$e;
	}
	return ($L,$R);
}

sub NOT
#function: translate align_position to not_align_position
#input:($length,@array) 
#output:(@array)
#for example (100,1,5,20,58) -> (6,19,59,100)
{
	my($length,@array) = @_;
	my $i = 0;
	my @output = ();
	my $start = 0;
	my $end = 0;
	
	@array = &cat(1,@array);
	for ($i=0;$i<@array;$i+=2) 
	{
		$array[$i] -= 1;
		$array[$i+1] += 1;
	}
	push(@array,$length);
	unshift(@array,1);
	@output = &cat(1,@array);
	return(@output);
}

sub dog
#function: translate position to the contrary strand
#input:($length,@array) 
#output:(@array)
#for example (100,1,5,20,58) -> (43,81,96,100)
{
	my($length,@array) = @_;
	my $i = 0;
	my @output = ();
	my $start = 0;
	my $end = 0;

	for ($i=0;$i<@array;$i+=2) 
	{
		$start = $length+1-$array[$i+1];
		$end = $length+1-$array[$i];
		unshift(@output,$end);
		unshift(@output,$start);
	}

	return(@output);
}

sub L
#function: compute total length by start and stop
#input:(@array)
#output:($length)
#for example (1,3,5,8)->(7)
{
	my(@input) = @_;
	my $i = 0;
	my $length = 0;

	for ($i=0;$i<@input;$i+=2)
	{
		$length += $input[$i+1] - $input[$i] + 1;
	}

	return($length);
}

################################################################################################################
sub aveStd{
	#calculate and return the average and variance of the arry.
	#计算数组的平均值和方差
	my @arr=@_;
	my $x=0;
	my $x2=0;
	my $count=scalar @arr;
	for (my $i=0;$i<$count ;$i++) {
		$x+=$arr[$i];
		$x2+=$arr[$i]**2;
	}
	my $aver=$x/$count;
	my $variance=($count*$x2-$x**2)/($count*($count-1));
	return $aver,$variance;
}

################################################################################################################
sub over
#function: find overlap between two groups
#input:($array1,$array2)
#output:(@array)
#for example (0,1,5,8) + (1,3,4,8) -> (1,1,5,8)
{
	my($array1,$array2) = @_;
	my @array1 = split(/\s+/,$array1);
	my @array2 = split(/\s+/,$array2);
	my $i = 0;
	my $j = 0;
	my $s = 0;
	my @output = ();
	my $start = 0;
	my $end = 0;

	for ($i=0;$i<@array1;$i+=2)
	{
		for ($j=$s;$j<@array2;$j+=2)
		{
			if($array1[$i+1] < $array2[$j]) { last; }
			elsif($array1[$i] > $array2[$j+1]) { $s = $j + 2; next; }
			else
			{
				if($array1[$i] < $array2[$j]) { $start = $array2[$j]; } else { $start = $array1[$i]; }
				if($array1[$i+1] < $array2[$j+1]) { $end = $array1[$i+1]; } else { $end = $array2[$j+1]; }
				if($end < $start) { die"error: END $end < START $start\n"; }
				push(@output,$start);
				push(@output,$end);
				if($array1[$i+1] > $array2[$j+1]) { $s = $j + 2; }
			}
		}
	}

	return(@output);
}

sub over_len
#function: find overlap between two groups
#input:($array1,$array2)
#output:(@array)
#for example (0,1,5,8) + (1,3,4,8) -> 5
{
	my($array1,$array2) = @_;
	my @array1 = split(/\s+/,$array1);
	my @array2 = split(/\s+/,$array2);
	my $i = 0;
	my $j = 0;
	my $s = 0;
	my @output = ();
	my $start = 0;
	my $end = 0;

	for ($i=0;$i<@array1;$i+=2)
	{
		for ($j=$s;$j<@array2;$j+=2)
		{
			if($array1[$i+1] < $array2[$j]) { last; }
			elsif($array1[$i] > $array2[$j+1]) { $s = $j + 2; next; }
			else
			{
				if($array1[$i] < $array2[$j]) { $start = $array2[$j]; } else { $start = $array1[$i]; }
				if($array1[$i+1] < $array2[$j+1]) { $end = $array1[$i+1]; } else { $end = $array2[$j+1]; }
				if($end < $start) { die"error: END $end < START $start\n"; }
				push(@output,$start);
				push(@output,$end);
				if($array1[$i+1] > $array2[$j+1]) { $s = $j + 2; }
			}
		}
	}

	my $len = 0;
	for (my $i=0; $i<@output; $i+=2) {
		$len += $output[$i+1]-$output[$i]+1;
	}

	return $len;
}


################################################################################################################
sub cat
#function:quit redundance
#input:($para,@array), para is the merge length
#output:(@array),
#for example (0,1,3,4,7,5,8)->(1,3,4,8) (1,1,3,4,7,5,8)->(1,8)
{
	my($merge,@input) = @_;
	my $i = 0;
	my @output = ();
	my %hash = ();
	my $each = 0;
	my $begin = "";
	my $end = 0;


	for ($i=0;$i<@input;$i+=2)
	{
		my $Qb = $input[$i];
		my $Qe = $input[$i+1];

		if($Qb > $Qe) { next; }
		if(defined($hash{$Qb}))	{ if($hash{$Qb} < $Qe) { $hash{$Qb} = $Qe; } }
		else { $hash{$Qb} = $Qe; }
		$Qb = 0;
	}

	foreach $each (sort {$a <=> $b} keys %hash)
	{
		if($begin eq "")
		{
			$begin = $each;
			$end = $hash{$each};
		}
		else
		{
			if($hash{$each} > $end)
			{
				if($each > $end + $merge)
				{
					push(@output,$begin);
					push(@output,$end);
					$begin = $each;
					$end = $hash{$each};
				}
				else { $end = $hash{$each}; }
			}
		}
	}
	if(keys %hash > 0)
	{
		push(@output,$begin);
		push(@output,$end);
	}

	%hash = ();

	return(@output);
}


sub overlen {
	my ($p1,$p2,$t1,$t2)=@_;
		($p1,$p2)=sort {$a <=> $b} ($p1,$p2);
		($t1,$t2)=sort {$a <=> $b} ($t1,$t2);

		my $m = $p1<$t1 ? $t1 : $p1;
		my $n = $p2<$t2 ? $p2 : $t2;
		if ($m > $n) {
			return 0;
		}
		else {
			return $n-$m+1;
		}
}

1;

