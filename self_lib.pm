## 打分
## score:满分
sub score_single{
	my ($v, $score, $minb, $maxb, $min, $max)=@_;
	return $score if($v eq "NULL");
	my $disl = $minb - $min;
	my $disr = $max - $maxb;
	if($disl<0 || $disr<0){
		print "Wrong Score set: $minb, $maxb, $min, $max\n";
		die;
	}

	my $s;
	if($v<$minb){
		$disl=$disl==0? 0.1: $disl;
		$s = int($score * (1 - ($minb-$v)/$disl));
	}elsif($v<=$maxb){
		$s = $score;
	}else{
		$disr= $disr==0? 0.1: $disr;
		$s = int($score * (1 - ($v-$maxb)/$disr));
	}
	return $s;
}

## 平均提取函数
## dis是距离，rfloat是浮动的比例（0.2）,apos是位置hash，ascore是打分hash，select:"UP",表示选取分值最高的，否则选取分值最低的。
sub average{
	my ($dis, $rfloat, $apos, $ascore,$select)=@_;
    my $ss = 0;
    my $es = 0;
	my $last_select=0;
	my @id_sortp = sort{$apos->{$a}<=>$apos->{$b}} keys %{$apos};	
    my @pos_sort = sort{$a<=>$b} values %{$apos};
	my $left = 0;
	my $right = $#pos_sort;
	my @final;
	while($es<=$pos_sort[-1]){
        my ($s, $e)=($ss-$rfloat*$dis, $es+$rfloat*$dis);
		my ($is, $ie)=&get_ids_index($s, $e, \@pos_sort, $left, $right);
        if($is==-1 || $ie==-1){
            if($s<=$last_select+$rfloat*$dis){
                ($ss, $es)=($last_select+$rfloat*$dis+1, $e);
            }else{
                ($ss, $es)=($s, $e);
            }
			($left, $right)=(0, $#pos_sort);
            next;
        }
		
        my @id = @id_sortp[$is..$ie];
		my ($idb, $endd)=&get_best_and_middle_id(\@id, $ascore, $apos, $s, $e, $select);
		if($e-$apos->{$idb} < $dis*$rfloat*0.2){ ## the best id is near to the end
            if($s<=$last_select+$rfloat*$dis){
                ($ss, $es)=($last_select+$rfloat*$dis+1, $e-$rfloat*$dis*0.8);
            }else{
                ($ss, $es)=($s, $e-$rfloat*$dis*0.8);
            }
			($left, $right)=(0, $#pos_sort);
            next;
		}
		push @final, $idb;
        $ss = $apos->{$idb}+$dis;
        $es = $ss;
        $last_select = $apos->{$idb};
		($left, $right)=($is+1, $#pos_sort);
    }

	return @final;
}

## get id of best score and max middle, and its distance to end
sub get_best_and_middle_id{
	my ($aid, $ascore, $apos, $poss, $pose, $select)=@_;
	my @id = @{$aid};
	my @id_sort;
	if($select eq "UP"){
    	@id_sort = sort{$ascore->{$b} <=> $ascore->{$a}} @id;
	}else{
    	@id_sort = sort{$ascore->{$a} <=> $ascore->{$b}} @id;
	}
	my $best_score = $ascore->{$id_sort[0]};
	my %best_dis;
	for(my $i=0; $i<@id_sort; $i++){
		if($ascore->{$id_sort[$i]} == $best_score){
			my $pos = $apos->{$id_sort[$i]};
			$best_dis{$id_sort[$i]}=&min($pos-$poss, $pose-$pos); #distance between the best index and the end
		}else{
			last;
		}
	}
	## get the most middle best id
	my @idb_sort = sort{$best_dis{$b} <=> $best_dis{$a}} keys %best_dis;
	my $idb = $idb_sort[0];
	return ($idb, $best_dis{$idb});
}
sub get_ids_index{
	my ($s, $e, $aarr, $left, $right)=@_;
	my $istart = &binarySearch($s, $aarr, ">=", $left, $right);
	return (-1, -1) if($istart==-1 || $aarr->[$istart]>$e);

	my $iend = &binarySearch($e, $aarr, "<=", $istart+1, $right);
	return ($istart, $iend);
}

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


sub explain_bam_flag{
	my ($flag)=@_;
	my $flag_bin=sprintf("%b", $flag);
	my @flag_bin = split //, $flag_bin;
	my $is_read1 = $flag_bin[-7];
	my $is_read2 = @flag_bin>=8? $flag_bin[-8]: 0;
	my $is_supplementary = @flag_bin>=12? $flag_bin[-12]: 0;
	my $is_proper_pair = $flag_bin[-2];
	my $is_reverse = $flag_bin[-5];
	my $is_unmap = $flag_bin[-3];
	my $is_munmap = $flag_bin[-4];
	my $dup = @flag_bin>=11? $flag_bin[-11]: 0;
	my $rnum = $is_read1==1? 1: 2;
	return($rnum, $is_proper_pair, $is_reverse, $is_unmap, $is_munmap, $is_supplementary);
}

sub SHOW_TIME {
	#显示当时时间函数，参数内容是时间前得提示信息，为字符串
	my ($str)=@_;
	my ($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst)=localtime(time());
	my $temp=sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
	print "$str:\t[".$temp."]\n";
}

################################################################################################################

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

################################################################################################################

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

sub revcom{#&revcom($ref_seq);
	#获取字符串序列的反向互补序列，以字符串形式返回。ATTCCC->GGGAAT
	my $seq=shift;
	$seq=~tr/ATCGatcg/TAGCtagc/;
	$seq=reverse $seq;
	return uc $seq;
}

################################################################################################################

sub subseq{
	my $seqId_seq=shift;
	my $seqId=shift;
	my $beg=shift;
	my $end=shift;
	my $strand=shift;

	my $subseq=substr($seqId_seq->{$seqId},$beg-1,$end-$beg+1);
	if ($strand eq "-") {
		$subseq=revcom($subseq);
	}
	return uc $subseq;
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

################################################################################################################

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

################################################################################################################
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
sub SortByNum {	###### 对那种以数字结尾的list按数字的大小排序, 排序类型用"+"表示升序，"-"表示降序。
		###### 例如：LG1，LG6，LG3，LG5 按升序排序后得到：LG1，LG3，LG5，LG6
	my ($sort_type,$ref_list,$ref_sorted)=@_;
	my %hash;
	my $sum=0;
	foreach my $name (@{$ref_list}) {
		$sum++;
		my ($word,$num)=$name=~/(\S*\D+)(\d+)/;
		$hash{$word}{$sum}=$num;
	}
	if ($sort_type eq "+") {
		foreach my $word (sort {$a cmp $b} keys %hash) {
			foreach my $sum (sort {$hash{$word}{$a} <=> $hash{$word}{$b}} keys %{$hash{$word}}) {
				push @{$ref_sorted},$word.$hash{$word}{$sum};
			}
		}
	}elsif($sort_type eq "-"){
		foreach my $word (sort {$a cmp $b} keys %hash) {
			foreach my $sum (sort {$hash{$word}{$b} <=> $hash{$word}{$a}} keys %{$hash{$word}}) {
				push @{$ref_sorted},$word.$hash{$word}{$sum};
			}
		}
	}
}
################################################################################################################
sub fileFormat()
{	#查询文件格式，目前支持fasta，fastq，GenBank,unkonw
	my $file =shift;
	open (IN,"<",$file) or die $!;
	my $type=<IN>;
	close (IN) ;
	chomp $type;
	if ($type=~/^>/) {
		$type="fasta";
	}elsif($type=~/^@/){
		$type="fastq";
	}elsif($type=~/LOCUS/){
		$type="GenBank";
	}else{
		$type="nuknow";
	}
	return $type;
}

################################################################################################################
sub indexFastQA {#
	#获取fasta,fastq的Index
	my $file=shift;
	$file=~s/fa$//;
	$file=~s/fq$//;
	$file=~s/fasta$//;
	$file=~s/fastq$//;
	$file=~s/.$//;
	$file=~s/-$//;
	$file=~s/_$//;
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

sub trim
{
		my $s=shift;
		$s=~/^(.*):.*$/;
		return $1;
}

sub sub_format_datetime #datetime subprogram
{
		my($sec , $min , $hour , $day , $mon , $year , $wday , $yday , $isdst) = @_;

		sprintf("%4d-%02d-%02d %02d:%02d:%02d" , ($year + 1900) , $mon , $day , $hour , $min , $sec);
}


sub TimeCount() {
		my ($start,$stop)=@_;
		my $min=int(($stop-$start)/60);
		my $second=sprintf("%02d",($stop-$start) % 60);
		my $return=$min.":".$second;
		return ($return);
}
sub define_axis () {
		 my $i=0;
		 my ($max)=@_;
		 my $time=1;
		 my @ret=();
		 my @limit=(1,2,3,4,5,6,8,10,12,15,16,20,24,25,30,40,50,60,80,100,120);
		 my @unlim=(0,1,2,3,4,5,6,8 ,10,11,14,15,18.5,21,23,29,37,47,56,76 ,92 ,110);

		 while ($max >$unlim[21]) {
				   $max=$max/10;
				   $time=$time*10;
		 }
		 for ($i=0;$i<=20 ;$i++) {
				   if ($max>$unlim[$i] && $max<=$unlim[$i+1]) {
							$ret[0]=$limit[$i]*$time;
							
							if ($i==2 || $i==5 || $i==9 || $i==14) {
									 $ret[1]=$ret[0]/3;
							}
							elsif ($i==4 || $i==7 || $i==13 || $i==16){
									 $ret[1]=$ret[0]/5;
							}
							else {
									 $ret[1]=$ret[0]/4;
							}

				   }
		 }
		 return @ret;
}


################################################################################################################
sub SSR {#返回序列的低复杂度
	my @arr=split "",uc shift;
	my %hash;
	$hash{$_}++ foreach (@arr);
	my ($first,$second)=sort {$b <=> $a} values %hash;
	return ($first+$second)/scalar @arr;
}

1;
