

sub md_split{
   my($md)=@_;
   my @unit = split //, $md;
   my @mds;
   my $temp;
   for(my $i=0; $i<@unit; $i++){
       if($unit[$i]=~/[ATCGN]/){
		   push @mds, $temp;
           push @mds, $unit[$i];
           $temp = "";
       }elsif($unit[$i] eq "^"){
		   push @mds, $temp;
           $temp = $unit[$i];
		   my $j;
           for($j=$i+1;$j<@unit;$j++){
               last if($unit[$j]!~/[ATCGN]/);
               $temp .= $unit[$j];
           }
		   $i=$j-1;
		   push @mds, $temp;
		   $temp="";
       }elsif($unit[$i]=~/[0-9]/){
		   $temp.=$unit[$i];
	   }
   }
   push @mds, $temp;
   return(\@mds);
}

sub cigar_split{
    my($cigar, $keepH)=@_;
    my @ucigar = split //, $cigar;
    my (@match_n,@match_str);
    my $nstr="";
    foreach my $i(0..$#ucigar){
        if($ucigar[$i] eq "M" || $ucigar[$i] eq "I" || $ucigar[$i] eq "D" || $ucigar[$i] eq "S"){
            push @match_str, $ucigar[$i];
            push @match_n, $nstr;
            $nstr = "";
        }elsif($ucigar[$i] eq "H"){
			if(defined $keepH){
				push @match_str, $ucigar[$i];
				push @match_n, $nstr;
			}
            $nstr = "";
        }else{
            $nstr .= $ucigar[$i];
        }
    }
    return(\@match_n,\@match_str);
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


sub GC_stat{
        my ($p)=@_;
        my @u = split //, $p;
        my %stat;
        $stat{"total"}{'G'}=0;
        $stat{"total"}{'C'}=0;
        $stat{5}{'G'}=0;
        $stat{5}{'C'}=0;
        $stat{3}{'G'}=0;
        $stat{3}{'C'}=0;
        my $total = scalar @u;
        for(my $i=0; $i<$total; $i++){
            if ($i<$GC5_num){
                $stat{5}{$u[$i]}++;
            }
            if ($i>$total-$GC3_num-1){
                $stat{3}{$u[$i]}++;
            }
            $stat{"total"}{$u[$i]}++;
        }
        my $GC = ($stat{"total"}{'G'}+$stat{"total"}{'C'})/$total;
        my $GC5= ($stat{5}{'G'}+$stat{5}{'C'})/$GC5_num;
        my $GC3 = ($stat{3}{'G'}+$stat{3}{'C'})/$GC3_num;
        return ($GC, $GC5, $GC3);
}

sub GC{
	my ($s)=@_;
   	my @u=split //, $s;
	my $total = 0;
	my $gc = 0;
	foreach $b (@u){
		$total++;
		if($b eq 'G' || $b eq 'C' || $b eq "g" || $b eq "c"){
			$gc++;
		}
	}
	return($gc/$total);
}

sub GC_info_stat{
    my ($p, $Wind_GC)=@_;
    my @u = split //, $p;

	my $GC = &GC($p);
	my ($GC5, $GC3);
	my $min_GC=1;
	my $max_GC=0;
	for(my $i=0; $i<@u-$Wind_GC+1; $i++){
		my @sub = @u[$i..$i+$Wind_GC-1];
		my $sub_gc = &GC(@sub);
		if($i==0){
			$GC5 = $sub_gc;
		}
		if($i==@u-$Wind_GC){
			$GC3 = $sub_gc;
		}
		if($sub_gc>$max_GC){
			$max_GC = $sub_gc;
		}
		if($sub_gc<$min_GC){
			$min_GC = $sub_gc;
		}
	}
	$GC = sprintf "%0.2f", $GC;
    return ($GC, $GC5, $GC3, $max_GC-$min_GC);
}


sub revcom{#&revcom($ref_seq);
	#获取字符串序列的反向互补序列，以字符串形式返回。ATTCCC->GGGAAT
	my $seq=shift;
	$seq=~tr/ATCGatcg/TAGCtagc/;
	$seq=reverse $seq;
	return $seq;
}

###############################################################################################################
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
sub SSR {#返回序列的低复杂度
	my @arr=split "",uc shift;
	my %hash;
	$hash{$_}++ foreach (@arr);
	my ($first,$second)=sort {$b <=> $a} values %hash;
	return ($first+$second)/scalar @arr;
}

sub AbsolutePath
{		#获取指定目录或文件的决定路径
		my ($type,$input) = @_;

		my $return;
	$/="\n";

		if ($type eq 'dir')
		{
				my $pwd = `pwd`;
				chomp $pwd;
				chdir($input);
				$return = `pwd`;
				chomp $return;
				chdir($pwd);
		}
		elsif($type eq 'file')
		{
				my $pwd = `pwd`;
				chomp $pwd;

				my $dir=dirname($input);
				my $file=basename($input);
				chdir($dir);
				$return = `pwd`;
				chomp $return;
				$return .="\/".$file;
				chdir($pwd);
		}
		return $return;
}

sub GetTime {
	my ($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst)=localtime(time());
	return sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}

sub SHOW_TIME {
	#显示当时时间函数，参数内容是时间前得提示信息，为字符串
	my ($str)=@_;
	my ($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst)=localtime(time());
	my $temp=sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
	print "$str:\t[".$temp."]\n";
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


1;


