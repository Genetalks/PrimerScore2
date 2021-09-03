## return: -1:timeout; 1:finished
sub Run_monitor_timeout{
	my ($time, $cmd, $sh)=@_;
	if(defined $sh){
		print $sh $cmd,"\n";
		print "###", $cmd,"\n";
	}

	eval {
        local $SIG{ALRM} = sub { die "alarm\n" };
        alarm $time;
		$? = `$cmd`;
        alarm 0;
    };
    if ($@) { # timed out
        die unless $@ eq "alarm\n";
		print "Run time out! $cmd\n";
		return -1;
    }else { # didn't timeout
		return 1;
    }
}

sub Run{
    my ($cmd, $sh, $nodie)=@_;
	if(defined $sh){
		print $sh $cmd,"\n";
		print "###", $cmd,"\n";
	}
    my $ret = system($cmd);
    if ($ret) {
        print STDERR "Run $cmd failed!\n";
		if(!defined $nodie){
			die;
		}
    }
}


sub convert_StartEndRefAlt{ ## copy from annovar
	my ($start, $ref, $alt)=@_;
	my ($head, $newstart, $newend, $newref, $newalt);

	if (length ($ref) == 1 and length ($alt) == 1) {   #SNV
	    ($newstart, $newend) = ($start, $start+length($ref)-1);
	    ($newref, $newalt) = ($ref, $alt);
	
	} elsif (length ($ref) > length ($alt)) {       #deletion or block substitution
	    $head = substr ($ref, 0, length ($alt));
	    if ($head eq $alt) {
	        ($newstart, $newend) = ($start+length ($head), $start + length ($ref)-1);
	        ($newref, $newalt) = (substr($ref, length($alt)), '-');
	    } else {
	        ($newstart, $newend) = ($start, $start+length($ref)-1);     #changed to length(ref) on 20130820
	        ($newref, $newalt) = ($ref, $alt);
	    }
	
	    ($newstart, $newend, $newref, $newalt) = adjustStartEndRefAlt ($newstart, $newend, $newref, $newalt);  
	
	} elsif (length ($ref) < length ($alt)) {       #insertion or block substitution
	    $head = substr ($alt, 0, length ($ref));
	    if ($head eq $ref) {
	        ($newstart, $newend) = ($start+length($ref)-1, $start+length($ref)-1);
	        ($newref, $newalt) = ('-', substr($alt, length($ref)));
	    } else {
	        ($newstart, $newend) = ($start, $start+length($ref)-1);
	        ($newref, $newalt) = ($ref, $alt);
	    }
	
	    ($newstart, $newend, $newref, $newalt) = adjustStartEndRefAlt ($newstart, $newend, $newref, $newalt);
	}
	return ($newstart, $newend, $newref, $newalt);	
	
}
sub adjustStartEndRefAlt {
    my ($newstart, $newend, $newref, $newalt) = @_;

    until (substr($newref,-1) ne substr($newalt,-1)) {
        chop $newref;
        chop $newalt;
        $newend--;
        if (not $newref) {
            $newref = '-';
            last;
        }
        if (not $newalt) {
            $newalt = '-';
            last;
        }
    }
    until (substr($newref,0,1) ne substr ($newalt,0,1)) {
        substr ($newref,0,1) = '';
        substr ($newalt,0,1) = '';
        $newstart++;
        if (not $newref) {
            $newref = '-';
            last;
        }
        if (not $newalt) {
            $newalt = '-';
            last;
        }
    }
    return ($newstart, $newend, $newref, $newalt);
}



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

sub GC_array{
	my @u=@_;
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
	my $high_GC=0.8;
	my $low_GC=0.2;
	my ($high_len, $low_len)=(0,0);
	my $GC = &GC($p);
	my $min_GC=1;
	my $max_GC=0;
	for(my $i=0; $i<@u-$Wind_GC+1; $i++){
		my @sub = @u[$i..$i+$Wind_GC-1];
		my $sub_gc = &GC_array(@sub);
		if($sub_gc>=$high_GC){
			$high_len++;
		}
		if($sub_gc<=$low_GC){
			$low_len++;
		}
		if($sub_gc>$max_GC){
			$max_GC = $sub_gc;
		}
		if($sub_gc<$min_GC){
			$min_GC = $sub_gc;
		}
	}
	$GC = sprintf "%0.2f", $GC;
    return ($GC, $high_len, $low_len, $max_GC, $min_GC);
}

sub revcom{#&revcom($ref_seq);
	#获取字符串序列的反向互补序列，以字符串形式返回。ATTCCC->GGGAAT
	my $seq=shift;
	$seq=~tr/ATCGatcg/TAGCtagc/;
	$seq=~tr/RYMKSWHDBVN/YRKMSWDHVBN/;
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


