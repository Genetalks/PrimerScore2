# $aarr: target info array with order
# $avalue, $ascore, $type: eg: 
#       (1, 0.9, 0), (1, 0, -1), "Number": when v>=1, s=1; 1>v>=0.9, s=0; 0.9>v>=0, s=-1;
#       ("|", "*", "-", "^"), (2,1,-2,-2), "String": when v="|", s=2; when v="-", s=-2
# $max_extend: max extend num, ie. max sequential num whose score<0
# return: ((is1,ie1), (is2, ie2)...)
sub segment{
	my ($aarr, $avalue, $ascore, $max_extend, $type, $is_back)=@_;
	my @arr=@{$aarr};
	my @value=@{$avalue};
	my @score=@{$ascore};
	my $stotal=0;
	my ($istart, $iend)=(-1, -1);
	my @result;
	for(my $i=0; $i<@arr; $i++){
		for(my $j=0; $j<@value; $j++){
			if(($type eq "Number" && $arr[$i]>=$value[$j]) || ($type eq "String" && $arr[$i] eq $value[$j])){
				if($istart==-1){## not start
					if($score[$j]>0){
						$istart=$i;
						$stotal+=$score[$j];
					}
				}else{
					$stotal+=$score[$j];
					if($score[$j]>0){
						$extend=0;
					}else{
						$extend++;
					}

					if($stotal<=0 || $extend>=$max_extend){## Interrupt
						$iend=$i;
					}
				}
				last;
			}
		}
		if($iend!=-1){
			push @result, [$istart, $iend];
			($istart, $iend)=(-1, -1);
			$stotal=0;
		}
	}
	if($istart!=-1 && $iend==-1){
		$iend=scalar @arr -1;
		push @result, [$istart, $iend];
	}
	if($is_back){
		for(my $x=0; $x<@result; $x++){
			my ($is, $ie)=@{$result[$x]};
			my $is_back=0;
			for(my $i=$ie; $i>=$is; $i--){
				for(my $j=0; $j<@value; $j++){
					if(($type eq "Number" && $arr[$i]>=$value[$j]) || ($type eq "String" && $arr[$i] eq $value[$j])){
						if($score[$j]>0){
							$result[$x][1]=$i;
							$is_back=1;
						}
						last;
					}
				}
				last if($is_back);
			}
		}
	}
	return (@result);
}

sub window_ratio{
	my ($ainfo, $wd, $type)=@_;
	my @info=@{$ainfo};
	my $len = scalar @info;
	my @ratio;
	for(my $i=0; $i<$len-$wd+1; $i++){
		my $r;
		if($type eq "CLW_mismatch"){
			$r = &mismatch_ratio(@info[$i..($i+$wd-1)]);
		}elsif($type eq "GC"){
			$r = &GC_stat_array(@info[$i..($i+$wd-1)]);
		}
		push @ratio, $r;
	}
	return @ratio;
}

sub mismatch_ratio{
	my @unit=@_;
	my ($mat, $mis)=(0,0);
	foreach my $u(@unit){
		if($u eq "*"){
			$mat++;
		}else{
			$mis++;
		}
	}
	return ($mat/($mat+$mis));
}

1;
