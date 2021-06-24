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

1;

