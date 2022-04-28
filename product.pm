require "$Bin/score.pm";
require "$Bin/math.pm";


#my $prob=$dis."/".join(",", $chr, $sd.$pos, $mvisual1,sprintf("%.2f",$tm1), $sd2.$pos2,$mvisual2,sprintf("%.2f",$tm2));
sub probe_bounds_on_products{
	my ($pb, $abound, $aprod, $aresult, $PCRsize, $opt_tm_probe, $min_tm_spec)=@_;
	my %bdid=%{$abound};
	my $n=0;
	foreach my $prod(keys %{$aprod}){
		my ($dis, $info)=split /\//, $prod;
		my ($chr, $p, $sdpos, $mv, $tm, $p2, $sdpos2, $mv2, $tm2)=split /,/, $info;
		my ($sd1, $pos1)=$sdpos=~/([+-])(\d+)/;
		my ($sd2, $pos2)=$sdpos2=~/([+-])(\d+)/;
		next if(!exists $bdid{$chr});
		my ($minp, $maxp)=sort($pos1, $pos2);
		foreach my $sd(keys %{$bdid{$chr}}){
			my @bds=@{$bdid{$chr}{$sd}};
			for(my $i=0; $i<@bds; $i++){
				my ($pos3, $pos5, $tm, $end3_base, $mvisual)=@{$bds[$i]};
				if($pos5>$minp-10 && $pos5<$maxp+10){
					my $dis=&min($pos5-$minp, $maxp-$pos5);
					my $eff=&efficiency_probe($tm, $dis, $PCRsize, $opt_tm_probe, $min_tm_spec) * $aprod->{$prod};
					my $bpd=$dis."/".join(",", $chr, $sd.$pos5,$mvisual,sprintf("%.2f",$tm)).":".$info;
					$aresult->{$pb}{$bpd}=$eff;
					$n++;
				}
			}
		}
	}
	
	return $n;
}

#push @{$bound{$id}{$chr}{$strand}}, [$chr, $strand, $pos3, $pos5, $tm, $end3_base, $mvisual, $seq, p1_p2, index, ef1, ef_tm, ef_end];
sub caculate_product{
	my ($tid1, $p1, $tid2, $p2, $abound1, $abound2, $ptype, $aprod, $arecord, $PCRsize, $opt_tm, $min_tm, $mind, $maxd, $min_eff, $Max_prodn)=@_;
	my $bdp1=$abound1; # an array of alignment pos
	my $bdp2=$abound2;
	my $prodn=0;

    for (@$bdp1){
        $_->[8] = 0;
    }    
    for (@$bdp2){
        $_->[8] = 1;
    }

	my @all_align_info = (@$bdp1, @$bdp2); #all alignment info of p1 and p2

	##                Nested		 face-to-face         back-to-back
##P1			 	  -->| 		   |-->                    |-->
##P2min		  -->|...      		    <--|             <--|...  
##P2max			      -->| 		       ...<--|                 ...<--|   
	my ($ixpos, $dmin, $dmax);
	if($ptype eq "Nested"){
		$ixpos=2; #pos3
		#Nested: P1 --->|   x        distance range: 5,30 (mind>0)
		#          P2 --->| x

		#Nested: P2 --->|   x        distance range: -30,-5 (mind<0)
		#          P1 --->| x
		if($mind>0){
			$dmin=0;
			$dmax=$PCRsize;
		}else{
			$dmin=-1*$PCRsize;
			$dmax=0;
		}
	}elsif($ptype eq "back-to-back"){
			$ixpos=3; #pos5
			$dmin=-1*$PCRsize+2*20; #20: primer min len
			$dmax=$PCRsize;
	}else{ ##face-to-face
			$ixpos=3; #pos5
			$dmin=30;
			$dmax=$PCRsize;
	}

    my $is_nested_mode = $ptype eq "Nested";

	#sort according to chr and pos
	@all_align_info = sort { $a->[0] cmp $b->[0] or $a->[$ixpos] <=> $b->[$ixpos]} @all_align_info;

	# iter all alignments and found effective products
	for (my $i = 0; $i < @all_align_info; ++$i){  # O(n+m) 
		my $a = $all_align_info[$i];
        my $sd = $a->[1];
        my $pos = $a->[$ixpos];
        my $sd2 = $sd;
        my $is_rev = $sd eq "-";
        if (!$is_nested_mode){
            $sd2 = ("+", "-")[1-$is_rev];
        }
        my $pmin = $is_rev ? $pos - $dmax : $pos + $dmin;;
        my $pmax = $is_rev ? $pos - $dmin : $pos + $dmax;
        #my ($sd2, $pmin, $pmax)=&primer2_scope(\$ptype, $a->[1], $a->[$ixpos], $dmin, $dmax);
        for (my $j = $i + 1; ; $j++){ # look forward
			last if ($j >= @all_align_info); # no more alignments
			my $b = $all_align_info[$j];
			last if ($b->[$ixpos] > $pmax);  # exceed valid dist range
			last if ($b->[0] ne $a->[0]);  # diff chr
			next if ($b->[1] ne $sd2); # not valid strand
			# $a and $b create effective product
			my $dis;
			if($a->[1] eq "+"){
				$dis = $ixpos==2? $b->[$ixpos]-$a->[$ixpos]: $b->[$ixpos]-$a->[$ixpos]+1;
			}else{
				$dis = $ixpos==2? $a->[$ixpos]-$b->[$ixpos]: $a->[$ixpos]-$b->[$ixpos]+1;
			}

            my $eff_dis = &score_single($dis, 1, $mind,$maxd, $dmin, $dmax);
		    if (!defined $a->[10]){
                ($a->[10], $a->[11], $a->[12]) = &efficiency($a->[4], $a->[6], $opt_tm, $min_tm, $a->[5]);
            }
		    if (!defined $b->[10]){
                ($b->[10], $b->[11], $b->[12]) = &efficiency($b->[4], $b->[6], $opt_tm, $min_tm, $b->[5]);
            }
			my $eff=$a->[10]*$b->[10]*$eff_dis;
			next if($eff<$min_eff);
			$prodn++;
			#my $prob=$dis."/".join(",", $a->[8], $b->[8], $a->[9], $b->[9]);
			my $prob=$dis."/".join(",", $a->[0], $a->[8] ? $p2 : $p1, $a->[1].$a->[3], $a->[6], sprintf("%.2f", $a->[4]), $b->[8] ? $p2 : $p1, $b->[1].$b->[3], $b->[6], sprintf("%.2f", $b->[4]));
			$aprod->{$tid1}{$tid2}{$prob}=$eff;
			if($prodn==$Max_prodn){ # p1_p2, p1_p1, p2_p2
				$aprod->{$tid1}{$tid2}{"Max"}=1; ## mark:reach max production num
				return $prodn;
			}
		}
		for (my $j = $i - 1; ; $j--){ # look backward
			last if ($j < 0);
			my $b = $all_align_info[$j];
			last if ($b->[$ixpos] < $pmin);  # exceed valid dist range
			last if ($b->[0] ne $a->[0]);
			next if ($b->[1] ne $sd2);

			# $a and $b create effective product
			my $dis = 0;
			if($a->[1] eq "+"){
				$dis = $ixpos==2? $b->[$ixpos]-$a->[$ixpos]: $b->[$ixpos]-$a->[$ixpos]+1;
			}else{
				$dis = $ixpos==2? $a->[$ixpos]-$b->[$ixpos]: $a->[$ixpos]-$b->[$ixpos]+1;
			}

			my $eff_dis=&score_single($dis, 1, $mind,$maxd, $dmin, $dmax);
		    if (!defined $a->[10]){
                ($a->[10], $a->[11], $a->[12]) = &efficiency($a->[4], $a->[6], $opt_tm, $min_tm, $a->[5]);
            }
		    if (!defined $b->[10]){
                ($b->[10], $b->[11], $b->[12]) = &efficiency($b->[4], $b->[6], $opt_tm, $min_tm, $b->[5]);
            }
            my $eff=$a->[10]*$b->[10]*$eff_dis;
			next if($eff<$min_eff);
			$prodn++;
			#my $prob=$dis."/".join(",", $a->[8], $b->[8], $a->[9], $b->[9]);
			my $prob=$dis."/".join(",", $a->[0], $a->[8] ? $p2 : $p1, $a->[1].$a->[3], $a->[6], sprintf("%.2f", $a->[4]), $b->[8] ? $p2 : $p1, $b->[1].$b->[3], $b->[6], sprintf("%.2f", $b->[4]));
			$aprod->{$tid1}{$tid2}{$prob}=$eff;
			if($prodn==$Max_prodn){ # p1_p2, p1_p1, p2_p2
				$aprod->{$tid1}{$tid2}{"Max"}=1; ## mark:reach max production num
				return $prodn;
			}

		}
	}
	
	return ($prodn);
}

# [$chr, $strand, $pos3, $pos5, $tm, $end3_base, $mvisual, $tid, $type, undef, undef, undef, undef]
sub caculate_product_evaluation{
	my ($all_align_info, $ptype, $aprod, $arecord, $PCRsize, $opt_tm, $min_tm, $mind, $maxd, $min_eff, $Max_prodn, $etype, $eval_all, $tp0)=@_;	
	my $prodn=0;

	##                Nested		 face-to-face         back-to-back
##P1			 	  -->| 		   |-->                    |-->
##P2min		  -->|...      		    <--|             <--|...  
##P2max			      -->| 		       ...<--|                 ...<--|   
	my ($ixpos, $dmin, $dmax);
	if($ptype eq "Nested"){
		$ixpos=2; #pos3
		#Nested: P1 --->|   x        distance range: 5,30 (mind>0)
		#          P2 --->| x

		#Nested: P2 --->|   x        distance range: -30,-5 (mind<0)
		#          P1 --->| x
		if($mind>0){
			$dmin=0;
			$dmax=$PCRsize;
		}else{
			$dmin=-1*$PCRsize;
			$dmax=0;
		}
	}elsif($ptype eq "back-to-back"){
			$ixpos=3; #pos5
			$dmin=-1*$PCRsize+2*20; #20: primer min len
			$dmax=$PCRsize;
	}else{ ##face-to-face
			$ixpos=3; #pos5
			$dmin=30;
			$dmax=$PCRsize;
	}

    my $is_nested_mode = $ptype eq "Nested";

	#sort according to chr and pos
	@$all_align_info = sort { $a->[0] cmp $b->[0] or $a->[$ixpos] <=> $b->[$ixpos]} @$all_align_info;

	# iter all alignments and found effective products
	for (my $i = 0; $i < @$all_align_info; ++$i){  # O(n+m) 
		my $a = $all_align_info->[$i];
		next if($a->[8] eq "P");
		next if($etype eq "SinglePlex" && !defined $eval_all && $a->[8] ne $tp0); ## only evalue one type 
		my $tid1 = $a->[7];
        my $sd = $a->[1];
        my $pos = $a->[$ixpos];
        my $sd2 = $sd;
        my $is_rev = $sd eq "-";
        if (!$is_nested_mode){
            $sd2 = ("+", "-")[1-$is_rev];
        }
        my $pmin = $is_rev ? $pos - $dmax : $pos + $dmin;;
        my $pmax = $is_rev ? $pos - $dmin : $pos + $dmax;
        #my ($sd2, $pmin, $pmax)=&primer2_scope(\$ptype, $a->[1], $a->[$ixpos], $dmin, $dmax);
        for (my $j = $i + 1; ; $j++){ # look forward
			last if ($j >= @$all_align_info); # no more alignments
			my $b = $all_align_info->[$j];
			last if ($b->[$ixpos] > $pmax);  # exceed valid dist range
			last if ($b->[0] ne $a->[0]);  # diff chr
			next if($b->[8] eq "P");
			next if ($b->[1] ne $sd2); # not valid strand
			next if(!defined $eval_all && $a->[8] eq $b->[8]);

			my $tid2 = $b->[7];
			next if($etype eq "SinglePlex" && $tid2 ne $tid1); ## SinglePlex: not evalue product between different tid
			# $a and $b create effective product
			my $dis;
			if($a->[1] eq "+"){
				$dis = $ixpos==2? $b->[$ixpos]-$a->[$ixpos]: $b->[$ixpos]-$a->[$ixpos]+1;
			}else{
				$dis = $ixpos==2? $a->[$ixpos]-$b->[$ixpos]: $a->[$ixpos]-$b->[$ixpos]+1;
			}

            my $eff_dis = &score_single($dis, 1, $mind,$maxd, $dmin, $dmax);
		    if (!defined $a->[10]){
                ($a->[10], $a->[11], $a->[12]) = &efficiency($a->[4], $a->[6], $opt_tm, $min_tm, $a->[5]);
            }
		    if (!defined $b->[10]){
                ($b->[10], $b->[11], $b->[12]) = &efficiency($b->[4], $b->[6], $opt_tm, $min_tm, $b->[5]);
            }
			my $eff=$a->[10]*$b->[10]*$eff_dis;
			next if($eff<$min_eff);
			$prodn++;
			
			if ($tid1 lt $tid2){
				my $prob=$dis."/".join(",", $a->[0], $tid1."-".$a->[8], $a->[1].$a->[3], $a->[6], sprintf("%.2f", $a->[4]), $tid2."-".$b->[8], $b->[1].$b->[3], $b->[6], sprintf("%.2f", $b->[4]));
				$aprod->{$tid1}{$tid2}{$prob}=$eff;
			}else{
				my $prob=$dis."/".join(",", $b->[0], $tid2."-".$b->[8], $b->[1].$b->[3], $b->[6], sprintf("%.2f", $b->[4]), $tid1."-".$a->[8], $a->[1].$a->[3], $a->[6], sprintf("%.2f", $a->[4]));
				$aprod->{$tid2}{$tid1}{$prob}=$eff;
			}
			
			if($prodn==$Max_prodn){ # p1_p2, p1_p1, p2_p2
				$aprod->{$tid1}{$tid2}{"Max"}=1; ## mark:reach max production num
				return $prodn;
			}
		}
		for (my $j = $i - 1; ; $j--){ # look backward
			last if ($j < 0);
			my $b = $all_align_info->[$j];
			last if ($b->[$ixpos] < $pmin);  # exceed valid dist range
			last if ($b->[0] ne $a->[0]);
			next if($b->[8] eq "P");
			next if ($b->[1] ne $sd2); # not valid strand
			next if(!defined $eval_all && $a->[8] eq $b->[8]);
			my $tid2 = $b->[7];
			next if($etype eq "SinglePlex" && $tid2 ne $tid1); ## SinglePlex: not evalue product between different tid

			# $a and $b create effective product
			my $dis = 0;
			if($a->[1] eq "+"){
				$dis = $ixpos==2? $b->[$ixpos]-$a->[$ixpos]: $b->[$ixpos]-$a->[$ixpos]+1;
			}else{
				$dis = $ixpos==2? $a->[$ixpos]-$b->[$ixpos]: $a->[$ixpos]-$b->[$ixpos]+1;
			}

			my $eff_dis=&score_single($dis, 1, $mind,$maxd, $dmin, $dmax);
		    if (!defined $a->[10]){
                ($a->[10], $a->[11], $a->[12]) = &efficiency($a->[4], $a->[6], $opt_tm, $min_tm, $a->[5]);
            }
		    if (!defined $b->[10]){
                ($b->[10], $b->[11], $b->[12]) = &efficiency($b->[4], $b->[6], $opt_tm, $min_tm, $b->[5]);
            }
            my $eff=$a->[10]*$b->[10]*$eff_dis;
			next if($eff<$min_eff);
			$prodn++;
			
			if ($tid1 lt $tid2){
				my $prob=$dis."/".join(",", $a->[0], $tid1."-".$a->[8], $a->[1].$a->[3], $a->[6], sprintf("%.2f", $a->[4]), $tid2."-".$b->[8], $b->[1].$b->[3], $b->[6], sprintf("%.2f", $b->[4]));
				$aprod->{$tid1}{$tid2}{$prob}=$eff;
			}else{
				my $prob=$dis."/".join(",", $b->[0], $tid2."-".$b->[8], $b->[1].$b->[3], $b->[6], sprintf("%.2f", $b->[4]), $tid1."-".$a->[8], $a->[1].$a->[3], $a->[6], sprintf("%.2f", $a->[4]));
				$aprod->{$tid2}{$tid1}{$prob}=$eff;
			}
			
			if($prodn==$Max_prodn){ # p1_p2, p1_p1, p2_p2
				$aprod->{$tid1}{$tid2}{"Max"}=1; ## mark:reach max production num
				return $prodn;
			}

		}
	}
	
	return ($prodn);
}


#my ($sd, $pmin, $pmax)=&primer2_scope($type, $strand,$pos, $dis_min, $dis_max);
sub primer2_scope{
	my ($type, $strand, $pos, $dis_min, $dis_max)=@_;
	my ($sd, $pmin, $pmax);
	if($$type ne "Nested"){ ## strand: opposite; Dis: from pos5
		if($strand eq "+"){
			$sd="-";
			$pmin=$pos+$dis_min;
			$pmax=$pos+$dis_max;
		}else{
			$sd="+";
			$pmin=$pos-$dis_max;
			$pmax=$pos-$dis_min;
		}
	}else{ ## Nested: Dis from pos3
		$sd=$strand;
		if($strand eq "+"){
			$pmin=$pos+$dis_min; #dis <0
			$pmax=$pos+$dis_max;
		}else{
			$pmin=$pos-$dis_max;
			$pmax=$pos-$dis_min;
		}
	}
	return ($sd, $pmin, $pmax);
}



sub efficiency_probe{
	my ($tm, $dis, $PCRsize, $opt_tm_probe, $min_tm)=@_;
	# tm eff
	my $opt=$opt_tm_probe;
	my @etm=($opt-1, $opt+100, $min_tm-0.1, $opt+100, $min_tm-0.1, $opt+100); 
	my $eff_tm = &score_growth_curve($tm, 1, @etm);
	$eff_tm = $eff_tm>0? $eff_tm: 0;
	
	# dis to primer 3end
	my @dis=(1,10,-5,$PCRsize);
	my $eff_dis=&score_single($dis, 1, @dis);
	$eff_dis=$eff_dis>0? $eff_dis: 0;
	return $eff_tm*$eff_dis;
}



sub efficiency{
	my ($tm, $mvisual, $opt_tm, $min_tm, $end3b)=@_;
	
	# tm eff
	my @etm=($opt_tm, $opt_tm+100, $min_tm-0.1, $opt_tm+100, $min_tm-0.1, $opt_tm+100); 
	my $eff_tm = &score_growth_curve($tm, 1, @etm);
	$eff_tm = $eff_tm>0? $eff_tm: 0;

	# mismatch pos to 3end
	my @mis_pos;
	my $tlen = &get_3end_mismatch($mvisual, \@mis_pos);
	my $eff_end = 1;
	if(scalar @mis_pos>0){
		my @mis_end=($tlen,100,1,100); #end: 0-$tlen/2 => score: 0-fulls
		if($mis_pos[0]==1 && ($end3b eq "C" || $end3b eq "G")){
			@mis_end=($tlen,100,0.95,100);
		}
		for(my $i=0; $i<@mis_pos; $i++){
			$eff_end *= &score_single($mis_pos[$i], 1, @mis_end);
		}
	}
	
	return ($eff_tm*$eff_end, $eff_tm, $eff_end);
}

sub get_3end_mismatch{
	my ($mv, $apos)=@_;
	
	$mv=~s/\^+/\*\*/; ## Del ==> 2 mismatch
	$mv=~s/\-+/\*\*/; ## Insert ==> 2 mismatch
	my @unit = split //, $mv;
	my $num=0;
	for(my $i=$#unit; $i>=0; $i--){
		$num++;
		if($unit[$i] eq "*" || $unit[$i] eq "#"){
			push @{$apos}, $num;
		}
	}
	return $num;
}

sub get_highest_bound{
	my ($abound, $maxn, $type)=@_;
	my %bound=%{$abound};
	my $is_max=0;
	if(exists $bound{"Max"}){
		$is_max=1;
		delete $bound{"Max"};
	}
	my @binfo=sort{$bound{$b} <=> $bound{$a}} keys %bound;
	my $bnum = scalar @binfo;
	my $n=0;
	my @bvalue;
	my @binfos;
	foreach $binfo(@binfo){
		push @binfos, $binfo;
		if($type eq "Eff"){
			push @bvalue, sprintf("%.5f", $bound{$binfo});
		}else{
			push @bvalue, sprintf("%.2f", $bound{$binfo});
		}
		$n++;
		last if($n==$maxn);
	}
	if($is_max==1){
		$bnum.="+";
	}
	return ($bnum, \@bvalue, \@binfos);
}



1;

