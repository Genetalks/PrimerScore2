require "$Bin/score.pm";
require "$Bin/math.pm";

my $end_len=10;
my @mis_end=(10,100,0,100); #end: 0-10 => score: 0-fulls
my $Max_prodn=50;

#my $prob=$dis."/".join(",", $chr, $sd.$pos, $mvisual1,sprintf("%.2f",$tm1), $sd2.$pos2,$mvisual2,sprintf("%.2f",$tm2));
sub probe_bounds_on_products{
	my ($pb, $abound, $aprod, $aresult, $PCRsize, $opt_tm_probe)=@_;
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
				my ($pos3, $pos5, $tm, $end_match, $mvisual)=@{$bds[$i]};
				if($pos5>$minp-10 && $pos5<$maxp+10){
					my $dis=&min($pos5-$minp, $maxp-$pos5);
					my $eff=&efficiency_probe($tm, $dis, $PCRsize, $opt_tm_probe) * $aprod->{$prod};
					my $bpd=$dis."/".join(",", $chr, $sd.$pos5,$mvisual,sprintf("%.2f",$tm)).":".$info;
					$aresult->{$pb}{$bpd}=$eff;
					$n++;
				}
			}
		}
	}
	
	return $n;
}

#push @{$bound{$id}{$chr}{$strand}}, [$pos3, $pos5, $tm, $end_match, $mvisual];
sub caculate_product{
	my ($tid1, $p1, $tid2, $p2, $abound1, $abound2, $ptype, $aprod, $arecord, $PCRsize, $opt_tm, $mind, $maxd, $min_eff)=@_;
	my %bdp1=%{$abound1};
	my %bdp2=%{$abound2};
	my $prodn=0;
	foreach my $chr(keys %bdp1){
		foreach my $sd(keys %{$bdp1{$chr}}){
			my @bds=@{$bdp1{$chr}{$sd}};
			for(my $i=0; $i<@bds; $i++){
				my ($tm1, $end_match1, $mvisual1)=@{$bds[$i]}[2..4];
##                Nested		 face-to-face         back-to-back
##P1			 	  -->| 		   |-->                    |-->
##P2min		  -->|...      		    <--|             <--|...  
##P2max			      -->| 		       ...<--|                 ...<--|   
				my ($ixpos, $dmin, $dmax);
				if($ptype eq "Nested"){
					$ixpos=0; #pos3
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
					$ixpos=1; #pos5
					$dmin=-1*$PCRsize+2*20; #20: primer min len
					$dmax=$PCRsize;
				}else{ ##face-to-face
					$ixpos=1; #pos5
					$dmin=30;
					$dmax=$PCRsize;
				}
				my $pos=$bds[$i][$ixpos];
				my ($sd2, $pmin, $pmax)=&primer2_scope($ptype,$sd,$pos,$dmin,$dmax);
				next if(!exists $bdp2{$chr}{$sd2});
				my @bds2=@{$bdp2{$chr}{$sd2}};
				for(my $j=0; $j<@bds2; $j++){
					my ($tm2, $end_match2, $mvisual2)=@{$bds2[$j]}[2..4];
					my $pos2 = $bds2[$j][$ixpos];
					if($pos2>=$pmin && $pos2<=$pmax){
						##
						my ($eff1, $eff2);
						my $b1=join(",", $p1, $chr, $sd, $bds[$i][0]);
						if(exists $arecord->{"eff"}->{$b1}){
							$eff1=$arecord->{"eff"}->{$b1};
						}else{
							$eff1=&efficiency($tm1, $mvisual1, $opt_tm);
							$arecord->{"eff"}->{$b1}=$eff1;
						}
						my $b2=join(",", $p2, $chr, $sd2, $bds2[$j][0]);
						if(exists $aeff->{$b2}){
							$eff2=$aeff->{$b2};
						}else{
							$eff2=&efficiency($tm2, $mvisual2, $opt_tm);
							$aeff->{$b2}=$eff2;
						}
						my $dis=$pos2-$pos;
						if($sd eq "-"){
							$dis=$dis*(-1);
						}
						my @rsize=($mind,$maxd, $dmin, $dmax);
						my $eff_dis=&score_single($dis, 1, @rsize);
						my $eff=$eff1*$eff2*$eff_dis;
						my $prob=$dis."/".join(",", $chr, $p1, $sd.$pos, $mvisual1,sprintf("%.2f",$tm1), $p2, $sd2.$pos2,$mvisual2,sprintf("%.2f",$tm2));
						next if($eff<$min_eff);
						$aprod->{$tid1}{$tid2}{$prob}=$eff;
					#	print join("\t", $eff, $prob),"\n";
					#	print Dumper %{$aprod};
						$prodn++;
						if($prodn==$Max_prodn){
							$aprod->{$tid1}{$tid2}{"Max"}=1; ## mark:reach max production num
							return $prodn;
						}
					}
				}
			}
		}
	}
	return ($prodn);
}

#my ($sd, $pmin, $pmax)=&primer2_scope($type, $strand,$pos, $dis_min, $dis_max);
sub primer2_scope{
	my ($type, $strand, $pos, $dis_min, $dis_max)=@_;
	my ($sd, $pmin, $pmax);
	if($type ne "Nested"){ ## strand: opposite; Dis: from pos5
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
	my ($tm, $dis, $PCRsize, $opt_tm_probe)=@_;
	# tm eff
	my $opt=$opt_tm_probe;
	my @etm=($opt-1, $opt+100, 45, $opt+100); 
	my $eff_tm = &score_single($tm, 1, @etm);
	$eff_tm = $eff_tm>0? $eff_tm: 0;
	
	# dis to primer 3end
	my @dis=(1,10,-5,$PCRsize);
	my $eff_dis=&score_single($dis, 1, @dis);
	$eff_dis=$eff_dis>0? $eff_dis: 0;
	return $eff_tm*$eff_dis;
}



sub efficiency{
	my ($tm, $mvisual, $opt_tm)=@_;
	# tm eff
	my @etm=($opt_tm-1, $opt_tm+100, 40, $opt_tm+100); 
	my $eff_tm = &score_single($tm, 1, @etm);
	$eff_tm = $eff_tm>0? $eff_tm: 0;

	# mismatch pos to 3end
	my @mis_pos;
	&get_3end_mismatch($mvisual, \@mis_pos, $end_len);
	my $eff_end = 1;
	for(my $i=0; $i<@mis_pos; $i++){
		$eff_end *= &score_single($mis_pos[$i], 1, @mis_end);
	}
	return $eff_tm*$eff_end;
}

sub get_3end_mismatch{
	my ($mv, $apos, $max_end)=@_;
	
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
}

sub get_highest_bound{
	my ($abound, $maxn)=@_;
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
		push @bvalue, sprintf("%.2f",$bound{$binfo});
		$n++;
		last if($n==$maxn);
	}
	if($is_max==1){
		$bnum.="+";
	}
	return ($bnum, \@bvalue, \@binfos);
}



1;

