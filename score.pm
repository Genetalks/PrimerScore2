sub probe_meth_score{
	my ($opt_tm, $len, $tm, $gc, $hairpin, $dimert, $dimers, undef, undef, $snp, $poly, $cpgs, $cs, $bnumtm, $is_G5, $CGd)=@_;
	my @tm = ($opt_tm, $opt_tm+1, $opt_tm-3, $opt_tm+5);
	my @gc = (0.48, 0.6, 0.4, 0.7);
	my @self = (-50, 40, -50, 55); ## self tm
	my @CGd = (0.1, 1, 0, 1);
	my @G5 = (0, 0, 0, 0.5);
	my @cpg=(4,100,2,100,0,100);
	my @cs=(8,100,4,100,0,100);
	my $fulls=10;
	
	my $stm=int(&score_single($tm, $fulls, @tm)+0.5);
	my $sgc=int(&score_single($gc, $fulls, @gc)+0.5);
	my $sself=int(&score_single($hairpin, $fulls, @self)+0.5);
	my ($snpv)=split /:/, $snp;	
	my $ssnp = int(&SNP_score($snpv, $len, "Probe")*$fulls +0.5);
	my $spoly = int(&poly_score_curve($poly, $len, "Probe")*$fulls +0.5);
	my $sCGd=int(&score_single($CGd, $fulls, @CGd)+0.5);
	my $sG5=int(&score_single($is_G5, $fulls, @G5)+0.5);
	#specificity: bound
	my $cpgn=scalar(split /,/, $cpgs);
	my $csn=scalar(split /,/, $cs);
	my $scpg=int(&score_growth_curve($cpgn, $fulls, @cpg)+0.5);
	my $scs=int(&score_growth_curve($csn, $fulls, @cs)+0.5);
	my $sbound=$fulls;
	if(defined $bnumtm){
		my ($bnum, $btm)=split /\|/, $bnumtm;
		$sbound = &bound_score($bnum, $btm, $fulls, "Probe_Tm");
	}
	
	my @score = ($stm, $sgc, $sself, $sCGd, $sG5, $ssnp, $spoly, $sbound, $scpg, $scs);
	my @weight =( 1.5,   1,     1,      1,    1,    0.8,    0.2,   0.5,  1.5,  1.5);
	my $sadd=0;
	for(my $i=0; $i<@score; $i++){
#		$score[$i]=$score[$i]<0? 0: $score[$i];
		$sadd+=$weight[$i]*$score[$i];
	}

	my $score_info=join(",", @score);
	return($sadd, $score_info);
}



sub probe_oligo_score{
	my ($opt_tm, $len, $tm, $gc, $hairpin, $dimert, $dimers, undef, undef, $snp, $poly, $bnumtm, $is_G5, $CGd)=@_;

	my @tm = ($opt_tm, $opt_tm+1, $opt_tm-3, $opt_tm+5);
	my @gc = (0.48, 0.6, 0.4, 0.7);
	my @self = (-50, 40, -50, 55); ## self tm
	my @CGd = (0.1, 1, 0, 1);
	my @G5 = (0, 0, 0, 0.5);
	my $fulls=10;
	
	my $stm=int(&score_single($tm, $fulls, @tm)+0.5);
	my $sgc=int(&score_single($gc, $fulls, @gc)+0.5);
	my $sself=int(&score_single($hairpin, $fulls, @self)+0.5);
	my ($snpv)=split /:/, $snp;	
	my $ssnp = int(&SNP_score($snpv, $len, "Probe")*$fulls +0.5);
	my $spoly = int(&poly_score_curve($poly, $len, "Probe")*$fulls +0.5);
	my $sCGd=int(&score_single($CGd, $fulls, @CGd)+0.5);
	my $sG5=int(&score_single($is_G5, $fulls, @G5)+0.5);
	#specificity: bound
	my $sbound=$fulls;
	if(defined $bnumtm){
		my ($bnum, $btm)=split /\|/, $bnumtm;
		$sbound = &bound_score($bnum, $btm, $fulls, "Probe_Tm");
	}
	my @score = ($stm, $sgc, $sself, $sCGd, $sG5, $ssnp, $spoly, $sbound);
	my @weight =( 2,   1.5,     1,      1,    1,    1.5,    0.5,      1.5);
	my $sadd=0;
	for(my $i=0; $i<@score; $i++){
#		$score[$i]=$score[$i]<0? 0: $score[$i];
		$sadd+=$weight[$i]*$score[$i];
	}

	my $score_info=join(",", @score);
	return($sadd, $score_info);
}

## score_growth_curve
sub primer_meth_score{
	my ($opt_tm, $len, $tm, $gc, $hairpin, $dimert, $dimers, $nendA, $enddG, $snp, $poly, $cpgs, $cs, $bnumtm)=@_;
	my $fulls = 10;
	my @nendA = (1, 1, -1, 3, -1, 8);
	my @enddG = (-9,-7,-12,-6.2,-14,-5);
	my @gc = (0.52, 0.6, 0.42, 0.7, 0.32, 0.8);
	my @tm = ($opt_tm, $opt_tm+1, $opt_tm-2, $opt_tm+6, $opt_tm-4, $opt_tm+11);
	my @self = (-50, 47, -50, 52, -50, 57); ## self tm
	my @cpg=(4,100,2,100,0,100);
	my @cs=(8,100,4,100,0,100);

	my $snendA=int(&score_growth_curve($nendA, $fulls, @nendA)+0.5);
	my $senddG=int(&score_single($enddG, $fulls, @enddG)+0.5);
	my $stm=int(&score_growth_curve($tm, $fulls, @tm)+0.5);
	my $sgc=int(&score_growth_curve($gc, $fulls, @gc)+0.5);
	my $sself=int(&score_growth_curve($hairpin, $fulls, @self)+0.5);
	my ($snpv)=split /:/, $snp; 
	my $ssnp = int(&SNP_score($snpv, $len, "Primer")*$fulls +0.5);
	my $spoly = int(&poly_score($poly, $len, "Primer")*$fulls +0.5);
#	int(&cpgs_score("0", 24) *$fulls +0.5);
#	int(&cpgs_score("3", 24) *$fulls +0.5);
#	int(&cpgs_score("5", 24) *$fulls +0.5);
#	int(&cpgs_score("10", 24) *$fulls +0.5);
#	int(&cpgs_score("15", 24) *$fulls +0.5);
#	int(&cpgs_score("20", 24) *$fulls +0.5);
#	int(&cpgs_score("0,1", 24) *$fulls +0.5);
#	int(&cpgs_score("5,8,11,13,15", $len) *$fulls +0.5);
#	int(&cpgs_score("1,11,17", $len) *$fulls +0.5);
#	die;
	my $scpgs = int(&cpgs_score($cpgs, $len) *$fulls +0.5);
	my $scs = int(&cpgs_score($cs, $len) *$fulls +0.5);
	my $sbound=10;
	if(defined $bnumtm){
		my ($bnum, $btm)=split /\|/, $bnumtm;
		$sbound=&bound_score($bnum, $btm, $fulls, "Primer_Tm");
	}
	my @score = ($stm, $sgc, $sself, $snendA, $senddG, $ssnp, $spoly, $sbound,$scpgs, $scs);
	my @weight =(1,    1.5,    1,    0.5,      1,       1,     0.5,   0.5,    1.5,   1.5);
	my $sadd=0;
	for(my $i=0; $i<@score; $i++){
#		$score[$i]=$score[$i]<0? 0: $score[$i];
		$sadd+=$weight[$i]*$score[$i];
	}
	my $score_info=join(",", @score);
	return ($sadd, $score_info);
}


## score_growth_curve
sub primer_oligo_score{
	my ($opt_tm, $len, $tm, $gc, $hairpin, $dimert, $dimers, $nendA, $enddG, $snp, $poly, $bnumtm)=@_;
	my $fulls = 10;
	my @nendA = (1, 1, -1, 3, -1, 8);
	my @enddG = (-9,-7,-12,-6.2,-14,-5);
	my @gc = (0.52, 0.6, 0.42, 0.7, 0.32, 0.8);
	my @tm = ($opt_tm, $opt_tm+1, $opt_tm-2, $opt_tm+6, $opt_tm-4, $opt_tm+11);
	my @self = (-50, 47, -50, 52, -50, 57); ## self tm

	my $snendA=int(&score_growth_curve($nendA, $fulls, @nendA)+0.5);
	my $senddG=int(&score_single($enddG, $fulls, @enddG)+0.5);
	my $stm=int(&score_growth_curve($tm, $fulls, @tm)+0.5);
	my $sgc=int(&score_growth_curve($gc, $fulls, @gc)+0.5);
	my $sself=int(&score_growth_curve($hairpin, $fulls, @self)+0.5);
	my ($snpv)=split /:/, $snp; 
	my $ssnp = int(&SNP_score($snpv, $len, "Primer")*$fulls +0.5);
	my $spoly = int(&poly_score($poly, $len, "Primer")*$fulls +0.5);
	my $sbound= $fulls;
	if(defined $bnumtm){
		my ($bnum, $btm)= split /\|/, $bnumtm;
		$sbound=&bound_score($bnum, $btm, $fulls, "Primer_Tm");
	}
	my @score = ($stm, $sgc, $sself, $snendA, $senddG, $ssnp, $spoly, $sbound);
	my @weight =(1.5,     2,    1.5,    0.5,      1,    1.5,   1.5,   0.5);
	my $sadd=0;
	for(my $i=0; $i<@score; $i++){
#		$score[$i]=$score[$i]<0? 0: $score[$i];
		$sadd+=$weight[$i]*$score[$i];
	}
	my $score_info=join(",", @score);
	return ($sadd, $score_info);
}



sub bound_score{
	my ($bnum, $bvalue, $fulls, $type)=@_;
	#bvalue: tm, Eff
	my $sbound;
	if($bnum==0){## too low tm
		return $fulls*(-1);
	}elsif($bnum==1){
		$sbound=$fulls;
	}else{
		my @value=split /,/, $bvalue;
		if(scalar @value==1){
			return $fulls;
		}
		my (@bvalue, $evalue, $enum, $etotal);
		if($type eq "Primer_Tm"){
			my @bnum = (1,5,1,500); #bound number
			@bvalue = (0,$value[0]*0.65,0,$value[0]); #bound sec value
			$evalue = &score_single($value[1], 0.5, @bvalue);
			$enum = &score_single($bnum, 0.5, @bnum);
			$etotal = $evalue+$enum;
		}elsif($type eq "Probe_Tm"){
			my @bnum = (1,5,1,100); #bound number
			@bvalue = (0,$value[0]*0.65,0,$value[0]); #bound sec value
			$evalue = &score_single($value[1], 0.7, @bvalue);
			$enum = &score_single($bnum, 0.3, @bnum);
			$etotal = $evalue+$enum;
		}else{ #Eff
			@bvalue = (0,0.00001,0,0.1);
			if(!defined $value[1]){
				print STDERR join("\t", $bnum, $bvalue),"\n";
			}
			$evalue = &score_single($value[1], 1, @bvalue);
			#$enum = &score_single($bnum, 1, @bnum);
			$etotal = $evalue;
		}
		$sbound = int($etotal*$fulls+0.5);
		$sbound = $sbound>0? $sbound: 0;
	}
	return $sbound;
}

#Primer: AAA=>10,AAAA=>8,AAAAA=>6,
sub poly_score_curve{
	my ($info, $len, $type)=@_;
	return 1 if($info eq "NA");
	my @polys = split /,/, $info;
	my $score=1;
	my @dprobe=($len*0.3, $len*0.5, 0, $len*0.5);
	my @dprimer=(8,$len,0,$len);
	my @polylen1=(0,2.5,0,4,0,10);
	my @polylen2=(0,2.5,0,6,0,15);
	for(my $i=0; $i<@polys; $i++){
		my ($d3, $t, $l)=$polys[$i]=~/(\d+)(\w)(\d+)/;
		my $s;
		my @polylen=@polyleno;
		if(($type eq "Probe" && $t=~/CG/) || ($type eq "Primer" && $t=~/AT/)){
			@polylen=@polylen1;
		}else{
			@polylen=@polylen2;
		}
		$slen = (1-&score_growth_curve($l,1,@polylen));
		if($type eq "Probe"){
			my $d5=$len-$d3-1;
			my $de=$d5>$d3?$d3:$d5;
			$sdis = (2 - &score_single($de, 1, @dprobe));
		}else{
			$sdis = (2 - &score_single($d3, 1, @dprimer));
		}
		$score-=$sdis*$slen;
	}
	return $score;
}


#Primer: AAA=>10,AAAA=>8,AAAAA=>6,
sub poly_score{
	my ($info, $len, $type)=@_;
	return 1 if($info eq "NA");
	my @polys = split /,/, $info;
	my $score=1;
	my @dprobe=($len*0.3, $len*0.5, 0, $len*0.5);
	my @dprimer=(8,$len,0,$len);
	my @polylen1=(0,2.5,0,4);
	my @polylen2=(0,2.5,0,6);
	for(my $i=0; $i<@polys; $i++){
		my ($d3, $t, $l)=$polys[$i]=~/(\d+)(\w)(\d+)/;
		my $s;
		my @polylen=@polyleno;
		if(($type eq "Probe" && $t=~/CG/) || ($type eq "Primer" && $t=~/AT/)){
			@polylen=@polylen1;
		}else{
			@polylen=@polylen2;
		}
		$slen = (1-&score_single($l,1,@polylen));
		if($type eq "Probe"){
			my $d5=$len-$d3-1;
			my $de=$d5>$d3?$d3:$d5;
			$sdis = (2 - &score_single($de, 1, @dprobe));
		}else{
			$sdis = (2 - &score_single($d3, 1, @dprimer));
		}
		$score-=$sdis*$slen;
	}
	return $score;
}

sub get_poly_value{
	my ($seq)=@_;
	$seq = reverse $seq;
	my @unit = split //, $seq;
	
	my $min_plen = 3; # min poly length
	my @polys;
	my $last = $unit[0];
	my $poly = $unit[0];
	my $dis = 0; #dis to the 3 end
	for(my $i=1; $i<@unit; $i++){
		if($unit[$i] eq $last){
			$poly.=$unit[$i];
		}else{
			if(length($poly)>=$min_plen){
				push @polys, [$poly, $dis];
			}
			$dis = $i;
			$last = $unit[$i];
			$poly=$last;
		}
	}
	if(length($poly)>=$min_plen){
		push @polys, [$poly, $dis];
	}
	
	my $value1 = 0; ## end poly impact
	my $pslen = 0;
	my $psnum = 0;
	for(my $i=0; $i<@polys; $i++){
		my $l = length($polys[$i][0]);
		my $d = $polys[$i][1];
		$pslen+=$l;
		$psnum++;
		if($i==0){
			#$value1+=($l-$min_plen+1)*(4-$d)*3; # end poly impact
			$value1+=($l-$min_plen+1)*(6-$d)*2; # end poly impact
		}
	}
	$value1 = $value1>=0? $value1: 0;

	## other impact
	my $value2 = 0.2*$pslen*(6-$psnum);
	$value2 = $value2>=0? $value2: 0;
	my $value = $value1 + $value2;
	return $value;
}

sub cpgs_score{
	my ($info, $len)=@_;
	my @units=split /,/, $info;
	my $qua=int($len/4);
	my $min=-2*$qua;
	my @pos=($min,$min,$min,$len,$min,$len);
	## first score 
	my $score0=0;
	for(my $i=0; $i<@units; $i++){
		$score0+=&score_growth_curve($units[$i], 100, @pos);
	}

	# score final
	@score=(150,1000,75,1000,0,1000);
	my $score=&score_growth_curve($score0, 1, @score);
	return $score;
}

sub SNP_score{
	my ($info, $len, $type)=@_;
	return 1 if($info eq "NA");
	my @snps = split /,/, $info;
	my $score=1;
	my @dprobe=($len*0.4, $len*0.5, 0, $len*0.5);
	my @dprimer=(10,$len,4,$len);
	my @snplen=(0,0,0,10); #snp length: SNP-1, Indel-indel len
	for(my $i=0; $i<@snps; $i++){
		my ($d3, $t, $l)=$snps[$i]=~/(\d+)(\w)(\d+)/;
		my $s;
		if($type eq "Probe"){
			my $d5=$len-$d3-1;
			my $de=$d5>$d3?$d3:$d5;
			$s = &score_single($de,1,@dprobe)*&score_single($l,1,@snplen);
		}else{
			$s = &score_single($d3, 1, @dprimer)*&score_single($l,1,@snplen);
		}
		$score*=$s;
	}
	return $score;
}


### score accord growth curve
## growth curve(a,b,K)=(1,1,1) is a symmetrical curve,  y--(0,1), y tend to 0 when x<-5, y tend to 1 when x>5.
## minl is the limit min value when y tend to 0
## maxl is the limit max value when y tend to 1
#    |<--curve1->|   |<--curve2->|
#----|------|----|---|-----|-----|---
# minl     min minb maxb  max   maxl

sub score_growth_curve{
	my ($v, $score, $minb, $maxb, $min, $max, $minl, $maxl)=@_;
	if(!defined $v){
		print STDERR join("\t", $v, $score, $minb, $maxb, $min, $max, $minl, $maxl),"\n";
	}
	return $score if($v eq "NULL");
	if(!defined $minb){
		print STDERR join("\t", $v, $score, $minb, $maxb, $min, $max),"\n";
	}
	my $disl = $minb - $min;
	my $disr = $max - $maxb;
	if($disl<0 || $disr<0){
		print "Wrong Score set: $minb, $maxb, $min, $max\n";
		die;
	}
	##  growth curve
	my ($K, $a, $b, $e)=(1,1,1,2.718);
	my $x;
	#y = $K/(1+$b*$e**(-1*$a*$x));

	my $s;
	if($v<$minb){
		#up curve1
		$x = ($v-$minl)/(($minb-$minl)/10)-5;
		my $y1 = $K/(1+$b*$e**(-1*$a*$x));
		$x = ($min-$minl)/(($minb-$minl)/10)-5; 
		my $y10 = $K/(1+$b*$e**(-1*$a*$x)); ## v=min, score=0
		$s = $score * ($y1-$y10)/(1-$y10);
	}elsif($v<=$maxb){
		$s = $score;
	}else{
		if($maxl-$maxb==0){
			print join("\t", ($v, $score, $minb, $maxb, $min, $max, $minl, $maxl)),"\n";
			die;
		}
		#down curve2
		$x = ($maxb-$v)/(($maxl-$maxb)/10)+5;
		my $y2 = $K/(1+$b*$e**(-1*$a*$x));
		$x = ($maxb-$max)/(($maxl-$maxb)/10)+5; ## v=max, score=0
		my $y20 = $K/(1+$b*$e**(-1*$a*$x));
		$s = $score * ($y2-$y20)/(1-$y20);
	}
	return $s;
}



## 打分
## score:满分
sub score_single{
	my ($v, $score, $minb, $maxb, $min, $max)=@_;
	if(!defined $v){
		print STDERR join("\t", $v, $score, $minb, $maxb, $min, $max),"\n";
	}
	return $score if($v eq "NULL");
	if(!defined $minb){
		print STDERR join("\t", $v, $score, $minb, $maxb, $min, $max),"\n";
	}
	my $disl = $minb - $min;
	my $disr = $max - $maxb;
	if($disl<0 || $disr<0){
		print "Wrong Score set: $minb, $maxb, $min, $max\n";
		die;
	}

	my $s;
	if($v<$minb){
		$disl=$disl==0? 0.1: $disl;
		$s = $score * (1 - ($minb-$v)/$disl);
	}elsif($v<=$maxb){
		$s = $score;
	}else{
		$disr= $disr==0? 0.1: $disr;
		$s = $score * (1 - ($v-$maxb)/$disr);
	}
	return $s;
}



1;
