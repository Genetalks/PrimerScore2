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
			my @bnum = (1,5,1,100); #bound number
			@bvalue = (0,$value[0]*0.65,0,$value[0]); #bound sec value
			$evalue = &score_single($value[1], 0.7, @bvalue);
			$enum = &score_single($bnum, 0.3, @bnum);
			$etotal = $evalue+$enum;
		}elsif($type eq "Probe_Tm"){
			my @bnum = (1,5,1,100); #bound number
			@bvalue = (0,$value[0]*0.65,0,$value[0]); #bound sec value
			$evalue = &score_single($value[1], 0.9, @bvalue);
			$enum = &score_single($bnum, 0.1, @bnum);
			$etotal = $evalue+$enum;
		}else{ #Eff
			@bvalue = (0,0.001,0,0.1);
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
sub poly_score{
	my ($info, $len, $type)=@_;
	return 1 if($info eq "NA");
	my @polys = split /,/, $info;
	my $score=1;
	my @dprobe=($len*0.3, $len*0.5, 0, $len*0.5);
	my @dprimer=(8,$len,0,$len);
	my @polyleno=(0,2.5,0,8);
	my @polylenG=(0,2.5,0,5);
	for(my $i=0; $i<@polys; $i++){
		my ($d3, $t, $l)=$polys[$i]=~/(\d+)(\w)(\d+)/;
		my $s;
		my @polylen=@polyleno;
		if($t eq "G"){
			@polylen=@polylenG;
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

1;
