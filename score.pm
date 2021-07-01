sub bound_score{
	my ($bnum, $bvalue, $fulls, $type)=@_;
	#bvalue: tm, Eff
	my $sbound;
	if($bnum==1){
		$sbound=$fulls;
	}else{
		my @value=split /,/, $bvalue;
		if(scalar @value==1){
			return $fulls;
		}
		my (@bvalue, $evalue, $enum, $etotal);
		if($type eq "Tm"){
			my @bnum = (1,5,1,100); #bound number
			@bvalue = (0,$value[0]*0.65,0,$value[0]); #bound sec value
			$evalue = &score_single($value[1], 0.7, @bvalue);
			$enum = &score_single($bnum, 0.3, @bnum);
			$etotal = $evalue+$enum;
		}else{ #Eff
			#my @bnum = (1,5,1,100); #bound number
			#my @bnum = (1,1,1,5); #bound number
			@bvalue = (0,0.001,0,0.1);
			#@bvalue = (0,0.001,0,0.1);
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
	my @dprobe=($len*0.4, $len*0.5, 0, $len*0.5);
	my @dprimer=(15,$len,0,$len);
	my @polylen=(0,2.8,0,8);
	for(my $i=0; $i<@polys; $i++){
		my ($d3, $t, $l)=$polys[$i]=~/(\d+)(\w)(\d+)/;
		my $s;
		if($type eq "Probe"){
			my $d5=$len-$d3-1;
			my $de=$d5>$d3?$d3:$d5;
			$s = &score_single($de, 1, @dprobe)*&score_single($l,1,@polylen);
		}else{
			$s = &score_single($d3, 1, @dprimer)*&score_single($l,1,@polylen);
		}
		$score*=$s;
	}
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

1;
