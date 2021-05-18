
sub poly_score{
	my ($info, $len, $type)=@_;
	return 1 if($info eq "NA");
	my @polys = split /,/, $info;
	my $score=1;
	my @dprobe=($len*0.4, $len*0.5, 0, $len*0.5);
	my @dprimer=(10,$len,0,$len);
	my @polylen=(0,2,0,10);
	for(my $i=0; $i<@polys; $i++){
		my ($d3, $t, $l)=$polys[$i]=~/(\d+)(\w)(\d+)/;
		my $s;
		if($type eq "Probe"){
			my $d5=$len-$d3-1;
			my $de=$d5>$d3?$d3:$d5;
			$s = 0.2+&score_single($de, 0.2, @dprobe)+&score_single($l,0.6,@polylen);
		}else{
			$s = 0.1+&score_single($d3, 0.3, @dprime)+&score_single($l,0.6,@polylen);
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
	my @dprimer=(10,$len,0,$len);
	my @snplen=(1,1,1,10); #snp length: SNP-1, Indel-indel len
	for(my $i=0; $i<@snps; $i++){
		my ($d3, $t, $l)=$snps[$i]=~/(\d+)(\w)(\d+)/;
		my $s;
		if($type eq "Probe"){
			my $d5=$len-$d3-1;
			my $de=$d5>$d3?$d3:$d5;
			$s = 0.4+&score_single($de,0.4,@dprobe)+&score_single($l,0.1,@snplen);
		}else{
			$s = 0.1+&score_single($d3, 0.7, @dprime)+&score_single($l,0.1,@snplen);
		}
		$score*=$s;
	}
	return $score;
}

1;
