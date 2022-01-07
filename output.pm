my $score_des_probe= "##ScoreInfo: scores of length, tm, self-complementary, CG content diff, 5endG, snp, poly, mapping\n";
my $score_des_probe_meth= "##ScoreInfo: scores of length, tm, self-complementary, CG content diff, 5endG, snp, poly, mapping, CpGs, Cs\n";

my $base_head="#ID\tSeq\tLen";
my $feature_head="Tm\tGC\tHairpin\tDimerType\tDimerSize\tEndANum\tEndStability\tSNP\tPoly";
my $bound_head="MapNum\tMapTm\tMapInfo";
my $meth_head="CpGNum\tCNum";
my $score_head="Score\tScoreInfo";

sub score_des{
	my ($type, $Methylation)=@_;
	my $des;
	if($type eq "Probe"){
		if(defined $Methylation){
			$des=$score_des_probe_meth;
		}else{
			$des=$score_des_probe;
		}
	}else{

	}
	return $des;
}

sub score_head{
	my ($Methylation, $NoSpecificity)=@_;
	my $head=join("\t", $base_head, $score_head, $feature_head);
	if(defined $Methylation){
		$head.="\t".$meth_head;
	}
	if(!defined $NoSpecificity){
		$head.="\t".$bound_head;
	}
	return $head;
}

sub evaluation_head{
	my ($Methylation, $NoSpecificity)=@_;
	my $head=$base_head."\t".$feature_head;;
	if(defined $Methylation){
		$head.="\t".$meth_head;
	}
	if(!defined $NoSpecificity){
		$head.="\t".$bound_head;
	}
	return $head;
}

