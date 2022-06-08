my $score_des_probe= "scores of Tm, GC, self-complementary, CG content diff, dupGNum, 5endG, SNPs, polys, bounds.";
my $score_des_probe_meth= "scores of Tm, GC, self-complementary, CG content diff, dupGNum, 5endG, SNPs, polys, bounds, CpGs, Cs.";
my $score_des_pair_primer= "scores of primer1 3'end to target(Probe 5'end), distance between pair, length diff, Tm diff, specificity of primer pair.";
my $score_des_pair_probe= "score of probe bounds.";
my $score_des_primer= "scores of Tm, GC, self-complementary, end3 A num, end3 stability, SNPs, polys, bounds.";
my $score_des_primer_meth= "scores of Tm, GC, self-complementary, end3 A num, end3 stability, SNPs, polys, bounds, CpGs, Cs.";

my @feature = ("Tm", "GC", "Hairpin", "DimerType", "DimerSize", "EndANum", "EndStability", "SNPs", "Polys");
my @meth = ("CpGs", "Cs");
my @bound = ("BoundsNum", "BoundsTm", "BoundsInfo");
my @bound_final = ("OligoBounds", "ProdsNum", "ProdsEff", "ProdsInfo");

my @feature_score_primer = ("Tm", "GC", "Hairpin", "EndANum", "EndStability", "SNPs", "Polys");
my @feature_score_probe  = ("Tm", "GC", "Hairpin", "CGdiff", "5endG", "SNPs", "Polys");
my @meth_score = @meth;
my @bound_score = ("Bounds");

my $base_head = "#ID\tSeq\tLen";
my $base_head_final = "#Chr\tStart\tStrand\tID\tSeq\tLen";
my $score_head = "Score\tScoreRel\tScoreInfo";
my $score_head_final = "Score\tScoreRel\tScoreOligo";

sub score_des{
	my ($type, $Methylation)=@_;
	my $des;
	if($type eq "Probe"){
		if(defined $Methylation){
			$des=$score_des_probe_meth;
		}else{
			$des=$score_des_probe;
		}
	}elsif($type eq "Primer"){
		if(defined $Methylation){
			$des=$score_des_primer_meth;
		}else{
			$des=$score_des_primer;
		}
	}elsif($type eq "PairPrimer"){
		$des=$score_des_pair_primer;
	}elsif($type eq "PairProbe"){
		$des=$score_des_pair_probe;
	}
	return $des;
}

sub read_evaluation_info{
	my ($idx, $line, $Methylation, $NoSpecificity)=@_;
	my @unit = split /\t/, $_;
	my @base=@unit[0..($idx+2)]; ## ID\tSeq\tLen
	my @feature=@unit[($idx+3)..($idx+11)];
	my $ix=12+$idx;
	my (@meth, @spec);
	if(defined $Methylation){
		@meth=@unit[($idx+12)..($idx+13)];
		$ix=14+$idx;
	}
	my $bnumtm;
	if(!defined $NoSpecificity){
		@spec=@unit[$ix..($ix+2)];
		$bnumtm=$unit[$ix]."|".$unit[$ix+1];
	}
	return (\@base, \@feature, \@meth, \@spec, $bnumtm);
}


sub final_evalue_head{
	my ($Methylation, $NoSpecificity)=@_;
	my @title = ($base_head);
	push @title, $score_head;
	push @title, &features("FinalHead", $Methylation, $NoSpecificity);
	return join("\t", @title);
}



sub final_head{
	my ($Methylation, $NoSpecificity, $range_pos, $ptype)=@_;
	my @title = ($base_head_final);
	if(defined $range_pos){
		push @title, "Dis2Target";
	}
	if($ptype=~/face-to-face/){
		push @title, "ProductSize";
	}else{
		push @title, "DisBetweenPairs";
	}
	push @title, $score_head_final;
	push @title, &features("FinalHead", $Methylation, $NoSpecificity);
	return join("\t", @title);
}

sub score_head{
	my ($Methylation, $NoSpecificity)=@_;
	my $head=join("\t", $base_head, $score_head, &features("Head", $Methylation, $NoSpecificity));
	return $head;
}

sub evaluation_head{
	my ($Methylation, $NoSpecificity)=@_;
	my $head=join("\t", $base_head, &features("Head", $Methylation, $NoSpecificity));
	return $head;
}

sub features{
	my ($type, $Methylation, $NoSpecificity)=@_;
	my @features;
	my (@fea, @met, @bou);
	if($type eq "Head"){
		@fea = @feature;
		@met = @meth;
		@bou = @bound;
	}elsif($type eq "FinalHead"){
		@fea = @feature;
		@met = @meth;
		@bou = @bound_final;
	}elsif($type eq "ScorePrimer"){
		@fea = @feature_score_primer;
		@met = @meth_score;
		@bou = @bound_score;
	}elsif($type eq "ScoreProbe"){
		@fea = @feature_score_probe;
		@met = @meth_score;
		@bou = @bound_score;
	}
	push @features, @fea;
	if(defined $Methylation){
		push @features, @met;
	}
	if(!defined $NoSpecificity){
		push @features, @bou;
	}
	return @features;
}

