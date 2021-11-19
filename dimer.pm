#Example: 
#mv:                               ##|||||||||||||*|||||**-|||
#SEQ                               TT             A     AT-   ----
#SEQ                                 TGGGTGGTGCTAC TCTTC   AAT
#STR                                 ACCCACCACGATG AGAAG   TTA
#STR     AACCGTCCCAACCCCCAACACACCCCCC             G     AGT   AGGA

#mv:                             ###|||*|*|||**^^^^||||||*||||||||||##
#SEQ                             ATA   A G   GGTTTC      C          AC---
#SEQ                                AAT A CCT      AACCCT TAATTCTTGC
#STR                                TTA T GGA      TTGGGA ATTAAGAACG
#STR     AACTGAAATCTTCTTCTCTCCGACAGG   G G   AA----      C          ATTCC

sub dimer_amplify{
	my ($info)=@_;
	if(!defined $info || $info eq ""){
		return(-100, -1, -1);
	}
	my ($head,@line) = split /\n/, $info;
	if(!defined $head){
		print $info,"\n";
		die;
	}
	my ($tm)=$head=~/t = (\S+)/;
	if(!defined $tm){
		print "yes:\n";
		print $info,"\n";
		print $head,"\n";
		die;
	}
	$line[0]=~s/SEQ\t//;
	$line[1]=~s/SEQ\t//;
	$line[3]=~s/STR\t//;
	my @punmap = split //, $line[0];
	my @pmap = split //, $line[1];
	my @tunmap = split //, $line[3];
	my ($tleft, $tright, $pleft, $pright)=(0,0,0,0);
	## trim end 
	for(my $i=0; $i<@punmap; $i++){
		last if($pmap[$i] ne " ");
		if($punmap[$i] ne " "){
			$pleft++;
		}
		if($tunmap[$i] ne " "){
			$tleft++;
		}
	}
	for(my $i=$#punmap; $i>0; $i--){
		last if($pmap[$i] ne " ");
		if($punmap[$i] ne "-"){
			$pright++;
		}
		if($tunmap[$i] ne "-"){
			$tright++;
		}
	}

	## end3 map length
	my @maps=split /\s+/, $line[1];
	my $mlen3=length($maps[-1]);
	return ($tm, $pright, $tleft, $mlen3);
}

#sub efficiency_dimer{
#	my ($tm, $pright, $tleft)=@_;
#	my @etm=(0, 100, -30, 100);
#	my $eff_tm = &score_single($tm, 1, @etm);
#	$eff_tm = $eff_tm>0? $eff_tm: 0;
#
#	my @
#
#}

1;
