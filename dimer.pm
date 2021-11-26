
require "$Bin/score.pm";
my $min_tm_omega = 17;
my $min_tm_amp = 29; ##Dis-7:26
my $min_meetlen=3; ## min end3 match len when Enddimer
my $min_amplen=15;


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
### return: 
## pright: primer end3 length of non-mapped
## tleft:  template end3 length of non-mapped
## amplen: primer amplified length, ie. template end5 length of non-mapped
## mlen3:  primer end3 mapped length
## len:    product length
## msum:   mapped crowd sum
sub dimer_amplify{
	my ($info)=@_;
	if(!defined $info || $info eq ""){
		return(-100, -1, -1,-1, -1, -1, -1);
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
	$line[2]=~s/STR\t//;
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

	## amplify length
	my $amplen=$tright;

	## end3 map length
	$line[1]=~s/^\s+//;
	my @maps=split /\s+/, $line[1];
	my $msum=scalar @maps;
	my $mlen3=length($maps[-1]);
	my @unmaps=split /\s+/, $line[0];
	if($mlen3<=2){## not really map
		my $unmap;
		if($amplen>0){## $unmaps[-1] is "---"
			$unmap=$unmaps[-2];
			$unmap=~s/\-//g;
			$pright=$mlen3+length($unmap);
		}else{
			$unmap=$unmaps[-1];
			$unmap=~s/\-//g;
			$pright=$mlen3+length($unmap);
		}
		if(defined $maps[-2]){
			$mlen3=length($maps[-2]);
		}

	}
	
	## indel
	my $indel=0;
	pop @unmaps;
	my @tunmaps=split /\s+/, $line[3];
	my @unmapall=(@unmaps, @tunmaps);
	for(my $i=0; $i<@unmapall; $i++){
		if($unmapall[$i]=~/\-\-/){
			$indel++;
		}
	}
	## 
	my $len = length $line[0];
	return ($tm, $pright, $tleft, $amplen, $mlen3, $len, $msum, $indel);
}

sub judge_amplify_endmeet{
	my ($tm, $end31, $end32, $mlen3)=@_;
	my $is_amplify=0;
	my $type = "NoAmplify";
	my $eff=0;
	my @tm;
	if($end31==0 && $end32==0 && $mlen3>=$min_meetlen){
		$is_amplify=1;
		@tm=(10,100,-50,100,-50,100); #$minb, $maxb, $min, $max, $minl, $maxl
		$eff=&score_growth_curve($tm, 1, @tm);
		$type = "AmpEndMeet";
	}
	return ($is_amplify, $type, $eff);
}

sub judge_amplify{
	my ($tm, $end31, $end32, $amplen, $mlen3, $msum, $indel)=@_;
	
	my $type="NoAmplify";
	my $eff=0;
	my @tm;
	### judge is amplify or not
	if($end31==0 && $end32==0){
		if($msum==1){
			@tm=(10,100,-50,100,-50,100); #$minb, $maxb, $min, $max, $minl, $maxl
			$eff=&score_growth_curve($tm, 1, @tm);
			$type = "AmpEndMeet";
		}else{
			if($indel==0){
				@tm=(20,100,0,100,0,100);
			}else{
				@tm=(30,100,10,100,10,100);
			}
			$eff=&score_growth_curve($tm, 1, @tm);
			$type = "AmpEndMap";
		}
	}elsif($end31<=1 || $end32<=1){
		if($amplen>=$min_amplen){
			if($tm>=$min_tm_omega){
				if($indel==0){
					@tm=(30,100,15,100,15,100); #$minb, $maxb, $min, $max, $minl, $maxl
				}else{
					@tm=(30,100,15,100,15,100); #$minb, $maxb, $min, $max, $minl, $maxl
				}
				my @alen=(30,100,$min_tm_omega-5,100,$min_tm_omega-5,100);
				$eff=&score_growth_curve($tm, 1, @tm) * &score_growth_curve($amplen, 1, @alen);
				$type = "AmpOmega";
			}
		}else{
			if($tm>=$min_tm_amp){
				if($indel==0){
					@tm=(50,100,25,100,25,100);
				}else{
					@tm=(60,100,40,100,40,100);
				}
				$eff=&score_growth_curve($tm, 1, @tm);
				$type = "AmpEachOther";
			}
		}
	}
	my $is_amplify=0;
	if($type ne "NoAmplify"){
		$is_amplify=1;
	}
	return($is_amplify, $type, $eff);
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
