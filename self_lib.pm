sub get_highest_bound{
	my ($abound, $maxn)=@_;
	my %bound=%{$abound};
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
	my $binfos=join(";", @binfos);
	my $bvalues=join(",", @bvalue);
	return ($bnum, $bvalues, $binfos);
}

#my @rank_end=(3,   5,   8); #  PCR efficiency when dis of mismatch pos to 3end <= @rank_end
#my @eff_end =(0.1, 0.4, 0.8 );
#my $max_end3_eff = $rank_end[-1];
#my $eff_times = 10;
#my $min_eff = 0.1;
sub efficienty_end3_mismatch{
	my ($mvisual,$max_end3, $arank_end, $aeff_end)=@_;
	my $mapr= reverse ($mvisual);
	my @unit = split //, $mapr;
	my @mis_pos;
	for(my $i=0; $i<$max_end3; $i++){
		if($unit[$i] eq "*"){# mismatch
			push @mis_pos, $i;
		}elsif($unit[$i] eq "-" || $unit[$i] eq "^"){# Del or Ins ==> 2 mismatch
			push @mis_pos, ($i, $i);
			for(my $j=$i+1; $j<$max_end3;$j++){
				if($unit[$j] ne $unit[$i]){
					$i=$j-1;
					last;
				}
			}
		}
	}
	my $eff_end = 1;
	for(my $i=0; $i<@mis_pos; $i++){
		$eff_end *= &get_eff_rank($mis_pos[$i], $arank_end, $aeff_end, "<=");
	}
	return $eff_end;
}

sub get_eff_rank{
	my ($v, $arank, $aeff, $compare)=@_;
	my @rank = @{$arank};
	my @eff = @{$aeff};
	
	my $eff;
	for(my $i=0; $i<@rank; $i++){
		if($compare eq ">="){
			if($v >= $rank[$i]){
				$eff = $eff[$i];
				last;
			}
		}else{
			if($v <= $rank[$i]){
				$eff = $eff[$i];
				last;
			}
		}
	}
	if(!defined $eff){
		$eff = 0;
	}
	return $eff;

}


sub poly_check{
	my ($seq)=@_;
	$seq = reverse $seq;
	my @unit = split //, $seq;
	
	my $min_plen = 3; # min poly length
	my @polys;
	my $last = $unit[0];
	my $poly = $unit[0];
	my $maxl=0; ## longest poly len
	my $maxb; ## longest poly base
	my $total=0; ## total poly len
	for(my $i=1; $i<@unit; $i++){
		if($unit[$i] eq $last){
			$poly.=$unit[$i];
		}else{
			my $len = length($poly);
			if($len>=$min_plen){
				$total+=$len;
				if($len>$maxl){
					$maxl=$len;
					$maxb=$last;
				}
				push @polys, ($i-$len).$last.$len;
			}
			$last = $unit[$i];
			$poly=$last;
		}
	}
	my $len = length $poly;
	if($len>=$min_plen){
		$total+=$len;
		if($len>$maxl){
			$maxl=$len;
			$maxb=$last;
		}
		push @polys, ($#unit+1-$len).$last.$len;
	}
	return ($total, $maxl, $maxb, join(",", @polys));
}

#calculate a poly factor according to poly len、num and position
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

sub get_end3_detaG{
	my ($primer, $endlen)=@_;
	my $len = length $primer;
	my $endseq=substr($primer, $len-$endlen);
	my %ix=("A",0, "C",1, "G",2, "T",3);
	my @matrix=([-1.9,-1.3,-1.6, -1.5],[-1.9,-3.1,-3.6,-1.6],[-1.6,-3.1,-3.1,-1.3],[-1.0,-1.6,-1.9,-1.9]);
	my @b = split //, $endseq;
	my $detaG=0;
	for(my $i=0; $i<@b-1; $i++){
		my $ix1=$ix{$b[$i]};
		my $ix2=$ix{$b[$i+1]};
		$detaG+=$matrix[$ix1][$ix2];
	}
	return $detaG;
}

## get endA number
sub get_end_A{
	my ($seq)=@_;
	my @unit = split //, $seq;
	my $Anum = 0;
	my $rnum = 0;
	my %class=(A=>"AT",T=>"AT",C=>"CG",G=>"CG",a=>"AT",t=>"AT",c=>"CG",g=>"CG");

	my $endc = $class{$unit[-1]};;
	for(my $i=@unit; $i>0; $i--){
		last if($class{$unit[$i-1]} ne $endc);
		$rnum++;
	}
	$Anum = $endc eq "AT"? $rnum: 0;
	return ($Anum);
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
## eg:(off=6)  
##:  #|||**-||||^^^|||||||||||||             =>      ||||^^^||||||||||||| 
##:  ||*||||^^|||----------||||||||||||      =>      ####||||||||||||
##:  #######||*|||||||||||  (off=4)          =>      ###||*|||||||||||
sub map_visual_trim{
	my ($mv, $off, $end3_match)=@_;
	if($off ==0){
		return($mv, $end3_match);
	}
	#trim from offset
	my @unit=split //, $mv;
	my $ntrim=0;
	while($ntrim<$off){
		if($unit[0] ne "-" && $unit[0] ne "^"){
			shift @unit;
			$ntrim++;
		}else{
			shift @unit;
		}
	}

	#rm indel 
	while($unit[0] eq "-" || $unit[0] eq "^"){
		shift @unit;
	}
	my $mvsub = join("", @unit);
	
	##
	if($unit[0] eq "#"){ ## no need adjust
		return($mvsub, $end3_match);
	}

	## adjust
	my $min_score=8;
	my ($wmatch, $wmis, $windel)=(2,1,2);
	my $irm=-1;
	my $score=0;
	for(my $i=0; $i<@unit; $i++){
		if($unit[$i] eq "|"){
			$score+=$wmatch;
		}elsif($unit[$i] eq "*"){
			$score-=$wmis;
		}elsif($unit[$i] eq "-" || $unit[$i] eq "^"){## InDel
			my $len=1;
			for(my $j=$i+1; $j<@unit; $j++){
				if($unit[$j] ne $unit[$i]){
					$score-=$windel * (int($len/5)+1);
					$i=$j-1;
					last;
				}
				$len++;
			}
		}
		if($score<=0){
			$irm=$i;
			$score=0;
		}elsif($score>=$min_score){
			last;
		}
	}
	my $mvaf = substr($mvsub, $irm+1);
	my $mvbf ="";
	for(my $i=0; $i<=$irm;$i++){
		if($unit[$i] eq "*" || $unit[$i] eq "|"){
			$mvbf .= "#";
		}
	}
	my $mvn = $mvbf.$mvaf;
	if($irm==-1 && $end3_match > length($mvsub)){
		$end3_match=length($mvsub);
	}
	return ($mvn, $end3_match);
}



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

sub map_visual_from_ntthal{
	my ($aline, $is_reverse, $start, $end)=@_;
	my @line = @{$aline};
	$line[0]=~s/SEQ\t//;
	$line[3]=~s/STR\t//;
	my @pmap = split //, $line[0];
	my @tmap = split //, $line[3];
	my ($tleft, $tright, $pleft, $pright)=(0,0,0,0);
	## trim end 
	while($pmap[0] eq " " && $tmap[0] ne " "){
		shift @pmap;
		shift @tmap;
		$tleft++;
	}
	while($pmap[-1] eq "-"){
		pop @pmap;
		pop @tmap;
		$tright++;
	}

	## get map visual
	my @mv;
	for(my $i=0; $i<@pmap; $i++){
		if($pmap[$i] eq " "){
			$mv[$i]="|";
		}elsif($pmap[$i] eq "-"){
			$mv[$i]="-";
		}else{
			if($tmap[$i] eq "-"){
				$mv[$i]="^";
			}else{
				$mv[$i]="*";
			}
		}
	}
	## * of end => #
	for(my $i=0; $i<@mv;$i++){
		last if($mv[$i] ne "*");
		$mv[$i]="#";
		$pleft++;
		$tleft++;
	}
	for(my $i=$#mv; $i>=0;$i--){
		last if($mv[$i] ne "*");
		$mv[$i]="#";
		$pright++;
		$tright++;
	}

	my $mv = join("", @mv);
	my ($pos3, $pos5);
	if($is_reverse){
		$pos5=$end-$tleft;
		$pos3=$start+$tright;
	}else{
		$pos5=$start+$tleft;
		$pos3=$end-$tright;
	}
	return ($mv, $pos3, $pos5);
}


sub end_match_length{
	my ($mv)=@_;
	my @mv=split //, $mv;
	my $len=1;
	my $ref=$mv[$#mv];
	for($i=$#mv-1; $i>=0; $i--){
		last if($mv[$i] ne $ref);
		$len++;
	}
	my $endm;
	if($ref eq "|"){
		$endm=$len;
	}elsif($ref eq "#"){
		$endm=$len*(-1);
	}else{ ## primer end is not covered by template
		print "Warn: end3 map visual info $mv\n";
		$endm="Fail";
	}
	return $endm;
}



#1S5M1D4M3I13M   3G0T0^G17       =>      #|||**-||||^^^|||||||||||||
#13M16D5M        13^GGGCGTCAGATGCAGG5    =>      |||||||||||||----------------|||||
#10M1D4M1I3M     10^T3C3 =>      ||||||||||-|||*^|||
sub map_visualation{
	my ($cigar, $md)=@_;
	my ($acigar_n, $acigar_str)=&cigar_split($cigar);
	my @numcg=@{$acigar_n};
	my @strcg=@{$acigar_str};

	my ($amds)=&md_split($md);
	my @mds=@{$amds};
	my $minfo;
	my $imd=0; ## md index
	for(my $i=0; $i<@strcg; $i++){
		if($strcg[$i] eq "S" || $strcg[$i] eq "H"){
			$minfo.="#" x $numcg[$i];
		}elsif($strcg[$i] eq "M"){
			my $Mlen=0;
			for(my $j=$imd; $j<@mds; $j++){
				#print join("\t", $i, $j, $mds[$j], $Mlen, $numcg[$i]),"\n";
				if($mds[$j]=~/[0-9]/){
					if($Mlen+$mds[$j]<=$numcg[$i]){
						$minfo.="|" x $mds[$j];
						$Mlen+=$mds[$j];
					}else{
						my $l=$numcg[$i]-$Mlen;
						$mds[$j]-=$l;
						$minfo.="|" x $l;
						$imd=$j;
						last;
					}
				}elsif($mds[$j]=~/\^/){
					die "Wrong: $md, $cigar\n";
				}elsif($mds[$j]=~/[ATCG]/){
					if(length $mds[$j] > 1){
						die "Wrong: $md\n";
					}
					$minfo.="*";
					$Mlen++;
				}

				if($Mlen==$numcg[$i]){
					$imd=$j+1;
					if($imd<$#mds && $mds[$imd] eq '0'){
						$imd++;
					}
					last;
				}

			}
		}elsif($strcg[$i] eq "I"){
			#$minfo.="-" x $numcg[$i];
			$minfo.="^" x $numcg[$i];
		}elsif($strcg[$i] eq "D"){
			#$minfo.="^" x $numcg[$i];
			$minfo.="-" x $numcg[$i];
			$imd++;
		}else{
			die "Wrong: $cigar\n";
		}
	}
	#print join("\t", $cigar, $md, "=>", $minfo),"\n";
	return $minfo;
}

#Example: 
#cigar: 2S21M1D3M
#md: 13C5T0C^A3
#SEQ                               TT             A     AT-   ----
#SEQ                                 TGGGTGGTGCTAC TCTTC   AAT
#STR                                 ACCCACCACGATG AGAAG   TTA
#STR     AACCGTCCCAACCCCCAACACACCCCCC             G     AGT   AGGA
sub get_match_cigar_from_ntthal{
	my ($aline, $is_reverse, $start, $end)=@_;
	my @line = @{$aline};
	$line[0]=~s/SEQ\t//;
	$line[1]=~s/SEQ\t//;
	$line[2]=~s/STR\t//;
	$line[3]=~s/STR\t//;
	my @punmap = split //, $line[0];
	my @pmap = split //, $line[1];
	my @tmap = split //, $line[2];
	my @tunmap = split //, $line[3];
	my $is_start = 1;
	my ($tleft, $tright, $pleft, $pright)=(0,0,0,0);
	my ($cigar, $md);
	my ($mlen, $Mlen_cigar)=(0, 0);
	for(my $i=0; $i<@punmap; $i++){
		if($pmap[$i] ne " " && $tmap[$i] ne " "){ ### match
			$mlen++;
			$Mlen_cigar++;
			$is_start=0;
		}else{### not match
			if($is_start){ ## start: soft
				if($punmap[$i] ne " "){
					$pleft++;
				}
				if($tunmap[$i] ne " "){
					$tleft++;
				}
				next;
			}
				
			if($punmap[$i] ne "-" && $tunmap[$i] ne "-"){ ## mismatch
				$Mlen_cigar++;
				my $b = $tunmap[$i];
				$b=~tr/ATCG/TAGC/;
				$md .= $mlen.$b;
				$mlen=0;
			}else{## Indel
				$cigar.=$Mlen_cigar."M";
				$Mlen_cigar=0;

				my $len = 0;
				my $str;
				if($punmap[$i] eq "-"){ # del
					$md.=$mlen;
					$mlen=0;
					$str=$tunmap[$i];
					$len++;
					my $j=$i+1;
					while($j<@punmap){
						last if($punmap[$j] ne "-");
						$len++;
						$str.=$tunmap[$j];
						$j++;
					}
					$str=~tr/ATCG/TAGC/;
					$md.="^".$str;
					$cigar.=$len."D";
					$i=$j-1;
				}else{ # insert
					$len++;
					my $j=$i+1;
					while($j<@punmap){
						last if($tunmap[$j] ne "-");
						$len++;
						$str.=$punmap[$j];
						$j++;
					}
					$cigar.=$len."I";
					$i=$j-1;
				}
			}
		}
	}
	if($mlen>0){
		$cigar.=$Mlen_cigar."M";
		$md.=$mlen;
	}

	## end handle D
	if($cigar=~/\d+D$/){
		($tright)=$cigar=~/(\d+)D$/;
		$cigar=~s/\d+D$//;
		$md=~s/\^[A-Z]+$//;
	}elsif($cigar=~/\d+I$/){
#eg: 10M12D8M2I
#
#                     ------------     TGACC
#           GCCACTGCCT            GCTGG
#           CGGTGACGGA            CGACC
#CGTACACCCAC          GAAACACTGTAT     GAC--
		$cigar=~s/\d+I//;
	}
	## end mismatch to Soft
	if($md=~/[A-Z]$/){ ##end mismatch
		$md.="0"; 
	}
	my $end_match;
	my $mislen=0;
	if($md=~/[ATCG]0$/){
		while($md=~/[ATCG]0$/){
			$md=~s/[ATCG]0$//;
			$mislen++;
		}
		my ($mlen)=$cigar=~/(\d+)M$/;
		$mlen-=$mislen;
		$cigar=~s/(\d+)M$//;
		$cigar.=$mlen."M".$mislen."S";
		$end_match = -1*$mislen;
	}else{
		($end_match)=$md=~/(\d+)$/;
	}
	if($pleft>0){
		$cigar = $pleft."S".$cigar;
	}

	my $pos_new = $is_reverse? $start+$tright+$mislen: $end-$tright-$mislen; ## position of  primer3
	#my $pos_new = $is_reverse? $start+$tright+$mislen: $start+$tleft;
	return($pos_new, $cigar, $md, $end_match);
}


