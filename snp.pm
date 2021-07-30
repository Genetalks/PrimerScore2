
our @Ins=("I", "J", "L", "L");
my @dg  =("R",   "Y",   "M",   "K",   "S",   "W",   "H",     "B",     "V",     "D",     "N");
my @atcg=("A/G", "C/T", "A/C", "G/T", "C/G", "A/T", "A/C/T", "C/G/T", "A/C/G", "A/G/T", "A/C/G/T");
sub snp_to_degenerate{
	my ($info)=@_;
	$info=~s/,$//;
	if($info=~/N/){
		return "N";
	}
	my @snp=split /,/, $info;
	for(my $i=0; $i<@snp; $i++){
		$snp[$i]=uc($snp[$i]);
		if($snp[$i]!~/[ATCGatcg]/){
			die "Wrong snps:", $info,"\n";
		}
	}
	my @snps=sort{ord($a) <=>ord($b)} @snp;
	my %hash;
	for(my $i=0; $i<@dg; $i++){
		$hash{$atcg[$i]}=$dg[$i];
	}
	return $hash{join("/", @snps)};
}
sub degenerate_to_snp{
	my ($v)=@_;
	my %hash;
	for(my $i=0; $i<@dg; $i++){
		$hash{$dg[$i]}=$atcg[$i];
	}
	return $hash{$v};
}

sub SNP_parse{
	my ($seq)=@_;
	$seq = reverse $seq;
	my @u=split //, $seq;
	my @snp;
	my %snp;
	for(my $i=0; $i<@dg; $i++){
		$snp{$dg[$i]}=1;
	}

	for(my $i=0; $i<@u; $i++){
		if(exists $snp{$u[$i]}){
			push @snp, $i."S1";
		}elsif($u[$i] eq "E"){
			my $l=1;
			while(1){
				if($u[$i+1] eq "E"){
					$l++;
					$i++;
				}else{
					last;
				}
			}
			push @snp, $i."D".$l;
		}elsif($u[$i]=~/IJL/){
			my $l=ord($u[$i])-ord('I')+1;
			push @snp, $i."I".$l;
		}
	}
	return (scalar @snp, join(",", @snp));
}



### add snp/indel info to sequence, convert as following:
#G       A         ==> R(degenerate)
#G       A,T       ==> D(degenerate)
#AAAG    A         ==> AEEE(D is degenerate)
#G       GA        ==> I
#G       GAA       ==> J
#G       GAAA      ==> L(K is degenerate)
#G       GATTTT    ==> L
#G       GA,GAAA   ==> L
## aseq: addr of sequence array
## p: polym(snp or indel) position, count from 1
## ref: ref base, ref in vcf
## alt: alt base, alt in vcf
sub sequence_convert_snp{
	my ($aseq, $p, $ref, $alt)=@_;
	my $seqlen = scalar(@{$aseq});
	
	my $type=&mutant_type($ref, $alt);
	my @types = split /,/, $type;
	my ($polym, $len, $off)=("S", 0, 0);
	my $num=0;
	#choose the longest mutant
	foreach my $tp(@types){##usually D will not be with S and I.
		$num++;
		my ($t, $l, $o)=split /_/, $tp;
		if($l>$len){
			$polym=$t;
			$len=$l;
			$off=$o;
		}
	}
	if($p-1>=0 && $p-1<=$seqlen-1){ ## the pos is on the sequence 
		if($polym eq "S"){
			$alt=~s/<[xX]>//;
			my $dg=&snp_to_degenerate($ref.",".$alt);
			$aseq->[$p-1]=$dg;
		}elsif($polym eq "I"){
			my $ix=($len-1)>=3? 3: ($len-1);
			$aseq->[$p-1]=$Ins[$ix];
		}
	}
	if($polym eq "D"){
		for(my $i=$p+$off; $i<($p+$off+$len); $i++){
			if($i-1>=0 && $i-1<=$seqlen-1){
				$aseq->[$i-1]="E";
			}
		}
	}
}


sub mutant_type{
	my ($ref, $alt)=@_;
	my $lenr=$ref eq "-"? 0: length $ref;
	my @alts=split /,/, $alt;
	my $lena=0;
	my @type;
	foreach my $at(@alts){
		my $type;
		my $lena = $at eq "-"? 0: length $at;
		if($lenr == $lena){ ## SNP
			$type="S_0_0";
		}else{
			my $dl=abs($lenr-$lena);
			if($lenr > $lena){## Del
				$type="D_".$dl."_".$lena;
			}else{ ## Ins
				$type="I_".$dl."_".$lenr;
			}
		}
		push @type, $type;
	}
	return (join(",", @type));
}



