
my @dg  =("R",   "Y",   "M",   "K",   "S",   "W",   "H",     "B",     "V",     "D",     "N");
my @atcg=("A/G", "C/T", "A/C", "G/T", "C/G", "A/T", "A/C/T", "C/G/T", "A/C/G", "A/G/T", "A/C/G/T");
sub snp_to_degenerate{
	my ($info)=@_;
	$info=~s/,$//;
	my @snp=split /,/, $info;
	foreach my $b(@snp){
		if($b!~/[ATCGatcg]/){
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
## p: poly(snp or indel) position, count from 1
## ref: ref base, ref in vcf
## alt: alt base, alt in vcf
sub sequence_convert_poly{
	my ($aseq, $p, $ref, $alt)=@_;
	my @Ins=("I", "J", "L", "L");

	my $type=&mutant_type($ref, $alt);
	my @types = split /,/, $type;
	my ($poly, $len, $off)=("S", 0, 0);
	my $num=0;
	foreach my $tp(@types){##usually D will not be with S and I.
		$num++;
		my ($t, $l, $o)=split /_/, $tp;
		$poly=$t;
		if($l>$len){
			$len=$l;
			$off=$o;
		}
	}
	if($p-1>=0){
		if($poly eq "S"){
			$alt=~s/<[xX]>//;
			my $dg=&snp_to_degenerate($ref.",".$alt);
			$aseq->[$p-1]=$dg;
		}elsif($poly eq "I"){
			my $ix=($len-1)>=3? 3: ($len-1);
			$aseq->[$p-1]=$Ins[$ix];
		}
	}
	if($poly eq "D"){
		for(my $i=$p+$off; $i<($p+$off+$len); $i++){
			if($i-1>=0){
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



