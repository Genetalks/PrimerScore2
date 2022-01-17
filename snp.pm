
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
			return "Error";
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
		if(exists $snp{$u[$i]}){ ##SNP 
			push @snp, $i."S1";
		}elsif($u[$i] eq "E"){
			my $l=1;
			while($i+1<@u){
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
##Ref   Alt        to   
#G       A         ==> R(degenerate)
#G       A,T       ==> D(degenerate)
#GA      AT        ==> RW(MNP is several SNPs)
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
		if($polym eq "S" || $polym eq "P"){ ##SNP MNP
			$alt=~s/<[xX]>//;
			my @uref=split //, $ref;
			my @ualt=split //, $alt;
			for(my $i=$off; $i<@uref; $i++){
				next if($uref[$i] eq $ualt[$i]);
				my $dg=&snp_to_degenerate($uref[$i].",".$ualt[$i]);
				if(!defined $dg || $dg eq "Error"){
					print STDERR join("\t", $uref[$i], $ualt[$i], $ref, $alt, $type, $off),"\n";
					return;
				}
				$aseq->[$p-1+$i]=$dg;
			}
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
		if($lenr == $lena){ 
			if($lenr==1){ ##SNP
				$type="S_0_0";
			}else{ ##MNP
				my @unitr=split //, $ref;
				my @unita=split //, $at;
				my $off=0;
				for(my $i=0; $i<@unitr; $i++){
					last if($unitr[$i] ne $unita[$i]);
					$off++;
				}
				my $l=$lenr-$off;
				$type="P_".$l."_".$off;
			}
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



