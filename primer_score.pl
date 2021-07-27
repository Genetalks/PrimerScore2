#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
require "$Bin/common.pm";
require "$Bin/score.pm";
require "$Bin/self_lib.pm";
require "$Bin/average.pm";
require "$Bin/math.pm";

my $BEGIN_TIME=time();
my $version="1.0.0";
#######################################################################################

# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
my ($foligo,$ftm,$fbound,$fprobe,$fkey,$outdir);
my $min_len=20;
my $max_len=30;
my $opt_tm=60;
my $rfloat = 0.2;
my $dis = 500;
my $range_dis="10,15,1,30"; ## dis score, pair dis range(best_min, best_max, min, max)
my $range_pos;
my $ptype = "face-to-face";
my $ctype = "Single";
my $onum = 3;
my $max_probe_num=6;
my $PCRsize=600;
my $NoFilter;
GetOptions(
				"help|?" =>\&USAGE,
				"io:s"=>\$foligo,
				"it:s"=>\$ftm,
				"ib:s"=>\$fbound,
				"ip:s"=>\$fprobe,
				"k:s"=>\$fkey,
				"NoFilter:s"=>\$NoFilter,
				"maxl:s"=>\$max_len,
				"minl:s"=>\$min_len,
				"opttm:s"=>\$opt_tm,
				"PCRsize:s"=>\$PCRsize,
				"tp:s"=>\$ptype,
				"on:s"=>\$onum,
				"pn:s"=>\$max_probe_num,
				"rd:s"=>\$range_dis,
				"rp:s"=>\$range_pos,
				"ct:s"=>\$ctype,
				"ds:s"=>\$dis,
				"rf:s"=>\$rfloat,
				"od:s"=>\$outdir,
				) or &USAGE;
&USAGE unless ($foligo and $fkey);
$outdir||="./";
`mkdir $outdir`	unless (-d $outdir);
$outdir=AbsolutePath("dir",$outdir);

my @nendA = (0, 0, 0, 3);
my @enddG = (-8,0,-10,0);
my $dlen=$max_len-$min_len;
my @len = ($min_len,$max_len-$dlen*0.5,$min_len-5, $max_len+1);
my @tm = ($opt_tm-1, $opt_tm+1, $opt_tm-5, $opt_tm+5);
my @self = (-50, 40, -50, 55); ## self tm
my $max_best_dis_primer_probe=3;
my $max_dis_primer_probe=20;
my $min_eff=0.01;
my $end_len=10;
my @mis_end=(10,100,0,100); #end: 0-10 => score: 0-fulls
my @lendif=(0,3,0,8); ## tm diff between F and R primer
my @tmdif=(0,1,0,2); ## tm diff between F and R primer
my ($fulls_pos, $fulls_dis, $fulls_lend, $fulls_tmd, $fulls_prod)=(25,15,10,20,30);
my $fulls=10;
my @rdis=split /,/, $range_dis;
my @rpos;
if(defined $range_pos){
	@rpos=split /,/, $range_pos;
}
my $ftype;
my $min_end_match=-1;
my %tlength;
my %tempos;
my %target;
&SHOW_TIME("#Read in template file");
open(F, $ftm) or die $!;
$/=">";
while(<F>){
	chomp;
	next if(/^$/);
	my ($head, @seq)=split /\n/, $_;
	my ($id) = split /\s+/, $head;
	my $seq = join("", @seq);
	my $len = length $seq;
	$tlength{$id}=$len;

	if($head=~/XP:Z:/){
		my (undef, $xs, $xe, $xp, $xi)=split /\t/, $head;
		$xs=~s/XS:i://;
		$xe=~s/XE:i://;
		$xi=~s/XI:Z://;
		my ($strand, $chr, $s, $e)=$xp=~/XP:Z:([+-])(\S+):(\d+)-(\d+)/; #>NF2-401-U      XS:i:1  XE:i:1  XP:Z:+chr22:30038027-30038227   XD:Z:NF2,+,401,COSM4414719;
		@{$tempos{$id}}=($chr, $s, $e, $strand, $xs, $xe);
		my @targets=split /;/,$xi;
		$id=~s/-[UD]//;
		@{$target{"tem"}{$id}}=@targets;
		foreach my $t(@targets){
			$target{"id"}{$t}=$id;
		}
		$ftype="SNP";
	}else{
		@{$tempos{$id}}=($id, 1, $len, "+", 0, 0);
		@{$target{"tem"}{$id}}=($id);
		$target{"id"}{$id}=$id;
		$ftype="fasta";
	}
}
$/="\n";
close(F);
open(O,">$outdir/$fkey.primer.score") or die $!;
print O "##Score: scores of length, tm, self-complementary, end3 A num, end3 stability, snp, poly, bounding\n";
my %oligo_info;
my %oligo_pos;
my %oligo_score;
open(P, $foligo) or die $!;
&SHOW_TIME("#Primer Score");
while(<P>){
	chomp;
	my ($id, $seq, $len, $tm, $gc, $hairpin, $END, $ANY, $nendA, $enddG, $snp, $poly, $bnum, $btm, $binfo)=split /\t/, $_;
	my ($tid, $dis, $chr, $pos3, $pos5, $strand)=&get_position_info($id, $len, \%tempos);
	@{$oligo_info{$id}}=($chr, $pos3, $pos5, $strand, $dis, $seq, $len, $tm, $gc, $hairpin, $END, $ANY, $nendA, $enddG, $snp, $poly, $bnum."|".$btm);
	next if(!defined $NoFilter && ($len>$max_len || $len<$min_len));
	next if(!defined $NoFilter && ($tm<$opt_tm-5 || $tm>$opt_tm+5));

	## score
	my $snendA=int(&score_single($nendA, $fulls, @nendA)+0.5);
	my $senddG=int(&score_single($enddG, $fulls, @enddG)+0.5);
	my $slen=int(&score_single($len, $fulls, @len)+0.5);## round: int(x+0.5)
	my $stm=int(&score_single($tm, $fulls, @tm)+0.5);
	my $self = &max($hairpin, $END, $ANY);
	my $sself=int(&score_single($self, $fulls, @self)+0.5);
	my ($snpv)=split /:/, $snp; 
	my $ssnp = int(&SNP_score($snpv, $len, "Primer")*$fulls +0.5);
	my $spoly = int(&poly_score($poly, $len, "Primer")*$fulls +0.5);
	my $sbound=&bound_score($bnum, $btm, $fulls, "Tm");
	my @score = ($slen, $stm, $sself,$snendA, $senddG, $ssnp, $spoly, $sbound);
	my @weight =(0.5,   2.5,     1,    1,       1,     2,    1.5,      0.5);
	my $sadd=0;
	for(my $i=0; $i<@score; $i++){
#		$score[$i]=$score[$i]<0? 0: $score[$i];
		$sadd+=$weight[$i]*$score[$i];
	}
	my $score_info=join(",", @score);
	@{$oligo_score{$id}}=($sadd, $score_info);
	if($ftype eq "SNP"){
		$tid=~s/-[UD]$//;
	}
	@{$oligo_pos{$tid}{$strand}{$id}}=($dis, $pos3, $pos5);
	print O join("\t", $id, $chr, $strand, $pos5, $seq, $len, $sadd, $score_info, $tm, $gc, $hairpin, $END, $ANY, $nendA, $enddG, $snp, $poly, $bnum, $btm, $binfo),"\n";
}
close(P);
close(O);

my %bound;
&SHOW_TIME("#Read in Bound file");
open(B, $fbound) or die $!;
while(<B>){
	chomp;
	my ($id, $strand, $chr, $pos3, $seq, $tm, $end_match, $mvisual)=split /\t/, $_;
	next if($end_match<$min_end_match);
	my $len = length $seq;
	my $pos5=$strand eq "+"? $pos3-$len+1: $pos3+$len-1;
	push @{$bound{$id}{$chr}{$strand}}, [$pos3, $pos5, $tm, $end_match, $mvisual];
}
close(B);

my %probe;
if(defined $fprobe){
	&SHOW_TIME("#Read in Probe file");
	open(B, $fprobe) or die $!;
	while(<B>){
		chomp;
		next if(/^#/);
		my ($id, $seq, $len, $score, $score_info)=split /\t/, $_;
		my $tid;
		if($id=~/-[UD]-/){
			($tid)=$id=~/^(\S+)-[UD]-/;
		}else{
			($tid)=$id=~/^(\S+)-[FR]/;
		}
		@{$probe{$tid}{$id}}=($score, $score_info);
	}
	close(B);
}

#### score pair and output
open(O,">$outdir/$fkey.final.result") or die $!;
print O "##ScorePair(Primer): scores of P1 position to target(Probe 5'end), distance between pair, length diff, tm diff, specificity of primer pair\n";
print O "##ScorePair(Probe) : score of probe boundings\n";
#($slen, $stm, $sself,$snendA, $senddG, $ssnp, $spoly, $sbound)
print O "##ScoreOligo(Primer): total score | scores of length, tm, self-complementary, end3 A num, end3 stability, snp, poly, bounding\n";
print O "##ScoreOligo(Probe) : total score | scores of length, tm, self-complementary, CG content diff, snp, poly, bounding\n";
my @title=("#Chr\tStart\tStrand\tID\tSeq\tLen");
if(defined $range_pos){
	push @title, "Dis2Target";
}
if($ptype=~/face-to-face/){
	push @title, "ProductSize";
}else{
	push @title, "DisBetweenPairs";
}
push @title, ("ScoreTotal\tScorePair\tScoreOligo");
push @title, ("Tm\tGC\tHairpin\tEND_Dimer\tANY_Dimer\tEndANum\tEndStability\tSNP\tPoly\tOligoBound\tBoundNum\tHighestTm\tHighestInfo");
print O join("\t", @title),"\n";

my %success;
&SHOW_TIME("#Primer Pair Score and Select");
foreach my $tid(sort {$a cmp $b} keys %{$target{"tem"}}){
	#primer conditions
	my @condv;
	if(defined $fprobe){
		my %probec=%{$probe{$tid}};
		my $n=0;
		my %record;
		foreach my $id (sort {$probec{$b}->[0]<=>$probec{$a}->[0]} keys %probec){
			my ($subid)=$id=~/(\S+)_\d+$/;
			next if(exists $record{$subid}); ## more probes with the same position will only keep one
			$record{$subid}=1;

			my ($chr,$pos3, $pos5, $strand)=@{$oligo_info{$id}};
			my ($pos3_min, $pos3_max, $pos3_bmin, $pos3_bmax);
			if($strand eq "+"){
				$pos3_bmin = $pos5-$max_best_dis_primer_probe;
				$pos3_bmax = $pos5-1;
				$pos3_min = $pos5-$max_dis_primer_probe;
				$pos3_max = $pos5-1;
			}else{
				$pos3_bmax = $pos5+$max_best_dis_primer_probe;
				$pos3_bmin = $pos5+1;
				$pos3_max = $pos5+$max_dis_primer_probe; ## loose region
				$pos3_min = $pos5+1;
			}
			if(!defined $range_pos){##SNP
				push @condv, [$strand, "Pos3", $pos3_min.",".$pos3_max, $pos3_bmin.",".$pos3_bmax, $id];
			}else{
				push @condv, [$strand, "Dis,Pos3", $rpos[-2].",".$rpos[-1].";".$pos3_min.",".$pos3_max, $pos3_bmin.",".$pos3_bmax, $id];
			}
			$n++;
	#		last if($n==100*$max_probe_num);
		}
	}else{
		if(defined $range_pos){ ## SNP
			@condv=(["+", "Dis", $rpos[-2].",".$rpos[-1]], ["-", "Dis", $rpos[-2].",".$rpos[-1]]);
		}else{ ## input file is fasta
			@condv=(["+", "No", 0]);
		}
	}
	
	my %primer_eff;
	my $pbnum=0;
	for(my $i=0; $i<@condv; $i++){
		## get candidate primer1
		my @primer1 = &get_candidate($condv[$i], $oligo_pos{$tid});
		my %score_pair;
		my %score_pair_info;
		my %pair_info;
		my %pos_pair;
		my %probe_final;
		foreach my $p1(@primer1){
			my ($chr, $pos3, $pos5, $strand, $dis_tg, $seq, $len, $tm)=@{$oligo_info{$p1}};
			my ($score, $score_info)=@{$oligo_score{$p1}};
			#score for pos
			my $spos = $fulls_pos;
			if($condv[$i][1] eq "Pos3"){ #probe
				my ($min, $max)=split /,/, $condv[$i][2];
				my ($bmin, $bmax)=split /,/, $condv[$i][3];
				$spos=int(&score_single($pos3, $fulls_pos, ($bmin, $bmax, $min, $max))+0.5);
			}elsif($condv[$i][1] eq "Dis"){ #SNP
				$spos=int(&score_single($dis_tg, $fulls_pos, @rpos)+0.5);
			}
			#candidate p2
			my ($pos, $pform);
			if($ptype eq "Nested"){
				$pform="Pos3";
				$pos=$pos3;
			}else{
				$pform="Pos5";
				$pos=$pos5;
			}
			my ($sd, $pmin, $pmax)=&primer2_scope($ptype, $strand, $pos, $rdis[-2], $rdis[-1]);
			my @condv2=($sd, $pform, $pmin.",".$pmax);
			
			my @primer2=&get_candidate(\@condv2, $oligo_pos{$tid});
			foreach my $p2(@primer2){
				my ($chr, $pos32, $pos52, $strand2, $dis_tg2, $seq2, $len2, $tm2)=@{$oligo_info{$p2}};
				my ($score2, $score_info2)=@{$oligo_score{$p2}};
				# score for tm diff 
				my $tmdif=abs($tm2-$tm);
				my $stmd=int(&score_single($tmdif, $fulls_tmd, @tmdif)+0.5);
				# score for len diff 
				my $lendif=abs($len2-$len);
				my $slend=int(&score_single($lendif, $fulls_lend, @lendif)+0.5);

				# score for dis
				my $dis;
				if($ptype eq "Nested"){
					$dis=$strand eq "+"? $pos32-$pos3: $pos3-$pos32;
				}else{
					$dis=$strand eq "+"? $pos52-$pos5: $pos5-$pos52;
				}
				my $sdis=int(&score_single($dis, $fulls_dis, @rdis)+0.5);
				
				# specificity, product 
				my %prod;
				my ($pnum) = &caculate_product($p1, $p2, \%bound, $ptype, \%prod, \%primer_eff);
				my @effs = sort{$b<=>$a} values %prod;
				my $sprod = &bound_score($pnum, join(",",@effs), $fulls_prod, "Eff");
				
				# score pair
				my $stotal=$score+$score2;
				if(defined $fprobe){#Probe
					# probe num on products
					my %pdr;
					my $pb=$condv[$i][4];
					my ($pdnum)= &probe_bounds_on_products($pb, \%bound, \%prod, \%pdr); #my ($id, $abound, $aprod, $aresult)=@_;
					my @pdeffs = sort{$b<=>$a} values %pdr;
					my $spdr = &bound_score($pdnum, join(",", @pdeffs), $fulls_prod, "Eff"); ## probe specificity
					$sprod=($sprod*1+$spdr*2)/3; ## specificity weight is primer:probe=1:2
					$stotal+=$probe{$tid}{$pb}->[0];
					@{$probe_final{$p1.",".$p2}}=($spdr, \%pdr);
				}
				$stotal +=$spos+$slend+$stmd+$sdis+$sprod;

				$score_pair{$p1.",".$p2}=$stotal;
				@{$score_pair_info{$p1.",".$p2}}=($spos, $sdis, $slend, $stmd, $sprod);
				$pos_pair{$p1.",".$p2}=$pos3;

				#my $prob=join(",", $chr, $eff_dis, $sd, $pos, $eff1,$tm1, $mvisual1, $sd2, $pos2, $eff2, $tm2, $mvisual2);
				my $effs_info;
				@{$pair_info{$p1.",".$p2}}=($dis, \%prod);
			}
		}
		next if(scalar keys %score_pair==0);
		## choose according score
		my @final;
		if($ctype eq "Single"){
			my @pair_sort = sort{$score_pair{$b}<=>$score_pair{$a}} keys %score_pair;
			my %record;
			my $num = 0;
			#LBR-1096-D-28-12,LBR-1096-D-26-32
			for(my $i=0; $i<@pair_sort; $i++){
				my ($p1, $p2)=split /,/, $pair_sort[$i];
				my $d1=$oligo_info{$p1}->[1]; 
				my $d2=$oligo_info{$p2}->[1]; 
				my $min_dis1=100;
				my $min_dis2=100;
				foreach my $ds(keys %record){
					my ($d1r, $d2r)=split /,/, $ds;
					my $dis1 = abs($d1-$d1r);
					my $dis2 = abs($d2-$d2r);
					if($dis1<$min_dis1){
						$min_dis1 = $dis1;
					}
					if($dis2<$min_dis2){
						$min_dis2 = $dis2;
					}
				}
				next if($min_dis1<3 || $min_dis2<3 );
				push @final, $pair_sort[$i];
				$record{$d1.",".$d2}=1;
				$num++;
				last if($num==$onum);
			}
		}else{
			@final = &average($dis, $rfloat, \%pos_pair,\%score_pair,"UP"); ##my ($len, $dis, $rfloat, $apos, $ascore,$select)
		}
	
		##output
		my @mk=("A","B","C","D","E","F","G");
		my $n=0;
		foreach my $pair (@final){
			my ($size, $aprod, $apdr)=@{$pair_info{$pair}};
			my $stotal=$score_pair{$pair};
			my ($p1, $p2)=split /,/, $pair;
			my ($chr, $pos3, $pos5, $strand, $dis, $seq, $len, @info)=@{$oligo_info{$p1}};
			my ($chr2, $pos32, $pos52, $strand2, $dis2, $seq2, $len2, @info2)=@{$oligo_info{$p2}};
			my ($score, $score_info)=@{$oligo_score{$p1}};
			my ($score2, $score_info2)=@{$oligo_score{$p2}};
			my ($target)=$tid;
			my $pname;
			if($ctype eq "Single"){
				my $UD=$strand eq "+"? "U":"D";
				$UD=defined $fprobe? "B".($pbnum+1): $UD; ##Probe i
				$pname=$target."-".$UD."-".$mk[$n];
			}else{
				$pname=$target."-".($n+1);
			}
			my $p1_new=$pname."-P1";
			my $p2_new=$pname."-P2";
				
			my ($pdnum, $pdeffs, $pdinfos)=&get_highest_bound($aprod, 3);
			my @opos=($chr, $pos5, $strand, $p1_new, $seq, $len);
			my @oinfo=($stotal, join(",", @{$score_pair_info{$pair}}), $score."|".$score_info, @info, $pdnum, $pdeffs, $pdinfos);
			my @opos2=($chr2, $pos52, $strand2, $p2_new, $seq2, $len2);
			my @oinfo2=($stotal, join(",", @{$score_pair_info{$pair}}), $score2."|".$score_info2, @info2, $pdnum, $pdeffs, $pdinfos);
			if(defined $range_pos){
				print O join("\t", @opos, $dis, $size, @oinfo),"\n";
			}else{
				print O join("\t", @opos, $size, @oinfo),"\n";
			}
			if(defined $fprobe){#probe
				my $pb=$condv[$i][4];
				my $pb_new = $target."-"."B".($pbnum+1)."-Probe";
				my ($pbs, $pbsi)=@{$probe{$tid}{$pb}};
				my ($chrp, $pos3p, $pos5p, $strandp, $disp, $seqp, $lenp, @infop)=@{$oligo_info{$pb}};
				my ($spdr, $apdr)=@{$probe_final{$pair}};
				my ($pdpnum, $pdpeffs, $pdpinfos)=&get_highest_bound($apdr, 3);
				$infop[6]=0; ##end stability

				my @oposp=($chrp, $pos5p, $strandp, $pb_new, $seqp, $lenp);
				my @oinfop=($stotal, $spdr, $pbs."|".$pbsi, @infop, $pdpnum, $pdpeffs, $pdpinfos);
				if(defined $range_pos){
					print O join("\t", @oposp, $disp, 0, @oinfop),"\n";
				}else{
					print O join("\t", @oposp, 0, @oinfop),"\n";
				}
			}
			if(defined $range_pos){
				print O join("\t", @opos2, $dis2, $size, @oinfo2),"\n";
			}else{
				print O join("\t", @opos2, $size, @oinfo2),"\n";
			}
			
			foreach my $t(@{$target{"tem"}{$tid}}){
				$success{$t}=1;
			}
			$n++;
		}
		if($n>0){
			$pbnum++;
		}
		last if($pbnum==$max_probe_num);
	}
}
close(O);


## check failed target
open(O, ">$outdir/$fkey.design.status") or die $!;
print O "## If Failed, you can try turning -maxl/-maxlp up or -opttm/-opttmp down, or choosing --Nofilter!\n";
print O "#TargetID\tPrimerID\tStatus\n";
foreach my $t(sort {$a cmp $b} keys %{$target{"id"}}){
	if(exists $success{$t}){
		print O join("\t", $t, $target{"id"}{$t}, "Successful"),"\n";
	}else{
		print O join("\t", $t, $target{"id"}{$t}, "Failed"),"\n";
		print "Warn: $t design failed! you can try turning the -maxl/-maxlp up or turn the -opttm/-opttmp down, or amplify -regions!\n";
	}
}
close(O);


#######################################################################################
print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################

# ------------------------------------------------------------------
# sub function
# ------------------------------------------------------------------

#my $prob=$dis."/".join(",", $chr, $sd.$pos, $mvisual1,sprintf("%.2f",$tm1), $sd2.$pos2,$mvisual2,sprintf("%.2f",$tm2));
sub probe_bounds_on_products{
	my ($id, $abound, $aprod, $aresult)=@_;
	my %bdid=%{$abound->{$id}};
	my $n=0;
	foreach my $prod(keys %{$aprod}){
		my ($dis, $info)=split /\//, $prod;
		my ($chr, $sdpos, $mv, $tm, $sdpos2, $mv2, $tm2)=split /,/, $info;
		my ($sd, $pos)=$sdpos=~/([+-])(\d+)/;
		my ($sd2, $pos2)=$sdpos2=~/([+-])(\d+)/;
		my @sd = ($sd, $sd2);
		my @pos= ($pos, $pos2);
		next if(!exists $bdid{$chr});
		for(my $i=0; $i<2; $i++){
			my ($sd0, $pos0)=($sd[$i], $pos[$i]);
			next if(!exists $bdid{$chr}{$sd0});
			my @bds=@{$bdid{$chr}{$sd0}};
			for(my $i=0; $i<@bds; $i++){
				my ($pos3, $pos5, $tm, $end_match, $mvisual)=@{$bds[$i]};
				if(($sd0 eq "+" && $pos5>$pos0-10 && $pos5<$pos0+$PCRsize) || ($sd0 eq "-" && $pos5<$pos0+10 && $pos5>$pos0-$PCRsize)){
					my $dis=abs($pos5-$pos0);
					my $eff=&efficiency_probe($tm, $dis);
					my $bpd=$dis."/".join(",", $chr, $sd.$pos5,$mvisual,sprintf("%.2f",$tm));
					$aresult->{$bpd}=$eff;
					$n++;
				}
			}
		}
	}
	return $n;
}

#push @{$bound{$id}{$chr}{$strand}}, [$pos3, $pos5, $tm, $end_match, $mvisual];
sub caculate_product{
	my ($p1, $p2, $abound, $ptype, $aprod, $aeff)=@_;
	my %bdp1=%{$abound->{$p1}};
	my %bdp2=%{$abound->{$p2}};
	my $prodn=0;
	foreach my $chr(keys %bdp1){
		foreach my $sd(keys %{$bdp1{$chr}}){
			my @bds=@{$bdp1{$chr}{$sd}};
			for(my $i=0; $i<@bds; $i++){
				my ($tm1, $end_match1, $mvisual1)=@{$bds[$i]}[2..4];
##              Nested      face-to-face         back-to-back
##P1			-->|          |-->                    |-->
##P2min			-->|           <--|             <--|...  
##P2max			  ...-->|        ...<--|                 ...<--|   
				my ($ixpos, $dmin, $dmax);
				if($ptype eq "Nested"){
					$ixpos=0; #pos3
					$dmin=$rdis[-2];
					$dmax=$PCRsize-$min_len;
				}elsif($ptype eq "back-to-back"){
					$ixpos=1; #pos5
					$dmin=-1*$PCRsize+2*$min_len; #min_len: primer min len
					$dmax=$PCRsize;
				}else{ ##face-to-face
					$ixpos=1; #pos5
					$dmin=30;
					$dmax=$PCRsize;
				}
				my $pos=$bds[$i][$ixpos];
				my ($sd2, $pmin, $pmax)=&primer2_scope($ptype,$sd,$pos,$dmin,$dmax);
				next if(!exists $bdp2{$chr}{$sd2});
				my @bds2=@{$bdp2{$chr}{$sd2}};
				for(my $j=0; $j<@bds2; $j++){
					my ($tm2, $end_match2, $mvisual2)=@{$bds2[$j]}[2..4];
					my $pos2 = $bds2[$j][$ixpos];
					if($pos2>=$pmin && $pos2<=$pmax){
						##
						my ($eff1, $eff2);
						my $b1=join(",", $p1, $chr, $sd, $bds[$i][0]);
						if(exists $aeff->{$b1}){
							$eff1=$aeff->{$b1};
						}else{
							$eff1=&efficiency($tm1, $mvisual1);
							$aeff->{$b1}=$eff1;
						}
						my $b2=join(",", $p2, $chr, $sd2, $bds2[$j][0]);
						if(exists $aeff->{$b2}){
							$eff2=$aeff->{$b2};
						}else{
							$eff2=&efficiency($tm2, $mvisual2);
							$aeff->{$b2}=$eff2;
						}
						my $dis=abs($pos2-$pos);
						my @rsize=($rdis[-2],$rdis[-1], $dmin, $dmax);
						my $eff_dis=&score_single($dis, 1, @rsize);
						my $eff=$eff1*$eff2*$eff_dis;
						my $prob=$dis."/".join(",", $chr, $sd.$pos, $mvisual1,sprintf("%.2f",$tm1), $sd2.$pos2,$mvisual2,sprintf("%.2f",$tm2));
						next if($eff<$min_eff);
						$aprod->{$prob}=$eff;
						$prodn++;
					}
				}
			}
		}
	}
	return ($prodn);
}

sub efficiency_probe{
	my ($tm, $dis)=@_;
	# tm eff
	my $opt=$opt_tm+10;
	my @etm=($opt-1, $opt+100, 45, $opt+100); 
	my $eff_tm = &score_single($tm, 1, @etm);
	$eff_tm = $eff_tm>0? $eff_tm: 0;
	
	# dis to primer 3end
	my @dis=(1,10,-5,$PCRsize);
	my $eff_dis=&score_single($dis, 1, @dis);
	$eff_dis=$eff_dis>0? $eff_dis: 0;
	return $eff_tm*$eff_dis;
}



sub efficiency{
	my ($tm, $mvisual)=@_;
	# tm eff
	my @etm=($opt_tm-1, $opt_tm+100, 40, $opt_tm+100); 
	my $eff_tm = &score_single($tm, 1, @etm);
	$eff_tm = $eff_tm>0? $eff_tm: 0;

	# mismatch pos to 3end
	my @mis_pos;
	&get_3end_mismatch($mvisual, \@mis_pos, $end_len);
	my $eff_end = 1;
	for(my $i=0; $i<@mis_pos; $i++){
		$eff_end *= &score_single($mis_pos[$i], 1, @mis_end);
	}
	return $eff_tm*$eff_end;
}

sub get_3end_mismatch{
	my ($mv, $apos, $max_end)=@_;
	
	$mv=~s/\^+/\*\*/; ## Del ==> 2 mismatch
	$mv=~s/\-+/\*\*/; ## Insert ==> 2 mismatch
	my @unit = split //, $mv;
	my $num=0;
	for(my $i=$#unit; $i>=0; $i--){
		$num++;
		if($unit[$i] eq "*" || $unit[$i] eq "#"){
			push @{$apos}, $num;
		}
	}
}


#my ($sd, $pmin, $pmax)=&primer2_scope($type, $strand,$pos, $dis_min, $dis_max);
sub primer2_scope{
	my ($type, $strand, $pos, $dis_min, $dis_max)=@_;
	my ($sd, $pmin, $pmax);
	if($type ne "Nested"){ ## strand: opposite; Dis: from pos5
		if($strand eq "+"){
			$sd="-";
			$pmin=$pos+$dis_min;
			$pmax=$pos+$dis_max;
		}else{
			$sd="+";
			$pmin=$pos-$dis_max;
			$pmax=$pos-$dis_min;
		}
	}else{ ## Nested: Dis from pos3
		$sd=$strand;
		if($strand eq "+"){
			$pmin=$pos+$dis_min; #dis <0
			$pmax=$pos+$dis_max;
		}else{
			$pmin=$pos-$dis_max;
			$pmax=$pos-$dis_min;
		}
	}
	return ($sd, $pmin, $pmax);
}

#my @primer1 = &get_candidate($tid,$condkey[$i], $condv[$i], \%oligo);
#@{$oligo_pos{$tid}{$strand}{$id}}=($dis, $pos3, $pos5);
sub get_candidate{
	my ($acondv, $aoligo)=@_;
	my ($strand, @conds)=@{$acondv};
	my @cond=split /,/, $conds[0];
	my @value=split /;/, $conds[1];
	my (@oligo, @final);
	for(my $i=0; $i<@cond; $i++){
		my $ix;
		if($cond[$i] eq "No"){
			return (keys %{$aoligo->{$strand}});
		}elsif($cond[$i] eq "Pos3"){
			$ix=1;
		}elsif($cond[$i] eq "Pos5"){
			$ix=2;
		}elsif($cond[$i] eq "Dis"){
			$ix=0;
		}
		if($i==0){
			@oligo = sort{$aoligo->{$strand}{$a}->[$ix] <=> $aoligo->{$strand}{$b}->[$ix]} keys %{$aoligo->{$strand}};
		}else{
			@oligo=@final;
		}
		my ($min, $max)=split /,/, $value[$i];
		@final=();
		for(my $i=0; $i<@oligo; $i++){
			next if($aoligo->{$strand}{$oligo[$i]}->[$ix]<$min);
			last if($aoligo->{$strand}{$oligo[$i]}->[$ix]>$max);
			push @final, $oligo[$i];
		}

	}
	return @final;
}

sub get_position_info{
	my ($id, $plen, $atpos)=@_;
	my ($tid, $tori, $startt, $endt, $pori, $off)=$id=~/^(\S+)-([FR])-(\d+)_(\d+)_([FR])_(\d+)$/; ## primer pos
	if(!defined $tid){
		die "Wrong oligo ID: $id\n";
	}
	my $strandp = $tori ne $pori? "-": "+";
	my ($tidt, $start, $end, $strandt)=@{$atpos->{$tid}}; ## templet pos
	my ($pos3t, $pos5t); # pos3/5 on tid
	if($strandp eq "+"){
		$pos3t = $endt;
		$pos5t = $pos3t-$plen+1;
	}else{
		$pos3t = $startt;
		$pos5t = $pos3t+$plen-1;
	}
	my $dis = $pos3t; #dis to target
	my ($pos3, $pos5, $strand);
	if($strandt eq "+"){
		$pos3=$start+$pos3t;
		$pos5=$start+$pos5t;
		$strand = $strandp;
	}else{
		$pos3=$end-$pos3t;
		$pos5=$end-$pos5t;
		$strand = $strandp eq "+"? "-": "+";
	}
	return ($tid, $dis, $tidt, $pos3, $pos5, $strand);
}


sub USAGE {#
	my $usage=<<"USAGE";
Program:
Version: $version
Contact:zeng huaping<huaping.zeng\@genetalks.com> 


-rd: distance range of pair primers, (best_min, best_max, min, max) separted by ",", example:
      face-to-face: |---> P1 x            dis_range: 180,220,150,250
         (SNP)               x  P2 <---|  

      face-to-face: P1 |--->        x     dis_range:100,150,70,200
        (Region)     x        <---| P2    
                     ________________

      back-to-back:  x <---| P2           dis_range:5,10,0,15
                      P1 |---> x          (Overlap between p1 and p2: dis > 0)
                     __________

      back-to-back:  x <---| P2           dis_range:-50,-40,-60,-30
                           P1 |---> x     (No overlap between p1 and p2: dis < 0)
                     _______________

            Nested: P2 --->|   x          dis_range(P2-P1):-15,-10,-30,-3
                      P1 --->| x                   (dis < 0)

Usage:
  Options:
  -io  <file>   Input oligo features file, forced
  -it  <file>   Input template file, forced
  -ib  <file>   Input oligo bounded info file, forced
  -ip  <file>   Input probe score file, optional
  -k   <str>	Key of output file, forced
  -tp  <str>    primer type, "face-to-face", "back-to-back", "Nested", ["face-to-face"]

  --NoFilter    Not filter any primers
  -minl      <int>  min len of primer, [$min_len]
  -maxl      <int>  max len of primer, [$max_len]
  -opttm     <int>  opt tm of primer, [$opt_tm]
  -PCRsize   <int>  max PCR fragment size, [$PCRsize]
  -rd    <str>	distance range of pair primers, (best_min, best_max, min, max) separted by ",", [$range_dis]
  -rp    <str>	position range(best_min, best_max, min, max), distance of p1 to the detected site when designing oligos for target spot, optional
  -ct    <str>	    primer covered type, "Single" or "Full-covered", ["Single"]
  	 -ds  <int>		distance when -ct "Full-covered", [500]
	 -rf  <float>	distance float ratio when -ct "Full-covered", [0.2]
	 -on  <int>     output num when -ct is "Single",[$onum]
	 -pn  <int>     output probe num when design probe,[$max_probe_num]

  -od <dir>	   Dir of output file, default ./
  -h		   Help

USAGE
	print $usage;
	exit;
}

