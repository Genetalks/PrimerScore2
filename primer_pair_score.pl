#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
require "$Bin/self_lib.pm";
require "$Bin/score.pm";

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
my $PCRsize=600;
GetOptions(
				"help|?" =>\&USAGE,
				"io:s"=>\$foligo,
				"it:s"=>\$ftm,
				"ib:s"=>\$fbound,
				"ip:s"=>\$fprobe,
				"k:s"=>\$fkey,
				"maxl:s"=>\$max_len,
				"minl:s"=>\$min_len,
				"opttm:s"=>\$opt_tm,
				"PCRsize:s"=>\$PCRsize,
				"tp:s"=>\$ptype,
				"on:s"=>\$onum,
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
my @enddG = (-9,0,-10,0);
my $dlen=$max_len-$min_len;
my @len = ($min_len,$max_len-$dlen*0.5,$min_len-5, $max_len+1);
my @tm = ($opt_tm-1, $opt_tm+1, $opt_tm-5, $opt_tm+5);
my @self = (-50, 40, -50, 55); ## self tm
my $max_best_dis_primer_probe=3;
my $max_dis_primer_probe=15;
my $min_eff=0.01;
my $end_len=10;
my @mis_end=(10,100,0,100); #end: 0-10 => score: 0-fulls
my $fulls=10;
my $fulls_dis = 10;
my $fulls_pos = 10;
my $max_probe_num=3;
my @rdis=split /,/, $range_dis;
my @rpos;
if(defined $range_pos){
	@rpos=split /,/, $range_pos;
}

my $min_end_match=-1;
my %tlength;
my %tempos;
my %target;
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
	}else{
		@{$tempos{$id}}=($id, 1, $len, "+", 0, 0);
		@{$target{"tem"}{$id}}=($id);
		$target{"id"}{$id}=$id;
	}
}
$/="\n";
close(F);

my %oligo_info;
my %oligo_pos;
my %oligo_score;
open(P, $foligo) or die $!;
while(<P>){
	chomp;
	my ($id, $seq, $len, $tm, $gc, $hairpin, $END, $ANY, $nendA, $enddG, $snp, $poly, $bnum, $btm)=split /\t/, $_;
	my ($tid, $dis, $chr, $pos3, $pos5, $strand)=&get_position_info($id, $len, \%tempos);
	@{$oligo_info{$id}}=($chr, $pos3, $pos5, $strand, $dis, $seq, $len, $tm, $gc, $hairpin, $END, $ANY, $nendA, $enddG, $snp, $poly, $bnum.":".$btm);
	next if($len>$max_len);
	next if($tm<$opt_tm-5 || $tm>$opt_tm+5);

	## score
	my $snendA=int(&score_single($nendA, $fulls, @nendA)+0.5);
	my $senddG=int(&score_single($enddG, $fulls, @enddG)+0.5);
	my $slen=int(&score_single($len, $fulls, @len)+0.5);## round: int(x+0.5)
	my $stm=int(&score_single($tm, $fulls, @tm)+0.5);
	my $self = &max($hairpin, $END, $ANY);
	my $sself=int(&score_single($self, $fulls, @self)+0.5);
	my $ssnp = int(&SNP_score($snp, $len, "Primer")*$fulls +0.5);
	my $spoly = int(&poly_score($poly, $len, "Primer")*$fulls +0.5);
	my $sbound=&bound_score($bnum, $btm, $fulls, "Tm");
	my @score = ($slen, $stm, $sself,$snendA, $senddG, $ssnp, $spoly, $sbound);
	my @weight =(0.5,   2.5,     1,    1,       0.5,     2,    2,      0.5);
	my $sadd=0;
	for(my $i=0; $i<@score; $i++){
#		$score[$i]=$score[$i]<0? 0: $score[$i];
		$sadd+=$weight[$i]*$score[$i];
	}
	my $score_info=join(",", @score);
	@{$oligo_score{$id}}=($sadd, $score_info);
	@{$oligo_pos{$tid}{$strand}{$id}}=($dis, $pos3, $pos5);
}
close(P);

my %bound;
open(B, $fbound) or die $!;
while(<B>){
	chomp;
	my ($id, $strand, $chr, $pos3, $seq, $tm, $end_match, $mvisual)=split /\t/, $_;
	next if($end_match<$min_end_match);
	my $len = length $seq;
	my $pos5=$strand eq "+"? $pos3-$len: $pos3+$len;
	push @{$bound{$id}{$chr}{$strand}}, [$pos3, $pos5, $tm, $end_match, $mvisual];
}
close(B);

my %probe;
if(defined $fprobe){
	open(B, $fprobe) or die $!;
	while(<B>){
		chomp;
		next if(/^#/);
		my ($id, $seq, $len, $score, $score_info)=split /\t/, $_;
		my ($tid)=$id=~/(\w+)-[FR]/;
		@{$probe{$tid}{$id}}=($score, $score_info);
	}
	close(B);
}


#### score pair and output
open(O,">$outdir/$fkey.final.result") or die $!;
my @title=("#Chr\tStart\tStrand\tID\tSeq\tLen");
if(defined $range_pos){
	push @title, "Dis2Target";
}
if($ptype=~/face-to-face/){
	push @title, "ProductSize";
}
push @title, ("ScoreTotal\tScore_PosDisSpec\tScore\tScoreInfo");
push @title, ("Tm\tGC\tHairpin\tEND_Dimer\tANY_Dimer\tEndANum\tEndStability\tSNP\tPoly\tOligoBound\tBoundNum\tHighestTm\tHighestInfo");
print O join("\t", @title),"\n";


foreach my $tid(keys %{$target{"tem"}}){
	#primer conditions
	my @condv;
	if(defined $fprobe){
		my %probec=%{$probe{$tid}};
		my $n=0;
		my %record;
		foreach my $id (sort {$probec{$b}->[0]<=>$probec{$a}->[0]} keys %probec){
			my ($subid)=$id=~/(\w+)_\d+$/;
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
				$pos3_max = $pos5+$max_dis_primer_probe;
				$pos3_min = $pos5+1;

			}
			push @condv, [$strand, "Pos3", $pos3_min.",".$pos3_max, $pos3_bmin.",".$pos3_bmax, $id];
			print join("\t", @{$oligo_info{$id}}, "Pos3", $pos3_min.",".$pos3_max, $pos3_bmin.",".$pos3_bmax),"\n";
			$n++;
			last if($n==$max_probe_num);
		}
	}else{
		if(defined $range_pos){ ## SNP
			@condv=(["+", "Dis", $rpos[-1].",".$rpos[-2]], ["-", "Dis", $rpos[-1].",".$rpos[-2]]);
		}else{ ## input file is fasta
			@condv=(["+", "No"], ["-", "No"]);
		}
	}
	
	my %primer_eff;
	for(my $i=0; $i<@condv; $i++){
		## get candidate primer1
		my @primer1 = &get_candidate($condv[$i], $oligo_pos{$tid});
		print "######", join("\t", @{$condv[$i]}),"\n";
		print Dumper @primer1;
		my %score_pair;
		my %score_pair_info;
		my %pair_info;
		my %pos_pair;
		foreach my $p1(@primer1){
			my ($chr, $pos3, $pos5, $strand, $dis, $seq)=@{$oligo_info{$p1}};
			my ($score, $score_info)=@{$oligo_score{$p1}};
			#score for pos
			my $spos = $fulls_pos;
			if($condv[$i][1] eq "Pos3"){ #probe
				my ($min, $max)=split /,/, $condv[$i][2];
				my ($bmin, $bmax)=split /,/, $condv[$i][3];
				$spos=int(&score_single($pos3, $fulls_pos, ($bmin, $bmax, $min, $max))+0.5);
			}elsif($condv[$i][1] eq "Dis"){ #SNP
				$spos=int(&score_single($dis, $fulls_pos, @rpos)+0.5);
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
			print "####", join("\t", $p1, $chr, $pos3, $pos5, $strand),"\n";
			print join(",", @condv2),"\n";
			print Dumper @primer2;
			foreach my $p2(@primer2){
				my ($chr, $pos32, $pos52, $strand2, $dis2, $seq2)=@{$oligo_info{$p2}};
				my ($score2, $score_info2)=@{$oligo_score{$p2}};
				print join("\t", $p2, @{$oligo_info{$p2}}),"\n";
				# score for dis
				my $dis;
				if($ptype eq "Nested"){
					$dis=$strand eq "+"? $pos32-$pos3: $pos3-$pos32;
				}else{
					$dis=$strand eq "+"? $pos52-$pos5: $pos5-$pos52;
				}
				my $sdis=int(&score_single($dis, $fulls, @rdis)+0.5);
				print "sdis:",join("\t", $dis, $sdis),"\n";
				# specificity, product 
				my %prod;
				my ($pnum) = &caculate_product($p1, $p2, \%bound, $ptype, \%prod, \%primer_eff);
				print "Product: $pnum\n";
				print Dumper %prod;
				my @effs = sort{$b<=>$a} values %prod;
				my $sprod = &bound_score($pnum, join(",",@effs), $fulls, "Eff");
				
				# score pair
				my $stotal = $score+$score2+$spos+$sdis+$sprod;
				if($condv[$i][1] eq "Pos3"){#Probe
					my $pb=$condv[$i][4];
					$stotal+=$probe{$tid}{$pb}->[0];
				}
				$score_pair{$p1.",".$p2}=$stotal;
				@{$score_pair_info{$p1.",".$p2}}=($stotal, $spos, $sdis, $sprod);
				
				$pos_pair{$p1.",".$p2}=$pos3;

				#my $prob=join(",", $chr, $eff_dis, $sd, $pos, $eff1,$tm1, $mvisual1, $sd2, $pos2, $eff2, $tm2, $mvisual2);
				my $effs_info;
				@{$pair_info{$p1.",".$p2}}=($dis, \%prod);
				print join("\t", $p1.",".$p2, $stotal, $score, $score_info, $score2, $score_info2, $spos, $sdis, $sprod),"\n";
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
		my $n=0;
		foreach my $pair (@final){
			my ($size, $aprod)=@{$pair_info{$pair}};
			my ($stotal, $spos, $sdis, $sprod)=@{$score_pair_info{$pair}};
			my ($p1, $p2)=split /,/, $pair;
			my ($chr, $pos3, $pos5, $strand, $dis, $seq, $len, @info)=@{$oligo_info{$p1}};
			my ($chr2, $pos32, $pos52, $strand2, $dis2, $seq2, $len2, @info2)=@{$oligo_info{$p2}};
			my ($score, $score_info)=@{$oligo_score{$p1}};
			my ($score2, $score_info2)=@{$oligo_score{$p2}};
			if($condv[$i][1] ne "Dis"){
				$dis="NA";
			}
			my ($target)=$p1=~/(\w+)-[FR]/;
			my $UD=$strand eq "+"? "U":"D";
			my $p1_new = $target."-".$UD."-".$pos3."-P1";
			my $p2_new = $target."-".$UD."-".$pos3."-P2";
				
			my ($pdnum, $pdeffs, $pdinfos)=&get_highest_bound($aprod, 3);
			print O join("\t", $chr, $pos3, $strand, $p1_new, $seq, $len, $size, $stotal, join(",", $spos, $sdis, $sprod), $score, $score_info, @info, $pdnum, $pdeffs, $pdinfos),"\n";
			if($condv[$i][1] eq "Pos3"){#probe
				my $pb=$condv[$i][4];
				my $pb_new = $target."-".$UD."-".$pos3."-Probe";
				my ($pbs, $pbsi)=@{$probe{$tid}{$pb}};
				my ($chrp, $pos3p, $pos5p, $strandp, $disp, $seqp, $lenp, @infop)=@{$oligo_info{$pb}};
				print O join("\t", $chrp, $pos3p, $strandp, $pb_new, $seqp, $lenp, "NA", $stotal, "NA", $pbs, $pbsi, @infop),"\n";
			}
			print O join("\t", $chr2, $pos32, $strand2, $p2_new, $seq2, $len2, $size, $stotal, join(",", $spos, $sdis, $sprod), $score2, $score_info2, @info2, $pdnum, $pdeffs, $pdinfos),"\n";
		}
	}
}
close(O);



#######################################################################################
print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################

# ------------------------------------------------------------------
# sub function
# ------------------------------------------------------------------
#push @{$bound{$id}{$chr}{$strand}}, [$pos3, $pos5, $tm, $end_match, $mvisual];
sub caculate_product{
	my ($p1, $p2, $abound, $ptype, $aprod, $aeff)=@_;
	my %bd=%{$abound};
	my $prodn=0;
	print ">>", join(",", $p1, $p2),"\n";
#	print Dumper %{$bd{$p1}};
	foreach my $chr(keys %{$bd{$p1}}){
		foreach my $sd(keys %{$bd{$p1}{$chr}}){
			my @bds=@{$bd{$p1}{$chr}{$sd}};
			for(my $i=0; $i<@bds; $i++){
				my ($tm1, $end_match1, $mvisual1)=@{$bds[$i]}[2..4];
##              Nested      face-to-face         back-to-back
##P1			-->|          |-->                    |-->
##P2min			-->|           <--|             <--|...  
##P2max			  ...-->|        ...<--|                 ...<--|   
				my ($ixpos, $dmin, $dmax);
				if($ptype eq "Nested"){
					$ixpos=0; #pos3
					$dmin=0;
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
				print "##", join("\t", $chr, $sd, $pos, "=>", $sd2, $pmin, $pmax),"\n";
				next if(!exists $bd{$p2}{$chr}{$sd2});
				my @bds2=@{$bd{$p2}{$chr}{$sd2}};
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
						print join("\t", $eff,$eff1,$eff2,$eff_dis,$dis, $prob),"\n";
						next if($eff<$min_eff);
						$prodn++;
						$aprod->{$prob}=$eff;
					}
				}
			}
		}
	}
	return ($prodn);
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
			$pmin=$pos-$dis_max;
			$pmax=$pos-$dis_min;
		}else{
			$pmin=$pos+$dis_min;
			$pmax=$pos+$dis_max;
		}
	}
	return ($sd, $pmin, $pmax);
}

#my @primer1 = &get_candidate($tid,$condkey[$i], $condv[$i], \%oligo);
#@{$oligo_pos{$tid}{$strand}{$id}}=($dis, $pos3, $pos5);
sub get_candidate{
	my ($acondv, $aoligo)=@_;
	my ($strand, @cond)=@{$acondv};
	my $ix;
	if($cond[0] eq "No"){
		return (keys %{$aoligo->{$strand}});
	}elsif($cond[0] eq "Pos3"){
		$ix=1;
	}elsif($cond[0] eq "Pos5"){
		$ix=2;
	}elsif($cond[0] eq "Dis"){
		$ix=0;
	}
	my @oligo = sort{$aoligo->{$strand}{$a}->[$ix] <=> $aoligo->{$strand}{$b}->[$ix]} keys %{$aoligo->{$strand}};

	my ($min, $max)=split /,/, $cond[1];
	my @final;
	for(my $i=0; $i<@oligo; $i++){
		next if($aoligo->{$strand}{$oligo[$i]}->[$ix]<$min);
		last if($aoligo->{$strand}{$oligo[$i]}->[$ix]>$max);
		push @final, $oligo[$i];
	}
	return @final;
}

sub get_position_info{
	my ($id, $plen, $atpos)=@_;
	my ($tid, $tori, $startt, $endt, $pori, $off)=$id=~/(\w+)-([FR])-(\d+)-(\d+)_([FR])_(\d+)/; ## primer pos
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
		$pos3=$start+$pos3t-1;
		$pos5=$start+$pos5t-1;
		$strand = $strandp;
	}else{
		$pos3=$end-$pos3t+1;
		$pos5=$end-$pos5t+1;
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
      face-to-face: |---> P1 x            dis_range(P2-P1): 180,220,150,250
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

            Nested: P2 --->|   x          dis_range(P2-P1):10,15,1,30
                      P1 --->| x

Usage:
  Options:
  -io  <file>   Input oligo features file, forced
  -it  <file>   Input template file, forced
  -ib  <file>   Input oligo bounded info file, forced
  -ip  <file>   Input probe score file, optional
  -k   <str>	Key of output file, forced
  -tp  <str>    primer type, "face-to-face", "back-to-back", "Nested", ["face-to-face"]

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

  -od <dir>	   Dir of output file, default ./
  -h		   Help

USAGE
	print $usage;
	exit;
}

