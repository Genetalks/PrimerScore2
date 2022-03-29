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
require "$Bin/product.pm";
require "$Bin/io.pm";


my $BEGIN_TIME=time();
my $version="1.0.0";
#######################################################################################

# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
my ($foligo,$ftm,$fbound,$fprobe,$fkey,$outdir);
my ($Methylation,$NoSpecificity,$NoFilter);
my $min_len=20;
my $max_len=30;
my $opt_tm=60;
my $rfloat = 0.2;
my $dis = 500;
my $range_dis="120,160,80,200"; ## dis score, pair dis range(best_min, best_max, min, max)
my $range_pos;
my $ptype = "face-to-face";
my $ctype = "Single";
my $onum = 3;
my $max_probe_num=6;
my $PCRsize=1000;
my $min_eff=0.00001;
my $max_prodn=50;
my $min_tm_spec=45;
my $max_bound_num = 10000;
GetOptions(
				"help|?" =>\&USAGE,
				"io:s"=>\$foligo,
				"it:s"=>\$ftm,
				"ib:s"=>\$fbound,
				"ip:s"=>\$fprobe,
				"k:s"=>\$fkey,
				"Probe:s"=>\$fprobe,
				"Methylation:s"=>\$Methylation,
				"NoSpecificity:s"=>\$NoSpecificity,
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
				"stm:s"=>\$min_tm_spec,
				"mine:s"=>\$min_eff,
				"maxp:s"=>\$max_prodn,
				"od:s"=>\$outdir,
				) or &USAGE;
&USAGE unless ($foligo and $fkey);
$outdir||="./";
`mkdir $outdir`	unless (-d $outdir);
$outdir=AbsolutePath("dir",$outdir);
my $min_cpgs=1;
my $min_cs=3;
my $max_best_dis_primer_probe=3;
my $max_dis_primer_probe=30;
my $end_len=10;
my @mis_end=(10,100,0,100); #end: 0-10 => score: 0-fulls
my @lendif=(0,4,0,8); ## tm diff between F and R primer
my @tmdif=(0,3,0,6); ## tm diff between F and R primer
#my ($fulls_pos, $fulls_dis, $fulls_lend, $fulls_tmd, $fulls_prod)=(25,25,10,10,30);
my ($fulls_pos, $fulls_dis, $fulls_lend, $fulls_tmd, $fulls_prod)=(20,30,10,10,30);
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
open(F,">$outdir/$fkey.primer.filter") or die $!;
print F "#Filter\tID\tTarget\tStrand\tPos5\tSeq\tLen\t",join("\t", &features("Head", $Methylation, $NoSpecificity)),"\n";
print O "##Score: ". &score_des("Primer", $Methylation)."\n";
print O "#ID\tTarget\tStrand\tPos5\tSeq\tLen\tScore\tScoreInfo\t",join("\t", &features("Head", $Methylation, $NoSpecificity)),"\n";
my %oligo_info;
my %oligo_pos;
my %oligo_score;
open(P, $foligo) or die $!;
&SHOW_TIME("#Primer Score");
while(<P>){
	chomp;
	#($id, $seq, $len), ($tm, $gc, $hairpin, $dimertype, $dimersize, $nendA, $enddG, $snp, $poly), ($CpGs, $Cs), ($bnum, $btm, $binfo);
	my ($abase, $afeature, $ameth, $aspec, $bnumtm)=&read_evaluation_info(0, $_, $Methylation, 1);
	my ($id, $seq, $len)=@{$abase};
	my $tm = $afeature->[0];
	my ($tid, $dis, $chr, $pos3, $pos5, $strand)=&get_position_info($id, \%tempos);
	@{$oligo_info{$id}}=($chr, $pos3, $pos5, $strand, $dis, $seq, $len, @{$afeature}, @{$ameth});
	if(defined $bnumtm){
		push @{$oligo_info{$id}}, $bnumtm;
	}
	my @out=($id, $chr, $strand, $pos5, $seq, $len, @{$afeature}, @{$ameth}, @{$aspec});
	if(!defined $NoFilter && ($len>$max_len || $len<$min_len)){
		print F join("\t", "Len", @out),"\n";
		next;
	}
	if(!defined $NoFilter && ($tm<$opt_tm-5 || $tm>$opt_tm+5)){
		print F join("\t", "Tm", @out),"\n";
		next;
	}
	if(!defined $NoFilter && defined $Methylation){
		my $cpgs=scalar(split /,/, $ameth->[0]);
		if($cpgs<$min_cpgs){
			print F join("\t", "CpGs", @out),"\n";
			next;
		}
		my $cs=scalar(split /,/, $ameth->[1]);
		if($cs<$min_cs){
			print F join("\t", "Cs", @out),"\n";
			next;
		}
	}

	if($ftype eq "SNP"){
		$tid=~s/-[UD]$//;
	}
	@{$oligo_pos{$tid}{$strand}{$id}}=($dis, $pos3, $pos5);
	
	## score
	my ($sadd, $score_info);
	if(defined $Methylation){
		($sadd, $score_info)=&primer_meth_score($opt_tm, $len, @{$afeature}, @{$ameth}, $bnumtm);
	}else{
		($sadd, $score_info)=&primer_oligo_score($opt_tm, $len, @{$afeature}, $bnumtm);
	}
	@{$oligo_score{$id}}=($sadd, $score_info);
	print O join("\t", $id, $chr, $strand, $pos5, $seq, $len, $sadd, $score_info, @{$afeature}, @{$ameth}, @{$aspec}),"\n";
}
close(P);
close(O);

my %bound;
if(!defined $NoSpecificity){
	&SHOW_TIME("#Read in Bound file");
	open(B, $fbound) or die $!;
	my $n=0;
	my $last_id = "NA";
	while(<B>){
		chomp;
		my ($id, $strand, $chr, $pos5, $seq, $tm, $end_match, $mvisual)=split /\t/, $_;
		if(!defined $NoFilter){
			if($id ne $last_id){
				if($n>$max_bound_num){
					print F join("\t", "BoundsTooMore", $id, $n),"\n";
					delete $bound{$last_id};
					delete $oligo_info{$last_id};
				}
				$last_id = $id;
				$n=0;
			}
			$n++;
		}
		my $len = length $seq;
		my $pos3=$strand eq "+"? $pos5+$len-1: $pos5-$len+1;
		push @{$bound{$id}{$chr}{$strand}}, [$pos3, $pos5, $tm, $end_match, $mvisual, $seq];
	}
	close(B);
	if($n>$max_bound_num){
		print F join("\t", "BoundsTooMore", $last_id, $n),"\n";
		delete $bound{$last_id};
		delete $oligo_info{$last_id};
	}
}
close(F);
	
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
print O "##ScorePair: " . &score_des("PairPrimer", $Methylation)."\n";
if(defined $fprobe){
	print O "##ScorePair(Probe): " . &score_des("PairProbe", $Methylation)."\n";
}
print O "##ScoreOligo: total score | ". &score_des("Primer", $Methylation)."\n";
if(defined $fprobe){
	print O "##ScoreOligo(Probe): total score | " . &score_des("Probe", $Methylation)."\n";
}
print O &final_head($Methylation, $NoSpecificity, $range_pos, $ptype),"\n";

my %success;
my %output;
&SHOW_TIME("#Primer Pair Score and Select");
foreach my $tid(sort {$a cmp $b} keys %{$target{"tem"}}){
	#primer conditions
	my @condv;
	if(defined $fprobe){
		my %probec=%{$probe{$tid}};
		my $n=0;
		my %recordp;
		foreach my $id (sort {$probec{$b}->[0]<=>$probec{$a}->[0]} keys %probec){
			my ($subid)=$id=~/(\S+)_\d+$/;
			next if(exists $recordp{$subid}); ## more probes with the same position will only keep one
			$recordp{$subid}=1;

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
	
	my %record;
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
			next if(!exists $oligo_info{$p1});
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
				next if(!exists $oligo_info{$p2});
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
				my $sprod=$fulls_prod;
				my %prod;
				if(!defined $NoSpecificity){
					if(!exists $bound{$p1} || !exists $bound{$p2}){
						print "Warn: $p1 or $p2 No bound info!\n";
						next;
					}
	
					&caculate_product($tid, "P1", $tid, "P2", $bound{$p1}, $bound{$p2}, $ptype, \%prod, \%record, $PCRsize, $opt_tm, $min_tm_spec, $rdis[-2], $rdis[-1], $min_eff, $max_prodn); ## 1<-->2
					$record{"pro"}{$p1}{$p2}=1;
					&caculate_product($tid, "P1", $tid, "P1", $bound{$p1}, $bound{$p1}, $ptype, \%prod, \%record, $PCRsize, $opt_tm, $min_tm_spec, $rdis[-2], $rdis[-1], $min_eff, $max_prodn);## 1<-->1
					$record{"pro"}{$p1}{$p1}=1;
					&caculate_product($tid, "P2", $tid, "P2", $bound{$p2}, $bound{$p2}, $ptype, \%prod, \%record, $PCRsize, $opt_tm, $min_tm_spec, $rdis[-2], $rdis[-1], $min_eff, $max_prodn);## 2<-->2
					$record{"pro"}{$p2}{$p2}=1;
					my @effs = sort{$b<=>$a} values %{$prod{$tid}{$tid}};
					my $pnum=scalar @effs;
					$sprod = &bound_score($pnum, join(",",@effs), $fulls_prod, "Eff");
				}
					
				# score pair
				my $stotal=$score+$score2;
				if(defined $fprobe){#Probe
					# probe num on products
					my %pdr;
					my $pb=$condv[$i][4];
					my $spdr=$fulls_prod;
					if(!defined $NoSpecificity){
						my ($pdnum)= &probe_bounds_on_products($pb, $bound{$pb}, $prod{$tid}{$tid}, \%pdr, $PCRsize, $opt_tm+8, $min_tm_spec+8); #my ($id, $abound, $aprod, $aresult)=@_;
						my @pdeffs = sort{$b<=>$a} values %{$pdr{$pb}};
						$spdr = &bound_score($pdnum, join(",", @pdeffs), $fulls_prod, "Eff"); ## probe specificity
						$sprod=($sprod*1+$spdr*2)/3; ## specificity weight is primer:probe=1:2
					}
					$stotal+=$probe{$tid}{$pb}->[0];
					@{$probe_final{$p1.",".$p2}}=($spdr, $pdr{$pb});
				}
				$stotal +=$spos+$slend+$stmd+$sdis+$sprod;

				$score_pair{$p1.",".$p2}=$stotal;
				@{$score_pair_info{$p1.",".$p2}}=($spos, $sdis, $slend, $stmd, $sprod);
				$pos_pair{$p1.",".$p2}=$pos3;

				#my $prob=join(",", $chr, $eff_dis, $sd, $pos, $eff1,$tm1, $mvisual1, $sd2, $pos2, $eff2, $tm2, $mvisual2);
				my $effs_info;
				@{$pair_info{$p1.",".$p2}}=($dis, $prod{$tid}{$tid});
			}
		}
		next if(scalar keys %score_pair==0);
		## choose according score
		my @final;
		if($ctype eq "Single"){
			my @pair_sort = sort{$score_pair{$b}<=>$score_pair{$a}} keys %score_pair;
			my %recordt;
			my $num = 0;
			#LBR-1096-D-28-12,LBR-1096-D-26-32
			for(my $i=0; $i<@pair_sort; $i++){
				my ($p1, $p2)=split /,/, $pair_sort[$i];
				my $d1=$oligo_info{$p1}->[1]; 
				my $d2=$oligo_info{$p2}->[1]; 
				my $min_dis1=100;
				my $min_dis2=100;
				foreach my $ds(keys %recordt){
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
				$recordt{$d1.",".$d2}=1;
				$num++;
				last if($num==$onum || defined $fprobe); ## probe: only one pair for one probe
			}
		}else{
			@final = &average($dis, $rfloat, \%pos_pair,\%score_pair,"UP"); ##my ($len, $dis, $rfloat, $apos, $ascore,$select)
		}
	
		##output
		my @mk=("A","B","C","D","E","F","G","H","I","J","K","L","M","N");
		my $n=0;
		foreach my $pair (@final){
			my ($size, $aprod)=@{$pair_info{$pair}};
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
				$pname=$target."-".$UD."-".$mk[$n];
				if(defined $fprobe){
					$pname=$target."-"."P"."-".$mk[$pbnum];
				}
			}else{
				$pname=$target."-".($n+1);
			}
			my $p1_new=$pname."-1";
			my $p2_new=$pname."-2";
			$output{$p1_new}=$p1;
			$output{$p2_new}=$p2;
			my @prod_info=();
			if(!defined $NoSpecificity){
				my ($pdnum, $apdeffs, $apdinfos)=&get_highest_bound($aprod, 3, "Eff");
				my $pdeffs=join(",", @{$apdeffs});
				my $pdinfos=join(";", @{$apdinfos});
				@prod_info=($pdnum, $pdeffs, $pdinfos);
			}
			my @opos=($chr, $pos5, $strand, $p1_new, $seq, $len);
			my @oinfo=($stotal, join(",", @{$score_pair_info{$pair}}), $score."|".$score_info, @info, @prod_info);
			my @opos2=($chr2, $pos52, $strand2, $p2_new, $seq2, $len2);
			my @oinfo2=($stotal, join(",", @{$score_pair_info{$pair}}), $score2."|".$score_info2, @info2, @prod_info);
			if(defined $range_pos){
				print O join("\t", @opos, $dis, $size, @oinfo),"\n";
			}else{
				print O join("\t", @opos, $size, @oinfo),"\n";
			}
			if(defined $fprobe){#probe
				my $pb=$condv[$i][4];
				my $pb_new = $target."-"."P"."-".$mk[$pbnum]."-P";
				$output{$pb_new}=$pb;
				my ($pbs, $pbsi)=@{$probe{$tid}{$pb}};
				my ($chrp, $pos3p, $pos5p, $strandp, $disp, $seqp, $lenp, @infop)=@{$oligo_info{$pb}};
				my ($spdr, $apdr)=@{$probe_final{$pair}};
				my @pdprod_info=();
				if(!defined $NoSpecificity){
					my ($pdpnum, $apdpeffs, $apdpinfos)=&get_highest_bound($apdr, 3, "Eff");
					my $pdpeffs=join(",", @{$apdpeffs});
					my $pdpinfos=join(";", @{$apdpinfos});
					@pdprod_info=($pdpnum, $pdpeffs, $pdpinfos);
				}
				$infop[6]=0; ##end stability

				my @oposp=($chrp, $pos5p, $strandp, $pb_new, $seqp, $lenp);
				my @oinfop=($stotal, $spdr, $pbs."|".$pbsi, @infop, @pdprod_info);
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

## output final bound info
open(B, ">$outdir/$fkey.final.bound.info") or die $!;
foreach my $idn(sort {$a cmp $b} keys %output){
	my $id=$output{$idn};
	foreach my $chr(sort {$a cmp $b} keys %{$bound{$id}}){
		foreach my $strand (sort {$a cmp $b} keys %{$bound{$id}{$chr}}){
			my @info = @{$bound{$id}{$chr}{$strand}};
			for(my $i=0; $i<@info; $i++){
				my ($pos3, $pos5, $tm, $end_match, $mvisual, $seq)=@{$info[$i]};
				print B join("\t", $idn, $strand, $chr, $pos3, $seq, $tm, $end_match, $mvisual),"\n";
			}
		}
	}
}
close(B);

## check failed target
open(O, ">$outdir/$fkey.design.status") or die $!;
print O "## If Failed, you can try turning -maxl/-maxlp up or -opttm/-opttmp down, or choosing --Nofilter!\n";
print O "#TargetID\tPrimerID\tStatus\n";
foreach my $t(sort {$a cmp $b} keys %{$target{"id"}}){
	if(exists $success{$t}){
		print O join("\t", $t, $target{"id"}{$t}, "Successful"),"\n";
	}else{
		print O join("\t", $t, $target{"id"}{$t}, "Failed"),"\n";
		print "Warn: $t design failed! you can try turning up the -pnum, or check file $outdir/design/$fkey.oligo.filter.list and adjust the -maxl/-maxlp/-minl/-minp or -opttm/-opttmp!\n";
	}
}
close(O);


#######################################################################################
print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################

# ------------------------------------------------------------------
# sub function
# ------------------------------------------------------------------
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


sub USAGE {#
	my $usage=<<"USAGE";
Program:
Version: $version
Contact:zeng huaping<huaping.zeng\@genetalks.com> 


-rd: distance range of pair primers, (best_min, best_max, min, max) separted by ",", example:
      face-to-face: |---> P1 x            dis_range: 120,160,80,200
         (SNP)               x  P2 <---|  

      face-to-face: P1 |--->        x     dis_range:120,160,80,200
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

  --Methylation     Design methylation oligos
  --NoSpecificity   Not evalue specificity
  --NoFilter        Not filter any primers
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
  -stm    <int>     min tm to be High_tm in specifity, [$min_tm_spec]
  -mine      <float> min efficiency to consider a product, [$min_eff]
  -maxp      <int>   maximum products number to be caculated, to reduce running time. [$max_prodn]

  -od <dir>	   Dir of output file, default ./
  -h		   Help

USAGE
	print $usage;
	exit;
}

