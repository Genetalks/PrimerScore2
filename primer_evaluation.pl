#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
require "$Bin/path.pm";
require "$Bin/self_lib.pm";

my $BEGIN_TIME=time();
my $version="1.0.0";
#######################################################################################

# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
my ($fprimer, $fkey,$detail,$outdir);
my $NoSpecificity;
my $pnum = 1;
my $min_tm_spec = 40; #when caculate specificity
my $nohead;
my $thread = 3;
our $PATH_PRIMER3;
our $REF_HG19;
my $fdatabase = $REF_HG19;
my $NoFilter;
my $len_map=20; ##bwa result is the most when 20bp
my $opt_tm = 65;
my $opt_tm_probe;
my $opt_size = 100;
my $dis_range="80,120,30,300"; ## pair dis range(best_min, best_max, min, max)
my $type = "face-to-face";
GetOptions(
				"help|?" =>\&USAGE,
				"p:s"=>\$fprimer,
				"d:s"=>\$fdatabase,
				"n:s"=>\$pnum,
				"k:s"=>\$fkey,
				"NoFilter:s"=>\$NoFilter,
				"NoSpecificity:s"=>\$NoSpecificity,
				"rdis:s"=>\$dis_range,
				"type:s"=>\$type,
				"nohead:s"=>\$nohead,
				"maplen:s"=>\$len_map,
				"opttm:s"=>\$opt_tm,
				"opttmp:s"=>\$opt_tm_probe,
				"stm:s"=>\$min_tm_spec,
				"Detail:s"=>\$detail,
				"thread:s"=>\$thread,
				"od:s"=>\$outdir,
				) or &USAGE;
&USAGE unless ($fprimer and $fkey);

$outdir||="./";
`mkdir $outdir`	unless (-d $outdir);
$outdir=AbsolutePath("dir",$outdir);
my $min_end_match = 6;
my $merge_len = 100;
my $Wind_GC = 8;
my $MAX_hairpin_tm = 55;
my $MAX_endA = 4;
my $MAX_poly = 20;
my @rank_end=(3,   5,   8); #  PCR efficiency when dis of mismatch pos to 3end <= @rank_end
my @eff_end =(0.1, 0.4, 0.8 );
my $eff_times = 10;
my $min_eff = 0.1;
my $MIN_tm = $opt_tm-5;
my $MAX_tm = defined $opt_tm_probe? $opt_tm_probe+5: $opt_tm+5;
my $MIN_gc = 0.15;
my $MAX_gc = 0.85;
my @dis = split /,/, $dis_range;
my @tm = ($opt_tm*0.8, $opt_tm*2, $opt_tm*0.6, $opt_tm*2);
my $extend = $dis[-1];
	
my $oligotm = "$PATH_PRIMER3/src/oligotm";
my $ntthal = "$PATH_PRIMER3/src/ntthal";
my $primer3_config = "$PATH_PRIMER3/src/primer3_config/";

## creat primer.fa
my @PN=();
for (my $i=0; $i<$pnum; $i++){
	open ($PN[$i], ">$outdir/$fkey.primer_$i.fa") or die $!;
}
if(!defined $NoFilter){
	open(F, ">$outdir/$fkey.filter.list") or die $!;
}
my %evalue;
my %seq;
my %map;
my %map_id;
open (P, $fprimer) or die $!;
while (<P>){
	chomp;
	next if(/^$/);
	$_=~s/\s+$//;
	my ($id, @seq)=split /\s+/, $_;
	@{$seq{$id}}=@seq;
	for (my $i=0; $i<$pnum; $i++){
		my $primer_seq = $seq[$i];
		## len Tm GC Hairpin, Dimer between itself
		if(!defined $primer_seq){
			print $_,"\n";
			die;
		}

		## filter endA and poly
		my $nendA = &get_end_A($primer_seq);
		my $vpoly = &get_poly_value($primer_seq);
		my $ftype;
		if(!defined $NoFilter && $nendA>$MAX_endA){
			$ftype = "EndA";
			print F join("\t",  $id, $primer_seq, $ftype, $nendA),"\n";
			next;
		}
		if(!defined $NoFilter && $vpoly >$MAX_poly){
			$ftype = "Poly";
			print F join("\t",  $id, $primer_seq, $ftype, $vpoly),"\n";
			next;
		}

		## filter TM and GC
		my $Tm = `$oligotm $primer_seq`;
		chomp $Tm;
		$Tm = sprintf "%0.2f", $Tm;
		my @GC_info = &GC_info_stat($primer_seq, $Wind_GC);
		if(!defined $NoFilter && ($Tm<$MIN_tm || $Tm>$MAX_tm)){
			$ftype = "Tm";
			print F join("\t",  $id, $primer_seq, $ftype, $Tm),"\n";
			next;
		}
		if(!defined $NoFilter && ($GC_info[0]<$MIN_gc || $GC_info[0]>$MAX_gc)){
			$ftype = "GC";
			print F join("\t",  $id, $primer_seq, $ftype, $GC_info[0]),"\n";
			next;
		}
		
		## filter hairpin
		my $len = length($primer_seq);
		my ($hairpin_dg, $hairpin_tm, $dimer_dg, $dimer_tm)=('','','','');
		
		my $hairpin_result = `$ntthal -path $primer3_config -a HAIRPIN -s1 $primer_seq`;
		($hairpin_dg, $hairpin_tm)=$hairpin_result=~/dG = ([\d\-\.]+)\tt = ([\d\-\.]+)/;
		if (!defined $hairpin_dg){
			$hairpin_dg = '';
			$hairpin_tm = '';
		}else{
			$hairpin_dg =  sprintf "%0.2f",$hairpin_dg;
			$hairpin_tm =  sprintf "%0.2f",$hairpin_tm;
			if(!defined $NoFilter && $hairpin_tm > $MAX_hairpin_tm){
				print F join("\t", $id, $primer_seq, "Hairpin:".$hairpin_tm),"\n";
				next;
			}
		}
		
		my $dimer_result = `$ntthal -path $primer3_config -a END1 -s1 $primer_seq -s2 $primer_seq`;
		($dimer_dg, $dimer_tm)=$dimer_result=~/dG = ([\d\-\.]+)\tt = ([\d\-\.]+)/;
		if (!defined $dimer_dg){
			$dimer_dg = '';
			$dimer_tm = '';
		}else{
			$dimer_dg =  sprintf "%0.2f",$dimer_dg;
			$dimer_tm = sprintf "%0.2f",$dimer_tm;
		}
		
		push @{$evalue{$id}},[$len, $Tm, @GC_info, $hairpin_dg, $hairpin_tm, $dimer_dg, $dimer_tm];
		
		my $id_new = "$id\_$i";
		my $off = (length $primer_seq)-$len_map;
		$off=$off>0? $off: 0;
		my $mseq = substr($primer_seq, $off);
		if(exists $map{$i}{$mseq}){
			$map_id{$id_new}=$map{$i}{$mseq};
			next;
		}
		$map{$i}{$mseq}=$id_new;
		$map_id{$id_new}=$id_new;
		print {$PN[$i]} ">$id_new\n";
		print {$PN[$i]} $mseq,"\n";
	}
}
for (my $i=0; $i<$pnum; $i++){
	close($PN[$i]);
}
if(!defined $NoFilter){
	close(F);
}
exit(0) if(scalar keys %evalue==0);

if(defined $detail){
	open (Detail, ">$outdir/$fkey.evaluation.detail") or die $!;
}
my $DB;
my %speci;
my %efficiency;
if(!defined $NoSpecificity){
	for (my $pn=0; $pn<$pnum; $pn++){
		if($pn!=$pnum-1){
			open ($DB, ">$outdir/$fkey.database_$pn") or die $!;
		}
		
		my %db_region;
		### bwa
		my $fa_primer = "$outdir/$fkey.primer_$pn.fa";
		
		Run("bwa mem -D 0 -k 9 -t $thread -c 5000000 -y 100000 -T 12 -B 1 -L 2,2 -h 200 -a  $fdatabase $fa_primer > $fa_primer.sam");
	#	Run("bwa mem -D 0 -k 9 -t 4 -c 5000000 -y 100000 -T 12 -B 1 -O 2,2 -L 1,1 -h 200 -a  $fdatabase $fa_primer > $fa_primer.sam");
		
		### read in sam
		my %mapping;
		open (I, "$fa_primer.sam") or die $!;
		while (<I>){
			chomp;
			next if(/^$/ || /^\@/ || /^\[/);
			my ($id, $flag, $chr, $pos, $score, $cigar, undef, undef, undef, $seq)=split /\s+/,$_;
			my ($is_unmap, $is_reverse)=&explain_bam_flag_unmap($flag);
			my ($md)=$_=~/MD:Z:(\S+)/;
			if ($is_unmap){
				push @{$speci{$id}}, [0,0,0,'']; #[$align_num, $high_tm_num, $end_match_num,$high_tm_info];
				print "Warn: primer $id is not mapped!\n";
				print $_,"\n";
				next;
			}
			#next if($pn!=0 && $is_reverse);
			my ($H3)=&get_3end1_mismatch($is_reverse, $cigar);
			next if($H3>1);
			push @{$mapping{$id}},[$is_reverse, $flag, $chr, $pos, $score, $cigar, $md];
		}
		close(I);
		### evaluate
		foreach my $ido (sort {$a cmp $b} keys %evalue){
			my ($primer_seq)=$seq{$ido}->[$pn];
			my $mid = $map_id{$ido."_".$pn};
			my $align_num = 0;
			my $map_num = scalar @{$mapping{$mid}};
			my @high_tm;
			my %high_info;
			for (my $i=0; $i<$map_num; $i++){
				my ($is_reverse, $flag, $chr, $pos, $score, $cigar, $md)=@{$mapping{$mid}->[$i]};

#				print join("\t", ("test", $ido, $is_reverse, $flag, $chr, $pos, $score, $cigar, $md)),"\n";
				## specificity
				my ($dG, $tm, $end_match, $ntthal_result, $pos3, $cigar_new, $md_new) = &ntthal_map($primer_seq, $chr, $pos, $is_reverse, $fdatabase);
				next if($end_match == -1); ## map to NNNN region, invalid
				($pos, $cigar, $md) = ($pos3, $cigar_new, $md_new);
#				print "new:", join(",",($tm, $pos, $cigar, $md)),"\n";

				## get pcr efficiency
				my ($last_eff, $dis);
				if($pn==0){
					$last_eff = 1;
					$dis = 0;
				}else{
					$last_eff = $efficiency{$chr}*$eff_times; ## eff_times: few DNA in the last step will be amplified, so multiply by eff_times
					$dis = $pos;
				}
				my ($eff)=&get_pcr_efficienty($tm, $dis, $cigar, $md, $min_tm_spec, $last_eff, \@rank_end, \@eff_end, $pn);

				$align_num++;
				if(($tm<$min_tm_spec && $align_num>200) || ($align_num>400)){
					$align_num.="+";
					last;
				}
				
				## store info
				my $strand = $is_reverse? "-": "+";
				if($pn==0){
					push @{$high_info{sprintf("%.2f",$tm)}},[$strand, $chr, $pos, $cigar, $md, $end_match, $eff];
				}else{
					push @{$high_info{sprintf("%.2f",$eff)}},[$strand, $chr, $pos, $cigar, $md, $end_match, $tm];
				}

				## out database
				if ($tm >= $min_tm_spec){
					if ($pn != $pnum-1){
						my $seq;
						my $spos;
						if($type eq "back-to-back"){
							$seq = &get_database_seq_back($chr, $pos, $is_reverse, $fdatabase, $extend); 
							@dis=($extend-$dis[0]/2, $extend, 0, $extend);
							my $len = length($primer_seq);
							$spos=$is_reverse? $pos+$extend-$len: $pos-$extend+$len; ## to be compatible with other type
						}else{
							$seq = &get_database_seq($primer_seq, $chr, $pos, $is_reverse, $fdatabase, $extend); #($primer, $fdatabae, $chr, $pos, $is_reverse)
							$spos=$pos;
						}
						$chr=$is_reverse? "-".$chr: "+".$chr;
						my $id_db = join(",", $ido."_".$pn, $chr, $spos); ## spos in $id_db is used to calculate primer positions in sub get_high_info
						$efficiency{$id_db}=$eff;
						print $DB ">$id_db\n$seq\n";
					}
				}
				if(defined $detail){
					print Detail join("\t",$primer_seq, join(",",$chr, $pos, $is_reverse),$dG, $tm, $end_match),"\n";
					print Detail $ntthal_result;
				}
			}

			my ($high_tm_num, $high_eff_num, $high_info);
			if($pn==0){
				($high_tm_num, $high_eff_num, $high_info)=&get_high_info(\%high_info, $min_tm_spec, $min_eff, $pn);
			}else{
				($high_eff_num, $high_tm_num, $high_info)=&get_high_info(\%high_info, $min_eff, $min_tm_spec, $pn);
			}

			if(defined $detail){
				print Detail "##",join("\t", $ido, $pn, $align_num, $high_tm_num, $high_eff_num, $high_info),"\n";
			}

			push @{$speci{$ido}}, [$align_num, $high_tm_num, $high_eff_num, $high_info];
		}
		
		if($pn !=$pnum-1){
			close($DB);
			$fdatabase = "$outdir/$fkey.database_$pn";
			if(-f "$fdatabase.amb"){ ## rm old bwa indexed files 
				`rm $fdatabase.*`;
			}
			Run("bwa index $fdatabase");
		}
		
	}
	if(defined $detail){
		close (Detail);
	}
}

### output
open (O, ">$outdir/$fkey.evaluation.out") or die $!;
if(!defined $nohead){
	print O "##High_Info(n=1): TM        : Flag/DatabaseID/Pos/Cigar/MD/End_match_num/Efficiency\n";
	if($pnum>1){
		print O "##High_Info(n>1): Efficiency: Flag/DatabaseID/Pos/Cigar/MD/End_match_num/TM\n";
	}
	print O "#ID";
	for (my $i=0; $i<$pnum; $i++){
		print O "\tSeq\tLen\tTm\tGC\tGC5\tGC3\tdGC\tdG_Hairpin\tTm_Hairpin\tdG_Dimer\tTm_Dimer\tAlign_Num\tHigh_Tm_Num\tHigh_Efficiency_Num\tHigh_Info";
	}
	print O "\n";
}

foreach my $id (sort {$a cmp $b} keys %evalue){
	print O $id;
	for (my $i=0; $i<$pnum; $i++){
		my @spe = defined $NoSpecificity? ("","","",""):@{$speci{$id}->[$i]};
		print O "\t",$seq{$id}->[$i],"\t",join("\t",@{$evalue{$id}->[$i]}),"\t",join("\t",@spe);
	}
	print O "\n";
}
close(O);
#######################################################################################
print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################

# ------------------------------------------------------------------
# sub function
# ------------------------------------------------------------------

sub get_high_info{
	my ($ahigh_info, $min_key, $min_other, $pn)=@_;
	my %high_info=%{$ahigh_info};
	
	my @high_info=();
	my @high_value=();
	my $n = 0;
	my $high_key_num = 0;
	my $high_other_num=0;
	
	my %record;
	foreach my $t (sort{$b<=>$a} keys %high_info){
		#last if($t<$min_key && $n>=3);
		last if($t<$min_key);
		for(my $i=0; $i<@{$high_info{$t}}; $i++){
			if($pn>0){
				my ($strand, $database_id, $pos, $cigar)=@{$high_info{$t}->[$i]}[0..3];
				#check strand 
				if(($type eq "Nested" && $strand eq "-") || ($type ne "Nested" && $strand eq "+")){
					print "Wrong strand: ", join(",", ($strand, $database_id, $pos, $cigar)),"\n";
					next;
				}
				
				my ($idt, $chrt_ori, $p3t)=split /,/, $database_id;
				my ($idto)=$idt=~/(\S+)\_\d+$/;
				my ($strandt,$chrt)=$chrt_ori=~/([+-])(\S+)/;

				my $plen = length ($seq{$idto}->[$pn-1]);
				my ($unmap_len5) = $cigar=~/^(\d+)H/; ## when primer 5end is not mapped, indicate that: the primer is connected to the last primer, it is not expected
				if(!defined $unmap_len5){
					$unmap_len5=0;
				}else{
					if($type eq "Nested" && $pos==1){
						print "Wrong: the last primer $database_id been connected on the 5end\n";
					}
				}
				if($type eq "Nested"){
					$p3t = $strandt eq "+"? $p3t+$pos-$plen-$unmap_len5: $p3t-$pos+$plen+$unmap_len5;
					#print "Pos3:",join("\t", ($p3t, $strand, $database_id, $pos, $cigar)),"\n";
				}else{
					$p3t = $strandt eq "+"? $p3t+$pos-$plen: $p3t-$pos+$plen;
				}
#				print join("\t","pos3:", $database_id, $p3t, $pos, $plen, $unmap_len5),"\n";
				if(exists $record{$chrt."_".$p3t}){
					my $ix = $record{$chrt."_".$p3t};
					$high_info[$ix].="@".$database_id; 
					next;
				}
				$record{$chrt."_".$p3t}=$n;
			}
			if($t >= $min_key){
				$high_key_num++;
			}
			my ($other) = @{$high_info{$t}->[$i]}[-1];
			if ($other >=$min_other){
				$high_other_num++;
			}
			if($n<3){
				push @high_value, $t;
				push @high_info, join("/",@{$high_info{$t}->[$i]});
			}
			$n++;
		}
	}
#	print Dumper %record;
	my $high_info = join(",",@high_value).":".join(";",@high_info);
	return($high_key_num, $high_other_num, $high_info);
}

sub get_pcr_efficienty{ 
	my ($tm, $dis, $cigar, $md, $min_tm, $last_eff, $arank_end, $aeff_end, $pn)=@_;
	if($tm < $min_tm){
		return 0;
	}

	# tm eff
	my $eff_tm = &score_single($tm, 1, @tm);
	$eff_tm = $eff_tm>0? $eff_tm: 0;
	# dis eff
	my $eff_dis = $pn>0? &score_single($dis, 1, @dis) : 1;
	$eff_dis = $eff_dis>0? $eff_dis: 0;

	# mismatch pos to 3end
	my @mis_pos;
	&get_3end_mismatch($cigar, $md, \@mis_pos, $arank_end->[-1]);
	my $eff_end = 1;
	for(my $i=0; $i<@mis_pos; $i++){
		$eff_end *= &get_eff_rank($mis_pos[$i], $arank_end, $aeff_end, "<=");
	}

	my $eff = sprintf("%.2f", $last_eff * $eff_tm * $eff_dis * $eff_end);
#	print "eff:\t", join("\t", $eff, $last_eff, $eff_tm, $eff_dis, $eff_end, ":", $tm, $dis, $md),"\n";

	#$eff=$eff>1? 1: $eff;
	return $eff;
}

sub get_3end_mismatch{
	my ($cigar, $md, $aresult, $end_len)=@_;
	## insert: 2 mismatch
	if($cigar=~/I/){
		$cigar=~s/\d+S//;
		{
			my @num=split /[MID]/, $cigar;
			my @minfo=split /\d+/, $cigar;
			my $sum=0;
			for(my $i=$#num; $i>=0; $i--){
				$sum+=!defined $num[$i]? 0: $num[$i];
				$sum++;
				if($sum<=$end_len && $minfo[$i] eq "I"){
					push @{$aresult}, ($sum, $sum);
				}
			}
		}
	}

	## del
	$md=~s/\^[ATCG]+/A0A/g; ## one del ==> 2 mismatch
		
	## mismatch
	my @num=split /[ATCG]+/, $md;
	my $sum=0;
	for(my $i=$#num; $i>=0; $i--){
		$sum+=!defined $num[$i]? 0: $num[$i];
		$sum++;
		if($sum<=$end_len){
			push @{$aresult}, $sum;
		}
	}
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

sub get_3end1_mismatch{
	my ($is_reverse, $cigar)=@_;
	my $H3;
	if($is_reverse==0){
		($H3)=$cigar=~/M(\d+)H$/;
	}else{
		($H3)=$cigar=~/^(\d+)H/;
	}
	if(!defined $H3){ ## 1H can be prc from experiment data
		$H3=0;
	}
	return $H3;
}

sub get_database_seq_back{
	my ($chr, $pos, $is_reverse, $fref, $extend)=@_;
	my ($s, $e); ## contain primer sequence
	if($is_reverse==0){
		$s=$pos-$extend;
		$e=$pos;
	}else{
		$s=$pos;
		$e=$s+$extend;
	}
	$s=$s>1? $s: 1;
	# get seq
	my $rout = `samtools faidx $fref $chr:$s-$e`;
	chomp $rout;
	my ($id, @seq)=split /\n/, $rout;
	my $seq = join("", @seq);
	if($is_reverse){
		$seq=~tr/ATCGatcg/TAGCtagc/;
		$seq=reverse $seq;
	}
	return ($seq);
}

sub get_database_seq_back_old{
	my ($adrange, $primer, $chr, $pos, $is_reverse, $fref)=@_;
	my $plen = length $primer;
	my $dextend = 30; ## fragments amplified by primers mapped on downstream sequence will be removed because of no cyclizing marks

#	#### downstream: not contain primer sequence
#	my ($s, $e);
#	if ($is_reverse == 0){
#		$s = $pos+1;
#		$e = $s+$dextend;
#	}else{
#		$e = $pos-1;
#		$s = $e - $dextend;
#	}
#	$s=$s>1? $s: 1;
#
#	# get seq
#	my $rout = `samtools faidx $fref $chr:$s-$e`;
#	chomp $rout;
#	my ($id, @seq)=split /\n/, $rout;
#	my $seq = join("", @seq);
#	if($is_reverse){
#		$seq=~tr/ATCGatcg/TAGCtagc/;
#		$seq=reverse $seq;
#	}
	
	#### upstream: contain primer revcom sequence
	my $bseq;
	{
		my ($bs, $be);
		if($is_reverse){
			$bs=$pos;
			$be=$bs+$extend;
		}else{
			$bs=$pos-$extend;
			$be=$pos;
		}
		$bs=$bs>1? $bs: 1;
		# get seq
		my $rout = `samtools faidx $fref $chr:$bs-$be`;
		chomp $rout;
		my ($id, @seq)=split /\n/, $rout;
		my $seq = join("", @seq);
		if($is_reverse){
			$seq=~tr/ATCGatcg/TAGCtagc/;
			$seq=reverse $seq;
		}
		$bseq = $seq; ## contain primer revcom sequence
	}

#	$seq = $seq.$bseq;

	my @dr=@{$adrange};
	@{$adrange} = ($extend+$dextend-2*$plen+$dr[2], $extend+$dextend, $dextend, $extend+$dextend); ## position of primer 3end
	return ($bseq);

}

## pos is postion of primer 3end
sub get_database_seq{
	my ($primer, $chr, $pos, $is_reverse, $fref, $extend)=@_;
	my $plen = length $primer;
	my ($s, $e);
	if ($is_reverse == 0){
		$s = $pos+1;
		$e = $s+$extend;
	}else{
		$e = $pos-1;
		$s = $e - $extend;
	}
	$s=$s>1? $s: 1;

	# get seq
	my $rout = `samtools faidx $fref $chr:$s-$e`;
	chomp $rout;
	my ($id, @seq)=split /\n/, $rout;
	my $seq = join("", @seq);
	if($is_reverse){
		$seq=~tr/ATCGatcg/TAGCtagc/;
		$seq=reverse $seq;
	}
	$seq=$primer.$seq;

	return ($seq);
}

sub Run{
    my ($cmd, $log_fh)=@_;

    if (defined $log_fh) {
        print $log_fh $cmd,"\n";
    }
    print STDERR $cmd, "\n";
    my $ret = system($cmd);
    if ($ret) {
        die "Run $cmd failed!\n";
    }
}

sub ntthal_map{
	my ($primer_seq, $chr, $pos, $is_reverse, $fdatabase)=@_;
	my $len = length($primer_seq);
	my $extend = 10;
	
	### get seq
	my ($start, $end);
	my $off = (length $primer_seq)-$len_map;
	if($is_reverse){# "-"
		$start = $pos-$extend;
		$end = $pos+$len+$extend;
	}else{
		$start = $pos-$off-$extend;
		$end = $pos-$off-1+$len+$extend;
	}
	#print join("\t", $pos, $is_reverse, $start, $end),"\n";
	$start=$start>1? $start: 1;
	my $seq_info = `samtools faidx $fdatabase $chr:$start-$end`;
	my @seq_info = split /\n/, $seq_info;
	shift @seq_info;
	my $seq = join("",@seq_info);
	if($seq!~/[ATCGatcg]+/){
		return (-1,-1,-1,-1);
	}
	$seq=uc($seq);
	if ($is_reverse == 0){
		$seq=~tr/atcg/tagc/;
		$seq=~tr/ATCG/TAGC/;
		$seq = reverse $seq;
	}
	### ntthal
	if(defined $detail){
		print Detail "$ntthal -path $primer3_config -a END1 -s1 $primer_seq -s2 $seq\n";
	}
	my $result = `$ntthal -path $primer3_config -a END1 -s1 $primer_seq -s2 $seq`;

	### get tm 
	my ($dG, $tm) = $result =~/dG = ([\d\+\-\.]+)\tt = ([\d\+\-\.]+)/;
	my @line = split /\n/, $result;
	shift @line;	
	if(!defined $line[1]){ ## just print and return
		print join("\t",$primer_seq, $chr, $pos, $is_reverse, $fdatabase),"\n";
		print "samtools faidx $fdatabase $chr:$start-$end\n";
		print "$ntthal -path $primer3_config -a END1 -s1 $primer_seq -s2 $seq\n";
		return (-1,-1,-1,-1);
	}
	## match position and cigar
	my ($pos_new, $cigar, $md)=&get_match_cigar_from_ntthal(\@line, $is_reverse, $start, $end);
	#print $result,"\n";
	#print "ntthal map:\n"; 
	#print join("\t", $pos_new, $cigar, $md,  $chr, $is_reverse, $start, $end),"\n";

	## get end_match
	my @aunit = split //, $line[1];
	my $end_match = 0;
	my $count_flag = 0;
	for (my $i=$#aunit; $i>=0; $i--){
		if ($aunit[$i] ne " "){
			$count_flag = 1;
		}
		if ($count_flag ==1){
			if($aunit[$i] ne " "){
				$end_match++;
			}else{
				last;
			}
		}
	}

	return($dG, $tm, $end_match, $result, $pos_new, $cigar, $md);
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

	## end handle
	if($cigar=~/(\d+)D$/){
		($tright)=$cigar=~/(\d+)D$/;
		$cigar=~s/\d+D$//;
		$md=~s/\^[ATCG]+$//;
	}elsif($cigar!~/M$/){
		print "Wrong:$cigar\n";
		print join("\n", @line),"\n";
		die;
	}
	if($pleft>0){
		$cigar = $pleft."S".$cigar;
	}

	my $pos_new = $is_reverse? $start+$tright: $end-$tright; ## position of  primer3
#	print "ntthal map:", join("\t", $cigar, $md, $pos_new, $tleft, $tright, $pleft, $pright),"\n";
	return($pos_new, $cigar, $md);
}

sub explain_bam_flag_unmap{
	my ($flag)=@_;
	my $flag_bin=sprintf("%b", $flag);
	my @flag_bin = split //, $flag_bin;
#	my $is_read1 = @flag_bin>=7? $flag_bin[-7]:0;
#	my $is_read2 = @flag_bin>=8? $flag_bin[-8]: 0;
#	my $is_supplementary = @flag_bin>=12? $flag_bin[-12]: 0;
#	my $is_proper_pair = @flag_bin>=2? $flag_bin[-2]:0;
	my $is_reverse = @flag_bin>=5? $flag_bin[-5]: 0;
	my $is_unmap = @flag_bin>=3? $flag_bin[-3]:0;
#	my $is_munmap = @flag_bin>=4? $flag_bin[-4]:0;
#	my $dup = @flag_bin>=11? $flag_bin[-11]: 0;
#	my @result = ($is_read1, $is_proper_pair, $is_reverse, $is_unmap, $is_munmap, $is_supplementary);
	return ($is_unmap, $is_reverse);
}


sub AbsolutePath
{		#获取指定目录或文件的决定路径
		my ($type,$input) = @_;

		my $return;
	$/="\n";

		if ($type eq 'dir')
		{
				my $pwd = `pwd`;
				chomp $pwd;
				chdir($input);
				$return = `pwd`;
				chomp $return;
				chdir($pwd);
		}
		elsif($type eq 'file')
		{
				my $pwd = `pwd`;
				chomp $pwd;

				my $dir=dirname($input);
				my $file=basename($input);
				chdir($dir);
				$return = `pwd`;
				chomp $return;
				$return .="\/".$file;
				chdir($pwd);
		}
		return $return;
}

sub GetTime {
	my ($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst)=localtime(time());
	return sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}


sub USAGE {#
	my $usage=<<"USAGE";
Program:
Version: $version
Contact:zeng huaping<huaping.zeng\@genetalks.com> 
	

Usage:
  Options:
  -p  <file>   Input primer list file, forced
  -d  <file>   Input database file, [$fdatabase]
  -n  <int>    combined primer num, single primer: 1, primer pair: 2, forced
  -k  <str>	Key of output file, forced

  --NoFilter             Not filter any primers
  --NoSpecificity        Not evalue specificity
  -type      <str>       primer type, "face-to-face", "back-to-back", "Nested", [$type]
  -rdis      <str>       distance range between pair primers when evaluate primer efficiency, (opt_min, opt_max, min, max) separted by ",", [$dis_range]
  -opttm    <int>       optimal tm of primer, [$opt_tm]
  -opttmp    <int>     optimal tm of probe, not design probe when not set the parameter, optional
  -maplen    <int>      length to map with bwa, [$len_map]
  -stm       <int>      min tm to be High_TM when caculate specificity, [$min_tm_spec]
  -thread    <int>      thread in bwa, [$thread]
  --Detail              Output Detail Info, optional
  -od        <dir>      Dir of output file, default ./
  -h		 Help

USAGE
	print $usage;
	exit;
}
