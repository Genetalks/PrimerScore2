#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
require "$Bin/path.pm";
require "$Bin/common.pm";
require "$Bin/self_lib.pm";
require "$Bin/snp.pm";
require "$Bin/product.pm";
require "$Bin/dimer.pm";

my $BEGIN_TIME=time();
my $version="1.0.0";
######################################################################################

our $PATH_PRIMER3;
our $REF_GRCh37;
our $SAMTOOLS;
our $BWA;
# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
my ($foligo, $fkey,$detail,$outdir);
my ($NoSpecificity,$NoFilter);
my $min_tm_spec = 45; #when caculate specificity
my $nohead;
my $thread = 1;
my $fdatabases = $REF_GRCh37;
my $len_map=20; ##bwa result is the most when 20bp
my $opt_tm = 60;
my $opt_tm_probe=70;
my $opt_size = 100;
my $Debug;
my $probe;
my ($mv, $dv, $dNTP, $dna, $tp, $sc)=(50, 1.5, 0.6, 50, 1, 1);
my $olens;
my $revcom;
my $max_time=120;
GetOptions(
				"help|?" =>\&USAGE,
				"p:s"=>\$foligo,
				"d:s"=>\$fdatabases,
				"k:s"=>\$fkey,
				"Revcom:s"=>\$revcom,
				"Probe:s"=>\$probe,
				"NoFilter:s"=>\$NoFilter,
				"maxtime:s"=>\$max_time,
				"NoSpecificity:s"=>\$NoSpecificity,
				"nohead:s"=>\$nohead,
				"maplen:s"=>\$len_map,
				"olens:s"=>\$olens,
				"opttm:s"=>\$opt_tm,
				"opttmp:s"=>\$opt_tm_probe,
				"stm:s"=>\$min_tm_spec,
				"mv:s"=>\$mv,
				"dv:s"=>\$dv,
				"dNTP:s"=>\$dNTP,
				"dna:s"=>\$dna,
				"tp:s"=>\$tp,
				"sc:s"=>\$sc,
				"Detail:s"=>\$detail,
				"thread:s"=>\$thread,
				"Debug:s"=>\$Debug,
				"od:s"=>\$outdir,
				) or &USAGE;
&USAGE unless ($foligo and $fkey);

$outdir||="./";
`mkdir $outdir`	unless (-d $outdir);
$outdir=AbsolutePath("dir",$outdir);
my $merge_len = 100;
my $Wind_GC = 8;
my $MAX_hairpin_tm = 60; ##  from experience
my $Max_Dimer_tm = 30;
my $Max_SNP_ratio = 0.071; ## 1/14
my $Max_poly_ratio = 0.8;
my $Max_poly_G = 5;
my $Max_poly_ATC = 8; 
my $MAX_endA = 4;
my $MAX_end_stability=9;
my $MAX_poly = 20;
my $Max_Bound_Num = 800; ## too high tm aligns will be filtered
my $MIN_tm = $opt_tm-5;
my $MAX_tm = defined $opt_tm_probe? $opt_tm_probe+5: $opt_tm+5;
my $MIN_gc = 0.2;
my $MAX_gc = 0.8;
my @tm = ($opt_tm*0.8, $opt_tm*2, $opt_tm*0.6, $opt_tm*2);
	
my $oligotm = "$PATH_PRIMER3/src/oligotm -mv $mv -dv $dv -n $dNTP -d $dna -tp $tp -sc $sc";
my $ntthal = "$PATH_PRIMER3/src/ntthal -mv $mv -dv $dv -n $dNTP -d $dna";

## creat oligo.fa
open(PN, ">$outdir/$fkey.oligo.fa") or die $!;
open(F, ">$outdir/$fkey.filter.list") or die $!;
open(L, ">$outdir/$fkey.oligo.list") or die $!;
my %evalue;
my %map;
my %map_id;
my %olen_oligo;
open (P, $foligo) or die $!;
#&SHSHOW_TIME("Routine analysis:");
while (<P>){
	chomp;
	next if(/^$/);
	$_=~s/\s+$//;
	my ($id0, $seq)=split /\s+/, $_;
	my ($oligo_seq0, $oligo_seq_snp0) = split /:/, $seq;
	$oligo_seq0=uc $oligo_seq0;
	if(!defined $oligo_seq_snp0){
		$oligo_seq_snp0=$oligo_seq0;
	}
	$oligo_seq_snp0=uc $oligo_seq_snp0;
	if(!defined $oligo_seq0){
		print $_,"\n";
		die;
	}
	my $id0new=defined $olens? $id0."_F_0": $id0;
	push @{$olen_oligo{$id0}}, [$id0new,"+", 0, $oligo_seq0, $oligo_seq_snp0];
	my ($oligo_seq0r, $oligo_seq_snp0r);
	if(defined $revcom){
		$oligo_seq0r=&revcom($oligo_seq0);
		$oligo_seq_snp0r=&revcom($oligo_seq_snp0);
		my $id0newR=defined $olens? $id0."_R_0": $id0;
		push @{$olen_oligo{$id0}}, [$id0newR, "-", 0, $oligo_seq0r, $oligo_seq_snp0r];
	}
	if(defined $olens){
		my ($min, $max, $scale)=split /,/, $olens;
		my $len0 = length($oligo_seq0);
		for(my $l=$min; $l<$max; $l+=$scale){
			last if($l>=$len0);
			my $off=$len0-$l;
			my $id=$id0."_F"."_".$off;
			my $seq = substr($oligo_seq0, $off);
			my $seq_snp = substr($oligo_seq_snp0, $off);
			push @{$olen_oligo{$id0}}, [$id, "+", $off, $seq, $seq_snp];
			if(defined $revcom){
				$id=$id0."_R"."_".$off;
				my $seq = substr($oligo_seq0r, $off);
				my $seq_snp = substr($oligo_seq_snp0r, $off);
				push @{$olen_oligo{$id0}}, [$id,"-", $off, $seq, $seq_snp];
			}
		}
	}

	## filter
	## oligo 3end evalue
	my $nendA0 = &get_end_A($oligo_seq0);
	my $dG_end30 = &get_end3_detaG($oligo_seq0, 5);
	my ($nendA0r, $dG_end30r);
	if(defined $revcom){
		$nendA0r = &get_end_A($oligo_seq0r);
		$dG_end30r = &get_end3_detaG($oligo_seq0r, 5);
	}
	my $is_all_filter=1;
	for(my $i=0; $i<@{$olen_oligo{$id0}}; $i++){
		my ($id, $ori, $off, $oligo_seq, $oligo_seq_snp)=@{$olen_oligo{$id0}->[$i]};
		my $len = length($oligo_seq);
		print L $id,"\n";
		my $ftype;
		## filter endA and end stability
		my ($nendA, $dG_end3) = $ori eq "+"? ($nendA0, $dG_end30): ($nendA0r, $dG_end30r);
		if(!defined $probe && !defined $NoFilter){## Primer filter
			$ftype="No";
			if($nendA>$MAX_endA){
				$ftype = "EndA";
			}
			if($dG_end3 > $MAX_end_stability){
				$ftype = "EndStability";
			}
			if($ftype ne "No"){
				print F join("\t",  $id, $oligo_seq, $ftype, $nendA),"\n";
				next;
			}
		}
		
		## filter TM and GC
		my $Tm = `$oligotm $oligo_seq `;
		chomp $Tm;
		$Tm = sprintf "%0.2f", $Tm;
		my $GC=&GC($oligo_seq);
		if(!defined $NoFilter && ($Tm<$MIN_tm || $Tm>$MAX_tm)){
			$ftype = "Tm";
			print F join("\t",  $id, $oligo_seq, $ftype, $Tm),"\n";
			next;
		}
		if(!defined $NoFilter && ($GC<$MIN_gc || $GC>$MAX_gc)){
			$ftype = "GC";
			print F join("\t",  $id, $oligo_seq, $ftype, $GC),"\n";
			next;
		}
		
		## filter Self_Complementarity
		my $hairpin_tm = `$ntthal -a HAIRPIN -s1 $oligo_seq -r`;
		chomp $hairpin_tm;
		if(!defined $NoFilter && $hairpin_tm > $MAX_hairpin_tm){
			print F join("\t", $id, $oligo_seq, "Hairpin", $hairpin_tm),"\n";
			next;
		}
		my $ENDinfo = `$ntthal -a END1 -s1 $oligo_seq -s2 $oligo_seq`;
		my ($END_tm, $END31, $END32)=&dimer_amplify($ENDinfo);
		if(!defined $NoFilter && ($END_tm > $Max_Dimer_tm && $END31+$END32==0)){
			print F join("\t", $id, $oligo_seq, "END_Dimer",join(",",$END_tm, $END31,$END32)),"\n";
			next;
		}
		my $ANYinfo = `$ntthal -a ANY -s1 $oligo_seq -s2 $oligo_seq`;
		my ($ANY_tm, $ANY31, $ANY32)=&dimer_amplify($ANYinfo);
		if(!defined $NoFilter && ($ANY_tm > $Max_Dimer_tm && $ANY31+$ANY32==0)){
			print F join("\t", $id, $oligo_seq, "ANY_Dimer",join(",", $ANY_tm, $ANY31, $ANY32)),"\n";
			next;
		}
		
		## filter SNP
		my ($SNP_num, $SNP_info)=&SNP_parse($oligo_seq_snp);
		$SNP_info=$SNP_num==0? "NA":$SNP_info.":".$oligo_seq_snp;
#		$SNP_info.=",".$oligo_seq_snp;
		if(!defined $NoFilter && $SNP_num/$len > $Max_SNP_ratio){
			print F join("\t", $id, $oligo_seq, "SNP",$SNP_info),"\n";
			next;
		}
	
		## filter poly
		my ($total, $max_len, $max_base, $poly_info)=&poly_check($oligo_seq);
		$poly_info = "NA" if($total==0);
		if(!defined $NoFilter && $total>0){
			if($total/$len > $Max_poly_ratio || ($max_base eq "G" && $max_len>$Max_poly_G) || ($max_base ne "G" && $max_len>$Max_poly_ATC)){
				print F join("\t", $id, $oligo_seq, "Poly",$poly_info),"\n";
				next;
			}
		}
		push @{$evalue{$id}},($oligo_seq, $len, $Tm, sprintf("%.3f",$GC), sprintf("%.2f",$hairpin_tm), sprintf("%.2f",$END_tm), sprintf("%.2f",$ANY_tm),$nendA,$dG_end3, $SNP_info, $poly_info);
		$is_all_filter=0;
	}
	
	
	if($is_all_filter==0){
		my $off = (length $oligo_seq0)-$len_map;
		$off=$off>0? $off: 0;
		my $mseq = substr($oligo_seq0, $off);
		print PN ">$id0\n";
		print PN $mseq,"\n";
		if(defined $revcom){
			my $mseqL = substr($oligo_seq0, 0, $len_map); ##
			print PN ">$id0\_L\n";
			print PN $mseqL,"\n";
		}
	}
	
}
close(PN);

exit(0) if(scalar keys %evalue==0);
foreach my $id(keys %evalue){
	print L "Final: $id\n";
}
close(L);
#&SHSHOW_TIME("Specificity analysis:");

if(defined $detail){
	open (Detail, ">$outdir/$fkey.evaluation.detail") or die $!;
}
my $DB;
my %bound;
if(!defined $NoSpecificity){
	my %db_region;
	### bwa
	my %mapping;
	my $fa_oligo = "$outdir/$fkey.oligo.fa";
	my @fdatabase=split /,/, $fdatabases;	
	my $sum=0;
	foreach my $fdatabase(@fdatabase){
		if(!-e "$fdatabase\.ann"){
			my $fdbname = basename($fdatabase);
			my $fdatabase_new = "$outdir/$fdbname";
			if(!-e "$fdatabase_new\.ann"){
				`ln -s $fdatabase $fdatabase_new`;
				`$BWA index $fdatabase_new`;
			}
			$fdatabase=$fdatabase_new;
		}
		my $dname = basename($fdatabase);
		my $cmd="$BWA mem -D 0 -k 9 -t $thread -c 1000000000 -y 1000000000 -T 12 -B 1 -L 2,2 -h 200 -a  $fdatabase $fa_oligo >$fa_oligo\_$dname.sam 2>$fa_oligo\_$dname.sam.log";
		&Run_monitor_timeout($max_time, $cmd);
		open(I, "$fa_oligo\_$dname.sam.log") or die $!;
		my $ret = `grep -aR Killed $fa_oligo\_$dname.sam.log`;
		chomp $ret;
		if($ret eq "Killed"){## time out
			open(I, $fa_oligo) or die $!;
			$/=">";
			while(<I>){
				chomp;
				next if(/^$/);
				my ($id0, $seq)=split /\n/, $_; 
				next if(!exists $olen_oligo{$id0});	
				for(my $i=0; $i<@{$olen_oligo{$id0}}; $i++){
					my ($idn, $ori, $off, $pseqn)=@{$olen_oligo{$id0}->[$i]};
					print F join("\t",  $idn, $pseqn, "BwaTimeout", $max_time."s"),"\n";
				}
			}
			close(I);
			print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
			exit(0); ## once time out, all primer is filtered with flag "BwaTimeout", and exit!
			$/="\n";
		}
		### read in sam
		open (I, "$SAMTOOLS view $fa_oligo\_$dname.sam|") or die $!;
		while (<I>){
			chomp;
			my ($id, $flag, $chr, $pos, $score, $cigar, undef, undef, undef, $seq)=split /\s+/,$_;
			my ($is_unmap, $is_reverse)=&explain_bam_flag_unmap($flag);
			my ($md)=$_=~/MD:Z:(\S+)/;
			next if ($is_unmap);
			my ($H3)=&get_3end1_mismatch($is_reverse, $cigar);
			#print "H3:", join("\t", $H3,$id, $is_reverse, $flag, $chr, $pos, $score, $cigar, $md),"\n";
			next if(!defined $probe && !defined $revcom && $H3>1);
			$sum++;
			push @{$mapping{$id}{$dname}},[$is_reverse, $flag, $chr, $pos, $score, $cigar, $md, $fdatabase];
		}
		close(I);
	}
	### evaluate
	open(O, ">$outdir/$fkey.bound.info") or die $!;
	foreach my $id0 (sort {$a cmp $b} keys %olen_oligo){
		my ($id,undef, undef, $oligo_seq)=@{$olen_oligo{$id0}->[0]};
		my %lowtm;
		my $filter_by_bound=0;
		foreach my $dname(keys %{$mapping{$id0}}){
			my @id0=($id0);
			if(defined $revcom){
				push @id0, $id0."_L";
			}
			foreach my $id0t(@id0){
				my $bound_num = 0;
				my $map_num = scalar @{$mapping{$id0t}{$dname}};
				for (my $i=0; $i<$map_num; $i++){
					my ($is_reverse, $flag, $chr, $pos, $score, $cigar, $md, $fdatabase)=@{$mapping{$id0t}{$dname}->[$i]};
					#filter lowtm
#					my $map_form=&map_form_standard($is_reverse, $cigar, $md);
					my $map_form=join(",", $is_reverse, $cigar, $md);
					next if(exists $lowtm{$map_form});## filter bound regions whose map info are same with where bound tm too low 
					my $strand=$is_reverse? "-": "+";
					if(defined $detail){
						print Detail "\nOriginal:",join("\t",$id0t, $is_reverse, $flag, $chr, $pos, $score),"\n";
					}
					######## specificity
					my $len = length($oligo_seq);
					my $extend = 40;
					
					### get seq
					my ($start, $end);
					my $off = $len-$len_map;
					if($is_reverse){# "-"
						$start = $pos-$extend;
						$end = $pos+$len+$extend;
					}else{
						$start = $pos-$off-$extend;
						$end = $pos-$off-1+$len+$extend;
					}
					#print join("\t", $pos, $is_reverse, $start, $end),"\n";
					$start=$start>1? $start: 1;
					my $seq_info = `$SAMTOOLS faidx $fdatabase $chr:$start-$end`;
					my @seq_info = split /\n/, $seq_info;
					shift @seq_info;
					my $seq = join("",@seq_info);
					if($seq!~/[ATCGatcg]+/){ ## seq is NNN...NNN, false mapping, because bwa will convert N to one of ATCG randomly
						print "Warn: extract aligned region sequence failed\n";
						print join("\t", $chr, $start, $end, $fdatabase),"\n";
						next;
					}
					$seq=uc($seq);
					$seq=&revcom($seq)if(!$is_reverse);
					### ntthal
					my $result = `$ntthal -a ANY -s1 $oligo_seq -s2 $seq`;
					if(defined $detail){
						print Detail "$ntthal -a ANY -s1 $oligo_seq -s2 $seq\n";
						print Detail $result;
					}
				
					### get tm 
					my ($dG, $tm) = $result =~/dG = ([\d\+\-\.]+)\tt = ([\d\+\-\.]+)/;
					next if($dG eq ""); ## map to NNNN region, invalid
					if($tm<$min_tm_spec){
						if($map_form=~/H/){
							$lowtm{$map_form}=1;
						}
						next;
					}
					my @line = split /\n/, $result;
					shift @line;	
					if(!defined $line[1]){ ## just print and return
						print join("\t",$oligo_seq, $chr, $pos, $is_reverse, $fdatabase),"\n";
						print "$SAMTOOLS faidx $fdatabase $chr:$start-$end\n";
						print "$ntthal -a ANY -s1 $oligo_seq -s2 $seq\n";
						return (-1,-1,-1,-1);
					}
					## match visual
					my ($mvisual, $pos3, $pos5)=&map_visual_from_ntthal(\@line, $is_reverse, $start, $end);
					my ($end_match3) = &end_match_length($mvisual);
					my $mvisualr=reverse $mvisual;
					my ($end_match5) = &end_match_length($mvisualr);
					next if($end_match3 eq "Fail" || $end_match5 eq "Fail");#oligo end is not covered by template 
					if(defined $detail){
						print Detail "ntthal map:"; 
						print Detail join("\t", ($mvisual, $pos3, $pos5, $end_match3, $end_match5)),"\n";
					}
					next if(!defined $probe && !defined $revcom && $end_match3<0);
					if(defined $detail){
						print Detail "New info:",join("\t",$id, $strand, $chr, $pos, $oligo_seq, $tm, $end_match3,$mvisual),"\n";
					}
					print O join("\t",$id, $strand, $chr, $pos3, $oligo_seq, $tm, $end_match3,$mvisual),"\n";
					$bound{$id}{$strand."/".$chr."/".$pos3.":".$mvisual}=$tm;
					$bound_num++;
					if(!defined $NoFilter && $bound_num>$Max_Bound_Num){
						print F join("\t",  $id, $oligo_seq, "BoundTooMore", $bound_num."+"),"\n";
						delete $evalue{$id};
						$filter_by_bound=1;
						for(my $i=1; $i<@{$olen_oligo{$id0}}; $i++){
							my ($idn, $ori, $off, $pseqn)=@{$olen_oligo{$id0}->[$i]};
							print F join("\t",  $idn, $pseqn, "BoundTooMore", $bound_num."+"),"\n";
							delete $evalue{$idn};
						}
						last;
					}
					## other len's oligos and revcom
					for(my $i=1; $i<@{$olen_oligo{$id0}}; $i++){
						my ($idn, $ori, $off, $pseqn)=@{$olen_oligo{$id0}->[$i]};
						next if(!exists $evalue{$idn});
						my $posn = $pos3;
						my $end_matchn=$end_match3;
						my $seqn=$seq;
						my $strandn=$strand;
						my $mvisualn=$mvisual;
						if($ori eq "-"){
							$mvisualn=reverse($mvisual);
							$seqn=&revcom($seq);
							$posn = $pos5;
							$end_matchn=$end_match5;
							$strandn=$strand eq "+"? "-": "+";
						}
						my $tmn = `$ntthal -a ANY -s1 $pseqn -s2 $seqn -r`;
						chomp $tmn;
						next if($tmn<$min_tm_spec);
						my ($mvn, $ematchn) = &map_visual_trim($mvisualn, $off, $end_matchn); 
						next if(!defined $probe && $ematchn<0);
						if(defined $detail){
							print Detail "New Sam:",join("\t",$idn, $strandn, $chr, $posn, $pseqn, $tmn,$ematchn,$mvn),"\n";
						}
						print O join("\t",$idn, $strandn, $chr, $posn, $pseqn, $tmn,$ematchn,$mvn),"\n";
						$bound{$idn}{$strandn."/".$chr."/".$posn.":".$mvn}=$tmn;
					}
				}
				last if($filter_by_bound==1);
			}
		}
	}
	close(O);
	if(defined $detail){
		close (Detail);
	}
}

close(F);

#&SHSHOW_TIME("Output:");
### output
open (O, ">$outdir/$fkey.evaluation.out") or die $!;
if(!defined $nohead){
	print O "#ID\tSeq\tLen\tTm\tGC\tHairpin\tEND_Dimer\tANY_Dimer\tEndANum\tEndStability\tSNP\tPoly\tBoundNum\tHighestTm\tHighestInfo\n";
}

foreach my $id (sort {$a cmp $b} keys %evalue){
	## get bounds info of the max tm
	my $bnum=0;
	my ($abtms, $abinfos);
	if(exists $bound{$id}){## when tm too low, not exists
		($bnum, $abtms, $abinfos)=&get_highest_bound($bound{$id}, 3);
	}
	my ($btms, $binfos)=("NA", "NA"); ## all bound tm < min_sepc_tm
	if(defined $abtms){
		$btms = join(",", @{$abtms});
		$binfos = join(";", @{$abinfos});
	}
	print O join("\t",$id, @{$evalue{$id}}, $bnum, $btms, $binfos), "\n";
}
close(O);
#######################################################################################
print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################

# ------------------------------------------------------------------
# sub function
# ------------------------------------------------------------------
sub map_form_standard{
	my ($is_reverse, $cigar, $md)=@_;
	#($is_reverse, $cigar, $md)=(0, "16M4H", 16);
	if($cigar=~/[ID]/){ ## indel -> fail
		return "Fail";
	}
	my @mds=split /[ATCG]/, $md;
	if(scalar @mds >2){ ## mismatch>1 -> fail
		return "Fail";
	}
	my ($acigar_n, $acigar_str)=&cigar_split($cigar, "keepH");
	my @numcg=@{$acigar_n};
	my @strcg=@{$acigar_str};
	if(scalar @numcg>2){
		die "cigar: $cigar\n";
	}
	my @maps;
	for(my $i=0; $i<@strcg; $i++){
		if($strcg[$i] eq "M"){
			push @maps, $mds[0];
			if(scalar @mds==2){
				push @maps, ("N", $mds[1]);
			}
		}elsif($strcg[$i] eq "H"){
			push @maps, $numcg[$i]."H";
		}else{
			die "cigar: $cigar\n";
		}
	}
	my $form = join(",", @maps);
	if($is_reverse){
		$form=$maps[-1];
		for(my $i=$#maps-1; $i>=0; $i--){
			$form.=",".$maps[$i];
		}
	}
	return $form;
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


sub USAGE {#
	my $usage=<<"USAGE";
Program:
Version: $version
Contact:zeng huaping<huaping.zeng\@genetalks.com> 
	

Usage:
  Options:
  -p  <file>             Input oligo list file, forced
  -d  <files>            Input database files separated by "," to evalue specificity, [$fdatabases]
  -k  <str>              Key of output file, forced

  --Revcom               Also evalue revcom oligos
  --Probe                Design probe and will consider mapping region where oligo 3end not matched exactly when caculate specificity
  -olen   <int,int,int>  Evalue other length's oligos, <min,max,scale> of length, optional
  -opttm     <int>       optimal tm of oligo, [$opt_tm]
  -opttmp    <int>       optimal tm of probe, [$opt_tm_probe]
  -maplen    <int>      length to map with bwa, [$len_map]
  -stm       <int>      min tm to be High_TM when caculate specificity, [$min_tm_spec]
  -mv        <int>      concentration of monovalent cations in mM, [$mv]
  -dv        <float>    concentration of divalent cations in mM, [$dv]
  -dNTP      <float>    concentration of deoxynycleotide triphosphate in mM, [$dNTP]
  -dna       <int>      concentration of DNA strands in nM, [$dna]
  -tp        [0|1]      Specifies the table of thermodynamic parameters and the method of melting temperature calculation,[$tp]
                        0   Breslauer et al., 1986 and Rychlik et al., 1990
						1   Use nearest neighbor parameters from SantaLucia 1998
  -sc        [0|1|2]    Specifies salt correction formula for the melting temperature calculation, [$sc]
                        0   Schildkraut and Lifson 1965, used by oligo3 up to and including release 1.1.0.
						1   SantaLucia 1998
						2   Owczarzy et al., 2004
  -thread    <int>      thread in bwa, [$thread]
  -maxtime  <int>       max bwa running time, killed when time out, not killed please set 10000000000, [$max_time]

  --NoFilter             Not filter any oligos
  --NoSpecificity        Not evalue specificity
  --Detail              Output Detail Info to xxx.evaluation.detail, optional
  -od        <dir>      Dir of output file, default ./
  -h		 Help

USAGE
	print $usage;
	exit;
}
