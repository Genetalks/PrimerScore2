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
my $min_tm_spec = 40; #when caculate specificity
my $nohead;
my $thread = 3;
our $PATH_PRIMER3;
our $REF_HG19;
our $SAMTOOLS;
our $BWA;
my $fdatabases = $REF_HG19;
my $NoFilter;
my $len_map=20; ##bwa result is the most when 20bp
my $opt_tm = 60;
my $opt_tm_probe=70;
my $opt_size = 100;
my $Debug;
my $probe;
my ($mv, $dv, $dNTP, $dna, $tp, $sc)=(50, 1.5, 0.6, 50, 1, 1);
my $olens;
my $revcom;
GetOptions(
				"help|?" =>\&USAGE,
				"p:s"=>\$fprimer,
				"d:s"=>\$fdatabases,
				"k:s"=>\$fkey,
				"Revcom:s"=>\$revcom,
				"Probe:s"=>\$probe,
				"NoFilter:s"=>\$NoFilter,
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
&USAGE unless ($fprimer and $fkey);

$outdir||="./";
`mkdir $outdir`	unless (-d $outdir);
$outdir=AbsolutePath("dir",$outdir);
my $merge_len = 100;
my $Wind_GC = 8;
my $MAX_hairpin_tm = 55; ## 55 from experience
my $Max_Dimer_tm = 55;
my $Max_SNP_num = 1;
my $Max_poly_total = 12;
my $Max_poly_G = 3;
my $Max_poly_ATC = 5; 
my $MAX_endA = 4;
my $MAX_poly = 20;
my $Max_Bound_Num = 300; ## too high tm aligns will be filtered
my $MIN_tm = $opt_tm-5;
my $MAX_tm = defined $opt_tm_probe? $opt_tm_probe+5: $opt_tm+5;
my $MIN_gc = 0.2;
my $MAX_gc = 0.8;
my @tm = ($opt_tm*0.8, $opt_tm*2, $opt_tm*0.6, $opt_tm*2);
	
my $oligotm = "$PATH_PRIMER3/src/oligotm";
my $ntthal = "$PATH_PRIMER3/src/ntthal";

## creat primer.fa
open(PN, ">$outdir/$fkey.primer.fa") or die $!;
if(!defined $NoFilter){
	open(F, ">$outdir/$fkey.filter.list") or die $!;
}
my %evalue;
my %map;
my %map_id;
my %olen_primer;
open (P, $fprimer) or die $!;
#&SHSHOW_TIME("Routine analysis:");
while (<P>){
	chomp;
	next if(/^$/);
	$_=~s/\s+$//;
	my ($id0, $seq)=split /\s+/, $_;
	my ($primer_seq0, $primer_seq_snp0) = split /:/, $seq;
	if(!defined $primer_seq_snp0){
		$primer_seq_snp0=$primer_seq0;
	}
	if(!defined $primer_seq0){
		print $_,"\n";
		die;
	}
	push @{$olen_primer{$id0}}, [$id0."_F_0","+", 0, $primer_seq0, $primer_seq_snp0];
	my ($primer_seq0r, $primer_seq_snp0r);
	if(defined $revcom){
		$primer_seq0r=&revcom($primer_seq0);
		$primer_seq_snp0r=&revcom($primer_seq_snp0);
		push @{$olen_primer{$id0}}, [$id0."_R_0", "-", 0, $primer_seq0r, $primer_seq_snp0r];
	}
	if(defined $olens){
		my ($min, $max, $scale)=split /,/, $olens;
		my $len0 = length($primer_seq0);
		for(my $l=$min; $l<$max; $l+=$scale){
			last if($l>=$len0);
			my $off=$len0-$l;
			my $id=$id0."_F"."_".$off;
			my $seq = substr($primer_seq0, $off);
			my $seq_snp = substr($primer_seq_snp0, $off);
			push @{$olen_primer{$id0}}, [$id, "+", $off, $seq, $seq_snp];
			if(defined $revcom){
				$id=$id0."_R"."_".$off;
				my $seq = substr($primer_seq0r, $off);
				my $seq_snp = substr($primer_seq_snp0r, $off);
				push @{$olen_primer{$id0}}, [$id,"-", $off, $seq, $seq_snp];
			}
		}
	}

	## filter
	my $nendA0 = &get_end_A($primer_seq0);
	my $nendA0r = &get_end_A($primer_seq0r);
	my $is_all_filter=1;
	for(my $i=0; $i<@{$olen_primer{$id0}}; $i++){
		my ($id, $ori, $off, $primer_seq, $primer_seq_snp)=@{$olen_primer{$id0}->[$i]};
		my $ftype;
		## filter endA
		my $nendA = $ori eq "+"? $nendA0: $nendA0r;
		if(!defined $NoFilter && $nendA>$MAX_endA){
			$ftype = "EndA";
			print F join("\t",  $id, $primer_seq, $ftype, $nendA),"\n";
			next;
		}
		
		## filter TM and GC
		my $Tm = `$oligotm -mv $mv -dv $dv -n $dNTP -d $dna -tp $tp -sc $sc $primer_seq `;
		chomp $Tm;
		$Tm = sprintf "%0.2f", $Tm;
		my $GC=&GC($primer_seq);
		if(!defined $NoFilter && ($Tm<$MIN_tm || $Tm>$MAX_tm)){
			$ftype = "Tm";
			print F join("\t",  $id, $primer_seq, $ftype, $Tm),"\n";
			next;
		}
		if(!defined $NoFilter && ($GC<$MIN_gc || $GC>$MAX_gc)){
			$ftype = "GC";
			print F join("\t",  $id, $primer_seq, $ftype, $GC),"\n";
			next;
		}
		
		## filter Self_Complementarity
		my $hairpin_tm = `$ntthal -a HAIRPIN -s1 $primer_seq -r`;
		chomp $hairpin_tm;
		if(!defined $NoFilter && $hairpin_tm > $MAX_hairpin_tm){
			print F join("\t", $id, $primer_seq, "Hairpin:".$hairpin_tm),"\n";
			next;
		}
		my $END_tm = `$ntthal -a END1 -s1 $primer_seq -s2 $primer_seq -r`;
		chomp $END_tm;
		if(!defined $NoFilter && $END_tm > $Max_Dimer_tm){
			print F join("\t", $id, $primer_seq, "END Dimer:".$END_tm),"\n";
			next;
		}
		my $ANY_tm = `$ntthal -a ANY -s1 $primer_seq -s2 $primer_seq -r`;
		chomp $ANY_tm;
		if(!defined $NoFilter && $ANY_tm > $Max_Dimer_tm){
			print F join("\t", $id, $primer_seq, "ANY Dimer:".$ANY_tm),"\n";
			next;
		}
		
		## filter SNP
		my ($SNP_num, $SNP_info)=&SNP_check($primer_seq_snp);
		$SNP_info=$SNP_num==0? "NA":$SNP_info;
#		$SNP_info.=",".$primer_seq_snp;
		if(!defined $NoFilter && $SNP_num > $Max_SNP_num){
			print F join("\t", $id, $primer_seq, "SNP:".$SNP_info),"\n";
			next;
		}
	
		## filter poly
		my ($total, $max_len, $max_base, $poly_info)=&poly_check($primer_seq);
		$poly_info = "NA" if($total==0);
		if(!defined $NoFilter && $total>0){
			if($total > $Max_poly_total || ($max_base eq "G" && $max_len>$Max_poly_G) || ($max_base ne "G" && $max_len>$Max_poly_ATC)){
				print F join("\t", $id, $primer_seq, "Poly:".$poly_info),"\n";
				next;
			}
		}
		my $len = length($primer_seq);
		push @{$evalue{$id}},($primer_seq, $len, $Tm, sprintf("%.3f",$GC), sprintf("%.2f",$hairpin_tm), sprintf("%.2f",$END_tm), sprintf("%.2f",$ANY_tm), $SNP_info, $poly_info);
		$is_all_filter=0;
	}
	
	
	if($is_all_filter==0){
		my $off = (length $primer_seq0)-$len_map;
		$off=$off>0? $off: 0;
		my $mseq = substr($primer_seq0, $off);
		print PN ">$id0\n";
		print PN $mseq,"\n";
	}
	
}
close(PN);
exit(0) if(scalar keys %evalue==0);

#&SHSHOW_TIME("Specificity analysis:");

if(defined $detail){
	open (Detail, ">$outdir/$fkey.evaluation.detail") or die $!;
}
my $DB;
my %bound;
my %bound_stat;
if(!defined $NoSpecificity){
	my %db_region;
	### bwa
	my %mapping;
	my $fa_primer = "$outdir/$fkey.primer.fa";
	my @fdatabase=split /,/, $fdatabases;	
	foreach my $fdatabase(@fdatabase){
		if(!-e "$fdatabase\.ann"){
			`bwa index $fdatabase`;
		}
		my $dname = basename($fdatabase);
		Run("$BWA mem -D 0 -k 9 -t $thread -c 5000000 -y 100000 -T 12 -B 1 -L 2,2 -h 200 -a  $fdatabase $fa_primer |samtools view -bS - >$fa_primer\_$dname.bam");
	#	Run("bwa mem -D 0 -k 9 -t 4 -c 5000000 -y 100000 -T 12 -B 1 -O 2,2 -L 1,1 -h 200 -a  $fdatabase $fa_primer > $fa_primer.sam");
		### read in sam
		open (I, "samtools view $fa_primer\_$dname.bam|") or die $!;
		while (<I>){
			chomp;
			my ($id, $flag, $chr, $pos, $score, $cigar, undef, undef, undef, $seq)=split /\s+/,$_;
			my ($is_unmap, $is_reverse)=&explain_bam_flag_unmap($flag);
			my ($md)=$_=~/MD:Z:(\S+)/;
			next if ($is_unmap);
			my ($H3)=&get_3end1_mismatch($is_reverse, $cigar);
			#print "H3:", join("\t", $H3,$id, $is_reverse, $flag, $chr, $pos, $score, $cigar, $md),"\n";
			next if(!defined $probe && !defined $revcom && $H3>1);
			push @{$mapping{$id}{$dname}},[$is_reverse, $flag, $chr, $pos, $score, $cigar, $md, $fdatabase];
		}
		close(I);
	}
	### evaluate
	open(O, ">$outdir/$fkey.specificity.sam") or die $!;
	foreach my $id0 (sort {$a cmp $b} keys %olen_primer){
		my ($id,undef, undef, $primer_seq)=@{$olen_primer{$id0}->[0]};
		my $bound_num = 0;
		my %lowtm;
		foreach my $dname(keys %{$mapping{$id0}}){
			my $map_num = scalar @{$mapping{$id0}{$dname}};
			for (my $i=0; $i<$map_num; $i++){
				my ($is_reverse, $flag, $chr, $pos, $score, $cigar, $md, $fdatabase)=@{$mapping{$id0}{$dname}->[$i]};
				my $strand=$is_reverse? "-": "+";
				next if(exists $lowtm{join(",",$is_reverse,$cigar, $md)});## filter bound regions whose map info are same with where bound tm too low 
				if(defined $detail){
					print Detail "\nOriginal:",join("\t",$id0, $is_reverse, $flag, $chr, $pos, $score),"\n";
				}
				######## specificity
				my $len = length($primer_seq);
				my $extend = 20;
				
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
				if($seq!~/[ATCGatcg]+/){
					print "Wrong: extract aligned region sequence failed\n";
					print $_,"\n";
					print join("\t", $chr, $start, $end, $fdatabase),"\n";
					die;
				}
				$seq=uc($seq);
				$seq=&revcom($seq)if(!$is_reverse);
				### ntthal
				my $result = `$ntthal -a ANY -s1 $primer_seq -s2 $seq`;
				if(defined $detail){
					print Detail "$ntthal -a ANY -s1 $primer_seq -s2 $seq\n";
					print Detail $result;
				}
			
				### get tm 
				my ($dG, $tm) = $result =~/dG = ([\d\+\-\.]+)\tt = ([\d\+\-\.]+)/;
				next if($dG eq ""); ## map to NNNN region, invalid
				if($tm<$min_tm_spec){
					$lowtm{join(",",$is_reverse,$cigar, $md)}=1;
					next;
				}
				my @line = split /\n/, $result;
				shift @line;	
				if(!defined $line[1]){ ## just print and return
					print join("\t",$primer_seq, $chr, $pos, $is_reverse, $fdatabase),"\n";
					print "$SAMTOOLS faidx $fdatabase $chr:$start-$end\n";
					print "$ntthal -a ANY -s1 $primer_seq -s2 $seq\n";
					return (-1,-1,-1,-1);
				}
				## match visual
				my ($mvisual, $pos3, $pos5)=&map_visual_from_ntthal(\@line, $is_reverse, $start, $end);
				my ($end_match3) = &end_match_length($mvisual);
				my $mvisualr=reverse $mvisual;
				my ($end_match5) = &end_match_length($mvisualr);
				if(defined $detail){
					print Detail "ntthal map:"; 
					print Detail join("\t", ($mvisual, $pos3, $pos5, $end_match3, $end_match5)),"\n";
				}
				next if(!defined $probe && !defined $revcom && $end_match3<0);
				$bound_num++;
				if(!defined $NoFilter && $bound_num>$Max_Bound_Num){
					for(my $x=0; $x<@{$olen_primer{$id0}}; $x++){
						my ($idt, undef, undef, $seqt) = @{$olen_primer{$id0}->[$x]};
						##here, Bound contain regions where end_match3<0 when defined $probe or $revcom; it maybe imprecise and will filter some oligos fit for primer, but the effect should be very small when parameter $Max_Bound_Num is very big
						print F join("\t", $idt, $seqt, "Bound_Too_More"),"\n";
						delete $evalue{$idt};
					}
					last;
				}
				if(defined $detail){
					print Detail "New info:",join("\t",$id, $strand, $chr, $pos, $primer_seq, $tm, $end_match3,$mvisual),"\n";
				}
				if(defined $probe || (!defined $probe && $end_match3>0)){
					$bound_stat{$id}++;
					print O join("\t",$id, $strand, $chr, $pos3, $primer_seq, $tm, $end_match3,$mvisual),"\n";
					push @{$bound{$id}{$tm}}, [$strand, $chr, $pos3, $end_match3,$mvisual];
				}
				## other len's primers and revcom
				for(my $i=1; $i<@{$olen_primer{$id0}}; $i++){
					my ($idn, $ori, $off, $pseqn)=@{$olen_primer{$id0}->[$i]};
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
					$bound_stat{$idn}++;
					print O join("\t",$idn, $strandn, $chr, $posn, $pseqn, $tmn,$ematchn,$mvn),"\n";
					push @{$bound{$idn}{$tmn}}, [$strandn, $chr, $posn, $ematchn,$mvn];
				}
			}
		}
	}
	close(O);
	if(defined $detail){
		close (Detail);
	}
}

if(!defined $NoFilter){
	close(F);
}

#&SHSHOW_TIME("Output:");
### output
open (O, ">$outdir/$fkey.evaluation.out") or die $!;
if(!defined $nohead){
	print O "#ID\tSeq\tLen\tTm\tGC\tHairpin\tEND_Dimer\tANY_Dimer\tSNP\tPoly\tBoundNum\tHighestTm\tHighestInfo\n";
}

foreach my $id (sort {$a cmp $b} keys %evalue){
	## get bounds info of the max tm
	my @tms=sort{$b <=> $a} keys %{$bound{$id}};
	my $maxn=3;
	my @tm;
	my @binfo;
	my $n=0;
	for(my $i=0; $i<@tms; $i++){
		for(my $j=0; $j<@{$bound{$id}{$tms[$i]}}; $j++){
			push @tm, sprintf("%.2f",$tms[$i]);
			my ($strand, $chr, $pos, $ematchn,$mvn)=@{$bound{$id}{$tms[$i]}->[$j]};
#			push @binfo, join("/", $strand, $chr, $pos, $mvn);
			push @binfo, $strand."/".$chr."/".$pos.":".$mvn;
			$n++;
			last if($n>=$maxn);
		}
		last if($n>=$maxn);
	}
	my $bnum=exists $bound_stat{$id}? $bound_stat{$id}: 0;
	print O join("\t",$id, @{$evalue{$id}}, $bnum, join(",", @tm), join(";", @binfo)), "\n";
}
close(O);
#######################################################################################
print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################

# ------------------------------------------------------------------
# sub function
# ------------------------------------------------------------------
sub SNP_check{
	my ($seq)=@_;
	$seq = reverse $seq;
	my @u=split //, $seq;
	my @snp;
	for(my $i=0; $i<@u; $i++){
		if($u[$i] eq "S"){
			push @snp, $i."S1";
		}elsif($u[$i] eq "D"){
			my $l=1;
			while(1){
				if($u[$i+1] eq "D"){
					$l++;
					$i++;
				}else{
					last;
				}
			}
			push @snp, $i."D".$l;
		}elsif($u[$i]=~/IJKL/){
			my $l=ord($u[$i])-ord('I')+1;
			push @snp, $i."I".$l;
		}
	}
	return (scalar @snp, join(",", @snp));
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
  -p  <file>             Input primer list file, forced
  -d  <files>            Input database files separated by ",", [$fdatabases]
  -k  <str>              Key of output file, forced

  --Revcom               evalue revcom primers
  --Probe                design probe and will consider mapping region where primer 3end not matched exactly when caculate specificity
  -olen   <int,int,int>  evalue other length's primers, <min,max,scale> of length, optional
  -opttm     <int>       optimal tm of primer, [$opt_tm]
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
                        0   Schildkraut and Lifson 1965, used by primer3 up to and including release 1.1.0.
						1   SantaLucia 1998
						2   Owczarzy et al., 2004
  -thread    <int>      thread in bwa, [$thread]
  --NoFilter             Not filter any primers
  --NoSpecificity        Not evalue specificity
  --Detail              Output Detail Info to xxx.evaluation.detail, optional
  -od        <dir>      Dir of output file, default ./
  -h		 Help

USAGE
	print $usage;
	exit;
}
