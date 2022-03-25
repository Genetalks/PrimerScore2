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
require "$Bin/io.pm";

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
my ($Methylation,$NoSpecificity,$NoFilter, $Precise);
my $min_tm_spec = 45; #when caculate specificity
my $nohead;
my $thread = 10;
my $fdatabases = $REF_GRCh37;
my $len_map=20; ##bwa result is the most when 20bp
my $opt_tm = 60;
my $opt_tm_probe=70;
my $opt_size = 100;
my $Debug;
my $probe;
my $ftemplate;
my ($mv, $dv, $dNTP, $dna, $tp, $sc)=(50, 1.5, 0.6, 50, 1, 1);
my $olens;
my $revcom;
my $sublen = 8; ## substr end3's seq to detect dimer, because primer3 always don't predict dimers with lowtm although end3 is matched exactly
my $max_time;
GetOptions(
				"help|?" =>\&USAGE,
				"p:s"=>\$foligo,
				"t:s"=>\$ftemplate,
				"d:s"=>\$fdatabases,
				"k:s"=>\$fkey,
				"Revcom:s"=>\$revcom,
				"Methylation:s"=>\$Methylation,
				"Probe:s"=>\$probe,
				"NoFilter:s"=>\$NoFilter,
				"Precise:s"=>\$Precise,
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
my $MinCpG = 1;
my $MinC = 3;
my @tm = ($opt_tm*0.8, $opt_tm*2, $opt_tm*0.6, $opt_tm*2);
	
my $oligotm = "$PATH_PRIMER3/src/oligotm -mv $mv -dv $dv -n $dNTP -d $dna -tp $tp -sc $sc";
my $ntthal = "$PATH_PRIMER3/src/ntthal -mv $mv -dv $dv -n $dNTP -d $dna";

if(defined $Precise){
	$max_time=defined $max_time? $max_time: 3000;
}else{
	$max_time=defined $max_time? $max_time: 120;
}

## check
if(!defined $NoSpecificity && !defined $ftemplate){
	print "Wrong: a template file(-t) is required when not choose --NoSpecificity! \n";
	exit(-1);
}
open(F, ">$outdir/$fkey.filter.list") or die $!;
open(L, ">$outdir/$fkey.oligo.list") or die $!;
my %evalue;
my %map;
my %map_id;
my %olen_oligo;
my %primer3_tm;
open (P, $foligo) or die $!;
#&SHSHOW_TIME("Routine analysis:");
while (<P>){
	chomp;
	next if(/^$/);
	$_=~s/\s+$//;
	my ($id0, $seq)=split /\s+/, $_;
	my ($oligo_seq0, $oligo_seq_snp0, $oligo_seq_mark0) = split /,/, $seq;
	$oligo_seq0=uc $oligo_seq0;
	if(!defined $oligo_seq_snp0 || $oligo_seq_snp0 eq "NA"){
		$oligo_seq_snp0=$oligo_seq0;
	}
	$oligo_seq_snp0=uc $oligo_seq_snp0;
	if(!defined $oligo_seq_mark0 || $oligo_seq_mark0 eq "NA"){
		$oligo_seq_mark0="|" x length($oligo_seq0);
	}
	if(!defined $oligo_seq0){
		print $_,"\n";
		die;
	}
	my $id0new=defined $olens? $id0."_F_0": $id0;
	push @{$olen_oligo{$id0}}, [$id0new,"+", 0, $oligo_seq0, $oligo_seq_snp0, $oligo_seq_mark0];
	my ($oligo_seq0r, $oligo_seq_snp0r, $oligo_seq_mark0r);
	if(defined $revcom){
		$oligo_seq0r=&revcom($oligo_seq0);
		$oligo_seq_snp0r=&revcom($oligo_seq_snp0);
		$oligo_seq_mark0r=reverse($oligo_seq_mark0);
		my $id0newR=defined $olens? $id0."_R_0": $id0;
		push @{$olen_oligo{$id0}}, [$id0newR, "-", 0, $oligo_seq0r, $oligo_seq_snp0r, $oligo_seq_mark0r];
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
			my $seq_mark = substr($oligo_seq_mark0, $off);
			push @{$olen_oligo{$id0}}, [$id, "+", $off, $seq, $seq_snp, $seq_mark];
			if(defined $revcom){
				$id=$id0."_R"."_".$off;
				my $seq = substr($oligo_seq0r, $off);
				my $seq_snp = substr($oligo_seq_snp0r, $off);
				my $seq_mark = substr($oligo_seq_mark0r, $off);
				push @{$olen_oligo{$id0}}, [$id,"-", $off, $seq, $seq_snp, $seq_mark];
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
		my ($id, $ori, $off, $oligo_seq, $oligo_seq_snp, $oligo_seq_mark)=@{$olen_oligo{$id0}->[$i]};
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
		$primer3_tm{$oligo_seq}=$Tm;
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
		
		## filter Self_Complementarity: Hairpin
		my $hairpin_tm = `$ntthal -a HAIRPIN -s1 $oligo_seq -r`;
		chomp $hairpin_tm;
		if(!defined $NoFilter && $hairpin_tm > $MAX_hairpin_tm){
			print F join("\t", $id, $oligo_seq, "Hairpin", $hairpin_tm),"\n";
			next;
		}

		## dimer check
		my ($is_amplify, @dtype, @eff, @dlen);
		{
		my $info = `$ntthal -a END1 -s1 $oligo_seq -s2 $oligo_seq`;
		chomp $info;
		my ($dtm, $end31, $end32, $amplen, $mlen3, $dlen, $msum, $indel)=&dimer_amplify($info);
		my ($is_amplify, $dtype, $eff)=&judge_amplify($dtm, $end31, $end32, $amplen, $mlen3, $msum, $indel); 
		if($eff>0){
			push @dtype, $dtype;
			push @eff, $eff;
			push @dlen, $dlen;
		}
		}
		
		if(scalar @dtype==0 || $dtype[0] ne "AmpEndMeet"){
		my $subseq = substr($oligo_seq, length($oligo_seq)-$sublen, $sublen);
		my $info = `$ntthal -a END1 -s1 $subseq -s2 $subseq`;
		chomp $info;
		my ($dtm, $end31, $end32, $amplen, $mlen3, $dlen, $msum, $indel)=&dimer_amplify($info);
		$dlen+=length($oligo_seq)-$sublen+length($oligo_seq)-$sublen;
		my ($is_amplify, $dtype, $eff)=&judge_amplify_endmeet($dtm, $end31, $end32, $mlen3);
		if($eff>0){
			push @dtype, $dtype;
			push @eff, $eff;
			push @dlen, $dlen;
		}
		}
		my $dimertype="NA";
		my $dimersize="NA";
		if(scalar @dtype>0){
			$dimertype=join(",", @dtype);
			$dimersize=join(",", @dlen);

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

		## Methylation oligos filter
		my ($CpGs, $nonCpG_Cs);
		if(defined $Methylation){
			($CpGs, $nonCpG_Cs)=&CpG_info($oligo_seq_mark);
			my @cpg=split /,/, $CpGs;
			my @cs=split /,/, $nonCpG_Cs;
			if(!defined $NoFilter){
				if(scalar @cpg < $MinCpG){
					print F join("\t", $id, $oligo_seq, "CpG",$CpGs),"\n";
					next;
				}
				if(scalar @cs < $MinC){
					print F join("\t", $id, $oligo_seq, "nonCpG_C",$nonCpG_Cs),"\n";
					next;
				}
			}
		}
		
		push @{$evalue{$id}},($oligo_seq, $len, $Tm, sprintf("%.3f",$GC), sprintf("%.2f",$hairpin_tm), $dimertype, $dimersize,$nendA,$dG_end3, $SNP_info, $poly_info);
		if(defined $Methylation){
			push @{$evalue{$id}}, ($CpGs, $nonCpG_Cs);
		}
		$is_all_filter=0;
	}
}

exit(0) if(scalar keys %evalue==0);
foreach my $id(keys %evalue){
	print L "Final: $id\n";
}
close(L);

### output
open (O, ">$outdir/$fkey.evaluation.out") or die $!;
if(!defined $nohead){
	print O &evaluation_head($Methylation, $NoSpecificity), "\n";
}

foreach my $id (sort {$a cmp $b} keys %evalue){
	print O join("\t",$id, @{$evalue{$id}});
	print O "\n";
}
close(O);


### Specificity, get bound info
if(!defined $NoSpecificity){
	&Run("perl $Bin/get_bound_info.pl -tm $min_tm_spec -io $outdir/$fkey.evaluation.out -it $ftemplate -id $fdatabases -k $fkey -od $outdir -t $thread");
}
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
	my ($is_reverse, $cigar, $md)=@_;
	my $H3;
	if($is_reverse==0){
		($H3)=$cigar=~/M(\d+)H$/;
	}else{
		($H3)=$cigar=~/^(\d+)H/;
	}
	if(!defined $H3){ ## 1H can be prc from experiment data
		if(($is_reverse==0 && $md=~/0[ATCG]0$/) || ($is_reverse==1 && $md=~/^0[ATCG]0/)){
			$H3=2;
		}elsif(($is_reverse==0 && $md=~/[ATCG]0$/) || ($is_reverse==1 && $md=~/^0[ATCG]/)){
			$H3=1;
		}else{
			$H3=0;
		}
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
  -t  <file>             Input template file to evalue specificity, not required when --NoSpecificity, optional
  -d  <files>            Input database files separated by "," to evalue specificity, [$fdatabases]
  -k  <str>              Key of output file, forced

  --Revcom               Also evalue revcom oligos
  --Methylation          Design methylation oligos
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
  -thread    <int>      thread in blast, [$thread]
  -maxtime  <int>       max bwa running time, killed and filtered when time out, [120]

  --NoFilter             Not filter any oligos
  --NoSpecificity        Not evalue specificity
  --Precise              Evalue specificity precisely, but will consume a long time, -maxtime sets 3000
  --Detail              Output Detail Info to xxx.evaluation.detail, optional
  -od        <dir>      Dir of output file, default ./
  -h		 Help

USAGE
	print $usage;
	exit;
}
