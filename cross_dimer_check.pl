#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
require "$Bin/path.pm";
require "$Bin/dimer.pm";
my $BEGIN_TIME=time();
my $version="1.0.0";
#######################################################################################

# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
my ($fIn,$fkey,$outdir);
my $type = "Result";
my ($mv, $dv, $dNTP, $dna, $tp, $sc)=(50, 1.5, 0.6, 50, 1, 1);
my $MultiPlex;
my $SelfComplementary;
my $high_tm = 45;
my $min_tm_omega = 17;
my $min_tm_amp = 29; ##Dis-7:26
my $min_meetlen=3; ## min end3 match len when Enddimer
my $min_amplen=15;
my $sublen = 8; ## substr end3's seq to detect dimer, because primer3 always don't predict dimers with lowtm although end3 is matched exactly
my $adapter1="AGATGTGTATAAGAGACAG";
#my $adapter2="GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT";
my $adapter2="CTGGAGTTCAGACGTGTGCTCTTCCGATCT"; ## remove 4 bp, for too long to ntthal
my $Nofilter;
GetOptions(
				"help|?" =>\&USAGE,
				"i:s"=>\$fIn,
				"t:s"=>\$type,
				"k:s"=>\$fkey,
				"mv:s"=>\$mv,
				"dv:s"=>\$dv,
				"dNTP:s"=>\$dNTP,
				"dna:s"=>\$dna,
				"adpt1:s"=>\$adapter1,
				"adpt2:s"=>\$adapter2,
				"MultiPlex:s"=>\$MultiPlex,
				"Nofilter:s"=>\$Nofilter,
				"SelfComplementary:s"=>\$SelfComplementary,
				"od:s"=>\$outdir,
				) or &USAGE;
&USAGE unless ($fIn and $fkey);

$outdir||="./";
`mkdir $outdir`	unless (-d $outdir);
$outdir=AbsolutePath("dir",$outdir);

our $PATH_PRIMER3;
my $ntthal = "$PATH_PRIMER3/src/ntthal -mv $mv -dv $dv -n $dNTP -d $dna";

my ($ix1, $ix2);
if($type eq "Result"){
	($ix1, $ix2)=(3,4);
}elsif($type eq "Primer"){
	($ix1, $ix2)=(0,1);
}

my %seq;
my %id;
open(I, $fIn) or die $!;
while(<I>){
	chomp;
	next if(/^$/ || /^#/);
	my ($id, $seq)=(split /\s+/, $_)[$ix1, $ix2];
	my ($tid, $LR)=$id=~/(\S+)-([12LRP])/;
	$LR=~s/L/1/;
	$LR=~s/R/2/;
	$LR=~s/P/3/;
	my $subseq = substr($seq, length($seq)-$sublen, $sublen);
	if($LR==1){
		$seq=$adapter1.$seq;
	}elsif($LR==2){
		$seq=$adapter2.$seq;
	}
	@{$seq{$tid}->[$LR-1]}=($seq, $subseq);
	$id{$tid}[$LR-1]=$id;
}
close(I);


open(O,">$outdir/$fkey.cross.dimer")or die $!;
foreach my $tid(sort {$a cmp $b} keys %seq){
	my @seq=@{$seq{$tid}};
	for(my $j=0; $j<@seq; $j++){
		my ($oligo_seq, $subseq) = @{$seq{$tid}->[$j]};
		my $id = $id{$tid}->[$j];
#		print O ">",$id,"\t",$oligo_seq,"\n";
		if(defined $SelfComplementary){
			my $HPinfo = `$ntthal -a HAIRPIN -s1 $oligo_seq`;
			chomp $HPinfo;
			if(!defined $HPinfo){
				print $HPinfo;
				my ($HPtm)=$HPinfo=~/t = (\S+)/;
				if($HPtm >= $high_tm){
					print O "#Hairpin\t$id\t$HPtm\n"; 
					print O $HPinfo, "\n";
				}
			}
			my ($is_amplify, $atype, $eff);
			{
			my $info = `$ntthal -a END1 -s1 $oligo_seq -s2 $oligo_seq`;
			chomp $info;
			my ($tm, $end31, $end32, $amplen, $mlen3, $len, $msum, $indel)=&dimer_amplify($info);
			($is_amplify, $atype, $eff)=&judge_amplify($tm, $end31, $end32, $amplen, $mlen3, $msum, $indel, $min_tm_omega, $min_tm_amp, $min_amplen, $min_meetlen); 
			if($is_amplify || defined $Nofilter){
				print O join("\t", "#Dimer", $id, $id, $atype, $eff, $len, $tm, $end31, $end32, $amplen, $mlen3, $msum),"\n";
				print O "$ntthal -a END1 -s1 $oligo_seq -s2 $oligo_seq\n";
				print O $info, "\n";
			}
			}

			if($atype ne "AmpEndMeet"){
			my $info = `$ntthal -a END1 -s1 $subseq -s2 $subseq`;
			chomp $info;
			my ($tm, $end31, $end32, $amplen, $mlen3, $len, $msum, $indel)=&dimer_amplify($info);
			next if($msum>1 && $atype eq "AmpEndMap");
			$len+=length($oligo_seq)-$sublen+length($oligo_seq)-$sublen;
			($is_amplify, $atype, $eff)=&judge_amplify_endmeet($tm, $end31, $end32, $mlen3, $min_meetlen);
			if($is_amplify || defined $Nofilter){
				$amplen="NA";
				print O join("\t", "#EndDimer", $id, $id, $atype, $eff, $len, $tm, $end31, $end32, $amplen, $mlen3, $msum),"\n";
				print O "$ntthal -a END1 -s1 $subseq -s2 $subseq\n";
				print O $info, "\n";
			}
			}

			
		}
		## cross primer of same tid
		for(my $k=$j+1; $k<@seq; $k++){
			my ($info, $tm, $umap1, $umap2);
			my ($seq, $subs) = @{$seq{$tid}->[$k]};
			my $id2=$id{$tid}->[$k];
			my ($is_amplify, $atype, $eff);
			{
			my $info = `$ntthal -a END1 -s1 $oligo_seq -s2 $seq`;
			chomp $info;
			my ($tm, $end31, $end32, $amplen, $mlen3, $len, $msum, $indel)=&dimer_amplify($info);
			($is_amplify, $atype, $eff)=&judge_amplify($tm, $end31, $end32, $amplen, $mlen3, $msum, $indel, $min_tm_omega, $min_tm_amp, $min_amplen, $min_meetlen); 
			if($is_amplify || defined $Nofilter){
				print O join("\t", "#Dimer", $id, $id2, $atype, $eff, $len, $tm, $end31, $end32, $amplen, $mlen3, $msum),"\n";
				print O "$ntthal -a END1 -s1 $oligo_seq -s2 $seq\n";
				print O $info, "\n";
			}
			}
			if(!$is_amplify || defined $Nofilter){
			my $info = `$ntthal -a END2 -s1 $oligo_seq -s2 $seq`;
			chomp $info;
			my ($tm, $end31, $end32, $amplen, $mlen3, $len, $msum, $indel)=&dimer_amplify($info);
			($is_amplify, $atype, $eff)=&judge_amplify($tm, $end31, $end32, $amplen, $mlen3, $msum, $indel, $min_tm_omega, $min_tm_amp, $min_amplen, $min_meetlen); 
			if($is_amplify || defined $Nofilter){
				print O join("\t", "#Dimer", $id, $id2, $atype, $eff, $len, $tm, $end31, $end32, $amplen, $mlen3, $msum),"\n";
				print O "$ntthal -a END2 -s1 $oligo_seq -s2 $seq\n";
				print O $info, "\n";
			}
			}
			
			if($atype ne "AmpEndMeet"){
			my $info = `$ntthal -a END1 -s1 $subseq -s2 $subs`;
			chomp $info;
			my ($tm, $end31, $end32, $amplen, $mlen3, $len, $msum, $indel)=&dimer_amplify($info);
			next if($msum>1 && $atype eq "AmpEndMap");
			$len+=length($seq)-$sublen+length($oligo_seq)-$sublen;
			
			($is_amplify, $atype, $eff)=&judge_amplify_endmeet($tm, $end31, $end32, $mlen3, $min_meetlen);
			if($is_amplify || defined $Nofilter){
				$amplen="NA";
				print O join("\t", "#Dimer", $id, $id2, $atype, $eff, $len, $tm, $end31, $end32, $amplen, $mlen3, $msum),"\n";
				print O "$ntthal -a END1 -s1 $subseq -s2 $subs\n";
				print O $info, "\n";
			}
			}
		}


		## cross tid 
		if(defined $MultiPlex){
			foreach my $tid2(keys %seq){
				next if($tid2 eq $tid);
				my @seq2=@{$seq{$tid2}};
				for(my $k=0; $k<@seq2; $k++){
					my ($seq, $subs) = @{$seq{$tid2}->[$k]};
					my $id2=$id{$tid2}->[$k];
					my ($is_amplify, $atype, $eff);
					{
					my $info = `$ntthal -a END1 -s1 $oligo_seq -s2 $seq`;
					chomp $info;
					my ($tm, $end31, $end32, $amplen, $mlen3, $len, $msum, $indel)=&dimer_amplify($info);
					($is_amplify, $atype, $eff)=&judge_amplify($tm, $end31, $end32, $amplen, $mlen3, $msum, $indel, $min_tm_omega, $min_tm_amp, $min_amplen, $min_meetlen); 
					if($is_amplify || defined $Nofilter){
						print O join("\t", "#Dimer", $id, $id2, $atype, $eff, $len, $tm, $end31, $end32, $amplen, $mlen3, $msum),"\n";
						print O "$ntthal -a END1 -s1 $oligo_seq -s2 $seq\n";
						print O $info, "\n";
					}
					}
					
					if(!$is_amplify || defined $Nofilter){
					my $info = `$ntthal -a END2 -s1 $oligo_seq -s2 $seq`;
					chomp $info;
					my ($tm, $end31, $end32, $amplen, $mlen3, $len, $msum, $indel)=&dimer_amplify($info);
					($is_amplify, $atype, $eff)=&judge_amplify($tm, $end31, $end32, $amplen, $mlen3, $msum, $indel, $min_tm_omega, $min_tm_amp, $min_amplen, $min_meetlen); 
					if($is_amplify || defined $Nofilter){
						
						print O join("\t", "#Dimer", $id, $id2, $atype, $eff, $len, $tm, $end31, $end32, $amplen, $mlen3, $msum),"\n";
						print O "$ntthal -a END2 -s1 $oligo_seq -s2 $seq\n";
						print O $info, "\n";
					}
					}

					if($atype ne "AmpEndMeet"){
					my $info = `$ntthal -a END1 -s1 $subseq -s2 $subs`;
					chomp $info;
					my ($tm, $end31, $end32, $amplen, $mlen3, $len, $msum, $indel)=&dimer_amplify($info);
					next if($msum>1 && $atype eq "AmpEndMap");
					$len+=length($seq)-$sublen+length($oligo_seq)-$sublen;
					($is_amplify, $atype, $eff)=&judge_amplify_endmeet($tm, $end31, $end32, $mlen3, $min_meetlen);
					if($is_amplify || defined $Nofilter){
						$amplen="NA";
						print O join("\t", "#EndDimer", $id, $id2, $atype, $eff, $len, $tm, $end31, $end32, $amplen, $mlen3, $msum),"\n";
						print O "$ntthal -a END1 -s1 $subseq -s2 $subs\n";
						print O $info, "\n";
					}
					}
				}
			}
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

	-t "Primer": 2 columns (IDs and Sequences), as following:
	xxx-1  AAAAAAAAAAAAAAAAA
	xxx-2  AAAAAAAAAAAAAAAAA
	xxx-P  AAAAAAAAAAAAAAAAA

Usage:
  Options:
  -i  <file>   Input primer file, forced
  -t  <str>    primer file type, "Result" (xxx.final.result) or "Primer", [$type]
  -k  <str>	   Key of output file, forced
  --SelfComplementary   Evalue self-complementary
  --MultiPlex           Evalue dimer among different primer pairs
  --Nofilter            No filter and output all dimer info.
  -mv        <int>      concentration of monovalent cations in mM, [$mv]
  -dv        <float>    concentration of divalent cations in mM, [$dv]
  -dNTP      <float>    concentration of deoxynycleotide triphosphate in mM, [$dNTP]
  -dna       <int>      concentration of DNA strands in nM, [$dna]

  -od <dir>	Dir of output file, default ./
  -h		 Help

USAGE
	print $usage;
	exit;
}

