#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
require "$Bin/score.pm";
require "$Bin/common.pm";
require "$Bin/self_lib.pm";
require "$Bin/math.pm";

my $BEGIN_TIME=time();
my $version="1.0.0";
#######################################################################################

# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
my ($foligo,$fkey,$outdir);
my $min_len=20;
my $max_len=36;
my $opt_tm=70;
my $NoFilter;
GetOptions(
				"help|?" =>\&USAGE,
				"i:s"=>\$foligo,
				"maxl:s"=>\$max_len,
				"minl:s"=>\$min_len,
				"opttm:s"=>\$opt_tm,
				"k:s"=>\$fkey,
				"NoFilter:s"=>\$NoFilter,
				"od:s"=>\$outdir,
				) or &USAGE;
&USAGE unless ($foligo and $fkey);

$outdir||="./";
`mkdir $outdir`	unless (-d $outdir);
$outdir=AbsolutePath("dir",$outdir);

my $max_bound_num=100;
open(F,">$outdir/$fkey.probe.filter") or die $!;
open(O, ">$outdir/$fkey.probe.score") or die $!;
#$slen, $stm, $sself, $ssnp, $spoly, $sbound, $sCGd
print O "##ScoreInfo: scores of length, tm, self-complementary, CG content diff, 5endG, snp, poly, bounding\n";
print O "#ID\tSeq\tLen\tScore\tScoreInfo\tTM\tGC\tHairpin\tDimerType\tDimerSize\tSNP\tPoly\tBoundNum\tBoundTM\tBoundInfo\n";
open(P, $foligo) or die $!;
while(<P>){
	chomp;
	my ($id, $seq, $len, $tm, $gc, $hairpin, $dimert, $dimers, undef, undef, $snp, $poly, $bnum, $btm)=split /\t/, $_;
	## filter
	if(!defined $NoFilter && ($tm<$opt_tm-8 || $tm>$opt_tm+8)){
		print F "TM\t$_\n";
		next;
	}
	my ($is_G5, $CGd) = &G_content($seq);
	if(!defined $NoFilter && $is_G5){
		print F "5endG\t$_\n";
		next;
	}
	if(!defined $NoFilter && $CGd<=0){
		print F "C<G\t$_\n";
		next;
	}
#	if($bnum>$max_bound_num){
#		print F "Bound\t$_\n";
#		next;
#	}
	## score
	my ($sadd, $score_info)=&probe_oligo_score($opt_tm, $len, $tm, $gc, $hairpin, $snp, $poly, $bnum, $btm, $is_G5, $CGd);
	print O join("\t", $id, $seq, $len, $sadd, $score_info, $tm, $gc, $hairpin, $dimert, $dimers, $snp, $poly, $bnum, $btm),"\n";
}
close(P);
close(O);
close(F);


#######################################################################################
print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################

# ------------------------------------------------------------------
# sub function
# ------------------------------------------------------------------

sub G_content{
	my ($seq)=@_;
	$seq = uc($seq);
	my @unit = split //, $seq;
	my $is_G5 = $unit[0] eq "G"? 1: 0;
	my ($n, $nc, $ng)=(0,0,0);

	for(my $i=0; $i<@unit; $i++){
		$n++;
		if($unit[$i] eq "C"){
			$nc++;
		}elsif($unit[$i] eq "G"){
			$ng++;
		}
	}
	return($is_G5, ($nc-$ng)/$n);
}



sub USAGE {#
	my $usage=<<"USAGE";
Program:
Version: $version
Contact:zeng huaping<huaping.zeng\@genetalks.com> 

Usage:
  Options:
  -i  <file>   Input oligo evaluation file, forced
  -k  <str>	   Key of output file, forced

  --NoFilter    Not filter any probes
  -minl  <int>  min len of probe, [$min_len]
  -maxl  <int>  max len of probe, [$max_len]
  -opttm <int>  opt tm of probe, [$opt_tm]
  -od <dir>	   Dir of output file, default ./
  -h		   Help

USAGE
	print $usage;
	exit;
}

