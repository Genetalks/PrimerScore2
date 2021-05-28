#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
require "$Bin/score.pm";
require "$Bin/self_lib.pm";

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
GetOptions(
				"help|?" =>\&USAGE,
				"i:s"=>\$foligo,
				"maxl:s"=>\$max_len,
				"minl:s"=>\$min_len,
				"opttm:s"=>\$opt_tm,
				"k:s"=>\$fkey,
				"od:s"=>\$outdir,
				) or &USAGE;
&USAGE unless ($foligo and $fkey);

$outdir||="./";
`mkdir $outdir`	unless (-d $outdir);
$outdir=AbsolutePath("dir",$outdir);
my $dlen=$max_len-$min_len;
my @endA = (0, 0, 0, 2);
my @len = ($min_len,$max_len-$dlen*0.5,$min_len, $max_len+1);
my @tm = ($opt_tm-1, $opt_tm+1, $opt_tm-5, $opt_tm+5);
my @self = (-50, 40, -50, 55); ## self tm
my @CGd = (0.1, 1, 0, 1);
my $max_self_tm=$self[3];
my $max_bound_num=100;

my $fulls=10;
open(F,">$outdir/$fkey.probe.filter") or die $!;
open(O, ">$outdir/$fkey.probe.score") or die $!;
print O "#ID\tSeq\tLen\tScore\tScoreInfo\tTM\tGC\tHairpin\tEND\tANY\tSNP\tPoly\tBoundNum\tBoundTM\tBoundInfo\n";
open(P, $foligo) or die $!;
while(<P>){
	chomp;
	my ($id, $seq, $len, $tm, $gc, $hairpin, $END, $ANY, undef, undef, $snp, $poly, $bnum, $btm)=split /\t/, $_;
	
	## filter
	if($tm<$opt_tm-5 || $tm>$opt_tm+5){
		print F "TM\t$_\n";
		next;
	}
	my ($is_G5, $CGd) = &G_content($seq);
	if($is_G5){
		print F "5endG\t$_\n";
		next;
	}
	if($CGd<=0){
		print F "GCdiff\t$_\n";
		next;
	}
	if($bnum>$max_bound_num){
		print F "Bound\t$_\n";
		next;
	}
	## score
	my $slen=int(&score_single($len, $fulls, @len)+0.5);## round: int(x+0.5)
	my $stm=int(&score_single($tm, $fulls, @tm)+0.5);
	my $self = &max($hairpin, $END, $ANY);
	my $sself=int(&score_single($self, $fulls, @self)+0.5);
	my $ssnp = int(&SNP_score($snp, $len, "Probe")*$fulls +0.5);
	my $spoly = int(&poly_score($poly, $len, "Probe")*$fulls +0.5);
	my $sCGd=int(&score_single($CGd, $fulls, @CGd)+0.5);
	#specificity: bound
	my $sbound=&bound_score($bnum, $btm, $fulls);
	my @score = ($slen, $stm, $sself, $ssnp, $spoly, $sbound, $sCGd);
	my @weight =(0.5,   3,     1,      2,    2,      0.5,     1);
	my $sadd=0;
	for(my $i=0; $i<@score; $i++){
		$score[$i]=$score[$i]<0? 0: $score[$i];
		$sadd+=$weight[$i]*$score[$i];
	}
	my $score_info=join(",", @score);
	print O join("\t", $id, $seq, $len, $sadd, $score_info, $tm, $gc, $hairpin, $END, $ANY, $snp, $poly, $bnum, $btm),"\n";
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

  -minl  <int>  min len of probe, [$min_len]
  -maxl  <int>  max len of probe, [$max_len]
  -opttm <int>  opt tm of probe, [$opt_tm]
  -od <dir>	   Dir of output file, default ./
  -h		   Help

USAGE
	print $usage;
	exit;
}

