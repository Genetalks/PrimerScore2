#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $BEGIN_TIME=time();
my $version="3.0.0";

# ********************************************************************
# Smith-Waterman Algorithm in Perl
# Author: litc
# Contact: litc\@biomarker.com.cn; litc\@qq.com
# ********************************************************************


# ------------------------------------------------------------------
# Get sequence
# ------------------------------------------------------------------
if (@ARGV<2) {
	USAGE();
}

my $seq1=shift;
my $seq2=shift;


# ------------------------------------------------------------------
# GetOptions and Glob score scheme
# ------------------------------------------------------------------
my $DEBGU_MODE=0;

#my $match_score     ||=  1;
#my $mismatch_score  ||=  -1;
#my $open_gap        ||=  -5;
#my $extend_gap      ||=  -3;
my $match_score     ||=  1;
my $mismatch_score  ||=  -0.8;
my $open_gap        ||=  -5;
my $extend_gap      ||=  -3;


GetOptions(
				"help|?" =>\&USAGE,
				"match_score:s"=>\$match_score,
				"mismatch_score:s"=>\$mismatch_score,
				"open_gap:s"=>\$open_gap,
				"extend_gap:s"=>\$extend_gap,
				"debug"=>\$DEBGU_MODE,
				) or &USAGE;

# guarantee sign of parameters
$match_score     =    abs($match_score);
$mismatch_score  =   -abs($mismatch_score);
$open_gap        =   -abs($open_gap);
$extend_gap      =   -abs($extend_gap);


# ------------------------------------------------------------------
# Smith-Waterman
# ------------------------------------------------------------------

#initialization
my @score_matrix;
my @trace_matrix;


my $seq1_length=length $seq1;
my $seq2_length=length $seq2;


#initialization top left cell.
$score_matrix[0][0]=0;
$trace_matrix[0][0]="done";

for (my $i=1;$i<=$seq1_length ;$i++) {
	$trace_matrix[0][$i]="left";
	$score_matrix[0][$i]=0;
}
for (my $i=1;$i<=$seq2_length ;$i++) {
	$trace_matrix[$i][0]="up";
	$score_matrix[$i][0]=0;
}


#split seq into array
my @seq1=split "",$seq1;
my @seq2=split "",$seq2;


#store the position of the end of alignment
my $max_i=1;
my $max_j=1;
my $max_score=0;


#calculate score matrix and traceback matrix
for (my $i=1;$i<=$seq2_length ;$i++) {
	for (my $j=1;$j<=$seq1_length ;$j++) {
		my $score=$seq1[$j-1] eq $seq2[$i-1] ? $match_score : $mismatch_score;
		my $Qdiag=$score_matrix[$i-1][$j-1]+$score;

		my $upGap=$trace_matrix[$i-1][$j] eq "up" ? $extend_gap : $open_gap;
		my $Qup=$score_matrix[$i-1][$j]+$upGap;

		my $leftGap=$trace_matrix[$i][$j-1] eq "left" ? $extend_gap : $open_gap;
		my $Qleft=$score_matrix[$i][$j-1]+$leftGap;

		if (0 > $Qdiag && 0 > $Qup && 0 > $Qleft) {
			$score_matrix[$i][$j]=0;
			$trace_matrix[$i][$j]="done";
		}elsif ($Qdiag > $Qup && $Qdiag > $Qleft) {
			$score_matrix[$i][$j]=$Qdiag;
			$trace_matrix[$i][$j]="diag";
		}elsif($Qup > $Qdiag && $Qup > $Qleft){
			$score_matrix[$i][$j]=$Qup;
			$trace_matrix[$i][$j]="up";
		}else{
			$score_matrix[$i][$j]=$Qleft;
			$trace_matrix[$i][$j]="left";
		}

		if ($score_matrix[$i][$j] > $max_score) {
			$max_score=$score_matrix[$i][$j];
			$max_i=$i;
			$max_j=$j;
		}
	}
}


#used for debug
if ($DEBGU_MODE) {
	for (my $i=0;$i<@score_matrix ;$i++) {
		for (my $j=0;$j<@{$score_matrix[$i]} ;$j++) {
			print $score_matrix[$i][$j],"\t";
		}
		print "\n";
	}
	print "\n********************************************************************\n";
	for (my $i=0;$i<@trace_matrix ;$i++) {
		for (my $j=0;$j<@{$trace_matrix[$i]} ;$j++) {
			print $trace_matrix[$i][$j],"\t";
		}
		print "\n";
	}
	print "\n********************************************************************\n";
}


#traceback for alignment sequence
my $str1="";
my $str2="";
my $alig;

my $i=$max_i;
my $j=$max_j;

my $match=0;
my $misMatch=0;
my $gapCount=0;
my $gapSize=0;


my ($index_s1, $index_s2, $index_e1, $index_e2) ;
while ($trace_matrix[$i][$j] ne "done" && $i>0 && $j>0) {
	$DEBGU_MODE && print $i,"\t",$j,"\t",$trace_matrix[$i][$j],"\t";
	if ($trace_matrix[$i][$j] eq "diag") {
		$str1.=$seq1[$j-1];
		$str2.=$seq2[$i-1];
		if ($seq1[$j-1] eq $seq2[$i-1]){
			if(length($str1)==1){$index_e1 = $j-1;}
			if(length($str2)==1){$index_e2 = $i-1;}
			$index_s1 = $j-1;
			$index_s2 = $i-1;
			$alig.="|";
			$match++;
		}else{
			$alig.="*";
			$misMatch++;
		}
		$i--;
		$j--;
		$DEBGU_MODE && print "diag","\n";
	}elsif($trace_matrix[$i][$j] eq "left"){
		$str1.=$seq1[$j-1];
		$str2.="-";
		$alig.="^";
		$j--;
		$gapSize++;
		$DEBGU_MODE && print "left","\n";
	}else{
		$str1.="-";
		$str2.=$seq2[$i-1];
		$alig.="-";
		$i--;
		$gapSize++;
		$DEBGU_MODE && print "up","\n";
	}
}
$DEBGU_MODE && print "\n********************************************************************\n";

#reverse traceback string
$str1=reverse $str1;
$str2=reverse $str2;
$alig=reverse $alig;
#print join("\n", $str1, $alig, $str2),"\n";
#print join("\t", $index_s1, $index_s2, $index_e1, $index_e2),"\n";

## complement seq
if($index_s1 !=0 || $index_s2 !=0){
	my $sub1 = substr($seq1, 0, $index_s1);
	my $sub2 = substr($seq2, 0, $index_s2);
	my $max_len;
	if($index_s1 > $index_s2){
		my $l = length $sub2;
		$str1 = $sub1.$str1;
		$alig = "#" x $l. $alig;
		$str2 = " " x ($index_s1 - $index_s2). $sub2. $str2;
		$alig = " " x ($index_s1 - $l) . $alig;
	}else{
		my $l = length $sub1;
		$str2 = $sub2. $str2;
		$alig = "#" x $l . $alig;
		$str1 = " " x ($index_s2 - $index_s1) . $sub1. $str1;
		$alig = " " x ($index_s2 - $l) . $alig;
	}
	#print join("\t", $index_s1, $index_s2, $sub1, $sub2),"\n";
	#print $alig,"\n";
}
my $imax1 = length($seq1)-1;
my $imax2 = length($seq2)-1;
if($index_e1 != $imax1 || $index_e2 != $imax2){
	my $sub1 = substr($seq1, $index_e1+1);
	my $sub2 = substr($seq2, $index_e2+1);
	my $alen1 = $imax1 - $index_e1;
	my $alen2 = $imax2 - $index_e2;
	my $dlen = abs($alen1 - $alen2);
#	print join("\t", $sub1, $sub2, $alen1, $alen2),"\n";
	if($alen1 > $alen2){
		$str1 .= $sub2;
		$str2 .= $sub2. " " x $dlen;
		$alig .= "#" x $alen2;
		$alig .= " " x ($alen1-$alen2);
	}else{
		$str2 .= $sub2;
		$str1 .= $sub1 . " " x $dlen;
		$alig .= "#" x $alen1;
		$alig .= " " x ($alen2-$alen1);
	}
	#print join("\t", $index_e1, $index_e2, $imax1, $imax2, $alen1, $alen2, $sub1, $sub2, $alig),"\n";
}




#exit if no alignment
if ($alig eq "") {
	exit(0);
}


# ------------------------------------------------------------------
# Calculate score, identity
# ------------------------------------------------------------------
if ($gapSize == 1) {
	$gapCount = 1;
}else{
	my $temp=$str1;
	$temp=~s/\-+/\-/g;
	foreach my $char (split "",$temp) {
		$gapCount++ if ($char eq "-") ;
	}
	$temp=$str2;
	$temp=~s/\-+/\-/g;
	foreach my $char (split "",$temp) {
		$gapCount++ if ($char eq "-") ;
	}
}

my $score=$match*$match_score + $misMatch*$mismatch_score + $gapCount*$open_gap + ($gapSize-$gapCount)*$extend_gap;
my $alignment_length=length $alig;
my $identity=sprintf "%.2f",$match/$alignment_length;


# ------------------------------------------------------------------
# Output result
# ------------------------------------------------------------------
print $str1,"\n";
print $alig,"\n";
print $str2,"\n";
print "\n";
print "Match     : $match","\n";
print "MisMatch  : $misMatch","\n";
print "Gap Count : $gapCount","\n";
print "Gap Size  : $gapSize","\n";
print "Score     : $score","\n";
print "Identity  : $identity","\n";
print "\n";



sub USAGE {#
	my $usage=<<"USAGE";
Program: sw.pl (Smith-Waterman alignment by litc)
Version: $version
Contact: Li Tiancheng <litc\@biomarker.com.cn> <ltc_gs\@qq.com>

Usage: perl $0 <seq1> <seq2> [option]
  Options:
  -match_score      <int>   match score, default 1
  -mismatch_score   <int>   mismatch score, default -1
  -open_gap         <int>   open gap score, default -5
  -extend_gap       <int>   extend gap score, default -3
  -debug                    Debug mode, default off
  -h                        Help

USAGE
	print $usage;
	exit;
}
