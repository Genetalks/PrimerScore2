#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
require "$Bin/common.pm";
require "$Bin/algorithm.pm";

my $BEGIN_TIME=time();
my $version="1.0.0";
#######################################################################################

# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
my ($fIn,$fkey,$outdir);
my $min_len=300;
my $type="Strict";
my $max_extend=50;
GetOptions(
				"help|?" =>\&USAGE,
				"i:s"=>\$fIn,
				"k:s"=>\$fkey,
				"me:s"=>\$max_extend,
				"ml:s"=>\$min_len,
				"tp:s"=>\$type,
				"od:s"=>\$outdir,
				) or &USAGE;
&USAGE unless ($fIn and $fkey);

$outdir||="./";
`mkdir $outdir`	unless (-d $outdir);
$outdir=AbsolutePath("dir",$outdir);

my ($windsize, @value, @score);
if($type eq "Strict"){
	$windsize=16;
	@value=(0.8,0.7,0.6,0.4,0.3,0.2,0);
	@score=(-3,-1,0,2,0,-1,-3);
}else{
	$windsize=16;
	@value=(0.8,0.7,0.6,0.4,0.3,0.2,0);
	@score=(-3,-1,0.5,2,0.5,-1,-3);
}

open(I, $fIn) or die $!;
open(O, ">$outdir/$fkey.filterGC.fasta") or die $!;
open(F, ">$outdir/$fkey.filterGC.info") or die $!;
$/=">";
while(<I>){
	chomp;
	next if(/^$/);
	my ($head, @seq)=split /\n/, $_;
	my ($id)=split /\s+/, $head;
	my $seq=join("", @seq);
	my @unit = split //, $seq;
	my @GC=&window_ratio(\@unit, $windsize, "GC");
	my @segs = &segment(\@GC, \@value, \@score, $max_extend, "Number", 1);
	my $n=1;
	my @result=((0) x scalar @unit);
	for(my $j=0; $j<@segs; $j++){
		my ($s, $e)=@{$segs[$j]};
		$e+=$windsize-1;
		@result[$s..$e]=((1) x ($e-$s+1));
		next if($e-$s+1<$min_len);
		print O ">$id\_$n\n";
		print O join("", @unit[$s..$e]),"\n";
		$n++;
	}
	print F join("", @unit),"\n";
	print F join("", @result),"\n";
}
close(I);
close(O);
close(F);




#######################################################################################
print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################

# ------------------------------------------------------------------
# sub function
# ------------------------------------------------------------------



sub USAGE {#
	my $usage=<<"USAGE";
Program:
Version: $version
Contact:zeng huaping<huaping.zeng\@genetalks.com> 

Usage:
  Options:
  -i  <file>   Input file, forced
  -k  <str>	   Key of output file, forced
  -me <int>    Max extend length, [$max_extend]
  -ml <int>    Min seq length to filter, [$min_len]
  -tp <str>    Segment handle type: "Strict", "Loose", [$type]
  -od <dir>	   Dir of output file, default ./
  -h		   Help

USAGE
	print $usage;
	exit;
}

