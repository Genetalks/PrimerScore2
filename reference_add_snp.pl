#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
require "$Bin/snp.pm";
require "$Bin/common.pm";

my $BEGIN_TIME=time();
my $version="1.0.0";
#######################################################################################

# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
my ($fref, $fvcf,$fkey,$outdir);
GetOptions(
				"help|?" =>\&USAGE,
				"ir:s"=>\$fref,
				"iv:s"=>\$fvcf,
				"k:s"=>\$fkey,
				"od:s"=>\$outdir,
				) or &USAGE;
&USAGE unless ($fref and $fvcf and $fkey);

$outdir||="./";
`mkdir $outdir`	unless (-d $outdir);
$outdir=AbsolutePath("dir",$outdir);


&SHOW_TIME("read reference file");
my @chr;
my %seq;
open(R, $fref) or die $!;
$/=">";
while(<R>){
	chomp;
	next if(/^$/);
	my ($head, @seq)=split /\n/, $_;
	my ($id)=split /\s+/, $head;
	print $id,"...\n";
	push @chr, $id;
	$id=~s/chr//;
	my $x = join("", @seq);
	@{$seq{$id}} = split //, $x;
}
close(R);

&SHOW_TIME("read in vcf file and add snp info");
if($fvcf=~/.gz/){
	open(V, "zcat $fvcf|") or die $!;
}else{
	open(V, $fvcf) or die $!;
}
$/="\n";
while(<V>){
	chomp;
	next if(/^$/ || /^#/);
	my ($c, $p, undef, $ref, $alt)=split /\s+/, $_;
	$c=~s/chr//;
	next if(!exists $seq{$c});
	&sequence_convert_snp($seq{$c}, $p, $ref, $alt);
}
close(V);

my $L=60;
open(O, ">$outdir/$fkey.add_snp.fasta") or die $!;
foreach my $c (@chr){
	print O ">$c\n";
	$c=~s/chr//;
	my $len = scalar @{$seq{$c}};
	for(my $i=0; $i<$len; $i+=$L){
		my $e = $i+$L-1;
		$e= $e >($len-1)? $len-1: $e;
#		for(my $j=$i; $j<=$e;$j++){
#			if(!defined $seq{$c}->[$j]){
#				print join("\t",$i,$e,$j),"\n";
#				print join("", @{$seq{$c}}[$i..$e]),"\n";
#				die;
#			}
#
#		}
		print O join("", @{$seq{$c}}[$i..$e]),"\n";
	}
}
close(O);

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
  -ir  <file>   Input reference file, forced
  -iv  <file>   Input vcf gz file, forced
  -k  <str>	   Key of output file, forced
  -od <dir>	   Dir of output file, default ./
  -h		   Help

USAGE
	print $usage;
	exit;
}

