#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
require "/data/bioit/biodata/zenghp/bin/perl_mylib/tool.lib.pl";

my $BEGIN_TIME=time();
my $version="1.0.0";
#######################################################################################

# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
my ($ftem, $fsam,$fkey,$outdir);
GetOptions(
				"help|?" =>\&USAGE,
				"if:s"=>\$ftem,
				"is:s"=>\$fsam,
				"k:s"=>\$fkey,
				"od:s"=>\$outdir,
				) or &USAGE;
&USAGE unless ($ftem and $fsam and $fkey);

$outdir||="./";
`mkdir $outdir`	unless (-d $outdir);
$outdir=AbsolutePath("dir",$outdir);

my $n=0;
my %query;
open(I, $ftem) or die $!;
$/=">";
while(<I>){
	chomp;
	next if(/^$/);
	my ($id)=split /\s+/, $_;
	$n++;
	$query{"Query_$n"}=$id;
}
close(I);

$/="\n";
open(O, "|samtools view -bS -> $outdir/$fkey.bam") or die $!;
open(S, $fsam) or die $!;
while(<S>){
	chomp;
	next if(/^$/);
	
	if(/^\@/){
		if(/^\@SQ/){
			my ($qid)=$_=~/SN:(\S+)/;
			my $id = $query{$qid};
			$_=~s/$qid/$id/;

		}
		print O $_,"\n";
		next;
	}

	my @unit = split /\t/, $_;
	$unit[2]=$query{$unit[2]};
	print O join("\t", @unit),"\n";
}
close(S);
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
  -if  <file>   Input fasta file, forced
  -is  <file>   Input blast(-m17) sam file, forced
  -k  <str>	   Key of output bam file, forced
  -od <dir>	   Dir of output file, default ./
  -h		   Help

USAGE
	print $usage;
	exit;
}

