#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
require "$Bin/common.pm";
require "$Bin/product.pm";
require "$Bin/path.pm";

my $BEGIN_TIME=time();
my $version="1.0.0";
#######################################################################################

# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
my ($foligo,$fbound,$fprobe,$fkey,$outdir);
our $REF_HG19;
my $fdatabases = $REF_HG19;
my $ptype = "face-to-face";
my $PCRsize=1000;
my $min_eff=0.01;
my $opt_tm = 60;
my $opt_tm_probe = 70;
my $range;
my $thread = 3;
my $etype="SinglePlex";
my $AllEvalue;
GetOptions(
				"help|?" =>\&USAGE,
				"io:s"=>\$foligo,
				"id:s"=>\$fdatabases,
				"k:s"=>\$fkey,
				"tp:s"=>\$ptype,
				"ep:s"=>\$etype,
				"AllEvalue:s"=>\$AllEvalue,
				"rd:s"=>\$range,
				"tm:s"=>\$opt_tm,
				"tmb:s"=>\$opt_tm_probe,
				"sz:s"=>\$PCRsize,
				"me:s"=>\$min_eff,
				"td:s"=>\$thread,
				"od:s"=>\$outdir,
				) or &USAGE;
&USAGE unless ($foligo and $fkey);
$outdir||="./";
`mkdir $outdir`	unless (-d $outdir);
$outdir=AbsolutePath("dir",$outdir);

if($ptype eq "face-to-face"){
	$range=defined $range? $range:"70,200";
}elsif($ptype eq "Nested"){
	$range=defined $range? $range:"5,30";
}elsif($ptype eq "back-to-back"){
	$range=defined $range? $range:"0,15";
}else{
	die "Wrong primer type! Must be face-to-face, back-to-back or Nested!\n";
}
my ($mind, $maxd)=split /,/, $range;

### evalue single oligo
&Run("perl $Bin/oligo_evaluation.pl --nohead -p $foligo -d $fdatabases -thread $thread -stm 45 --NoFilter -k $fkey -od $outdir");

### evalue combination oligos
my %bound;
my %boundm;
&SHOW_TIME("#Read in Bound file");
open(B, "$outdir/$fkey.bound.info") or die $!;
while(<B>){
	chomp;
	my ($id, $strand, $chr, $pos3, $seq, $tm, $end_match, $mvisual)=split /\t/, $_;
	my ($tid, $type)=$id=~/(\S+)-([12LRP])/;
	if(!defined $tid || !defined $type){
		die "Wrong oligo ID $id: Primer ID must end by -L/R/1/2, eg. xxx-L, xxx-R, xxx-1, xxx-2; Probe ID must end by -P, eg. xxx-P!\n";
	}
	my $len = length $seq;
	my $pos5=$strand eq "+"? $pos3-$len+1: $pos3+$len-1;
	push @{$bound{$type}{$tid}{$chr}{$strand}}, [$pos3, $pos5, $tm, $end_match, $mvisual, $tid];
}
close(B);

my %product;
my %record;
my @tps = sort {$a cmp $b} keys %bound;
&SHOW_TIME("#Evalue");
foreach my $tp1(@tps){
	next if($tp1 eq "P"); 
	## AllEvalue: L<->R L<->L R<->R
	next if($etype eq "SinglePlex" && !defined $AllEvalue && $tp1 ne $tps[0]); ## only evalue one type 
	foreach my $tid1 (sort {$a cmp $b} keys %{$bound{$tp1}}){
		foreach my $tp2(@tps){
			next if($tp2 eq "P");
			next if(!defined $AllEvalue && $tp2 eq $tp1); ## Only evalue product between pairs with different types, L<-->R; 1<-->2
			foreach my $tid2(keys %{$bound{$tp2}}){
				next if($etype eq "SinglePlex" && $tid2 ne $tid1); ## SinglePlex: not evalue product between different tid
				next if(exists $record{"pro"}{$tid2."-".$tp2}{$tid1."-".$tp1});
				&caculate_product($tid1, $tid1."-".$tp1, $tid2, $tid2."-".$tp2, $bound{$tp1}{$tid1}, $bound{$tp2}{$tid2}, $ptype,\%product, \%record, $PCRsize, $opt_tm, $mind, $maxd, $min_eff);
			}
		}
	}
}

my %productp;
if(exists $bound{"P"}){#Probe
	foreach my $tid(keys %{$bound{"P"}}){
		# probe num on products
		&probe_bounds_on_products($tid."-P", $bound{"P"}{$tid}, $product{$tid}{$tid}, \%productp, $opt_tm_probe); #my ($id, $abound, $aprod, $aresult)=@_;
	}
}

### output
open(I, "$outdir/$fkey.evaluation.out") or die $!;
open(O, ">$outdir/$fkey.final.evaluation") or die $!;
print O "#ID\tSeq\tLen\tTm\tGC\tHairpin\tEND_Dimer\tANY_Dimer\tEndANum\tEndStability\tSNP\tPoly\tOligoBound\tBoundNum\tHighestTm\tHighestInfo\n";
open(P, ">$outdir/$fkey.final.pair.product") or die $!;
while(<I>){
	chomp;
	next if(/^$/);
	my ($id, @info)=split /\t/, $_;
	my $sbinfo=pop @info;
	my $sbeff=pop @info;
	my $sbnum=pop @info;
	my $sbd = $sbnum."|".$sbeff;
	my ($tid,$tp)=$id=~/(\S+)-([12LRP])/;
	my ($pnum, $apeff, $apinfos);
	if($tp ne "P"){
		($pnum, $apeff, $apinfos)=&get_highest_bound($product{$tid}{$tid}, 1000000);
	}else{
		($pnum, $apeff, $apinfos)=&get_highest_bound($productp{$id}, 1000000);
	}
	my ($peffs, $pinfos);
	if($pnum!~/\+/ && $pnum<=3){
		$peffs = join(",", @{$apeff});
		$pinfos = join(";", @{$apinfos});
	}else{
		$peffs = join(",", @{$apeff}[0..2]);
		$pinfos = join(";", @{$apinfos}[0..2]);
	}
	print O join("\t", $id, @info, $sbd, $pnum, $peffs, $pinfos),"\n";
	print P ">$id\t$pnum\n";
	my @peff=@{$apeff};
	my @pinfo=@{$apinfos};
	for(my $i=0; $i<@{$apeff}; $i++){
		print P join("\t", $apinfos->[$i], $apeff->[$i]),"\n";
	}
}
close(I);
close(O);
close(P);


open(C, ">$outdir/$fkey.final.cross.product") or die $!;
foreach my $tid1(keys %product){
	foreach my $tid2(keys %{$product{$tid1}}){
		next if($tid1 eq $tid2);
		
		my ($pnum, $apeff, $apinfos)=&get_highest_bound($product{$tid1}{$tid2}, 1000000);
		print C ">$tid1\t$tid2\t$pnum\n";
		my @peff=@{$apeff};
		my @pinfo=@{$apinfos};
		for(my $i=0; $i<@{$apeff}; $i++){
			print C join("\t", $apinfos->[$i], $apeff->[$i]),"\n";
		}
	}
}
close(C);
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

	-rd: distance range between primer pair(P2-P1)

      face-to-face: |---> P1             distance range: 70,200
                          P2 <---|  

      back-to-back:  x <---| P2          distance range(>0): 0,15 
                      P1 |---> x           

      back-to-back:  x <---| P2         distance range(<0): -40,-10
                           P1 |---> x          

            Nested: P1 --->|   x        distance range: 5,30
                      P2 --->| x                

            Nested: P2 --->|   x        distance range: -30,-5
                      P1 --->| x                

Usage:
  Options:
  -io        <file>   Input oligo file(ID must like: xxx-L, xxx-R, xxx-P, xxx-1, xxx-2), forced
  -id        <files>  Input database files separated by "," to evalue specificity, [$fdatabases]
  -k         <str>    Key of output file, forced
  -tp        <str>    primer type, "face-to-face", "back-to-back", "Nested", ["face-to-face"]
  -ep        <str>    evalue type, "SinglePlex", "MultiPlex", ["SinglePlex"]
  --AllEvalue         All Evalue, not only L<->R, contain L<->L and R<->R.
  -td        <int>    thread in bwa, [$thread]
  -od        <dir>	  Dir of output file, default ./

  Caculating efficiency options:
  -sz        <int>    max PCR fragment size, [$PCRsize]
  -tm        <int>    optimal tm of primer, [$opt_tm]
  -tmb       <int>    optimal tm of probe, [$opt_tm_probe]
  -rd      <int,int>  distance range between primer pair, optional
                      default:
                      when -tp is "face-to-face", [70,200]
                                  "back-to-back", [0,15]
                                  "Nested",       [5,30]
  -me       <float>   min efficiency to consider a product, [$min_eff]

  -h		   Help

USAGE
	print $usage;
	exit;
}

