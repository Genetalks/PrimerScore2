#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
require "/data/bioit/biodata/zenghp/bin/perl_mylib/tool.lib.pl";
require "$Bin/snp.pm";
require "$Bin/algorithm.pm";

my $BEGIN_TIME=time();
my $version="1.0.0";
#######################################################################################

# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
my ($fIn,$fkey,$outdir);
my $samplen;
my $max_extend=80;
my $min_len=300;
my $type="Strict";
GetOptions(
				"help|?" =>\&USAGE,
				"i:s"=>\$fIn,
				"sn:s"=>\$samplen,
				"me:s"=>\$max_extend,
				"ml:s"=>\$min_len,
				"tp:s"=>\$type,
				"k:s"=>\$fkey,
				"od:s"=>\$outdir,
				) or &USAGE;
&USAGE unless ($fIn and $fkey and $samplen);

$outdir||="./";
`mkdir $outdir`	unless (-d $outdir);
$outdir=AbsolutePath("dir",$outdir);

my ($windsize, @value, @score);
if($type eq "Strict"){
	$windsize=25;
	@value=(1, 0.95, 0);
	@score=(2, 0, -2);
}elsif($type eq "Medium"){
	$windsize=10;
	@value=(1, 0.9, 0);
	@score=(2, 1, -1);
}else{
	$windsize=10;
	@value=(1, 0.8, 0);
	@score=(2, 1, -1);
}
open(I, $fIn) or die $!;
<I>;
<I>;
<I>;
#s1: r1, r2, r3, r4...
#s2: r1, r2, r3, r4...
#s3: r1, r2, r3, r4...
#...
my @seq; 
my @align;
my $n=0;
## first read: get sample list and spaces num
my @sample;
my $off; ## offset of seq and alignment
my $is_unit1=1;
my $is_start=0;
while(<I>){
	my @info;
	my $line=$_;
	for(my $i=0; $i<$samplen+1; $i++){
		chomp $line;
		if($is_unit1){
			my ($s, $seq)=split /\s+/, $line;
			push @sample, $s;
			if(!defined $off){
				$off=(length $line)-(length $seq);
			}
		}
		my $sq=substr($line, $off);
		push @info, $sq;
		$line=<I>;
	}
	pop @sample;
	$is_unit1=0;
	if($line!~/^$/){
		die "Wrong: non blank line!\n";
		print join("\n", @info, $line),"\n";
	}
	if($info[-1]!~/\*/){
		if($is_start){
			$n++;
			$is_start=0;
		}
		next;
	}else{
		$is_start=1;
		for(my $i=0; $i<$samplen; $i++){
			push @{$seq[$i][$n]}, split //, $info[$i];
		}
		push @{$align[$n]}, split //, $info[-1];
	}
}
close(I);

open(O, ">$outdir/$fkey.conserved.fa") or die $!;
open(D, ">$outdir/$fkey.conserved_addsnp.fa") or die $!;
open(F, ">$outdir/$fkey.conserved.info") or die $!;
my $regionn=scalar @align;
our @Ins;
$n=1;
for(my $i=0; $i<$regionn; $i++){
	my @ratio=&window_ratio($align[$i], $windsize,"CLW_mismatch");
	my @segs=&segment(\@ratio, \@value, \@score, $max_extend, "Number", 1); #my ($aarr, $avalue, $ascore, $extend, $type)=@_;
	my @result=((0) x scalar @{$align[$i]});
	for(my $j=0; $j<@segs; $j++){
		my ($s, $e)=@{$segs[$j]};
		for(my $x=$s; $x<=$e+$windsize-1; $x++){
			$result[$x]=1;
		}
		$e+=$windsize-1;
		next if($e-$s+1<$min_len);
		my $tseq;
		my $dseq;## seq with degenerated base
		my $len_insert=0;
		for(my $x=$s; $x<=$e; $x++){
			my %base;
			##current is insert or not? yes: len++ and next; 
			if($align[$i][$x] ne "*"){
				for(my $k=0; $k<$samplen; $k++){
					$base{$seq[$k][$i][$x]}=1;
				}
				if(exists $base{"-"}){
					if($seq[0][$i][$x] eq "-"){## ref -: Ins
						$len_insert++;
						next;
					}
				}
			}

			if($len_insert>0){## handle insert
				if(defined $dseq){
					chop $dseq;
				}
				my $ix=($len_insert-1)>=3? 3: ($len_insert-1);
				$dseq.=$Ins[$ix];
				$len_insert=0;
			}

			$tseq.=$seq[0][$i][$x];
			if($align[$i][$x] eq "*"){
				$dseq.=$seq[0][$i][$x];
			}else{
				if(exists $base{"-"}){## Del
					$dseq.="E";
				}else{
					my $dbase = &snp_to_degenerate(join(",",keys %base));
					if($dbase eq "Error"){
						die "Wrong snp: ", join(",",keys %base), "\n";
					}
					$dseq.=$dbase;
				}
			}
		}
		my $l=length $tseq;
		print O ">Region$n\t$l\n";
		print O $tseq,"\n";
		print D ">Region$n\t$l\n";
		print D $dseq,"\n";
		$n++;
	}
	print F ">\n";
	for(my $k=0; $k<$samplen; $k++){
		print F join("", @{$seq[$k][$i]}),"\n";
	}
	print F join("", @{$align[$i]}),"\n";
	print F join("", @result),"\n\n";
}
close(O);
close(D);
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
  -sn <int>    Sample number, forced
  -me <int>    Max extend length, [$max_extend]
  -ml <int>    Min seq length to filter, [$min_len]
  -tp <str>    Segment handle type: "Strict", "Medium", "Loose", [$type]
  -od <dir>	   Dir of output file, default ./
  -h		   Help

USAGE
	print $usage;
	exit;
}

