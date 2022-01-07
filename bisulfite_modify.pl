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
my ($fIn,$fkey,$outdir);
my $Unmeth;
GetOptions(
				"help|?" =>\&USAGE,
				"i:s"=>\$fIn,
				"k:s"=>\$fkey,
				"Unmeth:s"=>\$Unmeth,
				"od:s"=>\$outdir,
				) or &USAGE;
&USAGE unless ($fIn and $fkey);

$outdir||="./";
`mkdir $outdir`	unless (-d $outdir);
$outdir=AbsolutePath("dir",$outdir);
if(defined $Unmeth){
	$fkey=$fkey."_Unmeth";
}


$/=">";
open(O, ">$outdir/$fkey.bisulfite.fasta") or die $!;
open(M, ">$outdir/$fkey.bisulfite.mark.fasta") or die $!;
open(RO, ">$outdir/$fkey.rev.bisulfite.fasta") or die $!;
open(RM, ">$outdir/$fkey.rev.bisulfite.mark.fasta") or die $!;
open (I, $fIn) or die $!;
while(<I>){
	chomp;
	next if(/^$/);
	my ($head, @line)=split /\n/, $_;
	my ($id)=split /\s+/, $head;
	my $seq=join("", @line);
	$seq=uc($seq);
	my ($aseqct, $amarkct, $aseqga, $amarkga)=&bisulfite_modify($seq, length $line[0], $Unmeth);

	print O ">$id\n";
	print O join("\n", @{$aseqct}), "\n";
	print M ">$id\n";
	print M join("\n", @{$amarkct}), "\n";
	
	print RO ">$id\n";
	print RO join("\n", @{$aseqga}), "\n";
	print RM ">$id\n";
	print RM join("\n", @{$amarkga}), "\n";

	#print ">$id\n";
	#print $seq,"\n";
	#print join("", @{$amarkct}), "\n";
	#print  join("", @{$aseqct}), "\n";
	#print "\n";	
	#print $seq,"\n";
	#print  join("", @{$amarkga}), "\n";
	#print  join("", @{$aseqga}), "\n";
	#print "\n";	
	undef $seq;
	undef @{$aseqct};
	undef @{$amarkct};
	undef @{$aseqga};
	undef @{$amarkga};
}
close(I);

close(O);
close(M);
close(RO);
close(RM);



#######################################################################################
print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################

# ------------------------------------------------------------------
# sub function
# ------------------------------------------------------------------


sub bisulfite_modify{
	my ($seq, $w, $Unmeth)=@_;
	my @unit = split //, $seq;
	my @ct;
	my @ctm;
	my @ga;
	my @gam;
	for(my $i=0; $i<@unit; $i++){
		if(!defined $Unmeth && $unit[$i] eq "C" && $i+1<=$#unit && $unit[$i+1] eq "G"){ ## meth
			$ct[$i]="C";
			$ct[$i+1]="G";
			$ctm[$i]="+";
			$ctm[$i+1]="*";

			$ga[$i]="C";
			$ga[$i+1]="G";
			$gam[$i]="*";
			$gam[$i+1]="+";
			$i++;
			next;
		}

		if($unit[$i] eq "C"){
			$ct[$i]="T";
			$ctm[$i]=":";
		}else{
			$ct[$i]=$unit[$i];
			$ctm[$i]="|";
		}

		if($unit[$i] eq "G"){
			$ga[$i]="A";
			$gam[$i]=":";
		}else{
			$ga[$i]=$unit[$i];
			$gam[$i]="|";
		}
	}

	
	my @seqct;
	my @seqga;
	my @markct;
	my @markga;
	my $i=0;
	while($i<=$#unit){
		my $iend = $i+$w-1;
		if($iend>$#unit){
			$iend = $#unit;
		}
		push @seqct, join("", @ct[$i..$iend]);
		push @seqga, join("", @ga[$i..$iend]);
		push @markct, join("", @ctm[$i..$iend]);
		push @markga, join("", @gam[$i..$iend]);
		$i+=$w;
	}
	undef @ct;
	undef @ga;
	undef @ctm;
	undef @gam;
	return (\@seqct, \@markct, \@seqga, \@markga);
}


sub USAGE {#
	my $usage=<<"USAGE";
Program:
Version: $version
Contact:zeng huaping<huaping.zeng\@genetalks.com> 

Usage:
  Options:
  -i  <file>   Input fasta file, forced
  -k  <str>	   Key of output file, forced
  --Unmeth     Un-methylation
  -od <dir>	   Dir of output file, default ./
  -h		   Help

USAGE
	print $usage;
	exit;
}

