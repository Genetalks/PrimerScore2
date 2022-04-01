#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
require "$Bin/path.pm"; 
require "$Bin/common.pm"; 

my $BEGIN_TIME=time();
my $version="1.0.0";
#######################################################################################

# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
our $blastn;
our $makeblastdb;
our $SAMTOOLS;
our $REF_GRCh37;
my ($foligo, $ftemp, $fprefix, $outdir);
my $threads = 10;
my $stm = 45;
my $max_bounds_num = 10000;
my $fdatabases = $REF_GRCh37;
GetOptions(
				"help|?" =>\&USAGE,
				"io:s"=>\$foligo,
				"it:s"=>\$ftemp,
				"id:s"=>\$fdatabases,
				"k:s"=>\$fprefix,
				"tm:s"=>\$stm,
				"mn:s"=>\$max_bounds_num,
				"t:s"=>\$threads,
				"od:s"=>\$outdir,
				) or &USAGE;
&USAGE unless ($foligo and $ftemp and $fdatabases and $fprefix);

$outdir||="./";
`mkdir $outdir`	unless (-d $outdir);
$outdir=AbsolutePath("dir",$outdir);

my $blast_qid_convert = "perl $Bin/blast_qid_convert.pl";
my $ptools = "$Bin/primer_tools/build/src/ptools -p $Bin/primer_tools/src/primer3_config/ -M $max_bounds_num";
my $sh;
open($sh, ">$outdir/$fprefix.bound.sh") or die $!;

my @fdatabase=split /,/, $fdatabases;
my $fbound = "$outdir/$fprefix.bound.info";
if(-e $fbound){
	`rm $fbound`;
}
my $is_mv=1;
foreach my $fdb(@fdatabase){
	my $dname = basename($fdb);
	if(!-e "$fdb.db.nsq"){
		`ln -s $fdb $outdir`;
		$fdb = "$outdir/$dname";
		&Run("$makeblastdb -in $fdb -dbtype nucl -out $fdb\.db -parse_seqids", $sh);
	}
	my $fkey = $fprefix."_".$dname;
	my $len = &max_template_len($ftemp);
	my $evalue = 385.74*$len-10836+8000;
	&Run("$blastn -db $fdb\.db -query $ftemp -task blastn -num_threads $threads -evalue $evalue -word_size 7 -outfmt '17 SQ' -out $outdir/$fkey.sam", $sh);
	&Run("$blast_qid_convert -if $ftemp -is $outdir/$fkey.sam -k $fkey -od $outdir", $sh);
	&Run("$SAMTOOLS index $outdir/$fkey.bam", $sh);
	&Run("$ptools -b $stm -i $foligo -f $ftemp -a $outdir/$fkey.bam -o $outdir/$fkey.bound.info -t $threads", $sh);
	if($is_mv){
		`mv $outdir/$fkey.bound.info $fbound`;
		$is_mv=0;
	}else{
		`cat >>$fbound $outdir/$fkey.bound.info`;
		`rm $outdir/$fkey.bound.info`;
	}
}

close($sh);




#######################################################################################
print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################

# ------------------------------------------------------------------
# sub function
# ------------------------------------------------------------------
sub max_template_len{
	my ($fa)=@_;
	open(I, "$fa") or die $!;
	$/=">";
	my $max_len = 0;
	while(<I>){
		next if(/^$/);
		my ($head, @line)=split /\n/, $_;
		my $seq = join("", @line);
		my $l = length($seq);
		if($l>$max_len){
			$max_len = $l;
		}
	}
	$/="\n";

	return $max_len;
}

sub USAGE {#
	my $usage=<<"USAGE";
Program:
Version: $version
Contact:zeng huaping<huaping.zeng\@genetalks.com> 

Usage:
  Options:
  -io  <file>   Input oligo evaluate file, forced
  -it  <file>   Input template file, forced
  -id  <files>  Input database files separated by "," to evalue specificity, [$fdatabases]
  -k  <str>	   Key of output file, forced
  -tm <int>    Min tm to bind and amplify template, [$stm]
  -mn <int>    Max bounds num of a oligo, if exceeded, it will be filtered: its bounds info will not be output. -1 means no filtering, [$max_bounds_num]
  -t  <int>    threads num, [$threads]
  -od <dir>	   Dir of output file, default ./
  -h		   Help

USAGE
	print $usage;
	exit;
}

