#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
require "$Bin/path.pm";
my $BEGIN_TIME=time();
my $version="1.0.0";
#######################################################################################

# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
my ($fIn,$fkey,$outdir);
GetOptions(
				"help|?" =>\&USAGE,
				"i:s"=>\$fIn,
				"k:s"=>\$fkey,
				"od:s"=>\$outdir,
				) or &USAGE;
&USAGE unless ($fIn and $fkey);

$outdir||="./";
`mkdir $outdir`	unless (-d $outdir);
$outdir=AbsolutePath("dir",$outdir);

our $PATH_PRIMER3;
my $ntthal = "$PATH_PRIMER3/src/ntthal";
my $primer3_config = "$PATH_PRIMER3/src/primer3_config/";

my %seq;
my @group;
open(I, $fIn) or die $!;
while(<I>){
	chomp;
	next if(/^$/ || /^#/);
	my ($id, $seq)=(split /\s+/, $_)[3,4];
	$seq{$id}=$seq;
	push @group, $id;
}
close(I);


open(O,">$outdir/$fkey.dimer.check")or die $!;
for(my $j=0; $j<@group; $j++){
	my $primer_seq = $seq{$group[$j]};
	print O ">",$group[$j],"\t",$primer_seq,"\n";
	for(my $k=$j+1; $k<@group; $k++){
		my $seq = $seq{$group[$k]};
		my $dimer_result = `$ntthal -path $primer3_config -a END1 -s1 $primer_seq -s2 $seq`;
		my ($dG, $tm) = $dimer_result =~/dG = ([\d\+\-\.]+)\tt = ([\d\+\-\.]+)/;
		next if(!defined $tm || $tm<-30);
		print O "#Dimer\t",join("\t",$group[$j], $group[$k], $dG, $tm),"\n";
		print O $dimer_result,"\n";
	}
}

close(O);

#######################################################################################
print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################

# ------------------------------------------------------------------
# sub function
# ------------------------------------------------------------------
sub AbsolutePath
{		#获取指定目录或文件的决定路径
		my ($type,$input) = @_;

		my $return;
	$/="\n";

		if ($type eq 'dir')
		{
				my $pwd = `pwd`;
				chomp $pwd;
				chdir($input);
				$return = `pwd`;
				chomp $return;
				chdir($pwd);
		}
		elsif($type eq 'file')
		{
				my $pwd = `pwd`;
				chomp $pwd;

				my $dir=dirname($input);
				my $file=basename($input);
				chdir($dir);
				$return = `pwd`;
				chomp $return;
				$return .="\/".$file;
				chdir($pwd);
		}
		return $return;
}

sub GetTime {
	my ($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst)=localtime(time());
	return sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}


sub USAGE {#
	my $usage=<<"USAGE";
Program:
Version: $version
Contact:zeng huaping<huaping.zeng\@genetalks.com> 

Usage:
  Options:
  -i  <file>   Input primer file, forced
  -k  <str>	   Key of output file, forced
  -od <dir>	Dir of output file, default ./
  -h		 Help

USAGE
	print $usage;
	exit;
}

