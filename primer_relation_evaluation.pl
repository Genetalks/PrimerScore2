#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $BEGIN_TIME=time();
my $version="1.0.0";
#######################################################################################

# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
my ($fIn,$fkey,$outdir);
my $pnum = 1;
GetOptions(
				"help|?" =>\&USAGE,
				"i:s"=>\$fIn,
				"n:s"=>\$pnum,
				"k:s"=>\$fkey,
				"od:s"=>\$outdir,
				) or &USAGE;
&USAGE unless ($fIn and $fkey);

$outdir||="./";
`mkdir $outdir`	unless (-d $outdir);
$outdir=AbsolutePath("dir",$outdir);

my $ntthal = "/data/bioit/biodata/zenghp/software/primer3-2.4.0/src/ntthal";
my $primer3_config = "/data/bioit/biodata/zenghp/software/primer3-2.4.0/src/primer3_config/";

my %seq;
my @group;
open(I, $fIn) or die $!;
while(<I>){
	chomp;
	next if(/^$/);
	my ($id, @seq)=split /\s+/, $_;
	for(my $i=0; $i<$pnum; $i++){
		my $id_new = $id."_".$i;
		push @{$group[$i]}, $id_new;
		$seq{$id_new}=$seq[$i];
	}
}
close(I);


open(O,">$outdir/$fkey.relation")or die $!;
for(my $i=0; $i<$pnum; $i++){
	my $gn = scalar @{$group[$i]};
	for(my $j=0; $j<$gn; $j++){
		my $primer_seq = $seq{$group[$i][$j]};
		print O ">",$group[$i][$j],"\t",$primer_seq,"\n";
		#hairpin
		my $hairpin_result = `$ntthal -path $primer3_config -a HAIRPIN -s1 $primer_seq`;
		my ($hdG, $htm) = $hairpin_result =~/dG = ([\d\+\-\.]+)\tt = ([\d\+\-\.]+)/;
		if(defined $hdG){
			print O "#Hairpin\t",join("\t",$group[$i][$j],$group[$i][$j],$hdG, $htm),"\n";
			print O $hairpin_result,"\n";
		}
		my $dimer_result = `$ntthal -path $primer3_config -a END1 -s1 $primer_seq -s2 $primer_seq`;
		my ($ddG, $dtm) = $dimer_result =~/dG = ([\d\+\-\.]+)\tt = ([\d\+\-\.]+)/;
		if(defined $ddG ){
			print O "#Dimer1\t",join("\t",$group[$i][$j],$group[$i][$j],$ddG, $dtm),"\n";
			print O $dimer_result,"\n";
		}

		for(my $k=$j+1; $k<$gn; $k++){
			my $seq = $seq{$group[$i][$k]};
			my $dimer_result = `$ntthal -path $primer3_config -a END1 -s1 $primer_seq -s2 $seq`;
			my ($dG, $tm) = $dimer_result =~/dG = ([\d\+\-\.]+)\tt = ([\d\+\-\.]+)/;
			next if(!defined $tm || $tm<-30);
			print O "#Dimer2\t",join("\t",$group[$i][$j], $group[$i][$k], $dG, $tm),"\n";
			print O $dimer_result,"\n";
		}
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
  -i  <file>   Input primer list file, forced
  -k  <str>	   Key of output file, forced
  -n  <int>    primer num, [$pnum]
  -od <dir>	Dir of output file, default ./
  -h		 Help

USAGE
	print $usage;
	exit;
}

