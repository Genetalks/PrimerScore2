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
my ($ftxt, $flist,$fkey,$outdir, $col, $Head);
my $type = "Extract";
GetOptions(
				"help|?" =>\&USAGE,
				"i:s"=>\$ftxt,
				"l:s"=>\$flist,
				"k:s"=>\$fkey,
				"c:s"=>\$col,
				"Type:s"=>\$type,
				"Head:s"=>\$Head,
				"od:s"=>\$outdir,
				) or &USAGE;
&USAGE unless ($ftxt and $flist and $fkey);

$outdir||="./";
`mkdir $outdir`	unless (-d $outdir);
$outdir=AbsolutePath("dir",$outdir);
$col=defined $col? $col: 0;

my %list;
open (L, $flist) or die $!;
while (<L>) {
	chomp;
	next if(/^$/);
	my ($id)=split;
	$list{$id}=1;
}
close (L);

open (O, ">$outdir/$fkey.extract") or die $!;
open (I, "zcat $ftxt|") or die $!;
if (defined $Head) {
	my $h = <I>;
	print O $h;
}

my %search;
while (<I>) {
	chomp;
	next if(/^$/ || /^#/);
	my ($id)=(split /\t/, $_)[$col];
	if (exists $list{$id}) {
		if($type eq "Extract"){
			print O $_,"\n";
		}
		$search{$id}++;
	}else{
		if($type eq "Filter"){
			print O $_,"\n";
		}
	}
}
close (O);
close (I);

foreach my $id (keys %list) {
	if (!exists $search{$id}) {
		print "No Search ID: $id\n";
	}
}




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
Contact:zeng huaping<zenghp\@biomarker.com.cn> 

Usage:
  Options:
  -i  <file>   Input gzipped txt file(xx.gz), forced
  -l  <file>   Input ID list file, forced
  -k  <str>    Key of output file, forced
  -type <str>  Extract or Filter, [$type]

  -c  <int>    column of ID in txt file, [0]
  --Head       output txt head
  -od <dir>    Dir of output file, default ./
  -h         Help

USAGE
	print $usage;
	exit;
}
