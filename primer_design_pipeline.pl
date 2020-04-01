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
my ($ftarget, $fref,$fkey,$outdir);
my $para_num=10;
my $step =1;
my $min_tm = 65;
my $max_tm = 75;
GetOptions(
				"help|?" =>\&USAGE,
				"it:s"=>\$ftarget,
				"ir:s"=>\$fref,
				"k:s"=>\$fkey,
				"mintm:s"=>\$min_tm,
				"maxtm:s"=>\$max_tm,
				"s:s"=>\$step,
				"p:s"=>\$para_num,
				"od:s"=>\$outdir,
				) or &USAGE;
&USAGE unless ($ftarget and $fref and $fkey);

$outdir||="./";
`mkdir $outdir`	unless (-d $outdir);
$outdir=AbsolutePath("dir",$outdir);

my $sh;
open($sh, ">$outdir/$fkey.sh") or die $!;
if($step ==1){
	&Run("perl $Bin/get_template_for_ctbest.pl -i $ftarget -r $fref -k $fkey -od $outdir/ > $outdir/$fkey.get_template.log 2>&1", $sh);
	$step++;
}
if($step ==2){
	&Run("perl $Bin/primer_design.pl -i $outdir/$fkey.template.fa -r $fref -mintm $min_tm -maxtm $max_tm -para $para_num -k $fkey -od $outdir/design > $outdir/design.log 2>&1", $sh);
	&Run("perl $Bin/double_primer_select.pl -i $outdir/design/$fkey.primer.score -it $outdir/$fkey.template.fa -k $fkey -od $outdir", $sh);
	&Run("perl $Bin/output_order_format.pl -i $outdir/$fkey.final.primer -k $fkey.final.primer -od $outdir", $sh);
	$step ++;
}

if($step ==3){
	&Run("perl $Bin/check_primer_and_hotspot.pl -i1 $ftarget -i2 $outdir/$fkey.final.primer.txt -k $fkey.final -od $outdir", $sh);
	## double primer specificity evalue
	my @UD=("U", "D");
	my @ABC=("A", "B", "C");
	foreach my $ud(@UD){
		foreach my $abc(@ABC){
			`mkdir "$outdir/check_ud/"` unless(-d "$outdir/check_ud/");
			&Run("less $outdir/$fkey.final.primer |perl -ne '{chomp; \@a=split; if(\$_=~/SP1/){\$TN=\$a[4];}  if(\$_=~/SP2/){ print join(\"\\t\", \$a[3],\$a[4],\$TN),\"\\n\";}}'|grep '\\-$ud\\-$abc'|less >$outdir/check_ud/$fkey.final.primer.list_$ud$abc", $sh);
			&Run("perl $Bin/primer_evaluation.pl -p $outdir/check_ud/$fkey.final.primer.list_$ud$abc -d $fref -n 2 -k $fkey.final.primer.list_$ud$abc -od $outdir/check_ud >$outdir/check_ud/$fkey.final.primer.list_$ud$abc.log 2>&1", $sh);
		}
	}
	&Run("cat >$outdir/$fkey.final.primer.list.evaluation.out $outdir/check_ud/*evaluation.out", $sh);
}



#######################################################################################
print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################

# ------------------------------------------------------------------
# sub function
# ------------------------------------------------------------------

sub Run{
    my ($cmd, $sh, $nodie)=@_;
	print $sh $cmd,"\n";
    my $ret = system($cmd);
    if (!defined $nodie && $ret) {
        die "Run $cmd failed!\n";
    }
}

sub AbsolutePath
{		#获取指定目录或文件的决定路径
		my ($type,$input) = @_;

		my $return;
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

		## design RSQ-TN primer pipeline

Usage:
  Options:
  -it  <file>   Input target file, forced
  -ir  <file>   Input ref file, forced
  -k   <str>	Key of output file, forced

  -mintm <int>		min tm, [$min_tm]
  -maxtm <int>		max tm, [$max_tm]
  -p   <int>    parallel num, [$para_num]
  -s   <int>    step, 1:get template file; 2: design; 3: evaluation, [$step]
  -od  <dir>	Dir of output file, default ./
  -h		 Help

USAGE
	print $usage;
	exit;
}

