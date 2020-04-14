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
my $min_len = 20;
my $max_len = 35;
my $range_dis = "5,60,1";
my $min_tm = 65;
my $max_tm = 75;
my $min_gc = 0;
my $max_gc = 0.8;
my $step = 1;
my $recov;
my $para_num = 10;
my $stm = 45;
my $scale_len = 2;
my $scale_dis = 2;
my $line_num = 3;

my $rfloat = 0.2;
my $dis_aver = 500;
my $score_dis_range="10,10,15,1,30"; ## dis score, pair dis range(best_min, best_max, min, max)
my $score_pos_range="10,10,15,1,40"; ## pos score, pos range(best_min, best_max, min, max)
my $dtype = "Single";
my $onum = 3;
my $face_to_face;
my $back_to_back;
my $dimer_check;
GetOptions(
				"help|?" =>\&USAGE,
				"it:s"=>\$ftarget,
				"ir:s"=>\$fref,
				"k:s"=>\$fkey,
				"recov:s"=>\$recov,
				"dimer_check:s"=>\$dimer_check,
				## primer design
				"mintm:s"=>\$min_tm,
				"maxtm:s"=>\$max_tm,
				"mingc:s"=>\$min_gc,
				"maxgc:s"=>\$max_gc,
				"minl:s"=>\$min_len,
				"maxl:s"=>\$max_len,
				"scalel:s"=>\$scale_len,
				"rdis:s"=>\$range_dis,
				"stm:s"=>\$stm,
				"lnum:s"=>\$line_num,

				## double select
				"face-to-face:s"=>\$face_to_face,
				"back-to-back:s"=>\$back_to_back,
				"rd:s"=>\$score_dis_range,
				"rp:s"=>\$score_pos_range,
				"dt:s"=>\$dtype,
				"ds:s"=>\$dis_aver,
				"rf:s"=>\$rfloat,
				"on:s"=>\$onum,

				"para:s"=>\$para_num,
				"od:s"=>\$outdir,
				) or &USAGE;
&USAGE unless ($ftarget and $fref and $fkey);

$outdir||="./";
`mkdir $outdir`	unless (-d $outdir);
$outdir=AbsolutePath("dir",$outdir);

my $sh;
open($sh, ">$outdir/$fkey.sh") or die $!;

### get template file
my $ftemplate;
my $type;
my $head = `head -1 $ftarget`;
chomp $head;
if($head=~/^>/){
	$type = "fasta";
	$ftemplate = $ftarget;
}else{
	$type = "SNP";
}
if($type eq "SNP"){
	&Run("perl $Bin/get_template.pl -i $ftarget -r $fref -k $fkey -od $outdir/ > $outdir/$fkey.get_template.log 2>&1", $sh);
	$ftemplate = "$outdir/$fkey.template.fa";
}

### primer design
&Run("perl $Bin/primer_design.pl -i $ftemplate -r $fref -minl $min_len -maxl $max_len -mintm $min_tm -maxtm $max_tm -mingc $min_gc -maxgc $max_gc -scalel $scale_len -rdis $range_dis -lnum $line_num -stm $stm -para $para_num -k $fkey -od $outdir/design > $outdir/design.log 2>&1", $sh);
if(defined $recov){
	&Run("perl $Bin/primer_design.pl -i $ftemplate -r $fref --recov -minl $min_len -maxl $max_len -mintm $min_tm -maxtm $max_tm -mingc $min_gc -maxgc $max_gc -scalel $scale_len -rdis $range_dis -lnum $line_num -stm $stm -para $para_num -k $fkey\_rev -od $outdir/design_rev > $outdir/design_rev.log 2>&1", $sh);
}


### select primer pair
my $cmd = "perl $Bin/double_primer_select.pl -i $outdir/design/$fkey.primer.score -it $ftemplate -k $fkey -rd $score_dis_range -rp $score_pos_range  -od $outdir --dt $dtype";
if(defined $recov){
	$cmd .= " -ir $outdir/design_rev/$fkey\_rev.primer.score";
}
if(defined $face_to_face){
	$cmd .= " --face-to-face";
}
if(defined $back_to_back){
	$cmd .= " --back-to-back";
}
if($dtype eq "Region"){
	$cmd .= " -ds $dis_aver -rf $rfloat";
}else{
	$cmd .= " -on $onum";
}
&Run($cmd, $sh);

### primers dimer check
if(defined $dimer_check){
	&Run("perl $Bin/primers_dimer_check.pl -i $outdir/$fkey.final.primer -k $fkey -od $outdir", $sh);
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

    ############## PrimerScore pipeline  ################

    ###input file(target spot file) format:
     format 1(note: "ID=xxx;" must be given):
        chr7    107350577       107350577       A       G       ID=SLC26A4_2168A_G;GENE=SLC26A4;STRAND=+;CDS=c.2168A>G;AA=p.H723R;NM=NM_000441.1
        chr13   20763421        20763422        AT      -       ID=GJB2_299delAT;GENE=GJB2;STRAND=-;CDS=c.299_300del;AA=p.H100fs;NM=NM_004004.5
        chr19   45412079        45412079        C       T       ID=rs7412;GENE=Unkown
        chr19   45411941        45411941        T       C       ID=rs429358;GENE=Unkown
     format 2(first 5 columns of vcf format): 
        chr1    10247   rs796996180     T       C
        chr1    10248   rs148908337     A       T
        chr1    10249   rs774211241     AAC     A
        chr1    10250   rs199706086     A       C

    ### -rd: distance range of pair primers, (best_min, best_max, min, max) separted by ",", example:
       #P1 is F and just the one which is used to count distance to target(-rp: position range);
       semi-Nested: P2 --->|   x         dis_range(P2-P1):10,15,1,30
                      P1 --->| x
      
      face-to-face: --->| P1             dis_range(P1+P2):180,220,150,250
       (Single)              x  P2 |<--- 
    
      face-to-face: P1 --->|             dis_range(P2-P1):-550,-450,-800,-200
        (Region)             <---| P2
    
      back-to-back:   <---| P2           dis_range(P2-P1):20,30,10,50
                      P1 --->| x


Usage:
  Options:
  -it  <file>   Input target file(target spot file or template fasta file), forced
  -ir  <file>   Input ref file, forced
  -k   <str>    Key of output file, forced
  --recov          recov to creat reverse primer of template when back-to-back and face-to-face(Region), optional
  --dimer_check    dimer check among final designed primers, optional

  ### design parameters
  -mintm   <int>      min tm, [$min_tm]
  -maxtm   <int>      max tm, [$max_tm]
  -minl    <int>      min primer len, [$min_len]
  -maxl    <int>      max primer len, [$max_len]
  -mingc   <float>    min gc, [$min_gc]
  -maxgc   <float>    max gc, [$max_gc]
  -scalel  <int>      scale len(step size) [$scale_len]
  -rdis    <str>      template regions and scale(step size) to generate candidate primer(start,end,scale), start <= end, count from right to left and count from 0, separated by ",", [$range_dis]
                      Example:
                           sanger sequence primer: 100,150,2,400,500,5
                           ARMS PCR primer: 1,1,1,80,180,2
  -lnum    <int>      line num in one separated file, [$line_num]
  -stm     <int>      min tm to be High_tm in specifity, [$stm]

  ### select parameters
  -rd     <str>     distance range of pair primers, (dis_score, best_min, best_max, min, max) separted by ",", [$score_dis_range]
  -rp     <str>     position range(distance to the detected site) when -dt is Single, [$score_pos_range]
  --face-to-face    design face-to-face primer
  --back-to-back    design back-to-back primer
  --dt    <str>     primer design type, "Single" or "Region", ["Single"]
     -ds  <int>     distance when -dt Region, [500]
     -rf  <float>   distance float ratio when -dt Region, [0.2]
     -on  <int>     output num when -dt is Single,[$onum]

  -para  <int>      parallel num, [$para_num]
  -od    <dir>      Dir of output file, default ./
  -h                Help

USAGE
    print $usage;
    exit;
}

