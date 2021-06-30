#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $BEGIN_TIME=time();
my $version="1.0.0";
require "$Bin/path.pm";
require "$Bin/common.pm";
require "$Bin/math.pm";
#######################################################################################

# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
our ($VCF_dbSNP, $REF_HG19, $REF_HG19_SNP, $BLAT);
my ($ftarget, $fkey,$outdir, $NoFilter, $ComeFromRefer);
my $fref = $REF_HG19;
my $fref_snp;
my $fdatabase = $REF_HG19;
my $step = 1;
my $para_num = 10;
my $stm = 45;
my $opt_tm=60;
my $opt_tm_probe=70;
my $min_len=18;
my $max_len=28;
my $min_len_probe=18;
my $max_len_probe=36;
my $scale_len=2;
my $pcr_size=600;
my $pnum = 100; ## position num for one oligo, its candidate oligos num is roughly: pnum*(maxl-minl)/scalel.
my $choose_num = 10; ## for a primer, at least $choose_num primers can be selected as its pair with the best dis.
my $rfloat = 0.2;
my $dis_aver = 500;
my $dis_range="100,150,70,200"; ##distance between oligos range(best_min, best_max, min, max), that is product size range when face-to-face
my $pos_range; ## pos range(best_min, best_max, min, max)
my $type = "face-to-face";
my $ctype = "Single";
my $onum = 3;
my ($dimer_check, $homology_check, $SNP_check);
my $probe;
my ($regions, $regions_rev);
GetOptions(
				"help|?" =>\&USAGE,
				"it:s"=>\$ftarget,
				"ComeFromRefer:s"=>\$ComeFromRefer,
				"ir:s"=>\$fref,
				"is:s"=>\$fref_snp,
				"id:s"=>\$fdatabase,
				"p:s"=>\$fkey,
				"dimer_check:s"=>\$dimer_check,
				"homology_check:s"=>\$homology_check,
				"NoFilter:s"=>\$NoFilter,

				## oligo design
				"opttm:s"=>\$opt_tm,
				"opttmp:s"=>\$opt_tm_probe,
				"minl:s"=>\$min_len,
				"maxl:s"=>\$max_len,
				"minlp:s"=>\$min_len_probe,
				"maxlp:s"=>\$max_len_probe,
				"scalel:s"=>\$scale_len,

				"regions:s"=>\$regions,
				"rdis:s"=>\$dis_range,
				"rpos:s"=>\$pos_range,

				"type:s"=>\$type,
				"ctype:s"=>\$ctype,
				"probe:s"=>\$probe,
				"ds:s"=>\$dis_aver,
				"rf:s"=>\$rfloat,
				"on:s"=>\$onum,

				"stm:s"=>\$stm,
				"pnum:s"=>\$pnum,
				"size:s"=>\$pcr_size,
				"para:s"=>\$para_num,
				"step:s"=>\$step,
				"od:s"=>\$outdir,
				) or &USAGE;
&USAGE unless ($ftarget and $fdatabase and $fkey);
$outdir||="./";
`mkdir $outdir`	unless (-d $outdir);
$outdir=AbsolutePath("dir",$outdir);
$ftarget=AbsolutePath("file", $ftarget);
$fref=AbsolutePath("file", $fref);
$fdatabase=AbsolutePath("file", $fdatabase);

my ($score_dis, $score_pos)=(10,10);
my $sh;
open($sh, ">$outdir/$fkey.sh") or die $!;

### get template file
my ($ftemplate, $ftemplate_snp);
my $ftype;
my $head = `head -1 $ftarget`;
chomp $head;
if($head=~/^>/){
	$ftype = "fasta";
	$ftemplate = $ftarget;
	if(!defined $ComeFromRefer){
		$fdatabase.=",".$ftemplate;
	}
}else{
	$ftype = "SNP";
	my @unit = split /\s+/, $head;
	my $cnum = scalar @unit;
	if($cnum == 1 && $head=~/rs\d+/){
		&Run("perl $Bin/extract_by_id.pl -i $VCF_dbSNP -l $ftarget -k $fkey.SNP -c 2 -od $outdir", $sh);
		&Run("less $outdir/$fkey.SNP.extract|awk '{print \"chr\"\$_}' > $outdir/$fkey.target.txt", $sh);
		$ftarget = "$outdir/$fkey.target.txt";
	}elsif($cnum < 5){
		die "Wrong input target file: $ftarget\n";
	}
}

if($ftype eq "SNP"){
	my $extend_len;
	my @rdiss = split /,/, $dis_range;
	$extend_len = (int(($rdiss[-1]+$max_len+10)/100)+1) * 100;
	
	if(!defined $pos_range){
		die "Wrong: -rpos must be given when -it SNP file!\n";
	}
	my ($optmin, $optmax, $min, $max) = split /,/, $pos_range;
	my $md = int (($max+$optmax)/2) - $optmin;
	my $cmd = "perl $Bin/get_template.pl -i $ftarget -r $fref -k $fkey -et $extend_len -md $md -od $outdir/ --dieC";
	if(defined $fref_snp){
		$cmd .= " -s $fref_snp";
	}
	&Run($cmd, $sh);
	$ftemplate = "$outdir/$fkey.template.fa";
	$ftemplate_snp = "$outdir/$fkey.template_snp.fa";
}


### homology check
if($step==1){
	if(defined $homology_check){
		&Run("perl $Bin/homology_check.pl -it $ftemplate -ir $fdatabase -k $fkey -od $outdir/homology_check", $sh);
	}
	$step++;
}

### oligo design
if($step==2){
	my $odir="$outdir/design/";
	my $rregion;
	if(defined $regions){
		$rregion=$regions;
	}else{
		$rregion=&caculate_rregion($dis_range, $pos_range, $type, $min_len, $max_len, $probe);
	}
	print "region range to design oligos: ", $rregion,"\n";
	
	### design oligo
	my $range_len=join(",", $min_len, $max_len, $scale_len);
	if(defined $probe){
		my $minl=&min($min_len, $min_len_probe);
		my $maxl=&max($max_len, $max_len_probe);
		$range_len=join(",", $minl, $maxl, $scale_len);
	}
	my $dcmd = "perl $Bin/oligo_design.pl -i $ftemplate -d $fdatabase -k $fkey -ptype $type -opttm $opt_tm -rlen $range_len -regions $rregion -stm $stm -para $para_num -od $odir";
	if(defined $ftemplate_snp){
		$dcmd .= " -is $ftemplate_snp";
	}
	if(defined $probe){
		$dcmd.=" --Probe -opttmp $opt_tm_probe";
	}
	if(defined $NoFilter){
		$dcmd .= " --NoFilter";
	}
	&Run($dcmd, $sh);
	$step++;
}

if($step==3){
	my $odir="$outdir/design/";
	### probe score
	if(defined $probe){
		&Run("perl $Bin/probe_score.pl -i $odir/$fkey.oligo.evaluation.out -k $fkey -minl $min_len_probe -maxl $max_len_probe -opttm $opt_tm_probe -od $odir", $sh);
	}
	
	### primer score and select
	my $score_dis_range = $score_dis.",".$dis_range;
	my $cmd = "perl $Bin/primer_score.pl -io $odir/$fkey.oligo.evaluation.out -it $ftemplate -ib $odir/$fkey.oligo.bound.info -k $fkey -tp $type -minl $min_len -maxl $max_len -opttm $opt_tm -PCRsize $pcr_size -rd $dis_range -ct $ctype -od $outdir";
	if(defined $probe){
		$cmd .= " -ip $odir/$fkey.probe.score";
	}
	if(defined $pos_range){
		$cmd .= " -rp $pos_range";
	}
	if($ctype eq "Full-covered"){
		$cmd .= " -ds $dis_aver -rf $rfloat";
	}else{
		$cmd .= " -on $onum";
	}
	&Run($cmd, $sh);
	$step++;
}

if($step==4){
	### primers dimer check
	if(defined $dimer_check){
		&Run("perl $Bin/cross_dimer_check.pl -i $outdir/$fkey.final.result -k $fkey -od $outdir/dimer_check", $sh);
	}

	$step++;
}

#######################################################################################
print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################

# ------------------------------------------------------------------
# sub function
# ------------------------------------------------------------------
sub template_database_check{
	my ($ftemplate, $fdatabase, $outdir, $fkey)=@_;
	my $fpsl = "$outdir/$fkey.psl";
	if(!-e $fpsl){
		&Run("$BLAT $fref $ftemplate $fpsl", $sh);
	}

	open(I, $ftemplate) or die $!;
	$/=">";
	my $tnum=0;
	while(<I>){
		chomp;
		next if(/^$/);
		$tnum++;
	}
	close(I);


	open(I, "$outdir/$fkey.psl") or die $!;
#psl v1:
#149     14      0       0       1       28      1       29      +       PTS-0639-U      201     10      201     chrX    155270560       38777907        38778099        2       99,64,   10,137, 38777907,38778035,
#176     13      0       0       0       0       1       1       +       PTS-0639-U      201     12      201     chr8    146364022       49220714        49220904        2       124,65,  12,136, 49220714,49220839,
#172     19      0       0       0       0       0       0       +       PTS-0639-U      201     10      201     chr8    146364022       144837574       144837765       1       191,     10,     144837574,
#173     18      0       0       0       0       0       0       +       PTS-0639-U      201     10      201     chr7    159138663       98139617        98139808        1       191,     10,     98139617,

#psl v2:
#501     0       0       0       0       0       0       0       +       chr11   135006516       1269880 1270381 chr11:1270382-U 501     0       501     1       501,    1269880,
#499     2       0       0       0       0       0       0       +       chr11   135006516       1266109 1266610 chr11:1270382-U 501     0       501     1       501,    1266109,
#490     11      0       0       0       0       0       0       +       chr11   135006516       1268209 1268710 chr11:1270382-U 501     0       501     1       501,    1268209,
#487     14      0       0       0       0       0       0       +       chr11   135006516       1264438 1264939 chr11:1270382-U 501     0       501     1       501,    1264438,
	my %is_in;
	$/="\n";
	while(<I>){
		chomp;
		next if($_!~/^\d/);
		my ($match, $mismatch, $gap1,$gaplen1, $gap2, $gaplen2, $ori, $pid, $len, $ps, $pe, $chr, $s, $e)=(split /\s+/, $_)[0,1,4,5,6,7,8,9,10, 11,12,13,15,16]; #psl v1
		if($match/$len>0.99){
			$is_in{$pid}=1;
		}
	}
	close(I);
	my $is_in=0;
	if(scalar keys %is_in==$tnum){
		$is_in=1;
	}
	return $is_in;
}

sub caculate_rregion{
	my ($dis_range, $pos_range, $type, $minl, $maxl, $probe)=@_;
	## caculate -rregion parameter
	my ($bmind, $bmaxd, $mind, $maxd)=split /,/, $dis_range;
	my @region; ##  two-dimensional array, @{region[1]} is regions of revcom template sequence
	my $alen = int(($minl+$maxl)/2);
	my ($min, $max);
	my $step0 = int(($bmaxd-$bmind)/$choose_num+0.5); ## max step for distance range $dis_range
	if(defined $pos_range){
		my ($bminp, $bmaxp, $minp, $maxp)=split /,/, $pos_range;
		push @region, ($minp, $maxp, int(($maxp-$minp)/$pnum)+1, "R");
		if(defined $probe){
			push @region, (1,20,1, "F"); ## oligos for probe
		}
		if($type eq "face-to-face"){
			$min = $mind-$maxp-2*$alen;
			$max = $maxd-$minp-2*$alen;
		}elsif($type eq "back-to-back"){
			$min = $minp+$alen-$maxd;
			$max = $maxp+$alen-$mind;
		}elsif($type eq "Nested"){
			$min = $maxp+$mind;
			$max = $maxp+$maxd;
		}else{
			die "Wrong type when defined -srpos, must be one of (face-to-face, back-to-back, Nested)!\n";
		}
		my $posnum = int(($max-$min+1)/$step0);
		if($posnum > $pnum*2){ ## check step is small enough to keep $choose_num primers to be selected as its pairs for one primer
			die "Step size $step0 in region $min-$max produces too many primers, which will take too long to design! Please narrow -rpos $pos_range, or magnify -rdis $dis_range!\n";
		}
		$min=$min<$maxp? $maxp: $min;
		push @region, ($min, $max, int(($max-$min)/$pnum)+1, "R");
	}else{
		if($ctype eq "Single"){## usually is generic:Region
			$min = 1;
			$max = $min+$step0*$pnum;
		}elsif($ctype eq "Full-covered"){
			$min = 1;
			$max = $min+$step0*$pnum*10; ## region is also limited when Full-covered considering running time
			print "Warn: Region of templates to design oligos is $min-$max, it will not cover the whole templates if templates length are more than $max, then you can magnify -pnum $pnum or -rdis $dis_range!\n";
		}
		push @region, ($min, $max, int(($max-$min)/$pnum+0.5), "FR");
	}
	my $range_region = join(",", @region);
	return ($range_region);
}

sub Run{
    my ($cmd, $sh, $nodie)=@_;
	if(defined $sh){
		print $sh $cmd,"\n";
		print "###", $cmd,"\n";
	}
    my $ret = system($cmd);
    if (!defined $nodie && $ret) {
        die "Run $cmd failed!\n";
    }
}

sub USAGE {#
	my $usage=<<"USAGE";
Program:
Version: $version
Contact:zeng huaping<huaping.zeng\@genetalks.com> 

    ############## PrimerScore pipeline  ################

    ###input file(target spot file) format:
     format 1(rdID list):
		 rs144701072
		 rs10794507
		 rs2245135
		 rs13127915
		 rs4620658
     format 2(first 5 columns of vcf format):
		 chr14    94847286     rs121912714	T              A
		 chr15    44859635     rs312262781	CTCAA          C
		 chr17    40695468     rs104894596	C              T
		 chr12    102855239    PAH_603T_G	A              C
		 chr12    102912801    PAH_158G_A	C              T

    ###-rdis: distance range of pair primers, (best_min, best_max, min, max) separted by ",", example:
   		
      face-to-face: |---> P1 x            dis_range:100,150,70,200(qPCR); 530,570,500,600(Sanger)
         (SNP)               x  P2 <---|  

      face-to-face: P1 |--->        x     dis_range:100,150,70,200(qPCR)
        (Region)     x        <---| P2    

      back-to-back:  x <---| P2           dis_range:5,10,0,15
                      P1 |---> x          (Overlap between p1 and p2: dis > 0)

      back-to-back:  x <---| P2           dis_range:-50,-40,-60,-30
                           P1 |---> x     (No overlap between p1 and p2: dis < 0)

            Nested: P2 --->|   x          dis_range(P2-P1):10,15,5,30
                      P1 --->| x

Usage:
  Options:
  -it        <file>   Input target file(SNP file or template fasta file), forced
   --ComeFromRefer    Sequences in target file(-it) come from reference file(-ir) when -it is fasta file, optional
  -ir        <file>   Input reference file to extract template sequence of SNP, needed when target file(-it) is SNP file, [$fref]
  -is        <file>   Input reference file containing snps to check SNP of oligos when -it is SNP file, optional
  -id        <file>   Input database file to check specificity, [$fdatabase] 
  -p         <str>    prefix of output file, forced
  --probe             design probe when -type "face-to-face", optional
  --NoFilter          Not filter any primers
  --homology_check    check homologous sequence of template sequence when design for NGS primers, optional
  --dimer_check       check cross dimers among selected primers, optional

  ### design parameters
  -opttm    <int>     optimal tm of primer, [$opt_tm]
  -opttmp   <int>     optimal tm of probe, [$opt_tm_probe]
  -minl     <int>     minimum length of primer, [$min_len]
  -maxl     <int>     maximum length of primer, [$max_len]
  -minlp    <int>     minimum length of probe, [$min_len_probe]
  -maxlp    <int>     maximum length of probe, [$max_len_probe]
  -scalel   <str>     candidate oligo length scale, [$scale_len]
  -rpos     <str>     position range, distance of p1 to the detected site, (opt_min, opt_max, min, max) separted by ",", must be given when -it is SNP file
  -rdis     <str>     distance range between pair primers, that is product size range when -type is "face-to-face", (opt_min, opt_max, min, max) separted by ",", [$dis_range]
  -regions  <str>     interested regions of candidate primers walking on template, format is "start,end,scale,start2,end2,scale2...", if not given, will caculate automatically, optional

  ### 
  -type   <str>     primer type, "face-to-face", "back-to-back", "Nested", ["face-to-face"]
  -ctype  <str>     primer covered type, "Single" or "Full-covered", ["Single"]
     -ds  <int>     average distance between adjacent primers when -ctype "Full-covered", [500]
     -rf  <float>   ratio of distance between adjacent primers can float when -ctype "Full-covered", [0.2]
     -on  <int>     output num when -ctype "Single",[$onum]

  -stm    <int>      min tm to be High_tm in specifity, [$stm]
  -pnum   <int>      position num of candidate oligos, [$pnum]
  -size   <int>      max PCR size, [$pcr_size]
  -para   <int>      parallel num, [$para_num]
  -step   <int>      step, [$step]
                     1: homology check
					 2: primer pair design 
					 3: primer pair specificity re-evalue, cross-dimer check
					 4: probe design
  -od     <dir>      Dir of output file, default ./
  -h                Help

USAGE
    print $usage;
    exit;
}

