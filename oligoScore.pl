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
our ($VCF_dbSNP, $REF_GRCh37, $REF_GRCh37_SNP, $REF_GRCh38, $REF_GRCh38_SNP, $BLAT);
my ($ftarget, $fkey,$outdir, $NoFilter, $ComeFromRefer);
my $fref = $REF_GRCh37;
my $fref_snp;
my $fdatabases = $REF_GRCh37;
my $step = 1;
my $para_num = 10;
my $stm = 45;
my $opt_tm=60;
my $opt_tm_probe=70;
my $min_len=18;
my $max_len=28;
my $min_len_probe=18;
my $max_len_probe=36;
my $scale_len=1;
my $pcr_size=1000;
my $min_eff=0.01;
my $max_prodn=50;
my $pnum = 20; ## position num for one oligo, its candidate oligos num is roughly: pnum*(maxl-minl)/scalel.
my $choose_num = 10; ## for a primer, at least $choose_num primers can be selected as its pair with the best dis.
my $rfloat = 0.2;
my $dis_aver = 500;
my $dis_range="100,150,70,200"; ##distance between oligos range(best_min, best_max, min, max), that is product size range when face-to-face
my $pos_range; ## pos range(best_min, best_max, min, max)
my $type = "face-to-face";
my $plex = "SinglePlex";
my $ctype = "Single";
my $onum = 3;
my ($homology_check, $multiplex_check);
my $probe;
my ($regions, $regions_rev);
my $thread;
GetOptions(
				"help|?" =>\&USAGE,
				"it:s"=>\$ftarget,
				"ComeFromRefer:s"=>\$ComeFromRefer,
				"ir:s"=>\$fref,
				"is:s"=>\$fref_snp,
				"id:s"=>\$fdatabases,
				"p:s"=>\$fkey,
				"Probe:s"=>\$probe,
				"Homology_check:s"=>\$homology_check,
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
				"ptype:s"=>\$plex,
				"ctype:s"=>\$ctype,
				"ds:s"=>\$dis_aver,
				"rf:s"=>\$rfloat,
				"on:s"=>\$onum,

				"stm:s"=>\$stm,
				"pnum:s"=>\$pnum,
				"size:s"=>\$pcr_size,
				"meff:s"=>\$min_eff,
				"mpro:s"=>\$max_prodn,
				"para:s"=>\$para_num,
				"thrd:s"=>\$thread,
				"step:s"=>\$step,
				"od:s"=>\$outdir,
				) or &USAGE;
&USAGE unless ($ftarget and $fdatabases and $fkey);

$outdir||="./";
`mkdir $outdir`	unless (-d $outdir);
$outdir=AbsolutePath("dir",$outdir);
$ftarget=AbsolutePath("file", $ftarget);
$fref=AbsolutePath("file", $fref);
my @fdb=split /,/, $fdatabases;
for(my $i=0; $i<@fdb; $i++){
	$fdb[$i]=AbsolutePath("file",$fdb[$i]);
}
$fdatabases=join(",", @fdb);

if(!defined $fref_snp){
	my $fref_name = basename($fref);
	if($fref_name eq basename($REF_GRCh37)){
		$fref_snp = $REF_GRCh37_SNP;
	}elsif($fref_name eq basename($REF_GRCh38)){
		$fref_snp = $REF_GRCh38_SNP;
	}
}


my ($score_dis, $score_pos)=(10,10);
my $sh;
open($sh, ">$outdir/$fkey.sh") or die $!;

### get file type
my $ftype;
my ($ftemplate, $ftemplate_snp);
my $head = `head -1 $ftarget`;
chomp $head;
if($head=~/^>/){
	$ftype = "fasta";
	$ftemplate = $ftarget;
	if(!defined $ComeFromRefer){
		$fdatabases.=",".$ftemplate;
	}
}else{
	my @unit = split /\s+/, $head;
	my $cnum = scalar @unit;
	if($cnum == 1){
		if($head=~/rs\d+/){
			&Run("perl $Bin/extract_by_id.pl -i $VCF_dbSNP -l $ftarget -k $fkey.SNP -c 2 -od $outdir", $sh);
			&Run("less $outdir/$fkey.SNP.extract|awk '{print \"chr\"\$_}' > $outdir/$fkey.target.txt", $sh);
			$ftarget = "$outdir/$fkey.target.txt";
			$ftype = "SNP";
		}else{
			die "Wrong input file!\n";	
		}
	}elsif($cnum==2){
		if($unit[1]=~/^[a-zA-Z\d]+$/){
			$ftype = "Primer";
		}else{
			die "Wrong input file!\n";
		}
	}elsif($cnum==3 || $cnum==4){
		if($unit[1]=~/^\d+$/ && $unit[2]=~/^\d+$/ && $unit[2]-$unit[1]>80){
			$ftype = "Bed";
		}else{
			die "Wrong input file!\n";
		}
	}elsif($cnum==5 || $cnum==8){## vcffile is 8col
		if($unit[1]=~/^\d+$/ && $unit[2]!~/^\d+$/){ ## SNP
			$ftype = "SNP";
		}else{
			die "Wrong input file!\n";
		}
	}else{
		die "Wrong input file!\n";
	}
}

### get template file
my @rdiss = split /,/, $dis_range;
if($ftype eq "SNP"){
	my $extend_len;
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
}elsif($ftype eq "Bed"){
	open(I, $ftarget) or die $!;
	$ftemplate = "$outdir/$fkey.template.fa";
	$ftemplate_snp = "$outdir/$fkey.template_snp.fa";
	open(T, ">$ftemplate") or die $!;
	open(TS, ">$ftemplate_snp") or die $!;
	while(<I>){
		chomp;
		next if(/^$/);
		my ($c, $s, $e, $id)=split /\s+/, $_;
		if(!defined $id){
			$id=join(",", $c, $s, $e);
		}
		my $info=`samtools faidx $fref $c:$s-$e`;
		my (undef, @line)=split /\n/, $info;
		print T ">$id\n";
		print T join("", @line),"\n";

		if(defined $fref_snp){
			my $info=`samtools faidx $fref_snp $c:$s-$e`;
			my (undef, @line)=split /\n/, $info;
			print TS ">$id\n";
			print TS join("", @line),"\n";
		}else{
			print TS ">$id\n";
			print TS join("", @line),"\n";
		}
	}
	close(T);
	close(TS);
}elsif($ftype eq "Primer"){
	my $range = join(",", $rdiss[-2], $rdiss[-1]);
	$thread=defined $thread? $thread: 10;
	&Run("perl $Bin/primer_evaluation.pl -io $ftarget -id $fdatabases -k $fkey -tp $type -ep $plex --AllEvalue --OutAllProduct -td $thread -sz $pcr_size -tm $opt_tm -tmb $opt_tm_probe -rd $range -me $min_eff -mp $max_prodn -od $outdir", $sh);
	exit();
}


### homology check
if($step==1){
	if(defined $homology_check){
		&Run("perl $Bin/homology_check.pl -it $ftemplate -ir $fdatabases -k $fkey -od $outdir/homology_check", $sh);
	}
	$step++;
}

### oligo design
my $odir="$outdir/design/";
if($step==2){
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
	my $dcmd = "perl $Bin/oligo_design.pl -i $ftemplate -d $fdatabases -k $fkey -ptype $type -opttm $opt_tm -rlen $range_len -regions $rregion -stm $stm -para $para_num -od $odir";
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
	### probe score
	if(defined $probe){
		my $cmd = "perl $Bin/probe_score.pl -i $odir/$fkey.oligo.evaluation.out -k $fkey -minl $min_len_probe -maxl $max_len_probe -opttm $opt_tm_probe -od $outdir";
		if(defined $NoFilter){
			$cmd .= " --NoFilter";
		}
		&Run($cmd, $sh);
	}
	
	### primer score and select
	my $score_dis_range = $score_dis.",".$dis_range;
	my $cmd = "perl $Bin/primer_score.pl -io $odir/$fkey.oligo.evaluation.out -it $ftemplate -ib $odir/$fkey.oligo.bound.info -k $fkey -tp $type -minl $min_len -maxl $max_len -opttm $opt_tm -PCRsize $pcr_size -rd=$dis_range -ct $ctype -mine $min_eff -maxp $max_prodn -od $outdir";
	if(defined $NoFilter){
		$cmd .= " --NoFilter";
	}
	if(defined $probe){
		$cmd .= " -ip $outdir/$fkey.probe.score";
	}
	if(defined $pos_range){
		$cmd .= " -rp=$pos_range";
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
	### primers multiplex check
	if($plex eq "MultiPlex"){
		my $range = join(",", $rdiss[-2], $rdiss[-1]);
		&Run("perl $Bin/primer_evaluation.pl -io $outdir/$fkey.final.result -ib $outdir/$fkey.final.bound.info -k $fkey -tp $type -ep MultiPlex -mp $max_prodn -AllEvalue -sz $pcr_size -tm $opt_tm -tmb $opt_tm_probe -rd=$range -me $min_eff -od $outdir");
		&Run("perl $Bin/cross_dimer_check.pl -i $outdir/$fkey.final.result -k $fkey.final -od $outdir", $sh);
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
	my ($ftemplate, $fdatabases, $outdir, $fkey)=@_;
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
		my $step=int(($maxp-$minp)/$pnum+0.5);
		$step=$step>=1? $step: 1;
		push @region, ($minp, $maxp, $step, "R");
		if(defined $probe){
			push @region, (1,20,1, "F"); ## oligos for probe
		}
		my $ori="R";
		if($type eq "face-to-face"){
			$min = $mind-$maxp-2*$alen;
			$max = $maxd-$minp-2*$alen;
		}elsif($type eq "back-to-back"){
			$min = $minp+$alen-$maxd;
			$max = $maxp+$alen-$mind;
			$ori="F";
		}elsif($type eq "Nested"){
			$min = $minp+$alen-$maxd-$alen; #start on template of primer 
			$max = $maxp+$alen-$mind-$alen;
		}else{
			die "Wrong type when defined -srpos, must be one of (face-to-face, back-to-back, Nested)!\n";
		}
		if($ori eq "R"){
			$min=$min<$maxp? $maxp: $min;
		}
		$step=int(($max-$min)/$pnum+0.5);
		$step=$step>=1? $step: 1;
		push @region, ($min, $max, $step, $ori);
		if($step > $step0*2){ ## check step is small enough to keep $choose_num primers to be selected as its pairs for one primer
			print "Warn: Step size $step is too big to choose primers as its pairs within optimal distance $bmind-$bmaxd! You can narrow -rpos $pos_range, or magnify -rdis $dis_range!\n";
		}
	}else{
		if($ctype eq "Single"){## usually is generic:Region
			$min = 1;
			$max = $min+$step0*$pnum;
		}elsif($ctype eq "Full-covered"){
			$min = 1;
			$max="Len"; ## total length of template
		}
		push @region, ($min, $max, $step0, "FR");
	}
	my $range_region = join(",", @region);
	return ($range_region);
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

      face-to-face: P1 |--->              dis_range:100,150,70,200(qPCR)
        (Region)              <---| P2    

      back-to-back:  x <---| P2           dis_range:5,10,0,15
                      P1 |---> x          (Overlap between p1 and p2: dis > 0)

      back-to-back:  x <---| P2           dis_range:-50,-40,-60,-30
                           P1 |---> x     (No overlap between p1 and p2: dis < 0)

            Nested: P2 --->|   x          dis_range(P2-P1):-15,-10,-30,-5
                      P1 --->| x                (dis < 0)

Usage:
  Options:
  -it        <file>   Input target file(SNP file or template fasta file with no non-ATCGatcg), forced
   --ComeFromRefer    Sequences in target file(-it) come from reference file(-ir) when -it is fasta file, optional
  -ir        <file>   Input reference file to extract template sequence of SNP, needed when target file(-it) is SNP file, [$fref]
  -is        <file>   Input reference file containing snps to check SNP of oligos when -it is SNP file, optional
  -id       <files>   Input database files separated by "," to check specificity, [$fdatabases]
  -p         <str>    prefix of output file, forced
  --Probe             Design probe when -type "face-to-face", optional
  --NoFilter          Not filter any oligos
  --Homology_check    Check homologous sequence of template sequence when design for NGS primers, optional

  ### oligo parameters
  -opttm    <int>     optimal tm of primer, [$opt_tm]
  -opttmp   <int>     optimal tm of probe, [$opt_tm_probe]
  -minl     <int>     minimum length of primer, [$min_len]
  -maxl     <int>     maximum length of primer, [$max_len]
  -minlp    <int>     minimum length of probe, [$min_len_probe]
  -maxlp    <int>     maximum length of probe, [$max_len_probe]
  -scalel   <str>     candidate oligo length scale, [$scale_len]
  -rpos     <str>     position range, distance of p1 to the detected site, (opt_min, opt_max, min, max) separted by ",", must be given when -it is SNP file
  -rdis     <str>     distance range between pair primers, that is product size range when -type is "face-to-face", (opt_min, opt_max, min, max) separted by ",", [$dis_range]
  -regions  <str>     interested regions of candidate primers walking on template, format is "start,end,scale,fr,start2,end2,scale2,fr2...", if not given, will caculate automatically, fr:F/R/FR, optional

  ### design parameters
  -type   <str>     primer type, "face-to-face", "back-to-back", "Nested", ["face-to-face"]
  -ptype  <str>     plex type, "SinglePlex" or "MultiPlex", [$plex]
  -ctype  <str>     primer covered type, "Single" or "Full-covered", ["Single"]
     -ds  <int>     average distance between adjacent primers when -ctype "Full-covered", [500]
     -rf  <float>   ratio of distance between adjacent primers can float when -ctype "Full-covered", [0.2]
     -on  <int>     output num when -ctype "Single",[$onum]
  -pnum   <int>     position num of candidate oligos, [$pnum]
 
  ### specificity parameters
  -stm    <int>     min tm to be High_tm in specifity, [$stm]
  -size   <int>     max PCR size, [$pcr_size]
  -mpro   <int>     maximum products number to be caculated, to reduce running time. [$max_prodn]
  -meff   <float>   min efficiency to consider a product, [$min_eff]

  ### run parameters
  -para   <int>      parallel num, [$para_num]
  -thrd   <int>      thread in bwa, [$thread]
  -step   <int>      step, [$step]
                     1: homology check
                     2: creat and evalue candidate oligos 
                     3: score for probe and primer
                     4: multiplex check
  -od     <dir>      Dir of output file, default ./
  -h                Help

USAGE
    print $usage;
    exit;
}

