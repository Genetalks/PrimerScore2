#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $BEGIN_TIME=time();
my $version="3.0";
require "$Bin/self_lib.pl";
#######################################################################################

# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
my ($fIn,$ftm, $frev,$fkey,$outdir);
my $rfloat = 0.2;
my $dis = 500;
my $range_dis="10,10,15,1,30"; ## dis score, pair dis range(best_min, best_max, min, max)
my $range_pos="10,10,15,1,40"; ## pos score, pos range(best_min, best_max, min, max)
my $dtype = "Single";
my $onum = 3;
my $face_to_face;
my $back_to_back;
GetOptions(
				"help|?" =>\&USAGE,
				"i:s"=>\$fIn,
				"ir:s"=>\$frev,
				"it:s"=>\$ftm,
				"k:s"=>\$fkey,
				"face-to-face:s"=>\$face_to_face,
				"back-to-back:s"=>\$back_to_back,
				"on:s"=>\$onum,
				"rd:s"=>\$range_dis,
				"rp:s"=>\$range_pos,
				"dt:s"=>\$dtype,
				"ds:s"=>\$dis,
				"rf:s"=>\$rfloat,
				"od:s"=>\$outdir,
				) or &USAGE;
&USAGE unless ($fIn and $ftm and $fkey);

$outdir||="./";
`mkdir $outdir`	unless (-d $outdir);
$outdir=AbsolutePath("dir",$outdir);

my ($sdis, @range_dis) = split /,/, $range_dis;
my ($spos, @range_pos) = $dtype eq "Single"? split /,/, $range_pos: ();
my ($sdtm, @range_dtm) = (10,0,2,0,8);

my %score;
my %pos;
my %plen;
my %info;
my %tm;
my (@explain, @title_info);
open(I, $fIn) or die $!;
while(<I>){
	chomp;
	next if(/^$/);
	if(/^##/){
		push @explain, $_;
		next;
	}elsif(/^#/){
		(undef, undef, undef, undef, undef, undef, @title_info) = split;
		next;
	}

	my @unit = split /\t/, $_;
	my ($id, $tid, $len, $pos, $tm, $score)=@unit[0,1,2,3,7,5];
	$score{$id}=$score;
	$plen{$id}=$len;
	$tm{$id}=$tm;
	$pos{$tid}{$id}=$pos;
	@{$info{$id}}=@unit[1..$#unit];
}
close(I);



my %tlength;
my %XP;
open(F, $ftm) or die $!;
$/=">";
while(<F>){
	chomp;
	next if(/^$/);
	my ($head, @seq)=split /\n/, $_;
	my ($id) = split /\t/, $head;
	my $seq = join("", @seq);
	$tlength{$id}=length $seq;

	if($head=~/XP:Z:/){
		my (undef, $xs, $xe, $xp)=split /\t/, $head;
		$xs=~s/XS:i://;
		$xe=~s/XE:i://;
		my ($strand, $chr, $s, $e)=$xp=~/XP:Z:([+-])(\S+):(\d+)-(\d+)/; #>NF2-401-U      XS:i:1  XE:i:1  XP:Z:+chr22:30038027-30038227   XD:Z:NF2,+,401,COSM4414719;
		@{$XP{$id}}=($chr, $s, $e, $strand, $xs, $xe);
	}
}
$/="\n";
close(F);

my %pos_rev;
if(defined $frev){
	open(I, $frev) or die $!;
	while(<I>){
		chomp;
		next if(/^$/);
		my @unit = split;
		my ($id, $tid,$plen, $pos, $score)=@unit[0,1,2,3,5];
		$score{$id}=$score;
	#	$unit[3]=$tlength{$tid}-$pos-$plen; ## |<- change to <-|
		$pos_rev{$tid}{$id}=$unit[3]; ## <-|
		@{$info{$id}}=@unit[1..$#unit];
	}
	close(I);
}

open(O,">$outdir/$fkey.final.primer") or die $!;
print O join("\n", @explain),"\n";
print O "#Chr\tStart\tStrand\tID\tSeq\tLen\tDis\tScorePair\tScorePair_PosDis\tScore\t",join("\t", @title_info),"\n";
foreach my $tid(sort {$a cmp $b}keys %pos){
	my @id1_sort = sort{$pos{$tid}{$a} <=> $pos{$tid}{$b}} keys %{$pos{$tid}};
	my @pos1_sort = sort{$a<=>$b} values %{$pos{$tid}};
	my @id2_sort = @id1_sort;
	my @pos2_sort = @pos1_sort;
	my ($dstart, $dend)=(1,1);
	if(exists $XP{$tid}){
		($dstart, $dend)=@{$XP{$tid}}[4..5];
	}
	
	if(defined $face_to_face && $dtype eq "Single"){ ## target
		my $tid2 = $tid;
		if($tid=~/\-U$/){$tid2=~s/\-U/\-D/;}elsif($tid=~/\-D$/){$tid2=~s/\-D/\-U/;}
		if(!exists $pos{$tid2}){
			print "Wrong: No corresponding temlate ID $tid2 for face_to_face primer, please Check!\n";
			next;
		}

		@id2_sort = sort{$pos{$tid2}{$a} <=> $pos{$tid2}{$b}} keys %{$pos{$tid2}};
		@pos2_sort = sort{$a<=>$b} values %{$pos{$tid2}};
	}
	if((defined $face_to_face && $dtype eq "Region") || defined $back_to_back){
		if(!defined $frev){
			print STDERR "Wrong: no rev primer file!\n";
			die;
		}
		@id2_sort = sort{$pos_rev{$tid}{$a} <=> $pos_rev{$tid}{$b}} keys %{$pos_rev{$tid}};
		@pos2_sort = sort{$a<=>$b} values %{$pos_rev{$tid}};
	}
	
	my %score_pair;
	my %score_pair_info;
	my %pos_pair;
	my $last_pos;
	for(my $i=0; $i<@id1_sort; $i++){
		my $id1 = $id1_sort[$i];
		my $pos = $pos1_sort[$i];
		next if($dtype eq "Single" && ($pos<$range_pos[2] || $pos>$range_pos[3]));
		my $min_pos = $pos+$range_dis[2];
		my $max_pos = $pos+$range_dis[3];
		if(defined $face_to_face && $dtype eq "Single"){
			$min_pos = $range_dis[2]-$pos;
			$max_pos = $range_dis[3]-$pos;
		}
		my $indexs = &binarySearch($min_pos, \@pos2_sort, ">=", 0, $#pos2_sort);
		my $indexe = &binarySearch($max_pos, \@pos2_sort, "<=", 0, $#pos2_sort);
		next if($indexs==-1 || $indexe==-1);

		my $score1 = $score{$id1};
		## get score_pos
		my $score_pos =0;
		if($dtype eq "Single"){
			my $pos_end = $pos+$dend-$dstart; ## pos to the end target, maybe more targets, -->  |dstart-dend| <--
			next if($pos > $range_pos[-1] || $pos_end > $range_pos[-1]); ## pos out of the pos range
			my $score_pos1 = &score_single($pos, $spos, @range_pos);
			my $score_pos2 = &score_single($pos_end, $spos, @range_pos);
			$score_pos = $score_pos1<$score_pos2? $score_pos1: $score_pos2;
		}
		for(my $j=$indexs; $j<=$indexe; $j++){
			my $id2 = $id2_sort[$j];
			my $dis = (defined $face_to_face && $dtype eq "Single")? $pos2_sort[$j]+$pos+($dend-$dstart+1): $pos2_sort[$j]-$pos; ## face-to-face Single primer, -> (p1) |dend| (p2) <- 
#			my $dis = (defined $face_to_face && $dtype eq "Single")? $pos2_sort[$j]+$pos: $pos2_sort[$j]-$pos; ## face-to-face Single primer, -> (p1) |dend| (p2) <- 
			my $score_dis = &score_single($dis, $sdis, @range_dis);
			my $dtm = abs($tm{$id1}-$tm{$id2});
			my $score_tm = &score_single($dtm, $sdtm, @range_dtm);
			my $score = $score1+$score_pos+$score_tm+$score{$id2}+$score_dis;
			$score_pair{$id1.",".$id2}=$score;
			@{$score_pair_info{$id1.",".$id2}}=($score_pos, $score_dis, $score_tm);
			$pos_pair{$id1.",".$id2}=$pos;
		}
	}
	next if(scalar keys %score_pair==0);
	my @final;
	if($dtype eq "Single"){
		my @pair_sort = sort{$score_pair{$b}<=>$score_pair{$a}} keys %score_pair;
		my %record;
		my $num = 0;
		#LBR-1096-D-28-12,LBR-1096-D-26-32
		for(my $i=0; $i<@pair_sort; $i++){
			my ($d1, $d2) = $pair_sort[$i]=~/\-(\d+),\S+\-(\d+)$/;
			my $min_dis1=100;
			my $min_dis2=100;
			foreach my $ds(keys %record){
				my ($d1r, $d2r)=split /,/, $ds;
				my $dis1 = abs($d1-$d1r);
				my $dis2 = abs($d2-$d2r);
				if($dis1<$min_dis1){
					$min_dis1 = $dis1;
				}
				if($dis2<$min_dis2){
					$min_dis2 = $dis2;
				}
			}
			next if($min_dis1<3 || $min_dis2<3 );
			push @final, $pair_sort[$i];
			$record{$d1.",".$d2}=1;
			$num++;
			last if($num==$onum);
		}
	}else{
		@final = &average($dis, $rfloat, \%pos_pair,\%score_pair,"UP"); ##my ($len, $dis, $rfloat, $apos, $ascore,$select)
	}

	##output
	my $n=0;
	foreach my $pair (@final){
		my @id = split /,/, $pair;
		my $pos = $pos_pair{$pair};
		my ($id1_new, $id2_new);
		if($dtype eq "Single"){
			my $le = chr(65+$n); ## A is 65
			$id1_new = $tid."-".$le."-SP1";
			$id2_new = $tid."-".$le."-SP2";
			$n++;
		}else{
			$id1_new = $tid."-$pos-P1";
			$id2_new = $tid."-$pos-P2";
		}
		
		my @info1 = @{$info{$id[0]}};
		my @info2 = @{$info{$id[1]}};

		my ($chr, $pos1, $strand1, $pos2, $strand2);
		my ($chrp1, $plen1, $epos1,$seq1) = @info1[0..3]; ## chrp is target id
		my ($chrp2, $plen2, $epos2, $seq2) = @info2[0..3];

#		my ($chrm, $sm, $em, $strandm)=($chrp, 1, $tlength{$chrp}, "+");
#		if(exists $XP{$chrp}){
#			($chrm,$sm, $em,$strandm)=@{$XP{$chrp}};
#		}
#		my $strandp1 = "+"; #P1 (F and nearer to target)
#		my $strandp2 = defined $frev? "-": "+"; 
#		$chr = $chrm;
#		if($strandm eq "+"){
#			$strand1 = $strandp1;
#			$strand2 = $strandp2;
#			$pos1 = $strand1 eq "+"? $em-$epos1-$plen1+1: $em-$epos1;
#			$pos2 = $strand2 eq "+"? $em-$epos2-$plen2+1: $em-$epos2;
#		}else{
#			$strand1 = $strandp1 eq "+"? "-": "+";
#			$strand2 = $strandp2 eq "+"? "-": "+";
#			$pos1 = $strand1 eq "+"? $sm+$epos1: $sm+$epos1+$plen1-1;
#			$pos2 = $strand2 eq "+"? $sm+$epos2: $sm+$epos2+$plen2-1;
#		}
		
		my ($chrm1, $sm1, $em1, $strandm1)= exists $XP{$chrp1}? @{$XP{$chrp1}}: ($chrp1, 1, $tlength{$chrp1}, "+");
		my $strandp1 = "+"; #P1 (F and nearer to target)
		($pos1, $strand1)=&get_chr_info($sm1, $em1, $strandm1, $strandp1, $plen1, $epos1);
		$chr = $chrm1;

		my ($chrm2, $sm2, $em2, $strandm2)= exists $XP{$chrp2}? @{$XP{$chrp2}}: ($chrp2, 1, $tlength{$chrp2}, "+");
		my $strandp2 = defined $frev? "-": "+"; 
		($pos2, $strand2)=&get_chr_info($sm2, $em2, $strandm2, $strandp2, $plen2, $epos2);

		print O join("\t", $chr, $pos1, $strand1, $id1_new, $seq1, $plen1, $epos1,$score_pair{$pair}, join(",",@{$score_pair_info{$pair}}), @info1[4..$#info1]),"\n";
		print O join("\t", $chr, $pos2, $strand2, $id2_new, $seq2, $plen2, $epos2,$score_pair{$pair}, join(",",@{$score_pair_info{$pair}}), @info2[4..$#info2]),"\n";
	}
}
close(O);

#######################################################################################
print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################

# ------------------------------------------------------------------
# sub function
# ------------------------------------------------------------------
sub get_chr_info{
	my ($sm, $em, $strandm, $strandp, $plen, $epos)=@_;
	my ($strand, $pos);
	if($strandm eq "+"){
		$strand = $strandp;
		$pos = $strand eq "+"? $em-$epos-$plen+1: $em-$epos;
	}else{
		$strand = $strandp eq "+"? "-": "+";
		$pos = $strand eq "+"? $sm+$epos: $sm+$epos+$plen-1;
	}

	return($pos, $strand);
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

		-pr: distance range of pair primers, (best_min, best_max, min, max) separted by ",", example:
			Rsq_Tn: 10,15,1,30
			face-to-face: -550,-450,-800,-200
			

	###v2: change order of p1 and p2: p1 is F and the one which is used to count dis to target;
            Rsq_Tn: --->|  P2        dis_range(P2-P1):10,15,1,30
                      --->|  P1

      face-to-face: --->|  P1        dis_range(P2-P1):-550,-450,-800,-200
        (Region)         <---| P2

             planD: <---|  P2        dis_range(P2-P1):20,30,10,50
                       --->| P1

		note: face-to-face result looked messy!

	###v3: 1) add face-to-face and Single primer(--> x   <--)
	       2) add ABC letter to ID when Single primer
		
	  face-to-face: --->| P1         dis_range(P1+P2):180,220,150,250
	   (Single)              x  P2 |<--- 

      back-to-back:   <---| P2           dis_range(P2-P1):20,30,10,50
                      P1 --->| x
Usage:
  Options:
  -i   <file>   Input primer score file, forced
  -it  <file>   Input template file, forced
  -ir  <file>   Input rev primer score file to design face-to-face and Region primer, optional
  -k   <str>	Key of output file, forced
  -rd  <str>	distance range of pair primers, (dis_score, best_min, best_max, min, max) separted by ",", ["10,10,15,1,30"]
  -rp  <str>	position range(distance to the detected site) when -dt is Single, ["10,10,15,1,40"]
  
  --face-to-face    design face-to-face primer
  --back-to-back    design back-to-back primer
  --dt    <str>	    primer design type, "Single" or "Region", ["Single"]
  	 -ds  <int>		distance when -dt Region, [500]
	 -rf  <float>	distance float ratio when -dt Region, [0.2]
	 -on  <int>     output num when -dt is Single,[$onum]

  -od <dir>	Dir of output file, default ./
  -h		 Help

USAGE
	print $usage;
	exit;
}


