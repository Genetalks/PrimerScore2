#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
my $BEGIN_TIME=time();
my $version="3.0";
require "$Bin/self_lib.pm";
#######################################################################################

# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
my ($fIn,$ftm, $frev,$fkey,$outdir);
my $rfloat = 0.2;
my $dis = 500;
my $range_dis="10,10,15,1,30"; ## dis score, pair dis range(best_min, best_max, min, max)
my $range_pos;
my $type = "face-to-face:SNP";
my $ctype = "Single";
my $onum = 3;
GetOptions(
				"help|?" =>\&USAGE,
				"i:s"=>\$fIn,
				"ir:s"=>\$frev,
				"it:s"=>\$ftm,
				"k:s"=>\$fkey,
				"tp:s"=>\$type,
				"on:s"=>\$onum,
				"rd:s"=>\$range_dis,
				"rp:s"=>\$range_pos,
				"ct:s"=>\$ctype,
				"ds:s"=>\$dis,
				"rf:s"=>\$rfloat,
				"od:s"=>\$outdir,
				) or &USAGE;
&USAGE unless ($fIn and $ftm and $fkey);

$outdir||="./";
`mkdir $outdir`	unless (-d $outdir);
$outdir=AbsolutePath("dir",$outdir);

my ($sdis, @range_dis) = split /,/, $range_dis;
my ($spos, @range_pos);
if(defined $range_pos){
	($spos, @range_pos) = split /,/, $range_pos;
}
my ($sdtm, @range_dtm) = (10,0,1,0,3);

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
	my ($id) = split /\s+/, $head;
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
		my ($id, $tid, $len, $pos, $tm, $score)=@unit[0,1,2,3,7,5];
		$score{$id}=$score;
		$plen{$id}=$len;
		$tm{$id}=$tm;
	#	$unit[3]=$tlength{$tid}-$pos-$plen; ## |<- change to <-|
		@{$info{$id}}=@unit[1..$#unit];
		$tid =~s/_rev$//;
		$pos_rev{$tid}{$id}=$unit[3]; ## <-|
	}
	close(I);
}

open(O,">$outdir/$fkey.final.primer") or die $!;
print O join("\n", @explain),"\n";
if($type =~/face-to-face/){
	print O "#Chr\tStart\tStrand\tID\tSeq\tLen\tDis\tProductSize\tScorePair\tScorePair_PosDisTm\tScore\t",join("\t", @title_info),"\n";
}else{
	print O "#Chr\tStart\tStrand\tID\tSeq\tLen\tDis\tDisBetweenPair\tScorePair\tScorePair_PosDisTm\tScore\t",join("\t", @title_info),"\n";
}
foreach my $tid(sort {$a cmp $b}keys %pos){
	my @id1_sort = sort{$pos{$tid}{$a} <=> $pos{$tid}{$b}} keys %{$pos{$tid}};
	my @pos1_sort = sort{$a<=>$b} values %{$pos{$tid}};
	my @id2_sort = @id1_sort;
	my @pos2_sort = @pos1_sort;
	my ($dstart, $dend)=(1,1);
	my $dis_UD = 0; ## distance between U and D
	if(exists $XP{$tid}){ ## template by get_template.pl
		($dstart, $dend)=@{$XP{$tid}}[4..5];
		my ($tid_rmUD) = $tid=~/(\S+)-[UD]/;
		my ($chr1, $s1, $e1, $strand1) = @{$XP{$tid_rmUD."-U"}};
		my ($chr2, $s2, $e2, $strand2) = @{$XP{$tid_rmUD."-D"}};
		$dis_UD = $s2-$e1 -1;
	}
	
	if($type eq "face-to-face:SNP"){ ## target
		my $tid2 = $tid;
		if($tid=~/\-U$/){$tid2=~s/\-U/\-D/;}elsif($tid=~/\-D$/){$tid2=~s/\-D/\-U/;}
		if(!exists $pos{$tid2}){
			print "Wrong: No corresponding temlate ID $tid2 for face_to_face primer, please Check!\n";
			next;
		}
		@id2_sort = sort{$pos{$tid2}{$a} <=> $pos{$tid2}{$b}} keys %{$pos{$tid2}};
		@pos2_sort = sort{$a<=>$b} values %{$pos{$tid2}};
	}elsif(defined $frev){
		@id2_sort = sort{$pos_rev{$tid}{$a} <=> $pos_rev{$tid}{$b}} keys %{$pos_rev{$tid}};
		@pos2_sort = sort{$a<=>$b} values %{$pos_rev{$tid}};
	}
	
	my %score_pair;
	my %score_pair_info;
	my %pos_pair;
	my %dis_pair;
	my $last_pos;
	for(my $i=0; $i<@id1_sort; $i++){
		my $id1 = $id1_sort[$i];
		my $pos = $pos1_sort[$i];

		## get score of pos according to pos of id1
		my $score_pos =0;
		if(defined $range_pos){ # only valid when design for SNP
			next if($pos<$range_pos[2] || $pos>$range_pos[3]);
			my $pos_end = $pos+$dend-$dstart; ## pos to the last target, maybe more targets, -->  |x x|
			next if($pos_end > $range_pos[3]); ## pos out of the pos range
			my $score_pos1 = &score_single($pos, $spos, @range_pos);
			my $score_pos2 = &score_single($pos_end, $spos, @range_pos);
			$score_pos = $score_pos1<$score_pos2? $score_pos1: $score_pos2;
		}
		
		## get score of distance between id1 and id2

# note: P1/P2 is distance from primer1/primer2 3end to SNP or template 3end.
#		L1/L2 is length of primer1/primer2, L is length of template sequence.
#
#            Nested: P2 --->|   x          dis_range(P2-P1):10,15,1,30
#                      P1 --->| x
#
#      face-to-face: |---> P1 x            dis_range(P1+L1+D+P2+L2):180,220,150,250
#         (SNP)               x  P2 <---|  (D is distance between U and D)
#
#      face-to-face: P1 |--->        x     dis_range(P1+L1+P2+L2-L):100,150,70,200
#        (Region)     x        <---| P2    (L is length of template sequence)
#                     ________________
#
#      back-to-back:  x <---| P2           dis_range(P1+L1+P2+L2-L):5,10,0,15
#                      P1 |---> x          (Overlap between p1 and p2: dis > 0)
#                     __________
#
#      back-to-back:  x <---| P2           dis_range(P1+L1+P2+L2-L):-50,-40,-60,-30
#                           P1 |---> x     (No overlap between p1 and p2: dis < 0)
#                     _______________
#
		my ($min_pos, $max_pos); ## min_pos and max_pos of id2
		if($type eq "Nested"){ 
			$min_pos = $pos+$range_dis[2];
			$max_pos = $pos+$range_dis[3];
		}elsif($type eq "face-to-face:SNP"){
			$min_pos = $range_dis[2]-$pos-2*$plen{$id1}-$dis_UD; ## P2=dis-P1-L1-L2-D; use L1 instead L2 because of L2 unkown
			$max_pos = $range_dis[3]-$pos-2*$plen{$id1}-$dis_UD; 
		}elsif($type eq "face-to-face:Region" || $type eq "back-to-back"){
			$min_pos = $range_dis[2]+ $tlength{$tid} - $pos - 2*$plen{$id1};
			$max_pos = $range_dis[3]+ $tlength{$tid} - $pos - 2*$plen{$id1};
		}
		my $indexs = &binarySearch($min_pos, \@pos2_sort, ">=", 0, $#pos2_sort);
		my $indexe = &binarySearch($max_pos, \@pos2_sort, "<=", 0, $#pos2_sort);
		next if($indexs==-1 || $indexe==-1);

		my $score1 = $score{$id1};
		for(my $j=$indexs; $j<=$indexe; $j++){
			my $id2 = $id2_sort[$j];
			my $dis;
			if($type eq "Nested"){
				$dis = $pos2_sort[$j]-$pos;
			}elsif($type eq "face-to-face:SNP"){
				$dis = $pos2_sort[$j]+$plen{$id2}+$pos+$plen{$id1}+$dis_UD;
			}elsif($type eq "face-to-face:Region" || $type eq "back-to-back"){
				$dis = $pos2_sort[$j]+$plen{$id2}+$pos+$plen{$id1}-$tlength{$tid};
			}
			my $score_dis = &score_single($dis, $sdis, @range_dis);
			my $dtm = abs($tm{$id1}-$tm{$id2});
			my $score_tm = &score_single($dtm, $sdtm, @range_dtm);
			my $score = $score1+$score_pos+$score_tm+$score{$id2}+$score_dis;
			$score_pair{$id1.",".$id2}=$score;
			@{$score_pair_info{$id1.",".$id2}}=($score_pos, $score_dis, $score_tm);
			$pos_pair{$id1.",".$id2}=$pos;
			$dis_pair{$id1.",".$id2}=$dis;
		}
	}
	next if(scalar keys %score_pair==0);
	my @final;
	if($ctype eq "Single"){
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
		my $size = $dis_pair{$pair};
		my ($id1_new, $id2_new);
		if($ctype eq "Single"){
			my $le = chr(65+$n); ## A is 65
			$id1_new = $tid."-".$le."-P1";
			$id2_new = $tid."-".$le."-P2";
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

		my ($chrm1, $sm1, $em1, $strandm1)= exists $XP{$chrp1}? @{$XP{$chrp1}}: ($chrp1, 1, $tlength{$chrp1}, "+");
		my $strandp1 = "+"; #P1 (F and nearer to target)
		($pos1, $strand1)=&get_chr_info($sm1, $em1, $strandm1, $strandp1, $plen1, $epos1);
		$chr = $chrm1;
		
		my ($chrm2, $sm2, $em2, $strandm2);
		my $strandp2 = defined $frev? "-": "+"; ## when rev from primer, is "-"; not used in primerScore pipeline 
		if($chrp2=~/_rev$/){
			$chrp2=~s/_rev//;
			($chrm2, $sm2, $em2, $strandm2) = exists $XP{$chrp2}? @{$XP{$chrp2}}: ($chrp2, 1, $tlength{$chrp2}, "+");
			$strandm2 = $strandm2 eq "+"? "-":"+";
			$strandp2 = "+"; ## when rev from template, is "+"
		}else{
			($chrm2, $sm2, $em2, $strandm2)= exists $XP{$chrp2}? @{$XP{$chrp2}}: ($chrp2, 1, $tlength{$chrp2}, "+");
		}
		($pos2, $strand2)=&get_chr_info($sm2, $em2, $strandm2, $strandp2, $plen2, $epos2);

		print O join("\t", $chr, $pos1, $strand1, $id1_new, $seq1, $plen1, $epos1, $size, $score_pair{$pair}, join(",",@{$score_pair_info{$pair}}), @info1[4..$#info1]),"\n";
		print O join("\t", $chr, $pos2, $strand2, $id2_new, $seq2, $plen2, $epos2, $size, $score_pair{$pair}, join(",",@{$score_pair_info{$pair}}), @info2[4..$#info2]),"\n";
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
		if($strand eq "+"){ ## generic, SNP-U, template and primer are both "+"   ------->(-->) 
			$pos = $em-$epos-$plen+1;
			   ##                     <-|  |
		}else{ ## rev from primer, -------->, not used in primerScore pipeline
			$pos = $em-$epos;
		}
	}else{
		$strand = $strandp eq "+"? "-": "+";
		if($strand eq "-"){ ## generic, 1)SNP-D, 2)"face-to-face:Region" and "back-to-back"; template and primer are both "-". <------(<--)
			$pos = $sm+$epos+$plen-1;
			   ##                  |  |->
		}else{ ## rev from primer, <--------, not used in primerScore pipeline
			$pos = $sm+$epos;
		}
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

-rd: distance range of pair primers, (best_min, best_max, min, max) separted by ",", example:
   		
 note: P1/P2 is distance from primer1/primer2 3end to SNP or template 3end.
		L1/L2 is length of primer1/primer2, L is length of template sequence.

      face-to-face: |---> P1 x            dis_range(P1+L1+Lt+P2+L2):180,220,150,250
         (SNP)               x  P2 <---|  (Lt is scale length of target spots)

      face-to-face: P1 |--->        x     dis_range(P1+L1+P2+L2-L):100,150,70,200
        (Region)     x        <---| P2    (L is length of template sequence)
                     ________________

      back-to-back:  x <---| P2           dis_range(P1+L1+P2+L2-L):5,10,0,15
                      P1 |---> x          (Overlap between p1 and p2: dis > 0)
                     __________

      back-to-back:  x <---| P2           dis_range(P1+L1+P2+L2-L):-50,-40,-60,-30
                           P1 |---> x     (No overlap between p1 and p2: dis < 0)
                     _______________

            Nested: P2 --->|   x          dis_range(P2-P1):10,15,1,30
                      P1 --->| x


Usage:
  Options:
  -i   <file>   Input primer score file, forced
  -it  <file>   Input template file, forced
  -ir  <file>   Input rev primer score file to design face-to-face(Region) and back-to-back primer, optional
  -k   <str>	Key of output file, forced
  -tp  <str>    primer type, "face-to-face:SNP", "face-to-face:Region", "back-to-back", "Nested", ["face-to-face:SNP"]

  -rd  <str>	distance range of pair primers, (dis_score, best_min, best_max, min, max) separted by ",", ["10,10,15,1,30"]
  -rp  <str>	position range when score for distance of p1 to the detected site, can't choose when -ct is "Full-covered", optional
  -ct     <str>	    primer covered type, "Single" or "Full-covered", ["Single"]
  	 -ds  <int>		distance when -ct "Full-covered", [500]
	 -rf  <float>	distance float ratio when -ct "Full-covered", [0.2]
	 -on  <int>     output num when -ct is "Single",[$onum]

  -od <dir>	Dir of output file, default ./
  -h		 Help

USAGE
	print $usage;
	exit;
}


