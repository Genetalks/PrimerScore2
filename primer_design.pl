#!/usr/bin/perl -
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
require "$Bin/path.pm";
require "$Bin/self_lib.pm";

my $BEGIN_TIME=time();
my $version="1.0.0";
#######################################################################################

# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
my ($NoSpecificity, $FiterRepeat);
my $min_len = 20;
my $max_len = 35;
my $range_dis;
my $min_tm = 65;
my $max_tm = 75;
my $min_gc = 0;
my $max_gc = 0.8;
my ($fIn,$fkey,$outdir);
my $step = 1;
my $recov;
my $para_num = 10;
my $stm = 45;
my $scale_len = 2;
my $scale_dis = 2;
my $line_num = 3;
our $PATH_PRIMER3;
our $REF_HG19;
my $fref = $REF_HG19;
my $format_dis = 3;
GetOptions(
				"help|?" =>\&USAGE,
				"i:s"=>\$fIn,
				"r:s"=>\$fref,
				"k:s"=>\$fkey,
				"NoSpecificity:s"=>\$NoSpecificity,
				"FilterRepeat:s"=>\$FiterRepeat,
				"recov:s"=>\$recov,
				"mintm:s"=>\$min_tm,
				"maxtm:s"=>\$max_tm,
				"mingc:s"=>\$min_gc,
				"maxgc:s"=>\$max_gc,
				"minl:s"=>\$min_len,
				"maxl:s"=>\$max_len,
				"scalel:s"=>\$scale_len,
				"fdis:s"=>\$format_dis,
				"rdis:s"=>\$range_dis,
				"stm:s"=>\$stm,
				"lnum:s"=>\$line_num,
				"para:s"=>\$para_num,
				"step:s"=>\$step,
				"od:s"=>\$outdir,
				) or &USAGE;
&USAGE unless ($fIn and $fkey);

$outdir||="./";
`mkdir $outdir`	unless (-d $outdir);
$outdir=AbsolutePath("dir",$outdir);
my $oligotm = "$PATH_PRIMER3/src/oligotm";

## score value: min_best, max_best, min, max
my @endA = (0, 0, 0, 2);
my @len = ($min_len,$min_len+5,$min_len, $max_len);
my @tm = ($min_tm+int(($max_tm-$min_tm)/2), $max_tm-int(($max_tm-$min_tm)/2), $min_tm, $max_tm);
my @gc = (0.5, 0.65, 0.3, 0.8);
my @dgc = (-0.50, -0.25, -1, 0.75);
my @mdgc = (0, 0.5, 0, 0.75);
my @tmh = (-50, 40, -50, 55);
my @nalign = (0,40,0,55); #sec tm
my @poly = (0,5,0,20);

my $GC5_num = 8; #stat GC of primer start $GC5_num bp
my $GC3_num = 8;
my $min_tm_diff = 5;
my $plen = int(($max_len+$min_len)/2);


my %seq;

if ($step == 1){
	my @rdis;
	if(defined $range_dis){
		@rdis = &check_merge_rdis($range_dis);
	}
	open(P,">$outdir/$fkey.primer.list") or die $!;
	open(I, $fIn) or die $!;
	$/=">";
	while(<I>){
		chomp;
		next if(/^$/);
		my ($id_info, @line)=split /\n/, $_;
		my $seq = join("", @line);
		my ($id)=split /\s+/, $id_info;
		my ($dstart)=$id_info=~/XS:i:(\d+)/;
		my ($dend)=$id_info=~/XE:i:(\d+)/;
		if(!defined $dstart){
			$dstart = 0;
		}
		if(!defined $dend){
			$dend = 0;
		}
		
		## 
		my $tlen = length($seq);
		if(!defined $range_dis){
			@rdis = (1, $tlen, int(($tlen/50)+1));
		}
		if($format_dis == 5){ ## convert to format '3'
			for(my $i=0; $i<@rdis; $i+=3){
				my $e =$tlen-$rdis[$i]-$plen;
				my $s =$tlen-$rdis[$i+1]-$plen;
				$rdis[$i]=$s;
				$rdis[$i+1]=$e;
			}
		}
		for(my $l=$min_len; $l<=$max_len; $l+=$scale_len){
			my $rs=1;
			for(my $r=0; $r<@rdis; $r+=3){
				my ($min_dis, $max_dis) = ($rdis[$r], $rdis[$r+1]);
				my $min_p = $min_dis-$dstart>0? $min_dis-$dstart: 0; ## dstart=1
				my $max_p = $max_dis-$dend<$tlen-$l? $max_dis-$dend: $tlen-$l;
				if($max_p < $min_p){
					print "Wrong: max position < min position! maybe dend (XE:i:$dend) is too large, or -rdis range $range_dis is too narrow!\n";
					die;
				}
				my $sdis = $rdis[$r+2];
				for(my $p=$min_p; $p<=$max_p; $p+=$sdis){
					my $primer=&get_primer($id, $p, $l, $seq);
					my $id_new = $id."-".$l."-".$p;
					if(defined $FiterRepeat){
						my @match = ($primer=~/[atcg]/g);
						next if(scalar @match > length($primer)*0.4);
					}
					print P $id_new,"\t",$primer,"\n";
					$seq{$id_new}=$primer;
				}
			}
		}
	}
	close(P);
	$step ++;	
}else{
	open(I,"$outdir/$fkey.primer.list") or die $!;
	while(<I>){
		chomp;
		next if(/^$/);
		my ($id, $seq)=split;
		$seq{$id}=$seq;
	}
	close(I);
}

# evalue
if($step ==2){
	my $odir = "$outdir/split";
	`rm -r $odir` if(-d $odir);
	mkdir $odir;
	my ($rank, @fprimer) = &split_file("$outdir/$fkey.primer.list", $odir, $line_num, 100);

	open (SH, ">$outdir/$fkey.primer.evalue.sh") or die $!;
	foreach my $f (@fprimer){	
		my $fname = basename($f);
		my $dir_new = dirname($f);
		if(defined $NoSpecificity){
			print SH "perl $Bin/primer_evaluation.pl --nohead --NoSpecificity -p $f -d $fref -n 1 -thread $line_num -stm $stm -mintm $min_tm -maxtm $max_tm -mingc $min_gc -maxgc $max_gc -k $fname -od $dir_new >$dir_new/$fname.log 2>&1\n";
		}else{
			print SH "perl $Bin/primer_evaluation.pl --nohead -p $f -d $fref -n 1 -thread $line_num -stm $stm -mintm $min_tm -maxtm $max_tm -mingc $min_gc -maxgc $max_gc -k $fname -od $dir_new >$dir_new/$fname.log 2>&1\n";
		}
	}
	close (SH);
	my $timeout = 200*$line_num;
	#my $timeout = 1*$line_num;
	Run("parallel -j $para_num --timeout $timeout < $outdir/$fkey.primer.evalue.sh", 1);
	
	##cat
	if($rank==2){
		my @dirs = glob("$odir/dir_*");
		foreach my $dir (@dirs){
			Run("cat $dir/*.evaluation.out > $dir/evaluation.out", 1);
			Run("cat $dir/*.filtered_by_Tm_GC.list > $dir/filtered_by_Tm_GC.list", 1);
		}
		Run("cat $odir/*/evaluation.out >$outdir/$fkey.primer.evaluation.out", 1);
		Run("cat $odir/*/filtered_by_Tm_GC.list >$outdir/$fkey.primer.filtered_by_Tm_GC.list", 1);
	}else{
		Run("cat $odir/*evaluation.out > $outdir/$fkey.primer.evaluation.out", 1);
		Run("cat $odir/*filtered_by_Tm_GC.list > $outdir/$fkey.primer.filtered_by_Tm_GC.list", 1);
	}

	##get bwa fail list
	my %suc;
	my %filter;
	open(E, "$outdir/$fkey.primer.evaluation.out") or die $!;
	$/="\n";
	while(<E>){
		chomp;
		my ($id)=split;
		$suc{$id}=1;
	}
	close(E);
	open(F, "$outdir/$fkey.primer.filtered_by_Tm_GC.list") or die $!;
	while(<F>){
		chomp;
		my ($id)=split;
		$filter{$id}=1;
	}
	close(F);
	my %tid;
	open(FB, ">$outdir/$fkey.primer.filtered_by_timeout.list") or die $!;
	foreach my $id(keys %seq){
		my ($tid, $len, $dis)=$id=~/(\S+)-(\d+)-(\d+)$/;
		$tid{$tid}{"Total"}++;
		if(exists $suc{$id}){
			$tid{$tid}{"Suc"}++;
		}elsif(!exists $filter{$id}){
			print FB join("\t", $id, $seq{$id}),"\n";
		}
	}
	close(FB);
	open(FO, ">$outdir/$fkey.primer_design.summary") or die $!;
	print FO "TemplateID\tTotalCandidate\tEvaluatedSuccess\tPrimerDesignedSuccessOrNot\n";
	foreach my $tid(keys %tid){
		if(!exists $tid{$tid}{"Suc"}){
			$tid{$tid}{"Suc"} = 0;
		}
		my $type = "Success";
		if($tid{$tid}{"Suc"}<=2){
			$type = "Failure";
		}elsif($tid{$tid}{"Suc"}<10){
			$type = "Warning";
		}
		print FO join("\t", $tid, $tid{$tid}{"Total"}, $tid{$tid}{"Suc"}, $type), "\n";
	}
	close(FO);

	$step++;
}

my %info;
my %score;
# score and output
if($step==3){
	open(S,">$outdir/$fkey.primer.score") or die $!;
	print S "##Score_info: sendA, spoly, slen, stm, sgc, sdgc, smdgc, shairpin, snalign\n";
	print S "##High_Info(n=1): TM        : Flag/DatabaseID/Pos/Cigar/MD/End_match_num/Efficiency\n";
	print S "##High_Info(n>1): Efficiency: Flag/DatabaseID/Pos/Cigar/MD/End_match_num/TM\n";
#	open(O,">$outdir/$fkey.single.final.primer") or die $!;
#	print O "##Score_info: sendA, spoly, slen, stm, sgc, sdgc, smdgc, shairpin, snalign\n";
#	print O "##High_Info(n=1): TM        : Flag/DatabaseID/Pos/Cigar/MD/End_match_num/Efficiency\n";
#	print O "##High_Info(n>1): Efficiency: Flag/DatabaseID/Pos/Cigar/MD/End_match_num/TM\n";
	open (I, "$outdir/$fkey.primer.evaluation.out") or die $!;
	$/="\n";
	my @title_info=split /\t/, "Tm\tGC\tGC5\tGC3\tdGC\tdG_Hairpin\tTm_Hairpin\tdG_Dimer\tTm_Dimer\tAlign_Num\tHigh_Tm_Num\tHigh_Efficiency_Num\tHigh_Info";
	while(<I>){
		chomp;
		next if(/^$/);
		my ($id,$seq, $len, @info)=split/\t/, $_;
		my ($tm, $gc, $gc5, $gc3, $dgc, $dg_h, $tm_h, $dg_d, $tm_d, $nalign, $nhtm, $neff, $htm_info)=@info;
		#$nend=~s/\+//;
		my ($score, $score_info) = &score($seq{$id}, $len, $tm, $gc, $gc5, $gc3, $dgc, $dg_h, $tm_h, $nhtm, $htm_info);	
		my ($id_sub, $dis)=$id=~/(\S+)\-\d+\-(\d+)$/;
		if(defined $recov){
			$id_sub=~s/rev//;
		}
		push @{$score{$id_sub}{$score}},$id;
		@{$info{$id}}=($id_sub, $len,$dis, $seq{$id}, $score, $score_info, @info);
	}
	close(I);
	$step++;
	
	print S "#ID\tTarget\tLen\tDis\tSeq\tScore\tScore_info\t",join("\t", @title_info),"\n";
#	print O "#ID\tTarget\tLen\tDis\tSeq\tScore\tScore_info\t",join("\t", @title_info),"\n";
	foreach my $id_sub(sort {$a cmp $b} keys %score){
		my $n=0;
		my $flag = 0;
		foreach my $s(sort {$b<=>$a} keys %{$score{$id_sub}}){
			my @id = @{$score{$id_sub}{$s}};
			foreach my $id(@id){
				print S $id,"\t",join("\t",@{$info{$id}}),"\n";
				if($flag == 0){
					#print O $id,"\t",join("\t",@{$info{$id}}),"\n";
					$flag = 1;
				}
				$n++;
			}
#			if($n>$max_num){
#				last;
#			}
		}
	}
#	close(O);
	close(S);
}
#######################################################################################
print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################

# ------------------------------------------------------------------
# sub function
# ------------------------------------------------------------------

#dfnum: files num in one dir
#flnum: lines num in one file
#rnum: rank num of dir
sub split_file{
	my ($file, $odir, $flnum, $dfnum)=@_;
	my $total = `wc -l $file`;
	($total) = split /\s+/, $total;
	my $fname = basename($file);
	my $dlnum = $dfnum*$flnum;
	my @sfile;
	my $rnum;
	if($total/$dlnum < 1){ ## one rank
		Run("split -l $flnum $file $odir/$fname\_");
		@sfile = glob("$odir/$fname\_*");
		$rnum = 1;
	}else{ ## two rank
		Run("split -l $dlnum $file $odir/$fname\_");
		@sfile = glob("$odir/$fname\_*");
		foreach my $f (@sfile){
			my ($nid) = $f=~/$fname\_(\S+)/;
			my $dir_new = "$odir/dir_$nid";
			mkdir $dir_new unless(-d $dir_new);
			Run("split -l $flnum $f $dir_new/$fname\_$nid\_");
		}
		@sfile = glob("$odir/dir_*/$fname\_*");
		$rnum = 2;
	}
	return ($rnum, @sfile);
}


sub score_single{
	my ($v, $score, $minb, $maxb, $min, $max)=@_;
	return $score if($v eq "NULL");
	my $disl = $minb - $min;
	my $disr = $max - $maxb;
	if($disl<0 || $disr<0){
		print "Wrong Score set: $minb, $maxb, $min, $max\n";
		die;
	}

	my $s;
	if($v<$minb){
		$disl=$disl==0? 0.1: $disl;
		$s = int($score * (1 - ($minb-$v)/$disl));
	}elsif($v<=$maxb){
		$s = $score;
	}else{
		$disr= $disr==0? 0.1: $disr;
		$s = int($score * (1 - ($v-$maxb)/$disr));
	}
	return $s;
}
sub score{ #&score($seq{$id},$dis, $len, $tm, $gc, $gc5, $gc3, $dgc, $dg_h, $tm_h, $nhtm, $htm);
	my ($seq, $len, $tm, $gc, $gc5, $gc3, $mdgc, $dg_h, $tm_h, $nalign, $htm_info)=@_;
	my $s = 0;
	#my ($sendA, $spoly, $slen, $stm, $sgc, $sdgc, $smdgc, $snalign)  = (8, 15, 7, 22, 2, 6, 12,  28);
	my ($sendA, $spoly, $slen, $stm, $sgc, $sdgc, $smdgc, $shpin, $snalign)  = (6, 14, 4, 20, 2, 4, 8, 20, 22);
	my $nendA = &get_end_A($seq);
	my $vpoly = &get_poly_value($seq);

	
	my $dgc = $gc3-$gc5;
	my @score;
	push @score, &score_single($nendA, $sendA, @endA);
	push @score, &score_single($vpoly, $spoly, @poly);
	push @score, &score_single($len, $slen, @len);
	push @score, &score_single($tm, $stm, @tm);
	push @score, &score_single($gc, $sgc, @gc);
	push @score, &score_single($dgc, $sdgc, @dgc);
	push @score, &score_single($mdgc, $smdgc, @mdgc);
	if($tm_h ne ""){
		push @score, &score_single($tm_h, $shpin, @tmh);
	}else{
		push @score, $shpin;
	}
	
	my $s_nalign;
	if(defined $NoSpecificity){
		$s_nalign = $snalign;
	}else{
		$nalign=~s/\+//;
		my ($htm)=split /:/, $htm_info;
		my @htm = split /,/, $htm;
		if(@htm==1){
			$s_nalign = $snalign*1;
		}else{
			$s_nalign = &score_single($htm[1], $snalign, @nalign);
		}
	}
	push @score, $s_nalign;
	my $ssum = 0;
	$ssum += $_ foreach @score; 
	return ($ssum, join(",",@score));
}

sub get_poly_value{
	my ($seq)=@_;
	$seq = reverse $seq;
	my @unit = split //, $seq;
	
	my @polys;
	my $last = $unit[0];
	my $poly = $unit[0];
	my $dis = 0; #dis to the 3 end
	for(my $i=1; $i<@unit; $i++){
		if($unit[$i] eq $last){
			$poly.=$unit[$i];
		}else{
			if(length($poly)>=3){
				push @polys, [$poly, $dis];
			}
			$dis = $i;
			$last = $unit[$i];
			$poly=$last;
		}
	}
	if(length($poly)>=3){
		push @polys, [$poly, $dis];
	}
	
	my $value1 = 0;
	my $pslen = 0;
	my $psnum = 0;
	for(my $i=0; $i<@polys; $i++){
		my $l = length($polys[$i][0]);
		my $d = $polys[$i][1];
		$pslen+=$l;
		$psnum++;
		if($i==0){
			$value1+=$l*(6-$d); # end poly impact
		}
	}
	$value1 = $value1>=0? $value1: 0;

	## other impact
	my $value2 = 0.2*$pslen*(6-$psnum);
	$value2 = $value2>=0? $value2: 0;
	my $value = $value1 + $value2;
	return $value;
}

sub get_end_A{
	my ($seq)=@_;
	my @unit = split //, $seq;
	my $Anum = 0;
	my $rnum = 0;
	my %class=(A=>"AT",T=>"AT",C=>"CG",G=>"CG",a=>"AT",t=>"AT",c=>"CG",g=>"CG");

	my $endc = $class{$unit[-1]};;
	for(my $i=@unit; $i>0; $i--){
		last if($class{$unit[$i-1]} ne $endc);
		$rnum++;
	}
	$Anum = $endc eq "AT"? $rnum: 0;
	return ($Anum);
}


sub GC_stat{
        my ($p)=@_;
        my @u = split //, $p;
        my %stat;
        $stat{"total"}{'G'}=0;
        $stat{"total"}{'C'}=0;
        $stat{5}{'G'}=0;
        $stat{5}{'C'}=0;
        $stat{3}{'G'}=0;
        $stat{3}{'C'}=0;
        my $total = scalar @u;
        for(my $i=0; $i<$total; $i++){
            if ($i<$GC5_num){
                $stat{5}{$u[$i]}++;
            }
            if ($i>$total-$GC3_num-1){
                $stat{3}{$u[$i]}++;
            }
            $stat{"total"}{$u[$i]}++;
        }
        my $GC = ($stat{"total"}{'G'}+$stat{"total"}{'C'})/$total;
        my $GC5= ($stat{5}{'G'}+$stat{5}{'C'})/$GC5_num;
        my $GC3 = ($stat{3}{'G'}+$stat{3}{'C'})/$GC3_num;
        return ($GC, $GC5, $GC3);
}

sub get_primer{
	my ($id, $pos, $len, $seq)=@_;
	my $total_len = length $seq;
	my $primer = substr($seq, $total_len-$pos-$len, $len);
	if(defined $recov){
		$primer =~tr/ATGCatgc/TACGtacg/;
		$primer = reverse $primer;
		$id.="rev";
#		$pos = $total_len-$pos-$len;
	}
	return ($primer);
}

sub check_merge_rdis{
	my ($rdis)=$_;
	##check range_dis
	my @rdis = split /,/, $range_dis;
	my $nrdis = scalar @rdis;
	if($nrdis!=3 && $nrdis!=6){
		print "Wrong: number of -rdis must be 3 or 6!\n";
		die;
	}
	for(my $i=0; $i<@rdis; $i+=3){
		if($rdis[$i+1]-$rdis[$i] < 0){
			print "Wrong: -rdis region must be ascending ordered, eg:3,40,1,100,150,5\n";
			die;
		}
	}

	if($nrdis==3){
		return(@rdis);
	}else{ ##merge when overlap
		my ($s1, $e1, $b1, $s2, $e2, $b2) = @rdis;
		## sort two regions
		if($s1 > $s2){
			($s2, $e2, $b2) = @rdis[0..2];
			($s1, $e1, $b1) = @rdis[3..5];
		}
		if($e1>$s2){ ## overlap
			if($e1>$e2){ ## r1 include r2
				if($b1<$b2){ ## prefer min bin
					return($s1, $e1, $b1);
				}else{
					return ($s1, $s2, $b1, $s2, $e2, $b2, $e2, $e1, $b1);
				}
			}else{## intersect 
				if($b1<$b2){
					return ($s1, $e1, $b1, $e1, $e2, $b2);
				}else{
					return ($s1, $s2, $b1, $s2, $e2, $b2);
				}
			}
		}else{## no overlap
			return @rdis;
		}

	}
}


sub Run{
    my ($cmd, $nodie)=@_;

    my $ret = system($cmd);
    if (!$nodie && $ret) {
        die "Run $cmd failed!\n";
    }
}


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
  -i  	<file>   	Input template fa file, forced
  -r  	<file>   	Input ref fa file, [$fref]
  -k  	<str>		Key of output file, forced

  --recov           recov primer seq
  --NoSpecificity   not evalue specificity
  --FilterRepeat	filter primers with repeat region(lowercase in fref) more than 40%
  -minl	<int> 		min primer len, [$min_len]
  -maxl	<int>  		max primer len, [$max_len]
  -mintm <int>		min tm, [$min_tm]
  -maxtm <int>		max tm, [$max_tm]
  -mingc <float>	min gc, [$min_gc]
  -maxgc <float>	max gc, [$max_gc]
  -scalel <int>     scale len, [$scale_len]
  -fdis <5,3>       distance caculation format, 3: from primer right end to template 3'end; 5: from left right end to template 5'end, [$format_dis]
  -rdis <str>		region ranges(start,end,scale), start <= end, count format see -fdis, separated by ",", optional
  		            Example: 
			            sanger sequence primer: 100,150,2,400,500,5
			            ARMS PCR primer: 1,1,1,80,180,2

  -lnum  <int>      line num in one separated file, [$line_num]
  -stm  <int>		min tm to be High_tm in specifity, [$stm]
  -para  <int>		parallel num, [$para_num]
  -step	<int>		step, [$step]
  	1: get primer seq
	2: evalue primer
	3: score and output
  -od <dir>	Dir of output file, default ./
  -h		 Help

USAGE
	print $usage;
	exit;
}

