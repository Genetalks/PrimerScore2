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
my ($NoSpecificity, $FiterRepeat, $NoFilter, $NoCoverN);
my ($ftem,$ftem_snp,$fkey,$outdir);
my $step = 1;
my $para_num = 10;
my $stm = 40;
my $opt_tm=60;
my $opt_tm_probe=70;
our $PATH_PRIMER3;
our $REF_HG19;
my $fref = $REF_HG19;
my $format_dis = 3;
my $probe;
my $ptype = "face-to-face";
my $range_dis;
my $range_len="18,36,2";
GetOptions(
				"help|?" =>\&USAGE,
				"i:s"=>\$ftem,
				"is:s"=>\$ftem_snp,
				"r:s"=>\$fref,
				"k:s"=>\$fkey,
				"Probe:s"=>\$probe,
				"NoSpecificity:s"=>\$NoSpecificity,
				"FilterRepeat:s"=>\$FiterRepeat,
				"NoFilter:s"=>\$NoFilter,
				"NoCoverN:s"=>\$NoCoverN,
				"ptype:s"=>\$ptype,
				"rlen:s"=>\$range_len,
				"opttm:s"=>\$opt_tm,
				"opttmp:s"=>\$opt_tm_probe,
				"fdis:s"=>\$format_dis,
				"rdis:s"=>\$range_dis,
				"stm:s"=>\$stm,
				"para:s"=>\$para_num,
				"step:s"=>\$step,
				"od:s"=>\$outdir,
				) or &USAGE;
&USAGE unless ($ftem and $fkey);

$outdir||="./";
`mkdir $outdir`	unless (-d $outdir);
$outdir=AbsolutePath("dir",$outdir);
my $oligotm = "$PATH_PRIMER3/src/oligotm";

## score value: min_best, max_best, min, max
my ($min_len, $max_len, $scale_len)=split /,/, $range_len;
my @endA = (0, 0, 0, 2);
my @len = ($min_len,$min_len+5,$min_len, $min_len+10);
my @tm = ($opt_tm-1, $opt_tm+1, $opt_tm-5, $opt_tm+5);
my @gc = (0.45, 0.65, 0.3, 0.8);
my @dgc = (-0.50, -0.25, -1, 0.75);
my @mdgc = (0, 0.5, 0, 0.75);
my @tmh = (-50, 40, -50, 55); ## hairpin tm
my @nalign = (0,40,0,55); #sec tm
my @poly = (0,2,0,10);

my $GC5_num = 8; #stat GC of primer start $GC5_num bp
my $GC3_num = 8;
my $min_tm_diff = 5;
my %seq;

if ($step == 1){
	my $n=0;
	my @rdis;
	if(defined $range_dis){
		@rdis = &check_merge_rdis($range_dis);
	}

	## get template seq containing snp
	my %seq_snp;
	if(defined $ftem_snp){
		open(I, $ftem_snp) or die $!;
		$/=">";
		while(<I>){
			chomp;
			next if(/^$/);
			my ($id_info, @line)=split /\n/, $_;
			my $seq = join("", @line);
			my ($id)=split /\s+/, $id_info;
			$seq_snp{$id}=$seq;
		}
	}
	
	## get primer seq
	open(I, $ftem) or die $!;
	open(PT, ">$outdir/$fkey.primer.list") or die $!;
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
		my $seq_snp = $seq_snp{$id};
		
		## 
		my $tlen = length($seq);
		if(!defined $range_dis){
			@rdis = (1, $tlen, 1);
		}
		my @rdisc=@rdis;
		if($format_dis == 5){ ## convert to format '3'
			for(my $i=0; $i<@rdis; $i+=3){
				my $e =$tlen-$rdis[$i]-$min_len;
				my $s =$tlen-$rdis[$i+1]-$max_len;
				$rdisc[$i]=$s;
				$rdisc[$i+1]=$e;
			}
		}
		for(my $r=0; $r<@rdisc; $r+=3){
			my ($min_dis, $max_dis) = ($rdisc[$r], $rdisc[$r+1]);
			my $min_p = $min_dis-$dstart>0? $min_dis-$dstart: 0; ## dstart=1
			my $max_p = $max_dis-$dend<$tlen? $max_dis-$dend: $tlen;
			if($max_p < $min_p){
				print "Wrong: max position $max_p < min position $min_p! ($min_dis, $max_dis, $dstart, $dend) Maybe dend (XE:i:$dend) is too large, or -rdis range $range_dis is too narrow!\n";
				die;
			}
			my $sdis = $rdis[$r+2];
			for(my $p=$min_p; $p<=$max_p; $p+=$sdis){
				$n++;
				my $dn=int($n/1000);
				my $dir = "$outdir/split_$dn";
				`mkdir $dir` unless(-d $dir);
				my $fn=$n%1000;
				open(P, ">$dir/$fkey.primer.list_$fn") or die $!;
				for(my $l=$max_len; $l>=$min_len; $l-=$scale_len){
					next if($p+$l>$tlen);
					my ($primer)=&get_primer($id, $p, $l, $seq);
					next if(defined $NoCoverN && ($primer=~/N/ || $primer=~/n/));
					my $id_new = $id."-".$p; ##primers of different length are evalued in primer_evaluation.pl
					#my $id_new = $id."-".$p."-".$l;
					if(defined $FiterRepeat){
						my @match = ($primer=~/[atcg]/g);
						next if(scalar @match > length($primer)*0.4);
					}
					
					if(defined $ftem_snp){
						my ($primer_snp)=&get_primer($id, $p, $l, $seq_snp{$id});
						print P $id_new,"\t",$primer,":", $primer_snp, "\n";
						print PT $id_new,"\t",$primer,":", $primer_snp, "\n";
					}else{
						print P $id_new,"\t",$primer,"\n";
						print PT $id_new,"\t",$primer,"\n";
					}
					$seq{$id_new}=$primer;
					last; ##  primers of different length are evalued in primer_evaluation.pl
				}
				close(P);
			}
		}
	}
	close(PT);
	$step ++;
}

# evalue
if($step ==2){
	my @fprimer = glob("$outdir/split_*/$fkey.primer.list*");
	open (SH, ">$outdir/$fkey.primer.evalue.sh") or die $!;
	foreach my $f (@fprimer){	
		my $fname = basename($f);
		my $dir_new = dirname($f);
		my $olens=join(",", $min_len,$max_len, $scale_len); 
		my $cmd = "perl $Bin/primer_evaluation.pl --nohead -p $f -d $fref -thread 1 -stm $stm -k $fname -opttm $opt_tm -olen $olens -od $dir_new";
		if($ptype eq "face-to-face" || $ptype eq "back-to-back"){
			$cmd .= " --Revcom";
		}
		if(defined $probe){
			$cmd .= " --Probe -opttmp $opt_tm_probe";
		}
		if(defined $NoSpecificity){
			$cmd .= " --NoSpecificity";
		}
		if(defined $NoFilter){
			$cmd .= " --NoFilter";
		}
		$cmd .= " >$dir_new/$fname.log 2>&1";
		print SH $cmd,"\n";
	}
	close (SH);
	Run("parallel -j $para_num  < $outdir/$fkey.primer.evalue.sh", 1);
	
	##cat
	my @dirs = glob("$outdir/split_*");
	foreach my $dir (@dirs){
		Run("cat $dir/*.evaluation.out > $dir/evaluation.out", 1);
		Run("cat $dir/*.bound.info > $dir/bound.info", 1);
		Run("cat $dir/*.filter.list > $dir/filter.list", 1);
	}
	Run("cat $outdir/*/evaluation.out >$outdir/$fkey.primer.evaluation.out", 1);
	Run("cat $outdir/*/bound.info >$outdir/$fkey.primer.bound.info", 1);
	Run("cat $outdir/*/filter.list >$outdir/$fkey.primer.filter.list", 1);

	$step++;
}
die;
my %info;
my %score;
# score and output
if($step==3){
	open(S,">$outdir/$fkey.primer.score") or die $!;
	print S "##Score_info: sendA, spoly, slen, stm, sgc, sdgc, smdgc, shairpin, snalign\n";
	print S "##High_Info(n=1): TM        : Flag/DatabaseID/Pos/Cigar/MD/End_match_num/Efficiency\n";
	print S "##High_Info(n>1): Efficiency: Flag/DatabaseID/Pos/Cigar/MD/End_match_num/TM\n";
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
		my ($id_sub, $dis)=$id=~/(\S+)\-(\d+)\-\d+$/;
		push @{$score{$id_sub}{$score}},$id;
		@{$info{$id}}=($id_sub, $len,$dis, $seq{$id}, $score, $score_info, @info);
	}
	close(I);
	$step++;
	
	print S "#ID\tTarget\tLen\tDis\tSeq\tScore\tScore_info\t",join("\t", @title_info),"\n";
	foreach my $id_sub(sort {$a cmp $b} keys %score){
		my $n=0;
		my $flag = 0;
		foreach my $s(sort {$b<=>$a} keys %{$score{$id_sub}}){
			my @id = @{$score{$id_sub}{$s}};
			foreach my $id(@id){
				print S $id,"\t",join("\t",@{$info{$id}}),"\n";
				if($flag == 0){
					$flag = 1;
				}
				$n++;
			}
		}
	}
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



sub score{ #&score($seq{$id},$dis, $len, $tm, $gc, $gc5, $gc3, $dgc, $dg_h, $tm_h, $nhtm, $htm);
	my ($seq, $len, $tm, $gc, $gc5, $gc3, $mdgc, $dg_h, $tm_h, $nalign, $htm_info)=@_;
	my $s = 0;
	#my ($sendA, $spoly, $slen, $stm, $sgc, $sdgc, $smdgc, $snalign)  = (8, 15, 7, 22, 2, 6, 12,  28);
#	my ($sendA, $spoly, $slen, $stm, $sgc, $sdgc, $smdgc, $shpin, $snalign)  = (6, 14, 4, 20, 2, 4, 8, 20, 22);
	my ($sendA, $spoly, $slen, $stm, $sgc, $sdgc, $smdgc, $shpin, $snalign)  = (6, 12, 6, 20, 6, 4, 8, 20, 18);
	my $nendA = &get_end_A($seq);
	my $vpoly = &get_poly_value($seq);
	my $dgc = $gc3-$gc5;
	my @score;
	push @score, int(&score_single($nendA, $sendA, @endA));
	push @score, int(&score_single($vpoly, $spoly, @poly));
	push @score, int(&score_single($len, $slen, @len));
	push @score, int(&score_single($tm, $stm, @tm));
	push @score, int(&score_single($gc, $sgc, @gc));
	push @score, int(&score_single($dgc, $sdgc, @dgc));
	push @score, int(&score_single($mdgc, $smdgc, @mdgc));
	if($tm_h ne ""){
		push @score, int(&score_single($tm_h, $shpin, @tmh));
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
			$s_nalign = int(&score_single($htm[1], $snalign, @nalign));
		}
	}
	push @score, $s_nalign;
	my $ssum = 0;
	$ssum += $_ foreach @score; 
	return ($ssum, join(",",@score));
}

sub get_primer{
	my ($id, $pos, $len, $seq)=@_;
	my $total_len = length $seq;
	if($total_len-$pos-$len<0){
		print "Extract primer failed! $seq, $pos, $len, $total_len-$pos-$len\n";
		die;
	}
	my $primer = substr($seq, $total_len-$pos-$len, $len);
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



sub USAGE {#
	my $usage=<<"USAGE";
Program:
Version: $version
Contact:zeng huaping<huaping.zeng\@genetalks.com> 

Usage:
  Options:
  -i  	<file>   	Input template fa file, forced
  -is  	<file>   	Input template add snp fa file, optional
  -r  	<file>   	Input ref fa file, [$fref]
  -k  	<str>		Key of output file, forced

  --Probe                design probe
  --NoSpecificity   not evalue specificity
  --FilterRepeat	filter primers with repeat region(lowercase in fref) more than 40%
  --NoFilter             Not filter any primers
  --NoCoverN             candidate primer sequences can't contain N/n

  -ptype     <str>       primer type, "face-to-face", "back-to-back", "Nested", [$ptype]
  -opttm    <int>       optimal tm, [$opt_tm]
  -opttmp   <int>       optimal tm of probe, [$opt_tm_probe]
  -rlen     <str>       primer len ranges(start,end,scale), start <= end, [$range_len]
  -rdis     <str>       region ranges(start,end,scale), start <= end, count format see -fdis, separated by ",", (1,len,1) when not given, optional
  		                Example: 
			               sanger sequence primer: 100,150,2,400,500,5
			               ARMS PCR primer: 1,1,1,80,180,2
  -fdis    <5,3>        distance caculation format, 3: from primer right end to template 3'end; 5: from left right end to template 5'end, [$format_dis]
  -stm     <int>		min tm to be High_tm in specifity, [$stm]
  -para    <int>		parallel num, [$para_num]
  -step	   <int>		step, [$step]
  	1: get primer seq
	2: evalue primer
	3: score and output
  -od <dir>	Dir of output file, default ./
  -h		 Help

USAGE
	print $usage;
	exit;
}

