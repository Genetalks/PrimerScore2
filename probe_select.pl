#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
require "$Bin/self_lib.pm";

my $BEGIN_TIME=time();
my $version="1.0.0";
#######################################################################################

# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
my ($fkey,$outdir);
my ($fprimer, $fprobes, $ftemplate);
my $opt_tm = 70;
my $min_len = 18;
my $max_len = 45;
my $nprobe=1;
GetOptions(
				"help|?" =>\&USAGE,
				"ip:s"=>\$fprimer,
				"it:s"=>\$ftemplate,
				"io:s"=>\$fprobes,

				"opttm:s"=>\$opt_tm,
				"minl:s"=>\$min_len,
				"maxl:s"=>\$max_len,
				"num:s"=>\$nprobe,
				"k:s"=>\$fkey,
				"od:s"=>\$outdir,
				) or &USAGE;
&USAGE unless ($fprimer and $fprobes and $ftemplate and $fkey);

$outdir||="./";
`mkdir $outdir`	unless (-d $outdir);
$outdir=AbsolutePath("dir",$outdir);

my @tm = ($opt_tm-1, $opt_tm+1, $opt_tm-5, $opt_tm+5);
my $dlen = $max_len-$min_len;
my @len = ($min_len+$dlen/3, $max_len-$dlen/3, $min_len, $max_len);
my @dis = (5,8,1,50);
my @CGd = (0.1, 1, 0, 1);

my %template;
my %primer;
##Chr    Start   Strand  ID      Seq
#PPBP    684     +       PPBP-A-SP1      TGTTTCTGTGCTATTGCCTTATCAGAGA
#PPBP    801     -       PPBP-A-SP2      GGGTTATCTGTGGGTGTCCTGATC    
#PPBP    653     +       PPBP-B-SP1      ATTGCAACCAAGTCGAAGTGATGTAA  
#PPBP    814     -       PPBP-B-SP2      GCAAAAGAGAACAGGGTTATCTGTGG  
#PPBP    800     +       PPBP-C-SP1      CCTGTTCTCTTTTGCAGAGCCA  
#PPBP    935     -       PPBP-C-SP2      GCAGAAACAGAACAAATTAATCAGCA
my @head;
open(I, $fprimer) or die $!;
while(<I>){
	chomp;
	next if(/^$/);
	if(/^#/){
		push @head, $_;
		next;
	}
	my ($chr, $start, $strand, $id, $seq)=split /\t/, $_;
	$id=~s/-P[12]//;
	my $len = length $seq;
	my $p3=$strand eq "+"? $start+$len: $start-$len;
	$template{$id}=$chr;
	push @{$primer{$id}{"info"}}, $_;
	push @{$primer{$id}{"pos"}}, $p3;
}
close(I);

my %tinfo;
open(T, $ftemplate) or die $!;
$/=">";
while(<T>){
	chomp;
	next if(/^$/);
	my ($head, @seq)=split /\n/;
	my $seq = join("", @seq);
	my ($id)=split /\s+/, $head;
	if($head=~/XP:Z:/){
		my ($strand, $chr, $start, $end)=$head=~/XP:Z:([+-])(\S+):(\d+)-(\d+)/;
		@{$tinfo{$id}}=($chr, $start, $end, $strand);
	}else{
		@{$tinfo{$id}}=($id, 1, length $seq, "+");
	}
}
close(T);


my %probe;
my %probe_info;
my @fprobe=split /,/, $fprobes;
foreach my $f(@fprobe){
	open(P, $f) or die $!;
	$/="\n";
	while(<P>){
		chomp;
		next if(/^$/ || /^#/);
		my ($id, $tid, $len, $dis, $seq, undef, $score_info, $tm, @info)=split /\t/, $_;
		my ($chr, $pos, $strand); ## pos is 5end position
		if(!exists $tinfo{$tid}){ ## design generic primer, rev
			if($tid=~/_rev/){
				$tid=~s/_rev//;		
				$pos = $dis+$len;
				$strand = "-";
				$chr=$tid;
			}else{
				die "$tid can't be found in $ftemplate!\n";
			}
		}else{
			my ($tchr, $tstart, $tend, $tstrand)=@{$tinfo{$tid}};
			($pos, $strand) = &get_chr_info($tstart, $tend, $tstrand, "+", $len, $dis);
			$chr=$tchr;
		}
#		print join("\t","pos:", $id, $tid, $len, $dis,$seq,"=>", $chr, $pos, $strand),"\n";
		@{$probe{$chr}{$pos}{$id}}=($len, $strand, $seq, $tm, $score_info);
		@{$probe_info{$id}}=($chr, $pos, $strand, $seq, $len, $tm, @info);
	}
	close(P);
}

open(O,">$outdir/$fkey.probe") or die $!;
print O "##Score_info(Probe): stm, slen, sdis, sG5, sCGd, spoly, sgc, shairpin, snalign\n";
print O "##ProductSize(Probe): Distance to Primer 3end\n";
print O join("\n", @head),"\n";
foreach my $id (sort{$a cmp $b} keys %primer){
	my ($start, $end) = sort @{$primer{$id}{"pos"}};
	my $tid = $template{$id};
	my %score;
#	if($id eq "rs151242354-D-A"){
#		print "yes\n";
#		print join("\t",$tid, $start, $end),"\n";
#		print join(",", sort {$a<=>$b} keys %{$probe{$tid}}),"\n";
#		print Dumper %{$probe{$tid}};
#		die;
#	}
	foreach my $pos (sort{$a<=>$b} keys %{$probe{$tid}}){
		foreach my $proid (keys %{$probe{$tid}{$pos}}){
			my ($len, $strand, $seq, $tm, $score_info, @info) = @{$probe{$tid}{$pos}{$proid}};
			next if($pos<$start || $pos>$end);
			my @score;
			my ($tstm, $tslen, $tsdis, $tsG5, $tsCGd)=(15, 5, 10, 15, 15);
			my $dis = $strand eq "+"? abs($pos-$start+1): abs($pos-$end+1); ## distance to primer 3end
			my ($is_G5, $CGd) = &G_content($seq);
			push @score, int(&score_single($tm, $tstm, @tm));
			push @score, int(&score_single($len, $tslen, @len));
			push @score, int(&score_single($dis, $tsdis, @dis));
			push @score, $is_G5? 0: $tsG5;
			push @score, int(&score_single($CGd, $tsCGd, @CGd));
			
			
			my ($spoly, $sgc, $shairpin, $snalign)=(split /,/, $score_info)[1,4,7,8]; ## original weight is 14, 2, 20, 16
			my ($tspoly, $tsgc, $tshairpin, $tsnalign)=(5, 5, 20, 10);
			push @score, int($spoly/14*$tspoly);
			push @score, int($sgc/2*$tsgc);
			push @score, int($shairpin/20*$tshairpin);
			push @score, int($snalign/16*$tsnalign);
			
			my $score=$score[0];
			for(my $i=1; $i<@score; $i++){
#				$score[$i]=$score[$i]<0? 0: $score[$i];
				$score+=$score[$i];
			}
			push @{$score{$score}}, $proid;
			push @{$probe_info{$proid}},($dis, $score, join(",",@score)); 
		}
	}
	print O join("\n", @{$primer{$id}{"info"}}),"\n";
	my $n=0;
	foreach my $s(sort {$b<=>$a} keys %score){
		foreach my $proid(@{$score{$s}}){
			$n++;
			last if($n>$nprobe);
			my ($chr, $pos, $strand, $seq, $len, $tm, @info)=@{$probe_info{$proid}}; 
			my $pidnew = $id."-Probe";
			my ($dis, $score, $score_info)=@info[-3,-2,-1];
			pop @info;
			pop @info;
			pop @info;
			print O join("\t", $chr, $pos, $strand, $pidnew, $seq, $len, $dis, "NA", "NA", $score, $score_info, $tm, @info),"\n";
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
sub G_content{
	my ($seq)=@_;
	$seq = uc($seq);
	my @unit = split //, $seq;
	my $is_G5 = $unit[0] eq "G"? 1: 0;
	my ($n, $nc, $ng)=(0,0,0);

	for(my $i=0; $i<@unit; $i++){
		$n++;
		if($unit[$i] eq "C"){
			$nc++;
		}elsif($unit[$i] eq "G"){
			$ng++;
		}
	}
	return($is_G5, ($nc-$ng)/$n);
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


sub USAGE {#
	my $usage=<<"USAGE";
Program:
Version: $version
Contact:zeng huaping<huaping.zeng\@genetalks.com> 

Usage:
  Options:
  -ip  <file>   input primer file, forced
  -it  <file>   input template file, forced
  -io  <file,>  input oligo score files, seperate by ",", forced
  -k  <str>	   Key of output file, forced

  -opttm <int>  optimal tm of probe, [$opt_tm]
  -num   <int>  designed probes number, [$nprobe]
  -minl  <int>  min length of probe, [$min_len]
  -maxl  <int>  max length of probe, [$max_len]
  -od <dir>	   Dir of output file, default ./
  -h		   Help

USAGE
	print $usage;
	exit;
}

