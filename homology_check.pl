#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
require "$Bin/path.pm";
my $BEGIN_TIME=time();
my $version="1.0.0";
#######################################################################################

# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
our $REF_HG19;
our $SAMTOOLS;
my ($fIn,$fkey,$outdir);
my ($ftemplate, $fpsl);
my $fref = $REF_HG19;
my $min_ratio = 0.8;
my $max_gap_num = 3;
my $max_gap_len = 20;
GetOptions(
				"help|?" =>\&USAGE,
				"it:s"=>\$ftemplate,
				"ir:s"=>\$fref,
				"k:s"=>\$fkey,
				"mr:s"=>\$min_ratio,
				"mgn:s"=>\$max_gap_num,
				"mgl:s"=>\$max_gap_len,
				"od:s"=>\$outdir,
				) or &USAGE;
&USAGE unless ($ftemplate  and $fkey);

$outdir||="./";
`mkdir $outdir`	unless (-d $outdir);
$outdir=AbsolutePath("dir",$outdir);
our $BLAT;

my $sh = "$BLAT $fref $ftemplate $outdir/$fkey.psl";
print $sh,"\n";
`$sh`;

open(I, $ftemplate) or die $!;
my %pseq;
my %idout;
$/=">";
while(<I>){
	chomp;
	next if(/^$/);
	my ($head, @seq)=split /\n/,$_;
	my $seq=join("", @seq);
	my ($id)=split /\s+/, $head;
	$seq=uc($seq);
	$pseq{$id}=$seq;
	$idout{$id}=$id;
	if(/XP:Z:/){
		my ($ori, $chr, $s, $e)=$head=~/XP:Z:([+-])(\w+):(\d+)-(\d+)/;
		my $idn = $id."(".$ori.$chr.":".$s."-".$e.")";
		$idout{$id}=$idn;
	}

}
close(I);

open(O, ">$outdir/$fkey.homology.seq") or die $!;
print O "#templateID\ttemplateStart\ttemplateEnd\tDatabaseChrAndOrient\tDatabaseStart\tDatabaseEnd\tMatch\tMismatch\tGapNumT\tGapLengthT\tGapNumD\tGapLengthD\n";
open(I, "$outdir/$fkey.psl") or die $!;
$/="\n";
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

while(<I>){
	chomp;
	next if($_!~/^\d/);
	my ($match, $mismatch, $gap1,$gaplen1, $gap2, $gaplen2, $ori, $pid, $len, $ps, $pe, $chr, $s, $e)=(split /\s+/, $_)[0,1,4,5,6,7,8,9,10, 11,12,13,15,16]; #psl v1
#	my ($match, $mismatch, $gap1, $gap2, $ori, $pid, $len, $ps, $pe, $chr, $s, $e)=(split /\s+/, $_)[0,1,4,6,8,13,14,15,16,9,11,12]; #v2
#	next if($match == $len && $mismatch==0 && $gap1==0 && $gap2==0); 
	next if($match<$len*$min_ratio || $gap1>$max_gap_num || $gaplen1>$max_gap_len || $gap2>$max_gap_num || $gaplen2>$max_gap_len);
#	next if($ori eq $ori0 && $chr eq $chr0 && $s==$s0-1 && $e==$e0);
	
	$s++;
	my $fa=`$SAMTOOLS faidx $fref $chr:$s-$e`;
	my ($head, @line) = split /\n/, $fa;
	my $seq = join("", @line);
	$seq = uc($seq);
	if($ori eq "-"){
		$seq=~tr/ATCG/TAGC/;
		$seq=reverse $seq;
	}
	my $minfo = `perl $Bin/sw.pl $pseq{$pid} $seq -print_overall`;
	#print "perl /home/zenghp/bin/sw.pl -print_overall $pseq{$pid} $seq\n";
	chomp $minfo;
	my @minfo = split /\n/, $minfo;

	print O ">",join("\t",$idout{$pid},$ps,$pe, $chr.$ori, $s, $e, $match, $mismatch, $gap1, $gaplen1, $gap2, $gaplen2),"\n";
	print O join("\n", @minfo[0..2]),"\n";
}
close(I);
close(O);



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
  -it  <file>   Input template fa file, forced
  -ir  <file>   reference fasta file, [$fref]
  -k   <str>	Key of output file, forced
  -mr  <float>  min match ratio, [$min_ratio]
  -mgn <int>    max gap num, [$max_gap_num]
  -mgl <int>    max gap length, [$max_gap_len]
  -od <dir>	Dir of output file, default ./
  -h		 Help

USAGE
	print $usage;
	exit;
}

