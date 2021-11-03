#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
require "$Bin/path.pm";
require "$Bin/dimer.pm";
my $BEGIN_TIME=time();
my $version="1.0.0";
#######################################################################################

# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
my ($fIn,$fkey,$outdir);
my $type = "Result";
my ($mv, $dv, $dNTP, $dna, $tp, $sc)=(50, 1.5, 0.6, 50, 1, 1);
my $MultiPlex;
my $SelfComplementary;
my $high_tm = 30;
my $min_tm = -30;
my $max_unmap=20;
GetOptions(
				"help|?" =>\&USAGE,
				"i:s"=>\$fIn,
				"t:s"=>\$type,
				"k:s"=>\$fkey,
				"mv:s"=>\$mv,
				"dv:s"=>\$dv,
				"dNTP:s"=>\$dNTP,
				"dna:s"=>\$dna,
				"MultiPlex:s"=>\$MultiPlex,
				"SelfComplementary:s"=>\$SelfComplementary,
				"od:s"=>\$outdir,
				) or &USAGE;
&USAGE unless ($fIn and $fkey);

$outdir||="./";
`mkdir $outdir`	unless (-d $outdir);
$outdir=AbsolutePath("dir",$outdir);

our $PATH_PRIMER3;
my $ntthal = "$PATH_PRIMER3/src/ntthal -mv $mv -dv $dv -n $dNTP -d $dna";

my ($ix1, $ix2);
if($type eq "Result"){
	($ix1, $ix2)=(3,4);
}elsif($type eq "Primer"){
	($ix1, $ix2)=(0,1);
}

my %seq;
my %id;
open(I, $fIn) or die $!;
while(<I>){
	chomp;
	next if(/^$/ || /^#/);
	my ($id, $seq)=(split /\s+/, $_)[$ix1, $ix2];
	my ($tid, $LR)=$id=~/(\S+)-([12LRP])/;
	$LR=~s/L/1/;
	$LR=~s/R/2/;
	$LR=~s/P/3/;
	$seq{$tid}->[$LR-1]=$seq;
	$id{$tid}[$LR-1]=$id;
}
close(I);


open(O,">$outdir/$fkey.cross.dimer")or die $!;
foreach my $tid(sort {$a cmp $b} keys %seq){
	my @seq=@{$seq{$tid}};
	for(my $j=0; $j<@seq; $j++){
		my $oligo_seq = $seq{$tid}->[$j];
		my $id = $id{$tid}->[$j];
#		print O ">",$id,"\t",$oligo_seq,"\n";
		if(defined $SelfComplementary){
			my $HPinfo = `$ntthal -a HAIRPIN -s1 $oligo_seq`;
			chomp $HPinfo;
			if(!defined $HPinfo){
				print $HPinfo;
				my ($HPtm)=$HPinfo=~/t = (\S+)/;
				if($HPtm >= $high_tm){
					print O "#Hairpin\t$id\t$HPtm\n"; 
					print O $HPinfo, "\n";
				}
			}
			
			my $ENDinfo = `$ntthal -a END1 -s1 $oligo_seq -s2 $oligo_seq`;
			chomp $ENDinfo;
			#print "$ntthal -a END1 -s1 $oligo_seq -s2 $oligo_seq\n";
			my ($END_tm, $END31, $END32)=&dimer_amplify($ENDinfo);
			if($END_tm>=$high_tm || ($END_tm>=$min_tm && $END31+$END32<=$max_unmap)){
				print O join("\t", "#ENDDimer", $id, $id, $END_tm, $END31, $END32),"\n";
				print O $ENDinfo, "\n";
			}
			
			my $ANYinfo = `$ntthal -a ANY -s1 $oligo_seq -s2 $oligo_seq`;
			chomp $ANYinfo;
			#print "$ntthal -a ANY -s1 $oligo_seq -s2 $oligo_seq\n";
			my ($ANY_tm, $ANY31, $ANY32)=&dimer_amplify($ANYinfo);
			if($ANY_tm>=$high_tm || ($ANY_tm>=$min_tm && $ANY31+$ANY32<=$max_unmap)){
				print O join("\t", "#ANYDimer", $id, $id, $ANY_tm, $ANY31, $ANY32),"\n";
				print O $ANYinfo, "\n";
			}
			
		}
		## cross primer of same tid
		for(my $k=$j+1; $k<@seq; $k++){
			my ($info, $tm, $umap1, $umap2);
			my $seq = $seq{$tid}->[$k];
			my $id2=$id{$tid}->[$k];

			my $ANYinfo = `$ntthal -a ANY -s1 $oligo_seq -s2 $seq`;
			chomp $ANYinfo;
			#print "$ntthal -a ANY -s1 $oligo_seq -s2 $seq\n";
			my ($ANY_tm, $ANY31, $ANY32)=&dimer_amplify($ANYinfo);
			if($ANY_tm>=$high_tm || ($ANY_tm>=$min_tm && $ANY31+$ANY32<=$max_unmap)){
				print O join("\t", "#ANYDimer", $id, $id2, $ANY_tm, $ANY31, $ANY32),"\n";
				print O $ANYinfo, "\n";
			}
			
			my $END1info = `$ntthal -a END1 -s1 $oligo_seq -s2 $seq`;
			chomp $END1info;
			#print "$ntthal -a END1 -s1 $oligo_seq -s2 $seq\n";
			#print "result1:\n", $END1info,"\n";
			my ($END_tm, $END31, $END32)=&dimer_amplify($END1info);
			if($END_tm>=$high_tm || ($END_tm>=$min_tm && $END31+$END32<=$max_unmap)){
				print O join("\t", "#ENDDimer", $id, $id2, $END_tm, $END31, $END32),"\n";
				print O $END1info, "\n";
			}
			my $END2info = `$ntthal -a END2 -s1 $oligo_seq -s2 $seq`;
			chomp $END2info;
			#print "$ntthal -a END2 -s1 $oligo_seq -s2 $seq\n";
			{
				my ($END_tm, $END31, $END32)=&dimer_amplify($END2info);
				if($END_tm>=$high_tm || ($END_tm>=$min_tm && $END31+$END32<=$max_unmap)){
					print O join("\t", "#ENDDimer", $id, $id2, $END_tm, $END31, $END32),"\n";
					print O $END2info, "\n";
				}
			}
			
		}


		## cross tid 
		if(defined $MultiPlex){
			foreach my $tid2(keys %seq){
				next if($tid2 eq $tid);
				my @seq2=@{$seq{$tid2}};
				for(my $k=0; $k<@seq2; $k++){
					my $seq = $seq{$tid2}->[$k];
					my $id2=$id{$tid2}->[$k];
					if(!defined $seq){
						print $tid2,"\n";
						print Dumper @seq2;
						print $k,"\n";
						die;
					}

					my $ANYinfo = `$ntthal -a ANY -s1 $oligo_seq -s2 $seq`;
					chomp $ANYinfo;
				#	print "$ntthal -a ANY -s1 $oligo_seq -s2 $seq\n";
					my ($ANY_tm, $ANY31, $ANY32)=&dimer_amplify($ANYinfo);
					if(!defined $ANY_tm){
						print $ANYinfo,"\n";
						die;
					}
					if($ANY_tm>=$high_tm || ($ANY_tm>=$min_tm && $ANY31+$ANY32<=$max_unmap)){
						print O join("\t", "#ANYDimer", $id, $id2, $ANY_tm, $ANY31, $ANY32),"\n";
						print O $ANYinfo, "\n";
					}

					my $END1info = `$ntthal -a END1 -s1 $oligo_seq -s2 $seq`;
					chomp $END1info;
				#	print "$ntthal -a END1 -s1 $oligo_seq -s2 $seq\n";
					my ($END_tm, $END31, $END32)=&dimer_amplify($END1info);
					if($END_tm>=$high_tm || ($END_tm>=$min_tm && $END31+$END32<=$max_unmap)){
						print O join("\t", "#ENDDimer", $id, $id2, $END_tm, $END31, $END32),"\n";
						print O $END1info, "\n";
					}
					my $END2info = `$ntthal -a END2 -s1 $oligo_seq -s2 $seq`;
					chomp $END2info;
				#	print "$ntthal -a END2 -s1 $oligo_seq -s2 $seq\n";
					{
						my ($END_tm, $END31, $END32)=&dimer_amplify($END2info);
						if($END_tm>=$high_tm || ($END_tm>=$min_tm && $END31+$END32<=$max_unmap)){
							print O join("\t", "#ENDDimer", $id, $id2, $END_tm, $END31, $END32),"\n";
							print O $END2info, "\n";
						}
					}
				}
			}
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

	-t "Primer": 2 columns (IDs and Sequences), as following:
	xxx-1  AAAAAAAAAAAAAAAAA
	xxx-2  AAAAAAAAAAAAAAAAA
	xxx-P  AAAAAAAAAAAAAAAAA

Usage:
  Options:
  -i  <file>   Input primer file, forced
  -t  <str>    primer file type, "Result" (xxx.final.result) or "Primer", [$type]
  -k  <str>	   Key of output file, forced
  --SelfComplementary      Evalue self-complementary
  --MultiPlex              Evalue dimer among different primer pairs
  -mv        <int>      concentration of monovalent cations in mM, [$mv]
  -dv        <float>    concentration of divalent cations in mM, [$dv]
  -dNTP      <float>    concentration of deoxynycleotide triphosphate in mM, [$dNTP]
  -dna       <int>      concentration of DNA strands in nM, [$dna]

  -od <dir>	Dir of output file, default ./
  -h		 Help

USAGE
	print $usage;
	exit;
}

