#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
require "$Bin/path.pm";
my $BEGIN_TIME=time();
my $version="2.0";
#######################################################################################

# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
our $REF_HG19;
my ($ftarget,$fkey,$outdir);
my $Max_Dis = 20;
my $vcf;
my $fref = $REF_HG19;
my $Extend = 200;
my $UD_devide;
my $die_check;
GetOptions(
				"help|?" =>\&USAGE,
				"i:s"=>\$ftarget,
				"r:s"=>\$fref,
				"k:s"=>\$fkey,
				"dieCk:s"=>\$die_check,
				"UD_devide:s"=>\$UD_devide,
				"md:s"=>\$Max_Dis,
				"et:s"=>\$Extend,
				"od:s"=>\$outdir,
				) or &USAGE;
&USAGE unless ($ftarget and $fkey);

$outdir||="./";
`mkdir $outdir`	unless (-d $outdir);
$outdir=AbsolutePath("dir",$outdir);

my $Max_Group_Len = $Max_Dis;
my $Min_Dis_Detect = 0;

my %target_start;
my %target_end;
my $check_success=1;
my $ftype;
{
	my $head = `head -1 $ftarget`;
	chomp $head;
	my @unit = split /\s+/, $head;
	if($unit[2]=~/^\d+$/ && $unit[2]-$unit[1]<1000 && $unit[2]-$unit[1]>0){
		if($unit[5]=~/ID=/){
			$ftype = "Target";
		}else{
			die "Wrong SNP file!\n";
		}
	}else{
		$ftype ="VCF";
	}
}
open(T, $ftarget) or die $!;
while (<T>){
	chomp;
	next if(/^$/);
	my ($chr, $s, $e, $ref, $alt, $info);
	my ($id, $gene, $gstrand, $cpos);
	if($ftype eq "VCF"){
		#13      19685715        rs564281463     G       GA      100     PASS    AC=2263;AF=0.451877;AN=5008;
		($chr, $s, $id, $ref, $alt)=split;
		$e = $s+length($ref)-1;
	}else{
		($chr, $s, $e, $ref, $alt, $info)=split;
		#10      43609948        43609948        T       C       ID=COSM966;GENE=RET;STRAND=+;CDS=c.1900T>C;AA=p.C634R;CNT=13;SOURCE=50Gene
		($id)=$info=~/ID=(\S+?);/;
	}
	## check ref base
	if($ref ne "-"){
		my $rinfo = `samtools faidx $fref $chr:$s-$e`;
		my ($id, @seq)=split /\n/, $rinfo;
		my $seq = join("", @seq);
		$seq = uc($seq);
		if($ref ne $seq){
			print "Wrong: ref base $ref($chr:$s-$e $seq) is not right!\n";
			$check_success=0;
		}
	}

	push @{$target_start{$chr}{$s}},[$id];
	push @{$target_end{$chr}{$e}}, [$id];
}
close (T);

if($check_success==0 && $die_check){
	print "ref base check failed!\n";
	die;
}

my %out;
if(defined $UD_devide){
	open ($out{"U"}, ">$outdir/$fkey.template_U.fa") or die $!;
	open ($out{"D"}, ">$outdir/$fkey.template_D.fa") or die $!;
}else{
	open ($out{"U"}, ">$outdir/$fkey.template.fa") or die $!;
}
foreach my $chr(sort {$a cmp $b} keys %target_start){
	my @merge_region;
	my @merge_id = ();
	my $start;
	my $end;

	my $last;
	my @s_sort = sort {$a <=>$b} keys %{$target_start{$chr}};
	my @e_sort = sort {$a <=>$b} keys %{$target_end{$chr}};
	my @s_group = &grouping_by_pos($Max_Dis, $Max_Group_Len, \@s_sort);
	my @e_group = &grouping_by_pos($Max_Dis, $Max_Group_Len, \@e_sort);

	## get template seq .fa
	get_template_seq(\%out, \@s_group, $target_start{$chr}, $chr, "Start");
	get_template_seq(\%out, \@e_group, $target_end{$chr}, $chr, "End");
}

close($out{"U"});
if(defined $UD_devide){
	close($out{"D"});
}
#######################################################################################
print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################

# ------------------------------------------------------------------
# sub function
# ------------------------------------------------------------------

sub get_template_seq{
	my ($aout, $agroups, $atarget, $chr, $mode)=@_;
	my @groups = @{$agroups};
	
	for(my $i=0; $i<@groups; $i++){
		my $info;
		my $id_str;
		for(my $j=0; $j<@{$groups[$i]}; $j++){
			if(!exists $atarget->{$groups[$i][$j]}){
				print $groups[$i][$j],"\n";
				die;
			}
			my @target_info = @{$atarget->{$groups[$i][$j]}};
			for (my $k=0; $k<@target_info; $k++){
				$id_str.=$target_info[$k][0].";";
				$info .= join(",",@{$target_info[$k]}).";";
			}
		}

		## get id 
		my $pos = $groups[$i][0];
		my $pos_end = $groups[$i][-1];
		my ($ds, $de) = ($Min_Dis_Detect, $pos_end-$pos+$Min_Dis_Detect);

		my ($start, $end);
		my $flank;
		if ($mode eq "Start"){
			$end = $pos-$Min_Dis_Detect;
			$start = $end - $Extend;
			if($start <0){
				$start=0;
			}
			$flank = "U";
		}else{
			$start = $pos_end + $Min_Dis_Detect;
			$end = $start + $Extend;
			$flank = "D";
		}
			
		## get seq
		my $seq_info = `samtools faidx $fref $chr:$start-$end`;
		my @line = split /\n/, $seq_info;
		shift @line;
		my $seq = join("", @line);
		my $strand = "+";
		if($flank eq "D"){
			$seq =~tr/atcg/tagc/;
			$seq =~tr/ATCG/TAGC/;
			$seq = reverse $seq;
			$strand = "-";
		}
		my $pos_info = $strand."$chr:$start-$end";
		
		## output
		my $id = $atarget->{$groups[$i][0]}->[0]->[0];
		if(!defined $UD_devide){
			$id.="-".$flank;
			print {$aout->{"U"}} ">$id\tXS:i:$ds\tXE:i:$de\tXP:Z:$pos_info\tXI:Z:$id_str\tXD:Z:$info\n";
			print {$aout->{"U"}} $seq,"\n";
		}else{
			print {$aout->{$flank}} ">$id\tXS:i:$ds\tXE:i:$de\tXP:Z:$pos_info\tXI:Z:$id_str\tXD:Z:$info\n";
			print {$aout->{$flank}} $seq,"\n";
		}
	}
	
}

sub grouping_by_pos{
	my ($Max_Dis, $Max_Group_Len, $apos)=@_;
	my @pos=@{$apos};
	my $last;
	my @final_group;
	my @group;
	for(my $i=0; $i<@pos+1; $i++){
		if (!defined $last || ($i!=@pos && $pos[$i]-$last<=$Max_Dis)){ 
			push @group, $pos[$i];
			$last = $pos[$i];
			next;
		}

		### handle with @group
		my $len = $group[-1] - $group[0];
		if ($len <= $Max_Group_Len){
			push @final_group, [@group];
		}else{
			#print join("\t", @group),"\n";
			my @groups= &break_group(@group);
			#print Dumper @groups;
			push @final_group, @groups;
		}

		### init
		if($i!=@pos){
			@group = ();
			push @group, $pos[$i];
			$last = $pos[$i];
		}
	}

	return @final_group;
}

sub break_group{
	my @group = @_;
	my @final_group;

	my $max_dis = 0;
	for(my $k=0; $k<@group-1; $k++){
		if($group[$k+1] - $group[$k] > $max_dis){
			$max_dis = $group[$k+1] - $group[$k];
		}
	}
	my @group_new1;
	my @group_new2 = @group;
	for(my $k=0; $k<@group-1; $k++){
		if($group[$k+1] - $group[$k] < $max_dis){
			push @group_new1, $group[$k];
			shift @group_new2;
		}elsif($group[$k+1] - $group[$k] == $max_dis){
			push @group_new1, $group[$k];
			shift @group_new2;
			last;
		}
	}


	## break
	if($group_new1[-1]-$group_new1[0]>$Max_Group_Len){
		my @grs = &break_group(@group_new1);
		push @final_group, @grs;
	}else{
		push @final_group, [@group_new1];
	}
	if($group_new2[-1]-$group_new2[0]>$Max_Group_Len){
		my @grs = &break_group(@group_new2);
		push @final_group, @grs;
	}else{
		push @final_group, [@group_new2];
	}

	return @final_group;
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

	###input file format:
	 format 1:
		chr7       107350577       107350577       A       G       ID=SLC26A4_2168A_G;GENE=SLC26A4;STRAND=+;CDS=c.2168A>G;AA=p.H723R;NM=NM_000441.1
		chr7       107323898       107323898       A       G       ID=SLC26A4_919-2A_G;GENE=SLC26A4;STRAND=+;CDS=c.919-2A>G;AA=;NM=NM_000441.1
		chr13      20763691        20763691        -       C       ID=GJB2_35dupG;GENE=GJB2;STRAND=-;CDS=c.35dupG;AA=p.V13fs;NM=NM_004004.5
		chr13      20763421        20763422        AT      -       ID=GJB2_299delAT;GENE=GJB2;STRAND=-;CDS=c.299_300del;AA=p.H100fs;NM=NM_004004.5
	 format 2: like first 5 columns of vcf
		chr1    10247   rs796996180     T       C
		chr1    10248   rs148908337     A       T
		chr1    10249   rs774211241     AAC     A
		chr1    10250   rs199706086     A       C


Usage:
  Options:
  -i  <file>   	Input file, forced
  -r  <file>   	Input ref file, [$fref]
  -k  <str>     Key of output file, forced

  --dieC		die when ref base check failed
  --UD_devide   output U and D template file separately
  -md <int>    	max distance permissible in one target spots cluster between two adjacent target spots , [$Max_Dis]
  -et <int>     extend length of UD flank of one target spot or target spots cluster , [$Extend]
  -od <dir>		Dir of output file, default ./
  -h		 	Help

USAGE
	print $usage;
	exit;
}

