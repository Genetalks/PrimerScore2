#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
require "$Bin/path.pm";
require "$Bin/common.pm";
require "$Bin/snp.pm";
my $BEGIN_TIME=time();
my $version="2.0";
#######################################################################################

# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
our $REF_GRCh37;
our $REF_GRCh37_SNP;
our $SAMTOOLS;
my ($ftarget,$fkey,$outdir);
my $Max_Dis = 20;
my $vcf;
my $fref = $REF_GRCh37;
my $fref_snp;
my $Extend = 200;
my $UD_devide;
my $die_check;
GetOptions(
				"help|?" =>\&USAGE,
				"i:s"=>\$ftarget,
				"r:s"=>\$fref,
				"s:s"=>\$fref_snp,
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
$fref_snp = defined $fref_snp? $fref_snp: $fref;
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
	if($unit[2]=~/^\d+$/ && $unit[2]-$unit[1]<1000 && $unit[2]-$unit[1]>=0){
		if($unit[5]=~/ID=/){
			$ftype = "Target";
		}else{
			die "Wrong SNP file!\n";
		}
	}else{
		$ftype ="VCF";
	}
}
my %target;
my %idcheck;
open(T, $ftarget) or die $!;
while (<T>){
	chomp;
	next if(/^$/);
	my ($chr, $s, $e, $ref, $alt, $info);
	my ($id, $gene, $gstrand, $cpos);
	if($ftype eq "VCF"){
		#13      19685715        rs564281463     G       GA      100     PASS    AC=2263;AF=0.451877;AN=5008;
		($chr, $s, $id, $ref, $alt)=split;
#		$e=$s+length($ref)-1;
		($s, $e, $ref, $alt) = &convert_StartEndRefAlt($s, $ref, $alt);
	}else{
		($chr, $s, $e, $ref, $alt, $info)=split;
		#10      43609948        43609948        T       C       ID=COSM966;GENE=RET;STRAND=+;CDS=c.1900T>C;AA=p.C634R;CNT=13;SOURCE=50Gene
		($id)=$info=~/ID=(\S+?);/;
	}
	if(!defined $s || !defined $e || $e-$s>100){
		print "Wrong: ref base $ref($chr:$s-$e) is not right! Target: $_\n";
		$check_success=0;
		next;
	}

	if(!defined $id || $id eq "" || $id eq "-" || $id eq "NA"){
		$id=join("_", $chr, $s);
	}

	## check ref base
	if($ref ne "-"){
		my $rinfo = `$SAMTOOLS faidx $fref $chr:$s-$e`;
		my (undef, @seq)=split /\n/, $rinfo;
		my $seq = join("", @seq);
		$seq = uc($seq);
		if($ref ne $seq){
			print "Wrong: ref base $ref($chr:$s-$e $seq) of $id is not right!\n";
			$check_success=0;
			next;
		}
	}
	
	$idcheck{$id}++;
	push @{$target_start{$chr}{$s}},[$id];
	push @{$target_end{$chr}{$e}}, [$id];
	@{$target{$id}}=($s, $e, $chr, $ref, $alt);
}
close (T);

#### check id
#foreach my $id(keys %idcheck){
#	if($idcheck{$id}>1){
#		print STDERR "Repeat ID: $id, please correct and run again!\n";
#		die;
#	}
#}

if($check_success==0 && $die_check){
	print "ref base check failed!\n";
	die;
}

my %out;
if(defined $UD_devide){
	open ($out{"U"}, ">$outdir/$fkey.template_U.fa") or die $!;
	open ($out{"D"}, ">$outdir/$fkey.template_D.fa") or die $!;
	open ($out{"Us"}, ">$outdir/$fkey.template_U_snp.fa") or die $!;
	open ($out{"Ds"}, ">$outdir/$fkey.template_D_snp.fa") or die $!;
}else{
	open ($out{"U"}, ">$outdir/$fkey.template.fa") or die $!;
	open ($out{"Us"}, ">$outdir/$fkey.template_snp.fa") or die $!;
}
my %tem_target;
my %tseq;
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
	get_template_seq(\%tem_target, \%tseq, \@s_group, $target_start{$chr}, $chr, "Start");
	get_template_seq(\%tem_target, \%tseq, \@e_group, $target_end{$chr}, $chr, "End");
	
}

## check target template-U and template-D are consistant, maybe not when indel is too long
foreach my $oid(sort{$a cmp $b}keys %target){
	my $tidU = $tem_target{"tar"}{$oid}{"U"};
	my $tidD = $tem_target{"tar"}{$oid}{"D"};
	next if($tidU eq $tidD);
	
	my @idsU=split /\|/, $tidU;
	my @idsD=split /\|/, $tidD;
	my %idsU;
	foreach my $id(@idsU){
		$idsU{$id}=1;
	}
	my @gr;
	my $tidnew=$oid;
	foreach my $id(@idsD){
		if(exists $idsU{$id}){
			push @gr, $id;
			$tem_target{"tar"}{$id}{"U"}=$tidnew;
			$tem_target{"tar"}{$id}{"D"}=$tidnew;
			delete $tem_target{"tem"}{$tidU}{"U"}{$id};
			delete $tem_target{"tem"}{$tidD}{"D"}{$id};
			$tem_target{"tem"}{$tidnew}{"U"}{$id}=1;
			$tem_target{"tem"}{$tidnew}{"D"}{$id}=1;
		}
	}
	@{$tseq{$tidnew}{"U"}}=@{$tseq{$tidU}{"U"}};
	@{$tseq{$tidnew}{"D"}}=@{$tseq{$tidD}{"D"}};
}

foreach my $tid(sort {$a cmp $b} keys %{$tem_target{"tem"}}){
	my @oidU=sort{$a cmp $b} keys %{$tem_target{"tem"}{$tid}{"U"}};
	my @oidD=sort{$a cmp $b} keys %{$tem_target{"tem"}{$tid}{"D"}};
	next if(scalar @oidU==0);
	#my $tidU=join("|", @oidU);
	#my $tidD=join("|", @oidD);
	my $tidU=$oidU[0];
	my $tidD=$oidD[0];
	if($tidU ne $tidD){
		die "U target ids $tidU are not consitant with D target ids $tidD!\n";
	}

	my @s;
	my @e;
	foreach my $id(@oidU){
		push @s, $target{$id}->[0];
		push @e, $target{$id}->[1];
	}
	my @ssort = sort{$a<=>$b}@s;
	my @esort = sort{$a<=>$b}@e;
	my ($Us, $Ue) = ($Min_Dis_Detect, $ssort[-1]-$ssort[0]+$Min_Dis_Detect);
	my ($Ds, $De) = ($Min_Dis_Detect, $esort[-1]-$esort[0]+$Min_Dis_Detect);
	my ($pos_infoU, $seqU, $seqUs)=@{$tseq{$tid}{"U"}};
	my ($pos_infoD, $seqD, $seqDs)=@{$tseq{$tid}{"D"}};
	my $id_str = join(";", @oidU);
	print {$out{"U"}} ">$tidU-U\tXS:i:$Us\tXE:i:$Ue\tXP:Z:$pos_infoU\tXI:Z:$id_str\n";
	print {$out{"U"}} $seqU,"\n";
	print {$out{"U"}} ">$tidU-D\tXS:i:$Ds\tXE:i:$De\tXP:Z:$pos_infoD\tXI:Z:$id_str\n";
	print {$out{"U"}} $seqD,"\n";

	print {$out{"Us"}} ">$tidU-U\tXS:i:$Us\tXE:i:$Ue\tXP:Z:$pos_infoU\tXI:Z:$id_str\n";
	print {$out{"Us"}} $seqUs,"\n";
	print {$out{"Us"}} ">$tidU-D\tXS:i:$Ds\tXE:i:$De\tXP:Z:$pos_infoD\tXI:Z:$id_str\n";
	print {$out{"Us"}} $seqDs,"\n";

}
close($out{"U"});
close($out{"Us"});
if(defined $UD_devide){
	close($out{"D"});
	close($out{"Ds"});
}
#######################################################################################
print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################

# ------------------------------------------------------------------
# sub function
# ------------------------------------------------------------------

sub get_template_seq{
	my ($atemp_target, $aseq, $agroups, $atarget, $chr, $mode)=@_;
	my @groups = @{$agroups};
	
	for(my $i=0; $i<@groups; $i++){
		my @ids;
		for(my $j=0; $j<@{$groups[$i]}; $j++){
			if(!exists $atarget->{$groups[$i][$j]}){
				print $groups[$i][$j],"\n";
				die;
			}
			my @target_info = @{$atarget->{$groups[$i][$j]}};
			for (my $k=0; $k<@target_info; $k++){
				push @ids, $target_info[$k][0];
			}
		}

		## get id
		my $pos = $groups[$i][0];
		my $pos_end = $groups[$i][-1];
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
		my $seq=&faidx_seq($fref, $chr, $start, $end);
		my $seq_snp=&faidx_seq($fref_snp, $chr, $start, $end);
		$seq_snp=&seq_snp_add_target($seq_snp, $atarget, $chr, $start, $end); #add targetspot snp
		my $strand = "+";
		if($flank eq "U"){## targetspot is on the start of sequence
			$seq = &revcom($seq);
			$seq_snp = &revcom($seq_snp);
			$strand = "-";
		}
		my $pos_info = $strand."$chr:$start-$end";
	
		## store
		my $tid = join("|", @ids);
		foreach my $id(@ids){
			$atemp_target->{"tem"}{$tid}{$flank}{$id}=1;
			$atemp_target->{"tar"}{$id}{$flank}=$tid;
		}
		@{$aseq->{$tid}{$flank}}=($pos_info, $seq, $seq_snp);
	}
}


sub seq_snp_add_target{
	my ($seq_snp, $atarget, $chr, $start, $end)=@_;
	my @unit = split //, $seq_snp;
	my @pos;
	foreach my $id(keys %target){
		my ($s, $e, $c, $r, $a)=@{$target{$id}};
		if($c eq $chr && $s<=$end && $e>=$start){
			&sequence_convert_snp(\@unit, $s-$start+1, $r, $a); ## convert @unit directly
		}
	}
	my $seq_new = join("", @unit);
	return $seq_new;
}

sub faidx_seq{
	my ($fref, $chr, $start, $end, $flank)=@_;
	my $seq_info = `$SAMTOOLS faidx $fref $chr:$start-$end`;
	my @line = split /\n/, $seq_info;
	shift @line;
	my $seq = join("", @line);
	return ($seq);
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
  -s  <file>   	Input ref file containing snps, optional
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

