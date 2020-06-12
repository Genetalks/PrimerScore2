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
my ($fprimer, $fkey,$detail,$outdir);
my $NoSpecificity;
my $pnum = 1;
my $min_tm = 0;
my $max_tm = 100;
my $min_gc = 0;
my $max_gc = 1;
my $min_tm_spec = 40; #when caculate specificity
my $eff_end1 = 0.2; # PCR ratio on the base of eff_tm when the last 1 base on 3end is not matched
my @rank_end=(3,   5,   8); #  PCR efficiency when dis of mismatch pos to 3end <= @rank_end
my @eff_end =(0.1, 0.4, 0.8 );
my @rank_tm= (55, 50,  45,   40);
my @eff_tm = (1, 0.6, 0.1, 0.05); # PCR efficiency when tm >= @rank_tm
my @rank_dis= (20, 50, 100, 200, 400, 600);
my @eff_dis = (1, 0.8, 0.6, 0.4, 0.2, 0.1); # PCR efficiency when dis to the last PCR product <= (20, 50, 100, 200, 400, 600)
my $eff_times = 10;
my $min_eff = 0.1;
my $nohead;
my $face_to_face;
my $thread = 3;
our $PATH_PRIMER3;
our $REF_HG19;
my $fdatabase = $REF_HG19;
GetOptions(
				"help|?" =>\&USAGE,
				"p:s"=>\$fprimer,
				"d:s"=>\$fdatabase,
				"n:s"=>\$pnum,
				"k:s"=>\$fkey,
				"NoSpecificity:s"=>\$NoSpecificity,
				"nohead:s"=>\$nohead,
				"face_to_face:s"=>\$face_to_face,
				"mintm:s"=>\$min_tm,
				"maxtm:s"=>\$max_tm,
				"mingc:s"=>\$min_gc,
				"maxgc:s"=>\$max_gc,
				"stm:s"=>\$min_tm_spec,
				"Detail:s"=>\$detail,
				"thread:s"=>\$thread,
				"od:s"=>\$outdir,
				) or &USAGE;
&USAGE unless ($fprimer and $fkey);

$outdir||="./";
`mkdir $outdir`	unless (-d $outdir);
$outdir=AbsolutePath("dir",$outdir);
my $min_end_match = 6;
my $merge_len = 100;
my $Wind_GC = 8;

if(defined $face_to_face){
	@rank_dis = (10, 40,  70, 170, 270, 500, 600); # 100-200 is best
	@eff_dis = (0.1, 0.3, 0.8, 1,  0.8, 0.3, 0.1);
}

my $oligotm = "$PATH_PRIMER3/src/oligotm";
my $ntthal = "$PATH_PRIMER3/src/ntthal";
my $primer3_config = "$PATH_PRIMER3/src/primer3_config/";

## creat primer.fa
my @PN=();
for (my $i=0; $i<$pnum; $i++){
	open ($PN[$i], ">$outdir/$fkey.primer_$i.fa") or die $!;
}
open(F, ">$outdir/$fkey.filtered_by_Tm_GC.list") or die $!;
my %evalue;
my %seq;
open (P, $fprimer) or die $!;
while (<P>){
	chomp;
	next if(/^$/);
	$_=~s/\s+$//;
	my ($id, @seq)=split /\s+/, $_;
	@{$seq{$id}}=@seq;
	for (my $i=0; $i<$pnum; $i++){
		my $primer_seq = $seq[$i];
		## len Tm GC Hairpin, Dimer between itself
		if(!defined $primer_seq){
			print $_,"\n";
			die;
		}
		my $Tm = `$oligotm $primer_seq`;
		chomp $Tm;
		$Tm = sprintf "%0.2f", $Tm;
		my @GC_info = &GC_info_stat($primer_seq);
		if($Tm<$min_tm || $Tm>$max_tm){
			print F join("\t", $id, $primer_seq, $Tm, $GC_info[0]),"\n";
			next;
		}

		my $len = length($primer_seq);
		my ($hairpin_dg, $hairpin_tm, $dimer_dg, $dimer_tm)=('','','','');
		
		my $hairpin_result = `$ntthal -path $primer3_config -a HAIRPIN -s1 $primer_seq`;
		($hairpin_dg, $hairpin_tm)=$hairpin_result=~/dG = ([\d\-\.]+)\tt = ([\d\-\.]+)/;
		if (!defined $hairpin_dg){
			$hairpin_dg = '';
			$hairpin_tm = '';
		}else{
			$hairpin_dg =  sprintf "%0.2f",$hairpin_dg;
			$hairpin_tm =  sprintf "%0.2f",$hairpin_tm;
		}
		
		my $dimer_result = `$ntthal -path $primer3_config -a END1 -s1 $primer_seq -s2 $primer_seq`;
		($dimer_dg, $dimer_tm)=$dimer_result=~/dG = ([\d\-\.]+)\tt = ([\d\-\.]+)/;
		if (!defined $dimer_dg){
			$dimer_dg = '';
			$dimer_tm = '';
		}else{
			$dimer_dg =  sprintf "%0.2f",$dimer_dg;
			$dimer_tm = sprintf "%0.2f",$dimer_tm;
		}
		
		push @{$evalue{$id}},[$len, $Tm, @GC_info, $hairpin_dg, $hairpin_tm, $dimer_dg, $dimer_tm];
		
		my $id_new = "$id\_$i";
		print {$PN[$i]} ">$id_new\n";
		print {$PN[$i]} $seq[$i],"\n";
	}
}
for (my $i=0; $i<$pnum; $i++){
	close($PN[$i]);
}
close(F);
exit(0) if(scalar keys %evalue==0);

if(defined $detail){
	open (Detail, ">$outdir/$fkey.evaluation.detail") or die $!;
}
my $DB;
my %speci;
my %efficiency;
if(!defined $NoSpecificity){
	for (my $pn=0; $pn<$pnum; $pn++){
		if($pn!=$pnum-1){
			open ($DB, ">$outdir/$fkey.database_$pn") or die $!;
		}
		
		my %db_region;
		### bwa
		my $fa_primer = "$outdir/$fkey.primer_$pn.fa";
		
		Run("bwa mem -D 0 -k 9 -t $thread -c 5000000 -y 100000 -T 12 -B 1 -L 2,2 -h 200 -a  $fdatabase $fa_primer > $fa_primer.sam");
	#	Run("bwa mem -D 0 -k 9 -t 4 -c 5000000 -y 100000 -T 12 -B 1 -O 2,2 -L 1,1 -h 200 -a  $fdatabase $fa_primer > $fa_primer.sam");
		
		### read in sam
		my %mapping;
		open (I, "$fa_primer.sam") or die $!;
		while (<I>){
			chomp;
			next if(/^$/ || /^\@/ || /^\[/);
			my ($id, $flag, $chr, $pos, $score, $cigar, undef, undef, undef, $seq)=split /\s+/,$_;
			my ($is_unmap, $is_reverse)=&explain_bam_flag($flag);
			my ($md)=$_=~/MD:Z:(\S+)/;
			if ($is_unmap){
				push @{$speci{$id}}, [0,0,0,'']; #[$align_num, $high_tm_num, $end_match_num,$high_tm_info];
				print "Warn: primer $id is not mapped!\n";
				print $_,"\n";
				next;
			}
			#next if($pn!=0 && $is_reverse);
			my ($H3)=&get_3end1_mismatch($is_reverse, $cigar);
			next if($H3>1);
			push @{$mapping{$id}},[$is_reverse, $flag, $chr, $pos, $score, $cigar, $md];
		}
		close(I);
		### evaluate
		foreach my $id (sort {$a cmp $b} keys %mapping){
			my ($ido)=$id=~/(\S+)_$pn/;
			my $align_num = 0;
			my $map_num = scalar @{$mapping{$id}};
			my @high_tm;
			my %high_info;
			for (my $i=0; $i<$map_num; $i++){
				my ($is_reverse, $flag, $chr, $pos, $score, $cigar, $md)=@{$mapping{$id}->[$i]};
				my $primer_seq = $seq{$ido}->[$pn];

				## specificity
				my ($dG, $tm, $end_match, $ntthal_result) = &END_tm($primer_seq, $chr, $pos, $is_reverse, $fdatabase);
				next if($end_match == -1); ## map to NNNN region, invalid

				## get pcr efficiency
				my ($last_eff, $dis);
				if($pn==0){
					$last_eff = 1;
					$dis = 0;
				}else{
					$last_eff = $efficiency{$chr}*$eff_times; ## eff_times: few DNA in the last step will be amplified, so multiply by eff_times
					$dis = $pos;
				}
				my ($eff)=&get_pcr_efficienty($tm, $dis, $cigar, $md, $is_reverse, $min_tm_spec, $last_eff, $eff_end1, \@rank_tm, \@eff_tm, \@rank_dis, \@eff_dis, \@rank_end, \@eff_end);

				$align_num++;
				if(($tm<$min_tm_spec && $align_num>200) || ($align_num>400)){
					$align_num.="+";
					last;
				}
				
				## store info
				my $strand = $is_reverse? "-": "+";
				if($pn==0){
					push @{$high_info{sprintf("%.2f",$tm)}},[$strand, $chr, $pos, $cigar, $md, $end_match, $eff];
				}else{
					push @{$high_info{sprintf("%.2f",$eff)}},[$strand, $chr, $pos, $cigar, $md, $end_match, $tm];
				}

				## out database
				if ($tm >= $min_tm_spec){
					if ($pn != $pnum-1){
						my ($p3, $seq) = &get_database_seq($primer_seq, $chr, $pos, $is_reverse, $cigar, $fdatabase); #($primer, $fdatabae, $chr, $pos, $is_reverse)
						$chr=$is_reverse? "-".$chr: "+".$chr;
						my $id_db = join(",", $id, $chr, $p3);
						$efficiency{$id_db}=$eff;
						print $DB ">$id_db\n$seq\n";
					}
				}
				if(defined $detail){
					print Detail join("\t",$primer_seq, join(",",$chr, $pos, $is_reverse),$dG, $tm, $end_match),"\n";
					print Detail $ntthal_result;
				}
			}

			my ($high_tm_num, $high_eff_num, $high_info);
			if($pn==0){
				($high_tm_num, $high_eff_num, $high_info)=&get_high_info(\%high_info, $min_tm_spec, $min_eff, $pn);
			}else{
				($high_eff_num, $high_tm_num, $high_info)=&get_high_info(\%high_info, $min_eff, $min_tm_spec, $pn);
			}

			if(defined $detail){
				print Detail "##",join("\t", $id, $pn, $align_num, $high_tm_num, $high_eff_num, $high_info),"\n";
			}

			push @{$speci{$ido}}, [$align_num, $high_tm_num, $high_eff_num, $high_info];
		}
		
		if($pn !=$pnum-1){
			close($DB);
			$fdatabase = "$outdir/$fkey.database_$pn";
			if(-f "$fdatabase.amb"){ ## rm old bwa indexed files 
				`rm $fdatabase.*`;
			}
			Run("bwa index $fdatabase");
		}
		
	}
	if(defined $detail){
		close (Detail);
	}
}

### output
open (O, ">$outdir/$fkey.evaluation.out") or die $!;
if(!defined $nohead){
	print O "##High_Info(n=1): TM        : Flag/DatabaseID/Pos/Cigar/MD/End_match_num/Efficiency\n";
	if($pnum>1){
		print O "##High_Info(n>1): Efficiency: Flag/DatabaseID/Pos/Cigar/MD/End_match_num/TM\n";
	}
	print O "#ID";
	for (my $i=0; $i<$pnum; $i++){
		print O "\tSeq\tLen\tTm\tGC\tGC5\tGC3\tdGC\tdG_Hairpin\tTm_Hairpin\tdG_Dimer\tTm_Dimer\tAlign_Num\tHigh_Tm_Num\tHigh_Efficiency_Num\tHigh_Info";
	}
	print O "\n";
}

foreach my $id (sort {$a cmp $b} keys %evalue){
	print O $id;
	for (my $i=0; $i<$pnum; $i++){
		my @spe = defined $NoSpecificity? ("","","",""):@{$speci{$id}->[$i]};
		print O "\t",$seq{$id}->[$i],"\t",join("\t",@{$evalue{$id}->[$i]}),"\t",join("\t",@spe);
	}
	print O "\n";
}
close(O);
#######################################################################################
print STDOUT "\nDone. Total elapsed time : ",time()-$BEGIN_TIME,"s\n";
#######################################################################################

# ------------------------------------------------------------------
# sub function
# ------------------------------------------------------------------

sub get_high_info{
	my ($ahigh_info, $min_key, $min_other, $pn)=@_;
	my %high_info=%{$ahigh_info};
	
	my @high_info=();
	my @high_value=();
	my $n = 0;
	my $high_key_num = 0;
	my $high_other_num=0;
	
	my %record;
	foreach my $t (sort{$b<=>$a} keys %high_info){
		last if($t<$min_key && $n>=3);
		for(my $i=0; $i<@{$high_info{$t}}; $i++){
			if($pn>0){
				my ($strand, $database_id, $pos, $cigar)=@{$high_info{$t}->[$i]}[0..3];
				#check strand 
				next if((!defined $face_to_face && $strand eq "-") || (defined $face_to_face && $strand eq "+"));
				
				my ($idt, $chrt_ori, $p3t)=split /,/, $database_id;
				my ($idto)=$idt=~/(\S+)\_\d+$/;
				my ($strandt,$chrt)=$chrt_ori=~/([+-])(\S+)/;

				my $plen = length ($seq{$idto}->[$pn-1]);
				my ($unmap_len5) = $cigar=~/^(\d+)H/; ## when primer 5end is not mapped, indicate that: the primer is connected to the last primer, it is not expected
				if(!defined $unmap_len5){
					$unmap_len5=0;
				}else{
					if(!defined $face_to_face && $pos==1){
						print "Wrong: the last primer $database_id been connected on the 5end\n";
					}
				}
				if(!defined $face_to_face){
					$p3t = $strandt eq "+"? $p3t+$pos-$plen-$unmap_len5: $p3t-$pos+$plen+$unmap_len5;
				}else{
					$p3t = $strandt eq "+"? $p3t+$pos-$plen: $p3t-$pos+$plen;
				}
				#print join("\t","pos3:", $database_id, $p3t, $pos, $plen, $unmap_len5),"\n";
				next if($record{$chrt."_".$p3t});
				$record{$chrt."_".$p3t}=1;
			}
			if($t >= $min_key){
				$high_key_num++;
			}
			my ($other) = @{$high_info{$t}->[$i]}[-1];
			if ($other >=$min_other){
				$high_other_num++;
			}
			$n++;
			push @high_value, $t;
			push @high_info, join("/",@{$high_info{$t}->[$i]});
		}
	}
	my $high_info = join(",",@high_value).":".join(";",@high_info);
	return($high_key_num, $high_other_num, $high_info);
}

sub get_pcr_efficienty{ 
	my ($tm, $dis, $cigar, $md, $is_reverse, $min_tm, $last_eff, $eff_end1, $arank_tm, $aeff_tm, $arank_dis, $aeff_dis, $arank_end, $aeff_end)=@_;
	if($tm < $min_tm){
		return 0;
	}
	# tm eff
	my $eff_tm = &get_eff_rank($tm, $arank_tm, $aeff_tm, ">=");
	# dis eff
	my $eff_dis = &get_eff_rank($dis, $arank_dis, $aeff_dis, "<=");

	# mismatch in 3'end(the last 1 base)
	my $H3 = &get_3end1_mismatch($is_reverse, $cigar);
	my $eff1 = $H3==0? 1: $eff_end1;

	# mismatch pos to 3end
	my @mis_pos;
	&get_3end_mismatch($md, $is_reverse, \@mis_pos, $arank_end->[-1]);
	my $eff_end = 1;
	for(my $i=0; $i<@mis_pos; $i++){
		$eff_end *= &get_eff_rank($mis_pos[$i], $arank_end, $aeff_end, "<=");
	}

	my $eff = $last_eff * $eff_tm * $eff_dis * $eff1 * $eff_end;
#	print join("\t", $eff, $last_eff, $eff_tm, $eff_dis, $eff1, $eff_end, ":", $tm, $dis, $cigar, $md, $is_reverse),"\n";

	#$eff=$eff>1? 1: $eff;
	return $eff;
}

sub get_3end_mismatch{
	my ($md, $is_reverse, $aresult, $end_len)=@_;
	$md=~s/\^//g; ## indel, usually in the middle of primer seq, so not important
	my @num=split /[ATCG]+/, $md;
	
	if($is_reverse==0){
		my $sum=0;
		for(my $i=$#num; $i>=0; $i--){
			$sum+=!defined $num[$i]? 0: $num[$i];
			$sum++;
			if($sum<=$end_len){
				push @{$aresult}, $sum;
			}
		}
	}else{
		my $sum=0;
		for(my $i=0; $i<@num; $i++){
			$sum+=!defined $num[$i]? 0: $num[$i];
			$sum++;
			if($sum<=$end_len){
				push @{$aresult}, $sum;
			}
		}
	}
}

sub get_eff_rank{
	my ($v, $arank, $aeff, $compare)=@_;
	my @rank = @{$arank};
	my @eff = @{$aeff};
	
	my $eff;
	for(my $i=0; $i<@rank; $i++){
		if($compare eq ">="){
			if($v >= $rank[$i]){
				$eff = $eff[$i];
				last;
			}
		}else{
			if($v <= $rank[$i]){
				$eff = $eff[$i];
				last;
			}
		}
	}
	if(!defined $eff){
		$eff = 0;
	}
	return $eff;

}
sub get_3end1_mismatch{
	my ($is_reverse, $cigar)=@_;
	my $H3;
	if($is_reverse==0){
		($H3)=$cigar=~/M(\d+)H$/;
	}else{
		($H3)=$cigar=~/^(\d+)H/;
	}
	if(!defined $H3){ ## 1H can be prc from experiment data
		$H3=0;
	}
	return $H3;
}

sub get_database_seq{
	my ($primer, $chr, $pos, $is_reverse, $cigar, $fref)=@_;
	my $Extend = 600;
	my $H5len;
	my $plen = length $primer;
	my ($s, $e, $p3);
	if ($is_reverse == 0){
		($H5len)=$cigar=~/^(\d+)H/;
		$H5len = defined $H5len? $H5len: 0;
		$s = $pos+$plen-$H5len; # primer 3' postion
		$e = $s+$Extend;
		$p3=$s;
	}else{
		($H5len)=$cigar=~/(\d+)H$/;
		$H5len = defined $H5len? $H5len: 0;
		$e = $pos-1;
		$s = $e - $Extend;
		$p3 = $e;
	}
	$s=$s>1? $s: 1;

	# get seq
	my $rout = `samtools faidx $fref $chr:$s-$e`;
	chomp $rout;
	my ($id, @seq)=split /\n/, $rout;
	my $seq = join("", @seq);
	if($is_reverse){
		$seq=~tr/ATCGatcg/TAGCtagc/;
		$seq=reverse $seq;
	}
	$seq=$primer.$seq;

	return ($p3, $seq);
}
sub GC_stat{
   	my @u=@_;
	my $total = 0;
	my $gc = 0;
	foreach $b (@u){
		$total++;
		if($b eq 'G' || $b eq 'C' || $b eq "g" || $b eq "c"){
			$gc++;
		}
	}
	return($gc/$total);
}

sub GC_info_stat{
    my ($p)=@_;
    my @u = split //, $p;

	my $GC = &GC_stat(@u);
	my ($GC5, $GC3);
	my $min_GC=1;
	my $max_GC=0;
	for(my $i=0; $i<@u-$Wind_GC+1; $i++){
		my @sub = @u[$i..$i+$Wind_GC-1];
		my $sub_gc = &GC_stat(@sub);
		if($i==0){
			$GC5 = $sub_gc;
		}
		if($i==@u-$Wind_GC){
			$GC3 = $sub_gc;
		}
		if($sub_gc>$max_GC){
			$max_GC = $sub_gc;
		}
		if($sub_gc<$min_GC){
			$min_GC = $sub_gc;
		}
	}
	$GC = sprintf "%0.2f", $GC;
    return ($GC, $GC5, $GC3, $max_GC-$min_GC);
}

sub Run{
    my ($cmd, $log_fh)=@_;

    if (defined $log_fh) {
        print $log_fh $cmd,"\n";
    }
    print STDERR $cmd, "\n";
    my $ret = system($cmd);
    if ($ret) {
        die "Run $cmd failed!\n";
    }
}

sub END_tm{
	my ($primer_seq, $chr, $pos, $is_reverse, $fdatabase)=@_;
	my $len = length($primer_seq);
	my $extend = 10;
	
	### get seq
	#my $seq = join("", split /\n/,`samtools faidx $fdatabase $chr:($pos-$extend)\-($pos+$extend)`);
	my $start = $pos-$extend>1? $pos-$extend: 1;
	my $end = $pos+$len+$extend;
	my $seq_info = `samtools faidx $fdatabase $chr:$start-$end`;
	my @seq_info = split /\n/, $seq_info;
	shift @seq_info;
	my $seq = join("",@seq_info);
	if($seq!~/[ATCGatcg]+/){
		return (-1,-1,-1,-1);
	}
	$seq=uc($seq);
	if ($is_reverse == 0){
		$seq=~tr/atcg/tagc/;
		$seq=~tr/ATCG/TAGC/;
		$seq = reverse $seq;
	}
	### ntthal
	if(defined $detail){
		print Detail "$ntthal -path $primer3_config -a END1 -s1 $primer_seq -s2 $seq\n";
	}
	my $result = `$ntthal -path $primer3_config -a END1 -s1 $primer_seq -s2 $seq`;

	### get tm and end match num
	my ($dG, $tm) = $result =~/dG = ([\d\+\-\.]+)\tt = ([\d\+\-\.]+)/;
	my @line = split /\n/, $result;
	if(!defined $line[2]){ ## just print and return
		print join("\t",$primer_seq, $chr, $pos, $is_reverse, $fdatabase),"\n";
		print "samtools faidx $fdatabase $chr:$start-$end\n";
		print "$ntthal -path $primer3_config -a END1 -s1 $primer_seq -s2 $seq\n";
		return (-1,-1,-1,-1);
	}
	my @aunit = split //, $line[2];
	my $end_match = 0;
	my $count_flag = 0;
	for (my $i=$#aunit; $i>=0; $i--){
		if ($aunit[$i] ne " "){
			$count_flag = 1;
		}
		if ($count_flag ==1){
			if($aunit[$i] ne " "){
				$end_match++;
			}else{
				last;
			}
		}
	}

	return($dG, $tm, $end_match, $result);
}

sub explain_bam_flag{
	my ($flag)=@_;
	my $flag_bin=sprintf("%b", $flag);
	my @flag_bin = split //, $flag_bin;
#	my $is_read1 = @flag_bin>=7? $flag_bin[-7]:0;
#	my $is_read2 = @flag_bin>=8? $flag_bin[-8]: 0;
#	my $is_supplementary = @flag_bin>=12? $flag_bin[-12]: 0;
#	my $is_proper_pair = @flag_bin>=2? $flag_bin[-2]:0;
	my $is_reverse = @flag_bin>=5? $flag_bin[-5]: 0;
	my $is_unmap = @flag_bin>=3? $flag_bin[-3]:0;
#	my $is_munmap = @flag_bin>=4? $flag_bin[-4]:0;
#	my $dup = @flag_bin>=11? $flag_bin[-11]: 0;
#	my @result = ($is_read1, $is_proper_pair, $is_reverse, $is_unmap, $is_munmap, $is_supplementary);
	return ($is_unmap, $is_reverse);
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
	
	v2: 1)stm from 45 to 40
	    2)database seq(amplified from the previous primer) add the previous primer seq
		3)add efficiency 

Usage:
  Options:
  -p  <file>   Input primer list file, forced
  -d  <file>   Input database file, [$fdatabase]
  -n  <int>    combined primer num, single primer: 1, primer pair: 2, forced
  -k  <str>	Key of output file, forced
  
  --NoSpecificity   not evalue specificity
  --face_to_face	evalue face-to-face primers
  -mintm <int>		min tm to evalue, [0]
  -maxtm <int>		max tm to evalue, [100]
  -mingc <float>		min gc to evalue, [0]
  -maxgc <float>		max gc to evalue, [1]
  -stm  <int>   min tm to be High_TM when caculate specificity, [40]
  -thread  <int>   thread in bwa, [$thread]
  --Detail     Output Detail Info, optional
  -od <dir>	Dir of output file, default ./
  -h		 Help

USAGE
	print $usage;
	exit;
}
