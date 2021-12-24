# PrimerScore2: a robust high-throughput primer design tool that designs multiple types of primers in one click

## Introduction
PrimerScore2 is a robust high-throughput primer design tool that designs multiple types of primers in one click, including face-to-face, back-to-back, and nested primers, as well as evenly full-covered primers on target regions. PrimerScore2 precisely evaluates candidate oligos′features, containing specificity, SNPs, dimers, and other basic features. Then it scores the oligos and pairs according to the features using growth curve piecewise function and weighted sum model. Finally, it provides the highest-scoring primer pairs or probes and examines for cross-products and cross-dimers among the outputted primers.

## Docker use (RECOMMENDED!)
PrimerScore2 can be run as a [Docker](https://www.docker.com) image. In this way you only need to [install Docker](https://docs.docker.com/install/).

Now users have two options of PrimerScore2 use in docker:  

(1) with already uploaded human reference genome GRCh37 and GRCh38,  
(2) without any reference genomes.  

The 1st variant is idead for use with GRCh37 and GRCh38 genome, but you will have to download about 4 Gb of data. The 2nd is ideal for use with other reference genomes, including other organisms. In this case you will have to download about only 0.5 Gb, but you also have to download database and reference genome files manually.

To use 1st variant, download docker image of PrimerScore2 as follows:

1. Download docker image  
   `docker pull zenghp88/primerscore2:latest`

2. Into the virtual machine command line  
   `docker run -it --entrypoint 'bash' --name primerscore2_ref -v '<directory where you are going to design new primers>:<name of this directory in the container>' zenghp88/primerscore2:latest`  
   where -v option lets you to mount some of your local directory to the virtual machine (container), eg. '/data/primerscore2/design/:/design/'. This command will put you into the virtual machine command line.

3. Preprocessing  
   `cd bin/database/`  
   `gunzip *gz`  
   `../third-party/bwa index GRCh37.fasta` or `../third-party/bwa index GRCh38.fasta` based on your needs, this will consume some time.
   
To use the 2nd variant, download docker image of PrimerScore2 as follows:

1. Download docker image
   `docker pull zenghp88/primerscore2:native`

2. Prepare files  
   You need to prepare the following files and do some preprocessing.  
   1) a database file to check specificity named 'database.fata'.  
   `bin/third-party/bwa index database.fasta`  
   
   2) maybe a reference genome to get template sequence when inputing a target spot file or a target region file(bed formate), named 'genome.fasta'.
   
   3) maybe a common SNP file to check SNP covered by oligos, named 'common.vcf.gz'.  
      `perl bin/reference_add_snp.pl -ir genome.fasta -iv common.vcf.gz -k genome.fasta`.  
      Note: It will consume a large memory, if your server doesn't have enough memory, you can split the genome file to more files and run them respectively, and cat the results into a whole result file.  

## Source code install
1. Denpendencies  
   PrimerScore2 needs some dependencies which should be attached to file ./path.pm   
   perl 5.22+   
   BWA 0.7.11+  
   Samtools 1.5+  
   Blat v.36+  
   primer3 2.5.0+  

2. Install  
   `git clone git@github.com:zenghp88/primerScore.git`

3. Prepare  
   1) download the denpendencies and modify file ./path.pm with their paths  
   2) prepare files as in docker 2nd variant.

## Example Usage
design generic primers targeting on template sequence:  
perl oligoScore.pl -it demo.fasta --ComeFromRefer -p demo -rdis 120,160,80,200 -opttm 60 -opttmp 70 --Probe -od outdir

design generic primers targeting SNPs:  
perl oligoScore.pl -it demo.txt -p demo -type face-to-face -rpos 10,20,5,45 -rdis 120,160,80,200 -ptype MultiPlex -od outdir

design arms primers:  
perl oligoScore.pl -it demo.txt -ir genome.fasta -is genome_add_snp.fasta -id database.fasta -p demo --Probe -type face-to-face -rpos 0,0,0,0 -rdis 120,160,80,200 -maxlp 40 -od outdir

design Nested primers:  
perl oligoScore.pl -it demo.txt -p demo -type Nested -rpos 10,20,5,45 -rdis=-15,-10,-30,-5 -ptype MultiPlex -od outdir

design Full-coverd primers:  
perl oligoScore.pl -it demo.bed -ir genome.fasta -is genome_add_snp.fasta -id database.fasta -p demo -type face-to-face -rdis 120,160,80,200 -ptype MultiPlex -ctype Full-covered -ds 80 -rf 0.2 -od outdir  

## Parameters description  

```
  -it        <file>   Input target file(SNP file or template fasta file with no non-ATCGatcg), forced
   --ComeFromRefer    Sequences in target file(-it) come from reference file(-ir) when -it is fasta file, optional
  -ir        <file>   Input reference file to extract template sequence of SNP, needed when target file(-it) is SNP file, [/bin/database/GRCh37.fasta]
  -is        <file>   Input reference file containing snps to check SNP of oligos when -it is SNP file, optional
  -id       <files>   Input database files separated by "," to check specificity, [/bin/database/GRCh37.fasta]
  -p         <str>    prefix of output file, forced
  --Probe             Design probe when -type "face-to-face", optional
  --NoFilter          Not filter any oligos
  --Precise           Evalue specificity precisely, but will consume a long time(reach to 100 times longer)
  --Homology_check    Check homologous sequence of template sequence when design for NGS primers, optional

  ### oligo parameters
  -opttm    <int>     optimal tm of primer, [60]
  -opttmp   <int>     optimal tm of probe, [70]
  -minl     <int>     minimum length of primer, [18]
  -maxl     <int>     maximum length of primer, [28]
  -minlp    <int>     minimum length of probe, [18]
  -maxlp    <int>     maximum length of probe, [36]
  -scalel   <str>     candidate oligo length scale, [2]
  -rpos     <str>     position range, distance of p1 to the detected site, (opt_min, opt_max, min, max) separted by ",", must be given when -it is SNP file
  -rdis     <str>     distance range between pair primers, that is product size range when -type is "face-to-face", (opt_min, opt_max, min, max) separted by ",", [120,160,80,200]
  -regions  <str>     interested regions of candidate primers walking on template, format is "start,end,scale,fr,start2,end2,scale2,fr2...", if not given, will caculate automatically, fr:F/R/FR, optional

  ### design parameters
  -type   <str>     primer type, "face-to-face", "back-to-back", "Nested", ["face-to-face"]
  -ptype  <str>     plex type, "SinglePlex" or "MultiPlex", [SinglePlex]
  -ctype  <str>     primer covered type, "Single" or "Full-covered", ["Single"]
     -ds  <int>     average distance between adjacent primers when -ctype "Full-covered", [500]
     -rf  <float>   ratio of distance between adjacent primers can float when -ctype "Full-covered", [0.2]
     -on  <int>     output num when -ctype "Single",[3]
  -pnum   <int>     position num of candidate oligos, [20]

  ### specificity parameters
  -stm    <int>     min tm to be High_tm in specifity, [45]
  -size   <int>     max PCR size, [1000]
  -mpro   <int>     maximum products number to be caculated, to reduce running time. [50]
  -meff   <float>   min efficiency to consider a product, [1e-05]

  ### run parameters
  -para   <int>      parallel num, [10]
  -thrd   <int>      thread in bwa, [3]
  -step   <int>      step, [1]
                     1: homology check
                     2: creat and evalue candidate oligos
                     3: score for probe and primer
                     4: multiplex check
  -od     <dir>      Dir of output file, default ./
  -h                Help
```

## note
Because PrimerScore2 evaluates almost all candidate primers, it takes about 5-20 minutes to design one target.

## Online Webserver  
[PrimerScore2](http://primerscore.gtxlab.com/)

