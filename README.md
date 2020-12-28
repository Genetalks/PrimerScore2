# primerScore: a full-automatic multiplex primer design tool

## Introduction
PrimerScore is a universal primer design tool that accurately and automatically generates primers for multiplex PCR panels, and it also apply to singleplex. It can design several types of primers, including face-to-face, back-to-back and nested primers, and both single primer pair and full-covered pairs for each target template (Fig.1). It adopts comprehensive scoring strategy to identity the most appropriate (i.e., highest scoring) primer pairs even on hostile template so that it can design primers successfully for all targets with no failure in one time.

The PrimerScore workflow includes several steps. First, PrimerScore reads in an input file, if the input is a SNP file, PrimerScore extracts the upstream (U) and downstream (D) template sequences flanking the SNPs. Second, PrimerScore generates all possible candidate primers by "walking" along the template region of interest with a user-defined step size. Third, PrimerScore evaluates various features for each candidate primer, including melting temperature (Tm), GC content, and the Tm and dG values for hairpins and dimers and so on. PrimerScore also performs an initial primer specificity evaluation by mapping to a database. Fourth, PrimerScore scores every feature for single primers using a weighted linear function and adds these scores to generate a final score for each single primer. Then, PrimerScore scores primer pairs and outputs the three highest-scoring primer pairs. Finally, PrimerScore performs three additional checks for all primer pairs output: it re-evaluates primer specificity by calculating the binding efficiency for each possible binding region; it checks for cross-dimerizations among all output primer pairs when designing primers for multiplex panels; and it checks whether the primer binding region includes any common SNPs by searching against dbSNP.

## Online Webserver  
[PrimerScore](http://primerscore.gtxlab.com/)

## Requirements
PrimerScore needs some dependencies which should be attached to file ./path.pm  
perl 5.22+  
BWA  
Samtools  
Blat  
Primer3  
Database: usually a reference genome sequence such as HG19, should be indexed by bwa.  
dbSNP VCF file: optional, if your input SNP file is list of rs IDs.  

## Install
1. git clone git@github.com:zenghp88/primerScore.git  
2. download the denpendencies and modify file ./path.pm with their paths  
3. bwa index the Database sequence  

## Example Usage
design generic primers targeting on template sequence:  
perl primerScore.pl -it template.fasta -p demo -tlen 1700 -rdis 100,150,70,200 -type face-to-face:Region -od outdir

design generic primers targeting SNPs:  
perl primerScore.pl -it demo.txt -p demo --homology_check --dimer_check -mintm 65 -maxtm 75 -minl 20 -maxl 35 -rdis 130,210,70,270 -rpos 8,12,4,17 -type "face-to-face:SNP" -od outdir

design arms primers:  
perl primerScore.pl -it demo.txt -p demo -rdis 100,150,70,200 -rpos 0,0,0,0 -type face-to-face:SNP -od outdir

design Nested primers:  
perl primerScore.pl -it demo.txt -p demo -rdis 10,15,1,30 -rpos 10,20,5,45 -type Nested -od outdir

## note
Because PrimerScore evaluates almost all candidate primers, it takes a long time, such as maybe 10 minutes to design one SNP.

