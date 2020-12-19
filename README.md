# primerScore: a full-automatic multiplex primer design tool

## Introduction
PrimerScore is a universal primer design tool that accurately and automatically generates primers for multiplex PCR panels, and it also apply to singleplex. It can design several types of primers, including face-to-face, back-to-back and nested primers, and both single primer pair and full-covered pairs for each target template (Fig.1). It adopts comprehensive scoring strategy to identity the most appropriate (i.e., highest scoring) primer pairs even on hostile template so that it can design primers successfully for all targets with no failure in one time.

The PrimerScore workflow includes several steps. First, PrimerScore reads in an input file, if the input is a SNP file, PrimerScore extracts the upstream (U) and downstream (D) template sequences flanking the SNPs. Second, PrimerScore generates all possible candidate primers by "walking" along the template region of interest with a user-defined step size. Third, PrimerScore evaluates various features for each candidate primer, including melting temperature (Tm), GC content, and the Tm and dG values for hairpins and dimers and so on. PrimerScore also performs an initial primer specificity evaluation by mapping to a database. Fourth, PrimerScore scores every feature for single primers using a weighted linear function and adds these scores to generate a final score for each single primer. Then, PrimerScore scores primer pairs and outputs the three highest-scoring primer pairs. Finally, PrimerScore performs three additional checks for all primer pairs output: it re-evaluates primer specificity by calculating the binding efficiency for each possible binding region; it checks for cross-dimerizations among all output primer pairs when designing primers for multiplex panels; and it checks whether the primer binding region includes any common SNPs by searching against dbSNP.

## Requirements
PrimerScore needs some dependencies which should be attached to file ./path.pm
perl
BWA  
Samtools
Blat
Primer3
Database: usually a reference genome sequence such as HG19, should be indexed by bwa.
dbSNP VCF file: optional, if your input SNP file is list of rs IDs.

## Install
git clone xxx.git
modify file ./path.pm with your path of denpendencies
bwa index database sequence
perl primerScore.pl

## Usage
