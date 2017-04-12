# Sequence quality
Check for over-represented sequences, adapters (need all the sequences) 
```
/usr/local/fastqc/fastqc -o /4/caroline/Pipa_parva/Rad_seq/ /4/caroline/1016_S2_L002_R1_001.fastq.gz
```
# STACKS: process_radtags
## Installing locally the new version of stacks 
(need new version of gcc see [here](https://github.com/caro46/Hymenochirus/blob/master/some_commands.Rmd) for how to install gcc)
```
export CC=/home/caroline/programs/gcc-6.2.0/bin/gcc
export CXX=/home/caroline/programs/gcc-6.2.0/bin/g++
export CPP=/home/caroline/programs/gcc-6.2.0/bin/cpp
export LD=/home/caroline/programs/gcc-6.2.0/bin/gcc

./configure --prefix=/home/caroline/programs/stacks-1.45 CC=/home/caroline/programs/gcc-6.2.0/bin/gcc
make
make install
```
## Run on the radseq
Reads seem pretty good (trimmed using `-t 75` but will do trimmomatic after to remove bad quality reads and contaminated sequences)
```
export LD_LIBRARY_PATH=/home/caroline/programs/gcc-6.2.0/lib64
/home/caroline/programs/stacks-1.45/bin/process_radtags -f /4/caroline/1016_S2_L002_R1_001.fastq.gz -b /4/caroline/Pipa_parva/Rad_seq/demultiplex/demutiplex_barcode -o /4/caroline/Pipa_parva/Rad_seq/demultiplex/ -e sbfI -r -c -q --filter_illumina -t 75
```
```
Processing single-end data.
Using Phred+33 encoding for quality scores.
Reads will be truncated to 75bp
Discarding reads marked as 'failed' by Illumina's chastity/purity filters.
Found 1 input file(s).
Searching for single-end, inlined barcodes.
Loaded 96 barcodes (10bp).
Will attempt to recover barcodes with at most 1 mismatches.

Outputing details to log: '/4/caroline/Pipa_parva/Rad_seq/demultiplex/process_radtags.caroline.log'

411680905 total sequences
        0 failed Illumina filtered reads (0.0%)
 52662332 ambiguous barcode drops (12.8%)
   124035 low quality read drops (0.0%)
  4709566 ambiguous RAD-Tag drops (1.1%)
354184972 retained reads (86.0%)
```
Checking the size of the files after demultiplexing, pretty small for the parents...
```
#!/usr/bin/perl

use warnings;
use strict;

# This script will run trimmomatic on SE radtag reads  

my $status;
my $commandline;
my @files;
my $path = "/4/caroline/Pipa_parva/Rad_seq/demultiplex/";
@files = glob($path."*.fq.gz");
my $x;
my @files_no_extension;
my @temp;

foreach(@files){
    @temp=split(".fq.gz",$_);
    push(@files_no_extension,$temp[0]);
}

for($x =0; $x <= $#files_no_extension; $x ++){
        $commandline = "java -jar /usr/local/trimmomatic/trimmomatic-0.36.jar SE -trimlog ".$_."_log.txt ".$files_no_extension[$x].".fq.gz ".$files_no_extension[$x]."_single.fastq.gz ILLUMINACLIP:/home/caroline/programs/Trimmomatic-0.32/adapters/RADseq_repeats_seq_Pipa_Laevis.fa:2:30:10 SLIDINGWINDOW:4:15 MINLEN:36";
        $status = system($commandline);
}
```
