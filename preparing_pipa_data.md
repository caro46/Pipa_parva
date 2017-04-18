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
Reads seem pretty good (trimmed using `-t 75` but will do trimmomatic after to remove bad quality reads and contaminated sequences).

Need that each sample name in demultiplex_barcode is unique (ex. BJE4294_1 ...)
```
export LD_LIBRARY_PATH=/home/caroline/programs/gcc-6.2.0/lib64
/home/caroline/programs/stacks-1.45/bin/process_radtags -f /4/caroline/1016_S2_L002_R1_001.fastq.gz -b /4/caroline/Pipa_parva/Rad_seq/demultiplex/demutiplex_barcode -o /4/caroline/Pipa_parva/Rad_seq/demultiplex/ -e sbfI -r -c -q --filter_illumina -t 75
```
Results of demultiplexing
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
Concatenate the different files for each *Pipa* sample
```perl
#!/usr/bin/perl 

use warnings;
use strict;

my $path_to_data="/4/caroline/Pipa_parva/Rad_seq/demultiplex/";
my @individuals=("BJE4294",
"BJE4295",
"BJE4296",
"BJE4299",
"BJE4300",
"BJE4301",
"BJE4302",
"BJE4303",
"BJE4304",
"BJE4305",
"BJE4306",
"BJE4307",
"BJE4308",
"BJE4309");

my $commandline;
my $status;

##Concatenate files for each individual after demultiplex 
foreach my $individual (@individuals) {
$commandline = "mkdir ".$path_to_data."\/".$individual;
$status = system($commandline);
$commandline = "zcat ".$path_to_data.$individual."_*.fq.gz | gzip \> ".$path_to_data."\/".$individual."\/".$individual.".fq.gz";
$status = system($commandline);

}

```
Checking the size of the files after demultiplexing, pretty small for the parents...
```perl
#!/usr/bin/perl

use warnings;
use strict;

# This script will run trimmomatic on SE radtag reads  

my $status;
my $commandline;
my @files;
my $path = "/4/caroline/Pipa_parva/Rad_seq/demultiplex/";
my @individuals=("BJE4294",
"BJE4295",
"BJE4296",
"BJE4299",
"BJE4300",
"BJE4301",
"BJE4302",
"BJE4303",
"BJE4304",
"BJE4305",
"BJE4306",
"BJE4307",
"BJE4308",
"BJE4309");

foreach my $individual (@individuals){

@files = glob($path.$individual."\/".$individual.".fq.gz");
my $x;
my @files_no_extension;
my @temp;

foreach(@files){
    @temp=split(".fq.gz",$_);
    push(@files_no_extension,$temp[0]);
}

for($x =0; $x <= $#files_no_extension; $x ++){
        $commandline = "java -jar /usr/local/trimmomatic/trimmomatic-0.36.jar SE -trimlog ".$files_no_extension[$x]."_log.txt ".$files_no_extension[$x].".fq.gz ".$files_no_extension[$x]."_single_trimmed.fastq.gz ILLUMINACLIP:/home/caroline/programs/adapters_TruSeq2_3_SE_repeats.fa:2:30:10 SLIDINGWINDOW:4:15 MINLEN:36";
        $status = system($commandline);
}

}
```
