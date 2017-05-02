# RADseq
## Sequence quality
Check for over-represented sequences, adapters (need all the sequences) 
```
/usr/local/fastqc/fastqc -o /4/caroline/Pipa_parva/Rad_seq/ /4/caroline/1016_S2_L002_R1_001.fastq.gz
```
## STACKS: process_radtags
### Installing locally the new version of stacks 
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
### Run on the radseq
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
Checking again quality after trimming
```
/usr/local/fastqc/fastqc -o /4/caroline/Pipa_parva/Rad_seq/demultiplex/fastqc/ /4/caroline/Pipa_parva/Rad_seq/demultiplex/BJE*/*_single_trimmed.fastq.gz
```
Mapping with BWA
```perl
#!/usr/bin/perl 

use warnings;
use strict;

# This script will execute alignment functions for single end reads
# using bwa mem for reference contigs

my $path_to_data="/4/caroline/Pipa_parva/Rad_seq/demultiplex/";
my $path_to_genome="/4/caroline/tropicalis_genome/";
my $path_to_output="/4/caroline/Pipa_parva/Rad_seq/bwa/";
my $genome="Xtropicalis_v9_repeatMasked_HARD_MASK.fa";
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

##Index                                                                                                                      
#$commandline = "bwa index -a bwtsw ".$path_to_genome."\/".$genome;                                 
#$status = system($commandline);                                                                                             
##bwa mem 
foreach my $individual (@individuals) {

#print $path_to_data.$individual."_single_trimmed.fastq.gz";

$commandline = "bwa mem -t 5 -R \'\@RG\\tID:".$individual."\\tSM:".$individual."\\tLB:library1\\tPL:illumina\' ".$path_to_genome.$genome." ".$path_to_data.$individual."\/".$individual."_single_trimmed.fastq.gz \> ".$path_to_output.$individual."_trop.sam";
$status = system($commandline);

##Samtools : sam to bam                                                                                                      
$commandline="samtools view -bt ".$path_to_genome.$genome." -o ".$path_to_output.$individual."_trop.bam ".$path_to_output.$individual."_trop.sam";
$status = system($commandline);

## Samtools : bam to _sorted.bam
$commandline="samtools sort ".$path_to_output.$individual."_trop.bam -o".$path_to_output.$individual."_trop_sorted.bam";
$status = system($commandline);
$commandline= "samtools index ".$path_to_output.$individual."_trop_sorted.bam";
$status = system($commandline);

##Delete unecessary files
print "Done with individual ",$individual,"\n";
$commandline= "rm -f ".$path_to_output.$individual."_trop.sam ".$path_to_output.$individual."_trop.bam";
$status = system($commandline);

}

```
Genotype calls with samtools
```perl
#!/usr/bin/perl

use warnings;
use strict;

#Using samtools to call the genotypes
##RADseq previopusly aligned with bwa (script_bwa_GBS.pl)

my $path_to_genome="/4/caroline/tropicalis_genome/";
my $path_to_bwa_res="/4/caroline/Pipa_parva/Rad_seq/bwa/";
my $path_to_output = "/4/caroline/Pipa_parva/Rad_seq/samtools_genotypes/";
my $genome="Xtropicalis_v9_repeatMasked_HARD_MASK.fa";
my $path_to_vcftab = "/usr/local/vcftools/src/perl/";
my $path_to_tabix = "/usr/local/tabix/";

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

#Variant calling
$commandline = "samtools mpileup -d8000 -uf ".$path_to_genome.$genome." ";

foreach my $individual (@individuals) {
$commandline = $commandline.$path_to_bwa_res.$individual."_trop_sorted.bam ";
}

$commandline = $commandline."| bcftools call -mv -Oz \> ".$path_to_output."Pipa_trop_var.vcf.gz";
$status = system($commandline);

#Variant filtering
#$commandline = "bcftools filter -O z ".$path_to_output."Pipa_trop_var.vcf.gz -s LowQual -e '%QUAL<5 || DP>1' \> ".$path_to_output."Pipa_trop_var.flt.vcf.gz";
$commandline = "bcftools filter -O z ".$path_to_output."Pipa_trop_var.vcf.gz -s LowQual -e 'DP>1' \> ".$path_to_output."Pipa_trop_var.flt.vcf.gz";
$status = system($commandline);

#vcf to tab
$commandline = "tabix -f -p vcf ".$path_to_output."Pipa_trop_var.flt.vcf.gz";
$status = system($commandline);
$commandline = "zcat ".$path_to_output."Pipa_trop_var.flt.vcf.gz | ".$path_to_vcftab."vcf-to-tab \> ".$path_to_output."Pipa_trop_var.tab";
$status = system($commandline);

```
```
#CHROM  POS     REF     BJE4294 BJE4295 BJE4296 BJE4299 BJE4300 BJE4301 BJE4302 BJE4303 BJE4304 BJE4305 BJE4306 BJE4307 BJE4308 BJE4309
Only_daughters_ZW scaffold_51   996407  A       G/G     ./.     ./.     ./.     ./.     G/G     G/G     ./.     ./.     ./.     G/G     G/G
     ./.     G/G
Only_daughters_ZW scaffold_51   996419  G       A/A     ./.     ./.     ./.     ./.     A/A     A/A     ./.     ./.     ./.     A/A     A/A
     ./.     A/A
Only_daughters_ZW scaffold_51   996425  G       T/T     ./.     ./.     ./.     ./.     T/T     T/T     ./.     ./.     ./.     T/T     T/T
     ./.     T/T
```
*X. tropicalis* scaffold_51 996407..996425 = Chromosome 1 *X. laevis* (in intron of npffr2?)

# HiSeq
## Sequence quality 
```
/work/cauretc/programs/FastQC/fastqc -o /work/cauretc/2017_pipoidea/fastqc_quality /work/cauretc/2017_pipoidea/2017_Pipa_Rhino_genomes/*.fastq.gz

```
```perl
#!/usr/bin/perl                                                                                                                                                                                                                                   
use warnings;
use strict;

# This script will read in the *fastq.gz file names in a directory, and                                                                                                                                                                           
# run trimmomatic on each one.                                                                                                                                                                                                                    

my $trimmomatic_path = "/work/cauretc/programs/Trimmomatic-0.36/";
my $data_path = "/work/cauretc/2017_pipoidea/2017_Pipa_Rhino_genomes/";
my $status;
my @files;
my $commandline;
my @temp;
my @pairsR1;
my @pairsR2;
my @unique;



my @filesR1 = glob($data_path."*R1*fastq.gz");
my @filesR2 = glob($data_path."*R2*fastq.gz");

foreach(@filesR1){
    @temp=split(".fastq.gz",$_);
    push(@pairsR1,$temp[0]);
}

foreach(@filesR2){
    @temp=split(".fastq.gz",$_);
    push(@pairsR2,$temp[0]);
}

# make sure the names are in the same order

@filesR1 = sort @filesR1;
@filesR2 = sort @filesR2;

my $x;
my @replace;
my $on_off_switch=1; # this is a switch to tell trimmomatic to work (1 = on, 0 = off)

if($#filesR1 ne $#filesR2){
    print "There is a different number of forward and reverse reads\n";
}
else{
    for($x =0; $x <= $#pairsR1; $x ++){
        ($replace[$x] = $pairsR1[$x]) =~ s/R1/R2/;
        if($replace[$x] ne $pairsR2[$x]){
            print "Problem with filenames\n";
            $on_off_switch = 0;
        }
    }
}

#       if (the name are the same (compare $pairsR1[$x] and $pairsR2[$x]){           

if($on_off_switch == 1){
    for($x =0; $x <= $#pairsR1; $x ++){                                                              
        #print $pairsR1[$x]," ",$pairsR2[$x],"\n";
        $commandline = "java -Xmx1G -jar ".$trimmomatic_path."trimmomatic-0.36.jar PE -trimlog ";
        $commandline = $commandline.$pairsR1[$x]."_log.txt ".$pairsR1[$x].".fastq.gz ".$pairsR2[$x].".fastq.gz ".$pairsR1[$x]."_trim_paired.fastq.gz ".$pairsR1[$x]."_trim_single.fastq.gz ".$pairsR2[$x]."_trim_paired.fastq.gz ".$pairsR2[$x]."_trim_single.fastq.gz ";
        $commandline = $commandline."ILLUMINACLIP:".$trimmomatic_path."adapters/Pipoidea_TruSeqPE_fastqc_adapters.fa:2:30:10 SLIDINGWINDOW:4:15 MINLEN:36";
        print $commandline,"\n";
        $status = system($commandline);
    }
}
```
```
/work/cauretc/programs/FastQC/fastqc -o /work/cauretc/2017_pipoidea/fastqc_quality /work/cauretc/2017_pipoidea/2017_Pipa_Rhino_genomes/*_trim_paired.fastq.gz
```
