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
Error for BJE4295 during the trimming (and fastqc failed to produce a file). Run it again. If an error again, need to copy the original data again to sharcnet (maybe error when copy). So when it was run by itself, it went fine; perhaps an issue on sharcnet during the previous run...
Exact command of the re-run
```
java -Xmx1G -jar /work/cauretc/programs/Trimmomatic-0.36/trimmomatic-0.36.jar PE -trimlog /work/cauretc/2017_pipoidea/2017_Pipa_Rhino_genomes/BJE4295_S21_L003_R1_001_log.txt /work/cauretc/2017_pipoidea/2017_Pipa_Rhino_genomes/BJE4295_S21_L003_R1_001.fastq.gz /work/cauretc/2017_pipoidea/2017_Pipa_Rhino_genomes/BJE4295_S21_L003_R2_001.fastq.gz /work/cauretc/2017_pipoidea/2017_Pipa_Rhino_genomes/BJE4295_S21_L003_R1_001_trim_paired.fastq.gz /work/cauretc/2017_pipoidea/2017_Pipa_Rhino_genomes/BJE4295_S21_L003_R1_001_trim_single.fastq.gz /work/cauretc/2017_pipoidea/2017_Pipa_Rhino_genomes/BJE4295_S21_L003_R2_001_trim_paired.fastq.gz /work/cauretc/2017_pipoidea/2017_Pipa_Rhino_genomes/BJE4295_S21_L003_R2_001_trim_single.fastq.gz ILLUMINACLIP:/work/cauretc/programs/Trimmomatic-0.36/adapters/Pipoidea_TruSeqPE_fastqc_adapters.fa:2:30:10 SLIDINGWINDOW:4:15 MINLEN:36
```
FastQC
```
/work/cauretc/programs/FastQC/fastqc -o /work/cauretc/2017_pipoidea/fastqc_quality /work/cauretc/2017_pipoidea/2017_Pipa_Rhino_genomes/*_trim_paired.fastq.gz
```
Some issues with kmers. Look if we correctly trimmed. From [Usadellab](http://www.usadellab.org/cms/?page=trimmomatic)

*Naming of the sequences indicates how they should be used. For 'Palindrome' clipping, the sequence names should both start with 'Prefix', and end in '/1' for the forward adapter and '/2' for the reverse adapter. All other sequences are checked using 'simple' mode. Sequences with names ending in '/1' or '/2' will be checked only against the forward or reverse read. Sequences not ending in '/1' or '/2' will be checked against both the forward and reverse read. If you want to check for the reverse-complement of a specific sequence, you need to specifically include the reverse-complemented form of the sequence as well, with another name.*

So we did not specify any /1 or /2 so all the adapters should be tested as forward and reverse in a simple mode. However, the `Palindrome` mode seems a little bit weird for me. A better explanation of the `palindrome` mode than in the manual can be found [here](http://seqanswers.com/forums/archive/index.php/t-11186.html). From *alisrpp*:
```

About the palindrome clipping, some days ago i wrote one of the creators of Trimmomatic asking about an alternative explanation to the one in the web site (i couldn't understand it either).
Here is the answer, for me was useful:

Simple clipping is just finding a contaminant sequence somewhere within a read. Conceptually, you get contaminant and read, and you slide them across each other, until you get a perfect or close enough match. So, with R being read bases, and C being contaminant, you check

1)
RRRRRRRRRRR
CCCC

2)
RRRRRRRRRRR
CCCC ->

etc.

Palindrome clipping is a bit more complex - and related to actual palindromes only in a twisted mind like mine. In this case, you 'ligate' the presumed adapter sequence to the start of each read in a pair, and try sliding them over each other.

So with F being bases from the forward read, R being bases from the reverse read, and A being either adapter (technically the two adapters are different, but lets ignore that for now).

AAAAAAFFFFFFF ->
<- RRRRRRRAAAAAA

In this case, the aligning region is much longer, since it consists of the entire read length plus part of the adapter. This gives a very high confidence that an apparent 'read-though' is a true-positive.
```
## Scythe
```
scythe -a /work/cauretc/programs/Trimmomatic-0.36/adapters/Pipoidea_TruSeqPE_fastqc_adapters_Wilson.fa -p 0.1 -o /work/cauretc/2017_pipoidea/2017_Pipa_Rhino_genomes/BJE4294_S20_L003_R1_001_trim_paired_Scythe.fastq /work/cauretc/2017_pipoidea/2017_Pipa_Rhino_genomes/BJE4294_S20_L003_R1_001_trim_paired.fastq.gz

#prior: 0.100

#Adapter Trimming Complete
#contaminated: 1914742, uncontaminated: 141064880, total: 142979622
#contamination rate: 0.013392

./scythe -a /work/cauretc/programs/Trimmomatic-0.36/adapters/Pipoidea_TruSeqPE_fastqc_adapters_Wilson.fa -p 0.1 -o /work/cauretc/2017_pipoidea/2017_Pipa_Rhino_genomes/BJE4294_S20_L003_R2_001_trim_paired_Scythe.fastq /work/cauretc/2017_pipoidea/2017_Pipa_Rhino_genomes/BJE4294_S20_L003_R2_001_trim_paired.fastq.gz
#prior: 0.100

#Adapter Trimming Complete
#contaminated: 3226899, uncontaminated: 139752723, total: 142979622
#contamination rate: 0.022569
```
[AdapterRemoval v2](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4751634/), [github page](https://github.com/MikkelSchubert/adapterremoval/), seems to be a pretty good program to try if `scythe` does not work as well as expected.
## Jellyfish/quake
For some reason cannot install on sharnet whereas no issue on info...
```
tar -zxvf jellyfish-2.2.4\(1\).tar.gz 
cd jellyfish-2.2.4
./configure --prefix=/home/caroline/programs/jellyfish-2.2.4
make
make install
```
`export R_LIBS=$HOME/Rlibs:$R_LIBS` used to load the `VGAM` library (version `1.0-3` downloaded from [here](https://www.stat.auckland.ac.nz/~yee/VGAM/download.shtml)).
```
zcat /4/caroline/2017_Pipoidea_Hiseq/after_trimmomatic/BJE4294_*_trim_paired.fastq.gz | /home/caroline/programs/jellyfish-2.2.4/bin/jellyfish count /dev/fd/0 -m 19 -s 100M -t 16 -C -o /4/caroline/2017_Pipoidea_Hiseq/jellyfish/BJE4294_jelly_count_19mers

#When BJE4295 trimmo done
zcat /4/caroline/2017_Pipoidea_Hiseq/after_trimmomatic/BJE4295_*_trim_paired.fastq.gz | /home/caroline/programs/jellyfish-2.2.4/bin/jellyfish count /dev/fd/0 -m 19 -s 100M -t 16 -C -o /4/caroline/2017_Pipoidea_Hiseq/jellyfish/BJE4295_jelly_count_19mers

/home/caroline/programs/jellyfish-2.2.4/bin/jellyfish dump -c -t /4/caroline/2017_Pipoidea_Hiseq/jellyfish/BJE4294_jelly_count_19mers -o /4/caroline/2017_Pipoidea_Hiseq/jellyfish/BJE4294_jelly_dump_19mers

export R_LIBS=$HOME/Rlibs:$R_LIBS
/usr/local/quake/bin/cov_model.py --int /4/caroline/2017_Pipoidea_Hiseq/jellyfish/BJE4294_jelly_dump_19mers

#When BJE4295 trimmo done
/home/caroline/programs/jellyfish-2.2.4/bin/jellyfish dump -c -t /4/caroline/2017_Pipoidea_Hiseq/jellyfish/BJE4295_jelly_count_19mers -o /4/caroline/2017_Pipoidea_Hiseq/jellyfish/BJE4295_jelly_dump_19mers
/usr/local/quake/bin/cov_model.py --int /4/caroline/2017_Pipoidea_Hiseq/jellyfish/BJE4295_jelly_dump_19mers
```
Cutoff obtained: BJE4294:
```
/usr/local/quake/bin/correct -f /4/caroline/2017_Pipoidea_Hiseq/quake/filenames_quake_pipa_female.txt -z -k 19 -c 2 -m /home/evanslab/tetra_project/jellyfish_results/19mers/BJE4294_jelly_dump_19mers -p 4

zcat /4/caroline/2017_Pipoidea_Hiseq/after_trimmomatic/CSL6209_*_trim_paired.fastq.gz | /home/caroline/programs/jellyfish-2.2.4/bin/jellyfish count /dev/fd/0 -m 19 -s 100M -t 16 -C -o /4/caroline/2017_Pipoidea_Hiseq/jellyfish/CSL6209_jelly_count_19mers

/home/caroline/programs/jellyfish-2.2.4/bin/jellyfish dump -c -t /4/caroline/2017_Pipoidea_Hiseq/jellyfish/CSL6209_jelly_count_19mers -o /4/caroline/2017_Pipoidea_Hiseq/jellyfish/CSL6209_jelly_dump_19mers
/usr/local/quake/bin/cov_model.py --int /4/caroline/2017_Pipoidea_Hiseq/jellyfish/CSL6209_jelly_dump_19mers

```
For some reason, `cov_model.py` does not work whereas it was fine for *Hymenochirus*.
```R
histoM=read.table("kmers_BJE4294.hist",sep="\t")
histoF=read.table("kmers_BJE4295.hist",sep="\t")
histoR=read.table("kmers_CSL6209.hist",sep="\t")
pdf('19mers_distribution_quake_pipa_Mom_Dad_Rhyno.pdf')
#pdf('19mers_distribution_quake_pipa_Mom_Dad_Rhyno_5_200.pdf')
#pdf('19mers_distribution_quake_pipa_Mom_Dad_Rhyno_4_200.pdf')
#pdf('19mers_distribution_quake_pipa_Mom_Dad_Rhyno_3_200.pdf')
#pdf('19mers_distribution_quake_pipa_Mom_Dad_Rhyno_2_200.pdf')

plot(histoF,type="l",col="pink",ylab="count",xlab="coverage",main="19-mers distribution (Pipa Female & Male & Rhyno)")
lines (histoM, col="blue")
lines (histoR, col="green")
legend('topright',c("Pipa Female BJE3814","Pipa Male BJE3815", "Rhyno"),lty=c(1,1),lwd=c(2.5,2.5),col=c("pink","blue","green"))
dev.off()

#plot(histoF[5:200,],type="l",col="pink",ylab="count",xlab="coverage",main="19-mers distribution (Hymenochirus Female & Male)")
#lines (histoM[5:200,], col="blue")
#lines (histoM[5:200,], col="green")
#legend('topright',c("Female BJE3814","Male BJE3815"),lty=c(1,1),lwd=c(2.5,2.5),col=c("pink","blue"))
#dev.off()
```
We will need to set different values for the cutoff for Pipa and Rhyno
