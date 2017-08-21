## 1- De novo assemblies 

### [*Pipa parva*](https://github.com/caro46/Pipa_parva/blob/master/assembly.md)

### *Rhinophrynus*
```
/work/cauretc/programs/SOAPdenovo2-src-r240/SOAPdenovo-63mer all -s /work/cauretc/2017_pipoidea/rhyno.config -K 43 -R -V -p 10 -F -o /work/cauretc/2017_pipoidea/Assemblies/SOAP_rhyno_genome_43mers 1>scaff_rhyno.log 2>scaff_rhyno.err
#Segmentation fault

/work/cauretc/programs/SOAPdenovo2-src-r240/SOAPdenovo-63mer all -s /work/cauretc/2017_pipoidea/rhyno.config -K 43 -R -V -p 10 -o /work/cauretc/2017_pipoidea/Assemblies/SOAP_rhyno_genome_43mers 1>ass_rhyno.log 2>ass_rhyno.err
#Segmentation fault
```
Failed with `Cannot open /work/cauretc/2017_pipoidea/Assemblies/SOAP_rhyno_genome_43mers.peGrads. Now exit to system...` but the file is in the good directory... Since there was issue with Sharcnet yesterday, let's delete everything and start again (same command used). Might need to rerun it later today (Aug 4) or tomorrow because of issues on the server. Ended with `Bus error`. Need to wait the issue on the servor to be solved. Ended again with `Segmentation fault`; tried with a bigger k-mer size (simplify the graph)
```
/work/cauretc/programs/SOAPdenovo2-src-r240/SOAPdenovo-63mer all -s /work/cauretc/2017_pipoidea/rhyno.config -K 63 -R -V -p 10 -o /work/cauretc/2017_pipoidea/Assemblies/SOAP_rhyno_genome_63mers 1>ass_rhyno.log 2>ass_rhyno.err
#Bus error
```
Server failed, cannot access any file (`permission denied`on Aug 7). Need to be rerun later (issues because of unpredicted shutdown, can take a while to fix everything).
Again `Segmentation fault`.
```
/work/cauretc/programs/SOAPdenovo2-src-r240/SOAPdenovo-63mer scaff -g /work/cauretc/2017_pipoidea/Assemblies/SOAP_rhyno_genome_63mers -F -p 10 1 >scaff_rhyno.log 2>scaff_rhyno.err
```
Tried to run the SOAP de novo assembly multiple times since the server is stable (3 times since the server fixation) using different values of kmers to try to use less memory but keep having `Segmentation fault`. 

Trying with Abyss:
```
module unload intel mkl openmpi
module load gcc/4.9.2
module load openmpi/gcc492-std/1.8.7
module load boost/gcc492-openmpi187std/1.59.0

abyss-pe np=8 name=CSL6209 k=64 in='/work/cauretc/2017_pipoidea/2017_Pipa_Rhino_genomes/CSL6209_S22_L003_R1_001_trim_paired_Scythe.cor.fastq.gz /work/cauretc/2017_pipoidea/2017_Pipa_Rhino_genomes/CSL6209_S22_L003_R2_001_trim_paired_Scythe.cor.fastq.gz'
```
```
Building the suffix array...
Building the Burrows-Wheeler transform...
Building the character occurrence table...
sort: write failed: /tmp/sortjgKeIs: No space left on device
error: `CSL6209-3.hist': No such file or directory
make: *** [CSL6209-3.dist] Error 1
make: *** Deleting file `CSL6209-3.dist'
``` 
Suggestions from [here](https://groups.google.com/forum/#!topic/abyss-users/x8Wd0tEnyIw) to change the Environment Variables and redirect the temporary directory to somewhere with more space 
```
TMPDIR=/scratch/cauretc/temp
```
## 2- Mugsy - not used

With all the *de novo* assemblies from Hymenochirus/Pipa/Rhyno and the 2 references genomes *X. tropicalis* and *X. laevis*.
```
/usr/local/mugsy/mugsy -duplications 1 --directory /home/caroline/hymeno/MUGSY --prefix Mugsy_alignment_BJE3814_BJE3815 /4/caroline/tropicalis_genome/Xtropicalis_v9_repeatMasked_HARD_MASK.fa /4/caroline/laevis_genome/Xla.v91_repeatMasked_HARD_MASK.fa /home/caroline/hymeno/MUMMER_analysis/BJE3814-8.fa /home/caroline/hymeno/MUMMER_analysis/BJE3815-8.fa
```

Mugsy is for closely related species which is not really our case.

## 3- [Mauve](http://darlinglab.org/mauve/user-guide/introduction.html) and [mafft](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC135756/)

Mauve (`2.3.1` on info): local multiple alignments, *identifies conserved genomic regions, rearrangements and inversions in conserved regions, and the exact sequence breakpoints of such rearrangements across multiple genomes*, large-scale evolutionary events.

```
/usr/local/mauve/current/mauveAligner [options] <seq1 filename> ... --output=pipoidea.mauve --output-alignment=pipoidea.alignment --scratch-path=<path> 
```
`--scratch-path` for large genomes

`--alignment-output-format=phylip` to be input for MAFFT.

MAFFT (`v7.205` on info) might improved the alignment obtained from Mauve. Default MAFFT output is fasta so good for Gblock.

```
mafft --threats 8
```


Maybe considering MULAN or PECAN or TBA

## 4- [Gblocks](http://molevol.cmima.csic.es/castresana/Gblocks/Gblocks_documentation.html)

Version `0.91b` on info.

To remove bad aligned sequences, too diverged/non informative sequences. 

```
/usr/local/Gblocks [alignment_name]
```
