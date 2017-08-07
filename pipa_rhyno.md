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

/work/cauretc/programs/SOAPdenovo2-src-r240/SOAPdenovo-63mer scaff -g /work/cauretc/2017_pipoidea/Assemblies/SOAP_rhyno_genome_63mers -F -p 10 1 >scaff_rhyno.log 2>scaff_rhyno.err
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
