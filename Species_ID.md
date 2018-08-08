# Getting mitochondrial sequences from the HiSeq

Wanted to check on the species. Using mitochondrial sequences from NCBI, extracted scaffolds from assemblies made for each species using HiSeq reads (de novo assemblies).

## 16S

Sequences from *X. mellotropicalis* (`KY080434.1`)

```
blastn -evalue 1e-1 -query mellotrop_NCBI_16S.fa -db ../canu_assembly/BJE3652_canu_assembly.contigsANDunassembled_blastable -out /work/cauretc/2017_Mellotropicalis/blast_chr7_genes/16S_mellotrop_canu_e1.out -outfmt 6
```
Same command for *Hymenochirus* and *Pipa* using Abyss genome from the mother `BJE3814` and the father `BJE3815`, SOAPchimerical; and the SOAP assembly for the *Pipa* mom (43mer).

## Mitochondrial genome

Using mitochondrial sequences from [NCBI](https://www.ncbi.nlm.nih.gov/): `NC_015617.1` (*Pipa carvalhoi*) since no other *Pipa sp.* complete mitochondrial genomes, `NC_015615.1` (*Hymenochirus boettgeri*), `AY581663.1` (*X. mellotropicalis*, partial mitochondrial sequence).

```
module load blast/2.2.28+

blastn -evalue 1e-1 -query ../../2017_Mellotropicalis/blast_chr7_genes/mellotrop_NCBI_incomplete_mitoch.fa -db /work/cauretc/2017_pipoidea/Assemblies/SOAP_pipa_genome_dad_43mers_blastable -out /work/cauretc/2017_pipoidea/mitochondria/SOAP_pipa_genome_dad_43_mitochmello_e1.out -outfmt 6

gunzip -c Pipa_carvalhoi_mitochondria.fa.gz | blastn -evalue 1e-1 -query - -db /work/cauretc/2017_pipoidea/Assemblies/SOAP_pipa_genome_dad_43mers_blastable -out /work/cauretc/2017_pipoidea/mitochondria/SOAP_pipa_genome_dad_43_mitochpcarvalhoi_e1.out -outfmt 6

gunzip -c Hymenochirus_boettgeri_mitochondria.fa.gz | blastn -evalue 1e-1 -query - -db /work/ben/2016_Hymenochirus/BJE3815/BJE3815_genomeAbyss_blastable -out /work/cauretc/2017_pipoidea/mitochondria/Hymenochirus_boettgeri_mitochondria_hymeno_BJE3815_e1_no_maxtarget.out -outfmt 6

blastn -evalue 1e-1 -query mellotrop_NCBI_incomplete_mitoch.fa -db ../canu_assembly/BJE3652_canu_assembly.contigsANDunassembled_blastable -out /work/cauretc/2017_Mellotropicalis/blast_chr7_genes/mitoch_mellotrop_canu_e1.out -outfmt 6
```
Tried also on SOAP and Allpaths assembly for *X. mellotropicalis*: NO HIT.
```
blastn -evalue 1e-1 -query mellotrop_NCBI_incomplete_mitoch.fa -db /work/cauretc/2017_Mellotropicalis/pseudomolecules/allpaths/final.assembly_blastable -out /work/cauretc/2017_Mellotropicalis/blast_chr7_genes/mitoch_incomplete_mellotrop_allapths_e1.ou

 blastn -evalue 1e-1 -query mellotrop_NCBI_incomplete_mitoch.fa -db /work/cauretc/2017_Mellotropicalis/SOAP_assembly/SOAP_Mellotropicalis_BJE3652_genome_33_memory_scaf_blastable -out /work/cauretc/2017_Mellotropicalis/blast_chr7_genes/mitoch_incomplete_mellotrop_SOAP_e1.out -outfmt 6
```

Extracting the corresponding scaffolds then used them either as input to blast against NCBI database (default `blastn` parameters) or making a quick phylogeny using Geneious.

# Geneious phylogeny

Using geneious alignment to detect direction for *X. mellotropicalis* scaffold and multiple sequences extracted from NCBI. Realigned with MAFFT. Geneious Tree Builder. **VALIDATED**

Mapping to reference with the *X. mellotropicalis* scaffold as reference and a bunch of *Hymenochirus sp.* and *Pipa sp.* sequences. Can't do the re-alignment -- messed up everything. Same tree parameters as *X. mellotropicalis*. **KINDA VALIDATED**
