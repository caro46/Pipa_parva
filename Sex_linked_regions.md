vcf to tab
```
/usr/local/vcftools/src/perl/vcf-to-tab </4/caroline/Pipa_parva/Rad_seq/GATK/SOAP_Chimerical_recalibrated_allsites_Phred20_HardFilteringSNP_filtered.vcf>/4/caroline/Pipa_parva/Rad_seq/GATK/SOAP_Chimerical_recalibrated_allsites_Phred20_HardFilteringSNP_filtered.tab

/usr/local/vcftools/src/perl/vcf-to-tab </4/caroline/Pipa_parva/Rad_seq/GATK/SOAP_Chimerical_recalibrated_allsites_Phred20_HardFilteringSNP.vcf>/4/caroline/Pipa_parva/Rad_seq/GATK/SOAP_Chimerical_recalibrated_allsites_Phred20_HardFilteringSNP.tab
```
On info `BLAST 2.6.0+`
```
makeblastdb -in /4/caroline/tropicalis_genome/Xtropicalis_v9_repeatMasked_HARD_MASK.fa -dbtype nucl -title UCSC_xenTro9 -out /4/caroline/tropicalis_genome/Xtropicalis_v9_repeatMasked_HARD_MASK_blastable
```
```
blastn -evalue 1e-20 -query /4/caroline/Pipa_parva/Rad_seq/samtools_genotypes/Sex_linked/SOAP_chim_assembly/Pipa_putative_sex_linked_polym1ratio0_HF.fa -db /4/caroline/tropicalis_genome/Xtropicalis_v9_repeatMasked_HARD_MASK_blastable -out /4/caroline/Pipa_parva/Rad_seq/samtools_genotypes/Sex_linked/SOAP_chim_assembly/blast_results/Pipa_chimerical_putative_sex_linked_Phred20_genocaller_xenTro9_hard_mask_e20 -outfmt 6 -max_target_seqs 1
```
```
/usr/local/RepeatMasker/RepeatMasker -dir /4/caroline/Pipa_parva/Rad_seq/samtools_genotypes/Sex_linked/SOAP_chim_assembly/ -species "xenopus genus" -pa 4 -a /4/caroline/Pipa_parva/Rad_seq/samtools_genotypes/Sex_linked/SOAP_chim_assembly/Pipa_putative_sex_linked_polym1ratio0_HF.fa
```
```
blastn -evalue 1e-20 -query /4/caroline/Pipa_parva/Rad_seq/samtools_genotypes/Sex_linked/SOAP_chim_assembly/Pipa_putative_sex_linked_polym1ratio0_HF.fa.masked -db /4/caroline/tropicalis_genome/Xtropicalis_v9_repeatMasked_HARD_MASK_blastable -out /4/caroline/Pipa_parva/Rad_seq/samtools_genotypes/Sex_linked/SOAP_chim_assembly/blast_results/Pipa_chimerical_putative_sex_linked_Phred20_genocaller_masked_xenTro9_hard_mask_e20 -outfmt 6 -max_target_seqs 1
```
### Without filtering
```
/usr/local/RepeatMasker/RepeatMasker -dir /4/caroline/Pipa_parva/Rad_seq/samtools_genotypes/Sex_linked/SOAP_chim_assembly/no_filtered/ -species "xenopus genus" -pa 4 -a /4/caroline/Pipa_parva/Rad_seq/samtools_genotypes/Sex_linked/SOAP_chim_assembly/no_filtered/Pipa_putative_sex_linked_polym1ratio0_HF.fa
```
```
blastn -evalue 1e-20 -query /4/caroline/Pipa_parva/Rad_seq/samtools_genotypes/Sex_linked/SOAP_chim_assembly/no_filtered/Pipa_putative_sex_linked_polym1ratio0_HF.fa.masked -db /4/caroline/tropicalis_genome/Xtropicalis_v9_repeatMasked_HARD_MASK_blastable -out /4/caroline/Pipa_parva/Rad_seq/samtools_genotypes/Sex_linked/SOAP_chim_assembly/blast_results/no_filtered/Pipa_chimerical_putative_sex_linked_Phred20_genocaller_masked_xenTro9_hard_mask_e20 -outfmt 6 -max_target_seqs 1

blastn -evalue 1e-20 -query /4/caroline/Pipa_parva/Rad_seq/samtools_genotypes/Sex_linked/SOAP_chim_assembly/no_filtered/Pipa_putative_sex_linked_polym1ratio0_HF.fa.masked -db /4/caroline/tropicalis_genome/Xtropicalis_v9_repeatMasked_HARD_MASK_blastable -out /4/caroline/Pipa_parva/Rad_seq/samtools_genotypes/Sex_linked/SOAP_chim_assembly/blast_results/no_filtered/Pipa_chimerical_putative_sex_linked_Phred20_genocaller_masked_xenTro9_hard_mask_e20 -outfmt 6 -max_target_seqs 1
blastn -evalue 1e-5 -query /4/caroline/Pipa_parva/Rad_seq/samtools_genotypes/Sex_linked/SOAP_chim_assembly/no_filtered/Pipa_putative_sex_linked_polym1ratio0_HF.fa.masked -db /4/caroline/tropicalis_genome/Xtropicalis_v9_repeatMasked_HARD_MASK_blastable -out /4/caroline/Pipa_parva/Rad_seq/samtools_genotypes/Sex_linked/SOAP_chim_assembly/blast_results/no_filtered/Pipa_chimerical_putative_sex_linked_Phred20_genocaller_masked_xenTro9_hard_mask_e5 -outfmt 6 -max_target_seqs 1
blastn -evalue 1e-1 -query /4/caroline/Pipa_parva/Rad_seq/samtools_genotypes/Sex_linked/SOAP_chim_assembly/no_filtered/Pipa_putative_sex_linked_polym1ratio0_HF.fa.masked -db /4/caroline/tropicalis_genome/Xtropicalis_v9_repeatMasked_HARD_MASK_blastable -out /4/caroline/Pipa_parva/Rad_seq/samtools_genotypes/Sex_linked/SOAP_chim_assembly/blast_results/no_filtered/Pipa_chimerical_putative_sex_linked_Phred20_genocaller_masked_xenTro9_hard_mask_e1 -outfmt 6 -max_target_seqs 1
blastn -evalue 1e-20 -query /4/caroline/Pipa_parva/Rad_seq/samtools_genotypes/Sex_linked/SOAP_chim_assembly/no_filtered/Pipa_putative_sex_linked_polym1ratio0_HF_HARDMASKED.fa -db /4/caroline/tropicalis_genome/Xtropicalis_v9_repeatMasked_HARD_MASK_blastable -out /4/caroline/Pipa_parva/Rad_seq/samtools_genotypes/Sex_linked/SOAP_chim_assembly/blast_results/no_filtered/Pipa_chimerical_putative_sex_linked_Phred20_genocaller_hard_masked_xenTro9_hard_mask_e20 -outfmt 6 -max_target_seqs 1
blastn -evalue 1e-1 -query /4/caroline/Pipa_parva/Rad_seq/samtools_genotypes/Sex_linked/SOAP_chim_assembly/no_filtered/Pipa_putative_sex_linked_polym1ratio0_HF_HARDMASKED.fa -db /4/caroline/tropicalis_genome/Xtropicalis_v9_repeatMasked_HARD_MASK_blastable -out /4/caroline/Pipa_parva/Rad_seq/samtools_genotypes/Sex_linked/SOAP_chim_assembly/blast_results/no_filtered/Pipa_chimerical_putative_sex_linked_Phred20_genocaller_hard_masked_xenTro9_hard_mask_e1 -outfmt 6 -max_target_seqs 1
blastn -evalue 1e-1 -query /4/caroline/Pipa_parva/Rad_seq/samtools_genotypes/Sex_linked/SOAP_chim_assembly/no_filtered/Pipa_putative_sex_linked_polym1ratio0_HF_prop07_HARDMASKED.fa -db /4/caroline/tropicalis_genome/Xtropicalis_v9_repeatMasked_HARD_MASK_blastable -out /4/caroline/Pipa_parva/Rad_seq/samtools_genotypes/Sex_linked/SOAP_chim_assembly/blast_results/no_filtered/Pipa_chimerical_putative_sex_linked_Phred20_genocaller_prop07_hard_masked_xenTro9_hard_mask_e1 -outfmt 6 -max_target_seqs 1
```
Soft masked - hard masked
```
tr '[:lower:]' 'N' </4/caroline/Pipa_parva/Rad_seq/samtools_genotypes/Sex_linked/SOAP_chim_assembly/no_filtered/Pipa_putative_sex_linked_polym1ratio0_HF_prop07.fa.masked >/4/caroline/Pipa_parva/Rad_seq/samtools_genotypes/Sex_linked/SOAP_chim_assembly/no_filtered/Pipa_putative_sex_linked_polym1ratio0_HF_prop07_HARDMASKED.fa
sed 's/>NNNNNNNN/>scaffold_/' /4/caroline/Pipa_parva/Rad_seq/samtools_genotypes/Sex_linked/SOAP_chim_assembly/no_filtered/Pipa_putative_sex_linked_polym1ratio0_HF_prop07_HARDMASKED.fa >/4/caroline/Pipa_parva/Rad_seq/samtools_genotypes/Sex_linked/SOAP_chim_assembly/no_filtered/Pipa_putative_sex_linked_polym1ratio0_HF_prop07_HARDMASKED1.fa
rm -f /4/caroline/Pipa_parva/Rad_seq/samtools_genotypes/Sex_linked/SOAP_chim_assembly/no_filtered/Pipa_putative_sex_linked_polym1ratio0_HF_prop07_HARDMASKED.fa
mv /4/caroline/Pipa_parva/Rad_seq/samtools_genotypes/Sex_linked/SOAP_chim_assembly/no_filtered/Pipa_putative_sex_linked_polym1ratio0_HF_prop07_HARDMASKED1.fa /4/caroline/Pipa_parva/Rad_seq/samtools_genotypes/Sex_linked/SOAP_chim_assembly/no_filtered/Pipa_putative_sex_linked_polym1ratio0_HF_prop07_HARDMASKED.fa
```

### Genes
- Scaffolds 726 -> NCOA2 (transcriptional coactivator for stroid receptor), 12050 -> chd7. A CHD gene is used to determined the sex of some birds ([see](http://www.academia.edu/10024376/Avian_Sex_Determination_Based_on_Chromo_Helicase_DNA-binding_CHD_Genes_Using_Polymerase_Chain_Reaction_PCR)).

- On Chr.06: catenin beta 1 (ctnnb1) involve in wnt signaling pathway, SOX4, wnt9a ([ovarian development](http://www.sciencedirect.com/science/article/pii/S0303720706005843?via%3Dihub)).

- Newly identified as W specific in *X.laevis*: scanw, on chr.6 there is also a zinc finger protein like `znf16` 
```
module load blast/2.2.28+
blastn -evalue 1e-1 -query /work/cauretc/2017_pipoidea/blast_results/znf16.S_xenbase_Xlaev.fa -db /work/cauretc/2017_pipoidea/Assemblies/SOAP_pipa_genome_mom_43mers_blastable -out /work/cauretc/2017_pipoidea/blast_results/Pipa_mom_Xlaev_znf16 -outfmt 6
blastn -evalue 1e-1 -query /work/cauretc/2017_pipoidea/blast_results/znf16.S_xenbase_Xlaev.fa -db /work/cauretc/2017_pipoidea/Assemblies/SOAP_pipa_genome_mom_43mers_blastable -out /work/cauretc/2017_pipoidea/blast_results/Pipa_mom_Xlaev_znf16_maxtarget1 -outfmt 6 -max_target_seqs 1
blastn -evalue 1e-1 -query /work/cauretc/2017_pipoidea/blast_results/znf16.S_xenbase_Xlaev.fa -db /work/cauretc/2017_pipoidea/Assemblies/SOAP_pipa_genome_mom_43mers_blastable -out /work/cauretc/2017_pipoidea/blast_results/Pipa_mom_Xlaev_znf16_maxtarget1 -outfmt 6 -max_target_seqs 1
blastn -evalue 1e-1 -query /work/cauretc/2017_pipoidea/blast_results/znf16.S_xenbase_Xlaev.fa -db /work/cauretc/2017_pipoidea/Assemblies/SOAP_pipa_genome_dad_43mers_blastable -out /work/cauretc/2017_pipoidea/blast_results/Pipa_dad_Xlaev_znf16_maxtarget1 -outfmt 6 -max_target_seqs 1
```
# HiSeq genotypes - chr.06
- Calling genotypes from HiSeq data from dad and mom (using supercontigs of chimerical assembly as ref)

- Identifying if SNP in mom different from dad in wnt9a

```
blastn -evalue 1e-20 -query /4/caroline/Pipa_parva/blast_genes/wnt9a_xenbase_Xtrop.fa -db /4/caroline/2017_Pipoidea_Hiseq/Assemblies/SOAP_pipa_genome_chimerical_43mers_blastable -out /4/caroline/Pipa_parva/blast_genes/Pipa_chimerical_Xtrop_wnt9a_e20 -outfmt 6 -max_target_seqs 1
blastn -evalue 1e-20 -query /4/caroline/Pipa_parva/blast_genes/wnt9a_xenbase_Xtrop.fa -db /4/caroline/2017_Pipoidea_Hiseq/Assemblies/SOAP_pipa_genome_chimerical_43mers_blastable -out /4/caroline/Pipa_parva/blast_genes/Pipa_chimerical_Xtrop_wnt9a_e20_pairwise -outfmt 0 -max_target_seqs 1
```
Mapped against `scaffold173038` of chimerical assembly. But cannot find the exon 1, even with:
```
blastn -evalue 1e-1 -query /4/caroline/Pipa_parva/blast_genes/wnt9a_xenbase_Xtrop_e1_only.fa -db /4/caroline/2017_Pipoidea_Hiseq/Assemblies/SOAP_pipa_genome_chimerical_43mers_blastable -out /4/caroline/Pipa_parva/blast_genes/Pipa_chimerical_Xtrop_wnt9a_exon1_only_e1 -outfmt 6 -max_target_seqs 1

#Trying on sex-specific (mom and dad) assemblies on Iqaluk. Still nothing for any of the assembly.
module load blast/2.2.28+
gunzip -c /work/cauretc/2017_pipoidea/Assemblies/SOAP_pipa_genome_mom_43mers.scafSeq.gz | makeblastdb -in - -dbtype nucl -title UCSC_xenTro9 -out /work/cauretc/2017_pipoidea/Assemblies/SOAP_pipa_genome_mom_43mers_blastable
blastn -evalue 1e-1 -query /work/cauretc/2017_pipoidea/blast_results/wnt9a_xenbase_Xtrop_e1_only.fa -db /work/cauretc/2017_pipoidea/Assemblies/SOAP_pipa_genome_mom_43mers_blastable -out /work/cauretc/2017_pipoidea/blast_results/Pipa_mom_Xtrop_wnt9a_exon1_only_e1 -outfmt 6
```
```
grep "scaffold173038" /4/caroline/2017_Pipoidea_Hiseq/Assemblies/SOAP_pipa_genome_chimerical_43mers_supercontigs.index
#6	scaffold173038	34782481	34793964
```
Extract from vcf file the genotypes corresponding to chimerical `scaffold173038`
```
vcftools --vcf /4/caroline/Pipa_parva/HiSeq_analysis/Pipa_chimerical_recalibrated_round1.vcf --chr "supercontig_6" --from-bp 34782481 --to-bp 34793964 --recode --recode-INFO-all --out /4/caroline/Pipa_parva/HiSeq_analysis/Pipa_chimerical_recalibrated_round1_scaffold173038
```
```
/usr/local/vcftools/src/perl/vcf-to-tab < /4/caroline/Pipa_parva/HiSeq_analysis/Pipa_chimerical_recalibrated_round1_scaffold173038.recode.vcf > /4/caroline/Pipa_parva/HiSeq_analysis/Pipa_chimerical_recalibrated_round1_scaffold173038.recode.tab
```
```
#CHROM  POS     REF     BJE4294 BJE4295
...
supercontig_6   34791518        AT      AT/A    A/A
supercontig_6   34791529        AT      AT/A    A/A
supercontig_6   34791536        A       A/AG    AG/AG
supercontig_6   34791537        T       T/TTCCA TTCCA/TTCCA
supercontig_6   34791615        GTTTA   GTTTA/G G/G
```
Positions on `scaffold173038`:
```
34791518-34782481=9037
34791529-34782481=9048
34791536-34782481=9055
34791537-34782481=9056
34791615-34782481=9134
```
Location of `wnt9a`: `scaffold173038 4743-4256`.

**Note:** `vcf` files zipped using `gzip /4/caroline/Pipa_parva/HiSeq_analysis/*.vcf`. To directly obtain genotypes from a region using `vcf.gz` files: `vcftools --gzvcf FILE [...]`.

- Identifying if SNP in mom different from dad in `sox4`
```
blastn -evalue 1e-20 -query /4/caroline/Pipa_parva/blast_genes/Sox4_xenbase_Xtrop.fa -db /4/caroline/2017_Pipoidea_Hiseq/Assemblies/SOAP_pipa_genome_chimerical_43mers_blastable -out /4/caroline/Pipa_parva/blast_genes/Pipa_chimerical_Xtrop_sox4_e20 -outfmt 6 -max_target_seqs 1
grep "scaffold51301" /4/caroline/2017_Pipoidea_Hiseq/Assemblies/SOAP_pipa_genome_chimerical_43mers_supercontigs.index
#2	scaffold51301	194665281	194689488
vcftools --vcf /4/caroline/Pipa_parva/HiSeq_analysis/Pipa_chimerical_recalibrated_round1.vcf --chr "supercontig_2" --from-bp 194665281 --to-bp 194689488 --recode --recode-INFO-all --out /4/caroline/Pipa_parva/HiSeq_analysis/Pipa_chimerical_recalibrated_round1_scaffold51301
vcftools --vcf /4/caroline/Pipa_parva/HiSeq_analysis/Pipa_chimerical_nonrecal_varonly.vcf --chr "supercontig_2" --from-bp 194665281 --to-bp 194689488 --recode --recode-INFO-all --out /4/caroline/Pipa_parva/HiSeq_analysis/Pipa_chimerical_non_recal_scaffold51301
/usr/local/vcftools/src/perl/vcf-to-tab </4/caroline/Pipa_parva/HiSeq_analysis/Pipa_chimerical_recalibrated_round1_scaffold51301.recode.vcf >/4/caroline/Pipa_parva/HiSeq_analysis/Pipa_chimerical_recalibrated_round1_scaffold51301.recode.tab
```
- Identifying if SNP in mom different from dad in `rab2a`
```
blastn -evalue 1e-20 -query /4/caroline/Pipa_parva/blast_genes/rab2a_xenbase_Xtrop.fa -db /4/caroline/2017_Pipoidea_Hiseq/Assemblies/SOAP_pipa_genome_chimerical_43mers_blastable -out /4/caroline/Pipa_parva/blast_genes/Pipa_chimerical_Xtrop_rab2a_e20 -outfmt 6 -max_target_seqs 1
grep "scaffold210848" /4/caroline/2017_Pipoidea_Hiseq/Assemblies/SOAP_pipa_genome_chimerical_43mers_supercontigs.index
#6	scaffold210848	166408241	166409321
vcftools --vcf /4/caroline/Pipa_parva/HiSeq_analysis/Pipa_chimerical_recalibrated_round1.vcf --chr "supercontig_6" --from-bp 166408241 --to-bp 166409321 --recode --recode-INFO-all --out /4/caroline/Pipa_parva/HiSeq_analysis/Pipa_chimerical_recalibrated_round1_scaffold210848

```
- Identifying SNP near exons for some genes important in development

Looked for SNP near exons for: `ctnnb1`, `hoxa11`, `sox4`, `NCOA2`, `wnt9a`, `cdk6`, `cyp51a`, `evx1`, `cbx3`, `cdh7`.

`Scaffold1250`: first exons of `NCOA2`. Using genotypes from HiSeq: position 19401 in intron 1 interesting (supercontig_1: 11067322): heterozygous `G/A` for dad, homozygous `G/G` for mom. Can amplify exon 2 (~180bp) with part of intron 1 containing the SNP: `scaffold1250` ~ 19200..20100 (~900bp). Position 18885 (supercontig_1	11066806, mom:T/T dad:T/TG) and 17755 (supercontig_1: 11065676, mom: G/G, dad: G/GT). Exon 1 from 18609..18694, should amplify from 18400..19300 (~900bp).

- Identifying if SNP in mom different from dad in `mmp16` (Chr06:110830741..110898268)

`Scaffold11527` contains a SNP with a sex inheritance pattern identified by RADseq (position 17118). It maps to an intronic region of `mmp16` of *X. tropicalis* v9.0 on xenbase. 
```
blastn -evalue 1e-1 -query /4/caroline/Pipa_parva/blast_genes/mmp16_cds_utr_xenbase_Xtrop.fa -db /4/caroline/2017_Pipoidea_Hiseq/Assemblies/SOAP_pipa_genome_chimerical_43mers_blastable -out /4/caroline/Pipa_parva/blast_genes/Pipa_chimerical_Xtrop_mmp16_cds_utr_e1_nomaxtarget -outfmt 6

#mmp16	scaffold177758	83.492	1363	151	40	1363	2675	12112	10774	0.0	1203
#mmp16	C50588817	94.156	308	18	0	274	581	629	322	1.58e-129	470
#mmp16	scaffold90826	88.444	225	25	1	734	958	3538	3761	1.72e-69	270
#mmp16	scaffold156786	91.946	149	12	0	1	149	910	762	3.81e-51	209
#mmp16	C48796834	94.915	118	6	0	1247	1364	110	227	6.42e-44	185
#mmp16	scaffold64654	86.164	159	20	2	944	1100	6123	6281	1.80e-39	171
```
We blast the cds of this gene onto our genome to try to identify sex-specific SNP using HiSeq data in the exons. 
```
grep "scaffold177758" /4/caroline/2017_Pipoidea_Hiseq/Assemblies/SOAP_pipa_genome_chimerical_43mers_supercontigs.index
#6	scaffold177758	52788721	52802656
vcftools --gzvcf /4/caroline/Pipa_parva/HiSeq_analysis/Pipa_chimerical_nonrecal_varonly.vcf.gz --chr "supercontig_6" --from-bp 52788721 --to-bp 52802656 --recode --recode-INFO-all --out /4/caroline/Pipa_parva/HiSeq_analysis/Pipa_chimerical_non_recal_scaffold177758
/usr/local/vcftools/src/perl/vcf-to-tab < /4/caroline/Pipa_parva/HiSeq_analysis/Pipa_chimerical_non_recal_scaffold177758.recode.vcf > /4/caroline/Pipa_parva/HiSeq_analysis/Pipa_chimerical_non_recal_scaffold177758.recode.tab
```
The `scaffold177758` maps to the last exon (~335bp long) of the gene from xenbase. A SNP 
```
supercontig_6	52799693	G	G/G	G/A
```
at position `10972` (52799693-52788721=10972) of the `scaffold177758` has been identified.

Need to be carefull when designing the primers. At position 10694, heterozygous with deletions. 

#### Checking with the sex-specific assembly
```
module load blast/2.2.28+
gunzip -c /work/cauretc/2017_pipoidea/Assemblies/scaffold177758.fa.gz| blastn -evalue 1e-60 -query - -db /work/cauretc/2017_pipoidea/Assemblies/SOAP_pipa_genome_dad_43mers_blastable -out /work/cauretc/2017_pipoidea/blast_results/Pipa_dad_chim_scaffold177758_maxtarget1 -outfmt 6 -max_target_seqs 1
gunzip -c /work/cauretc/2017_pipoidea/Assemblies/scaffold177758.fa.gz| blastn -evalue 1e-60 -query - -db /work/cauretc/2017_pipoidea/Assemblies/SOAP_pipa_genome_mom_43mers_blastable -out /work/cauretc/2017_pipoidea/blast_results/Pipa_mom_chim_scaffold177758_maxtarget1 -outfmt 6 -max_target_seqs 1
gunzip -c /work/cauretc/2017_pipoidea/Assemblies/scaffold177758.fa.gz| blastn -evalue 1e-60 -query - -db /work/cauretc/2017_pipoidea/Assemblies/SOAP_pipa_genome_mom_43mers_blastable -out /work/cauretc/2017_pipoidea/blast_results/Pipa_mom_chim_scaffold177758_nomaxtarget -outfmt 6

gunzip -c /work/cauretc/2017_pipoidea/Assemblies/SOAP_pipa_genome_dad_43mers.scafSeq.gz | awk -v seq="scaffold178622" -v RS='>' '$1 == seq {print RS $0}' - >/work/cauretc/2017_pipoidea/Assemblies/SOAP_pipa_genome_dad_scaffold178622.fa
gunzip -c /work/cauretc/2017_pipoidea/Assemblies/SOAP_pipa_genome_mom_43mers.scafSeq.gz | awk -v seq="scaffold174148" -v RS='>' '$1 == seq {print RS $0}' - >/work/cauretc/2017_pipoidea/Assemblies/SOAP_pipa_genome_mom_scaffold174148.fa
cat /work/cauretc/2017_pipoidea/Assemblies/SOAP_pipa_genome_dad_scaffold178622.fa /work/cauretc/2017_pipoidea/Assemblies/SOAP_pipa_genome_mom_scaffold174148.fa /work/cauretc/2017_pipoidea/Assemblies/scaffold177758.fa >/work/cauretc/2017_pipoidea/Assemblies/SOAP_pipa_genome_dad_mom_chim_scaffold177758.fa
```
Then used `MAFFT version 7` [online](mafft.cbrc.jp/alignment/software/).
The alignment confirm the absence of the repeat of Ts in the mom. In order to use the SNP identified in the mom and dad and check for sex inheritance pattern we need to not include this region of T repeats.

- Identifying if SNP in mom different from dad near `kctd1` (Chr06:94888462..94915523)

`scaffold96475` had a SNP with a XY sex inheritance pattern. It blasted in a region surrounded by `kcdtd1`, `taf4b`, `ss18`, `znf521-like`. 
```
blastn -evalue 1e-20 -query /4/caroline/Pipa_parva/blast_genes/kctd1_cds_xenbase_Xtrop.fa -db /4/caroline/2017_Pipoidea_Hiseq/Assemblies/SOAP_pipa_genome_chimerical_43mers_blastable -out /4/caroline/Pipa_parva/blast_genes/Pipa_chimerical_Xtrop_kctd1_cds_e20_nomaxtarget -outfmt 6
```
```
kctd1	scaffold297564	84.342	1437	205	15	309	1738	653	2076	0.0	1389
kctd1	scaffold68613	87.658	316	38	1	1892	2207	4861	4547	1.36e-98	366
kctd1	scaffold68613	91.195	159	14	0	2204	2362	3670	3512	1.45e-53	217
```
`scaffold68613` contains the 2 last exons of `kctd1`.
```
vcftools --gzvcf /4/caroline/Pipa_parva/HiSeq_analysis/Pipa_chimerical_nonrecal_varonly.vcf.gz --chr "supercontig_7" --from-bp 160739841 --to-bp 160745253 --recode --recode-INFO-all --out /4/caroline/Pipa_parva/HiSeq_analysis/Pipa_chimerical_non_recal_scaffold297564
/usr/local/vcftools/src/perl/vcf-to-tab < /4/caroline/Pipa_parva/HiSeq_analysis/Pipa_chimerical_non_recal_scaffold297564.recode.vcf > /4/caroline/Pipa_parva/HiSeq_analysis/Pipa_chimerical_non_recal_scaffold297564.recode.tab
vcftools --gzvcf /4/caroline/Pipa_parva/HiSeq_analysis/Pipa_chimerical_nonrecal_varonly.vcf.gz --chr "supercontig_3" --from-bp 111060401 --to-bp 111069555 --recode --recode-INFO-all --out /4/caroline/Pipa_parva/HiSeq_analysis/Pipa_chimerical_non_recal_scaffold68613
/usr/local/vcftools/src/perl/vcf-to-tab < /4/caroline/Pipa_parva/HiSeq_analysis/Pipa_chimerical_non_recal_scaffold68613.recode.vcf > /4/caroline/Pipa_parva/HiSeq_analysis/Pipa_chimerical_non_recal_scaffold68613.recode.tab
```
`scaffold68613` at position 3755 (supercontig_3:111064156 - 111060401=3755): mom=C/C, dad=C/T; position, 3214 (supercontig_3:111063615): C/C vs C/A; 4437 (supercontig_3:111064838): CATTAA/CATTAA vs C/CATTAA. 
## Primers

### NCOA2
*NCOA2 (nuclear receptor coactivator 2)*

- NCOA2 (exon2 - 180bp - from 19822 to 19992 of `Scaffold1250`). 

SNP at position 19401. Primer3 on NCBI: Forward primer from 19200 - Reverse primer to 20500 ; PCR product size: Min 500, Max 1500 ; Primer Pair Specificity Checking Parameters: Xenopus (taxid:8353) ; other parameters = default values.
```
Primer pair 2
	Sequence (5'->3')	Template strand	Length	Start	Stop	Tm	GC%	Self complementarity	Self 3' complementarity
Forward primer	CTGTCTGGCAAATTCACACCC	Plus	21	19327	19347	59.73	52.38	4.00	1.00
Reverse primer	GTGTGGCCTAAAGCACCAAC	Minus	20	20131	20112	59.69	55.00	4.00	2.00
Product length	805
```

Designed a pair further from the SNP. Forward primer from 19100 - Reverse primer to 20200. PCR product size: Min 800, Max 1100.
```
Primer pair 11
        Sequence (5'->3')       Template strand Length  Start   Stop    Tm      GC%     Self complementarity    Self 3' complementarity
Forward primer  TGTAAAATCGCAACCACGACA   Plus    21      19253   19273   58.79   42.86   3.00    0.00
Reverse primer  TGGCACTGATGAAAAAGCACTG  Minus   22      20127   20106   59.97   45.45   3.00    3.00
Product length  875
```

- NCOA2 (exon1 - 85bp - from 18609 to 18694 of `Scaffold1250`). 

SNP at 18885. Forward primer from 18200 - Reverse primer to 19500. Other: same as previously.
```
Primer pair 1
	Sequence (5'->3')	Template strand	Length	Start	Stop	Tm	GC%	Self complementarity	Self 3' complementarity
Forward primer	GCAGGGTTAAAGGCACCAGA	Plus	20	18272	18291	60.25	55.00	4.00	0.00
Reverse primer	AGAAAGGAGCCCACATGCAA	Minus	20	19087	19068	59.89	50.00	4.00	2.00
Product length	816

Primer pair 8
	Sequence (5'->3')	Template strand	Length	Start	Stop	Tm	GC%	Self complementarity	Self 3' complementarity
Forward primer	CACCTACGGCTGACTGCTC	Plus	19	18210	18228	60.15	63.16	3.00	1.00
Reverse primer	ACTGGTGGGAACAAAGGTGT	Minus	20	18991	18972	59.37	50.00	3.00	2.00
Product length	782
```
Best primer pair: Primer pair 8 (because of self complementarity).

### mmp16
*mmp16 (matrix metallopeptidase 16): [expressed in the rat testis](http://www.reproduction-online.org/content/134/3/455.full.pdf), [one of DMRT1-bound genes](http://www.pnas.org/content/107/30/13360.full.pdf)*

Need to include position 10972 of `scaffold177758`. and try to amplify the whole exon (~ position 12112..10774) and avoid position 10694. 

Forward primer from 10700 - Reverse primer to 12500. Other: same as previously.
```
Primer pair 4
	Sequence (5'->3')	Template strand	Length	Start	Stop	Tm	GC%	Self complementarity	Self 3' complementarity
Forward primer	CAGACAGACCAGCATAAGGGG	Plus	21	10917	10937	60.13	57.14	2.00	0.00
Reverse primer	TGCAAGAGTGGGTGTGATGT	Minus	20	11743	11724	59.53	50.00	4.00	0.00
Product length	827
```
The SNP we are interesting in is a little bit close to the primer (~35bp). Looking for primers further from the SNP, the issue is that we will include a potential deletion SNP which can cause not very clean sequences after sequencing. Need to ask Ben boss for his opinion.

Forward primer from 10600 - Reverse primer to 12500. Other: same as previously.
```
Primer pair 4
	Sequence (5'->3')	Template strand	Length	Start	Stop	Tm	GC%	Self complementarity	Self 3' complementarity
Forward primer	AGTCGCTCCCACTGTCATTTT	Plus	21	10692	10712	59.93	47.62	3.00	0.00
Reverse primer	AGTGGGTGTGATGTAGGGCT	Minus	20	11737	11718	60.55	55.00	2.00	2.00
Product length	1046
Primer pair 5
	Sequence (5'->3')	Template strand	Length	Start	Stop	Tm	GC%	Self complementarity	Self 3' complementarity
Forward primer	AGTTGAAAACAGTCGCTCCCA	Plus	21	10682	10702	60.13	47.62	3.00	0.00
Reverse primer	GCCTTTCCTCTGCACATCGT	Minus	20	11495	11476	60.67	55.00	4.00	2.00
Product length	814
Primer pair 6
	Sequence (5'->3')	Template strand	Length	Start	Stop	Tm	GC%	Self complementarity	Self 3' complementarity
Forward primer	CCAGACAGACCAGCATAAGGG	Plus	21	10916	10936	60.13	57.14	2.00	0.00
Reverse primer	CAACACAGCAAGCACTGTGA	Minus	20	11883	11864	59.27	50.00	5.00	3.00
Product length	968
Primer pair 7
	Sequence (5'->3')	Template strand	Length	Start	Stop	Tm	GC%	Self complementarity	Self 3' complementarity
Forward primer	GACAGACCAGCATAAGGGGA	Plus	20	10919	10938	58.80	55.00	2.00	0.00
Reverse primer	TGCCTTTGAGTGGGGATTCA	Minus	20	11593	11574	59.22	50.00	3.00	2.00
Product length	675
```
From the last primer pairs: the Primer pair 4.

New primers further from SNP. Forward primer from 10700 - Reverse primer to 12100. PCR product size: Min 800, Max 1500.
```
Primer pair 4
	Sequence (5'->3')	Template strand	Length	Start	Stop	Tm	GC%	Self complementarity	Self 3' complementarity
Forward primer	AGTCGCTCCCACTGTCATTTT	Plus	21	10742	10762	59.93	47.62	3.00	0.00
Reverse primer	TGGACAACACAGCAAGCACT	Minus	20	11937	11918	60.39	50.00	2.00	2.00
Product length	1196
```
### kctd1
*kctd1 (potassium channel tetramerization domain containing 1)*

Need to include position 3755 and 3214 of `scaffold68613`.

Forward primer from 3100 - Reverse primer to 4000 ; PCR product size: Min 600, Max 1500. Other: same as previously.
```
Primer pair 4
	Sequence (5'->3')	Template strand	Length	Start	Stop	Tm	GC%	Self complementarity	Self 3' complementarity
Forward primer	TGGAAGGATTTTACCAGCGGC	Plus	21	3129	3149	60.95	52.38	4.00	3.00
Reverse primer	CACTTGCTGTGTGCAGTGAAA	Minus	21	3890	3870	59.87	47.62	4.00	2.00
Product length	762

Primer pair 8
	Sequence (5'->3')	Template strand	Length	Start	Stop	Tm	GC%	Self complementarity	Self 3' complementarity
Forward primer	TTTTACCAGCGGCTAAGAGG	Plus	20	3137	3156	57.60	50.00	6.00	2.00
Reverse primer	AATCACTTGCTGTGTGCAGT	Minus	20	3893	3874	58.32	45.00	4.00	1.00
Product length	757

Primer pair 10
	Sequence (5'->3')	Template strand	Length	Start	Stop	Tm	GC%	Self complementarity	Self 3' complementarity
Forward primer	CCATCCTTAAGTTGCCGTCCT	Plus	21	3299	3319	60.07	52.38	6.00	0.00
Reverse primer	CACCCTTTAGTCAATAAGGCCC	Minus	22	3953	3932	58.71	50.00	5.00	3.00
Product length	655
```
Best primer pair: Primer pair 4.

Checked with geneious, the forward primer is too close from the 1st SNP. Seems like the fasta file needs to not contained the name in the fasta sequence for the NCBI version otherwise the location of the primers is wrong. Designed a new pair (without `>scaffold68613 28.4` in the fasta file on NCBI) that seems better:
```
Primer pair 9
	Sequence (5'->3')	Template strand	Length	Start	Stop	Tm	GC%	Self complementarity	Self 3' complementarity
Forward primer	ATGTAGTGGTTGCCTTGTGCC	Plus	21	3142	3162	61.16	52.38	2.00	1.00
Reverse primer	GCTCCAGTCCCAAGTCATAGT	Minus	21	3866	3846	59.17	52.38	2.00	1.00
Product length	725
```

#### Amplifications
Sent for sequencing on Sept.18 samples from the whole *Pipa* family except one individual for `mmp16_e9` (PCR 23/08), half of the family for `NCOA2_e2` (PCR 22/08, gradient from 09/08 and PCR 16/09) and `KCTD1_e4` (PCR 23/08, gradient from 09/08 and PCR 16/09). Should be ready for checking the quality and looking for the SNP inheritance pattern by Tomorrow (Sept.20). In case of bad sequences, should 1st try to amplify some individuals at a higher annealing temperature (+1 to increase the specificity). If good sequences and differences between males and females: need to amplify the 7 other individuals for `NCOA2` and `KCTD1`.

The sequences came back on Sept.9. They look pretty good for `NCOA2` and `KCTD1`. 2 snps show up for `NCOA2`, 1 SNP for `KCTD1` (and an indel but the sequence from the forward is very clean around the SNP). For `mmp16`, I think we are in a region with multiple indels, additional investigations are necessary and nested primers would be necessary. For now I'll focus on amplifying the whole family (7 individuals left) for `NCOA2` and `KCTD1`.

**28/09:** Ok so I actually think that for mmp16, the amplification was not specific enough. So 2 solutions: amplify at higher temperature (during the gradient we also had a band at 64.2); or having a nested primer to make sure we have a very good quality for our region of interest.
```
Primer pair 1
	Sequence (5'->3')	Template strand	Length	Start	Stop	Tm	GC%	Self complementarity	Self 3' complementarity
Forward primer	AGTCGCTCCCACTGTCATTTT	Plus	21	10742	10762	59.93	47.62	3.00	0.00
Reverse primer	TTTTCCCATCAAACAAGGGGC	Minus	21	11208	11188	59.30	47.62	3.00	2.00
Product length	467
```
Compared to other new reverse tested, this one has a closer Tm with the forward and not too big self complementarity values. The fact that the expected size is <500bp, it will most likely allow us to have both sequences (F&R) overlapping.

If BE agrees: I will amplify at a higher temperature and send for sequencing only the forward.

**4/10: NCOA2/KCTD1:** Gel from 3/10 of PCR 2/10: Last individuals for `NCOA2` and `KCTD1` amplified (in addition to the ones from PCR 24/09) - need to be exosap + sequence reaction today

**4/10: MMP16:** Gel from 3/10 of PCR 2/10 (PCR at 2 higher Tannealing 58 and 59): most individuals amplified and need to be clean + sequence reaction today

**3/10: Summary of PCR**

##### NCOA2_e2 (`scaffold1250`, `position19402`)
All the individuals except `BJE4300` (son) amplified and got OK sequences. Mother: G/G, father: A/G, 6 sons: A/G, 5 daughters: G/G.

29/03 - update: `BJE4300` sequenced: A/G

##### KCTD1_e4 (`scaffold68613`, `position3284`)
All the individuals amplified and got OK sequences. Mother: C/C, father: C/T, 7 sons: C/T, 5 daughters: C/C.

##### mmp16_e9 (`scaffold177758`, `position10973`)
All the individuals except `BJE4308` and `BJE4309` (1 son, 1 daughter) amplified and got OK sequences. Mother: G/G, father: A/G, 6 sons: A/G, 4 daughters: G/G.

29/03 - update: `BJE4308`: A/G, `BJE4309` G/G.
