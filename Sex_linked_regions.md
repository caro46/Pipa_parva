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
Looked for SNP near exons for: `ctnnb1`, `hoxa11`, `sox4`, `NCOA2`, `wnt9a`, `cdk6`, `cyp51a`, `evx1`, `cbx3`, `cdh7`.

`Scaffold1250`: first exons of `NCOA2`. Using genotypes from HiSeq: position 19401 in intron 1 interesting (supercontig_1: 11067322): heterozygous `G/A` for dad, homozygous `G/G` for mom. Can amplify exon 2 (~180bp) with part of intron 1 containing the SNP: `scaffold1250` ~ 19200..20100 (~900bp). Position 18885 (supercontig_1	11066806, mom:T/T dad:T/TG) and 17755 (supercontig_1: 11065676, mom: G/G, dad: G/GT). Exon 1 from 18609..18694, should amplify from 18400..19300 (~900bp).

## Primers
NCOA2 (exon2 - 180bp - from 19822 to 19992 of Scaffold1250). SNP at position 18885. Primer3 on NCBI: Forward primer from 19200 - Reverse primer to 20500 ; PCR product size: Min 500, Max 1500 ; Primer Pair Specificity Checking Parameters: Xenopus (taxid:8353) ; other parameters = default values.
```
Primer pair 2
	Sequence (5'->3')	Template strand	Length	Start	Stop	Tm	GC%	Self complementarity	Self 3' complementarity
Forward primer	CTGTCTGGCAAATTCACACCC	Plus	21	19327	19347	59.73	52.38	4.00	1.00
Reverse primer	GTGTGGCCTAAAGCACCAAC	Minus	20	20131	20112	59.69	55.00	4.00	2.00
Product length	805
```
NCOA2 (exon2 - 85bp - from 18609 to 18694 of Scaffold1250). SNP at 18885. Forward primer from 18200 - Reverse primer to 19500. Other: same as previously.
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
