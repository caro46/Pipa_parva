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
