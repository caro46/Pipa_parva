# Sequence quality
Check for over-represented sequences, adapters (need all the sequences) 
```
/usr/local/fastqc/fastqc -o /4/caroline/Pipa_parva/Rad_seq/ /4/caroline/1016_S2_L002_R1_001.fastq.gz
```
# STACKS: process_radtags

```
/usr/local/bin/process_radtags -f /4/caroline/563_S7_L007_R1_001.fastq.gz -b /4/caroline/Pipa_parva/Rad_seq/demultiplex/demutiplex_barcode -o /4/caroline/Pipa_parva/Rad_seq/demultiplex/ -e sbfI -t 75 -r -c -q --filter_illumina
```
