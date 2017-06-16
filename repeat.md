## RepARK 
`version 1.3.0`
```
/home/caroline/programs/RepARK-master/RepARK.pl -l /4/caroline/2017_Pipoidea_Hiseq/after_scythe/BJE4294_S20_L003_R1_001_trim_paired_Scythe.cor.fastq.gz -l /4/caroline/2017_Pipoidea_Hiseq/after_scythe/BJE4294_S20_L003_R2_001_trim_paired_Scythe.cor.fastq.gz -o /4/caroline/Pipa_parva/RepARK_output/BJE4294_repeats
#Couldn't calculate left threshold. Try adjusting the coverage or check your histogram file!
#Couldn't calculate left threshold. Try adjusting the coverage or check your histogram file!
#Theshold used: 56

/home/caroline/programs/RepARK-master/RepARK.pl -l /4/caroline/2017_Pipoidea_Hiseq/after_scythe/BJE4295_S21_L003_R1_001_trim_paired_Scythe.cor.fastq.gz -l /4/caroline/2017_Pipoidea_Hiseq/after_scythe/BJE4295_S21_L003_R2_001_trim_paired_Scythe.cor.fastq.gz -o /4/caroline/Pipa_parva/RepARK_output/BJE4295_repeats
#Theshold used: 62
#Check the directory /4/caroline/Pipa_parva/RepARK_output/BJE4295_repeats for results.
```

## TEclass
`version 2.1.3`
```
/home/caroline/programs/TEclass-2.1.3/TEclassTest.pl -r /4/caroline/Pipa_parva/RepARK_output/BJE4295_repeats/repeat_lib.fasta -o /4/caroline/Pipa_parva/RepARK_output/BJE4295_repeats/TEclass
/home/caroline/programs/TEclass-2.1.3/TEclassTest.pl -r /4/caroline/Pipa_parva/RepARK_output/BJE4294_repeats/repeat_lib.fasta -o /4/caroline/Pipa_parva/RepARK_output/BJE4294_repeats/TEclass
```
If when you try to run you got `Some of your classifiers appear to be missing.` even though you have the classifiers and the good paths in the scripts, try installing again with:
```
sh Compile_dependencies.sh
perl Configure.pl 
sh Install.sh /home/caroline/programs/TEclass-2.1.3
```
If `short forward_vs_reverse ... sh: line 1: 28134 Segmentation fault`, try to re-install everything from downloading.
