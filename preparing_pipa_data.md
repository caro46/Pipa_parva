# Sequence quality
Check for over-represented sequences, adapters (need all the sequences) 
```
/usr/local/fastqc/fastqc -o /4/caroline/Pipa_parva/Rad_seq/ /4/caroline/1016_S2_L002_R1_001.fastq.gz
```
# STACKS: process_radtags
Installing locally the new version of stacks (need new version of gcc see [here](https://github.com/caro46/Hymenochirus/blob/master/some_commands.Rmd) for how to install gcc)
```
export CC=/home/caroline/programs/gcc-6.2.0/bin/gcc
export CXX=/home/caroline/programs/gcc-6.2.0/bin/g++
export CPP=/home/caroline/programs/gcc-6.2.0/bin/cpp
export LD=/home/caroline/programs/gcc-6.2.0/bin/gcc

./configure --prefix=/home/caroline/programs/stacks-1.45 LDFLAGS /home/caroline/programs/gcc-6.2.0/lib64
make
make install
```

```
/usr/local/bin/process_radtags -f /4/caroline/563_S7_L007_R1_001.fastq.gz -b /4/caroline/Pipa_parva/Rad_seq/demultiplex/demutiplex_barcode -o /4/caroline/Pipa_parva/Rad_seq/demultiplex/ -e sbfI -t 75 -r -c -q --filter_illumina
```
