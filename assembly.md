 # SOAPdenovo - chimerical
 ### 63mers (K=63)
 ```
 /work/cauretc/programs/SOAPdenovo2-src-r240/SOAPdenovo-63mer all -s /work/cauretc/2017_pipoidea/pipa_chimerical.config -K 63 -R -V -p 10 -o /work/cauretc/2017_pipoidea/Assemblies/SOAP_pipa_genome_chimerical_63mers
 ```
```
<-- Information for assembly Scaffold '/work/cauretc/2017_pipoidea/Assemblies/SOAP_pipa_genome_chimerical_63mers.scafSeq'.(cut_off_length < 100bp) -->

Size_includeN   1836392386
Size_withoutN   1763334449
Scaffold_Num    2805070
Mean_Size       654
Median_Size     127
Longest_Seq     116786
Shortest_Seq    100
Singleton_Num   2509806
Average_length_of_break(N)_in_scaffold  26

Known_genome_size       NaN
Total_scaffold_length_as_percentage_of_known_genome_size        NaN

scaffolds>100   2759077 98.36%
scaffolds>500   313882  11.19%
scaffolds>1K    226995  8.09%
scaffolds>10K   41039   1.46%
scaffolds>100K  7       0.00%
scaffolds>1M    0       0.00%

Nucleotide_A    546977957       29.79%
Nucleotide_C    342179253       18.63%
Nucleotide_G    336184712       18.31%
Nucleotide_T    537992527       29.30%
GapContent_N    73057937        3.98%
Non_ACGTN       0       0.00%
GC_Content      38.47%          (G+C)/(A+C+G+T)

N10     27209   4962
N20     18611   13244
N30     13669   24818
N40     10094   40497
N50     7213    62052
N60     4780    93277
N70     2557    145074
N80     579     291921
N90     127     1221745

NG50    NaN     NaN
N50_scaffold-NG50_scaffold_length_difference    NaN
```

#### Improving: filling the gaps

In `/scratch/cauretc/SOAP_assemblies/` where I copied the results from the previous command. We redo the last step of scaffolding, this time trying `-F` to fill in the gaps (for *Hymenochirus* always failed because of the memory)
```
/work/cauretc/programs/SOAPdenovo2-src-r240/SOAPdenovo-63mer scaff -g SOAP_pipa_genome_chimerical_63mers -F -p 10 1 >scaff.log 2>scaff.err
```
```
<-- Information for assembly Scaffold 'SOAP_pipa_genome_chimerical_63mers.scafSeq'.(cut_off_length < 100bp) -->

Size_includeN   1836116521
Size_withoutN   1800826343
Scaffold_Num    2766152
Mean_Size       663
Median_Size     127
Longest_Seq     117088
Shortest_Seq    100
Singleton_Num   2470888
Average_length_of_break(N)_in_scaffold  12

Known_genome_size       NaN
Total_scaffold_length_as_percentage_of_known_genome_size        NaN

scaffolds>100   2721469 98.38%
scaffolds>500   314265  11.36%
scaffolds>1K    227334  8.22%
scaffolds>10K   41213   1.49%
scaffolds>100K  7       0.00%
scaffolds>1M    0       0.00%

Nucleotide_A    557881575       30.38%
Nucleotide_C    350128959       19.07%
Nucleotide_G    343941812       18.73%
Nucleotide_T    548873997       29.89%
GapContent_N    35290178        1.92%
Non_ACGTN       0       0.00%
GC_Content      38.54%          (G+C)/(A+C+G+T)

N10     27355   4941
N20     18707   13180
N30     13762   24685
N40     10168   40257
N50     7280    61634
N60     4838    92533
N70     2611    143512
N80     612     284285
N90     127     1186311

NG50    NaN     NaN
N50_scaffold-NG50_scaffold_length_difference    NaN
```
### 53mers (K=53)
```
/work/cauretc/programs/SOAPdenovo2-src-r240/SOAPdenovo-63mer all -s /work/cauretc/2017_pipoidea/pipa_chimerical.config -K 53 -R -V -p 10 -F -o /work/cauretc/2017_pipoidea/Assemblies/SOAP_pipa_genome_chimerical_53mers

The final rank

*******************************
 Scaffold number                  304566
 In-scaffold contig number        4259013
 Total scaffold length            1404611141
 Average scaffold length          4611
 Filled gap number                554653
 Longest scaffold                 140357
 Scaffold and singleton number    2259624
 Scaffold and singleton length    1594594292
 Average length                   705
 N50                              8160
 N90                              444
 Weak points                      0

*******************************

#Issue with gap filling

*** glibc detected *** /work/cauretc/programs/SOAPdenovo2-src-r240/SOAPdenovo-63mer: double free or corruption (!prev): 0x000000036424b9e0 ***

#Rerun the scaffolding step

/work/cauretc/programs/SOAPdenovo2-src-r240/SOAPdenovo-63mer scaff -g SOAP_pipa_genome_chimerical_53mers -F -p 10 1 >scaff.log 2>scaff.err

<-- Information for assembly Scaffold 'SOAP_pipa_genome_chimerical_53mers.scafSeq'.(cut_off_length < 100bp) -->

Size_includeN   1703163115
Size_withoutN   1662799375
Scaffold_Num    2238461
Mean_Size       760
Median_Size     107
Longest_Seq     140466
Shortest_Seq    100
Singleton_Num   1933895
Average_length_of_break(N)_in_scaffold  18

Known_genome_size       NaN
Total_scaffold_length_as_percentage_of_known_genome_size        NaN

scaffolds>100   2175265 97.18%
scaffolds>500   307367  13.73%
scaffolds>1K    228359  10.20%
scaffolds>10K   39121   1.75%
scaffolds>100K  6       0.00%
scaffolds>1M    0       0.00%

Nucleotide_A    515135733       30.25%
Nucleotide_C    321717614       18.89%
Nucleotide_G    317140165       18.62%
Nucleotide_T    508805863       29.87%
GapContent_N    40363740        2.37%
Non_ACGTN       0       0.00%
GC_Content      38.42%          (G+C)/(A+C+G+T)

N10     25825   4880
N20     17961   12901
N30     13410   23921
N40     10091   38611
N50     7471    58269
N60     5229    85487
N70     3226    126600
N80     1247    208372
N90     152     673284
```
### 43mers (K=43)
```
/work/cauretc/programs/SOAPdenovo2-src-r240/SOAPdenovo-63mer all -s /work/cauretc/2017_pipoidea/pipa_chimerical.config -K 43 -R -V -p 10 -F -o /work/cauretc/2017_pipoidea/Assemblies/SOAP_pipa_genome_chimerical_43mers
```
