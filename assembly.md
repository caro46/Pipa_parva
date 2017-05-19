 # SOAPdenovo - chimerical
 
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
## Improving: filling the gaps

In `/scratch/cauretc/SOAP_assemblies/` where I copied the results from the previous command. We redo the last step of scaffolding, this time trying `-F` to fill in the gaps (for *Hymenochirus* always failed because of the memory)
```
/work/cauretc/programs/SOAPdenovo2-src-r240/SOAPdenovo-63mer scaff -g SOAP_pipa_genome_chimerical_63mers -F -p 10 1 >scaff.log 2>scaff.err
```
