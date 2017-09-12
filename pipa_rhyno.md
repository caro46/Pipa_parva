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
#Bus error
```
Server failed, cannot access any file (`permission denied`on Aug 7). Need to be rerun later (issues because of unpredicted shutdown, can take a while to fix everything).
Again `Segmentation fault`.
```
/work/cauretc/programs/SOAPdenovo2-src-r240/SOAPdenovo-63mer scaff -g /work/cauretc/2017_pipoidea/Assemblies/SOAP_rhyno_genome_63mers -F -p 10 1 >scaff_rhyno.log 2>scaff_rhyno.err
```
Tried to run the SOAP de novo assembly multiple times since the server is stable (3 times since the server fixation) using different values of kmers to try to use less memory but keep having `Segmentation fault`. 

Trying with Abyss:
```
module unload intel mkl openmpi
module load gcc/4.9.2
module load openmpi/gcc492-std/1.8.7
module load boost/gcc492-openmpi187std/1.59.0
export PATH=/work/ben/abyss/bin:$PATH

abyss-pe np=8 name=CSL6209 k=64 in='/work/cauretc/2017_pipoidea/2017_Pipa_Rhino_genomes/CSL6209_S22_L003_R1_001_trim_paired_Scythe.cor.fastq.gz /work/cauretc/2017_pipoidea/2017_Pipa_Rhino_genomes/CSL6209_S22_L003_R2_001_trim_paired_Scythe.cor.fastq.gz'
```
```
Building the suffix array...
Building the Burrows-Wheeler transform...
Building the character occurrence table...
sort: write failed: /tmp/sortjgKeIs: No space left on device
error: `CSL6209-3.hist': No such file or directory
make: *** [CSL6209-3.dist] Error 1
make: *** Deleting file `CSL6209-3.dist'
``` 
Suggestions from [here](https://groups.google.com/forum/#!topic/abyss-users/x8Wd0tEnyIw) to change the Environment Variables and redirect the temporary directory to somewhere with more space 
```
export TMPDIR=/scratch/cauretc/temp
```
```
libgomp: Thread creation failed: Resource temporarily unavailable
Mateless   0
Unaligned  0
Singleton  0
FR         0
RF         0
FF         0
Different  0
Total      0
abyss-fixmate: error: All reads are mateless. This can happen when first and second read IDs do not match.
error: `CSL6209-3.hist': No such file or directory
make: *** [CSL6209-3.dist] Error 1
make: *** Deleting file `CSL6209-3.dist'
```

Deleted everything from the folder `/work/cauretc/2017_pipoidea/Assemblies/abyss/` to start again from scratch (Sept. 5) and using `np=4` but ended (Sept. 6) with
```
[iqaluk:03683] *** Process received signal ***
[iqaluk:03683] Signal: Bus error (7)
[iqaluk:03683] Signal code: Non-existant physical address (2)
[iqaluk:03683] Failing at address: 0x7fbdf4ebe006
[iqaluk:03683] [ 0] /lib64/libc.so.6[0x3d40432660]
[iqaluk:03683] [ 1] /lib64/libc.so.6(memcpy+0x15b)[0x3d4048995b]
[iqaluk:03683] [ 2] /opt/sharcnet/openmpi/1.8.7/gcc-4.9.2/std/lib/openmpi/mca_btl_vader.so(+0x2639)[0x7fbdf75ed639]
[iqaluk:03683] [ 3] /opt/sharcnet/openmpi/1.8.7/gcc-4.9.2/std/lib/openmpi/mca_pml_ob1.so(mca_pml_ob1_send_request_start_prepare+0x43)[0x7fbdf69222e3]
[iqaluk:03683] [ 4] /opt/sharcnet/openmpi/1.8.7/gcc-4.9.2/std/lib/openmpi/mca_pml_ob1.so(mca_pml_ob1_send+0x848)[0x7fbdf6917948]
[iqaluk:03683] [ 5] /opt/sharcnet/openmpi/1.8.7/gcc-4.9.2/std/lib/libmpi.so.1(MPI_Send+0x153)[0x7fc06199d9e3]
[iqaluk:03683] [ 6] ABYSS-P[0x414553]
[iqaluk:03683] [ 7] ABYSS-P[0x4323bc]
[iqaluk:03683] [ 8] ABYSS-P[0x432249]
[iqaluk:03683] [ 9] ABYSS-P[0x431f6f]
[iqaluk:03683] [10] ABYSS-P[0x41d763]
[iqaluk:03683] [11] ABYSS-P[0x427e30]
[iqaluk:03683] [12] ABYSS-P[0x42306d]
[iqaluk:03683] [13] ABYSS-P[0x417fcb]
[iqaluk:03683] [14] ABYSS-P[0x4182c9]
[iqaluk:03683] [15] ABYSS-P[0x405f64]
[iqaluk:03683] [16] /lib64/libc.so.6(__libc_start_main+0xfd)[0x3d4041ed1d]
[iqaluk:03683] [17] ABYSS-P[0x405769]
[iqaluk:03683] *** End of error message ***
--------------------------------------------------------------------------
mpirun noticed that process rank 1 with PID 3683 on node iqaluk exited on signal 7 (Bus error).
--------------------------------------------------------------------------
make: *** [CSL6209-1.fa] Error 135
```
Then tried to run again (no file in my `/scratch/cauretc/temp/`) and directly stopped with
```
[iqaluk.sharcnet.ca:14769] opal_os_dirpath_create: Error: Unable to create the sub-directory (/scratch/cauretc/temp/openmpi-sessions-cauretc@iqaluk_0) of (/scratch/cauretc/temp/openmpi-sessions-cauretc@iqaluk_0/14457/0/0), mkdir failed [1]
[iqaluk.sharcnet.ca:14769] [[14457,0],0] ORTE_ERROR_LOG: Error in file ../../openmpi-1.8.7/orte/util/session_dir.c at line 107
[iqaluk.sharcnet.ca:14769] [[14457,0],0] ORTE_ERROR_LOG: Error in file ../../openmpi-1.8.7/orte/util/session_dir.c at line 402
[iqaluk.sharcnet.ca:14769] [[14457,0],0] ORTE_ERROR_LOG: Error in file ../../../../../openmpi-1.8.7/orte/mca/ess/hnp/ess_hnp_module.c at line 638
--------------------------------------------------------------------------
It looks like orte_init failed for some reason; your parallel process is
likely to abort.  There are many reasons that a parallel process can
fail during orte_init; some of which are due to configuration or
environment problems.  This failure appears to be an internal failure;
here's some additional information (which may only be relevant to an
Open MPI developer):

  orte_session_dir failed
  --> Returned value Error (-1) instead of ORTE_SUCCESS
--------------------------------------------------------------------------
make: *** [CSL6209-1.fa] Error 213
```
Not sure what is going on and a lot of recent issues on multiple servers so started again on `wob101` (Sept. 6) - before I was running on `iqaluk`.
So `wobbie` is still on unstable set up. The run seems to take forever and right now (Sept.7) no file has been produced. Maybe need to switch again to another server. Killed the job on `wobbie` (nothing was happening and other people also switched because they couldn't run their jobs) and run again on `iqaluk` (Sept.8).

Interesting discussion [here](https://github.com/bcgsc/abyss/issues/104) about empty `.hist` files. 

Iqaluk (Sept.12): assembly seems to continue running pretty well
```
Mateless           0
Unaligned    2288333  1.56%
Singleton   18519513  12.6%
FR          54510728  37.2%
RF              7373  0.00503%
FF             11599  0.00792%
Different   71106123  48.6%
Total      146443669
```
From a [forum](https://groups.google.com/forum/#!topic/abyss-users/bCBAHHn5mQY) (answer from Shaun Jackman to an inquiry about abyss mapping statistics):

*The alignment stats look okay. They indicate that ~11% of the pairs align to the same contig and ~80% align to different contigs, which is likely a result of having small initial contigs due to low coverage.*

This person obtained ~80% Different. Our statistics seem pretty OK too compared to hers (by using the k-mer distribution from `Quake` I already knew about the low coverage). 

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
