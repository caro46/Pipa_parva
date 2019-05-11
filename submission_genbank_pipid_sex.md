# tbl2asn



# Renaming the files

`renaming_fasta_from_table.pl`
```perl

#!/usr/local/perl5.24/perl-5.24.0/perl

use strict;
use warnings;

my ($fasta_file, $table_name, $new_name_NCBI_fasta) = @ARGV;

if (not defined $fasta_file) {
  die "Need a fasta file \n";
}

if (not defined $table_name) {
  die "Need a file with the new sequence names\n";
}

if (not defined $new_name_NCBI_fasta) {
  die "Need an output name \n";
}

unless (open DATAINPUT, $fasta_file) {
                print "Can not find the input file $fasta_file.\n";
                exit;
}

unless (open DATAINPUT1, $table_name) {
                print "Can not find the input file $table_name.\n";
                exit;
}

unless (open(OUTFILE, ">$new_name_NCBI_fasta")){
        print "I can\'t write to $new_name_NCBI_fasta $!\n\n";
        exit;
}

#Usage:
#perl renaming_fasta_from_table.pl All_DMRT1_exon3_align_for_analysis_nobadstuff_nice.nex.fasta renaming_dmrt1_dmw_e3_alignment.tsv Evans_dmrt1_dmw_e3_alignment_xenopus.fsa

my %names;
my @columnF1;
my @temp;
my $switch=0;
my $complete_sequence=();
my $sequence_name=();

while (my $line1 = <DATAINPUT1>){
        chomp($line1);
        @columnF1=split("\t",$line1);
                $names{$columnF1[0]}=$columnF1[1];
}

while (my $line=<DATAINPUT>) {
    chomp($line);
    @temp=split(/[>\s]/,$line);

    if($switch==1) {
        if ($line!~/^>/){
            $complete_sequence=$complete_sequence.$temp[0];
        }
        else {
            $switch=0;
            print OUTFILE ">",$sequence_name,"\n";
            until(length($complete_sequence) < 80){
            print OUTFILE substr($complete_sequence, 0, 80),"\n";
            $complete_sequence = substr($complete_sequence,80);
        }
        print OUTFILE $complete_sequence,"\n";
        $complete_sequence=();
        }
    }

    if(defined($temp[1])){
        foreach my $old_name (keys %names) {
                if(grep $_ eq $temp[1], $old_name){
                    $sequence_name=$names{$old_name};
                    #@temp2 = split(' ',$sequence_name);
                    #$isolate = $temp2[0];
                    $switch=1;
                }
        }
    }
}
print OUTFILE ">",$sequence_name,"\n";
until(length($complete_sequence) < 80){
    print OUTFILE substr($complete_sequence, 0, 80),"\n";
    $complete_sequence = substr($complete_sequence,80);
}
print OUTFILE $complete_sequence,"\n";
$complete_sequence=();

close DATAINPUT;
print "Done with input file \n";
close DATAINPUT1;
print "Done with input file 1\n";
close OUTFILE;
print "Done with output file \n";

```

`geneious_to_tbl.pl`
```perl
#!/usr/local/perl5.24/perl-5.24.0/perl

use strict;
use warnings;

my ($geneious_table, $table_name, $NCBI_tbl_format) = @ARGV;

if (not defined $geneious_table) {
  die "Need a geneious annotation file \n";
}

if (not defined $table_name) {
  die "Need a file with the new sequence names\n";
}

if (not defined $NCBI_tbl_format) {
  die "Need an output name \n";
}

unless (open DATAINPUT, $geneious_table) {
                print "Can not find the input file $geneious_table.\n";
                exit;
}

unless (open DATAINPUT1, $table_name) {
                print "Can not find the input file $table_name.\n";
                exit;
}

unless (open(OUTFILE, ">$NCBI_tbl_format")){
        print "I can\'t write to $NCBI_tbl_format $!\n\n";
        exit;
}

#Usage:
#perl geneious_to_tbl.pl geneious_table_annotation_DMRT1_DMW_algt_e3.csv renaming_dmrt1_dmw_e3_alignment.tsv Xenopus_DMRT1_DMW_exon3.tbl

my %names;
my @columnF1;
my @temp;
my @temp2;
my @temp3;
my $seq_name;
my $annotation_name;
my $type_annotation;
my $start_annotation;
my $end_annotation;
my $max_length;
my $length;
my $sequence_name=();
my $annotation;
my $gene_name;

while (my $line1 = <DATAINPUT1>){
        chomp($line1);
        @columnF1=split("\t",$line1);
        @temp2=split(" ",$columnF1[1]);
                $names{$columnF1[0]}=$temp2[0];
}

while (my $line=<DATAINPUT>) {
        chomp($line);
        @temp=split(",",$line);
#       @temp2=split(" ",$temp[0]);
        $seq_name=$temp[0];
        $annotation_name=$temp[1];
        @temp3=split("_",$annotation_name);
        $gene_name=$temp3[0];
        $type_annotation=$temp[2];
        $start_annotation=$temp[3];
        $end_annotation=$temp[4];
        $max_length=$temp[5];
        $length=$temp[6];

        if(defined($seq_name)){
                foreach my $old_name (keys %names) {
                        if(grep $_ eq $seq_name, $old_name){
                                $sequence_name=$names{$old_name};
                                        if ($type_annotation eq "gene" ){
                                        print OUTFILE ">Feature \t", $sequence_name, "\n",
                                        "<",$start_annotation,"\t",">",$end_annotation,"\t",$type_annotation,"\n";
                                                if ($gene_name eq "dmrt1"){
                                                print OUTFILE "\t\t\t","gene\tdmrt1","\n";
                                                }
                                                else{
                                                print OUTFILE "\t\t\t","gene\tx-dmw","\n";
                                                }
                                        }
                                        elsif ($type_annotation eq "exon"){
                                        print OUTFILE $start_annotation,"\t",$end_annotation,"\t",$type_annotation,"\n";
                                        }
                                        else{print OUTFILE
                                        "<",$start_annotation,"\t",">",$end_annotation,"\t",$type_annotation,"\n";
                                        if ($type_annotation eq "CDS" ){
                                                if ($gene_name eq "dmrt1"){
                                                print OUTFILE "\t\t\t","product\tdoublesex- and mab-3-related transcription factor 1A","\n",
                                                "\t\t\t","codon_start\t1\n";
                                                        if ( grep $_ eq "pseudo", @temp3 ) {
                                                        print OUTFILE "\t\t\t","pseudo","\n",
                                                                "\t\t\t","note","\t","stop codon inside predicted CDS","\n";
                                                        }
                                                        else{
                                                        next}
                                                }else{
                                                print OUTFILE "\t\t\t","product\tdoublesex- and mab-3-related transcription factor DM-W","\n",
                                                "\t\t\t","codon_start\t1\n";
                                                }
                                        }
                                        }
                        }
                }
        }
}
#print OUTFILE ">Feature \t", $sequence_name, "\n", $start_annotation,"\t",$end_annotation,"\t",$type_annotation,"\n";


close DATAINPUT;
print "Done with input file \n";
close DATAINPUT1;
print "Done with input file 1\n";
close OUTFILE;
print "Done with output file \n";

```

# Batch

```
./mac.tbl2asn -t /Users/evanslab/Documents/caroline/Publi_pipid_SC/NCBI/submission/template.sbt -p /Users/evanslab/Documents/caroline/Publi_pipid_SC/NCBI/submission/ -r /Users/evanslab/Documents/caroline/Publi_pipid_SC/NCBI/submission/ -V v -a s
```

# Alignment

```
./mac.tbl2asn -t /Users/evanslab/Documents/caroline/Publi_pipid_SC/NCBI/submission/template.sbt -p /Users/evanslab/Documents/caroline/Publi_pipid_SC/NCBI/submission/alignments/ -r /Users/evanslab/Documents/caroline/Publi_pipid_SC/NCBI/submission/alignments/ -V v -a l1
```
