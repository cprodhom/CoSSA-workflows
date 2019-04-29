# Comparative Subsequence Sets Analysis (CoSSA) README

CoSSA (Prodhomme *et al.*, 2019) is a set of workflows used to identify quickly and efficiently haplotype specific SNPs linked to a trait of interest from Whole Genome Sequencing data.

CoSSA can be used for Bulk Segregant Analysis (BSA) using a reference genome. DNA pools of individuals with contrasting phenotypes are sequenced, reads are trimmed and *k*-mer tables are produced for each of the pools. *K*-mers specific to the haplotype of interest can be selected using set algebra between the different *k*-mer sets. The selected *k*-mers can then be mapped to the reference genome and the number of *k*-mers per chromosome bin is counted to identify the locus responsible for the trait variation.

CoSSA can be used for BSA when no reference genome is available as well. The *k*-mers specific to the haplotype of interest can be used to retrieve the paired-end reads from this haplotype. These paired-end reads are then assembled *de novo* to reconstruct the haplotype of interest. The specific *k*-mers are mapped back to the *de novo* assembly to identify haplotype specific SNPs. 

CoSSA can be used to refine the haplotype specificity of the identified SNPs linked with the trait of interest by including extra genotypes in the workflow which do not share the trait of interest.

CoSSA can be used for pedigree analysis. For instance, it can be used to verify if several genotypes sharing the same trait share the same locus involved in the variation of this trait.

## Prerequisites
### Bedtools
Download from https://bedtools.readthedocs.io/en/latest/ 
### GenomeTester4
Download from https://github.com/bioinfo-ut/GenomeTester4
### Samtools
Download from http://samtools.sourceforge.net/ 
### A read trimming software
You can use your favorite reads trimming software to clean your sequencing data. For Illumina data, we used [Trimmomatic](https://github.com/timflutre/trimmomatic).
### A read alignment software
Use your favorite alignment software. It should be adapted to the *k*-mer size you will use in your analysis. As GenomeTester4 creates *k*-mer tables for *k*-mers up to 32 nucleotides, we used [BWA-ALN]( http://bio-bwa.sourceforge.net/) which is designed for Illumina reads of up to 100bp.   
### A *de novo* assembly software
Use your favorite assembler to assemble the PE reads coming from your haplotype of interest. For Illumina 151bp PE reads, we used [SPAdes]( http://cab.spbu.ru/software/spades/). 

## Preprocessing
### Raw data cleaning
Use Trimmomatic or another software to clean the raw data of your pools and the other genotypes included in your CoSSA.
Example:
```
java -jar trimmomatic-0.32.jar PE -phred33 read1.fq read2.fq read1.trim.fq read1.trim.unpaired.fq read2.fq read2.trim.unpaired.fq ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:70
```
### *K*-mer tables building
Use the glistmaker tool from GenomeTester4 to create *k*-mer tables for each of your bulks and genotypes included in your CoSSA. With Illumina 151bp PE reads, we used a *k*-mer size of 31 nucleotides.
```
GenomeTester4/bin/glistmaker read1.trim.fq -w 31 -o sample.all.read1
GenomeTester4/bin/glistmaker read2.trim.fq -w 31 -o sample.all.read2
```
Use the union function of the glistcompare tool from GenomeTester4 to combine the *k*-mer tables from the forward and reverse reads.
```
GenomeTester4/bin/glistcompare sample.all.read1_31.list sample.all.read2_31.list -u -o sample.all
```
Use again glistcompare to remove the *k*-mers with a frequency (~ depth) of 1 as they are most likely due to sequencing errors and would increase the computational time. You need to use the intersection and the cutoff functions of glistcompare.
```
GenomeTester4/bin/glistcompare sample.all_31_union.list sample.all_31_union.list -i -c 2 -o sample.all.cutoff2
```
## CoSSA: BSA using a reference genome
### Set algebra with the *k*-mer sets
Use the set algebra functions available in glistmaker to select for the *k*-mers specific to your haplotype of interest.
We will use as an example the BSA performed by Prodhomme *et al.*, 2019. One resistant bulk (R-bulk), one susceptible bulk (S-bulk), a resistant parent (PR) and a susceptible parent (PS) have been sequences. *K*-mer tables have been built for the four samples.
To select for the R-bulk specific *k*-mers, use the difference function of glistcompare.
```
GenomeTester4/bin/glistcompare R-bulk.all.cutoff2_31_intersect.list S-bulk.all.cutoff2_31_intersect.list -d -o R-bulk-specific-kmers.cutoff2
```
To select for the R-bulk specific *k*-mers with a specific frequency (depending on the expected sequencing depth of your haplotype of interest), use the intersection and the difference functions of glistcompare. For instance for a lower cutoff of 10 and an upper cutoff of 20:
```
GenomeTester4/bin/glistcompare R-bulk-specific-kmers.cutoff2_31_intersect.list R-bulk-specific-kmers.cutoff2_31_intersect.list -i -c 10 -o R-bulk-specific-kmers.cutoff10
GenomeTester4/bin/glistcompare R-bulk-specific-kmers.cutoff2_31_intersect.list R-bulk-specific-kmers.cutoff2_31_intersect.list -i -c 21 -o R-bulk-specific-kmers.cutoff21
GenomeTester4/bin/glistcompare R-bulk-specific-kmers.cutoff10_intersect.list R-bulk-specific-kmers.cutoff21_intersect.list -d -o R-bulk-specific-kmers.cutoff10to20
```
In the example study chosen, the parents of the population have been sequenced. The R-bulk specific *k*-mers can be divided in function of their inheritance: inherited from PR, inherited from PS, inherited from both parents, inherited from none of the parents. This last set allows you to discard the contamination in your data. You need to use the different set algebra functions of glistcompare for this. There are different ways of doing it.
Example:
```
GenomeTester4/bin/glistcompare R-bulk-specific-kmers.cutoff10to20_31_0_diff1.list PR.all.cutoff2_31_ intrsec.list -d -o 1.R-bulk-specific-kmers.cutoff10to20-not-in-PR
GenomeTester4/bin/glistcompare R-bulk-specific-kmers.cutoff10to20_31_0_diff1.list 1.R-bulk-specific-kmers.cutoff10to20-not-in-PR_31_0_diff1.list -d -o 2.R-bulk-specific-kmers.cutoff10to20-in-common-with-PR
GenomeTester4/bin/glistcompare R-bulk-specific-kmers.cutoff10to20_31_0_diff1.list PS.all.cutoff2_31_ intrsec.list -d -o 3.R-bulk-specific-kmers.cutoff10to20-not-in-PS
GenomeTester4/bin/glistcompare R-bulk-specific-kmers.cutoff10to20_31_0_diff1.list 3.R-bulk-specific-kmers.cutoff10to20-not-in-PS_31_0_diff1.list -d -o 4.R-bulk-specific-kmers.cutoff10to20-in-common-with-PS
GenomeTester4/bin/glistcompare 1.R-bulk-specific-kmers.cutoff10to20-not-in-PR31_0_diff1.list 3.R-bulk-specific-kmers.cutoff10to20-not-in-PS_31_0_diff1.list -d -o 5.R-bulk-specific-kmers.cutoff10to20-from-PS
GenomeTester4/bin/glistcompare 3.R-bulk-specific-kmers.cutoff10to20-not-in-PS_31_0_diff1.list 1.R-bulk-specific-kmers.cutoff10to20-not-in-PR_31_0_diff1.list -d -o 6.R-bulk-specific-kmers.cutoff10to20-from-PR
GenomeTester4/bin/glistcompare 3.R-bulk-specific-kmers.cutoff10to20-not-in-PS_31_0_diff1.list 2.R-bulk-specific-kmers.cutoff10to20-in-common-with-PR_31_0_diff1.list -d -o 7.R-bulk-specific-kmers.cutoff10to20-from-none-parents
GenomeTester4/bin/glistcompare 2.R-bulk-specific-kmers.cutoff10to20-in-common-with-PR_31_0_diff1.list 4.R-bulk-specific-kmers.cutoff10to20-in-common-with-PS_31_0_diff1.list -i -o 8.R-bulk-specific-kmers.cutoff10to20-from-both-parents
```
When adding up the unique and total number of *k*-mers of the 4 last sets (5, 6, 7 and 8), you should find back the unique and total number of *k*-mers contained in the R-bulk-specific-kmers.cutoff10to20 original file. 
To check how many *k*-mers are contained in a set, use the glistquery tool of Genometester4.
```
GenomeTester4/bin/glistquery kmer.file.list -stat
```

Transform the binary *k*-mer tables files in text files using glistquery.
```
GenomeTester4/bin/glistquery kmer.file.list > kmer.file.kmer
```

### Mapping the *k*-mers to a reference

Transform the *k*-mer file in a fastq file using the python script kmer_to_fastq.py.
```
python kmer_to_fastq.py kmer.file.kmer > kmer.file.fq
```
Use an alignment software to align the selected *k*-mers to a reference genome. 

Example using [BWA]( http://bio-bwa.sourceforge.net/):
```
bwa index reference.fasta
bwa aln reference.fasta kmer.file.fq > kmer.file.reference.sai
bwa samse reference.fasta kmer.file.reference.sai kmer.file.fq | samtools view -Sbh - | samtools sort -T kmer.file.reference.samse.tmp - > kmer.file.reference.samse.srt.bam
samtools index kmer.file.reference.samse.srt.bam
```
Now the *k*-mers are mapped to the reference, you can count how many *k*-mers mapped to each chromosome bin.
Example for 1Mb bins:
```
bedtools makewindows -b chr01.bed -w 1000000 | bedtools intersect -a - - kmer.file.reference.samse.srt.bam -c -bed > kmer.file.reference.chr01.txt
bedtools makewindows -b chr02.bed -w 1000000 | bedtools intersect -a - - kmer.file.reference.samse.srt.bam -c -bed > kmer.file.reference.chr02.txt
....
```
## CoSSA: BSA without a reference genome
### Extract the PE reads containing the *k*-mers of interest

To extract the PE reads containing at least x *k*-mers from your *k*-mer set, use the shell script extract_reads_from_kmers.sh (it requires the *k*-mers list to be sorted but the output files of GenomeTester4 are already sorted). This shell script reformats the sorted *k*-mer table and extracts the reads containing at least x *k*-mers from the list using the kfastqfilter tool.
Example:
```
qsub extract_reads_from_kmers.sh kmer.file.sorted.kmer
```
*De novo* assemble these PE reads. 
For example using [SPAdes]( http://cab.spbu.ru/software/spades/):
```
python SPAdes-3.11.1/bin/spades.py -t 15 -k 21,33,55,77 --careful --pe1-1 extracted_read1.fq --pe1-2 extracted_read2.fq -o denovo.assembly
```
To identify the haplotype specific SNPs in your *de novo* assembly, map back the *k*-mers of interest to the assembly (for example using BWA-ALN).
## CoSSA: eliminate common SNPs
If your aim is to design haplotype specific markers associated with your trait of interest, it may be useful to remove common SNPs from your analysis which might not be diagnostic in a broader panel of genotypes.
After building the *k*-mer tables for genotypes which do not hold the trait of interest, make the union of their *k*-mers using glistcompare and subtract them from your selected *k*-mers using the difference function of glistcompare.
 ```
GenomeTester4/bin/glistcompare genotype1.all.cutoff2_31_ intrsec.list genotype2.all.cutoff2_31_ intrsec.list -u -o genotype1.genotype2.all.cutoff2
GenomeTester4/bin/glistcompare kmer.file.list genotype1.genotype2.all.cutoff2_31_union.list -d -o kmer.file.difference.genotype1.genotype2
```
## CoSSA: pedigree analysis
To inquire if two varieties sharing the same trait hold the same causal locus, use the intersection tool of glistcompare.
```
GenomeTester4/bin/glistcompare kmer.file.list genotype3 all.cutoff2_31_ intrsec.list -i -o kmer.file.inter.genotype3
```
Then map the output *k*-mers to the reference genome or the *de novo* assembly.

## Tips
### Manual validation
To design markers using the CoSSA outputs, load the bam file of the *k*-mers alignment to a reference genome or to a *de novo* assembly in a Genome Browser, e.g. [IGV]( http://software.broadinstitute.org/software/igv/) and look for stacks of 31 *k*-mers. They pinpoint SNPs without flanking SNPs. 
