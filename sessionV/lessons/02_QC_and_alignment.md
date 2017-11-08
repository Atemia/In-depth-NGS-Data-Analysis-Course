---
title: "FastQC and alignment"
author: "Mary Piper, Radhika Khetani (HSPH), Chris Fields (UIUC)"
date: "November 7, 2017"
output: 
  revealjs::revealjs_presentation:
    theme: solarized
    highlight: pygments
    transition: slide
    self_contained: true
    slide_level: 1
    css: styles.css
    reveal_options:
      slideNumber: true
---

Contributors: Mary Piper, Radhika Khetani, Chris Fields

Approximate time: 60 minutes

---

# Learning Objectives

* to use previous knowledge of quality control steps to perform FastQC
* to learn how to use Trimmomatic to perform quality trimming
* to understand parameters and perform alignment using Bowtie2

# Quality control of sequence reads

<center><img src="../img/chip_workflow_june2017_step1.png" width=400></center>

* Now that we have our files and directory structure, we are ready to begin our ChIP-Seq analysis. 
* For any NGS analysis method, our first step in the workflow is to explore the quality of our reads prior to aligning them to the reference genome and proceeding with downstream analyses. 

---

We will use FastQC to get a good idea of the overall quality of our data. We will use FastQC to identify whether any samples appear to be outliers, to examine our data for contamination, and to determine a trimming strategy.

>**NOTE:** We will trim poor quality bases and/or adapters prior to alignment because that was the workflow previously used by ENCODE. However, we do not need to trim as the downstream alignment tool, Bowtie2, has an option for soft-clipping.

---

# FASTQC

Let's run FastQC on all of our files. 

---

**Start an interactive session with 2 cores** if you don't have one going, and change directories to the `raw_data` folder.

```bash
$ cd ~/ngs_course/chipseq/raw_data 
```

---

First, what does the data look like?

```bash
$ head -n 8 H1hesc_Input_Rep1_chr12.fastq
```

---

This is [FASTQ format](https://en.wikipedia.org/wiki/FASTQ_format).

```
[instru03@compute-1-0 raw_data]$ head -n 8 H1hesc_Input_Rep1_chr12.fastq

@ILLUMINA-EAS295:72:70H93AAXX:3:17:5807:13712
CAATTGAGAATATAATTGCCTTGAAAATAAAAATGT
+
DGGEDEGEDDEGDD=GGGDG;BB;BGDGDEEDDDBG
@ILLUMINA-EAS295:72:70H93AAXX:3:18:18252:17693
TTTATTTTCAAGGCAATTATATTCTCAATTGGCTCT
+
IIIIIIIIIIFHIIIIIIIIIIIIIIDIIIIDIIID
```

This contains the actual base calls for each reads from the sequencing run.

---

Next, let's run FASTQC on this:

```
# Hint: use tab here
$ module load FastQC/0.11.5-IGB-gcc-4.9.4-Java-1.8.0_121

$ fastqc H1hesc_Input_Rep1_chr12.fastq 
```

---

What was generated?

```
[instru03@compute-1-0 raw_data]$ ls -l

total 191566
-rw-rw-r-- 1 instru03 instru03 59738626 Nov  6 15:55 H1hesc_Input_Rep1_chr12.fastq
-rw-rw-r-- 1 instru03 instru03   242724 Nov  6 16:08 H1hesc_Input_Rep1_chr12_fastqc.html
-rw-rw-r-- 1 instru03 instru03   292563 Nov  6 16:08 H1hesc_Input_Rep1_chr12_fastqc.zip
-rw-rw-r-- 1 instru03 instru03 36809331 Nov  6 15:55 H1hesc_Input_Rep2_chr12.fastq
-rw-rw-r-- 1 instru03 instru03 16481860 Nov  6 15:55 H1hesc_Nanog_Rep1_chr12.fastq
-rw-rw-r-- 1 instru03 instru03 32164197 Nov  6 15:55 H1hesc_Nanog_Rep2_chr12.fastq
-rw-rw-r-- 1 instru03 instru03 16895456 Nov  6 15:55 H1hesc_Pou5f1_Rep1_chr12.fastq
-rw-rw-r-- 1 instru03 instru03 33535661 Nov  6 15:55 H1hesc_Pou5f1_Rep2_chr12.fastq
```

---

Now, move all of the `fastqc` files to the `results/untrimmed_fastqc` directory:

```bash
$ mv *fastqc* ../results/untrimmed_fastqc/
```

---

Transfer the FastQC zip file for Input replicate 1 to your local machine using [your favorite SFTP transfer tool](https://help.igb.illinois.edu/File_Server_Access) and view the report.

Here, we will try [Cyberduck](https://help.igb.illinois.edu/File_Server_Access#Connect_From_OSX_Using_CyberDuck_.28Very_Secure.29).  In the instructions, instead of using `file-server.igb.illinois.edu`, we will enter `biologin.igb.illinois.edu`.  

---

<center><img src="../img/fastqc_input_rep1.png"></center>

Based on the sequence quality plot, we see across the length of the read the quality drops into the low range. Trimming should be performed from both ends of the sequences. 

---

# Question

What else is in the report?

---

# Trimmomatic

[*Trimmomatic*](http://www.usadellab.org/cms/?page=trimmomatic) can be used to trim away adapters and filter out poor quality score reads. *Trimmomatic* is a Java-based program that can remove sequencer specific reads and nucleotides that fall below a certain quality threshold. *Trimmomatic* offers the option to trim reads using a hard crop, sliding window or base-by-base methods. It can also trim adapter sequences and remove reads if below a minimum length. In addition, *Trimmomatic* can be multi-threaded to run quickly using a single, complex command. 

---

We will use Trimmomatic to trim the reads from both ends of the sequence.

So, let's do a quick refresher on modules.  First, let's check for the *Trimmomatic* module and load it:

```bash
# this will run a search for partial matches
$ module avail trim
```

---

Now let's load the module:

```
$ module load Trimmomatic/0.36-Java-1.8.0_121
To execute Trimmomatic run: java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.36.jar
```

Note how it says to run the tool.

---
  
What modules are loaded:

```
$ module list
```

You can also see how modules make our (command-line) lives a little easier:

```
$ echo $PATH
```

---

Okay, back to work.

By loading the *Trimmomatic* module, Java and Trimmomatic are loaded and appear in our PATH. 

We will run Trimmomatic using the following parameters:

* `SE`: Single End reads
* `-threads`: number of threads / cores
* `-phred33`: quality score format
* `LEADING`: cut bases off the start of a read, if below a threshold quality
* `TRAILING`: cut bases off the end of a read, if below a threshold quality
* `MINLEN`: drop an entire read if it is below a specified length

<small>
*NOTE:* We have to specify the `-threads` parameter because *Trimmomatic* uses all threads on a node by default.

</small>
---

*Trimmomatic* has a variety of other options and parameters:

* **_SLIDINGWINDOW_** Perform sliding window trimming, cutting once the average quality within the window falls below a threshold.
* **_CROP_** Cut the read to a specified length.
* **_HEADCROP_** Cut the specified number of bases from the start of the read.
* **_ILLUMINACLIP_** Cut adapter and other illumina-specific sequences from the read
* **_TOPHRED33_** Convert quality scores to Phred-33.
* **_TOPHRED64_** Convert quality scores to Phred-64.

---

Now that we know what parameters  we can set up our command. Since we are only trimming a single file, we will run the command in the interactive session rather than creating a script. 

Because *Trimmomatic* is Java based, it is run using the `java -jar` command. In addition to the options as described above, we have two arguments specifying our input file and output file names. 
<small>

>*NOTE:* `java -jar` calls the Java program, which is needed to run *Trimmomatic*, which is a 'jar' file (`trimmomatic-0.33.jar`). A 'jar' file is a special kind of java archive that is often used for programs written in the Java programming language.  If you see a new program that ends in '.jar', you will know it is a java program that is executed `java -jar` <*location of program .jar file*>. Even though *Trimmomatic* is in our PATH, we still need to specify the full path to the `.jar` file in the command.

</small>

---

Remember what the module prompt said about how to run Trimmomatic?

```bash
$ java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.36.jar SE \
  -threads $SLURM_NTASKS \
  -phred33 \
  H1hesc_Input_Rep1_chr12.fastq \
  ../results/trimmed/H1hesc_Input_Rep1_chr12.qualtrim20.minlen36.fq \
  LEADING:20 \
  TRAILING:20 \
  MINLEN:36
```

Note the line where we are pushing data to the `results/trimmed` folder.

---

You should see something like this:

```text
TrimmomaticSE: Started with arguments:
 -threads 1 -phred33 H1hesc_Input_Rep1_chr12.fastq ../results/trimmed/H1hesc_Input_Rep1_chr12.qualtrim20.minlen36.fq LEADING:20 TRAILING:20 MINLEN:36
Input Reads: 489620 Surviving: 418480 (85.47%) Dropped: 71140 (14.53%)
TrimmomaticSE: Completed successfully
```

---

Let's see how much trimming improved our reads by running FastQC again:

```bash
$ fastqc ../results/trimmed/H1hesc_Input_Rep1_chr12.qualtrim20.minlen36.fq
```

---

Move the FastQC folders to the results directory for trimmed FastQC results:

```bash
$ mv ../results/trimmed/*fastqc* ../results/trimmed_fastqc/
```

---

Using Cyberduck, transfer the file for the trimmed Input replicate 1 FastQC to your computer.

<center><img src="../img/chipseq_trimmed_fastqc.png"></center>

# Alignment

Now that we have removed the poor quality sequences from our data, we are ready to align the reads to the reference genome. 

Most ChIP-seq experiments do not require gapped alignments because the sequenced reads do not contain them, unlike exon junctions in RNA-seq analyses 

Therefore, we do not need a splice-aware aligner. We can use a traditional short-read aligner to quickly and accurately align reads to the genome.

---

[Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml) is a fast and accurate alignment tool that indexes the genome with an FM Index based on the Burrows-Wheeler Transform to keep memory requirements low for the alignment process. 

By default, Bowtie2 will perform a global end-to-end read alignment, which is best for quality-trimmed reads. However, it also has a local alignment mode, which will perform soft-clipping for the removal of poor quality bases or adapters from untrimmed reads.

---

*Bowtie2* supports gapped, local and paired-end alignment modes and works best for reads that are at least 50 bp (shorter read lengths should use Bowtie1). 

> _**NOTE:** Our reads are only 36 bp, so technically we should explore alignment Bowtie1 to see if it is better. However, since it is rare that you will have sequencing reads with less than 50 bp, we will show you how to perform alignment using Bowtie2._

# Question: Finding Bowtie2

How would you look for this on the cluster?

# Creating Bowtie2 index

To perform the Bowtie2 alignment, a genome index is required. **We previously generated the genome indices for you**, and they exist in the `reference_data` directory.

However, if you needed to create a genome index yourself, you would use the following command:

```bash
# DO NOT RUN

bowtie2-build <path_to_reference_genome.fa> <prefix_to_name_indexes>

# Though we don't always recommend it's use, you can find pre-built indexes for 
# the entire human genome (and other genomes) on biocluster using following path: 
# /home/mirror/igenome/Homo_sapiens/UCSC/hg19/Sequence/Bowtie2Index/
```

---

It is worth quickly checking the format for the reference.  Here we are using the chromosome 12 sequence from the UCSC hg19 (human genome release 37, UCSC release 19).

```
$ head -n 10 
```
---

This is [FASTA format](https://en.wikipedia.org/wiki/FASTA_format).

```
[instru03@compute-1-0 raw_data]$ head -n 10 ../reference_data/chr12.fa
>chr12
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
```

---

# Aligning reads with Bowtie2

Since we have our indices already created, we can get started with read alignment. Change directories to the `bowtie2` folder:

```bash
$ cd ~/ngs_course/chipseq/results/bowtie2
```

---

We will perform alignment on our single trimmed sample, `H1hesc_Input_Rep1_chr12.qualtrim20.minlen36.fq`. 

Details on Bowtie2 and its functionality can be found in the [user manual](http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml); we encourage you to peruse through to get familiar with all available options.

---

The basic options for aligning reads to the genome using Bowtie2 are:

* `-p`: number of processors / cores
* `-q`: reads are in FASTQ format
* `-x`: /path/to/genome_indices_directory
* `-U`: /path/to/FASTQ_file
* `-S`: /path/to/output/SAM_file

---

```bash
$ module load Bowtie2/2.3.2-IGB-gcc-4.9.4

$ bowtie2 -p $SLURM_NTASKS -q \
-x ~/ngs_course/chipseq/reference_data/chr12 \
-U ~/ngs_course/chipseq/results/trimmed/H1hesc_Input_Rep1_chr12.qualtrim20.minlen36.fq \
-S ~/ngs_course/chipseq/results/bowtie2/H1hesc_Input_Rep1_chr12_aln_unsorted.sam

```
> **NOTE:** If you had untrimmed fastq files, you would want to use local alignment to perform soft-clipping by including the option `--local`.

---

You should see output like this:

```text
418480 reads; of these:
  418480 (100.00%) were unpaired; of these:
    45644 (10.91%) aligned 0 times
    302306 (72.24%) aligned exactly 1 time
    70530 (16.85%) aligned >1 times
89.09% overall alignment rate
```

---

What does the alignment file look like?

```
$ head -n 20 ~/ngs_course/chipseq/results/bowtie2/H1hesc_Input_Rep1_chr12_aln_unsorted.sam
```

This is [SAM format](https://en.wikipedia.org/wiki/SAM_(file_format)).

```text
@HD     VN:1.0  SO:unsorted
@SQ     SN:chr12        LN:133851895
@PG     ID:bowtie2      PN:bowtie2      VN:2.3.2        CL:"/home/apps/software/Bowtie2/2.3.2-IGB-gcc-4.9.4/bin/bowtie2-align-s --wrapper basic-0 -p 1 -q -x /home/a-m/instru03/ngs_course/chipseq/reference_data/chr12 -S /home/a-m/instru03/ngs_course/chipseq/results/bowtie2/H1hesc_Input_Rep1_chr12_aln_unsorted.sam -U /home/a-m/instru03/ngs_course/chipseq/results/trimmed/H1hesc_Input_Rep1_chr12.qualtrim20.minlen36.fq"
ILLUMINA-EAS295:72:70H93AAXX:3:17:5807:13712    16      chr12   1000084 42      36M     *       0       0       ACATTTTTATTTTCAAGGCAATTATATTCTCAATTG    GBDDDEEDGDGB;BB;GDGGG=DDGEDDEGEDEGGD    AS:i:0  XN:i:0  XM:i:0  XO:i:0  XG:i:0  NM:i:0  MD:Z:36 YT:Z:UU
ILLUMINA-EAS295:72:70H93AAXX:3:18:18252:17693   0       chr12   1000089 42      36M     *       0       0       TTTATTTTCAAGGCAATTATATTCTCAATTGGCTCT    IIIIIIIIIIFHIIIIIIIIIIIIIIDIIIIDIIID    AS:i:0  XN:i:0  XM:i:0  XO:i:0  XG:i:0  NM:i:0  MD:Z:36 YT:Z:UU
...
```

These text files can be very large, about double the size from the trimmed FASTQ file.  These can be compressed though...

---

# Filtering reads

An important issue concerns the inclusion of multiple mapped reads (reads mapped to multiple loci on the reference genome). 

**Allowing for multiple mapped reads increases the number of usable reads and the sensitivity of peak detection; however, the number of false positives may also increase** [[1]](https://www.ncbi.nlm.nih.gov/pubmed/21779159/). 

Therefore we need to filter our alignment files to **contain only uniquely mapping reads** in order to increase confidence in site discovery and improve reproducibility. 

---

Since there is no parameter in Bowtie2 to keep only uniquely mapping reads, we will need to perform the following steps to generate alignment files containing only the uniquely mapping reads:

>1. Change alignment file format from SAM to BAM
>2. Sort BAM file by read coordinate locations
>3. Filter to keep only uniquely mapping reads (this will also remove any unmapped reads)

# 1. Changing file format from SAM to BAM

While the SAM alignment file output by Bowtie2 is human readable, we need a BAM alignment file for downstream tools. Therefore, we will use [Samtools](http://samtools.github.io) to convert the file formats. The command we will use is `samtools view` with the following parameters:

* `-h`: include header in output
* `-S`: input is in SAM format
* `-b`: output BAM format
* `-o`: /path/to/output/file

---

```bash
$ module load SAMtools/1.5-IGB-gcc-4.9.4

$ samtools view -h -S -b \
-o H1hesc_Input_Rep1_chr12_aln_unsorted.bam \
H1hesc_Input_Rep1_chr12_aln_unsorted.sam
```

The output is a [BAM file](https://en.wikipedia.org/wiki/Binary_Alignment_Map):

```bash
$ ls -lh 
total 101797
-rw-rw-r-- 1 instru03 instru03 19255572 Nov  6 22:26 H1hesc_Input_Rep1_chr12_aln_unsorted.bam
-rw-rw-r-- 1 instru03 instru03 84984011 Nov  6 22:16 H1hesc_Input_Rep1_chr12_aln_unsorted.sam
```

Note the BAM file size is significantly smaller than the original SAM output; BAM is a binary compressed format.  It is **not** human readable, though; this format is optimized for compression and data analysis.

---

You can find additional parameters for the samtools functions in the [manual](http://www.htslib.org/doc/samtools-1.2.html).

# 2. Sorting BAM files by genomic coordinates

Before we can filter to keep the uniquely mapping reads, we need to sort our BAM alignment files by genomic coordinates. To perform this sort, we will use [Sambamba](http://lomereiter.github.io/sambamba/index.html), which is a tool that quickly processes BAM and SAM files. It is similar to SAMtools, but has unique functionality.
The command we will use is `sambamba sort` with the following parameters:

* `-t`: number of threads / cores
* `-o`: /path/to/output/file

---

```bash
$ module load sambamba/0.6.6
$ sambamba sort -t $SLURM_NTASKS \
-o H1hesc_Input_Rep1_chr12_aln_sorted.bam \
H1hesc_Input_Rep1_chr12_aln_unsorted.bam 
```

---

Note the outputs:

```bash
[instru03@compute-1-0 bowtie2]$ ls -l
total 120426
-rw-rw-r-- 1 instru03 instru03 18821020 Nov  6 22:31 H1hesc_Input_Rep1_chr12_aln_sorted.bam
-rw-rw-r-- 1 instru03 instru03   254128 Nov  6 22:31 H1hesc_Input_Rep1_chr12_aln_sorted.bam.bai
-rw-rw-r-- 1 instru03 instru03 19255572 Nov  6 22:26 H1hesc_Input_Rep1_chr12_aln_unsorted.bam
-rw-rw-r-- 1 instru03 instru03 84984011 Nov  6 22:16 H1hesc_Input_Rep1_chr12_aln_unsorted.sam
```

There is a `.bai` file that is generated; this is a BAM index file, which is very useful for processing sortd alignment data.  `sambamba` generates this for you automatically when you sort a BAM file; if you use `samtools` you would need to do this in a second step.

# 3. Filtering uniquely mapping reads

Finally, we can filter the uniquely mapped reads. We will use the `sambamba view` command with the following parameters:

* `-t`: number of threads / cores
* `-h`: print SAM header before reads
* `-f`: format of output file (default is SAM)
* `-F`: set [custom filter](https://github.com/lomereiter/sambamba/wiki/%5Bsambamba-view%5D-Filter-expression-syntax) - we will be using the filter to remove multimappers and unmapped reads.

---

```bash
$ sambamba view -h -t $SLURM_NTASKS -f bam \
-F "[XS] == null and not unmapped " \
H1hesc_Input_Rep1_chr12_aln_sorted.bam > H1hesc_Input_Rep1_chr12_aln.bam
```

We filtered out unmapped reads by specifying in the filter `not unmapped`. 

Also, among the reads that were aligned, we filtered out multimappers by specifying `[XS] == null`. `XS` is a tag generated by Bowtie2 that gives an alignment score for the second-best alignment, and it is only present if the read is aligned and more than one alignment was found for the read.

---

Now that the alignment files contain only uniquely mapping reads, we are ready to perform peak calling.

> _**NOTE:** After performing read alignment, it's often useful to generate QC metrics for the alignment using tools such as [QualiMap](http://qualimap.bioinfo.cipf.es/doc_html/analysis.html#bam-qc) or [MultiQC](http://multiqc.info) prior to moving on to the next steps of the analysis._ 

***
*This lesson has been developed by members of the teaching team at the [Harvard Chan Bioinformatics Core (HBC)](http://bioinformatics.sph.harvard.edu/). These are open access materials distributed under the terms of the [Creative Commons Attribution license](https://creativecommons.org/licenses/by/4.0/) (CC BY 4.0), which permits unrestricted use, distribution, and reproduction in any medium, provided the original author and source are credited.*

