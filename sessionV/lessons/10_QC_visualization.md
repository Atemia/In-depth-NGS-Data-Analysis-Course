---
title: "QC visualization of peaks with IGV"
author: "Meeta Mistry"
date: "November 15, 2017"
---

Approximate time: 45 minutes

## Learning Objectives
* Generate bigWig files
* Use IGV to visualize BigWig, BED and data from ENCODE

## Visualization of ChIP-seq data

* We are going to generate some bigWig files to visualize our data using IGV.
  Instead of using BAM files that can be large and cannot be normalized, we
  generate bigWig files that have been normalized for read count relative to the
  input control. We will be using the `bamCompare` command within
  [deepTools](https://deeptools.github.io/) for this: [instructions to generate
  normalized bigwigs with
  deepTools](https://github.com/fidelram/deepTools/wiki/Normalizations).

> "**The bigWig format** is useful for dense, continuous data that will be
> displayed in the Genome Browser as a graph. BigWig files are created from
> [wiggle (wig)](https://genome.ucsc.edu/goldenpath/help/wiggle.html) type files
> using the program wigToBigWig.
>
> The bigWig files are in an indexed binary format. The main advantage of this
> format is that only those portions of the file needed to display a particular
> region are transferred to the Genome Browser server. Because of this, bigWig
> files have considerably faster display performance than regular wiggle files
> when working with large data sets. The bigWig file remains on your local
> web-accessible server (http, https or ftp), not on the UCSC server, and only
> the portion needed for the currently displayed chromosomal position is locally
> cached as a "sparse file"."
>
> -[https://genome.ucsc.edu/goldenpath/help/bigWig.html](https://genome.ucsc.edu/goldenpath/help/bigWig.html)

If you have logged out of your prior interactive session, restart one now:

```bash
$ srun -n 2 --mem 2000 -p classroom --pty bash
````

Switch back into the work directory, and create a new directory for generating the bigWig files:

```bash
$ cd ~/ngs_course/chipseq/results/

$ mkdir visualization
```

Load in deeptools:

```bash
$ module load deepTools/2.5.3-IGB-gcc-4.9.4-Python-2.7.13
```

We'll use the `bamCompare` script to generate some comparative tracks.

> For comprehensive help the deeptools
> [`bamCompare`](https://deeptools.readthedocs.io/en/latest/content/tools/bamCompare.html)
> manual page is actually quite detailed and gives common example use cases

Now, let's run a few quantitative tracks that compare the inputs and IP for the Nanog data:

```bash
$ bamCompare -b1 bowtie2/H1hesc_Nanog_Rep1_chr12_aln.bam \
  -b2 bowtie2/H1hesc_Input_Rep1_chr12_aln.bam \
  -o visualization/Nanog_Rep1_chr12_compare.bw --verbose 2> visualization/Nanog_Rep1_bamcompare.log

$ bamCompare -b1 bowtie2/H1hesc_Nanog_Rep2_chr12_aln.bam \
  -b2 bowtie2/H1hesc_Input_Rep2_chr12_aln.bam \
  -o visualization/Nanog_Rep2_chr12_compare.bw --verbose 2> visualization/Nanog_Rep2_bamcompare.log
```

Let's do the same for Pou5f1.  Use the above to do this on your own first (no peeking)!

<details>
<p>

```bash
$ bamCompare -b1 bowtie2/H1hesc_Pou5f1_Rep1_chr12_aln.bam \
  -b2 bowtie2/H1hesc_Input_Rep1_chr12_aln.bam \
  -o visualization/Pou5f1_Rep1_chr12_compare.bw --verbose 2> visualization/Pou5f1_Rep1_bamcompare.log

$ bamCompare -b1 bowtie2/H1hesc_Pou5f1_Rep2_chr12_aln.bam \
  -b2 bowtie2/H1hesc_Input_Rep2_chr12_aln.bam \
  -o visualization/Pou5f1_Rep2_chr12_compare.bw --verbose 2> visualization/Pou5f1_Rep2_bamcompare.log
```

</p>
</details>
<br><br>

Now, let's generate simple coverage data for all of our files. Here we
will use a `for` loop, and will use the `bamCoverage` tool (which works on
single files). We'll keep it simple and not perform any normalization since we
primarily are concerned with qualitatively looking at coverage across samples,
but it's worth keeping in mind you can normalize coverage using `bamCoverage` if
you want a more quantitative comparison.

> The deeptools help page for 
> [`bamCoverage`](https://deeptools.readthedocs.io/en/latest/content/tools/bamCoverage.html)
> also gives more information about options and example use cases

```bash
$ for bam in bowtie2/*aln.bam
do
bw=`basename $bam _aln.bam`
bamCoverage --binSize 10 -b $bam -o visualization/${bw}.bw -v > visualization/${bw}_bamcoverage.log
done
```

Copy over the bigWig files using Cyberduck or MobaXTerm; place them in the
`chipseq-project` folder in a new directory under `data`, called `bigWigs`.
Also, copy over the BEDtools overlap/intersect files to your computer and place
them in a separate folder under `data` called `intersects`.

## IGV

Start IGV, choose 'Human hg19' as the genome. Load the 2 Rep1 bamCompare bigWig
files, as well as the overlap BED files, and zoom in to chr12 where the coverage is.

You will notice that there are positive and negative values on the track, what
do you think this denotes in the context of normalization?

> You can generate a simple, non-normalized bigWig with `bamCoverage` and you
> won't see any negative values.

Now, load in the Rep1 bigWig bamCoverage files. These all have simple positive
values denoting read coverage for those regions. Note the input.

You can also load in other files, including BAM and narrowPeak files.

Let's open just two BAMs (Nanog Rep1, and Input Rep1). You could hunt and peck
read coverage; let's look at a pretty significant one reported for MACS2. In the
location box type (or paste) 'chr12:11761770-11762413', which should be the top
Nanog Rep1 peak call based on qvalue. Note the bimodal read distribution, and
notice there are overlapping calls in the BED data (including Pou5f1).

Now load in the Pou5f1 Rep1 BAM file. Note the enrichment again, but the
difference in peak profile.

Finally, we are going to visually compare our output to the output from the full
dataset from ENCODE, by loading that data from the IGV server.

***
*This lesson has been developed by members of the teaching team at the [Harvard
Chan Bioinformatics Core (HBC)](http://bioinformatics.sph.harvard.edu/). These
are open access materials distributed under the terms of the [Creative Commons
Attribution license](https://creativecommons.org/licenses/by/4.0/) (CC BY 4.0),
which permits unrestricted use, distribution, and reproduction in any medium,
provided the original author and source are credited.*
