---
title: "Handling replicates"
author: "Meeta Mistry"
date: "Monday, June 27th, 2017"
---

Contributors: Meeta Mistry, Radhika Khetani, Chris Fields

Approximate time: 75 minutes

## Learning Objectives

* Combining replicates using simple overlap with Bedtools
* Combining replicates to only get the highly reproducible peaks using the IDR
  method

## Handling replicates in ChIP-Seq

As with any high-throughput experiment, a single assay is often subject to a
substantial amount of variability. Thus, it is highly recommended to setup your
experimental design with a minimum of 2-3 biological replicates. Presumably, two
replicates measuring the same underlying biology should have high consistency
but that is not always the case. In order to evaluate consistency between
replicates **we require metrics that objectively assess the reproducibility of
high-throughput assays**.

<img src="../img/chip_workflow_june2017_step4.png" width=700>

Since we have 2 replicates in this example, we want to consider only those peaks
that are present in both replicates before we compare the peaks from the two
transcription factors to one another.

<img src=../img/idr_samples.png width=500>

Common methods for handling replicates includes taking overlapping peak calls
across replicates and then assessing differences in binding regions.
Additionally, there are more complex methods that employ statistical testing and
evaluate the reproducibility between replicates. In this lesson we will cover
both methods.

> **NOTE: A recent talk on "Accessing and using ENCODE data" [linked here](https://hstalks.com/t/2679/accessing-and-using-encode-data/?biosci)
> where they talk about handling replicates and the similarities and differences
> when using an overlap versus IDR analysis.**

## Overlapping peaks

In this section, our goal is to determine what peaks are in common between the
the two replicates for each factor (Nanog and Pou5f1). To perform this task we
are going to use a suite of tools called `bedtools`.

### `bedtools`

The idea is that genome coordinate information can be used to perform relatively
simple arithmetic, like combining, subsetting, intersecting, etc., to obtain all
sorts of information.
[bedtools](http://bedtools.readthedocs.org/en/latest/index.html) from [Aaron
Quinlan's group](http://quinlanlab.org/) at University of Utah is easy to use,
and an extremely versatile tool that performs tasks of this nature.

<img src="../img/bedtools.png" width="700">

As the name implies, this suite of tools works with [**BED** files](https://genome.ucsc.edu/FAQ/FAQformat.html#format1), but it also
works with other file formats that have genome coordinate information.

<img src="../img/bedtools-basic.png" width="600">

> **NOTE:** When working with multiple files to perform arithmetic on genomic
> coordinates, it is essential that all files have coordinate information for
> the same exact version of the genome and the same coordinate system (0-based
> or 1-based)!

### Setting up

Let's start an interactive session and change directories and set up a space for
the resulting overlaps.

```bash
$ srun -n 2 --mem 2000 -p classroom --pty bash

$ cd ~/ngs_course/chipseq/results/

$ mkdir bedtools
```

Load the modules for `bedtools` and `samtools`:

```bash
$ module load BEDTools/2.26.0-IGB-gcc-4.9.4

$ module load SAMtools/1.5-IGB-gcc-4.9.4
```

### Finding overlapping peaks between replicates

The `bedtools intersect` command only reports back the peaks that are
overlapping with respect to the file defined as `a` in the command.

<img src="../img/bedtools_intersect.png" width="600">

To find out more information on the parameters available when intersecting, use
the help flag:

```bash
$ bedtools intersect -h
```

The intersect tool evaluates A (file 1) and finds regions that overlap in B
(file 2). We will add the `-wo` which indicates to write the original A (file 1)
and B (file 2) entries plus the number of base pairs of overlap between the two
features.

Let's start with the Nanog replicates:

```bash
$ bedtools intersect \
-a macs2/Nanog-rep1_peaks.narrowPeak \
-b macs2/Nanog-rep2_peaks.narrowPeak \
-wo > bedtools/Nanog-overlaps.bed
```


#### Next steps:

1. **How many overlapping peaks did we get for the Nanog samples?** (hint hidden below)

<details>
  <p>

```bash
$ wc -l bedtools/Nanog-overlaps.bed
```

</p>
</details>
<br>

2. Do the same for the Pou5f1 replicates.  Try this on your own first.  No peeking!

<details>
  <p>

```bash
$ bedtools intersect \
-a macs2/Pou5f1-rep1_peaks.narrowPeak \
-b macs2/Pou5f1-rep2_peaks.narrowPeak \
-wo > bedtools/Pou5f1-overlaps.bed

$ wc -l bedtools/Pou5f1-overlaps.bed
```

</p>
</details>
<br><br>

Note that we are working with subsetted data and so our list of peaks for each
replicate is small. Thus, the overlapping peak set will be small as we found
with both Nanog and Pou5f1. What is interesting though, is that even though the
individual peak lists are smaller for Pou5f1 samples, the overlapping replicates
represent a higher proportion of overlap with respect to each replicate.

> **_Historical Note_:** A simpler heuristic for establishing reproducibility
> was previously used as a standard for depositing ENCODE data and was in effect
> when much of the currently available data was submitted. According to this
> standard, either 80% of the top 40% of the peaks identified from one replicate
> using an acceptable scoring method should overlap the list of peaks from the
> other replicate, OR peak lists scored using all available reads from each
> replicate should share more than 75% of regions in common. As with the current
> standards, this was developed based on experience with accumulated ENCODE
> ChIP-seq data, albeit with a much smaller sample size.

## Irreproducibility Discovery Rate (IDR)

[IDR](https://github.com/kundajelab/idr) is a framework developed by Qunhua Li
and Peter Bickel's group that **compares a pair of ranked lists of regions/peaks
and assigns values that reflect its reproducibility.**

<img src=../img/idr_figure.png>

It is extensively used by the ENCODE and modENCODE projects and is part of their
[ChIP-seq guidelines and
standards](http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3431496/). It has been
established for submission of ChIP-seq data sets and have been constructed based
on the historical experiences of ENCODE ChIP-seq data production groups.

### Why IDR?

* IDR avoids choices of initial cutoffs, which are not comparable for different
  callers
* IDR does not depend on arbitrary thresholds and so all regions/peaks are
  considered.
* It is based on ranks, so does not require the input signals to be calibrated
  or with a specific fixed scale (only order matters).

### Components of IDR

The IDR approach creates a curve, from which it then quantitatively assesses
when the ﬁndings are no longer consistent across replicates. There are three
main components:

1) A **correspondence curve**: a graphical representation of matched peaks as
you go down the ranked list. Qualitative, not adequate for selecting signals.

<img src=../img/corr_curve.png width=400>

2) An **inference procedure**: summarizes the proportion of reproducible and
irreproducible signals. Quantitative, using a [copula mixture model](https://en.wikipedia.org/wiki/Copula_(probability_theory)).

> What proportion of identifications have a poor correspondence, i.e. falling into ”noise”?
> How consistent are the identifications before reaching breakdown?

3) **Irreproducible Discovery Rate (IDR)**: Derive a significance value from the
inference procedure (#2) in a fashion similar to FDR, and can be used to control
the level of irreproducibility rate when selecting signals. i.e. 0.05 IDR means
that peak has a 5% chance of being an irreproducible discovery

### The IDR pipeline

There are three main steps to the IDR pipeline:

1. Evaluate peak consistency between **true replicates**
2. Evaluate peak consistency between **pooled pseudo-replicates**
3. Evaluate **self-consistency** for each individual replicate

<img src=../img/idr_pipeline.png>

> This figure is taken from the [ENCODE ChIP-Seq Guidelines](http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3431496/#box3).

_We will only be running Step 1 in this lesson, but will discuss steps 2 and 3
in a bit more detail._

## Running IDR

To run IDR *we should be using the full dataset*, and if you are interested the
full BAM files can be downloaded from ENCODE, however for timeliness and
consistency we will continue with our subsetted data for this lesson.

To run IDR the recommendation is to **run MACS2 less stringently** to allow a
larger set of peaks to be identified for each replicate. In addition the
narrowPeak files have to be **sorted by the `-log10(p-value)` column**. We have
already performed these 2 steps for the samples from both groups, and the
commands we used for it as listed below. To reiterate, _you do NOT NEED TO RUN
the following lines of code, we have already generated narrowPeak files for
you!_

```bash
###DO NOT RUN THIS###

# Call peaks using a liberal p-value cutoff
macs2 callpeak -t treatFile.bam -c inputFile.bam -f BAM -g 1.3e+8 -n macs/NAME_FOR_OUPUT -B -p 1e-3  2> macs/NAME_FOR_OUTPUT_macs2.log

#Sort peak by -log10(p-value)
sort -k8,8nr NAME_OF_INPUT_peaks.narrowPeak > macs/NAME_FOR_OUPUT_peaks.narrowPeak
```

> IDR will work with many different Peak callers, the following have been tested:
>
> * SPP - Works out of the box
> * MACS1.4 - DO NOT use with IDR
> * MACS2 - Works well with IDR with occasional problems of too many ties in
>   ranks for low quality ChIP-seq data.
> * HOMER - developers have a detailed pipeline and code (in beta) for IDR
>   analysis with HOMER at https://github.com/karmel/homer-idr
> * PeakSeq - Run with modified PeakSeq parameters to obtain large number of
>   peaks
> * HotSpot, MOSAiCS, GPS/GEM, …

### Setting up

IDR is an open source tool available on [GitHub](https://github.com/kundajelab/idr).
It is a Python program that has already been installed on Orchestra. The first
thing we need to do is load the module to run IDR:

```bash

$ module load IDR/2.0.4-IGB-gcc-4.9.4-Python-3.6.1

```

> *NOTE:* After loading the module, if your run the command `module list` you
> will notice that it has many dependencies which have also been loaded for you.

Now let's move into the `chipseq/results` directory and create a new directory
for the results of our IDR analysis.

```bash
$ cd ~/ngs_course/chipseq/results
$ mkdir IDR
```

Copy over the sorted narrowPeak files for each replicate for Nanog and Pou5f1:

```bash
$ cp /home/classroom/hpcbio/chip-seq/idr/macs2/*sorted_peaks* IDR/
```

### Peak consistency between true replicates

The first step is taking our replicates and evaluating how consistent they are
with one another.

<img src="../img/idr-idr.png" width=500>

To run IDR we use the `idr` command followed by any necessary parameters. To see
[what parameters](https://github.com/kundajelab/idr#command-line-arguments) we
have available to us, we can use:

```bash

# you can optionally run this, but it may take a little time
$ idr -h

```

For our run we will change only parameters associated with input and output
files. We will also change the field on which we want to create ranks since the
defaults are set for SPP peak calls. Parameters that pertain to the inference
procedure are left as is, since the chosen defaults are what the developers
believe are reasonable in the vast majority of cases.

Move into the IDR directory:

```bash
$ cd IDR
```

Let's start with the Nanog replicates:

```bash
$ idr --samples Nanog-rep1_sorted_peaks.narrowPeak Nanog-rep2_sorted_peaks.narrowPeak \
--input-file-type narrowPeak \
--rank p.value \
--output-file Nanog-idr \
--plot \
--log-output-file nanog.idr.log
```

And now with the Pou5f1 replicates (try this on your own first)

<details>

```bash
$ idr --samples Pou5f1-rep1_sorted_peaks.narrowPeak Pou5f1-rep2_sorted_peaks.narrowPeak \
--input-file-type narrowPeak \
--rank p.value \
--output-file Pou5f1-idr \
--plot \
--log-output-file pou5f1.idr.log
```

</details>
<br>

#### Output files

The output file format mimics the input file type, with some additional fields.
Note that the **first 10 columns are a standard narrowPeak file**, pertaining to
the merged peak across the two replicates.

**Column 5 contains the scaled IDR value, `min(int(log2(-125IDR), 1000)`** For
example, peaks with an IDR of 0 have a score of 1000, peaks with an IDR of 0.05
have a score of int(-125log2(0.05)) = 540, and IDR of 1.0 has a score of 0.
**Columns 11 and 12 correspond to the local and global IDR value,
respectively.** The global IDR is the value used to calculate the scaled IDR
number in column 5, it _is analogous to a multiple hypothesis correction on a
p-value to compute an FDR_. The local IDR is akin to the posterior probability
of a peak belonging to the irreproducible noise component. You can read
[this paper](http://projecteuclid.org/euclid.aoas/1318514284) for more details.

The next four columns correspond to Replicate 1 peak data and the following four
columns with Replicate 2 peak data.

More detail on the output can be [found in the user manual](https://github.com/kundajelab/idr#output-file-format).
Also, if you have any unanswered questions check out posts in the [Google groups forum](https://groups.google.com/forum/#!forum/idr-discuss).

Let's take a look at our output files. _How many common peaks are considered for
each TF?_

<details>

```bash
$ wc -l *-idr
```

</details>
<br>

To find out how may of those shared regions have an IDR < 0.05, we can take a
look at the log files. Alternatively, since we requested all peaks and their IDR
value as output we can also filter the file using an `awk` command.

```bash
$ awk '{if($5 >= 540) print $0}' Nanog-idr | wc -l
$ awk '{if($5 >= 540) print $0}' Pou5f1-idr | wc -l
```

_Which of the two TFs show better reproducibility between replicates? How does
this compare to the `bedtools` overlaps?_

#### Output plots

There is a single image file output for each IDR analyses (`.png` files). Within
each image you should see four plots. **Since we are working with such a small
subset of data, the plots are not as meaningful. Therefore, below we have
provided the images generated for Pou5f1 full dataset below**.

<img src=../img/Pou5f1-idr.png width=500>

The plot for each quadrant is described below:

**Upper Left**: Replicate 1 peak ranks versus Replicate 2 peak ranks - peaks
that do not pass the specified idr threshold are colored red.

**Upper Right**: Replicate 1 log10 peak scores versus Replicate 2 log10 peak
scores - peaks that do not pass the specified idr threshold are colored red.

**Bottom Row**: Peak rank versus IDR scores are plotted in black. The overlayed
boxplots display the distribution of idr values in each 5% quantile. The IDR
values are thresholded at the optimization precision - 1e-6 by default.

### Peak consistency between pooled pseudoreplicates

Once you have IDR values for true replicates, you want to see how this compares
to pooled replicates. This is a bit more involved, as it requires you to go back
to the BAM files, merge the reads and randomly split them into two
pseudo-replicates. If the original replicates are highly concordant, then
shuffling and splitting them should result in pseudo-replicates that the reflect
the originals. **Therefore, if IDR analysis on the pooled pseudo-replicates
results in a number of peaks that are similar (within a factor of 2) these are
truly good replicates.**

<img src="../img/idr-pool.png" width=500>

We provide [a bash script](https://github.com/HPCBio/In-depth-NGS-Data-Analysis-Course/blob/HPCBio-Fall2017/sessionV/scripts/pseudorep_idr_biocluster.sh) on the cluster for you to try this out. To run this
script you will need to:

* Provide BAM files and run it for each TF separately. These are located at
  `/home/classroom/hpcbio/chip-seq/bowtie2/`. Or you can
  point to the BAM files generated from Bowtie2 in the home directory.
* Be sure to also ask for enough memory in your `srun` command.
* Change the paths for output to the directories that are relevant to you

> _NOTE: For the paths and directories we are using `/scratch`. This script
> generates fairly large intermediate files which can quickly fill up your home
> directory. To avoid this please make use of the scratch space and once the
> analysis is complete move over only the relevant files._

Let's set up a new results directory:

```bash
$ cd ~/ngs_course/chipseq/results

$ mkdir IDR-pooled

$ cd IDR-pooled
```

Now, let's copy the script over:

```bash
$ cp /home/classroom/hpcbio/chip-seq/scripts/pseudorep_idr.sh .
```

Run on the Nanog samples:

```bash
$ ./pseudorep_idr.sh \
  ../bowtie2/H1hesc_Input_Rep1_chr12_aln.bam \
  ../bowtie2/H1hesc_Nanog_Rep1_chr12_aln.bam \
  ../bowtie2/H1hesc_Input_Rep2_chr12_aln.bam \
  ../bowtie2/H1hesc_Nanog_Rep2_chr12_aln.bam \
  pooled_Nanog
```

You should have output that looks like this:

```
Thu Nov 16 22:36:31 CST 2017
Merging BAM files for pseudoreplicates...
Merging input BAM files for pseudoreplicates...
Calling peaks for pseudoreplicate1
Calling peaks for pseudoreplicate2
Sorting peaks...
Running IDR on pseudoreplicates...
Initial parameter values: [0.10 1.00 0.20 0.50]
Final parameter values: [1.87 0.62 0.84 0.35]
Number of reported peaks - 174/174 (100.0%)

Number of peaks passing IDR cutoff of 0.05 - 55/174 (31.6%)
```

Your numbers (# peaks, peaks passing IDR, etc) will likely not match this, since
there are some steps in the process (shuffling) which lead to different results.
In general, this would need to be run many times to determine whether the #
significant peaks start to converge.

Now, try running this on the Pou5f1 files (note the order of the files in the above script).

<details>

```bash
$ ./pseudorep_idr.sh \
  ../bowtie2/H1hesc_Input_Rep1_chr12_aln.bam \
  ../bowtie2/H1hesc_Pou5f1_Rep1_chr12_aln.bam \
  ../bowtie2/H1hesc_Input_Rep2_chr12_aln.bam \
  ../bowtie2/H1hesc_Pou5f1_Rep2_chr12_aln.bam \
  pooled_Pou5f1
```

</details>
<br>

### Self-consistency analysis

An _optional step_ is to create pseudo-replicates for each replicate by randomly
splitting the reads and running them through the same workflow. Again, **if IDR
analysis on the self-replicates for Replicate 1 results in a number of peaks
that are similar (within a factor of 2) to self-replicates for Replicate 2 these
are truly good replicates.**

<img src="../img/idr-rep1-rep2.png" width=500>

### Threshold guidelines

The user manual provides [guidelines on IDR thresholds](https://sites.google.com/site/anshulkundaje/projects/idr#TOC-GETTING-THRESHOLDS-TO-TRUNCATE-PEAK-LISTS)
which are recommended for the different types of IDR analyses. Depending on the
organism you are studying and the total number of peaks you are starting with
you will want to modify the thresholds accordingly.

An example for our analysis is described below:

* If starting with < 100K pre-IDR peaks for large genomes (human/mouse): For
  true replicates and self-consistency replicates an IDR threshold of 0.05 is
  more appropriate
* Use a tighter threshold for pooled-consistency since pooling and subsampling
  equalizes the pseudo-replicates in terms of data quality. Err on the side of
  caution and use more stringent IDR threshold of 0.01

***

*This lesson has been developed by members of the teaching team at the [Harvard
Chan Bioinformatics Core (HBC)](http://bioinformatics.sph.harvard.edu/). These
are open access materials distributed under the terms of the [Creative Commons
Attribution license](https://creativecommons.org/licenses/by/4.0/) (CC BY 4.0),
which permits unrestricted use, distribution, and reproduction in any medium,
provided the original author and source are credited.*
