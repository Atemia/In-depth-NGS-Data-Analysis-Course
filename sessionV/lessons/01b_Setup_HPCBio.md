---
title: "Introduction to ChIP-Seq: Setup"
author: "Mary Piper, Radhika Khetani, Meeta Mistry (HSPH), Chris Fields (UIUC)"
date: "October 30, 2017"
output: revealjs::revealjs_presentation
---

## Set-up

Before we get started with the analysis, we need to set up our directory structure.

Login to biocluster2 and start an interactive session with two cores:

```bash
$ srun -c 2 --mem 2000 --pty bash
```

---

Create a new work directory

```bash
$ mkdir ngs_course 
$ cd ~/ngs_course
```

---

Create a `chipseq` directory and change directories into it:

```bash
$ mkdir chipseq
$ cd chipseq
```

---

Now let's setup the directory structure, we are looking for the following structure within the chipseq directory:

```text
chipseq/
├── logs/
├── meta/
├── raw_data/
├── reference_data/
├── results/
│   ├── bowtie2/
│   ├── trimmed/
│   ├── trimmed_fastqc/
│   └── untrimmed_fastqc/
└── scripts/
```

---

Here's how we create this directory structure

```bash
$ mkdir -p raw_data reference_data scripts logs meta

$ mkdir -p results/untrimmed_fastqc results/trimmed results/trimmed_fastqc results/bowtie2
```

---

Now that we have the directory structure created, let's copy over the data to perform our quality control and alignment, including our FASTQ files and reference data files:

```bash
$ cp /home/classroom/hpcbio/chip-seq/raw_fastq/*fastq raw_data/
$ cp /home/classroom/hpcbio/chip-seq/chr12* reference_data/
```

***
*This lesson has been developed by members of the teaching team at the [Harvard Chan Bioinformatics Core (HBC)](http://bioinformatics.sph.harvard.edu/). These are open access materials distributed under the terms of the [Creative Commons Attribution license](https://creativecommons.org/licenses/by/4.0/) (CC BY 4.0), which permits unrestricted use, distribution, and reproduction in any medium, provided the original author and source are credited.*
