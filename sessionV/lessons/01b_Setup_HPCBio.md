## Set-up

Before we get started with the analysis, we need to set up our directory structure.

Login to Orchestra and start an interactive session with two cores:

```bash
$ bsub -Is -n 2 -q interactive bash
```

---

Change directories to the `ngs_course` directory:

```bash
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

```bash
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

```bash
$ mkdir -p raw_data reference_data scripts logs meta

$ mkdir -p results/untrimmed_fastqc results/trimmed results/trimmed_fastqc results/bowtie2
```

---

Now that we have the directory structure created, let's copy over the data to perform our quality control and alignment, including our FASTQ files and reference data files:

```bash
$ cp /groups/hbctraining/ngs-data-analysis-longcourse/chipseq/raw_fastq/*fastq raw_data/

$ cp /groups/hbctraining/ngs-data-analysis-longcourse/chipseq/reference_data/chr12* reference_data/
```

---

You should have bcbio in you path, but please check that it is:

```bash
$ echo $PATH
```

---

If `/opt/bcbio/centos/bin` is not part of `$PATH`, add it by adding the following line within your `~/.bashrc` file and then run `source ~/.bashrc`:

```bash
export PATH=/opt/bcbio/centos/bin:$PATH
```

---

***
*This lesson has been developed by members of the teaching team at the [Harvard Chan Bioinformatics Core (HBC)](http://bioinformatics.sph.harvard.edu/). These are open access materials distributed under the terms of the [Creative Commons Attribution license](https://creativecommons.org/licenses/by/4.0/) (CC BY 4.0), which permits unrestricted use, distribution, and reproduction in any medium, provided the original author and source are credited.*
