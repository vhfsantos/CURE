<p align="center"><img src="misc/logo.png" alt="CURE" width="60%"></p>

**CURE** is an automated and parallel pipeline for the curation of ultraconserved elements (UCEs) for species-tree reconstruction. It is an automation/adaptation of the strategy proposed by [Van Dam et al. 2021. Genomic Characterization and Curation of UCEs Improves Species Tree Reconstruction](https://academic.oup.com/sysbio/article/70/2/307/5880562#227740768). 

It concatenates UCEs in two different ways, according to the available annotation:

* _by gene_: concatenates all UCEs from the same gene and treats different regions (exons and introns) as different partitions;
* _by region_: concatenates all UCEs from the same exons or introns of the same gene. 

By default, **CURE** runs both approaches, but this can be changed. The input files for the pipeline are the baits file used for UCE sequencing, the reference genome and annotation file, and the UCE alignments produced by [phyluce](https://phyluce.readthedocs.io/en/latest/).

# Table of contents

* [Installation](#installation)
* [How CURE works](#how-cure-works)
* [Quick usage examples](#quick-usage-examples)
* [Downstream analysis](#downstream-analysis)
    * [CURE output files](#cure-output-files)
    * [Estimating trees from output files](#estimating-trees-from-output-files)
    * [Summary analysis of estimated trees](#summary-analysis-of-estimated-trees)
* [License](#license)

# Installation

We recommend the installation of **CURE** with [conda](https://conda.io/) for the automatic installation of all dependencies.
First, clone this repo in your local machine and enter the created directory:

```
git clone https://github.com/vhfsantos/CURE.git
cd CURE
```
Then, create a conda environment for **CURE** using the `cure.yml`, and ensure all scripts have execution permission:

```
conda env create -n cure --file misc/cure.yml
chmod +x CURE
chmod +x scripts/*
```

After done all installations, activate the cure environment and run **CURE** with no arguments. **CURE** will tell you if any dependencies are missing.

```
conda activate cure
./CURE
```

# How CURE works

The main inputs of **CURE** are the UCE alignments and an annotated reference genome (note that **CURE** also needs to be provided with the baits file used for the UCE sequencing)

<p align="center"><img src="misc/img/input.png" alt="input" width="100%"></p>

The first step of **CURE** is running a custom version of the `uce_type` tool described by [Van Dam et al. 2021](https://academic.oup.com/sysbio/article/70/2/307/5880562#227740768) and available at the [Cal Academy's repository](https://github.com/calacademy-research/ccgutils).
Briefly, this step assigns each UCE to an exon, intron, or intergenic region of the given reference genome

<p align="center"><img src="misc/img/output1.png" alt="input" width="80%"></p>

Then **CURE** parses the results and merges the UCEs in two different ways: by gene and by region. 

When concatenating *by gene*, **CURE** merges all UCEs from the same gene and treats different regions (exons and introns) as different partitions (Note that different introns are placed under the same partition). It stores the results in `phylip` format inside the `concatenated_by_gene/` directory. 
Further phylogenetic analysis of UCEs merged with this approach would yield a phylogenetic tree for each gene. 

<p align="center"><img src="misc/img/cat_by_gene.png" alt="input" width="80%"></p>

When concatenating *by region*, **CURE** merges only UCEs from the region of the same gene.
It stores the results in `nexus` format inside the `concatenated_by_region/` directory. 
Further phylogenetic analysis of UCEs merged with this approach would yield several phylogenetic trees, one originated from each region. 

<p align="center"><img src="misc/img/cat_by_region.png" alt="input" width="80%"></p>

For any of the two concatenating approaches, **CURE** leaves unmerged UCEs in intergenic regions. 
These UCEs are just copied to the `intergenic_regions/` directory. 

<p align="center"><img src="misc/img/intergenic.png" alt="input" width="80%"></p>


# Quick usage examples


# Downstream analysis

## CURE output files

The main output files produced by **CURE** are the alignments of concatenated UCEs.
If you run **CURE** without `--only-by-gene` or `--only-by-region`, both concatenation approaches will be done.
In this case, your output dir will contain `concatenated-by-region/` and `concatenated-by-gene/` dirs.
If you raised any of these flags, only the corresponding dir will be created.
Besides, **CURE** creates the `intergenic-regions/` dir containing unmerged UCEs assigned to intergenic regions.

Alignments in `concatenated-by-region/` and `intergenic-regions/` dir are in NEXUS format.
Alignments in `concatenated-by-gene/` are in PHYLIP format, and its charsets are in NEXUS format.

Secondary outputs of **CURE** include `CURE-exons.txt`, `CURE-introns.txt`, and `CURE-intergenic.txt`, which contains the UCE names assigned to each region, as well as the region ID (for exons) and gene ID (for exons and introns).
The intergenic file contains only the UCE names.
**CURE** also maintain in the output directory the files produces by uce_kit pipeline (`uce_kit_output/` dir)

## Estimating trees from output files

## Summary analysis of estimated trees

# License

[GNU General Public License, version 3](https://www.gnu.org/licenses/gpl-3.0.html)
