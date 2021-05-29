<p align="center"><img src="misc/logo.png" alt="CURE" width="40%"></p>

**CURE** is an automated and parallel pipeline for curation of ultraconserved elements (UCEs) for species-tree reconstruction. It is an automation/adaptation of the strategy proposed by [Van Dam et al. 2021. Genomic Characterization and Curation of UCEs Improves Species Tree Reconstruction](https://academic.oup.com/sysbio/article/70/2/307/5880562#227740768). 

It concatenates UCEs in two different ways, according to the available annotation:

* _by gene_: concatenates all UCEs from the same gene and treats different regions (exons ans introns) as different partitions;
* _by region_: concatenates all UCEs from the same exons or from introns of the same gene. 

By default, **CURE** runs both approaches, but it can be changed. The input files for the pipeline are the baits file used for UCE sequencing, the reference genome and annotation file, and the UCE alignments produced by [phyluce](https://phyluce.readthedocs.io/en/latest/).

# Table of contents

* [Requirements](#requirements)
* [Installation](#installation)
* [How CURE works](#how-cure-works)
* [Quick usage examples](#quick-usage-examples)
* [Downstream analysis](#downstream-analysis)
    * [CURE output files](#cure-output-files)
    * [Estimating trees from output](#estimating-trees-from-output)
    * [Summary analysis of estimated trees](#summary-analysis-of-estimated-trees)
* [Known issues](#known-issues)
* [Acknowledgements](#acknowledgements)
* [License](#license)

# Installation

We recommend the installation of **CURE** with [conda](https://conda.io/) for the automatic installation of all dependencies.
First, clone this repo in your local machine and enter the created directory:

```
git clone https://github.com/vhfsantos/CURE.git
cd CURE
```
Then, create a conda environment for **CURE** using the `cure.yml`:

```
conda env create -n cure --file misc/cure.yml
```

After all installation, activate the cure environment and run **CURE** with no arguments. **CURE** will tell you if any dependency are missing.

```
conda activate cure
./CURE
```

# Downstream analysis

## CURE output files

The main output files produced by **CURE** are the alignments of concatenated UCEs.
If you run **CURE** without `--only-by-gene` or `--only-by-region`, both concatenation approaches will be done.
In this case, your output dir will contain `concatenated-by-region/` and `concatenated-by-gene/` dirs.
If you raised any of this flags, only the coresponding dir will be created.
Secondary outputs of **CURE** include `CURE-exons.txt`, `CURE-introns.txt`, and `CURE-intergenic.txt`, which contains the UCE names assigned to each region, as well as the region ID (for exons) and gene ID (for exons and introns).
The intergenic file contains only the UCE names.
**CURE** also maintain in the output directory the files produces by uce_kit pipeline (`uce_kit_output/` dir)

# License

[GNU General Public License, version 3](https://www.gnu.org/licenses/gpl-3.0.html)
