<p align="center"><img src="misc/logo.png" alt="CURE" width="60%"></p>

**CURE** is an automated and parallel pipeline for the curation of ultraconserved elements (UCEs) for species-tree reconstruction. It is an automation/adaptation of the strategies proposed by [Van Dam et al. 2021. Genomic Characterization and Curation of UCEs Improves Species Tree Reconstruction](https://academic.oup.com/sysbio/article/70/2/307/5880562#227740768) (named **GeneRegion** strategy), and [Freitas et al. 2021 Partitioned Gene-Tree Analyses and Gene-Based Topology Testing Help Resolve Incongruence in a Phylogenomic Study of Host-Specialist Bees (Apidae: Eucerinae)](https://academic.oup.com/mbe/article/38/3/1090/5976982) (name **UCERegion** strategy).

In the **GeneRegion** strategy (Van Dam et al. 2021) it concatenates UCEs in two different ways, according to the available annotation:

* _by gene_: concatenates all UCEs from the same gene and treats different regions (exons and introns) as different partitions;
* _by region_: concatenates all UCEs from the same exons or introns of the same gene.

By default, **CURE** runs both approaches for **GeneRegion** strategy, but this can be changed. The input files for this pipeline are the baits file used for UCE sequencing, the reference genome, and annotation file, and the UCE alignments produced by [phyluce](https://phyluce.readthedocs.io/en/latest/).

In the **UCERegion** strategy (Freitas et al. 2021) it runs SWSC-EN ([Tagliacollo & Lanfear 2018](https://academic.oup.com/mbe/article-abstract/35/7/1798/4969532)) in a parallelized way, that speeds up the process a lot, and creates charsets considering the left flank, right flank and core for each locus in a dataset for gene-tree estimation.

# Table of contents

* [Installation](#installation)
* [How CURE works](#how-cure-works)
* [Quick usage examples](#quick-usage-examples)
* [Output files](#output-files)
* [Estimating trees from output files](#estimating-trees-from-output-files)
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

## **GeneRegion** Strategy

The main inputs for this strategy are the UCE alignments and an annotated reference genome (note that **CURE** also needs to be provided with the baits file used for the UCE sequencing)

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
Further phylogenetic analysis of UCEs merged with this approach would yield several phylogenetic trees, one originating from each region.

<p align="center"><img src="misc/img/cat_by_region.png" alt="input" width="80%"></p>

For any of the two concatenating approaches, **CURE** leaves unmerged UCEs in intergenic regions.
These UCEs are just copied to the `intergenic_regions/` directory.

<p align="center"><img src="misc/img/intergenic.png" alt="input" width="80%"></p>


# Quick usage examples

## Running the test dataset for **Gene regions**

You can test **CURE** with the test dataset.
It usually takes about two minutes to run with 10 threads.
With the command line below, **CURE** will run the two concatenating approaches.

```sh
CURE GeneRegion --baits test_data/baits.fasta  \
                --reference test_data/ref.fa \
                --gff test_data/ref.gff \
                --phyluce-nexus test_data/uce_nexus/ \
                --output ./CURE-GeneRegion-output
```

### Running only one of the concatenating approaches

By default, **CURE** runs both concatenating approaches.
However, you can raise the `--only-by-gene` or `--only-by-region` flag to select only a single approach

#### Only by gene

```sh
CURE GeneRegion --baits test_data/baits.fasta  \
                --reference test_data/ref.fa \
                --gff test_data/ref.gff \
                --phyluce-nexus test_data/uce_nexus/ \
                --output ./CURE-GeneRegion-output \
                --only-by-gene
```

#### Only by gene region

```sh
CURE GeneRegion --baits test_data/baits.fasta  \
                --reference test_data/ref.fa \
                --gff test_data/ref.gff \
                --phyluce-nexus test_data/uce_nexus/ \
                --output ./CURE-GeneRegion-output \
                --only-by-region
```

## **UCERegion** Strategy

For this strategy you will need to provide a folder with all the alignments you want to use in nexus format (could be all your alignemnts or a subset).
You will ned to have [SWSC-EN](https://github.com/Tagliacollo/PFinderUCE-SWSC-EN) included in yout PATH environment or use --swsc to provide the path to the script.
What CURE will do is basically running SWSC-EN in paralell and use PHYLUCE to split the alignments according to regions identified by SWSC and re-concatenate them creating a charset file to be used in phylogenetic analyses to generate your gene-trees. For this last step check our [script](#estimating-trees-from-output-files) to do that using GNU parallel.

## Running the test dataset for UCE regions

```sh
CURE UCERegion  --phyluce-nexus test_data/uce_nexus/ \
                --output ./CURE-UCERegion-output
                

```

# Output files

## UCERegion

**CURE** UCERegion will generate three subfolders
1- logfiles: containing all the log files generated
2- partitioned-uces: with all the alignemnts and their respective charsets files
3- PF2-input: input file for PartitionFinder2 analysis

## GeneRegion

The main output files produced by **CURE** GeneRegion approach are the alignments of concatenated UCEs.
If you run **CURE** without `--only-by-gene` or `--only-by-region`, both of the concatenating approaches will be done.
In this case, your output-dir will contain `concatenated-by-region/` and `concatenated-by-gene/` dirs.
If you raised any of these flags, only the corresponding dir will be created.
Besides, **CURE** creates the `intergenic-regions/` dir containing unmerged UCEs assigned to intergenic regions.

Alignments in `concatenated-by-region/` and `intergenic-regions/` dir are in NEXUS format.
Alignments in `concatenated-by-gene/` are in PHYLIP format, and its charsets are in NEXUS format.

> To avoid troubles with further phylogenetic analysis, **CURE** replaces "-" with "_" in the gene and exon ID.

Secondary outputs of **CURE** include `CURE-exons.txt`, `CURE-introns.txt`, and `CURE-intergenic.txt`, which contains the UCE names assigned to each region, as well as the region ID (for exons) and gene ID (for exons and introns).
The `CURE-intergenic.txt` file contains only the UCE names.
Besides, **CURE** outputs the `CURE-summary.csv` file containing the UCE count assigned to exons, to introns, to both exon and intron, to intergenic regions, and unassigned UCEs.
UCEs assigned to both exon and intron are accounted for exons, and unassigned UCEs are accounted for intergenic regions.

**CURE** also maintains in the output directory the files produces by the uce_kit pipeline (`uce_kit_output/` dir)

# Estimating trees from output files

**CURE** provides the wrapper script `estimate-trees.sh` for the estimation of gene trees from the output alignments with [IQ-tree](http://www.iqtree.org/), and further summary analysis with [ASTRAL](https://github.com/smirarab/ASTRAL).
This script runs IQ-tree in parallel using [GNU Parallel](https://www.gnu.org/software/parallel/) following the structure of the **CURE** output-dir.
Then it prepares all inputs needed for a summary analysis with ASTRAL.
For instance, if you run **CURE** setting `CURE-output` as output directory, you can call `estimate-trees.sh` as the following:

```sh
scripts/estimate-trees.sh \
    --cure-out CURE-output/ \
    --estimated-trees estimated-trees
```

If you raised `--only-by-gene` or `--only-by-region` while running **CURE**, you can raise it here as well:

```sh
scripts/estimate-trees.sh \
    --cure-out CURE-output/ \
    --estimated-trees estimated-trees \
    --only-by-gene
```

or

```sh
scripts/estimate-trees.sh \
    --cure-out CURE-output/ \
    --estimated-trees estimated-trees \
    --only-by-region
```

Moreover, `estimate-trees.sh` can be used to estimate trees from alignments from any other source; not necessarily those produced by **CURE**.
In this case, you only need to use the parameter `--custom-alignments` instead of `--cure-out`.
So if you have a set of alignments (in Phylip, Fasta, or Nexus format) in a directory called `input-alignments`, and want to run IQ-tree on them, you can call `estimate-trees.sh` as the following:

```sh
scripts/estimate-trees.sh \
    --custom-alignments input-alignments/ \
    --estimated-trees estimated-trees
```

# Citation

If you use **CURE** in your research, please cite:

Felipe V. Freitas, Michael G. Branstetter, Vinicius H. Franceschini-Santos, Achik Dorchin, Karen Wright, Margarita Lopez-Uribe, Terry Griswold, Fernando A. Silveira, Eduardo A.B. Almeida, (xxxx). Phylogenomics, biogeography, and classification of long-horned bees (Apidae: Eucerini), with insights on the use of specimens with extremely degraded DNA for UCE phylogenomics. xxxxxxx

# License

[GNU General Public License, version 3](https://www.gnu.org/licenses/gpl-3.0.html)
