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


# License

[GNU General Public License, version 3](https://www.gnu.org/licenses/gpl-3.0.html)
