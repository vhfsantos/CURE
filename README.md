<p align="center"><img src="misc/logo.png" alt="CURE" width="40%"></p>

**CURE** is an automated and parallel pipeline for curation of phylogenomic analysis of ultraconserved elements (UCEs).
Basically, **CURE** concatenates UCEs with two different approaches (i) by gene: all UCEs bellonging to a gene are concatenated, and the different regions (exons, introns, and intergenic regions) are partinioned; 
(ii) by region: all UCEs belloging to the same exons or to introns of same gene are concatenated.
For such, **CURE** takes as input the baits file used for UCE sequencing, the reference genome and annotation file, and the UCE alignments produced by [phyluce](https://phyluce.readthedocs.io/en/latest/).

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
