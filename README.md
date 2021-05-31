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
   * [IQ-tree error while estimating trees](#IQ-tree-error-while-estimating-trees)
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
Besides, **CURE** creates the `intergenic-regions/` dir containing unmerged UCEs assigned to intergenic regions.

Alignments in `concatenated-by-region/` and `intergenic-regions/` dir are in NEXUS format.
Alignments in `concatenated-by-gene/` are in PHYLIP format, and its charsets are in NEXUS format.

Secondary outputs of **CURE** include `CURE-exons.txt`, `CURE-introns.txt`, and `CURE-intergenic.txt`, which contains the UCE names assigned to each region, as well as the region ID (for exons) and gene ID (for exons and introns).
The intergenic file contains only the UCE names.
**CURE** also maintain in the output directory the files produces by uce_kit pipeline (`uce_kit_output/` dir)

# Kwown issues

## IQ-tree error while estimating trees

Depending on the computer you are running on, you may face some trouble with IQ-tree while running the `estimating-trees.sh`. I am not sure why this happens, but it does. Frequently the log of `estimate-trees.sh` shows that `IQ-TREE CRASHES WITH SIGNAL ABORTED` during the model finder stage. 
To face this, you can create a list of the alignment files that have not been run and rerun IQ-tree for them.  
Say for instance that this happened while estimating trees for the alignments of `concatenated-by-gene/` dir.  Let's open our terminal and define some variables with the alignment dir and tree dir before starting solving this.

```
ALI_DIR=CURE-output/concatenated-by-gene/
TREE_DIR=iqtree-output/concatenated-by-gene/
```

Now our commands might look the same.

Let's write a file called `trees_already_done.txt` containing tha name of all files that did not had trouble with IQ-tree:

```
find $TREE_DIR -name *treefile -printf "%f\n" | cut -d '.' -f 1 > trees_already_done.txt
```
Now use `grep` to get the trees that were not done:

```
find $ALI_DIR -name *charsets -printf "%f\n" | cut -d '.' -f 1 | grep -v -f trees_already_done.txt > trees_not_done.txt
```

Now you can rerun IQ-tree for each tree in the `tree_not_done.txt`. 
The following code will do it in parallel, running 10 trees at a time, with 2 threads each (make sure you have activate the cure environment):

```
for alignment in $(cat trees_not_done.txt); do 
   sem --will-cite --max-procs 10 iqtree -s "$ALI_DIR"/"$alignment".phylip \
   -spp "$ALI_DIR"/"$alignment".charsets --prefix "$TREE_DIR"/"$alignment" \
   -bb 1000 -alrt 1000 --threads-max 2
done; sem --will-cite --wait
```

If you had this issue for alignments in the `concatenated-by-region/` or `intergenic-regions` dir, the steps are almost the same.
First difference is that you need to tell `find` command to look for `.nexus` files to write `trees_not_done.txt`:

```
find $ALI_DIR -name *nexus -printf "%f\n" | cut -d '.' -f 1 | grep -v -f trees_already_done.txt > trees_not_done.txt
```

Finally, IQ-tree command line needs to be changed:

```
for alignment in $(cat trees_not_done.txt); do 
   sem --will-cite --max-procs 10 iqtree -s "$ALI_DIR"/"$alignment".nexus \
   --prefix "$TREE_DIR"/"$alignment" \
   -bb 1000 -alrt 1000 --threads-max 2
done; sem --will-cite --wait
```

# License

[GNU General Public License, version 3](https://www.gnu.org/licenses/gpl-3.0.html)
