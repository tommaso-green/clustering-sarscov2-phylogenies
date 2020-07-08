# Clustering SARS-CoV2 phylogenetic trees and entropy-based validation using country of origin

This code was developed for a project for the Bioinformatics course of the University of Padua.
After having produced 4 different phylogenies using neighbor-joining
on the results of alignment and alignment-free procedures (provided in the tree_files folder)
we used [TreeCluster](https://github.com/niemasd/TreeCluster) to perform clustering.
Validation was implemented by using entropy on the country of origin of each
sequence.
Before executing the program, we converted our Newick tree Files in one-line Newick files so
that TreeCluster was able to parse them.
## Usage

If you need to convert your Newick tree files, just put them inside tree_files folder
and execute (the results will be saved in tree_files/converted):

```bash
python3 tree_converter.py
```
To perform the clustering and the entropy validation:

```bash
python3 clustering.py  -t THRESHOLD_FOR_ALIGNMENT_PHYLOGENY -td2 THRESHOLD_FOR_EPSIMD2_PHYLOGENY -td2* THRESHOLD_FOR_EPSIMD2*_PHYLOGENY [-v verbose] [-tf THRESHOLD_FREE]
```
The mandatory parameters are:
* -t which specifies the threshold for the trees derived using alignment (pairwise or multiple);
* -td2 for the threshold of the tree built using D2 alignment-free procedure;
* -td2* for the threshold of the tree built using D2* alignment-free procedure.

In outputs folder you can find the results for the execution of:

```bash
python3 clustering.py -t 7.12E-4 -td2 0.4 -td2* 0.65 -m sum_branch
```
## Requirements
* [NiemaDS](https://github.com/niemasd/NiemaDS)
* [TreeSwift](https://github.com/niemasd/TreeSwift)
* [Astropy](https://www.astropy.org/)
* [Matplotlib](https://matplotlib.org/) (with support to LaTeX)