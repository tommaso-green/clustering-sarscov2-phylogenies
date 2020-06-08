# Clustering SARS-CoV2 phylogenetic trees and entropy-based validation using country of origin

This code was developed for a project for the Bioinformatics course of the University of Padua, under the supervision of Prof. Matteo Comin.
After having produced 4 different phylogenies using neighbor-joining
on the results of alignment and alignment-free procedures (provided in the tree_files folder)
we used [TreeCluster](https://github.com/niemasd/TreeCluster) to perform clustering.
Validation was implemented by using entropy on the country of origin of each
sequence.
Before executing the program, we converted our Newick tree Files in one-line Newick files so
that TreeCluster was able to parse them.
## Usage
```bash
usage: python clustering.py  -t THRESHOLD_FOR_ALIGNMENT_PHYLOGENY -td2 THRESHOLD_FOR_EPSIMD2_PHYLOGENY -td2* THRESHOLD_FOR_EPSIMD2*_PHYLOGENY [-v verbose] [-tf THRESHOLD_FREE]
```


## Requirements
* [NiemaDS](https://github.com/niemasd/NiemaDS)
* [TreeSwift](https://github.com/niemasd/TreeSwift)
* [Astropy](https://www.astropy.org/)
* [Matplotlib](https://matplotlib.org/)