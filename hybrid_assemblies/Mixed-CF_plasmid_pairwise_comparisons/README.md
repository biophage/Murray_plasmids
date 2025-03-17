# PAE vs modern plasmid pairwise comparisons
## Protein-level comparison of PAE and modern plasmids of the same Mixed-CF
This directory contains the output of the pairwise comparisons we carried out across close families (CFs) that contain both pre-antibiotic era (PAE) and modern plasmids and in which we identified AMR genes (aka Mixed-CF AMR+; see our paper for details). The compared plasmids were selected randomly; the pMUR prefix identifies the PAE plasmids from the Murray collection.

Each of the 29 sub-directories corresponds to a different Mixed-CF AMR+. ___Note that the Mixed-CFs listed here correspond to those generated from the analysis of PAE plasmid hybrid assemblies (i.e. produced by combining short- and long-read sequencing data)___.

All sub-directories contain the following files:

- An _HTML_ file containing the results of the pairwise comparisons. This file was produced by clinker (https://github.com/gamcil/clinker) and it provides an interactive visualisation of the comparison: download the file and open it with a web browser to access the data (follow the embedded instructions to navigate and manipulate the genome maps). Genes encoding AMR genes and transposases are represented by red and blue arrows, respectively. _NOTE_: To access the predicted function of any gene in the genome maps, right-click on the corresponding arrow.
- A _TSV_ file describing the AMR genes identified in the plasmids featured in the comparison. This table includes several metadata; e.g. gene name, product, strand and coordinates of the gene in the corresponding plasmid (see ID in the "SEQUENCE" column to match with the genome map).

Additionally, the sub-directories of the five CFs selected as exemplars for our paper also include the _SVG_ file of the comparison presented in our figure.
