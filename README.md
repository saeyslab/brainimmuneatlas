# Brain immune atlas
### About
The brain immune atlas provides a resource and visualization tool for assessing various single-cell RNA sequencing datasets that together capture the diversity of the brain immune compartment, as published in Nature Neuroscience (*[A single-cell atlas of mouse brain macrophages reveals unique transcriptional identities shaped by ontogeny and tissue environment]*). Using the 10x genomics chromium platform, we analyzed more than 61.000 CD45+ immune cells that were obtained from whole brains or isolated border regions, including the dura mater, subdural meninges and choroid plexus. Datasets were obtained from WT homeostatic mice, conditional Irf8 knockout animals and aged APP/PS1 transgenics.

Click [here] for the full manuscript.

### Overview scripts
Here you can find the scripts used to analyse all the data:

+ script_bulkRNAseq.R = script to analyse bulk RNA-seq data
+ script_scRNAseq.R = consensus script to analyse single cell RNA-seq data
+ script_scRNAseq_choroidPlexusWtK11.R = script to analyse the WT choroid plexus sample. Here we did a correction for dissociation
+ script_scenic.R = consensus script to run __[SCENIC](https://github.com/aertslab/SCENIC)__ on single cell RNA-seq data
+ script_scorpius_durWt_K22.R = trajectory analysis of the WT dura sample using __[SCORPIUS](https://github.com/rcannood/SCORPIUS)__

### Webtool
Browse through all the data via our webtool: __[http://www.brainimmuneatlas.org/](http://www.brainimmuneatlas.org/)__

### Citation
Hannah Van Hove, Liesbet Martens, Isabelle Scheyltjens, Karen De Vlaminck, Ana Rita Pombo Antunes, Sofie De Prijck, Niels Vandamme, Sebastiaan De Schepper, Gert Van Isterdael, Charlotte L. Scott, Jeroen Aerts, Geert Berx, Guy E. Boeckxstaens, Roosmarijn E. Vandenbroucke, Lars Vereecke, Diederik Moechars, Martin Guilliams, Jo A. Van Ginderachter, Yvan Saeys and Kiavash Movahedi.*A single-cell atlas of mouse brain macrophages reveals unique transcriptional identities shaped by ontogeny and tissue environment.* Nature Neuroscience, 2019.
