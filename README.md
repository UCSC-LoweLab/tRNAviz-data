# F-451

Sequence feature analysis of 451 fungi

This repository is nowhere near production ready! Brave souls only!

## Pipeline overview

0. Prerequisites
	- Python 3: pandas, Biopython
	- R: plyr, dplyr, tidyr, ggplot2, Biostrings, RColorBrewer

1. Find taxonomy information for list of species along with relevant filepaths
	- Prerequisites: List of NCBI taxonomic IDs (`taxids.tsv`)
	- Scripts:
		- `get-taxonomy.R` (`taxids.tsv` -> `genomes.tsv`)
			- Output list of species, plus taxonomic info from NCBI. 
		- `generate_newick_tree.py` (`genomes.tsv` -> `newick-tree.txt`)
			- Output a Newick tree that can be used for visualizing phylogenetic trees.

2. Run tRNAscan-SE on all genomes. We need `.out`, `.iso`, and `.ss` files, and to run with detailed output.

3. Parse tRNAscan-SE output
  - Prerequisites:
  	- tRNAscan-SE output files (`.out`, `.ss`, `.iso`)
  	- A table containing paths for `.out`, `.iso`, and `.ss` files. Can be combined with `genomes.tsv` as additional columns. (`genomes.tsv`)
  	- Covariance model specialized for alignment and numbering (`numbering.sto`)
  	- Confidence set (if applicable, see notes) (`confidence-set.txt`)
  - Scripts:
  	- `parse-tRNAs.py`: main driver script (`genomes.tsv`, `.out`, `.iso`, `.ss`, `numbering.sto`, `confidence-set.txt` -> `trna-df.tsv`)
  	- `tRNA_position.py`: helper library that resolves positions
  	- `sstofa3` (`.ss` -> `.fa`)
  		-  Removes introns by parsing the `.ss` file. Note that the default output is the same name as tRNAscan-SE `.fa` output, which contains introns! This is not okay! 


4. Resolve consensus features
	- Prerequisites:
		- Parsed data frame with each position as a column (`trna-df.tsv`)
	- Scripts:
		- `resolve-consensus.py` (`trna-df.tsv` -> `consensus.tsv`, `isotype-specific.tsv`, `isotype-clade-specific.tsv`)
			- `consensus.tsv`: Lists positions with consensus features.
			- `isotype-specific.tsv` and `isotype-clade-specific.tsv`: Resolved features by top-level clade and isotype. Extending to resolve by a lower level clade is left as an exercise to the user.


## Tuning the pipeline

Depending on your genome set, you will certainly need to make some changes to your workflow. Here's some examples:

- Adding a couple of lower scoring tRNAs is helpful for identifying minor conserved variants. For non-mammalian vertebrates though, I've locked in the tRNA set to the high confidence set, due to the sheer number of highly amplified tRNA genes.
- Fungi have more diverged tRNAs and score a bit lower using the eukaryotic covariance model, so a score threshold needs to be tuned to maximize the quality of your tRNA set.