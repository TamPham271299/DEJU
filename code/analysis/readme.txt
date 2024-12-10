This directory contains scripts used to generate figures for,
- manuscript
- supplementary materials.

## Directory structure

analysis/
├── edgeR_diffSpliceDGE_simple.sh                  # Functions to run edgeR::diffSplice
├── limma_diffSpliceDGE_simple.sh              # Functions to run limma::diffSpliceDGE
├── main_figs_codes.R               # Script to generate figures for manuscript
├── supp_figs_codes.R                # Script to generate figures for supplementary data
├── plotJunc3.R                # Script to generate schematic exon-junction plots for limma::diffSplice
├── plotJunc3_diffSpliceDGE.R                # Script to generate schematic exon-junction plots for edgeR::diffSpliceDGE
├── make_bam_for_visualization                # prepare truncated BAM files for specific genes
└── readme.txt                  # Main documentation for the repository
