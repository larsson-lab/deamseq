This repository complements our study "Light-based footprinting of a eukaryotic genome" (DOI: 10.1126/sciadv.aec0059).

Programs and versions used to run all scripts are:
BWA 0.7.17-r1188
Samtools 1.9
UMI-tools 1.1.4
VarScan 2.3
Matlab R2024b


Variant processing and analysis in MATLAB, including all figures in the paper, is done as described below (adaptations to scripts, e.g. filenames, are needed if processing other Deam-seq datasets).
- setup_workspace.m: Load Deam-seq calls and perform basic processing. This script stores all key Deam-seq data variables, including normalized cell/naked ratios, in workspace.mat.
- setup_extra_data.m: Load and process some extra datasets, including yeast genomic data, reused in several analyses and needed to generate most figures. Stored in extra_data.mat.

To generate Figures 1-5 and Supplementary Figures S1-S10:
- plot/*.m: Scripts for each individual figure. Will load workspace mat  files generated above as needed.

To generate WIG files:
- make_wig.m: Generate wig files, written to the wig/ subdirectory.

To generate stats enabling comparison of deduplicated and non-duplicated samples:
- nondedup.stats.m: Saves nondedup_stats.mat (data used for Figure 1).


The subdirectory 'shell_scripts' contains the steps for generating the mutation calls.
- shell_scripts/align.sh: Align using BWA.
- shell_scripts/dedup.sh: Perform UMI-based deduplication (if applicable).
- shell_scripts/mutation_call.sh: Call mutations using VarScan.

Also in the same directory:
- shell_scripts/process_chip.sh: Used to align and get coverage for the Reb1 ChIP-exo data.
