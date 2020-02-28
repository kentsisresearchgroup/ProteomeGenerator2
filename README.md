# ProteomePersonalizer

ProteomePersonalizer is the first soup-to-nuts solution for contemporary proteogenomics, building upon the conceptual framework of its predecessor, ProteomeGenerator. It is a computational suite that, via comprehensive integration of genome sequencing and RNA-seq data (starting with raw FASTQ or BAM files), captures a diverse breadth of sequencing variation, structural variation, and aberrant epi/genetic events in order to generate a complete and maximally accurate, personalized protein sequence database (in FASTA format) that can be used in LC/MS proteomics searches, among other things. ProteomePersonalizer is effectively plug-n-play, only requiring installation of the common Python package manager Anaconda, Python-based pipeline manager Snakemake, and (optionally) the container/virtualization software Singularity. Finally, ProteomePersonalizer (via Snakemake) fully integrates task scheduling managers (such as LSF) for deployment to high performance clusters. It runs entirely on Linux-based systems.
