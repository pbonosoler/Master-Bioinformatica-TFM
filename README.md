This repository contains all scripts used in the analysis for my master’s thesis. Below is a description of each script and its purpose.
Scripts:

    Data Selection
        dataparser.py: Selects and extracts data from the server.

    Phenotypic Data Assignment
        pdata.py: Adds phenotypic data to each sample.

    Quality Control of Reads
        check_quality.sh: Checks the quality of reads from FASTQ files.

    Alignment to Reference Genome
        align_paired.sh / align_single.sh: Aligns paired-end or single-end reads to the reference genome.

    Count Table Generation
        generate_count_table.sh / generate_count_table_single.sh: Generates count tables from the aligned reads.

    Interaction Analysis
        genetic_interactions.py: Gets the number of genetic interactions for each gene.
        protein_interactions.py: Gets the number of PPIs for each identifier.

    Preprocessing, Differential Expression, Gene Origin and Pathway Analysis:
        DEAedgeR.R: Analyzes differential expression, gene origin, pathways, interaction counts, and more.
