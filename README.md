Scripts for paper "Multi-omics analysis provides insights into musk secretion in muskrat and musk deer".

These scripts are divided into 4 main sections: genome_sequencing_data_analysis, RNAseq_data_analysis, scRNAseq_data_analysis, and Hi-C_data_analysis.

1 Genome_sequencing_data_analysis:
The folder contains Assembly_Annotation.sh and Comparative_Genomics.sh shell scripts, there are step-by-step code comments in these two scripts, please refer to the code comments to run the code.

2 RNAseq_data_analysis:
This section of the code contains 17 plotting scripts, download the whole folder directly, please install the R software and the R packages loaded in these codes in Visual Studio Code in advance. Open the folder (downloaded) in Visual Studio Code and run the R plotting scripts directly to generate plots and tables.

3 scRNAseq_data_analysis:
This section of the code contains 4 plotting scripts, download the whole folder directly, please install the R software and the R packages loaded in these codes in Visual Studio Code in advance. Open the folder (downloaded) in Visual Studio Code and run the R plotting scripts directly to generate plots and tables. The integrated_seurat.rds in the 01_rds_file directory was not included due to its size, so please obtain it from the NCBI repository corresponding to the article.

4 Hi-C_data_analysis:
This section involves the juicer software analysis process, please install the relevant software in advance, and then follow the .sh suffix script to run the code to produce results.
