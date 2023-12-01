# Description 
DeepBlood is a genomic analysis tool that translates genotype data into blood group alleles. The software operates by accessing the ISBT database dump, containing a detailed description of blood group alleles and their underlying genomic variations. Eventually, it is a transcript of the ISBT pdf documentation (ISBT Blood Group Allele Tables). A second, self-made configuration file connects the ISBT variation annotation with the corresponding VCF file entries. The ISBT entries are annotated on specific cDNA sequences, which can differ from the human reference genome sequence. In addition, the VCF file is related to the human genome reference sequence and follows the left-aligned variation annotation, which is not necessarily the case for the ISBT entries. The mandatory input files are VCF, and genome wide coverage, all stored efficiently, like the configuration files as well, as C++ objects in-memory for swift accessibility. The coverage information could be generated from a BAM file, but we decided to do it with UCSC Genome browser tools and provide it as a separate input, as it is just a little extra step during secondary data analysis and the size is relatively small. That would it make easier to provide DeepBlood as a web service.

## Key Features:

The software recombines measured genotypes/haplotypes and compares the derived haplotypes with the ISBT DB entries.

### Real-time Genotype Data Access: 
All configuration data (ISBT related input files) is stored as C++ objects in-memory, accessible through dedicated functions. The same counts for the individual genome data, but the dataset is reduced to the data points required for the ISBT lookup.

### Automated Analysis and Allele Determination:
DeepBlood scores the matches between the recombined haplotypes and the ISBT DB entries. Implementing a custom scoring algorithm is straightforward.

### Efficiency: 
Quick, efficient, and multi-threading capable.

### Flexibility:
Object-oriented design enables easy adaptation to personal preferences. 
