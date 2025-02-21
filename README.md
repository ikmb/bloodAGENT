# Description 
bloodAGENT is a genomic analysis tool that translates genotype data into blood group alleles. The software operates by accessing the ISBT database dump, containing a detailed description of blood group alleles and their underlying genomic variations. Eventually, it is a transcript of the ISBT pdf documentation (ISBT Blood Group Allele Tables). A second, self-made configuration file connects the ISBT variation annotation with the corresponding VCF file entries. The ISBT entries are annotated on specific cDNA sequences, which can differ from the human reference genome sequence. In addition, the VCF file is related to the human genome reference sequence and follows the left-aligned variation annotation, which is not necessarily the case for the ISBT entries. The mandatory input files are VCF, and genome wide coverage, all stored efficiently, like the configuration files as well, as C++ objects in-memory for swift accessibility. The coverage information could be generated from a BAM file, but we decided to do it with UCSC Genome browser tools and provide it as a separate input, as it is just a little extra step during secondary data analysis and the size is relatively small. That would it make easier to provide bloodAGENT as a web service.

## Key Features:

The software recombines measured genotypes/haplotypes and compares the derived haplotypes with the ISBT DB entries.

### Real-time Genotype Data Access: 
All configuration data (ISBT related input files) is stored as C++ objects in-memory, accessible through dedicated functions. The same counts for the individual genome data, but the dataset is reduced to the data points required for the ISBT lookup.

### Automated Analysis and Allele Determination:
bloodAGENT scores the matches between the recombined haplotypes and the ISBT DB entries. Implementing a custom scoring algorithm is straightforward.

### Efficiency: 
Quick, efficient, and multi-threading capable.

### Flexibility:
Object-oriented design enables easy adaptation to personal preferences. 

# Installation
There are some dependencies that I solve in a special way. So I do not implement everything by submodules but clone repositories next to each other. 

## Installing the tool at the ikmbhead.uni-kiel.de

Setup environment first

`module load gcc/8.3.0`

To install this pipeline, clone the following repositories and build the libraries and executable by copy&paste the following commands:

`git clone git@github.com:ikmb/BfxCppClasses.git MyTools`

`cd MyTools`

`make all`

`cd ..`

MyTools libs created

`git clone git@github.com:samtools/htslib.git`

`cd htslib`

`git submodule update --init --recursive`

`make`

`cd ..`

htslib created

`git clone git@github.com:dpryan79/libBigWig.git`

`cd libBigWig`

`make`

`cd ..`

libBigWig created

`git clone --recursive git@github.com:ikmb/bloodAGENT.git`

`cd bloodAGENT`

Optional: `git submodule update --init --recursive`

If needed, create library symlinks. MyTools library is in the subfolder dist/[Release|Debug]/GNU-Linux/. htslib and libbigwig in their root folder: `export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/work_ifs/sukko545/haemo/tool/`

`make all`

`cd ..`

ngsblood created


# setup environment and run

`export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:<PATH to libhts.so.2> #/ifs/data/nfs_share/sukko545/software/htslib`

`./dist/Release/GNU-Linux/bloodAGENT`

run blood typing example:


## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details
Components of this project are licensed under different licenses



