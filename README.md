# bloodAGENT

## Introduction
**bloodAGENT** (Blood Antigen GENo Typer) is an open-source software tool designed for the determination of blood group alleles based on genetic markers. By analyzing genomic data from Next-Generation Sequencing (NGS) and Third-Generation Sequencing (TGS), bloodAGENT resolves blood group alleles and provides insights into genomic variations.

## Key Features
- **High accuracy** in allele determination under typical conditions.
- **Modular and flexible architecture**, allowing for future adaptations.
- **Uses cosine similarity scoring** to determine the best haplotype match.
- **Supports VCF and BigWig file formats** for variant and coverage data.
- **Open-source** for transparency and community collaboration.

## System Requirements
**Supported platforms:** Compatible with Windows, macOS, and Linux through the Singularity image bloodagent.sif
**Dependencies:** 
- GCC or Clang compiler
- `libhts` library: `https://github.com/samtools/htslib` (for vcf file reading)
- `libBigWig` library https://github.com/dpryan79/libBigWig.git` (for coverage file reading)
- `MyTools` library `https://github.com/ikmb/BfxCppClasses` (for file parsing functionality)
- Python (for parsing output files)
- `https://github.com/mirror/tclap.git` (for command-line argument parsing)
- `https://github.com/nlohmann/json` (for JSON output generation)

## Installation
1. Install all dependencies:
   ```sh
   sudo apt-get update && sudo apt-get install -y \
    g++ \
    make \
    zlib1g-dev \
    libbz2-dev \
    git \
    liblzma-dev \
    libcurl4-openssl-dev \
   ```
2. Clone the repository:
   ```sh
   git clone --recurse-submodules https://github.com/ikmb/bloodAGENT.git
   cd bloodAGENT
   git submodule update --init --recursive
   ```
3. Navigate to the project folder (if not already there):
   ```sh
   cd bloodAGENT
   ```
4. Build the software:
   ```sh
   cd external/htslib
   make
   cd ../libBigWig
   make
   cd ../
   make
   ```
5. Setup environment
   ```
   # find the two external libraries
   find . -iname "libhts.so" -o -iname "libBigWig.so"
   # add them to the LD_LIBRARY_PATH variable
   export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:<PATH to libhts.so.>:<PATH to libBigWig.so.>
   ```
6. Verify installation (e.g.):
   ```sh
   ./dist/Release/GNU-Linux/bloodAGENT --help
   ```

## Input Data Format
bloodAGENT requires two main input files:
- **VCF files**: Represent genomic variants (compatible with hg19 and hg38 reference genomes).
- **BigWig files**: Provide sequencing coverage information to determine the sequencing depth of SNVs that are not listed in the VCF file

Additionally, three configuration files are needed:
- **./data/config/exonic_annotation.${build}.BGStarget.txt**: Transcript annotation for blood group targets. A separate file for each genome build.
- **./data/config/variation_annotation_${Sec.Analysis.Pipeline}.dat**: Variant annotation for different pipelines. Pipelines means the combination of read aligner and variant caller.
- **./data/config/genotype_to_phenotype_annotation_${Sec.Analysis.Pipeline}.dat**: Genotype-to-phenotype mapping for different pipelines. Pipelines means the combination of read aligner and variant caller.
Different secondary analysis pipelines may produce varying VCF file entries. However, it is of critical importance that the representation of ISBT variants in the VCF is correctly annotated. Currently, the differences are limited to the representation of the 109bp insertion of RHCE*02, but additional discrepancies cannot be ruled out.

Pipeline settings:
- **HGDP**    The original HGDP Project secondary analysis pipeline
- **TGSGATK** For third generation sequencing using pbmm2 and GATK
- **TGSPBSV** For third generation sequencing using pbmm2 and assuming we use pbsv for detecting the RHCE 109bp insertion
- **Dragen**  For data coming out of the Dragen platform
- **TGS**     For third generation sequencing using pbmm2 and deepVariant
Our recommendation is TGSGATK where possible ...
All differences in these annotation files currently stem from the varying ways in which the 109bp insertion in RHCE is represented by different variant callers in the VCF files. One exception is HGDP, where all alleles for RHCE have been removed and replaced with the antigens C/c and E/e. In this case, identification of the C/c antigen requires the output from detect_RHCplusminus.py, while the E/e antigen is determined conventionally via the tagging SNV 676G>C.

This structure is expected to change in the near future, as the current setup results in highly redundant annotation storage and requires corrections to be applied separately for each individual file across multiple systems — an inefficient and impractical solution.

### Testdata
Data for testing can be found under ./data/testdata/. HGDP samples 001, 002 and 003. The complete HGDP dataset used for benchmarking in our original publication can be downloaded at: https://www.internationalgenome.org/data-portal/data-collection/hgdp

## Cosine Similarity Scoring
bloodAGENT uses **cosine similarity** to measure the similarity between observed haplotypes and reference haplotypes from the International Society of Blood Transfusion (ISBT). The score ranges from **0 to 2**, where:
- **1 per haplotype** is the theoretical maximum.
- **2** is the best possible match for a diploid genome.

However, due to varying numbers of relevant SNPs across blood groups and alleles, **scores between different blood groups or individuals are not directly comparable**. A score of **1.9 vs. 1.8 does not necessarily indicate a better result** unless both results refer to the same blood group and allele.

## Running bloodAGENT
### Job Type: Phenotype Analysis
A typical command:
```sh
bloodAGENT --job phenotype \
  --target ./data/config/exonic_annotation.hg38.BGStarget.txt \
  --variants ./data/config/variation_annotation_GATK.dat \
  --gt2pt ./data/config/genotype_to_phenotype_annotation_GATK.dat \
  --vcf ./data/testdata/HGDP00001/HGDP00001.phased.vcf.gz \
  --bigwig ./data/testdata/HGDP00001/HGDP00001.BGStarget.bw \
  --coverage 12 --verbose 2 --scoreRange 1 \
  --out HGDP00001.json \
  --build hg38 -k --id "HGDP00001"
```
```sh
### Singularity:
singularity exec bloodagent.sif /app/bloodAGENT --job phenotype \
  --target ./data/config/exonic_annotation.hg38.BGStarget.txt \
  --variants ./data/config/variation_annotation_GATK.dat \
  --gt2pt ./data/config/genotype_to_phenotype_annotation_GATK.dat \
  --vcf ./data/testdata/HGDP00001/HGDP00001.phased.vcf.gz \
  --bigwig ./data/testdata/HGDP00001/HGDP00001.BGStarget.bw \
  --coverage 12 --verbose 2 --scoreRange 1 \
  --out HGDP00001.json \
  --build hg38 -k --id "HGDP00001"
```

#### Parameters for `phenotype` Job
## Command-Line Parameters

| Short Code | Long Code | Description | Data Type | Required | Default Value |
|------------|------------|--------------|----------|----------|--------------|
| `-j` | `--job phenotype` | Runs bloodAGENT to determine blood group phenotypes. | String | Yes | - |
| `-t` | `--target <file>` | Annotation file containing transcript information for blood group targets. | File | Yes | - |
| `-s` | `--variants <file>` | Variant annotation file for ISBT blood group typing. | File | Yes | - |
| `-g` | `--gt2pt <file>` | Mapping file from genotype to phenotype. | File | Yes | - |
| `-v` | `--vcf <file>` | VCF file containing phased genetic variants. | File | Yes | - |
| `-b` | `--bigwig <file>` | BigWig file for coverage data. | File | No | - |
| `-c` | `--coverage <int>` | Minimum sequencing coverage required for reliable results. | Integer | No | `10` |
| `-d` | `--verbose <int>` | Level of verbosity (0: none, 1: warnings, 2: status, 3: detailed logs). | Integer | No | `1` |
| `-r` | `--scoreRange <float>` | Score threshold multiplier for reporting matches. | Float | No | - |
| `-o` | `--out <file>` | Output file in JSON format. | File | No | "bloodAGENT.json" |
| `-u` | `--build <hg19|hg38>` | Specifies genome reference build. | String | Yes | - |
| `-k` | -trick | Enables coverage-based typing of RhD instead of variant-based typing. | Boolean (Flag) | No | `false` |
| `-f` | `--id <string>` | Sample identifier. | String | No | `unknown` |
>>>>>>> develop


### Job Type: Simulated Data Generation
A typical command:
```sh
bloodAGENT --job vcf \
  --variants ./data/config/variation_annotation_TGS.dat \
  --gt2pt ./data/config/genotype_to_phenotype_annotation_TGS.dat -a "ABO*A1.01" -b "ABO*O.01.01" \
  --phased \
  --dropout 1 --crack 5
```
```sh 
### Singularity:
singularity exec bloodagent.sif /app/bloodAGENT --job vcf \
  --variants ./data/config/variation_annotation_TGS.dat \
  --gt2pt ./data/config/genotype_to_phenotype_annotation_TGS.dat -a "ABO*A1.01" -b "ABO*O.01.01" \
  --phased \
  --dropout 1 --crack 5
```

#### Parameters for `vcf` Job
## Command-Line Parameters
| Short Code | Long Code | Description | Data Type | Required | Default Value |
|------------|------------|--------------|----------|----------|--------------|
| `-j` | `--job vcf` | Generates simulated genetic data in VCF format. | String | Yes | - |
| `-s` | `--variants <file>` | Variant annotation file for TGS-based analysis. | File | Yes | - |
| `-g` | `--gt2pt <file>` | Genotype-to-phenotype mapping file. | File | Yes | - |
| `-a` | `--alleleA <string>` | First allele for in silico simulation. | String | Yes | - |
| `-b` | `--alleleB <string>` | Second allele for in silico simulation. | String | Yes | - |
| `-p` | `--phased` | Ensures output includes phased haplotypes. | Boolean (Flag) | No | `false` |
| `-o` | `--dropout <float>` | Probability of SNP dropout (0.0–1.0). | Float | No | - |
| `-x` | `--crack <float>` | Probability of haplotype breakage at heterozygous sites (0.0–1.0). | Float | No | - |


## Output Format
bloodAGENT produces results in JSON format
For simulated data (`--job vcf`), the output is the data lines of an vcf file written to stdout.



### JSON File Structure

#### General
- **genome**: Specifies the genome build (e.g., hg38).
- **sample_id**: Specifies the sample id given by the user.
- **version**: bloodAGENT version

#### Parameter Section
- **command line parameters** List of all command line parameters and their values.

#### Data Section
- **loci**: Contains blood group systems.
  - **System Name (e.g., ABO)**:
    - **calls**: List of genotype and phenotype determinations.
      - **alleles**: List of detected alleles.
        - **names**: Names of identified alleles.
        - **issues**: Quality and coverage issues related to each allele.
    - **haplotypes**: Genotype data.
    - **phenotypes**: Predicted blood group phenotypes.
    - **score**: Cosine similarity score: If ISBT-relevant SNV positions do not meet the coverage requirements (as defined by the --coverage parameter), this score is set to zero.
    - **weak_score**: Cosine similarity score: This score remains unaffected if ISBT-relevant SNV positions do not meet the coverage requirements (as defined by the --coverage parameter). It represents the default score, disregarding coverage-failed variants.
    - **coverage_failed_variants**: List of genetic variations with insufficient coverage.
    - **mean_coverage**: Coverage statistics for different genomic regions.
    - **relevant_variations**: List of all ISBT variants of this blood group system.
      - **chrom**: Chromosome location.
      - **position**: Position in the genome.
      - **reference / alternative**: Reference and detected allele.
      - **high_impact**: Whether the variation is of high impact.
      - **is_covered**: Whether the variant has sufficient read coverage.
      - **depth**: Sequencing depth at the variant position.

This structure provides an overview of the key components within the JSON file, organizing metadata and data into a hierarchical format.


```

## Limitations
- **Dropout effects**: Missing variants significantly impact allele determination. Accuracy drops below **50% at a 50% dropout rate**.
- **Phasing information**: While its effect on ambiguity is minor, it remains important for resolving certain blood group systems (e.g., KEL, ABO, Duffy).
- **Paralogous regions**: Some blood group alleles (e.g., RHCE) may be misclassified due to read alignment issues, variant calling issues, annotation issues or any issue we are not aware of

```

## Licensing

Third party licenses can be found at Third_Party_Licenses.md
The docker/singularity container(s) include third-party components under various open source licenses.
See the `/licenses` directory inside the image for details.

The source code of this project is available at: 

The official [GitHub repository](https://github.com/ikmb/bloodAGENT).


