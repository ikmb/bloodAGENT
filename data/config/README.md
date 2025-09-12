# bloodAGENT - Pipelines

For the analysis of different datasets, several configuration files are provided. Currently, there are four groups:

- **TGS - PacBio**
- **ONT - Oxford Nanopore**
- **HGDP - The publicly available HGDP dataset as analyzed in Bergström A. et al. 2020**
- **Dragen - NGS data generated with the Illumina Dragen pipeline**

For each analysis, two groups of configuration files must be loaded:

- `variation_annotation.dat` — contains the list of ISBT-relevant variants (`--variants`)
- `genotype_to_phenotype_annotation.dat` — contains the list of ISBT alleles with their haplotypes (`--gt2pt`)

Different analysis pipelines (bwa, minimap2, GATK, DeepVariant, Clair3, Sniffles, etc.) sometimes report certain variants differently. This applies especially to the RHCE*C/c-relevant 109 bp insertion. To address this, pipeline-specific annotation files exist that must be loaded in addition to the main annotation. These are also provided under the parameters `--variants` or `--gt2pt`, but specified as a comma-separated list.  
Example:  
--variants /path/variation_annotation.dat,/path/variation_annotation_GATK.dat

---

## Long-read analyses (PacBio and ONT)

### PacBio
- **GATK secondary analysis**  
```sh
--variants bloodAGENT/data/config/variation_annotation.dat,bloodAGENT/data/config/PacBio/variation_annotation_TGSGATK.dat
--gt2pt bloodAGENT/data/config/genotype_to_phenotype_annotation.dat
```
example:
```sh
# assuming you are at the bllodAGENT root directory
./dist/Release/GNU-Linux/bloodAGENT  --job phenotype \
  --target ./data/config/exonic_annotation.hg38.BGStarget.txt \
  --variants ./data/config/variation_annotation.dat,./data/config/PacBio/variation_annotation_TGSGATK.dat \
  --gt2pt ./data/config/genotype_to_phenotype_annotation.dat \
  --vcf ./data/testdata/NA24143/NA24143.GATK.phased.vcf.gz \
  --bigwig ./data/testdata/NA24143/NA24143.BGStarget.bw \
  --coverage 12 --verbose 2 --scoreRange 1 \
  --out PacBioSample.json \
  --build hg38 -k --id "PacBioSample"
```
- **DeepVariant secondary analysis with pbsv for larger insertions/deletions**  
```sh
--variants bloodAGENT/data/config/variation_annotation.dat,bloodAGENT/data/config/PacBio/variation_annotation_TGSPBSV.dat
--gt2pt bloodAGENT/data/config/genotype_to_phenotype_annotation.dat
```

### Oxford Nanopore (ONT)
We analyzed the publicly available dataset using minimap2, Clair3, WhatsHap, and Sniffles2.  
```sh
--variants bloodAGENT/data/config/variation_annotation.dat,bloodAGENT/data/config/ONT/variation_annotation_MINIMAP2SNIFFLES.dat
--gt2pt bloodAGENT/data/config/genotype_to_phenotype_annotation.dat
```
example:
```sh
# currently just an example, real testdata will we included soon
bloodAGENT --job phenotype \
  --target ./data/config/exonic_annotation.hg38.BGStarget.txt \
  --variants ./data/config/variation_annotation.dat,./data/config/ONT/variation_annotation_MINIMAP2SNIFFLES.dat \
  --gt2pt ./data/config/genotype_to_phenotype_annotation.dat \
  --vcf ./your/data/Sample.phased.vcf.gz \
  --bigwig ./your/data/Sample.BGStarget.bw \
  --coverage 12 --verbose 2 --scoreRange 1 \
  --out Sample.json \
  --build hg38 -k --id "Sample"
```

## Short-read NGS analyses

NGS data often struggles to reliably detect insertions/deletions, especially the critical 109 bp RHCE insertion. Therefore, different strategies are required:

### HGDP
This dataset is publicly available and already processed by secondary analysis:  
Bergström A, et al. *Insights into human genetic variation and population history from 929 diverse genomes.* **Science. 2020;367:eaay5012.**

The 109 bp insertion is not reliably detected. To compensate, we use the script `detect_RHCplusminus.py`, which determines the insertion based on RHCE exon 2 coverage. The script generates a VCF file that must be provided to bloodAGENT in addition to the main VCF file (comma-separated under the `--vcf` parameter).  
The required `variation_annotation` and `genotype_to_phenotype_annotation` files are in:  bloodAGENT/data/config/HGDP
example:
```sh
# assuming you are at the bllodAGENT root directory
./dist/Release/GNU-Linux/bloodAGENT --job phenotype \
  --target ./data/config/exonic_annotation.hg38.BGStarget.txt \
  --variants ./data/config/HGDP/variation_annotation_HGDP.dat \
  --gt2pt ./data/config/HGDP/genotype_to_phenotype_annotation_HGDP.dat \
  --vcf ./data/testdata/HGDP00001/HGDP00001.phased.vcf.gz,./data/testdata/HGDP00001/HGDP00001.RHC.vcf \
  --bigwig ./data/testdata/HGDP00001/HGDP00001.BGStarget.bw \
  --coverage 12 --verbose 2 --scoreRange 1 \
  --out Sample.json \
  --build hg38 -k --id "Sample"
```


### Dragen
The Dragen pipeline can detect the RHCE 109 bp insertion, but RHD+/– status cannot be reliably determined via coverage. Instead, CNV analysis results can be used.  
The script `parse_cnv_vcf.py` must be applied to the CNV VCF file from the secondary analysis. The resulting re-encoded VCF file must be provided to bloodAGENT in addition to the genome-wide regular VCF (comma-separated under the `--vcf` parameter).  
The required `variation_annotation` and `genotype_to_phenotype_annotation` files are in: bloodAGENT/data/config/Dragen
example:
```sh
# currently just an example, real testdata will we included soon
bloodAGENT --job phenotype \
  --target ./data/config/exonic_annotation.hg38.BGStarget.txt \
  --variants bloodAGENT/data/config/Dragen/variation_annotation_Dragen.dat \
  --gt2pt bloodAGENT/data/config/Dragen/genotype_to_phenotype_annotation_Dragen.dat \
  --vcf ./your/data/Sample.phased.vcf.gz,./your/data/Sample.cnv.recode.vcf \
  --bigwig ./your/data/Sample.BGStarget.bw \
  --coverage 12 --verbose 2 --scoreRange 1 \
  --out Sample.json \
  --build hg38 -k --id "Sample"
```
