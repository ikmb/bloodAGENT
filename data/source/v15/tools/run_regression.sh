#!/usr/bin/env bash
# Run bloodAGENT (V15 image) against the 4 bundled test samples and write
# the V15 phenotype JSONs into data/source/v15/derived/regression/.
#
# Usage:
#   data/source/v15/tools/run_regression.sh [docker_image]
#
# Default image: bloodagent:v15
set -eo pipefail

IMG="${1:-bloodagent:v15}"
HERE="$(cd "$(dirname "$0")" && pwd)"
ROOT="$(cd "${HERE}/../../../.." && pwd)"
OUT="${ROOT}/data/source/v15/derived/regression"
mkdir -p "${OUT}"

run_sample() {
  local sample=$1
  local vcf=$2
  echo "=== ${sample} ==="
  docker run --rm \
    -v "${ROOT}":/work \
    -w /work \
    "${IMG}" \
      --job phenotype \
      --target  /work/data/config/exonic_annotation.hg38.BGStarget.txt \
      --variants /work/data/config/HGDP/variation_annotation_HGDP.dat \
      --gt2pt   /work/data/config/HGDP/genotype_to_phenotype_annotation_HGDP.dat \
      --vcf     "/work/${vcf}" \
      --bigwig  "/work/data/testdata/${sample}/${sample}.BGStarget.bw" \
      --coverage 12 --verbose 2 --scoreRange 1 \
      --out "/work/data/source/v15/derived/regression/${sample}.v15.json" \
      --build hg38 -k --id "${sample}" \
      || echo "  !! ${sample} FAILED with exit $?"
}

run_sample HGDP00001 "data/testdata/HGDP00001/HGDP00001.phased.vcf.gz,data/testdata/HGDP00001/HGDP00001.RHC.vcf"
run_sample HGDP00003 "data/testdata/HGDP00003/HGDP00003.phased.vcf.gz,data/testdata/HGDP00003/HGDP00003.RHC.vcf"
run_sample HGDP00005 "data/testdata/HGDP00005/HGDP00005.phased.vcf.gz,data/testdata/HGDP00005/HGDP00005.RHC.vcf"
run_sample NA24143   "data/testdata/NA24143/NA24143.GATK.phased.vcf.gz"

echo
echo "Regression outputs:"
ls -la "${OUT}"
