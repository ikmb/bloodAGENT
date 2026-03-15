#!/usr/bin/env bash
set -euo pipefail

# Build bloodAGENT and extract pre-built binaries + config for platform-services.
#
# Usage:
#   ./build_for_platform.sh <path-to-platform-services>
#
# Example:
#   ./build_for_platform.sh ~/Desktop/platform-services
#
# This builds the Docker image, then copies only the necessary artifacts
# into the processor package: binary, shared libs, Dragen config, and licenses.

PLATFORM_DIR="${1:?Usage: $0 <path-to-platform-services>}"
PROCESSOR_DIR="${PLATFORM_DIR}/packages/genomics_blood_type_processor"

if [[ ! -d "$PROCESSOR_DIR" ]]; then
    echo "ERROR: Processor package not found at ${PROCESSOR_DIR}" >&2
    exit 1
fi

IMAGE_NAME="bloodagent-builder:local"
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"

echo "Building bloodAGENT Docker image..."
docker build -t "${IMAGE_NAME}" "${SCRIPT_DIR}"

echo "Extracting artifacts..."
CONTAINER_ID=$(docker create "${IMAGE_NAME}")

# Clean previous artifacts
rm -rf "${PROCESSOR_DIR}/bin" "${PROCESSOR_DIR}/bloodagent_data" "${PROCESSOR_DIR}/licenses"

# Binary and shared libraries
mkdir -p "${PROCESSOR_DIR}/bin"
docker cp "${CONTAINER_ID}:/usr/local/bin/bloodAGENT" "${PROCESSOR_DIR}/bin/bloodAGENT"
docker cp "${CONTAINER_ID}:/usr/local/lib/libhts.so" "${PROCESSOR_DIR}/bin/libhts.so"
docker cp "${CONTAINER_ID}:/usr/local/lib/libBigWig.so" "${PROCESSOR_DIR}/bin/libBigWig.so"

# Only the config files we need (Dragen pipeline + shared hg38 target annotation)
mkdir -p "${PROCESSOR_DIR}/bloodagent_data/config/Dragen"
docker cp "${CONTAINER_ID}:/data/config/exonic_annotation.hg38.BGStarget.txt" \
    "${PROCESSOR_DIR}/bloodagent_data/config/exonic_annotation.hg38.BGStarget.txt"
docker cp "${CONTAINER_ID}:/data/config/Dragen/variation_annotation_Dragen.dat" \
    "${PROCESSOR_DIR}/bloodagent_data/config/Dragen/variation_annotation_Dragen.dat"
docker cp "${CONTAINER_ID}:/data/config/Dragen/genotype_to_phenotype_annotation_Dragen.dat" \
    "${PROCESSOR_DIR}/bloodagent_data/config/Dragen/genotype_to_phenotype_annotation_Dragen.dat"

# Licenses (BSD 2-Clause compliance)
mkdir -p "${PROCESSOR_DIR}/licenses"
docker cp "${CONTAINER_ID}:/licenses/bloodAGENT-LICENSE" "${PROCESSOR_DIR}/licenses/bloodAGENT-LICENSE"
docker cp "${CONTAINER_ID}:/licenses/bloodAGENT-Third_Party_Licenses.md" "${PROCESSOR_DIR}/licenses/bloodAGENT-Third_Party_Licenses.md"

docker rm "${CONTAINER_ID}" > /dev/null

echo ""
echo "Done. Files extracted to ${PROCESSOR_DIR}:"
find "${PROCESSOR_DIR}/bin" "${PROCESSOR_DIR}/bloodagent_data" "${PROCESSOR_DIR}/licenses" -type f | sort | sed "s|${PROCESSOR_DIR}/|  |"
