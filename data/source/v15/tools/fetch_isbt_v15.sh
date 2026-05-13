#!/usr/bin/env bash
# Fetch full ISBT Blood Group Database V15 snapshot via the public REST API.
#
# Usage: ./fetch_isbt_v15.sh [RAW_DIR]
#   RAW_DIR defaults to ../raw relative to this script.
#
# What it does:
#   1. Dumps the 7 list endpoints (system, gene, antigen, allele, variant, phenotype, release)
#   2. Dumps /api/release/15 for the V15 release object incl. updatedAlleleIds
#   3. For every allele id in allele.json, dumps /api/allele/<id> into alleles/<id>.json
#      (this is required because the bulk allele endpoint does NOT include variants[])
#   4. Retries throttled responses (HTTP 429 / "Too Many Requests") with backoff
#
# Note: The API enforces a per-IP throttler. Concurrency is intentionally low.
set -euo pipefail

API="https://blooddatabase.isbtweb.org/api"
HERE="$(cd "$(dirname "$0")" && pwd)"
RAW_DIR="${1:-${HERE}/../raw}"
mkdir -p "${RAW_DIR}/alleles"

echo "[1/3] fetching list endpoints ..."
for ep in system gene antigen allele variant phenotype release publication genbank; do
  echo "  - ${ep}"
  curl -sL --retry 3 --retry-delay 2 "${API}/${ep}" -o "${RAW_DIR}/${ep}.json"
done

echo "[2/3] fetching release/15 ..."
curl -sL --retry 3 "${API}/release/15" -o "${RAW_DIR}/release_15.json"

echo "[3/3] fetching per-allele detail (with variants[]) ..."
python3 -c "
import json
print('\n'.join(str(x['id']) for x in json.load(open('${RAW_DIR}/allele.json'))))
" > "${RAW_DIR}/allele_ids.txt"

fetch_one() {
  local id=$1
  curl -sL "${API}/allele/${id}" -o "${RAW_DIR}/alleles/${id}.json"
}
export -f fetch_one
export RAW_DIR API

# First pass: 4-way parallel with small delay
xargs -P 4 -I {} bash -c 'fetch_one {}; sleep 0.15' < "${RAW_DIR}/allele_ids.txt"

# Retry pass: re-fetch anything <200 bytes (rate-limited "429" stubs)
for pass in 1 2 3; do
  count=$(find "${RAW_DIR}/alleles" -name '*.json' -size -200c | wc -l | tr -d ' ')
  echo "  retry pass ${pass}: ${count} files still throttled"
  [ "$count" -eq 0 ] && break
  find "${RAW_DIR}/alleles" -name '*.json' -size -200c -print0 \
    | xargs -0 -n1 basename \
    | sed 's/\.json//' \
    | while read id; do
        curl -sL "${API}/allele/${id}" -o "${RAW_DIR}/alleles/${id}.json"
        sleep 0.6
      done
done

final=$(find "${RAW_DIR}/alleles" -name '*.json' -size -200c | wc -l | tr -d ' ')
echo "Done. ${final} files still throttled (try again later if > 0)."
