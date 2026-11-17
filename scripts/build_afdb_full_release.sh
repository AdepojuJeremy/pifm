#!/usr/bin/env bash
# Full weak-label release: no taxon restriction, no sample cap.
set -euo pipefail
cd "$(dirname "$0")/.."
raw=data/raw/afdb
dist="$raw/full_release"
mkdir -p "$raw" "$dist"
if [[ $# -eq 0 ]]; then
  echo "Usage: $0 RELEASE_MATCHED_FASTA [RELEASE_MATCHED_FASTA ...]" >&2
  echo "Provide UniProt 2025_04 FASTA files covering all metadata accessions." >&2
  exit 2
fi
if [[ -e "$dist/COMPLETE" || -e "$raw/afdb_full.csv" || -e "$raw/afdb_full.audit.json" ]]; then
  echo "An existing build is present; move it before rebuilding." >&2
  exit 2
fi
metadata="$raw/heterodimer_metadata.csv"
url=https://ftp.ebi.ac.uk/pub/databases/alphafold/collaborations/nvda/heterodimer_metadata.csv
if [[ ! -s "$metadata" ]]; then
  curl --fail --location --retry 5 --continue-at - "$url" -o "$metadata.part"
  mv "$metadata.part" "$metadata"
fi
args=()
for fasta in "$@"; do
  [[ -s "$fasta" ]] || { echo "Missing FASTA: $fasta" >&2; exit 2; }
  args+=(--fasta "$fasta")
done
python -m unittest discover -s tests -p 'test_prepare_afdb_pseudolabels.py'
python scripts/prepare_afdb_pseudolabels.py \
  --metadata "$metadata" "${args[@]}" \
  --sequence-release 2025_04 --min-length 15 --max-length 1500 \
  --positive-min 0.8 --negative-max 0.4 --max-per-class 0 \
  --out "$raw/afdb_full.csv"
sha256sum "$metadata" "$@" > "$dist/input-checksums.sha256"
cp "$raw/afdb_full.audit.json" "$dist/"
git rev-parse HEAD > "$dist/source-commit.txt"
cp docs/AFDB_FULL_RELEASE.md "$dist/DATASET_CARD.md"
# Split by CSV rows before compression. Each shard has its own header,
# so each can be independently decompressed and ingested.
python scripts/shard_afdb_release.py "$raw/afdb_full.csv" "$dist"
(
  cd "$dist"
  sha256sum *.csv.gz *.json *.md *.txt input-checksums.sha256 > SHA256SUMS
)
touch "$dist/COMPLETE"
echo "Full export: $dist"
