#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
dist=data/raw/afdb/full_release
tag=afdb-heterodimers-iptm-v1
repo=AdepojuJeremy/pifm
[[ -f "$dist/COMPLETE" ]] || { echo "No completed dataset build." >&2; exit 1; }
(cd "$dist" && sha256sum -c SHA256SUMS)
# Do not overwrite a previous release. A failed upload remains a draft.
if gh release view "$tag" --repo "$repo" >/dev/null 2>&1; then
  echo "Release exists; inspect the existing release/draft before retrying." >&2
  exit 1
fi
gh release create "$tag" --repo "$repo" --target "$(git rev-parse HEAD)" \
  --draft --title "AFDB heterodimers — full ipTM weak-label dataset v1" \
  --notes-file "$dist/DATASET_CARD.md"
for asset in "$dist"/*.csv.gz "$dist"/*.json "$dist"/*.md "$dist"/*.txt \
  "$dist"/SHA256SUMS "$dist"/input-checksums.sha256; do
  gh release upload "$tag" "$asset" --repo "$repo"
done
gh release edit "$tag" --repo "$repo" --draft=false
gh release view "$tag" --repo "$repo" --json url,assets
