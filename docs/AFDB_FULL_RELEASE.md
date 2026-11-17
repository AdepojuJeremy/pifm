# SnapPPI AFDB full heterodimer weak-label dataset

Scope: all organisms in the NVIDIA/AFDB heterodimer metadata. No sample cap or
species restriction. This dataset is derived from computational structure
predictions, not experimental PPI labels. It is not an affinity dataset.
The number of released pairs is established by the actual build manifest.

## Source and sequence provenance

- Metadata: https://ftp.ebi.ac.uk/pub/databases/alphafold/collaborations/nvda/
- AFDB release description: https://www.alphafold.ebi.ac.uk/faq
- UniProt release 2025_04, matching the NVIDIA source release.
  Historical archives: https://ftp.uniprot.org/pub/databases/uniprot/previous_releases/release-2025_04/
- Supply release-matched FASTA covering all metadata accessions. Do not use
  the latest UniProt release as a silent fallback; accession mapping alone
  does not prove the sequence matches the modeled chain.

Thresholds are configurable research choices: ipTM >= 0.8 yields candidate
positives; ipTM <= 0.4 yields candidate negatives. Reject intermediate scores,
invalid scores, unresolved or ambiguous sequences, identical chains and
discordant repeated sequence pairs. Chain lengths are 15–1500 residues
(the source paper's monomer range); standard amino-acid alphabet only.
Low ipTM means low interface confidence, not experimentally proven absence
of binding. All predicted labels must remain distinct from curated labels.

## Full build

On a network-enabled machine with adequate disk and RAM, run from the
repository root:

    git switch research/afdb-pseudolabels
    bash scripts/build_afdb_full_release.sh

The no-argument build attempts to discover and stream the official UniProt
2025_04 archives, retaining only the accessions in AFDB metadata. This may
require substantial network traffic. Archive layout and accessibility have
not yet been validated in a live run; discovery fails explicitly if it
cannot identify the correct archives. Alternatively pass one or more
verified release-matched local FASTA files as arguments.
The builder automatically retrieves the complete metadata CSV, runs the
offline tests, processes the input without a row cap, and packages the
result. It does not silently fetch current UniProt sequences.

Output: data/raw/afdb/full_release/
- afdb-full-*.csv.gz: independent CSV shards containing sequences, labels,
  accession identifiers, representative model ID and ipTM ranges.
- afdb_full.audit.json: filtering and class counts, accession coverage.
- manifest.json: shard inventory, exact row totals and schema.
- sequence-source.json and unresolved-accessions.txt: archive provenance and
  unmatched accessions when automatic sequence resolution was used.
- SHA256SUMS, input-checksums.sha256, source-commit.txt, DATASET_CARD.md.
- COMPLETE: written only after successful packaging.

Check the accession coverage audit and verify model-chain sequences before
publication. A successful script run alone does not establish completeness
of the sequence mapping. The source metadata includes predictions
preselected from known or putative associations, so low scores are not a
random negative sample. Source selection bias must be considered.

## Hosting

After successful real-data processing and review of the audit, publish from
a GitHub-authenticated terminal (requires release creation permission):

    bash scripts/publish_afdb_release.sh

The publication script uploads to a draft release, then publishes only
after all assets have uploaded. The actual counts and checksums travel with
the data. Do not create a public release claiming a dataset exists until
the assets are present.

## Training

The files are candidates for training ONLY. Use a separately frozen
experimental validation/test set. Exclude any pseudo-labeled pair containing
held-out proteins or homologous sequence families (at least MMseqs2 50%
identity / 80% coverage), and resolve collisions in favor of experimental
labels. The reference holdout and its actual sequences are not in this
fork, so this audit is a separate prerequisite to model integration.
Use the current upstream SnapPPI classifier rather than the older
multitask fork. Do not use these scores as experimentally validated labels.

## Attribution

Source data: AlphaFold Protein Structure Database, EMBL-EBI, Google
DeepMind, NVIDIA, and the Steinegger Lab; Han et al., AlphaFold Database
expands to proteome-scale quaternary structures (2026). UniProt Consortium,
UniProt 2025_04. Derived dataset CC BY 4.0:
https://creativecommons.org/licenses/by/4.0/
