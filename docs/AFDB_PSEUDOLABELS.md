# AFDB weak-label dataset pilot

This opt-in experiment prepares protein sequence pairs from the AFDB/NVIDIA
heterodimer release. It does not alter SnapPPI training or aggregation.
This fork's main branch predates the current classification-only upstream.

## Inputs

- Official metadata and schema:
  https://ftp.ebi.ac.uk/pub/databases/alphafold/collaborations/nvda/
- Release description: https://alphafold.ebi.ac.uk/faq
- The metadata CSV is approximately 2.3 GB. No structure archives are needed
  for initial accession-based preprocessing.
- Supply local UniProt FASTA files corresponding to the modeled chains.
  The AFDB release documentation identifies UniProt 2025_04; use matching
  sequence snapshots and verify exact sequence identity against a sample of
  modeled chains before treating this as a release dataset. Do not silently
  substitute current sequences or strip isoform suffixes.

## Run locally (Python 3.10+; standard library only)

    git clone https://github.com/AdepojuJeremy/pifm.git
    cd pifm
    git switch research/afdb-pseudolabels
    python -m unittest discover -s tests -p 'test_prepare_afdb_pseudolabels.py'

    mkdir -p data/raw/afdb
    curl -fL -C - \
      https://ftp.ebi.ac.uk/pub/databases/alphafold/collaborations/nvda/heterodimer_metadata.csv \
      -o data/raw/afdb/heterodimer_metadata.csv

Obtain a release-matched human proteome FASTA from UniProt (reference proteome
UP000005640), or provide your own verified accession-to-modeled-sequence
FASTA. UniProt downloads: https://www.uniprot.org/downloads .
Save it as data/raw/afdb/human_2025_04.fasta.gz (or adjust the command and
release argument to describe the actual file). Then:

    python scripts/prepare_afdb_pseudolabels.py \
      --metadata data/raw/afdb/heterodimer_metadata.csv \
      --fasta data/raw/afdb/human_2025_04.fasta.gz \
      --sequence-release 2025_04 \
      --taxid 9606 \
      --positive-min 0.8 --negative-max 0.4 \
      --max-per-class 10000 --seed 1 \
      --out data/raw/afdb/human_pilot.csv

Repeat --fasta for multiple proteomes. The input sequences are loaded into
memory, while pair deduplication is disk-backed. Start with one proteome;
allow spare disk space for the source CSV and temporary SQLite database.
The class cap is applied after deduplication and is a maximum, not a promise
of a balanced class count. Output selection is deterministic by sequence-pair
hash and seed. For a small smoke test use local fixture inputs, not a partial
production download.

## Outputs

- human_pilot.csv: seq1, seq2, interaction_label, blank affinity_nm,
  accessions, representative model ID, ipTM range, model row count and
  weak-label provenance.
- human_pilot.audit.json: counts, policy, exclusions, and source information.

Record input checksums (for example, sha256sum on Linux), FASTA release,
script commit, thresholds and the audit report along with the final export.
Keep raw data and exports untracked; do not commit them to Git.

## Label semantics and next steps

1. ipTM >= 0.8 is pseudo-positive; ipTM <= 0.4 is candidate pseudo-negative.
   Intermediate and inconsistent repeated predictions are excluded. These
   thresholds are experimental, not established non-interaction criteria.
2. Canonicalize sequence pairs irrespective of partner order; do not label
   identical chains in this heterodimer pilot. Conflicting FASTA records,
   unresolved accessions, invalid scores and unsupported sequence characters
   are excluded.
3. Before combining with training, compare labels to experimental positives
   and negatives. Consider stricter ipSAE/pDockQ-based thresholds, check
   disorder and fold confidence, and examine species/length/degree bias.
4. Freeze an independent experimental validation/test split. Remove any
   pseudo-label pair containing a held-out sequence family via a joint
   homology audit (at least the upstream MMseqs2 50% identity/80% coverage
   policy). Exact deduplication alone is insufficient. Keep experimental
   labels authoritative when there is disagreement.
5. Sync/rebase with the current upstream before adding a dedicated optional
   afdb_pseudolabel source loader. Keep these rows out of validation/test.
   Compare experimental-only control with pseudo-positive-only and
   positive+negative weak supervision using the same holdout and seeds.
   Do not just drop the CSV into the literature source: this would collapse
   predicted and curated sources into one label and risk split leakage.

This pilot does not prove that low-confidence complexes do not bind and
does not estimate binding affinity. It is an AFDB teacher-distillation
experiment, not a new experimental interaction database.

Credit: AlphaFold Database, EMBL-EBI, Google DeepMind, NVIDIA and the
Steinegger Lab; cite Han et al., *AlphaFold Database Expands to Proteome-Scale
Quaternary Structures* (2026), DOI 10.64898/2026.03.27.714458, and the
relevant UniProt release. AFDB data are distributed under CC BY 4.0.
