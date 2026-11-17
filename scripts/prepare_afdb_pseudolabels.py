#!/usr/bin/env python3
"""Prepare weak PPI labels from AFDB/NVIDIA heterodimer metadata and local FASTA.

This is a preprocessing experiment, not a modification to model training.
Uses only the Python standard library. Does not download structure files.
"""
import argparse
import csv
import gzip
import hashlib
import json
import math
import sqlite3
import tempfile
from collections import Counter
from pathlib import Path

REQUIRED = {"modelEntityId", "uniprot_ac_1", "uniprot_ac_2",
            "tax_id_1", "tax_id_2", "ipTM"}
AMINO_ACIDS = set("ACDEFGHIKLMNPQRSTVWY")


def open_text(path):
    path = Path(path)
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8-sig", newline="")
    return path.open(encoding="utf-8-sig", newline="")


def fasta_records(path):
    """Accept UniProt sp|ACC|NAME/tr|ACC|NAME or bare accession headers."""
    acc, parts = None, []
    with open_text(path) as handle:
        for line in handle:
            line = line.strip()
            if line.startswith(">"):
                if acc is not None:
                    yield acc, "".join(parts).upper()
                fields = line[1:].split()
                if not fields:
                    raise ValueError(f"Empty FASTA header in {path}")
                token = fields[0].split("|")
                acc = token[1] if len(token) >= 3 and token[0] in ("sp", "tr") else token[0]
                parts = []
            elif line:
                if acc is None:
                    raise ValueError(f"Sequence before FASTA header in {path}")
                parts.append(line)
        if acc is not None:
            yield acc, "".join(parts).upper()


def load_sequences(paths, min_length, max_length, wanted=None):
    sequences, ambiguous = {}, set()
    # Do not silently strip isoform suffixes or resolve accessions to different
    # full-length proteins. Conflicting sequence records are excluded.
    for path in paths:
        for acc, seq in fasta_records(path):
            if wanted is not None and acc not in wanted:
                continue
            if acc in ambiguous:
                continue
            if acc in sequences and sequences[acc] != seq:
                ambiguous.add(acc)
                del sequences[acc]
            else:
                sequences[acc] = seq
    return {acc: seq for acc, seq in sequences.items()
            if min_length <= len(seq) <= max_length and set(seq) <= AMINO_ACIDS}, len(ambiguous)


def pair_key(a, b, seed):
    left, right = sorted((a, b))
    return hashlib.sha256(f"{seed}\0{left}\0{right}".encode()).hexdigest()


def score_label(low, high, negative_max, positive_min):
    """Conservative agreement across every scored model of the sequence pair."""
    if low >= positive_min:
        return 1
    if high <= negative_max:
        return 0
    return None


def prepare(args):
    if not 0 <= args.negative_max < args.positive_min <= 1:
        raise ValueError("Require 0 <= negative-max < positive-min <= 1")
    if not 1 <= args.min_length <= args.max_length or args.max_per_class < 0:
        raise ValueError("Invalid length or class cap")
    output = Path(args.out)
    audit = output.with_suffix(".audit.json")
    if output.suffix.lower() != ".csv":
        raise ValueError("--out must end in .csv")
    if output.exists() or audit.exists():
        raise FileExistsError("Output or audit already exists; use a new filename.")
    # First pass limits sequence RAM to accessions actually in the metadata.
    wanted = set()
    with open_text(args.metadata) as handle:
        reader = csv.DictReader(handle)
        missing = REQUIRED - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"Missing AFDB columns: {sorted(missing)}")
        for row in reader:
            if args.taxid and (row['tax_id_1'] != args.taxid or row['tax_id_2'] != args.taxid):
                continue
            wanted.update((row['uniprot_ac_1'].strip(), row['uniprot_ac_2'].strip()))
    wanted.discard('')
    sequences, n_ambiguous = load_sequences(args.fasta, args.min_length, args.max_length, wanted)
    if not sequences:
        raise ValueError("No usable FASTA records.")
    counts = Counter()
    output.parent.mkdir(parents=True, exist_ok=True)
    # Disk-backed deduplication: retain score ranges, not replicated sequences.
    with tempfile.TemporaryDirectory(prefix="afdb-", dir=output.parent) as tmp:
        with sqlite3.connect(str(Path(tmp) / "pairs.sqlite")) as db:
            db.execute("""CREATE TABLE pairs (
                key TEXT PRIMARY KEY, acc1 TEXT NOT NULL, acc2 TEXT NOT NULL,
                model_id TEXT NOT NULL, score_min REAL NOT NULL, score_max REAL NOT NULL,
                models INTEGER NOT NULL) WITHOUT ROWID""")
            with open_text(args.metadata) as handle:
                reader = csv.DictReader(handle)
                missing = REQUIRED - set(reader.fieldnames or [])
                if missing:
                    raise ValueError(f"Missing AFDB columns: {sorted(missing)}")
                for row in reader:
                    counts["metadata_rows"] += 1
                    if args.taxid and (row["tax_id_1"] != args.taxid or row["tax_id_2"] != args.taxid):
                        counts["other_taxon"] += 1
                        continue
                    try:
                        score = float(row["ipTM"])
                    except (ValueError, TypeError):
                        counts["invalid_iptm"] += 1
                        continue
                    if not math.isfinite(score) or not 0 <= score <= 1:
                        counts["invalid_iptm"] += 1
                        continue
                    a, b = row["uniprot_ac_1"].strip(), row["uniprot_ac_2"].strip()
                    s1, s2 = sequences.get(a), sequences.get(b)
                    if s1 is None or s2 is None:
                        counts["unresolved_or_filtered_sequence"] += 1
                        continue
                    if a == b or s1 == s2:
                        counts["identical_chain"] += 1
                        continue
                    if s1 > s2:
                        a, b, s1, s2 = b, a, s2, s1
                    model = row["modelEntityId"].strip()
                    if not model:
                        counts["missing_model_id"] += 1
                        continue
                    db.execute("""INSERT INTO pairs VALUES (?, ?, ?, ?, ?, ?, 1)
                        ON CONFLICT(key) DO UPDATE SET
                        score_min = min(score_min, excluded.score_min),
                        score_max = max(score_max, excluded.score_max),
                        models = models + 1""",
                        (pair_key(s1, s2, args.seed), a, b, model, score, score))
                    counts["eligible_model_rows"] += 1
                    if counts["eligible_model_rows"] % 50000 == 0:
                        db.commit()
            db.commit()
            counts["unique_sequence_pairs"] = db.execute("SELECT count(*) FROM pairs").fetchone()[0]
            counts["discordant_high_low_pairs"] = db.execute(
                "SELECT count(*) FROM pairs WHERE score_min <= ? AND score_max >= ?",
                (args.negative_max, args.positive_min)).fetchone()[0]
            counts["positive_candidates"] = db.execute(
                "SELECT count(*) FROM pairs WHERE score_min >= ?", (args.positive_min,)).fetchone()[0]
            counts["negative_candidates"] = db.execute(
                "SELECT count(*) FROM pairs WHERE score_max <= ?", (args.negative_max,)).fetchone()[0]
            counts["excluded_uncertain_pairs"] = (counts["unique_sequence_pairs"]
                - counts["positive_candidates"] - counts["negative_candidates"])
            if not counts["positive_candidates"] or not counts["negative_candidates"]:
                raise ValueError(f"No candidates in at least one class: {dict(counts)}")
            staged = Path(tmp) / "dataset.csv"
            with staged.open("w", newline="", encoding="utf-8") as handle:
                writer = csv.writer(handle)
                writer.writerow(["seq1", "seq2", "interaction_label", "affinity_nm",
                                 "uniprot_ac_1", "uniprot_ac_2", "representative_model_id",
                                 "iptm_min", "iptm_max", "metadata_rows_for_pair", "label_origin"])
                for label, condition, threshold in [
                    (1, "score_min >= ?", args.positive_min),
                    (0, "score_max <= ?", args.negative_max),
                ]:
                    query = "SELECT acc1, acc2, model_id, score_min, score_max, models FROM pairs WHERE " + condition + " ORDER BY key LIMIT ?"
                    for a, b, model, lo, hi, n in db.execute(query, (threshold, args.max_per_class or -1)):
                        writer.writerow([sequences[a], sequences[b], label, "", a, b, model,
                                         lo, hi, n, "afdb_iptm_pseudolabel"])
                        counts[f"exported_label_{label}"] += 1
            report = {
                "source": "AFDB/NVIDIA heterodimers (weak labels, not experimental)",
                "source_url": "https://ftp.ebi.ac.uk/pub/databases/alphafold/collaborations/nvda/",
                "license": "CC BY 4.0",
                "metadata": str(args.metadata),
                "fasta": [str(p) for p in args.fasta],
                "sequence_release": args.sequence_release,
                "policy": {"taxid": args.taxid, "negative_max": args.negative_max,
                           "positive_min": args.positive_min, "min_length": args.min_length,
                           "max_length": args.max_length, "seed": args.seed,
                           "max_per_class": args.max_per_class,
                           "duplicate_rule": "unanimous threshold across valid scores"},
                "counts": dict(counts), "usable_accessions": len(sequences),
                "requested_accessions": len(wanted),
                "unresolved_or_filtered_accessions": len(wanted - sequences.keys()),
                "ambiguous_fasta_accessions": n_ambiguous,
                "limitations": ["Low ipTM does not establish non-interaction.",
                    "Accession-matched FASTA must correspond to modeled sequences.",
                    "No experimental holdout/homology audit has been performed."],
            }
            staged_audit = Path(tmp) / "audit.json"
            staged_audit.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
            staged.replace(output)
            staged_audit.replace(audit)
    print(json.dumps(report, indent=2))
    return report


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--metadata", type=Path, required=True)
    p.add_argument("--fasta", type=Path, action="append", required=True)
    p.add_argument("--sequence-release", required=True, help="Record the FASTA release/version.")
    p.add_argument("--out", type=Path, required=True)
    p.add_argument("--taxid", help="Require both chains to have this taxonomy ID.")
    p.add_argument("--positive-min", type=float, default=0.8)
    p.add_argument("--negative-max", type=float, default=0.4)
    p.add_argument("--min-length", type=int, default=15)
    p.add_argument("--max-length", type=int, default=1500)
    p.add_argument("--max-per-class", type=int, default=0, help="0 disables cap (full release).")
    p.add_argument("--seed", type=int, default=1)
    prepare(p.parse_args())


if __name__ == "__main__":
    main()
