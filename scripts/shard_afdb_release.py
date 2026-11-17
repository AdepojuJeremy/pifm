"""Export independently readable compressed shards of a full AFDB CSV."""
import csv
import gzip
import json
import sys
from collections import Counter
from pathlib import Path

def shard(source, dest, shard_rows=100_000):
    dest = Path(dest)
    if list(dest.glob("afdb-full-*.csv.gz")):
        raise FileExistsError("Existing shards; use an empty release directory.")
    counts = Counter()
    shards = []
    out = None
    try:
        with Path(source).open(newline="", encoding="utf-8") as handle:
            reader = csv.reader(handle)
            header = next(reader)
            label_index = header.index("interaction_label")
            for i, row in enumerate(reader):
                if len(row) != len(header) or row[label_index] not in ("0", "1"):
                    raise ValueError(f"Invalid row {i+2}")
                if i % shard_rows == 0:
                    if out is not None:
                        out.close()
                    name = f"afdb-full-{len(shards):04d}.csv.gz"
                    out = gzip.open(dest / name, "wt", encoding="utf-8", newline="", compresslevel=5)
                    writer = csv.writer(out)
                    writer.writerow(header)
                    shards.append({"file": name, "rows": 0})
                writer.writerow(row)
                shards[-1]["rows"] += 1
                counts[row[label_index]] += 1
    finally:
        if out is not None:
            out.close()
    audit = json.loads((dest / "afdb_full.audit.json").read_text())["counts"]
    for label in ("0", "1"):
        if counts[label] != audit[f"exported_label_{label}"]:
            raise ValueError("Export row counts do not match preprocessing audit")
    manifest = {"rows": sum(counts.values()), "label_counts": dict(counts),
                "shards": shards, "schema": header, "label_semantics": "weak supervision"}
    (dest / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    return manifest

if __name__ == "__main__":
    print(json.dumps(shard(sys.argv[1], sys.argv[2]), indent=2))
