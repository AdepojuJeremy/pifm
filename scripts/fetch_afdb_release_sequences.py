#!/usr/bin/env python3
"""Extract only AFDB-needed FASTA sequences from UniProt 2025_04 archives.

Archive discovery is fail-closed: never substitutes current_release. Full
UniProt archives can be very large. Run on a network-enabled data server;
streaming avoids storing the full archives on disk.
"""
import argparse
import csv
import gzip
import io
import json
import re
import shutil
import sys
import tarfile
import urllib.parse
import urllib.request
from pathlib import Path
from html.parser import HTMLParser

BASE = "https://ftp.uniprot.org/pub/databases/uniprot/previous_releases/release-2025_04/knowledgebase/"

class Links(HTMLParser):
    def __init__(self):
        super().__init__()
        self.hrefs = []
    def handle_starttag(self, tag, attrs):
        if tag == "a":
            self.hrefs.extend(v for k, v in attrs if k == "href" and v)

def discover():
    with urllib.request.urlopen(BASE, timeout=120) as response:
        links = Links()
        links.feed(response.read().decode("utf-8", errors="replace"))
    entries = sorted({urllib.parse.urljoin(BASE, h) for h in links.hrefs
                      if urllib.parse.urljoin(BASE, h).startswith(BASE)})
    files = [url for url in entries if re.search(r"\.(?:tar(?:\.gz)?|tgz|fasta(?:\.gz)?)$", url, re.I)]
    direct = [u for u in files if re.search(r"uniprot_(?:sprot|trembl)\.fasta(?:\.gz)?$", u, re.I)]
    if any("sprot" in u for u in direct) and any("trembl" in u for u in direct):
        return direct
    separate = [u for u in files if re.search(r"uniprot_(?:sprot|trembl)[-_]only.*(?:\.tar(?:\.gz)?|\.tgz)$", u, re.I)]
    if any("sprot" in u for u in separate) and any("trembl" in u for u in separate):
        return separate
    combined = [u for u in files if re.search(r"/(?:uniprot|knowledgebase)[-_]?2025_04\.(?:tar(?:\.gz)?|tgz)$", u, re.I)]
    if len(combined) == 1:
        return combined
    raise RuntimeError(f"Could not safely identify complete Swiss-Prot + TrEMBL archives at {BASE}; discovered {files}. Provide explicit archive URLs or local FASTAs.")

def fasta(handle, wanted, output, observed):
    acc, seq = None, []
    def emit():
        if acc in wanted and seq:
            protein = "".join(seq).upper()
            if protein not in observed.setdefault(acc, set()):
                output.write(f">{acc}\n{protein}\n")
                observed[acc].add(protein)
    for line in handle:
        if line.startswith(">"):
            emit()
            tokens = line[1:].split()
            first = tokens[0].split("|") if tokens else []
            acc = first[1] if len(first) >= 3 and first[0] in ("sp", "tr") else (first[0] if first else None)
            seq = []
        elif acc in wanted:
            seq.append(line.strip())
    emit()

def scan_file(stream, name, wanted, output, observed, members, depth=0):
    name = name.lower()
    if depth > 4:
        raise ValueError("Unexpected nested archive depth")
    if re.search(r"\.(?:fasta|faa)(?:\.gz)?$", name):
        # Only protein sequence FASTA, not nucleotide or unrelated exports.
        members.append(name)
        if name.endswith(".gz"):
            with gzip.GzipFile(fileobj=stream) as uncompressed:
                with io.TextIOWrapper(uncompressed, encoding="utf-8") as txt:
                    fasta(txt, wanted, output, observed)
        else:
            txt = io.TextIOWrapper(stream, encoding="utf-8")
            try:
                fasta(txt, wanted, output, observed)
            finally:
                txt.detach()
    elif re.search(r"\.(?:tar(?:\.gz)?|tgz)$", name):
        with tarfile.open(fileobj=stream, mode="r|*") as archive:
            for member in archive:
                if not member.isfile() or not re.search(r"\.(?:fasta|faa|tar)(?:\.gz)?$|\.tgz$", member.name, re.I):
                    continue
                extracted = archive.extractfile(member)
                if extracted is not None:
                    with extracted:
                        scan_file(extracted, member.name, wanted, output, observed, members, depth+1)

def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--metadata", type=Path, required=True)
    p.add_argument("--out", type=Path, required=True)
    p.add_argument("--archive", action="append", default=[])
    args = p.parse_args()
    if args.out.suffix != ".gz" or args.out.exists():
        raise ValueError("Choose a new .fasta.gz output")
    with args.metadata.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        if not {"uniprot_ac_1", "uniprot_ac_2"} <= set(reader.fieldnames or []):
            raise ValueError("Unrecognized AFDB accession schema")
        wanted = set()
        for row in reader:
            wanted.update((row["uniprot_ac_1"].strip(), row["uniprot_ac_2"].strip()))
    wanted.discard("")
    if not wanted:
        raise ValueError("No accessions in source metadata")
    urls = args.archive or discover()
    observed, members = {}, []
    partial = args.out.with_name(args.out.name + ".part")
    with gzip.open(partial, "wt", encoding="utf-8") as out:
        for url in urls:
            # Never silently accept a different release.
            if not url.startswith("https://") or "2025_04" not in url or "/current_release/" in url:
                raise ValueError("Expected explicit HTTPS UniProt 2025_04 archive URL")
            print(f"Reading archive {url}", flush=True)
            before = len(members)
            with urllib.request.urlopen(url, timeout=180) as stream:
                scan_file(stream, urllib.parse.urlparse(url).path, wanted, out, observed, members)
            if len(members) == before:
                raise RuntimeError(f"No FASTA members found in {url}")
    missing = sorted(wanted - observed.keys())
    report = {"release": "2025_04", "archives": urls, "members": members,
              "accessions_requested": len(wanted), "accessions_found": len(observed),
              "accessions_missing": len(missing),
              "ambiguous_accessions": sum(len(v)>1 for v in observed.values())}
    args.out.with_suffix(".source.json").write_text(json.dumps(report, indent=2)+"\n")
    args.out.with_suffix(".missing.txt").write_text("\n".join(missing)+"\n")
    if not observed:
        raise RuntimeError("No metadata accessions found in archive.")
    partial.replace(args.out)
    print(json.dumps(report, indent=2))

if __name__ == "__main__":
    main()
