"""Small offline tests for AFDB pseudo-label preparation."""
import argparse
import csv
import tempfile
import unittest
from pathlib import Path

from scripts.prepare_afdb_pseudolabels import prepare, pair_key, score_label

class TestAFDB(unittest.TestCase):
    def test_thresholds(self):
        self.assertEqual(score_label(.8, .9, .4, .8), 1)
        self.assertEqual(score_label(.1, .4, .4, .8), 0)
        self.assertIsNone(score_label(.4, .8, .4, .8))
        self.assertIsNone(score_label(.5, .9, .4, .8))
        self.assertEqual(pair_key("AAA", "CCC", 1), pair_key("CCC", "AAA", 1))

    def test_pipeline(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            fasta = root / "proteins.fasta"
            fasta.write_text(">sp|A|one\nAAAA\n>B\nCCCC\n>C\nDDDD\n>D\nEEEE\n"
                             ">E\nFFFF\n>F\nGGGG\n>G\nHHHH\n>aliasA\nAAAA\n")
            metadata = root / "metadata.csv"
            with metadata.open("w", newline="") as handle:
                w = csv.writer(handle)
                w.writerow(["modelEntityId", "uniprot_ac_1", "uniprot_ac_2",
                            "tax_id_1", "tax_id_2", "ipTM"])
                for idx, (a,b,score,taxon) in enumerate([
                    ("A","B",.8,"9606"), ("B","aliasA",.95,"9606"),
                    ("C","D",.4,"9606"), ("D","C",.1,"9606"),
                    ("E","F",.2,"9606"), ("F","E",.9,"9606"),
                    ("A","C",.9,"9606"), ("A","C",.5,"9606"),
                    ("A","G",float("nan"),"9606"), ("A","unknown",.9,"9606"),
                    ("A","aliasA",.9,"9606"), ("G","B",.85,"10090"),
                ]):
                    w.writerow([f"model{idx}",a,b,taxon,taxon,score])
            args = argparse.Namespace(metadata=metadata, fasta=[fasta],
                sequence_release="fixture", out=root/"pilot.csv", taxid="9606",
                positive_min=.8, negative_max=.4, min_length=1, max_length=20,
                max_per_class=100, seed=1)
            result = prepare(args)
            with args.out.open() as handle:
                rows = list(csv.DictReader(handle))
            self.assertEqual(len(rows), 2)
            self.assertEqual([r["interaction_label"] for r in rows], ["1","0"])
            self.assertEqual(rows[0]["metadata_rows_for_pair"], "2")
            self.assertEqual(rows[1]["affinity_nm"], "")
            counts = result["counts"]
            self.assertEqual(counts["discordant_high_low_pairs"], 1)
            self.assertEqual(counts["excluded_uncertain_pairs"], 2)
            self.assertEqual(counts["invalid_iptm"], 1)
            self.assertEqual(counts["identical_chain"], 1)
            self.assertEqual(counts["unresolved_or_filtered_sequence"], 1)
            self.assertEqual(counts["other_taxon"], 1)
            self.assertTrue(args.out.with_suffix(".audit.json").exists())
            with self.assertRaises(FileExistsError):
                prepare(args)

if __name__ == "__main__":
    unittest.main()
