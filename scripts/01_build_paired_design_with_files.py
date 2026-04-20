"""
01_build_paired_design_with_files.py
Join the paired design table (which lists D1/D2 GSM pairs per subject)
with the fixed manifest to attach raw file paths for each time point.

Input:
  meta/sample_manifest_paired_with_files.fixed.tsv  — GSM → raw_file mapping
  meta/paired_design.tsv                            — subject-level D1/D2 GSM pairs

Output:
  meta/paired_design_with_files.fixed.tsv           — paired design with raw file paths
"""
import os
import csv

os.environ["PANDAS_USE_NUMEXPR"] = "0"

MANIFEST   = "meta/sample_manifest_paired_with_files.fixed.tsv"
PAIRED     = "meta/paired_design.tsv"
OUT        = "meta/paired_design_with_files.fixed.tsv"


def main():
    # Build GSM → raw_file lookup from manifest
    lookup = {}
    with open(MANIFEST, "r", encoding="utf-8") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            lookup[r["gsm"]] = r["raw_file"]

    # Attach raw file paths to paired design
    rows = []
    with open(PAIRED, "r", encoding="utf-8") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            r["D1_raw_file"] = lookup.get(r["D1_gsm"], "")
            r["D2_raw_file"] = lookup.get(r["D2_gsm"], "")
            rows.append(r)

    with open(OUT, "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(rows)

    print(f"Wrote {len(rows)} rows to {OUT}")
    missing_d1 = sum(1 for r in rows if not r["D1_raw_file"])
    missing_d2 = sum(1 for r in rows if not r["D2_raw_file"])
    if missing_d1 or missing_d2:
        print(f"Warning: {missing_d1} missing D1 paths, {missing_d2} missing D2 paths")


if __name__ == "__main__":
    main()
