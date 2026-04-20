"""
00_fix_manifest_paths.py
Prepend 'raw/' to raw_file paths in the paired manifest if not already present.
"""
import os
import csv

os.environ["PANDAS_USE_NUMEXPR"] = "0"

INP = "meta/sample_manifest_paired_with_files.tsv"
OUT = "meta/sample_manifest_paired_with_files.fixed.tsv"


def main():
    rows = []
    with open(INP, "r", encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        fieldnames = reader.fieldnames[:]
        for r in reader:
            p = r["raw_file"]
            if p and not p.startswith("raw/"):
                r["raw_file"] = os.path.join("raw", p)
            rows.append(r)

    with open(OUT, "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        w.writeheader()
        w.writerows(rows)

    print(f"Wrote {len(rows)} rows to {OUT}")


if __name__ == "__main__":
    main()
