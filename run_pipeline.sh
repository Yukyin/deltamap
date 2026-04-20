#!/bin/bash
# run_pipeline.sh
# Runs the full deltamap pipeline end to end.
# Execute from the project root directory after completing Data Preparation.
#
# Usage:
#   bash run_pipeline.sh
#
# Expected runtime: 30-90 minutes depending on machine (steps 02 and 03 are
# the slowest because they read all 116 .h5 files from disk).

set -e  # stop on first error

echo "[1/9] Fixing manifest paths..."
python scripts/00_fix_manifest_paths.py

echo "[2/9] Sample-level QC..."
python scripts/02_sample_level_qc.py

echo "[3/9] Building pseudobulk matrix..."
python scripts/03_build_pseudobulk.py

echo "[4/9] Computing paired delta (D2 minus D1)..."
python scripts/04_build_delta.py

echo "[5/9] Scoring evidence modules..."
python scripts/07_score_evidence_modules.py

echo "[6/9] PCA of evidence space..."
python scripts/08_plot_evidence_space_pca.py

echo "[7/9] Neighbor purity audit..."
python scripts/09_neighbor_purity_outlier_audit.py

echo "[8/9] Generating figures..."
python scripts/10_plot_module_effect_bar.py
python scripts/11_plot_neighbor_purity_bar.py
python scripts/12_plot_module_boxplot.py

echo "[9/9] Optional downstream analyses..."
python scripts/run_clustering.py
python scripts/run_bootstrap_stability.py

echo ""
echo "Pipeline complete. Results are in results/"
echo "Key figures:"
echo "  results/evidence_modules_v3_case_minus_control_barplot.png"
echo "  results/evidence_modules_v3_subject_pca.png"
echo "  results/evidence_modules_v3_neighbor_purity_k5.png"
echo "  results/evidence_modules_v3_case_control_boxplot.png"
