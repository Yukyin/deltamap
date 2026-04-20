"""
run_loso_classification.py
Leave-one-subject-out (LOSO) logistic regression on the full evidence feature set.
Reports accuracy, balanced accuracy, and ROC-AUC.

Features used:
  - stable_high / stable_low signature scores
  - signature_axis_score
  - champion module scores (humoral_bcell, innate_apc)
  - champion_axis_score
  - distance-to-prototype features
  - boundary_flag

Note: This validation uses the enriched feature table built by
build_evidencefm_ready_table_v1.py (evidence_fm_ready_subject_table_v1.tsv).
If you only have the basic module matrix, use the simpler 4-feature version below.
"""
import os
import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score, balanced_accuracy_score

os.environ["PANDAS_USE_NUMEXPR"] = "0"

IN_TABLE  = "results/evidence_fm_ready_subject_table_v1.tsv"
OUT_PRED  = "results/evidencefm_ready_loso_logreg_v1_subject_predictions.tsv"
OUT_SUM   = "results/evidencefm_ready_loso_logreg_v1_summary.tsv"
OUT_COEF  = "results/evidencefm_ready_loso_logreg_v1_feature_coefficients.tsv"

FEATURES = [
    "stable_high_signature_score",
    "stable_low_signature_score",
    "signature_axis_score",
    "humoral_bcell_champion",
    "innate_apc_champion",
    "champion_axis_score",
    "distance_to_stable_high_proto",
    "distance_to_stable_low_proto",
    "proto_margin_low_minus_high",
    "boundary_flag",
]


def loso_logreg(df: pd.DataFrame, features: list) -> tuple:
    df     = df.copy()
    df["y"] = (df["label"] == "case").astype(int)
    df["boundary_flag"] = df["boundary_flag"].astype(int)

    pred_rows, coef_rows = [], []
    for held_out in df["subject_id"].tolist():
        train = df[df["subject_id"] != held_out]
        test  = df[df["subject_id"] == held_out]

        X_tr = train[features].to_numpy(dtype=float)
        X_te = test[features].to_numpy(dtype=float)
        y_tr = train["y"].to_numpy(dtype=int)
        y_te = test["y"].to_numpy(dtype=int)

        mu, sd  = X_tr.mean(axis=0), X_tr.std(axis=0, ddof=0)
        sd[sd == 0] = 1.0
        X_tr_z  = (X_tr - mu) / sd
        X_te_z  = (X_te - mu) / sd

        clf = LogisticRegression(
            penalty="l2", C=1.0, solver="liblinear",
            max_iter=5000, class_weight="balanced", random_state=0,
        )
        clf.fit(X_tr_z, y_tr)

        prob_case = float(clf.predict_proba(X_te_z)[0, 1])
        pred_case = int(prob_case >= 0.5)
        pred_rows.append({
            "subject_id": held_out,
            "label":      test["label"].iloc[0],
            "y_true":     int(y_te[0]),
            "prob_case":  prob_case,
            "pred_case":  pred_case,
            "correct":    int(pred_case == int(y_te[0])),
        })
        for feat, coef in zip(features, clf.coef_[0]):
            coef_rows.append({"held_out_subject": held_out, "feature": feat, "coef": float(coef)})

    pred = pd.DataFrame(pred_rows)
    coef = pd.DataFrame(coef_rows)
    return pred, coef


def main():
    df   = pd.read_csv(IN_TABLE, sep="\t")
    keep = df["label"].isin(["case", "control"]).copy()
    for c in FEATURES:
        keep &= df[c].notna()
    df = df.loc[keep].copy()
    print(f"Usable subjects: {len(df)}")

    pred, coef = loso_logreg(df, FEATURES)

    y_true = pred["y_true"].to_numpy(dtype=int)
    y_prob = pred["prob_case"].to_numpy(dtype=float)
    y_pred = pred["pred_case"].to_numpy(dtype=int)

    summary = pd.DataFrame([
        {"metric": "n_subjects",        "value": int(pred.shape[0])},
        {"metric": "n_case",            "value": int((pred["label"] == "case").sum())},
        {"metric": "n_control",         "value": int((pred["label"] == "control").sum())},
        {"metric": "accuracy",          "value": float((y_true == y_pred).mean())},
        {"metric": "balanced_accuracy", "value": float(balanced_accuracy_score(y_true, y_pred))},
        {"metric": "roc_auc",           "value": float(roc_auc_score(y_true, y_prob))},
    ])
    coef_summary = (
        coef.groupby("feature")["coef"]
        .agg(coef_mean="mean", coef_std="std", coef_median="median")
        .reset_index()
    )

    pred.to_csv(OUT_PRED, sep="\t", index=False)
    summary.to_csv(OUT_SUM, sep="\t", index=False)
    coef_summary.to_csv(OUT_COEF, sep="\t", index=False)

    print("\nLOSO summary:")
    print(summary.to_string(index=False))
    print("\nMean feature coefficients:")
    print(coef_summary.sort_values("coef_mean", ascending=False).to_string(index=False))


if __name__ == "__main__":
    main()
