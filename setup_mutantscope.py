#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os, sys, zipfile, pathlib, textwrap

ROOT = pathlib.Path(__file__).resolve().parent
PROJECT = ROOT / "mutantscope_demo"

def write_text(path: pathlib.Path, content: str):
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8", newline="\n") as f:
        f.write(content)

def zipdir(folder: pathlib.Path, zip_path: pathlib.Path):
    with zipfile.ZipFile(zip_path, "w", zipfile.ZIP_DEFLATED) as zf:
        for p in folder.rglob("*"):
            if p.is_file():
                zf.write(p, p.relative_to(folder))

README_MD = textwrap.dedent("""\
    # MutantScope Demo

    Quick start:
    1) python -m venv .venv
    2) .venv\\Scripts\\Activate     (PowerShell: .\\.venv\\Scripts\\Activate)
    3) pip install -r requirements.txt
    4) streamlit run app/app.py
""")

REQUIREMENTS_TXT = "\n".join([
    "streamlit>=1.36",
    "pandas>=2.0",
    "numpy>=1.25",
    "scikit-learn>=1.3",
    "joblib>=1.3",
    "reportlab>=4.0"
]) + "\n"

SAMPLE_CSV = textwrap.dedent("""\
    protein,position,ref_aa,alt_aa,label
    EGFR,790,L,M,Up
    BRCA1,1203,S,F,Down
    TP53,175,R,H,Neutral
""")

AA_PROPS_PY = textwrap.dedent("""\
    _hydro = {
        'A': 1.8,'R': -4.5,'N': -3.5,'D': -3.5,'C': 2.5,'Q': -3.5,'E': -3.5,'G': -0.4,
        'H': -3.2,'I': 4.5,'L': 3.8,'K': -3.9,'M': 1.9,'F': 2.8,'P': -1.6,'S': -0.8,
        'T': -0.7,'W': -0.9,'Y': -1.3,'V': 4.2
    }
    _charge = {'D':-1,'E':-1,'K':1,'R':1,'H':0.1}

    def aa_hydro(a): return _hydro.get(a, 0.0)
    def aa_charge(a): return _charge.get(a, 0.0)

    def add_basic_features(df, colmap):
        df = df.copy()
        r = df[colmap['ref_aa']].str.upper()
        a = df[colmap['alt_aa']].str.upper()
        df['delta_hydro'] = r.map(aa_hydro) - a.map(aa_hydro)
        df['delta_charge'] = r.map(aa_charge) - a.map(aa_charge)
        df['is_conservative'] = (r == a).astype(int)
        return df
""")

BUILD_DATASET_PY = textwrap.dedent("""\
    from aa_props import add_basic_features

    def build_features(df, colmap):
        X = add_basic_features(df, colmap)
        cols = ['delta_hydro','delta_charge','is_conservative']
        if 'label' in colmap and colmap['label'] in df.columns:
            return X[cols + [colmap['label']]]
        return X[cols]
""")

TRAIN_MODEL_PY = textwrap.dedent("""\
    import joblib
    from sklearn.model_selection import train_test_split
    from sklearn.preprocessing import StandardScaler
    from sklearn.pipeline import Pipeline
    from sklearn.linear_model import LogisticRegression

    def train_and_save(Xy, label_col: str, model_path: str):
        X = Xy.drop(columns=[label_col])
        y = Xy[label_col]
        pipe = Pipeline([
            ('scaler', StandardScaler()),
            ('clf', LogisticRegression(max_iter=500, class_weight='balanced', multi_class='auto'))
        ])
        Xtr, Xte, ytr, yte = train_test_split(X, y, test_size=0.25, stratify=y, random_state=42)
        pipe.fit(Xtr, ytr)
        joblib.dump(pipe, model_path)
        return pipe.score(Xte, yte)
""")

SCORE_API_PY = textwrap.dedent("""\
    import joblib, pathlib
    MODEL = pathlib.Path(__file__).resolve().parent / 'models' / 'mutantscope.joblib'

    def score_rows(X):
        model = joblib.load(MODEL)
        return model.predict_proba(X)
""")

REPORT_PY = textwrap.dedent("""\
    from reportlab.lib.pagesizes import A4
    from reportlab.pdfgen import canvas

    def write_report(path, title="MutantScope Report", lines=None):
        c = canvas.Canvas(path, pagesize=A4)
        w, h = A4
        y = h - 72
        c.setFont("Helvetica-Bold", 16)
        c.drawString(72, y, title)
        y -= 24
        c.setFont("Helvetica", 11)
        if lines:
            for line in lines:
                c.drawString(72, y, str(line))
                y -= 14
                if y < 72:
                    c.showPage(); y = h - 72
        c.showPage()
        c.save()
""")

APP_STREAMLIT = textwrap.dedent("""\
    import streamlit as st
    import pandas as pd, joblib, pathlib
    from build_dataset import build_features
    from train_model import train_and_save
    from score_api import score_rows

    ROOT = pathlib.Path(__file__).resolve().parents[1]
    MODELS = ROOT / "models"
    DATA = ROOT / "data"
    MODELS.mkdir(exist_ok=True)
    DATA.mkdir(exist_ok=True)
    MODEL_PATH = MODELS / "mutantscope.joblib"

    st.set_page_config(page_title="MutantScope", layout="wide")
    st.title("🧬 MutantScope — mutation impact demo")

    tab1, tab2 = st.tabs(["Train from file", "Score trained model"])

    with tab1:
        st.subheader("Upload mutations CSV and train")
        up = st.file_uploader("CSV with columns: protein, position, ref_aa, alt_aa, label", type=["csv"])
        if st.button("Use bundled sample"):
            up = open(DATA / "sample_mutations.csv", "rb")
        if up:
            df = pd.read_csv(up)
            colmap = {'protein':'protein','position':'position','ref_aa':'ref_aa','alt_aa':'alt_aa','label':'label'}
            Xy = build_features(df, colmap)
            if 'label' not in Xy.columns:
                st.error("No label column found after featurization.")
            else:
                acc = train_and_save(Xy, 'label', str(MODEL_PATH))
                st.success(f"Model trained and saved to {MODEL_PATH.name}. Holdout accuracy ~ {acc:.3f}")

    with tab2:
        st.subheader("Score with saved model")
        if not MODEL_PATH.exists():
            st.info("Train a model first in the previous tab.")
        else:
            up2 = st.file_uploader("Feature table CSV (delta_hydro, delta_charge, is_conservative)", type=["csv"], key="sc")
            if st.button("Score sample (auto-build)"):
                df = pd.read_csv(DATA / "sample_mutations.csv")
                colmap = {'protein':'protein','position':'position','ref_aa':'ref_aa','alt_aa':'alt_aa','label':'label'}
                Xy = build_features(df, colmap)
                X = Xy.drop(columns=['label'])
                proba = score_rows(X)
                out = X.copy()
                if proba.shape[1] == 2:
                    out['prob_up'] = proba[:,1]
                else:
                    out[['prob_down','prob_neutral','prob_up']] = proba
                out.to_csv(ROOT / "scored_variants.csv", index=False)
                st.success("Wrote scored_variants.csv in project root.")
                st.dataframe(out.head())
            if up2:
                X = pd.read_csv(up2)
                proba = score_rows(X)
                out = X.copy()
                if proba.shape[1] == 2:
                    out['prob_up'] = proba[:,1]
                else:
                    out[['prob_down','prob_neutral','prob_up']] = proba
                st.dataframe(out.head())
""")

def main():
    print(f"[info] Working dir: {ROOT}")
    print("[info] Writing project files...")

    # structure & files
    write_text(PROJECT / "README.md", README_MD)
    write_text(PROJECT / "requirements.txt", REQUIREMENTS_TXT)
    write_text(PROJECT / "data" / "sample_mutations.csv", SAMPLE_CSV)
    write_text(PROJECT / "aa_props.py", AA_PROPS_PY)
    write_text(PROJECT / "build_dataset.py", BUILD_DATASET_PY)
    write_text(PROJECT / "train_model.py", TRAIN_MODEL_PY)
    write_text(PROJECT / "score_api.py", SCORE_API_PY)
    write_text(PROJECT / "report.py", REPORT_PY)
    write_text(PROJECT / "app" / "app.py", APP_STREAMLIT)
    (PROJECT / "models").mkdir(parents=True, exist_ok=True)

    zip_path = ROOT / "mutantscope_demo.zip"
    zipdir(PROJECT, zip_path)

    print(f"[ok] Created folder: {PROJECT}")
    print(f"[ok] Zipped to     : {zip_path}")
    print("\n[next] To run the demo:")
    print("  1) cd mutantscope_demo")
    print("  2) py -m venv .venv")
    print("  3) .\\.venv\\Scripts\\Activate")
    print("  4) pip install -r requirements.txt")
    print("  5) streamlit run app/app.py")

if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        print(f"[error] {e}", file=sys.stderr)
        sys.exit(1)
