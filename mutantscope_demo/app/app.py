# --- ensure imports work no matter how streamlit runs this file ---
import sys, pathlib
APP_DIR = pathlib.Path(__file__).resolve().parent
if str(APP_DIR) not in sys.path:
    sys.path.insert(0, str(APP_DIR))
# -----------------------------------------------------------------

from vep_client import fetch_vep_annotation
from interpret import summarize_variant

import streamlit as st

st.caption(f"Loaded app.py from: {APP_DIR}")
st.header("Variant Intelligence (Rare Disease Oriented)")

# ---------- small helpers ----------
def safe_str(x): 
    try: return str(x)
    except: return ""

def safe_join(items, sep=", "):
    return sep.join(safe_str(i) for i in (items or []) if i not in (None, "", [], {}))

def note(text: str):
    st.caption(text)

def explain_cadd(v):
    note("**CADD** combines many genomic features; higher = more likely damaging.")
    try: x = float(v)
    except: 
        return "Common rule of thumb: PHRED ≥20 (top ~1%) is considered high."
    if x >= 30: return "Very high (top 0.1% deleterious); concerning."
    if x >= 20: return "High (top ~1% deleterious); concerning."
    return "Below common pathogenic threshold (20); weigh against other evidence."

def explain_revel(v):
    note("**REVEL** is an ensemble for missense variants; 0→1 (higher = worse).")
    try: x = float(v)
    except: 
        return "Cutoffs often used: ≥0.5 concerning, ≥0.75 stronger."
    if x >= 0.85: return "Very strong damaging signal."
    if x >= 0.75: return "Strong; often supportive pathogenic evidence."
    if x >= 0.5:  return "Suggestive of damage."
    return "Below common damaging cutoffs."

def explain_sift(pred, score):
    note("**SIFT** predicts tolerance of the amino-acid change; <0.05 = deleterious.")
    try: s = float(score) if score is not None else None
    except: s = None
    if pred and str(pred).lower().startswith("del"): 
        return "SIFT says not tolerated (deleterious)."
    if s is not None and s < 0.05:
        return "Deleterious by numeric cutoff (<0.05)."
    return "Not strongly deleterious by SIFT."

def explain_polyphen(pred, score):
    note("**PolyPhen-2** predicts protein impact; >0.85 ≈ probably damaging.")
    try: s = float(score) if score is not None else None
    except: s = None
    if pred and "prob" in str(pred).lower(): 
        return "Predicted probably damaging."
    if s is not None and s > 0.85:
        return "Score >0.85: probably damaging."
    return "Not strongly damaging by PolyPhen-2."

def explain_alphamissense(am_class, am_score):
    note("**AlphaMissense** (AI) classifies missense pathogenicity; higher = worse.")
    try: s = float(am_score) if am_score is not None else None
    except: s = None
    if am_class and str(am_class).lower() in ("likely_pathogenic", "pathogenic"):
        return "Class indicates likely pathogenicity; supportive."
    if s is not None:
        if s >= 0.8: return "Score ≥0.8: concerning."
        if s >= 0.6: return "Moderate concern; combine with other lines."
    return "No strong pathogenic signal alone."

def explain_gnomad_af(af, popmax):
    note("**gnomAD AF** shows rarity in large populations; pathogenic variants are usually extremely rare.")
    def fmt(x):
        try: return f"{float(x):.2e}"
        except: return safe_str(x)
    if af is None and popmax is None:
        return "No population frequency returned."
    tip = f"AF={fmt(af)}; PopMax={fmt(popmax)}. "
    try:
        v = float(popmax if popmax is not None else af)
        if v <= 1e-6: return tip + "Ultra-rare; consistent with severe rare disease."
        if v <= 1e-5: return tip + "Extremely rare; consistent with rare disease."
        if v <= 1e-4: return tip + "Rare; could be compatible depending on phenotype."
        return tip + "Too common for most severe Mendelian disorders."
    except:
        return tip + "Rarity supports pathogenicity if phenotype matches."

def explain_conservation(gerp, phastcons):
    note("**Conservation**: important sites evolve slowly; high scores imply functional constraint.")
    parts = []
    try:
        if gerp is not None:
            g = float(gerp)
            if g > 4: parts.append("GERP++ >4: very strong constraint.")
            elif g > 2: parts.append("GERP++ >2: conserved.")
            else: parts.append("GERP++ low.")
    except: pass
    try:
        if phastcons is not None:
            p = float(phastcons)
            if p >= 0.99: parts.append("PhastCons ≈1: highly conserved across vertebrates.")
            elif p >= 0.8: parts.append("PhastCons high.")
            else: parts.append("PhastCons low.")
    except: pass
    return " ".join(parts) or "No conservation signal provided."

# ---------- sidebar ----------
st.sidebar.markdown("### VEP / HGVS Lookup")
hgvs = st.sidebar.text_input(
    "HGVS notation",
    value="NC_000023.11:g.153866826C>T",
    help="Paste an HGVS genomic notation (GRCh38). Example: NC_000023.11:g.153866826C>T",
)
species = st.sidebar.selectbox("Species", ["human"], index=0)
fetch_btn = st.sidebar.button("Fetch from VEP")

# ---------- main ----------
if fetch_btn:
    with st.spinner("Contacting VEP and interpreting..."):
        try:
            data = fetch_vep_annotation(hgvs, species)
            if not isinstance(data, list):
                st.error("Unexpected VEP response format.")
            else:
                summary = summarize_variant(data)
                if "error" in summary:
                    st.error(summary["error"])
                else:
                    st.subheader("Overview")
                    c1, c2, c3 = st.columns(3)
                    c1.metric("Gene", summary.get("gene") or "—")
                    c2.metric("Consequence", summary.get("consequence") or "—")
                    c3.metric("Impact (VEP)", summary.get("impact") or "—")
                    note("VEP consequence/impact are heuristic severities; combine with clinical & population evidence.")
                    st.write(f"**Variant (HGVSg):** {summary.get('hgvsg') or '—'}")
                    st.write(f"**Protein change:** {summary.get('hgvsp') or '—'}")

                    tabs = st.tabs([
                        "Clinical",
                        "Computational",
                        "Population & Conservation",
                        "Domains & GO",
                        "Literature",
                        "Actionable next steps",
                    ])

                    # Clinical
                    with tabs[0]:
                        st.markdown("#### Clinical evidence")
                        sig = safe_join(summary.get("clinvar_significance"))
                        rev = safe_join(summary.get("clinvar_review"))
                        st.write(f"- **ClinVar significance:** {sig or '—'}")
                        st.write(f"- **ClinVar review status:** {rev or '—'}")
                        phs = summary.get("phenotypes") or []
                        if phs:
                            st.write("**Associated phenotypes/diseases:**")
                            for p in phs[:15]:
                                st.write(f"- {safe_str(p)}")
                            if len(phs) > 15:
                                st.write(f"... and {len(phs)-15} more")
                        else:
                            st.info("No phenotype annotations returned in this response.")

                    # Computational
                    with tabs[1]:
                        st.markdown("#### In-silico / AI predictors")
                        am = summary.get("alphamissense") or {}
                        st.write(f"- **AlphaMissense:** {safe_str(am.get('class') or '—')} (score={safe_str(am.get('score'))})")
                        st.caption(explain_alphamissense(am.get('class'), am.get('score')))

                        sift = summary.get("sift", {}) or {}
                        st.write(f"- **SIFT:** {safe_str(sift.get('pred')) or '—'} (score={safe_str(sift.get('score'))})")
                        st.caption(explain_sift(sift.get('pred'), sift.get('score')))

                        poly = summary.get("polyphen", {}) or {}
                        st.write(f"- **PolyPhen-2:** {safe_str(poly.get('pred')) or '—'} (score={safe_str(poly.get('score'))})")
                        st.caption(explain_polyphen(poly.get('pred'), poly.get('score')))

                        cadd = summary.get("cadd_phred")
                        st.write(f"- **CADD (PHRED):** {safe_str(cadd)}")
                        st.caption(explain_cadd(cadd))

                        revel = summary.get("revel")
                        st.write(f"- **REVEL:** {safe_str(revel)}")
                        st.caption(explain_revel(revel))

                    # Population & Conservation
                    with tabs[2]:
                        st.markdown("#### Population & Conservation")
                        af = summary.get("gnomad_af"); pmax = summary.get("gnomad_popmax_af")
                        st.write(f"- **gnomAD AF (v4.1 joint):** {safe_str(af)}")
                        st.write(f"- **gnomAD PopMax AF:** {safe_str(pmax)}")
                        st.caption(explain_gnomad_af(af, pmax))

                        cons = summary.get("conservation") or {}
                        gerp = cons.get("gerp"); phc = cons.get("phastcons100way")
                        st.write(f"- **GERP++:** {safe_str(gerp)}")
                        st.write(f"- **PhastCons 100-way:** {safe_str(phc)}")
                        st.caption(explain_conservation(gerp, phc))

                    # Domains & GO
                    with tabs[3]:
                        st.markdown("#### Protein domains & GO")
                        doms = summary.get("domains") or []
                        if doms:
                            st.write("**Domains (from InterPro/Pfam/structure mappings):**")
                            st.write(safe_join(doms))
                        else:
                            st.info("No domain annotations returned.")
                        gos = summary.get("go_terms") or []
                        if gos:
                            st.write("**Gene Ontology terms (selected):**")
                            for g in gos[:20]:
                                st.write(f"- {safe_str(g)}")

                    # Literature
                    with tabs[4]:
                        pmids = summary.get("pubmed_ids") or []
                        if pmids:
                            st.markdown("#### Publications")
                            for pid in pmids[:20]:
                                st.write(f"- PubMed: https://pubmed.ncbi.nlm.nih.gov/{safe_str(pid)}/")
                            if len(pmids) > 20:
                                st.write(f"... and {len(pmids)-20} more")
                        else:
                            st.info("No PubMed identifiers returned in this payload.")

                        ints = summary.get("interactors") or []
                        if ints:
                            st.markdown("#### Known/Reported interactors (IntAct)")
                            st.write(safe_join(ints))

                    # Next steps
                    with tabs[5]:
                        st.markdown("#### Suggested next steps (non-medical)")
                        st.write(
                            "- Validate the variant context (GRCh38, HGVS form).\n"
                            "- Aggregate evidence: ClinVar, AI predictors (AlphaMissense/REVEL/CADD), conservation, and **extreme rarity** in gnomAD.\n"
                            "- Use PubMed links above; prioritize phenotype-matched reports.\n"
                            "- For classification or testing decisions, consult clinical genetics.\n"
                            "- Functional work: prioritize assays tied to the listed **domains** and **GO** processes."
                        )
                        st.caption("This tool provides research information only and is **not medical advice**.")
        except Exception as e:
            st.error(f"VEP fetch failed: {e}")
