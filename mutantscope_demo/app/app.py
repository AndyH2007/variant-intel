# app/app.py
import sys, pathlib, json, re, datetime as dt
from typing import List, Dict, Any, Optional

APP_DIR = pathlib.Path(__file__).resolve().parent
if str(APP_DIR) not in sys.path:
    sys.path.insert(0, str(APP_DIR))

from vep_client import fetch_vep_annotation
from interpret import summarize_variant

import streamlit as st

# ---------------- Page & global style ----------------
st.set_page_config(page_title="Variant Intelligence", layout="wide", initial_sidebar_state="expanded")
st.markdown("""
<style>
/* Layout */
section[data-testid="stSidebar"] { width: 320px !important; }
div.block-container { padding-top: .75rem; padding-bottom: 2rem; }

/* Cards */
.card { border-radius: 16px; padding: 16px 18px; background: #11151b; border: 1px solid #1f2630; }
.card + .card { margin-top: 10px; }

/* Hero cards */
.hero-grid { display:grid; grid-template-columns: repeat(4, minmax(0, 1fr)); gap:14px; margin: 6px 0 14px 0; }
.hero { border-radius:16px; border:1px solid #222a33; background:linear-gradient(180deg,#0f141a 0%, #0c1116 100%); padding:14px 16px; }
.hero .k { opacity:.85; font-size:.9rem; margin-bottom:6px; }
.hero .v { font-weight:700; font-size:1.05rem; }

/* Badge */
.badge { display:inline-flex; gap:8px; align-items:center; padding:6px 10px; border-radius:999px; border:1px solid var(--bd,#2b343f); background: var(--bg,#141a21); font-weight:600; }
.badge.dot::before { content:""; width:8px; height:8px; border-radius:999px; background: var(--dot,#64748b); display:inline-block; }

/* Gauges */
.gbar { height: 10px; background:#212834; border-radius:6px; position:relative; overflow:hidden; }
.gbar > div { height:10px; background:#3b82f6; width:0%; }

/* Chip grid */
.tag { display:inline-block; margin:6px 8px 0 0; padding:6px 12px;
       border-radius:999px; border:1px solid #2b2f36; background:#14171c; font-size:.92rem;
       max-width:100%; overflow:hidden; text-overflow:ellipsis; white-space:nowrap; }

/* Row layout */
.grid2 { display:grid; grid-template-columns: 1fr 1fr; gap:18px; }
hr { border:none; height:1px; background:#1d222a; margin:10px 0 14px 0; }

/* evidence table tweaks */
thead tr th { position: sticky; top: 0; background: #0e1318; }
.help { opacity:.75; font-size:.85rem; }
</style>
""", unsafe_allow_html=True)

# ---------------- Session state scaffolding ----------------
if "bookmarks" not in st.session_state:
    st.session_state.bookmarks: List[Dict[str, Any]] = []
if "last_summary" not in st.session_state:
    st.session_state.last_summary: Optional[Dict[str, Any]] = None
if "priors" not in st.session_state:
    st.session_state.priors = {
        "disease_context": "Rare, severe childhood-onset",
        "af_cutoff": 1e-5,
        "cadd_cutoff": 20.0,
        "revel_cutoff": 0.5,
        "revel_strong": 0.75,
        "sift_cutoff": 0.05
    }

# ---------------- Sidebar: Nav + Priors ----------------
st.sidebar.markdown("## Navigation")
mode = st.sidebar.radio("Mode", ["Single", "Batch", "Bookmarks", "Compare"], index=0, label_visibility="collapsed")

st.sidebar.markdown("## VEP / HGVS Lookup")
hgvs = st.sidebar.text_input("HGVS notation (GRCh38)",
    value="NC_000023.11:g.153866826C>T", key="input_hgvs",
    help="Example: NC_000023.11:g.153866826C>T")
species = st.sidebar.selectbox("Species", ["human"], index=0, key="input_species")
fetch_btn = st.sidebar.button("Fetch from VEP", use_container_width=True)

st.sidebar.markdown("---")
st.sidebar.markdown("## Profile / Priors")
st.sidebar.caption("Used to tune explanations (does not change raw data).")
pri = st.session_state.priors
pri["disease_context"] = st.sidebar.text_input("Clinical context", pri["disease_context"])
pri["af_cutoff"] = float(st.sidebar.text_input("Max AF consistent with phenotype", str(pri["af_cutoff"])))
c1, c2, c3 = st.sidebar.columns(3)
with c1: pri["cadd_cutoff"]  = float(st.text_input("CADD",  str(pri["cadd_cutoff"])))
with c2: pri["revel_cutoff"] = float(st.text_input("REVEL", str(pri["revel_cutoff"])))
with c3: pri["revel_strong"] = float(st.text_input("REVEL+", str(pri["revel_strong"])))
pri["sift_cutoff"] = float(st.sidebar.text_input("SIFT (< is deleterious)", str(pri["sift_cutoff"])))
st.sidebar.markdown("---")
st.sidebar.caption("Tip: Bookmark variants in *Single* mode and view them on the *Bookmarks* page.")

# ---------------- Helpers ----------------
def safe_float(x) -> Optional[float]:
    try: return float(x)
    except: return None

def severity_badge(score: float) -> str:
    """
    Map a 0-5 score to a qualitative badge.
    """
    score = max(0.0, min(5.0, score))
    if score >= 4.0:
        bg, dot, bd, txt = "#10231c", "#22c55e", "#1f422e", "likely pathogenic"
    elif score >= 3.0:
        bg, dot, bd, txt = "#2b1d00", "#f59e0b", "#433214", "uncertain — leaning pathogenic"
    elif score >= 2.0:
        bg, dot, bd, txt = "#1e1e27", "#a5b4fc", "#303046", "uncertain significance"
    else:
        bg, dot, bd, txt = "#22191a", "#ef4444", "#402728", "likely benign / conflicting"
    return f"""<span class="badge dot" style="--bg:{bg};--dot:{dot};--bd:{bd}">Pathogenicity: {txt} (index {score:.1f}/5)</span>"""

def interpret_ranges(value, cut_low, cut_high, good_low=True):
    """Return a short interpretation string vs cutoffs."""
    try:
        v = float(value)
    except:
        return ""
    if good_low:
        if v <= cut_low: return "Ultra-rare; supportive."
        if v <= cut_high: return "Rare; compatible."
        return "Common for severe Mendelian disease."
    else:
        if v >= cut_high: return "High; concerning."
        if v >= cut_low:  return "Moderate; supportive."
        return "Low; less concerning."

def explain_af(af, popmax):
    cutoff = pri["af_cutoff"]
    def fmt(x):
        try: return f"{float(x):.2e}"
        except: return str(x)
    txt = f"AF={fmt(af)}; PopMax={fmt(popmax)}. "
    try:
        v = float(popmax if popmax is not None else af)
    except:
        return txt + "Rarity can support pathogenicity if clinical context fits."
    if v <= cutoff/10: return txt + "Ultra-rare vs prior → supportive."
    if v <= cutoff:    return txt + "Rare enough for the current phenotype."
    if v <= cutoff*10: return txt + "Somewhat common; consider inheritance/penetrance."
    return txt + "Too common for most severe Mendelian disorders."

def chip_grid(items: List[str], cols=3, key_prefix="chips", placeholder="Filter…"):
    items = [x for x in (items or []) if x]
    flt = st.text_input(placeholder, key=f"{key_prefix}_filter")
    if flt:
        q = flt.lower()
        items = sorted([x for x in items if q in x.lower()], key=lambda s: s.lower().find(q))
    cols_list = st.columns(cols)
    for i, s in enumerate(items):
        with cols_list[i % cols]:
            st.markdown(f"<span class='tag' title='{s}'>{s}</span>", unsafe_allow_html=True)

def gbar(value: Optional[float], min_val: float, max_val: float, suffix: str = "", key: str = ""):
    if value is None:
        st.caption("—")
        return
    try:
        ratio = (float(value) - min_val) / (max_val - min_val)
        ratio = max(0.0, min(1.0, ratio))
    except:
        ratio = 0.0
    st.markdown(f"""
    <div class="gbar"><div style="width:{ratio*100:.1f}%"></div></div>
    <div class="help">{value}{suffix}</div>
    """, unsafe_allow_html=True)

def compute_pathogenicity_index(s: Dict[str, Any]) -> (float, List[Dict[str, Any]]):
    """
    Simple transparent heuristic (0-5) used ONLY to render the badge.
    Weights are deliberately small and additive; this is NOT ACMG classification.
    """
    w = []  # rows for the evidence table (Source, Value, Interpretation, Weight)
    score = 0.0

    # ClinVar
    clin = ", ".join(s.get("clinvar_significance") or []) or ""
    if "pathogenic" in clin:
        score += 3.0; w.append({"Source":"ClinVar","Value":clin,"Interpretation":"Reported pathogenic","Weight":"+3.0"})
    elif "likely_pathogenic" in clin:
        score += 2.5; w.append({"Source":"ClinVar","Value":clin,"Interpretation":"Likely pathogenic","Weight":"+2.5"})
    elif clin:
        w.append({"Source":"ClinVar","Value":clin,"Interpretation":"Other/conflicting","Weight":"+0.0"})

    # AlphaMissense
    am = s.get("alphamissense") or {}
    amc = (am.get("class") or "").lower()
    ams = safe_float(am.get("score"))
    if amc in ("pathogenic","likely_pathogenic"):
        score += 1.0; w.append({"Source":"AlphaMissense","Value":f"{am.get('class')} ({ams})","Interpretation":"AI pathogenic class","Weight":"+1.0"})
    elif ams and ams >= 0.8:
        score += 0.6; w.append({"Source":"AlphaMissense","Value":f"{ams}","Interpretation":"High score (≥0.8)","Weight":"+0.6"})
    elif ams:
        w.append({"Source":"AlphaMissense","Value":f"{ams}","Interpretation":"Moderate/low","Weight":"+0.0"})

    # REVEL
    revel = safe_float(s.get("revel"))
    if revel is not None:
        if revel >= pri["revel_strong"]: score += 0.6; w.append({"Source":"REVEL","Value":revel,"Interpretation":"Strong (≥0.75)","Weight":"+0.6"})
        elif revel >= pri["revel_cutoff"]: score += 0.3; w.append({"Source":"REVEL","Value":revel,"Interpretation":"Suggestive (≥0.5)","Weight":"+0.3"})
        else: w.append({"Source":"REVEL","Value":revel,"Interpretation":"Below damaging cutoffs","Weight":"+0.0"})

    # CADD
    cadd = safe_float(s.get("cadd_phred"))
    if cadd is not None:
        if cadd >= 30: score += 0.6; w.append({"Source":"CADD","Value":cadd,"Interpretation":"Very high (≥30)","Weight":"+0.6"})
        elif cadd >= pri["cadd_cutoff"]: score += 0.3; w.append({"Source":"CADD","Value":cadd,"Interpretation":f"High (≥{pri['cadd_cutoff']})","Weight":"+0.3"})
        else: w.append({"Source":"CADD","Value":cadd,"Interpretation":"Low/Moderate","Weight":"+0.0"})

    # SIFT / PolyPhen
    sp = ((s.get("sift") or {}).get("pred") or "").lower()
    if "del" in sp: score += 0.2; w.append({"Source":"SIFT","Value":sp,"Interpretation":"Deleterious","Weight":"+0.2"})
    pp = ((s.get("polyphen") or {}).get("pred") or "").lower()
    if "prob" in pp: score += 0.2; w.append({"Source":"PolyPhen-2","Value":pp,"Interpretation":"Probably damaging","Weight":"+0.2"})

    # Population
    pmax = safe_float(s.get("gnomad_popmax_af"))
    if pmax is not None:
        if pmax <= pri["af_cutoff"]: score += 1.0; w.append({"Source":"Population","Value":pmax,"Interpretation":"Rare enough for phenotype","Weight":"+1.0"})
        elif pmax <= pri["af_cutoff"]*10: score += 0.3; w.append({"Source":"Population","Value":pmax,"Interpretation":"Somewhat rare","Weight":"+0.3"})
        else: w.append({"Source":"Population","Value":pmax,"Interpretation":"Too common for severe phenotype","Weight":"0.0"})

    # Conservation
    gerp = safe_float((s.get("conservation") or {}).get("gerp"))
    if gerp and gerp > 4: score += 0.3; w.append({"Source":"Conservation","Value":gerp,"Interpretation":"GERP++ >4","Weight":"+0.3"})
    elif gerp and gerp > 2: score += 0.15; w.append({"Source":"Conservation","Value":gerp,"Interpretation":"GERP++ >2","Weight":"+0.15"})
    else: w.append({"Source":"Conservation","Value":gerp,"Interpretation":"Low/NA","Weight":"+0.0"})

    return min(5.0, score), w

def painlessly_download_button(label, content, filename):
    st.download_button(label, data=content, file_name=filename, use_container_width=True)

def bookmark_current(summary: Dict[str, Any], note: str = ""):
    if not summary: return
    hgvsg = summary.get("hgvsg")
    if not hgvsg: return
    st.session_state.bookmarks = [b for b in st.session_state.bookmarks if b.get("hgvsg") != hgvsg]
    st.session_state.bookmarks.append({
        "hgvsg": hgvsg,
        "gene": summary.get("gene"),
        "consequence": summary.get("consequence"),
        "impact": summary.get("impact"),
        "revel": summary.get("revel"),
        "cadd": summary.get("cadd_phred"),
        "am_class": (summary.get("alphamissense") or {}).get("class"),
        "am_score": (summary.get("alphamissense") or {}).get("score"),
        "timestamp": dt.datetime.utcnow().isoformat(timespec="seconds") + "Z",
        "note": note.strip()
    })

# ---------------- Single mode ----------------
def render_single():
    st.title("Variant Intelligence — Single")
    st.caption("Research-only; not medical advice.")

    # Fetch
    if fetch_btn:
        with st.spinner("Fetching from Ensembl VEP…"):
            try:
                raw = fetch_vep_annotation(hgvs, species)
                if not isinstance(raw, list):
                    st.error("Unexpected VEP response."); return
                s = summarize_variant(raw)
                if "error" in s: st.error(s["error"]); return
                st.session_state.last_summary = s
                st.session_state.last_raw = raw
            except Exception as e:
                st.error(f"VEP fetch failed: {e}"); return

    s = st.session_state.last_summary
    if not s:
        st.info("Enter an HGVS variant on the left and click **Fetch from VEP**.")
        return

    # ---- HERO OVERVIEW ----
    score, weights = compute_pathogenicity_index(s)
    st.markdown('<div class="hero-grid">', unsafe_allow_html=True)
    st.markdown(f'''
      <div class="hero"><div class="k">Gene</div><div class="v">{s.get("gene") or "—"}</div></div>
      <div class="hero"><div class="k">Variant (HGVSg)</div><div class="v">{s.get("hgvsg") or "—"}</div></div>
      <div class="hero"><div class="k">Consequence</div><div class="v">{s.get("consequence") or "—"}</div></div>
      <div class="hero"><div class="k">Impact (VEP)</div><div class="v">{s.get("impact") or "—"}</div></div>
    ''', unsafe_allow_html=True)
    st.markdown('</div>', unsafe_allow_html=True)

    st.markdown(severity_badge(score), unsafe_allow_html=True)
    st.caption(f"Context: {pri['disease_context']}")

    # Actions row
    a1, a2, a3, a4 = st.columns([1,1,1,2])
    with a1:
        note = st.text_input("Optional note for bookmark", "", key="bm_note")
        if st.button("🔖 Bookmark this variant", use_container_width=True, key="bm_btn"):
            bookmark_current(s, note)
            st.success("Bookmarked.")
    with a2:
        painlessly_download_button("⬇️ VEP JSON", json.dumps(st.session_state.last_raw, indent=2), f"{s.get('hgvsg','variant')}.json")
    with a3:
        md = f"""# Variant Summary

- Gene: **{s.get('gene','—')}**
- Variant: **{s.get('hgvsg','—')}**
- Protein: **{s.get('hgvsp','—')}**
- Consequence: **{s.get('consequence','—')}**, Impact: **{s.get('impact','—')}**
- REVEL: **{s.get('revel')}**, CADD: **{s.get('cadd_phred')}**, AlphaMissense: **{(s.get('alphamissense') or {}).get('class')}** ({(s.get('alphamissense') or {}).get('score')})
- AF: **{s.get('gnomad_af')}**, PopMax: **{s.get('gnomad_popmax_af')}**

> Heuristic pathogenicity index: {score:.1f}/5 (NOT ACMG).
"""
        painlessly_download_button("📝 Export Markdown", md, f"{s.get('hgvsg','variant')}.md")
    with a4:
        st.caption("Protein change")
        st.code(s.get("hgvsp") or "—", language="text")

    st.markdown("")

    # ---- CLINICAL / IN-SILICO (two columns) ----
    L, R = st.columns([1,1])

    with L:
        st.subheader("Clinical")
        st.write("**ClinVar significance:**", ", ".join(s.get("clinvar_significance") or []) or "—")
        st.caption("Curation by clinical labs; strongest single source when available.")
        st.write("**ClinVar review:**", ", ".join(s.get("clinvar_review") or []) or "—")

        # Phenotype matcher
        phs = s.get("phenotypes") or []
        with st.expander("Associated phenotypes (filterable & match)", expanded=False):
            query = st.text_input("Enter a patient phenotype phrase (optional)", key="phenomatch")
            if query:
                toks = {t for t in re.split(r"[^a-zA-Z0-9]+", query.lower()) if t}
                def score_row(p): 
                    pp = set(re.split(r"[^a-zA-Z0-9]+", p.lower()))
                    inter = len(toks & pp); return -inter  # sort desc
                phs_sorted = sorted(phs, key=score_row)
                chip_grid(phs_sorted, cols=2, key_prefix=f"ph_{s.get('hgvsg','curr')}")
            else:
                chip_grid(phs, cols=2, key_prefix=f"ph_{s.get('hgvsg','curr')}")
        st.markdown("")

        st.subheader("Population & Conservation")
        st.write("**gnomAD AF (v4.1):**")
        st.caption(explain_af(s.get("gnomad_af"), s.get("gnomad_popmax_af")))
        gbar(safe_float(s.get("gnomad_af")), 0.0, max(pri["af_cutoff"]*20, 1e-3), "", key="afbar")
        st.write("**PopMax AF:**", s.get("gnomad_popmax_af"))
        st.write("**GERP++:**", (s.get("conservation") or {}).get("gerp"),)
        st.caption("Higher means stronger evolutionary constraint.")
        st.write("**PhastCons 100-way:**", (s.get("conservation") or {}).get("phastcons100way"))
        st.caption("Probability of conservation across vertebrates.")

    with R:
        st.subheader("In-silico / AI")
        am = s.get("alphamissense") or {}
        st.write("**AlphaMissense:**", f"{am.get('class') or '—'} (score={am.get('score')})")
        st.caption("AI classifier trained on missense pathogenicity.")
        st.write("**SIFT:**", f"{(s.get('sift') or {}).get('pred','—')} (score={(s.get('sift') or {}).get('score')})")
        st.caption(f"Deleterious < {pri['sift_cutoff']}.")
        st.write("**PolyPhen-2:**", f"{(s.get('polyphen') or {}).get('pred','—')} (score={(s.get('polyphen') or {}).get('score')})")
        st.caption("Probably damaging typically > 0.85.")
        st.write("**CADD (PHRED):**")
        gbar(safe_float(s.get("cadd_phred")), 0.0, 35.0, "", key="caddbar")
        st.caption(f"≥{pri['cadd_cutoff']} high; ≥30 very high.")
        st.write("**REVEL:**")
        gbar(safe_float(s.get("revel")), 0.0, 1.0, "", key="revelbar")
        st.caption(f"≥{pri['revel_cutoff']} suggestive; ≥{pri['revel_strong']} stronger.")

    # ---- Evidence table (structured overview) ----
    st.subheader("Structured evidence overview")
    import pandas as pd
    df = pd.DataFrame(weights, columns=["Source","Value","Interpretation","Weight"])
    st.dataframe(df, use_container_width=True, hide_index=True)
    st.caption("Heuristic evidence and weights used for the badge above (NOT ACMG classification).")

    # ---- Domains & GO ----
    st.subheader("Domains & GO")
    doms = s.get("domains") or []
    if doms:
        with st.expander("Domains (click to expand)", expanded=False):
            chip_grid(doms, cols=4, key_prefix=f"dom_{s.get('hgvsg','curr')}")
    else:
        st.info("No domain annotations in this response.")
    gos = s.get("go_terms") or []
    if gos:
        with st.expander("GO terms (click to expand)", expanded=False):
            chip_grid(gos, cols=3, key_prefix=f"go_{s.get('hgvsg','curr')}")
    else:
        st.info("No GO terms in this response.")

    # ---- Protein ruler ----
    st.subheader("Protein position")
    with st.container(border=True):
        if s.get("protein_pos") is None:
            st.caption("Position unavailable in this payload.")
        else:
            pos = int(s["protein_pos"]); approx_len = s.get("protein_len") or max(1000, pos+50)
            frac = min(max(pos/approx_len, 0.0), 1.0)
            st.markdown(
                f"""
                <div style="height:12px;background:#1f2630;border-radius:6px;position:relative;margin:8px 6px 2px 6px;">
                  <div style="height:12px;width:{frac*100:.2f}%;background:#3b82f6;border-radius:6px;"></div>
                  <div style="position:absolute;top:-8px;left:calc({frac*100:.2f}% - 6px);width:0;height:0;border-left:6px solid transparent;border-right:6px solid transparent;border-top:8px solid #3b82f6;"></div>
                </div>
                <div class="help" style="margin-left:6px;">aa {pos} of ~{approx_len}</div>
                """, unsafe_allow_html=True
            )

    # ---- Literature ----
    st.subheader("Literature")
    pmids = s.get("pubmed_ids") or []
    if pmids:
        for pid in pmids[:20]:
            st.write(f"- PubMed: https://pubmed.ncbi.nlm.nih.gov/{pid}/")
        if len(pmids) > 20: st.caption(f"...and {len(pmids)-20} more")
    else:
        st.info("No PubMed IDs found in this payload.")

    # ---- Next steps ----
    st.subheader("Suggested next steps")
    st.write(
        "- Validate phenotype fit and inheritance.\n"
        "- Combine ClinVar, AI predictors, population rarity and conservation.\n"
        "- Use PubMed links for phenotype-matched reports.\n"
        "- For functional follow-up, target domains/GO terms.\n"
        "- For clinical questions, consult a medical genetics professional (ACMG/AMP)."
    )
    st.caption("This tool provides research information only and is **not** medical advice.")

# ---------------- Batch mode ----------------
def render_batch():
    st.title("Variant Intelligence — Batch")
    text = st.text_area("Paste multiple HGVS lines (one per line)", height=180, placeholder="NC_000023.11:g.153866826C>T\n...")
    run = st.button("Run batch", type="primary")
    if not run:
        st.info("Paste variants and click **Run batch**."); return
    lines = [l.strip() for l in text.splitlines() if l.strip()]
    if not lines:
        st.warning("No inputs given."); return
    rows = []
    with st.spinner("Fetching from VEP…"):
        for hg in lines:
            try:
                raw = fetch_vep_annotation(hg, "human")
                s = summarize_variant(raw) if isinstance(raw, list) else {"error":"bad response"}
            except Exception as e:
                s = {"error": str(e)}
            rows.append({
                "HGVS": hg, "Gene": s.get("gene"), "Consq": s.get("consequence"),
                "Impact": s.get("impact"), "REVEL": s.get("revel"),
                "CADD": s.get("cadd_phred"), "AF": s.get("gnomad_af"),
                "_summary": s
            })
    import pandas as pd
    df = pd.DataFrame([{k:v for k,v in r.items() if k != "_summary"} for r in rows])
    st.dataframe(df, use_container_width=True, hide_index=True)
    options = [f"{r['HGVS']} ({r.get('Gene') or '—'})" for r in rows]
    pick = st.multiselect("Bookmark selected", options)
    if st.button("Add to bookmarks"):
        for label, r in zip(options, rows):
            if label in pick and isinstance(r["_summary"], dict) and "error" not in r["_summary"]:
                bookmark_current(r["_summary"])
        st.success("Added selected to bookmarks.")

# ---------------- Bookmarks mode ----------------
def render_bookmarks():
    st.title("Bookmarks")
    if not st.session_state.bookmarks:
        st.info("No bookmarks yet. Use **Single** or **Batch** to add."); return
    q = st.text_input("Filter bookmarks (gene / HGVS / consequence)", "")
    data = st.session_state.bookmarks
    if q:
        ql = q.lower()
        data = [b for b in data if ql in (b.get("hgvsg","").lower()+b.get("gene","").lower()+str(b.get("consequence","")).lower())]
    for i, b in enumerate(data):
        with st.container(border=True):
            c1, c2, c3, c4 = st.columns([2,1,1,2])
            c1.write(f"**{b.get('hgvsg','—')}**")
            c1.caption(f"{b.get('gene','—')} • {b.get('consequence','—')} • {b.get('impact','—')}")
            c2.write(f"REVEL: {b.get('revel')}")
            c3.write(f"CADD: {b.get('cadd')}")
            c4.write(f"AlphaMissense: {b.get('am_class')} ({b.get('am_score')})")
            note = st.text_input("Note", b.get("note",""), key=f"note_{i}")
            if st.button("Save note", key=f"save_{i}"):
                b["note"] = note
                st.success("Saved.")
            d1, d2 = st.columns(2)
            with d1:
                if st.button("Remove", key=f"rm_{i}"):
                    st.session_state.bookmarks = [x for x in st.session_state.bookmarks if x.get("hgvsg") != b.get("hgvsg")]
                    st.experimental_rerun()
            with d2:
                st.caption(f"Saved {b.get('timestamp')}")

# ---------------- Compare mode ----------------
def render_compare():
    st.title("Compare")
    if len(st.session_state.bookmarks) < 2:
        st.info("Bookmark at least two variants to compare."); return
    labels = [f"{b['hgvsg']} ({b.get('gene') or '—'})" for b in st.session_state.bookmarks]
    c1, c2 = st.columns(2)
    with c1: A = st.selectbox("Left", labels, key="cmp_left")
    with c2: B = st.selectbox("Right", labels, index=1, key="cmp_right")
    if A == B:
        st.warning("Pick two different variants."); return
    a = st.session_state.bookmarks[labels.index(A)]
    b = st.session_state.bookmarks[labels.index(B)]
    import pandas as pd
    def row(name, va, vb): return {"Field": name, "Left": va, "Right": vb}
    rows = [
        row("HGVS", a["hgvsg"], b["hgvsg"]),
        row("Gene", a.get("gene"), b.get("gene")),
        row("Consequence", a.get("consequence"), b.get("consequence")),
        row("Impact", a.get("impact"), b.get("impact")),
        row("REVEL", a.get("revel"), b.get("revel")),
        row("CADD", a.get("cadd"), b.get("cadd")),
        row("AlphaMissense", a.get("am_class"), b.get("am_class")),
        row("AlphaMissense score", a.get("am_score"), b.get("am_score")),
        row("Note", a.get("note",""), b.get("note","")),
    ]
    st.dataframe(pd.DataFrame(rows), use_container_width=True, hide_index=True)

# ---------------- Router ----------------
if mode == "Single":
    render_single()
elif mode == "Batch":
    render_batch()
elif mode == "Bookmarks":
    render_bookmarks()
else:
    render_compare()
