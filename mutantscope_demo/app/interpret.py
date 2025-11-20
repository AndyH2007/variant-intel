# app/interpret.py
from typing import Any, Dict, List, Optional
import re

def _get_first(d: dict, *keys, default=None):
    for k in keys:
        if k in d and d[k] not in (None, "", [], {}):
            return d[k]
    return default

def _to_float(x) -> Optional[float]:
    try:
        return float(str(x))
    except Exception:
        return None

def _uniq(xs):
    seen = set(); out = []
    for x in xs or []:
        if not x: continue
        if x in seen: continue
        seen.add(x); out.append(x)
    return out

def _parse_protein_pos(hgvsp: Optional[str]) -> Optional[int]:
    if not hgvsp: return None
    m = re.search(r'(\d{1,6})', str(hgvsp))
    try:
        return int(m.group(1)) if m else None
    except Exception:
        return None

def summarize_variant(vep_json: List[Dict[str, Any]]) -> Dict[str, Any]:
    if not vep_json:
        return {"error": "Empty VEP response."}
    v = vep_json[0]

    impact = _get_first(v, "impact")
    consequence = _get_first(v, "most_severe_consequence")
    hgvsg = _get_first(v, "input", "id")

    gene = hgvsp = prot_id = None
    phenotypes: List[str] = []
    clinvar_signifs: List[str] = []
    clinvar_reviews: List[str] = []
    pubmed_ids: List[str] = []
    domains: List[str] = []
    go_terms: List[str] = []

    tcs = v.get("transcript_consequences", []) or []
    tc0 = tcs[0] if tcs else {}

    if tcs:
        tc = tc0
        gene = _get_first(tc, "gene_symbol", "genename")
        hgvsp = _get_first(tc, "hgvsp", "hgvsp_vep", "hgvsp_snpeff")
        prot_id = _get_first(tc, "protein_id", "ensembl_proteinid")

        for p in tc.get("phenotypes", []) or []:
            ph = _get_first(p, "phenotype", "description")
            if ph: phenotypes.append(ph)
            pm = _get_first(p, "pubmed_id")
            if pm:
                for tok in str(pm).split(","):
                    tok = tok.strip()
                    if tok: pubmed_ids.append(tok)

        for ddom in tc.get("domains", []) or []:
            nm = _get_first(ddom, "name")
            if nm: domains.append(str(nm))

        for g in tc.get("go", []) or []:
            desc = _get_first(g, "description", "go_term")
            if desc: go_terms.append(str(desc).replace("_", " "))

    for cv in v.get("colocated_variants", []) or []:
        if "clin_sig" in cv and isinstance(cv["clin_sig"], list):
            clinvar_signifs.extend([s for s in cv["clin_sig"] if s])
        for pmid in cv.get("pubmed", []) or []:
            pubmed_ids.append(str(pmid))

    if tc0:
        cv_sig = _get_first(tc0, "clinvar_clnsig", "clinvar_clin_sig")
        if cv_sig: clinvar_signifs.append(str(cv_sig))
        cv_rev = _get_first(tc0, "clinvar_review", "review_status")
        if cv_rev: clinvar_reviews.append(str(cv_rev))

    # In-silico / AI
    am = (tc0.get("alphamissense") if tc0 else None) or {}
    alphamissense_class = am.get("am_class")
    alphamissense_score = _to_float(am.get("am_pathogenicity"))

    sift_pred = _get_first(tc0, "sift_prediction", "sift_pred")
    sift_score = _to_float(_get_first(tc0, "sift_score", "sift4g_score"))

    polyphen_pred = _get_first(tc0, "polyphen_prediction")
    polyphen_score = _to_float(_get_first(tc0, "polyphen_score"))

    cadd_phred = _to_float(_get_first(tc0, "cadd_phred"))
    revel_score = _to_float(_get_first(tc0, "revel", "revel_score"))

    # Population
    af = _to_float(_get_first(tc0, "gnomad4.1_joint_af"))
    popmax_af = _to_float(_get_first(tc0, "gnomad4.1_joint_popmax_af"))

    # Conservation
    gerp = _to_float(_get_first(tc0, "gerp++_nr", "gerp_91_mammals"))
    phastcons = _to_float(_get_first(tc0, "phastcons100way_vertebrate"))

    # Position — prefer explicit protein_start, else parse from HGVSp
    protein_pos = None
    ps = _get_first(tc0, "protein_start", "protein_end")
    if isinstance(ps, (int, float)) and ps > 0:
        protein_pos = int(ps)
    else:
        protein_pos = _parse_protein_pos(hgvsp)
    protein_len = None  # not provided in this payload

    return {
        "gene": gene,
        "hgvsg": hgvsg,
        "hgvsp": hgvsp,
        "consequence": consequence,
        "impact": impact,
        "protein_id": prot_id,
        "protein_pos": protein_pos,
        "protein_len": protein_len,

        "clinvar_significance": _uniq([s.lower() for s in clinvar_signifs]),
        "clinvar_review": _uniq(clinvar_reviews),
        "phenotypes": _uniq(phenotypes),

        "alphamissense": {"class": alphamissense_class, "score": alphamissense_score},
        "sift": {"pred": sift_pred, "score": sift_score},
        "polyphen": {"pred": polyphen_pred, "score": polyphen_score},
        "cadd_phred": cadd_phred,
        "revel": revel_score,

        "gnomad_af": af,
        "gnomad_popmax_af": popmax_af,

        "conservation": {"gerp": gerp, "phastcons100way": phastcons},
        "domains": _uniq(domains),
        "go_terms": _uniq(go_terms),
        "pubmed_ids": _uniq(pubmed_ids),
    }
