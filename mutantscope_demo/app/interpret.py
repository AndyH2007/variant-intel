from typing import Any, Dict, List, Optional

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

def _to_str_list(items):
    out = []
    for it in items or []:
        if it is None:
            continue
        if isinstance(it, dict):
            cand = it.get("name") or it.get("description") or it.get("go_term")
            if cand is None:
                continue
            out.append(str(cand))
        else:
            out.append(str(it))
    # de-dupe, keep order
    seen, uniq = set(), []
    for s in out:
        if s not in seen:
            seen.add(s)
            uniq.append(s)
    return uniq

def _split_pubmed_values(values):
    ids = []
    for v in values or []:
        if isinstance(v, int):
            ids.append(str(v))
        elif isinstance(v, str):
            ids.extend([x.strip() for x in v.split(",") if x.strip()])
    # de-dupe, keep order
    seen, out = set(), []
    for x in ids:
        if x not in seen:
            seen.add(x)
            out.append(x)
    return out

def summarize_variant(vep_json: List[Dict[str, Any]]) -> Dict[str, Any]:
    if not vep_json:
        return {"error": "Empty VEP response."}
    v = vep_json[0]

    impact = _get_first(v, "impact")
    consequence = _get_first(v, "most_severe_consequence")
    hgvsg = _get_first(v, "input", "id")

    gene = None
    hgvsp = None
    prot_id = None
    phenotypes: List[str] = []
    clinvar_signifs: List[str] = []
    clinvar_reviews: List[str] = []
    pubmed_ids: List[str] = []
    domains: List[str] = []
    go_terms: List[str] = []

    tcs = v.get("transcript_consequences", []) or []
    if tcs:
        tc = tcs[0]
        gene = _get_first(tc, "gene_symbol", "genename")
        hgvsp = _get_first(tc, "hgvsp", "hgvsp_vep", "hgvsp_snpeff")
        prot_id = _get_first(tc, "protein_id", "ensembl_proteinid")

        for p in tc.get("phenotypes", []) or []:
            ph = _get_first(p, "phenotype", "description")
            if ph:
                phenotypes.append(str(ph))
            pm = _get_first(p, "pubmed_id")
            if pm:
                for token in str(pm).split(","):
                    tid = token.strip()
                    if tid:
                        pubmed_ids.append(tid)

        domains = _to_str_list(tc.get("domains"))
        go_terms = [s.replace("_", " ") for s in _to_str_list(tc.get("go"))]

    for cv in v.get("colocated_variants", []) or []:
        if "clin_sig" in cv and isinstance(cv["clin_sig"], list):
            for s in cv["clin_sig"]:
                if s:
                    clinvar_signifs.append(str(s).lower())
        for pmid in cv.get("pubmed", []) or []:
            pubmed_ids.append(str(pmid))

    if tcs:
        tc0 = tcs[0]
        cv_sig = _get_first(tc0, "clinvar_clnsig", "clinvar_clin_sig")
        if cv_sig:
            clinvar_signifs.append(str(cv_sig).lower())
        cv_rev = _get_first(tc0, "clinvar_review", "review_status")
        if cv_rev:
            clinvar_reviews.append(str(cv_rev))

    am = (tcs[0].get("alphamissense") if tcs else None) or {}
    alphamissense_class = am.get("am_class")
    alphamissense_score = _to_float(am.get("am_pathogenicity"))

    sift_pred = _get_first(tcs[0] if tcs else {}, "sift_prediction", "sift_pred")
    sift_score = _to_float(_get_first(tcs[0] if tcs else {}, "sift_score", "sift4g_score"))

    polyphen_pred = _get_first(tcs[0] if tcs else {}, "polyphen_prediction")
    polyphen_score = _to_float(_get_first(tcs[0] if tcs else {}, "polyphen_score"))

    cadd_phred = _to_float(_get_first(tcs[0] if tcs else {}, "cadd_phred"))
    revel_score = _to_float(_get_first(tcs[0] if tcs else {}, "revel", "revel_score"))

    af = _to_float(_get_first(tcs[0] if tcs else {}, "gnomad4.1_joint_af"))
    popmax_af = _to_float(_get_first(tcs[0] if tcs else {}, "gnomad4.1_joint_popmax_af"))

    gerp = _to_float(_get_first(tcs[0] if tcs else {}, "gerp++_nr", "gerp_91_mammals"))
    phastcons = _to_float(_get_first(tcs[0] if tcs else {}, "phastcons100way_vertebrate"))

    interactors: List[str] = []
    intact = v.get("intact") or v.get("IntAct")
    if intact and isinstance(intact, list):
        for it in intact[:20]:
            nm = _get_first(it, "interactor_xref", "interactor", "label", "name", "id")
            if nm:
                interactors.append(str(nm))

    def uniq(xs): return list(dict.fromkeys([x for x in xs if x]))

    return {
        "gene": gene,
        "hgvsg": hgvsg,
        "hgvsp": hgvsp,
        "consequence": consequence,
        "impact": impact,
        "protein_id": prot_id,

        "clinvar_significance": uniq(clinvar_signifs),
        "clinvar_review": uniq(clinvar_reviews),
        "phenotypes": uniq(phenotypes),

        "alphamissense": {"class": alphamissense_class, "score": alphamissense_score},
        "sift": {"pred": sift_pred, "score": sift_score},
        "polyphen": {"pred": polyphen_pred, "score": polyphen_score},
        "cadd_phred": cadd_phred,
        "revel": revel_score,

        "gnomad_af": af,
        "gnomad_popmax_af": popmax_af,

        "conservation": {"gerp": gerp, "phastcons100way": phastcons},
        "domains": uniq(domains),
        "go_terms": uniq(go_terms),
        "pubmed_ids": uniq(_split_pubmed_values(pubmed_ids)),
        "interactors": uniq(interactors),
    }
