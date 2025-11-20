# mutantscope_demo/app/acmg_helper.py
# Simple evidence grouping hints. NOT ACMG adjudication.
def acmg_flags_for_summary(s):
    flags = []
    # Pathogenic-leaning hints
    try:
        af = s.get("gnomad_af")
        if af is not None and float(af) <= 1e-6:
            flags.append(("P_PM2", "Very rare in population (PM2-like)"))
    except: pass
    am = s.get("alphamissense") or {}
    if (am.get("class") or "").lower() in {"pathogenic","likely_pathogenic"}:
        flags.append(("P_PP3", "In-silico support (AlphaMissense) (PP3-like)"))
    try:
        cadd = s.get("cadd_phred")
        if cadd is not None and float(cadd) >= 20:
            flags.append(("P_PP3b", "High CADD PHRED (PP3-like)"))
    except: pass
    try:
        revel = s.get("revel")
        if revel is not None and float(revel) >= 0.75:
            flags.append(("P_PP3c", "High REVEL (PP3-like)"))
    except: pass

    # Benign-leaning hints
    try:
        af = s.get("gnomad_af")
        if af is not None and float(af) >= 1e-3:
            flags.append(("B_BA1", "Too common for severe Mendelian (BA1-like)"))
    except: pass

    if not flags:
        flags.append(("INFO", "No strong heuristic flags — combine multiple evidence lines."))
    return flags
