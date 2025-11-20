# mutantscope_demo/app/protein_track.py
# Returns a lightweight SVG ruler with the variant's AA position marked.
def protein_track_svg(summary):
    hgvsp = summary.get("hgvsp") or ""
    # Expect something like "ENSP...:p.Val752Met" or "p.Val752Met"
    aa_pos = None
    for token in [hgvsp]:
        if not token: continue
        if "p." in token:
            part = token.split("p.",1)[1]
            # find the first integer sequence
            num = ""
            for ch in part:
                if ch.isdigit(): num += ch
                elif num: break
            if num:
                try: aa_pos = int(num)
                except: pass
    if aa_pos is None:
        return None

    # simple ruler: 0..max (normalize to 1200 for layout)
    maxlen = max(aa_pos+50, 100)  # unknown protein length; show context
    width = 1000; pad = 40
    x = int((aa_pos / maxlen) * width)

    # Domains (if any) as chips rendered below (not placed to exact coords)
    doms = (summary.get("domains") or [])[:8]
    chip_row = " ".join([f'<span style="border:1px solid #ccc;border-radius:10px;padding:2px 6px;margin-right:6px;font-size:12px">{d}</span>' for d in doms])

    svg = f"""
    <div style="font-family:Inter,Arial; color:#e6edf3; background:#161b22; padding:10px;border-radius:12px;border:1px solid rgba(255,255,255,.08);">
      <div style="font-size:12px;opacity:.8;margin-bottom:6px">Protein ruler (approximate context)</div>
      <svg width="{width+pad*2}" height="60" xmlns="http://www.w3.org/2000/svg">
        <rect x="{pad}" y="25" width="{width}" height="8" rx="4" fill="#2a303b" stroke="rgba(255,255,255,.08)"/>
        <line x1="{pad+x}" y1="18" x2="{pad+x}" y2="40" stroke="#58a6ff" stroke-width="2"/>
        <text x="{pad+x+6}" y="16" fill="#e6edf3" font-size="12">aa {aa_pos}</text>
      </svg>
      <div style="margin-top:6px">{chip_row}</div>
    </div>
    """
    return svg
