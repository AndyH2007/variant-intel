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
