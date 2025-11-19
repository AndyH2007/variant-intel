from aa_props import add_basic_features

def build_features(df, colmap):
    X = add_basic_features(df, colmap)
    cols = ['delta_hydro','delta_charge','is_conservative']
    if 'label' in colmap and colmap['label'] in df.columns:
        return X[cols + [colmap['label']]]
    return X[cols]
