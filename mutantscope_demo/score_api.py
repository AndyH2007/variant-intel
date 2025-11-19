import joblib, pathlib
MODEL = pathlib.Path(__file__).resolve().parent / 'models' / 'mutantscope.joblib'

def score_rows(X):
    model = joblib.load(MODEL)
    return model.predict_proba(X)
