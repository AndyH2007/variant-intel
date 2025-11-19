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
