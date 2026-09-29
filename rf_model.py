"""
Synthetic features and the Random Forest procedure used in the manuscript.

For each pair (f1, f2) two synthetic features are created:
  LR_f1__f2  log-ratio        log((f1 + eps) / (f2 + eps))
  GM_f1__f2  geometric mean   sqrt((f1 + eps) * (f2 + eps))

Random Forest: the number of features k is tuned on an internal 75/25 split of
the training data, the top-k features by RF importance are selected, and a
final RF trained on them is scored on the test data.
"""

import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import accuracy_score, roc_auc_score
from sklearn.model_selection import train_test_split

TOP_K_VALUES = [50, 100, 150, 200, 250, 300, 350, 400, 450, 500]
MODEL_PARAMS = {
    'n_estimators': 1000, 'max_features': 'sqrt', 'min_samples_leaf': 1,
    'class_weight': 'balanced', 'random_state': 42, 'n_jobs': -1
}
SYNTHETIC_PREFIXES = ('LR_', 'GM_')


def generate_synthetic_features(X, feature_pairs, epsilon=1e-9):
    """Return a DataFrame of LR_ and GM_ features for each (feature1, feature2) pair."""
    synthetic = {}
    for pair in feature_pairs:
        f1, f2 = sorted(pair)
        missing = [f for f in (f1, f2) if f not in X.columns]
        if missing:
            raise KeyError(f"MDFS pair feature(s) not found in the feature matrix: {missing}")
        synthetic[f"LR_{f1}__{f2}"] = np.log((X[f1] + epsilon) / (X[f2] + epsilon))
        synthetic[f"GM_{f1}__{f2}"] = np.sqrt((X[f1] + epsilon) * (X[f2] + epsilon))
    return pd.DataFrame(synthetic, index=X.index)


def find_best_k(X_train, y_train, k_values, model_params):
    """Choose the number of top-importance features k using an internal 75/25 split."""
    if X_train.shape[1] < min(k_values):
        return X_train.shape[1]

    X_int, X_val, y_int, y_val = train_test_split(
        X_train, y_train, test_size=0.25, random_state=42, stratify=y_train)

    ranking_model = RandomForestClassifier(**model_params).fit(X_int, y_int)
    sorted_features = pd.Series(ranking_model.feature_importances_,
                                index=X_int.columns).sort_values(ascending=False).index

    best_auc, best_k = -1, k_values[0]
    for k in sorted({min(k, len(sorted_features)) for k in k_values}):
        top_k = sorted_features[:k]
        model = RandomForestClassifier(**model_params).fit(X_int[top_k], y_int)
        auc = roc_auc_score(y_val, model.predict_proba(X_val[top_k])[:, 1])
        if auc > best_auc:
            best_auc, best_k = auc, k
    return best_k


def train_and_evaluate(X_train, y_train, X_test, y_test, k_values=None, model_params=None):
    """
    Tune k, select the top-k features by RF importance, train a final RF and score it.
    Returns auc, accuracy, best_k, selected_features and y_score (P(label = 1) on X_test).
    """
    k_values = k_values or TOP_K_VALUES
    model_params = model_params or MODEL_PARAMS

    best_k = find_best_k(X_train, y_train, k_values, model_params)

    ranking_model = RandomForestClassifier(**model_params).fit(X_train, y_train)
    importances = pd.Series(ranking_model.feature_importances_, index=X_train.columns)
    top_features = importances.nlargest(best_k).index.tolist()

    final_model = RandomForestClassifier(**model_params).fit(X_train[top_features], y_train)
    y_score = final_model.predict_proba(X_test[top_features])[:, 1]

    return {
        'auc': roc_auc_score(y_test, y_score),
        'accuracy': accuracy_score(y_test, final_model.predict(X_test[top_features])),
        'best_k': best_k,
        'selected_features': top_features,
        'y_score': y_score,
    }
