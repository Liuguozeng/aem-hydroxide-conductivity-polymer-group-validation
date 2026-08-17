#!/usr/bin/env python3
"""Schema-compatible public example for reported AEM hydroxide conductivity modeling.

This file demonstrates public column names, feature-family construction,
grouped validation, row-wise optimistic-control validation, and leakage checks
using a synthetic toy CSV. It does not expose the restricted curated dataset,
source-derived curation materials, or exact manuscript metrics.
"""

from __future__ import annotations

import argparse
import json
import re
from collections.abc import Iterable
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.impute import SimpleImputer
from sklearn.linear_model import Ridge
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score
from sklearn.model_selection import GroupKFold, KFold
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler


TARGET = "reported_hydroxide_conductivity_ms_cm"
GROUP = "polymer_group"

REQUIRED_COLUMNS = [
    "polymer_group",
    "chemistry_family_group",
    "source_context_group",
    "reported_temperature_c",
    "iec_mmol_g",
    "thickness_um",
    "water_uptake_wt",
    "swelling_ratio_pct",
    "hydration_number",
    "reported_hydroxide_conductivity_ms_cm",
    "synthetic_structure_token_1",
    "synthetic_structure_token_2",
    "ratio_1",
    "ratio_2",
    "aem_cation_to_anchor_descriptor",
    "free_volume_descriptor",
    "hydration_connectivity_descriptor",
    "transport_balance_descriptor",
]

BASE_COLUMNS = [
    "reported_temperature_c",
    "iec_mmol_g",
    "thickness_um",
    "water_uptake_wt",
    "swelling_ratio_pct",
    "hydration_number",
    "inv_temperature_k",
]

GENERIC_COLUMNS = [
    "structure_token_length_proxy",
    "lowercase_token_proxy",
    "uppercase_token_proxy",
    "separator_token_proxy",
    "composition_entropy_proxy",
]

AEM_DESCRIPTOR_COLUMNS = [
    "aem_cation_to_anchor_descriptor",
    "free_volume_descriptor",
    "hydration_connectivity_descriptor",
    "transport_balance_descriptor",
]

FEATURE_SETS = {
    "BASE": BASE_COLUMNS,
    "BASE_GENERIC": BASE_COLUMNS + GENERIC_COLUMNS,
    "BASE_GENERIC_AEM": BASE_COLUMNS + GENERIC_COLUMNS + AEM_DESCRIPTOR_COLUMNS,
}


def validate_schema(df: pd.DataFrame) -> None:
    missing = [column for column in REQUIRED_COLUMNS if column not in df.columns]
    if missing:
        raise ValueError(f"Input CSV is missing required public-schema columns: {missing}")
    target = pd.to_numeric(df[TARGET], errors="coerce")
    if target.isna().any() or (target <= 0).any():
        raise ValueError("reported_hydroxide_conductivity_ms_cm must be positive and numeric for log10 transformation.")
    n_groups = df[GROUP].astype(str).nunique()
    if n_groups < 2:
        raise ValueError("Grouped validation requires at least two polymer_group values.")


def schema_check(df: pd.DataFrame) -> dict[str, object]:
    missing = [column for column in REQUIRED_COLUMNS if column not in df.columns]
    target = pd.to_numeric(df[TARGET], errors="coerce") if TARGET in df.columns else pd.Series(dtype=float)
    group_counts = df[GROUP].astype(str).value_counts().sort_index() if GROUP in df.columns else pd.Series(dtype=int)
    return {
        "rows": int(len(df)),
        "required_column_count": len(REQUIRED_COLUMNS),
        "present_required_columns": [column for column in REQUIRED_COLUMNS if column in df.columns],
        "missing_required_columns": missing,
        "numeric_target_rows": int(target.notna().sum()),
        "positive_target": bool(len(target) > 0 and target.notna().all() and (target > 0).all()),
        "polymer_groups": int(group_counts.size),
        "repeated_observation_groups": int((group_counts >= 2).sum()),
        "observations_per_polymer_group": {str(key): int(value) for key, value in group_counts.items()},
        "chemistry_family_groups": int(df["chemistry_family_group"].astype(str).nunique()) if "chemistry_family_group" in df.columns else 0,
        "source_context_groups": int(df["source_context_group"].astype(str).nunique()) if "source_context_group" in df.columns else 0,
    }


def weighted_structure_tokens(df: pd.DataFrame) -> list[str]:
    token_1 = df["synthetic_structure_token_1"].fillna("").astype(str)
    token_2 = df["synthetic_structure_token_2"].fillna("").astype(str)
    ratio_2 = pd.to_numeric(df["ratio_2"], errors="coerce").fillna(0.0)
    return [f"{a}.{b}" if r2 > 0 and b else a for a, b, r2 in zip(token_1, token_2, ratio_2, strict=True)]


def count_regex(values: Iterable[str], pattern: str) -> list[int]:
    rx = re.compile(pattern)
    return [len(rx.findall(value)) for value in values]


def composition_entropy(r1: pd.Series, r2: pd.Series) -> np.ndarray:
    ratios = np.vstack([
        pd.to_numeric(r1, errors="coerce").fillna(0.0).to_numpy(float),
        pd.to_numeric(r2, errors="coerce").fillna(0.0).to_numpy(float),
    ]).T
    totals = ratios.sum(axis=1, keepdims=True)
    totals[totals == 0] = 1.0
    probabilities = ratios / totals
    with np.errstate(divide="ignore", invalid="ignore"):
        entropy = -(probabilities * np.where(probabilities > 0, np.log(probabilities), 0.0)).sum(axis=1)
    return entropy


def build_public_features(df: pd.DataFrame) -> pd.DataFrame:
    features = pd.DataFrame(index=df.index)
    for column in [
        "reported_temperature_c",
        "iec_mmol_g",
        "thickness_um",
        "water_uptake_wt",
        "swelling_ratio_pct",
        "hydration_number",
        *AEM_DESCRIPTOR_COLUMNS,
    ]:
        features[column] = pd.to_numeric(df[column], errors="coerce")
    features["inv_temperature_k"] = 1.0 / (features["reported_temperature_c"] + 273.15)

    structure_tokens = weighted_structure_tokens(df)
    features["structure_token_length_proxy"] = [len(value) for value in structure_tokens]
    features["lowercase_token_proxy"] = count_regex(structure_tokens, r"[a-z]")
    features["uppercase_token_proxy"] = count_regex(structure_tokens, r"[A-Z]")
    features["separator_token_proxy"] = count_regex(structure_tokens, r"[_-]")
    features["composition_entropy_proxy"] = composition_entropy(df["ratio_1"], df["ratio_2"])
    return features


def leakage_summary(groups: np.ndarray, splitter: GroupKFold, X: pd.DataFrame, y: np.ndarray) -> dict[str, object]:
    overlaps: list[int] = []
    fold_sizes: list[dict[str, int]] = []
    for train, test in splitter.split(X, y, groups):
        train_groups = set(groups[train])
        test_groups = set(groups[test])
        overlap = train_groups.intersection(test_groups)
        overlaps.append(len(overlap))
        fold_sizes.append({"train_rows": int(len(train)), "test_rows": int(len(test)), "overlap_groups": int(len(overlap))})
    return {"folds": len(overlaps), "max_overlap": int(max(overlaps)), "fold_sizes": fold_sizes}


def rowwise_leakage_summary(groups: np.ndarray, y: np.ndarray) -> dict[str, object]:
    splitter = KFold(n_splits=min(5, len(y)), shuffle=True, random_state=42)
    overlaps: list[int] = []
    fold_sizes: list[dict[str, int]] = []
    for train, test in splitter.split(y):
        train_groups = set(groups[train])
        test_groups = set(groups[test])
        overlap = train_groups.intersection(test_groups)
        overlaps.append(len(overlap))
        fold_sizes.append({"train_rows": int(len(train)), "test_rows": int(len(test)), "overlap_groups": int(len(overlap))})
    return {
        "folds": len(overlaps),
        "max_overlap": int(max(overlaps)),
        "fold_sizes": fold_sizes,
        "may_mix_repeated_polymer_rows": bool(max(overlaps) > 0),
        "role": "optimistic leakage-prone control path",
    }


def feature_family_check(features: pd.DataFrame) -> dict[str, object]:
    required_features = sorted({column for columns in FEATURE_SETS.values() for column in columns})
    missing = [column for column in required_features if column not in features.columns]
    return {
        "nested_order": list(FEATURE_SETS),
        "feature_counts": {name: len(columns) for name, columns in FEATURE_SETS.items()},
        "missing_feature_columns": missing,
        "base_columns": BASE_COLUMNS,
        "generic_columns": GENERIC_COLUMNS,
        "aem_descriptor_columns": AEM_DESCRIPTOR_COLUMNS,
        "full_feature_count": len(required_features),
    }


def package_display_path(path: Path) -> str:
    package_root = Path(__file__).resolve().parents[1]
    try:
        return path.resolve().relative_to(package_root.resolve()).as_posix()
    except ValueError:
        return str(path)


def model_template() -> object:
    return make_pipeline(
        SimpleImputer(strategy="median"),
        StandardScaler(),
        Ridge(alpha=0.1),
    )


def score_predictions(y: np.ndarray, preds: np.ndarray) -> dict[str, float]:
    mask = np.isfinite(preds)
    return {
        "r2": float(r2_score(y[mask], preds[mask])),
        "rmse": float(np.sqrt(mean_squared_error(y[mask], preds[mask]))),
        "mae": float(mean_absolute_error(y[mask], preds[mask])),
    }


def load_metric_index(metric_index_path: str | None) -> dict[str, object]:
    default_path = Path(__file__).resolve().parents[1] / "manuscript_metric_index_public.json"
    path = Path(metric_index_path) if metric_index_path else default_path
    if not path.exists():
        return {
            "available": False,
            "expected_path": str(path),
            "message": "Reported summary metrics are provided separately in manuscript_metric_index_public.json and Supplementary Data 1.",
        }
    data = json.loads(path.read_text(encoding="utf-8"))
    reported_metrics = data.get("reported_metrics", {})
    reported_metric_values = {
        key: {subkey: subvalue for subkey, subvalue in value.items() if subkey != "location"}
        for key, value in reported_metrics.items()
        if isinstance(value, dict)
    }
    return {
        "available": True,
        "path": package_display_path(path),
        "scope": data.get("scope", {}),
        "reported_metric_keys": list(reported_metrics),
        "reported_metric_values": reported_metric_values,
        "reported_metric_locations": {
            key: value.get("location", "") for key, value in reported_metrics.items() if isinstance(value, dict)
        },
        "metric_table_links": data.get("metric_table_links", {}),
    }


def load_claim_map(claim_map_path: str | None) -> dict[str, object]:
    default_path = Path(__file__).resolve().parents[1] / "claim_to_artifact_map_public.json"
    path = Path(claim_map_path) if claim_map_path else default_path
    if not path.exists():
        return {
            "available": False,
            "expected_path": str(path),
            "message": "Public claim-to-file map is provided separately in claim_to_artifact_map_public.json.",
        }
    data = json.loads(path.read_text(encoding="utf-8"))
    claims = data.get("claims", [])
    return {
        "available": True,
        "path": package_display_path(path),
        "scope": data.get("scope", {}),
        "claim_ids": [claim.get("claim_id", "") for claim in claims if isinstance(claim, dict)],
        "claim_count": len(claims),
    }


def load_reader_question_map(question_map_path: str | None) -> dict[str, object]:
    default_path = Path(__file__).resolve().parents[1] / "reader_question_map_public.json"
    path = Path(question_map_path) if question_map_path else default_path
    if not path.exists():
        return {
            "available": False,
            "expected_path": str(path),
            "message": "Public question file is provided separately in reader_question_map_public.json.",
        }
    data = json.loads(path.read_text(encoding="utf-8"))
    questions = data.get("questions", [])
    return {
        "available": True,
        "path": package_display_path(path),
        "scope": data.get("scope", {}),
        "question_ids": [question.get("question_id", "") for question in questions if isinstance(question, dict)],
        "question_count": len(questions),
    }


def load_figure_source_map(figure_source_map_path: str | None) -> dict[str, object]:
    default_path = Path(__file__).resolve().parents[1] / "figure_source_map_public.json"
    path = Path(figure_source_map_path) if figure_source_map_path else default_path
    if not path.exists():
        return {
            "available": False,
            "expected_path": str(path),
            "message": "Public figure file is provided separately in figure_source_map_public.json.",
        }
    data = json.loads(path.read_text(encoding="utf-8"))
    figures = data.get("figures", [])
    figure_ids = [figure.get("figure_id", "") for figure in figures if isinstance(figure, dict)]
    main_ids = [figure_id for figure_id in figure_ids if figure_id.startswith("Fig. ") and not figure_id.startswith("Fig. S")]
    supplementary_ids = [figure_id for figure_id in figure_ids if figure_id.startswith("Fig. S")]
    journal_upload_artifact_route_count = sum(
        len(figure.get("journal_upload_artifacts", []))
        for figure in figures
        if isinstance(figure, dict)
    )
    return {
        "available": True,
        "path": package_display_path(path),
        "scope": data.get("scope", {}),
        "figure_ids": figure_ids,
        "figure_count": len(figures),
        "main_figure_count": len(main_ids),
        "supplementary_figure_count": len(supplementary_ids),
        "all_entries_have_journal_upload_artifacts": all(
            bool(figure.get("journal_upload_artifacts")) for figure in figures if isinstance(figure, dict)
        ),
        "journal_upload_artifact_route_count": journal_upload_artifact_route_count,
        "summary": data.get("summary", {}),
    }


def load_schema_dictionary(schema_dictionary_path: str | None) -> dict[str, object]:
    default_path = Path(__file__).resolve().parents[1] / "public_schema_dictionary.csv"
    path = Path(schema_dictionary_path) if schema_dictionary_path else default_path
    if not path.exists():
        return {
            "available": False,
            "expected_path": str(path),
            "message": "Public schema dictionary is provided separately in public_schema_dictionary.csv.",
        }
    data = pd.read_csv(path)
    required = {"column", "public_role", "feature_family", "description", "manuscript_use", "boundary"}
    missing = sorted(required - set(data.columns))
    if missing:
        raise ValueError(f"public_schema_dictionary.csv is missing required columns: {missing}")
    clean_data = data.fillna("").astype(str)
    feature_family_to_columns = {
        family: clean_data.loc[clean_data["feature_family"] == family, "column"].tolist()
        for family in sorted(clean_data["feature_family"].unique().tolist())
    }
    public_role_to_columns = {
        role: clean_data.loc[clean_data["public_role"] == role, "column"].tolist()
        for role in sorted(clean_data["public_role"].unique().tolist())
    }
    return {
        "available": True,
        "path": package_display_path(path),
        "entry_count": int(len(data)),
        "columns": clean_data["column"].tolist(),
        "public_roles": sorted(clean_data["public_role"].unique().tolist()),
        "feature_families": sorted(clean_data["feature_family"].unique().tolist()),
        "feature_family_to_columns": feature_family_to_columns,
        "public_role_to_columns": public_role_to_columns,
        "records": clean_data.to_dict(orient="records"),
        "non_public_table_note_rows": int(clean_data["boundary"].str.contains("non-public observation table", case=False).sum()),
    }


def schema_feature_crosswalk(df: pd.DataFrame, features: pd.DataFrame, schema_dictionary: dict[str, object]) -> dict[str, object]:
    feature_set_membership = {name: list(columns) for name, columns in FEATURE_SETS.items()}
    feature_columns = list(dict.fromkeys(column for columns in FEATURE_SETS.values() for column in columns))
    dictionary_columns = [str(column) for column in schema_dictionary.get("columns", [])]
    dictionary_column_set = set(dictionary_columns)
    toy_input_or_feature_columns = list(dict.fromkeys([*df.columns.astype(str).tolist(), *features.columns.astype(str).tolist()]))
    toy_input_or_feature_set = set(toy_input_or_feature_columns)
    feature_columns_missing_from_dictionary = [
        column for column in feature_columns if column not in dictionary_column_set
    ]
    columns_missing_from_toy_input_or_features = [
        column for column in dictionary_columns if column not in toy_input_or_feature_set
    ]
    return {
        "available": schema_dictionary.get("available") is True,
        "feature_set_membership": feature_set_membership,
        "feature_set_column_counts": {name: len(columns) for name, columns in FEATURE_SETS.items()},
        "feature_columns": feature_columns,
        "dictionary_columns": dictionary_columns,
        "dictionary_covers_feature_columns": feature_columns_missing_from_dictionary == [],
        "feature_columns_missing_from_dictionary": feature_columns_missing_from_dictionary,
        "columns_missing_from_toy_input_or_features": columns_missing_from_toy_input_or_features,
        "dictionary_columns_not_in_feature_sets": [
            column for column in dictionary_columns if column not in set(feature_columns)
        ],
        "feature_family_to_columns": schema_dictionary.get("feature_family_to_columns", {}),
        "public_role_to_columns": schema_dictionary.get("public_role_to_columns", {}),
    }


def evaluate_split(X: pd.DataFrame, y: np.ndarray, groups: np.ndarray, split: str, feature_columns: list[str]) -> dict[str, float]:
    if split == "grouped_cv":
        splitter = GroupKFold(n_splits=min(5, len(np.unique(groups))))
        iterator = splitter.split(X[feature_columns], y, groups)
    elif split == "rowwise_control":
        splitter = KFold(n_splits=min(5, len(y)), shuffle=True, random_state=42)
        iterator = splitter.split(X[feature_columns], y)
    else:
        raise ValueError(f"Unknown split: {split}")

    preds = np.full(len(y), np.nan)
    for train, test in iterator:
        model = model_template()
        model.fit(X.iloc[train][feature_columns], y[train])
        preds[test] = model.predict(X.iloc[test][feature_columns])
    return score_predictions(y, preds)


def evaluate(
    path: str,
    metric_index_path: str | None = None,
    claim_map_path: str | None = None,
    reader_question_map_path: str | None = None,
    figure_source_map_path: str | None = None,
    schema_dictionary_path: str | None = None,
) -> dict[str, object]:
    df = pd.read_csv(path)
    validate_schema(df)
    public_schema_check = schema_check(df)
    y = np.log10(pd.to_numeric(df[TARGET], errors="coerce").to_numpy(float))
    X = build_public_features(df)
    groups = df[GROUP].astype(str).to_numpy()
    schema_dictionary = load_schema_dictionary(schema_dictionary_path)
    group_splitter = GroupKFold(n_splits=min(5, len(np.unique(groups))))
    group_leakage = leakage_summary(groups, group_splitter, X, y)

    metrics: dict[str, dict[str, dict[str, float]]] = {"grouped_cv": {}, "rowwise_control": {}}
    for name, columns in FEATURE_SETS.items():
        metrics["grouped_cv"][name] = evaluate_split(X, y, groups, "grouped_cv", columns)
        metrics["rowwise_control"][name] = evaluate_split(X, y, groups, "rowwise_control", columns)

    primary = metrics["grouped_cv"]["BASE_GENERIC_AEM"]
    return {
        "rows": int(len(df)),
        "polymer_groups": int(len(np.unique(groups))),
        "target": "log10(reported_hydroxide_conductivity_ms_cm)",
        "scope": {
            "toy_example_only": True,
            "exact_manuscript_metrics_use_non_public_observation_table": True,
            "not_public_contents": [
                "non-public literature-derived rows",
                "internal mapping metadata",
                "full literature-extraction tables",
                "trained models",
                "row-level prediction files",
            ],
        },
        "public_output_boundary": {
            "do_not_cite_toy_metrics_as_manuscript_performance": True,
            "toy_metrics_have_no_membrane_meaning": True,
            "reported_metric_source": "Supplementary Data 1 Tables S4, S6, S8, S11, and S13 plus manuscript_metric_index_public.json and claim_to_artifact_map_public.json",
            "public_curation_boundary": "Supplementary Data 1 Table S12 records the public-curation protocol and reader-facing curation rationale.",
            "public_route_maps": [
                "public_schema_dictionary.csv",
                "manuscript_metric_index_public.json",
                "claim_to_artifact_map_public.json",
                "reader_question_map_public.json",
                "figure_source_map_public.json",
            ],
            "figure_source_map": "figure_source_map_public.json maps Fig. 1-Fig. 5 and Fig. S1-Fig. S8 to public source and upload-artifact routes.",
            "public_example_purpose": "schema, split, leakage-check, row-wise-control, and feature-family code-path inspection",
        },
        "metric_role": {
            "metrics_field": "toy example metrics only",
            "toy_primary_metrics_field": "toy example summary only",
            "manuscript_metric_source": "manuscript_metric_index.reported_metric_values and Supplementary Data 1",
            "exact_metric_regeneration": "uses non-public observation table",
        },
        "feature_sets": list(FEATURE_SETS),
        "estimator_scope": {
            "public_template_estimator": "Ridge(alpha=0.1) with median imputation and standard scaling",
            "manuscript_estimator": "CatBoost on the non-public observation table",
            "reason": "The public toy example demonstrates schema, split, leakage-check, row-wise-control, and feature-family logic; exact CatBoost metric recomputation uses the non-public observation table.",
        },
        "schema_check": public_schema_check,
        "split_check": {
            "grouped_cv": group_leakage,
            "rowwise_control": rowwise_leakage_summary(groups, y),
        },
        "feature_family_check": feature_family_check(X),
        "group_leakage": group_leakage,
        "metrics": metrics,
        "manuscript_metric_index": load_metric_index(metric_index_path),
        "claim_to_artifact_map": load_claim_map(claim_map_path),
        "reader_question_map": load_reader_question_map(reader_question_map_path),
        "figure_source_map": load_figure_source_map(figure_source_map_path),
        "schema_dictionary": schema_dictionary,
        "schema_feature_crosswalk": schema_feature_crosswalk(df, X, schema_dictionary),
        "reader_links": {
            "recommended_tables": ["S4", "S6", "S8", "S11", "S12", "S13"],
            "public_checks": [
                "validate_schema",
                "group_leakage.max_overlap",
                "metrics.grouped_cv",
                "metrics.rowwise_control",
                "feature_sets",
                "claim_to_artifact_map.claim_ids",
                "reader_question_map.question_ids",
                "figure_source_map.figure_ids",
                "schema_dictionary.public_roles",
                "schema_feature_crosswalk.feature_set_membership",
                "schema_check",
                "split_check",
                "feature_family_check",
            ],
        },
        "manuscript_claim_checks": {
            "polymer_group_validation": "Use group_leakage.max_overlap only to inspect the public grouped-split code path; use Supplementary Data 1 Table S4 for the manuscript held-out-key logic.",
            "feature_family_sensitivity": "Use metrics.grouped_cv only to inspect that the toy feature-family loop runs; use Supplementary Data 1 Table S6 and manuscript_metric_index.reported_metric_values for manuscript feature-family results.",
            "rowwise_optimistic_control": "Use metrics.rowwise_control only to inspect the toy row-wise-control path; use Supplementary Data 1 Table S11 and manuscript_metric_index.reported_metric_values for manuscript row-wise control values.",
            "reported_manuscript_metrics": "Use Supplementary Data 1 Tables S4, S6, S8, S11, and S13 plus manuscript_metric_index.reported_metric_values for manuscript metrics and aggregate uncertainty checks; this public toy run is not a metric-regeneration source.",
            "public_curation_boundary": "Use Supplementary Data 1 Table S12 for the public-curation protocol, reader-facing curation rationale, and aggregate-only release logic.",
            "figure_source_routing": "Use figure_source_map_public.json and figure_source_map.figure_ids to route each figure to artwork files, public tables, entries, and boundary notes.",
        },
        "primary_split": "grouped_cv",
        "primary_feature_set": "BASE_GENERIC_AEM",
        "toy_primary_metrics": primary,
    }


def write_check_json(result: dict[str, object], path: str) -> None:
    Path(path).write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("csv")
    parser.add_argument(
        "--metric-index",
        default=None,
        help="Optional path to manuscript_metric_index_public.json; defaults to the copy beside this script package.",
    )
    parser.add_argument(
        "--claim-map",
        default=None,
        help="Optional path to claim_to_artifact_map_public.json; defaults to the copy beside this script package.",
    )
    parser.add_argument(
        "--reader-question-map",
        default=None,
        help="Optional path to reader_question_map_public.json; defaults to the copy beside this script package.",
    )
    parser.add_argument(
        "--figure-source-map",
        default=None,
        help="Optional path to figure_source_map_public.json; defaults to the copy beside this script package.",
    )
    parser.add_argument(
        "--schema-dictionary",
        default=None,
        help="Optional path to public_schema_dictionary.csv; defaults to the copy beside this script package.",
    )
    parser.add_argument(
        "--write-check",
        dest="write_check",
        default=None,
        help="Optional path where the example JSON should be written in addition to stdout.",
    )
    args = parser.parse_args()
    result = evaluate(
        args.csv,
        args.metric_index,
        args.claim_map,
        args.reader_question_map,
        args.figure_source_map,
        args.schema_dictionary,
    )
    output = json.dumps(result, indent=2, sort_keys=True)
    if args.write_check:
        write_check_json(result, args.write_check)
    print(output)


if __name__ == "__main__":
    main()
