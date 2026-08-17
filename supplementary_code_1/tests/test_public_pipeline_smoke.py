#!/usr/bin/env python3
"""Smoke test for Supplementary Code 1.

Run from the root of the unzipped Supplementary Code 1 package:

    python tests/test_public_pipeline_smoke.py

The test uses only the Python standard library. It checks that the public
example runs, loads the public metric/claim/question/figure/schema maps, and
keeps the public-curation scope explicit.
"""

from __future__ import annotations

import json
import math
import subprocess
import sys
import tempfile
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "src" / "public_pipeline_template.py"
EXAMPLE = ROOT / "example_data" / "schema_example.csv"
EXPECTED_OUTPUT = ROOT / "example_output" / "public_example_output.json"


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)



def equivalent(actual, expected, *, path="root"):
    """Compare schema exactly while tolerating platform-level float noise."""
    if isinstance(expected, bool) or expected is None or isinstance(expected, str):
        return actual == expected
    if isinstance(expected, int) and not isinstance(expected, bool):
        return type(actual) is int and actual == expected
    if isinstance(expected, float):
        return isinstance(actual, (int, float)) and not isinstance(actual, bool) and math.isclose(
            float(actual), expected, rel_tol=1e-12, abs_tol=1e-12
        )
    if isinstance(expected, list):
        return isinstance(actual, list) and len(actual) == len(expected) and all(
            equivalent(a, e, path=f"{path}[{index}]")
            for index, (a, e) in enumerate(zip(actual, expected))
        )
    if isinstance(expected, dict):
        return isinstance(actual, dict) and actual.keys() == expected.keys() and all(
            equivalent(actual[key], expected[key], path=f"{path}.{key}")
            for key in expected
        )
    return actual == expected


def main() -> None:
    with tempfile.TemporaryDirectory() as tmp:
        check_path = Path(tmp) / "public_check_output.json"
        result = subprocess.run(
            [
                sys.executable,
                str(SCRIPT),
                str(EXAMPLE),
                "--write-check",
                str(check_path),
            ],
            cwd=ROOT,
            check=True,
            capture_output=True,
            text=True,
        )
        output = json.loads(result.stdout)
        written = json.loads(check_path.read_text(encoding="utf-8"))

    require(equivalent(output, written), "stdout JSON and --write-check JSON should match")
    expected = json.loads(EXPECTED_OUTPUT.read_text(encoding="utf-8"))
    require(equivalent(output, expected), "fresh example JSON should match example_output/public_example_output.json")
    require(output.get("rows") == 10, "expected 10 toy rows")
    require(output.get("polymer_groups") == 5, "expected 5 toy polymer groups")
    require(output.get("target") == "log10(reported_hydroxide_conductivity_ms_cm)", "unexpected public target")

    scope = output.get("scope", {})
    require(scope.get("toy_example_only") is True, "toy-only scope flag missing")
    require(
        scope.get("exact_manuscript_metrics_use_non_public_observation_table") is True,
        "non-public observation-table metric route missing",
    )
    require(
        "non-public literature-derived rows" in scope.get("not_public_contents", []),
        "non-public literature-derived rows should be listed as absent",
    )

    boundary = output.get("public_output_boundary", {})
    require(
        boundary.get("do_not_cite_toy_metrics_as_manuscript_performance") is True,
        "toy metrics must be marked as non-manuscript performance",
    )
    require(boundary.get("toy_metrics_have_no_membrane_meaning") is True, "toy metric boundary missing")
    require(
        boundary.get("public_route_maps") == [
            "public_schema_dictionary.csv",
            "manuscript_metric_index_public.json",
            "claim_to_artifact_map_public.json",
            "reader_question_map_public.json",
            "figure_source_map_public.json",
        ],
        "public output boundary should list the five route-map files",
    )
    require(
        "Supplementary Data 1 Table S12" in boundary.get("public_curation_boundary", ""),
        "public curation boundary should point to Supplementary Data 1 Table S12",
    )
    require(
        "figure_source_map_public.json" in boundary.get("figure_source_map", ""),
        "public output boundary should point to the figure file",
    )
    metric_role = output.get("metric_role", {})
    require(metric_role.get("metrics_field") == "toy example metrics only", "metrics field role missing")
    require(
        metric_role.get("toy_primary_metrics_field") == "toy example summary only",
        "toy_primary_metrics role missing",
    )
    require(
        metric_role.get("manuscript_metric_source") == "manuscript_metric_index.reported_metric_values and Supplementary Data 1",
        "manuscript metric source role missing",
    )
    require(
        metric_role.get("exact_metric_regeneration") == "uses non-public observation table",
        "exact metric regeneration boundary missing",
    )

    grouped = output.get("split_check", {}).get("grouped_cv", {})
    rowwise = output.get("split_check", {}).get("rowwise_control", {})
    require(grouped.get("max_overlap") == 0, "grouped validation should have zero polymer-group overlap")
    require(
        rowwise.get("may_mix_repeated_polymer_rows") is True,
        "row-wise control should be flagged as potentially mixing repeated polymer rows",
    )

    expected_sets = ["BASE", "BASE_GENERIC", "BASE_GENERIC_AEM"]
    require(output.get("feature_sets") == expected_sets, "unexpected feature-set order")
    feature_check = output.get("feature_family_check", {})
    require(feature_check.get("nested_order") == expected_sets, "unexpected nested feature-family order")
    expected_full_features = [
        "reported_temperature_c",
        "iec_mmol_g",
        "thickness_um",
        "water_uptake_wt",
        "swelling_ratio_pct",
        "hydration_number",
        "inv_temperature_k",
        "structure_token_length_proxy",
        "lowercase_token_proxy",
        "uppercase_token_proxy",
        "separator_token_proxy",
        "composition_entropy_proxy",
        "aem_cation_to_anchor_descriptor",
        "free_volume_descriptor",
        "hydration_connectivity_descriptor",
        "transport_balance_descriptor",
    ]
    schema_feature_crosswalk = output.get("schema_feature_crosswalk", {})
    require(schema_feature_crosswalk.get("available") is True, "schema-feature crosswalk missing")
    require(
        schema_feature_crosswalk.get("feature_set_membership", {}).get("BASE_GENERIC_AEM") == expected_full_features,
        "schema-feature crosswalk should list the full feature-set columns",
    )
    require(
        schema_feature_crosswalk.get("dictionary_covers_feature_columns") is True,
        "schema dictionary should cover every toy feature column",
    )
    require(
        schema_feature_crosswalk.get("feature_columns_missing_from_dictionary") == [],
        "no toy feature column should be absent from the schema dictionary",
    )
    require(
        schema_feature_crosswalk.get("columns_missing_from_toy_input_or_features") == [],
        "schema dictionary should not list columns absent from the toy input/features",
    )

    metrics = output.get("metrics", {})
    for split_name in ["grouped_cv", "rowwise_control"]:
        require(split_name in metrics, f"missing {split_name} metrics")
        for feature_set in expected_sets:
            require(feature_set in metrics[split_name], f"missing {split_name}/{feature_set}")
            for key in ["r2", "rmse", "mae"]:
                require(key in metrics[split_name][feature_set], f"missing {split_name}/{feature_set}/{key}")

    reported_index = output.get("manuscript_metric_index", {})
    require(reported_index.get("available") is True, "manuscript metric index did not load")
    require(
        reported_index.get("scope", {}).get("copied_from_reported_manuscript_tables") is True,
        "manuscript metric index must identify values as copied from manuscript tables",
    )
    reported_values = reported_index.get("reported_metric_values", {})
    require(reported_values.get("polymer_group_full_descriptor_model", {}).get("r2") == 0.805, "reported full-descriptor-model R² missing")
    require(reported_values.get("polymer_group_full_descriptor_model", {}).get("rmse") == 0.128, "reported full-descriptor-model RMSE missing")
    require(reported_values.get("polymer_group_full_descriptor_model", {}).get("mae") == 0.096, "reported full-descriptor-model MAE missing")
    require(reported_values.get("polymer_group_base_generic", {}).get("r2") == 0.801, "reported BASE_GENERIC R² missing")
    require(
        reported_values.get("random_rowwise_controls", {}).get("base_generic_aem_descriptor_set_r2") == 0.941,
        "reported row-wise control R² missing",
    )
    require(
        reported_values.get("repeated_polymer_holdout_aem_descriptor_set", {}).get("win_rate_percent") == 75,
        "reported repeated-holdout win rate missing",
    )
    reader_links = output.get("reader_links", {})
    require(
        reader_links.get("recommended_tables") == ["S4", "S6", "S8", "S11", "S12", "S13"],
        "recommended table list should include S12 for the public-curation boundary and S13 for aggregate uncertainty checks",
    )
    claim_checks = output.get("manuscript_claim_checks", {})
    require(
        "Supplementary Data 1 Table S12" in claim_checks.get("public_curation_boundary", ""),
        "public-curation boundary should route to Supplementary Data 1 Table S12",
    )
    require(output.get("claim_to_artifact_map", {}).get("claim_count") == 9, "claim file should contain 9 entries")
    require(output.get("reader_question_map", {}).get("question_count") == 8, "question file should contain 8 entries")
    figure_source_map = output.get("figure_source_map", {})
    require(figure_source_map.get("available") is True, "figure file did not load")
    require(figure_source_map.get("figure_count") == 13, "figure file should cover five main and eight supplementary figures")
    require(figure_source_map.get("main_figure_count") == 5, "figure file should cover five main figures")
    require(figure_source_map.get("supplementary_figure_count") == 8, "figure file should cover eight supplementary figures")
    require(
        figure_source_map.get("all_entries_have_journal_upload_artifacts") is True,
        "each figure file entry should route to a journal upload artifact",
    )
    require(
        figure_source_map.get("journal_upload_artifact_route_count") == 14,
        "figure file should expose 14 journal upload artifact routes",
    )
    require("Fig. 5" in figure_source_map.get("figure_ids", []), "figure file should include Fig. 5")
    require("Fig. S8" in figure_source_map.get("figure_ids", []), "figure file should include Fig. S8")
    require(output.get("schema_dictionary", {}).get("entry_count", 0) >= 20, "schema dictionary is too small")

    for forbidden in ["r2", "rmse", "mae"]:
        require(forbidden not in output, f"{forbidden} should not appear as a top-level manuscript metric")

    print("PUBLIC_SMOKE_TEST_PASS rows=10 polymer_groups=5 group_leakage.max_overlap=0")


if __name__ == "__main__":
    main()
