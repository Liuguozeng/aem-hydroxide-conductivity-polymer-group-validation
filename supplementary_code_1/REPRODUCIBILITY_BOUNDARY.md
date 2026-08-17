# Reproducibility boundary for Supplementary Code 1

Supplementary Code 1 is a public code release for inspecting the field schema, grouped-validation logic, row-wise-control path, feature-family comparison pattern, and claim-to-file links. The literature-derived observation table follows the manuscript data-access route.

## What these public files show

- The public target definition: `log10(reported_hydroxide_conductivity_ms_cm)`.
- The grouped-validation key: `polymer_group`.
- The synthetic repeated-observation structure: 10 rows and 5 polymer groups, with two observations per group.
- The grouped split result: `group_leakage.max_overlap=0`.
- The row-wise optimistic-control path and its ability to mix repeated polymer rows.
- The nested feature-family order: `BASE`, `BASE_GENERIC`, and `BASE_GENERIC_AEM`.
- The public schema roles listed in `public_schema_dictionary.csv`.
- The field-to-feature-set list showing exact feature columns and public-dictionary coverage.
- The reported summary manuscript metrics indexed in `manuscript_metric_index_public.json`.
- The exact fixed benchmark settings and training-fold preprocessing contract in `fixed_model_configuration_public.json`.
- The claim file lists public inspection points for each major manuscript claim in `claim_to_artifact_map_public.json`.
- The question file links common methods or data questions to public files and manuscript locations in `reader_question_map_public.json`.
- The figure-level source guide from each main and supplementary figure to artwork files, public tables, map files, and boundary notes in `figure_source_map_public.json`.

## What remains outside this release

- Exact CatBoost metric recomputation from the observation-level table.
- Observation-level rows from the curated literature table.
- Source-derived curation materials and complete source-group identifiers.
- Full literature-extraction tables, source article files, local reference-manager files, caches, trained models, or row-level prediction files.

## Why the release is limited

The non-public observation table is derived from article-level extraction records and internal mapping metadata. Those rows contain extraction records and source-group assignments that are handled through the manuscript data-access route. The manuscript provides public summary metrics, curation rules, descriptor definitions, Supplementary Data summary tables, and this synthetic-schema code package.

## How to read the public files

1. Run `python src/public_pipeline_template.py example_data/schema_example.csv --write-check public_check_output.json`.
2. Run `python tests/test_public_pipeline_smoke.py`.
3. Inspect `example_output/public_example_output.json` to see the expected example JSON structure.
4. Inspect `manuscript_metric_index_public.json` for the reported values and their Supplementary Data table locations.
5. Use `claim_to_artifact_map_public.json` when starting from a manuscript claim.
6. Use `reader_question_map_public.json` when starting from a common methods or data question.
7. Use `figure_source_map_public.json` when starting from a figure.
8. Use Supplementary Data 1 Tables S4, S6, S8, S11, S12, and S13 for the locked validation, descriptor-family sensitivity, repeated-holdout, row-wise-control, public-curation, and uncertainty and split-size summaries.

The figure-level source guide links each figure to its public files and reported summary tables.

Toy metrics printed by the public example are deterministic example outputs. The JSON `metric_role` block labels them as toy-only outputs and points readers to `manuscript_metric_index.reported_metric_values` plus Supplementary Data 1 for manuscript metrics. They have no membrane-performance meaning and are not manuscript evidence.
