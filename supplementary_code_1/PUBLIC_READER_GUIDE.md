# Public guide for Supplementary Code 1

This guide is for readers who want to inspect the public code release quickly without treating the synthetic example as manuscript data. The files show the public field layout and validation paths, then link reported manuscript values to the relevant public table or text location.

## Five-minute guide

1. Run the synthetic example:

```bash
python src/public_pipeline_template.py example_data/schema_example.csv --write-check public_check_output.json
```

2. Read the scope fields in `public_check_output.json`:

- `scope.toy_example_only` is `true`.
- `scope.exact_metrics_use_non_public_observation_table` is `true`.
- `public_output_boundary.toy_metrics_are_not_manuscript_performance` is `true`.
- `public_output_boundary.public_files` lists the public support files.
- `public_output_boundary.public_curation_table` points to Supplementary Data 1 Table S12.
- `metric_role.metrics_field` is `toy example metrics only`.
- `metric_role.manuscript_values_source` is `manuscript_metric_index.reported_metric_values and Supplementary Data 1`.
- `figure_source_map.figure_count` is 13.

3. Read the validation fields:

- `schema_check.rows` is 10 and `schema_check.polymer_groups` is 5.
- `split_check.grouped_cv.max_overlap` is 0, showing that grouped validation keeps polymer groups separated.
- `split_check.rowwise_control.may_mix_repeated_polymer_rows` is `true`, showing why the row-wise path is an optimistic control.
- `feature_family_check.nested_order` is `BASE`, `BASE_GENERIC`, and `BASE_GENERIC_AEM`.
- `schema_feature_crosswalk.feature_set_membership` lists the exact columns in the BASE, BASE_GENERIC, and BASE_GENERIC_AEM toy feature sets and confirms dictionary coverage.

4. Inspect the reported summary metrics:

- `manuscript_metric_index_public.json` lists the reported manuscript metrics and their Supplementary Data table locations.
- The same values are in Supplementary Data 1, especially Tables S4, S6, S8, S11, and S13.
- For source-held-out transfer, Table S6 reports raw performance and post hoc offset-size checks, while Table S9 and Fig. S8 report aggregate baseline-shift, rank-retention, and exploratory polymer-theme diagnostics.
- `example_output/public_example_output.json` shows the expected example JSON from the synthetic run, including the same scope and metric-role fields.

5. Use the claim and question files:

- `public_schema_dictionary.csv` identifies the public role and feature family of each exposed field.
- `claim_to_artifact_map_public.json` starts from a manuscript claim and lists the public files that inspect it.
- `reader_question_map_public.json` starts from a common methods or data question and lists the public files, manuscript locations, and interpretation limits that answer it.
- `figure_source_map_public.json` maps Fig. 1-Fig. 5 and Fig. S1-Fig. S8 to public files.

## What these public files show

- Public field layout and target definition for reported AEM hydroxide conductivity.
- `public_schema_dictionary.csv` as a compact column-to-role map for the public schema and derived toy features.
- `schema_feature_crosswalk` as a machine-readable check that the public schema dictionary covers every column used by the toy feature-family loop.
- Polymer-group validation path and zero-overlap group check on the synthetic example.
- Separate row-wise-control path used to illustrate why row-wise validation is optimistic for repeated polymer observations.
- Nested feature-family comparison path for BASE, BASE_GENERIC, and BASE_GENERIC_AEM.
- Comparison of reported summary manuscript metrics against Supplementary Data 1.
- Links for the source-held-out offset claim to Supplementary Data 1 Table S6, Fig. S8, `claim_to_artifact_map_public.json`, and `reader_question_map_public.json`.
- Figure-level guide from every main and supplementary figure to artwork files, manuscript/SI locations, Supplementary Data tables, public files, and scope notes.
- Summary of the paper's main interpretation limits: reported hydroxide-conductivity target, non-public observation-level table, small AEM descriptor increment, source-held-out offset check, and interpretive AEM response profiles.

## What remains outside the public files

- Non-public literature-derived rows, internal mapping metadata, full extraction tables, trained models, row-level predictions, and restricted article files remain outside the public release.
- CatBoost metric regeneration uses the observation-level table described in the manuscript Data availability statement.
- The synthetic toy metrics are code-example outputs, not membrane-science evidence.
- AEM response profiles remain fitted summaries, not causal morphology, free-volume, or transport-pathway measurements.
- Post hoc source-held-out offset checks remain offset-size summaries, not prospective performance estimates.

## Question file

Use `reader_question_map_public.json` as the shortest path through eight expected questions:

- dataset_not_public
- exact_metric_regeneration
- polymer_group_leakage
- rowwise_split_optimism
- small_aem_descriptor_increment
- support_distance_and_source_context
- response_profiles_not_causal
- supplementary_code_scope

Each entry gives a short answer, public files, and a scope note. Each entry identifies the relevant public files and the limit of the synthetic example.
