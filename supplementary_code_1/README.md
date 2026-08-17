# Supplementary Code 1

Supplementary Code 1 demonstrates the public field layout, polymer-group split, row-wise control, and feature-family loop on a synthetic 10-row dataset. The manuscript Data availability statement describes access to the observation-level literature table.

## Contents

- `README.md`: this run guide.
- `PUBLIC_READER_GUIDE.md`: a short guide to the runnable example, manuscript metric index, and claim, question, and figure maps.
- `REPRODUCIBILITY_BOUNDARY.md`: a concise statement of what the public package can check and which manuscript values use the non-public observation-level table.
- `environment_report_public.json`: index of the separate public-toy, historical analysis-reference, submission-build, and current validation-runtime records.
- `example_data/schema_example.csv`: a 10-row, 5-polymer-group toy dataset with repeated observations per polymer.
- `summary_data/aem_response_profile_coordinates_public.csv`: aggregate fitted coordinates for the eight displayed AEM response profiles; no observation-level rows are included.
The response summary reports `training_support_n` from the reference profile and `plotted_support_n` for the valid observations inside the plotted 5th–95th percentile range.
- `fixed_model_configuration_public.json`: machine-readable outer-fold preprocessing, feature-budget, weighting, hyperparameter-grid, and fixed-model settings for all eleven benchmark families.
- `aem_descriptor_specification_public.csv` and `.json`: formula-level definitions for all 78 computed AEM descriptors.
- `src/aem_descriptor_reference.py`: deterministic RDKit reference implementation that accepts user-supplied repeat-unit SMILES and composition weights. The archive does not include the manuscript's restricted polymer structures.
- `public_schema_dictionary.csv`: a machine-readable dictionary of public input fields, derived toy features, feature-family roles, manuscript use, and release-scope notes.
- `requirements.txt`: version-bounded minimal packages needed to run the public toy example.
- `LICENSE.md`: license terms for the public code, public summary tables, and synthetic example data.
- `manuscript_metric_index_public.json`: public reported summary metrics and the Supplementary Data locations where they can be checked.
- `claim_to_artifact_map_public.json`: a public claim file linking manuscript claims to the code example, manuscript locations, and Supplementary Data tables that can be inspected.
- `reader_question_map_public.json`: a question file that links common methods and data-availability questions to the public files that can be checked.
- `figure_source_map_public.json`: a figure-level source map for Fig. 1-Fig. 5 and Fig. S1-Fig. S8; it maps each figure to artwork files, `journal_upload_artifacts`, manuscript/SI captions, Supplementary Data tables, entries, and boundary notes.
- `example_output/public_example_output.json`: deterministic output from the public toy command for readers who want to inspect the expected JSON structure before running the script.
- `src/public_pipeline_template.py`: a grouped-validation template that validates the public schema, checks polymer-group leakage, runs a row-wise control split, and reports toy descriptor-family sensitivity metrics.
- `tests/test_public_pipeline_smoke.py`: a standard-library smoke test that runs the example and checks the expected example fields.

## Toy schema

The example CSV contains the public-facing column pattern used by the example. For a compact machine-readable version of these roles, including derived toy features, see `public_schema_dictionary.csv`.

- `polymer_group`: grouped-validation key used for polymer-level holdout.
- `chemistry_family_group`: example chemistry-family grouping field.
- `source_context_group`: synthetic source-context grouping field.
- `reported_temperature_c`: reported measurement temperature.
- `iec_mmol_g`: ion exchange capacity.
- `thickness_um`: membrane thickness.
- `water_uptake_wt`, `swelling_ratio_pct`, `hydration_number`: example membrane/test covariates.
- `reported_hydroxide_conductivity_ms_cm`: reported hydroxide-conductivity target for the toy run.
- `synthetic_structure_token_1`, `synthetic_structure_token_2`, `ratio_1`, `ratio_2`: synthetic structure-composition fields.
- `aem_cation_to_anchor_descriptor`, `free_volume_descriptor`, `hydration_connectivity_descriptor`, `transport_balance_descriptor`: synthetic AEM descriptor fields used only to demonstrate the public feature-family pattern.

Run the toy example:

```bash
python src/public_pipeline_template.py example_data/schema_example.csv
```

The expected toy output includes `rows=10`, `polymer_groups=5`, `group_leakage.max_overlap=0`, and `toy_primary_metrics` for the full toy feature set. The nested `metrics` block reports toy grouped-CV and row-wise-control results for `BASE`, `BASE_GENERIC`, and `BASE_GENERIC_AEM`. The output reports the toy schema, grouped split, row-wise control, feature-family order, manuscript metric index, and public links to the manuscript tables. Ridge is the lightweight example estimator; CatBoost is the manuscript model.

To save the printed example JSON:

```bash
python src/public_pipeline_template.py example_data/schema_example.csv --write-check public_check_output.json
```

The saved JSON repeats the same public-only fields printed to the terminal; manuscript performance values are listed in the manuscript and Supplementary Data 1.

Optional package smoke test:

```bash
python tests/test_public_pipeline_smoke.py
```

The smoke test runs the same public example, verifies that the toy dataset has 10 rows and 5 polymer groups, checks zero overlap in grouped validation, confirms that row-wise validation is flagged as an optimistic control, and confirms that the manuscript metric index, claim file, question file, figure file, and schema dictionary load from the package. It is included to check public-package consistency; manuscript CatBoost metrics are listed in Supplementary Data 1.
It also compares the freshly generated example JSON with `example_output/public_example_output.json`, so the expected public output structure is inspectable even before running the code.

## What to inspect after running

Start with `rows=10`, `polymer_groups=5`, and `group_leakage.max_overlap=0` to confirm the synthetic grouped split. Then compare `metrics.grouped_cv` with `metrics.rowwise_control` to see that the two validation paths are separate. The `metric_role` and `scope` fields mark the toy outputs as example-only and point readers to Supplementary Data 1 and `manuscript_metric_index_public.json` for the reported manuscript values.

## Reported manuscript values

These values are copied from the manuscript tables, not generated by the toy run:

- Polymer-group full descriptor model: R² = 0.805, RMSE = 0.128, MAE = 0.096 (Supplementary Data 1, Table S6).
- BASE + GENERIC polymer-group model: R² = 0.801 (Supplementary Data 1, Table S6).
- Random row-wise controls: R² = 0.940 and R² = 0.941 (Supplementary Data 1, Table S11).
- AEM descriptors repeated-holdout summary: median ΔR² = +0.005 and win rate = 75% (Supplementary Data 1, Table S8).
- Source-held-out raw performance and baseline-shift diagnostics are reported in Tables S6 and S9; the post hoc mean-shift and affine values in Table S6 are offset-size checks, not prospective performance.
- Uncertainty and split-size summary: primary-model R² = 0.805 with 95% group-bootstrap CI 0.772 to 0.835; field-availability strata and split-size ranges are summarized in Supplementary Data 1, Table S13.

`manuscript_metric_index_public.json` repeats these reported values in machine-readable form. `claim_to_artifact_map_public.json`, `reader_question_map_public.json`, and `figure_source_map_public.json` let readers start from a claim, a methods question, or a figure and then find the public file, manuscript/SI location, Supplementary Data table, and scope note.

## How to read the package

- The public schema contains the fields needed to define the reported hydroxide-conductivity target, polymer-group key, source-context key, membrane/test covariates, generic structure descriptors, and AEM descriptor fields. The same roles are listed in `public_schema_dictionary.csv` so readers can inspect the column-to-claim mapping without reading the Python script.
- Polymer-group validation keeps each polymer group entirely on one side of every fold. The toy run reports `group_leakage.max_overlap=0`; the corresponding manuscript logic is summarized in Table S4, and aggregate fold-size ranges are summarized in Table S13.
- The row-wise control split is implemented separately from polymer-group validation. The real row-wise optimistic check is reported as summary metrics in Table S11; the toy data only demonstrate the code path.
- The feature-family loop evaluates `BASE`, `BASE_GENERIC`, and `BASE_GENERIC_AEM` in the same nested order used for the public manuscript narrative. The toy metrics confirm only that the comparison path runs; the manuscript values are reported in Table S6.
- The public script uses Ridge as a deterministic template estimator so it can run with only `numpy`, `pandas`, and `scikit-learn`. The manuscript benchmark uses CatBoost; the reported CatBoost metrics are indexed in `manuscript_metric_index_public.json` and Supplementary Data 1.
- The toy grouped-CV ordering is illustrative. It is not evidence for the manuscript AEM descriptor interpretation; the manuscript evidence is the reported validation, descriptor-family sensitivity, repeated-holdout, and row-wise-control metrics reported in Supplementary Data 1.
- The `metric_role` field labels `metrics` and `toy_primary_metrics` as toy example outputs. The top-level `toy_primary_metrics` field is included only to confirm that the example command runs. It must not be cited as manuscript performance.
- The `figure_source_map_public.json` map covers all five main figures and eight supplementary figures. It maps each figure to artwork files and public files while noting that exact row-level plotting data use the observation-level table.
- The observation-level literature extraction table is handled through the manuscript data-access route. The package also excludes internal mapping metadata, literature-extraction tables, trained models, and row-level prediction files.

## License

Supplementary Code 1 is released under MIT License (code) and CC BY 4.0 (public summary tables and synthetic example data). The license file also states that observation-level extraction records, source-derived curation materials, complete identifiers, trained models, row-level predictions, and third-party article full texts are outside the public release.

## Dependencies

The template requires Python with `numpy`, `pandas`, and `scikit-learn`; compatible version ranges for these minimal run dependencies are listed in `requirements.txt`. The public package records the toy runtime, the locked analysis package versions, the submission-build runtime, and the current independent validation runtime in separate files. The last two are build and verification records, not substitutes for the historical model-generation environment.

## Data release note

This package supports inspection of the field schema, polymer-group validation implementation, row-wise control implementation, and descriptor-family sensitivity pattern. Exact reported metrics use the non-public observation-level table described in the manuscript Data availability statement and Supplementary Data 1 Table S12.
