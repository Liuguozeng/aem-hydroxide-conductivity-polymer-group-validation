# Polymer-group validation of hydroxide-conductivity models for anion exchange membranes: data leakage, descriptor interpretation, and source heterogeneity

This release contains Supplementary Code 1, Supplementary Data 1, the manuscript-companion figures, and the 48-profile AEM response-atlas extension.

Run the synthetic grouped-validation example from `supplementary_code_1`:

```bash
python src/public_pipeline_template.py example_data/schema_example.csv
python tests/test_public_pipeline_smoke.py
```

The synthetic example checks the public schema and validation code path. Exact manuscript metrics are indexed in `supplementary_code_1/manuscript_metric_index_public.json` and Supplementary Data 1; observation-level verification follows the article's Data availability statement. See `supplementary_code_1/REPRODUCIBILITY_BOUNDARY.md` for the release boundary.
