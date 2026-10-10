Osteoarthritis Transcriptomics Analysis (v1.0)

Associated Manuscript

Manuscript title: ANTXR1 drives ECM dysregulation in osteoarthritis by blocking integrin β1 mechanosensing

Journal: Signal Transduction and Targeted Therapy

Manuscript status: Submitted

Publication DOI: [To be added upon publication]

Note: The repository description does not, by itself, establish the analyses or results of a specific publication.

Release Contents
This release contains ten R scripts supporting osteoarthritis (OA) transcriptomic analysis, exploratory candidate-gene identification, diagnostic-model development, and model interpretation. The scripts are provided for methodological transparency and research reference. They are not packaged as a turnkey, one-command pipeline.

Repository Structure
Osteoarthritis-Transcriptomics-R-Scripts/
├── README.md
├── .gitignore
└── scripts/
├── 01_GEO_preprocessing.R
├── 02_differential_expression.R
├── 03_WGCNA.R
├── 04_gene_intersection.R
├── 05_GO_enrichment.R
├── 06_KEGG_enrichment.R
├── 07_LASSO_feature_selection.R
├── 08_machine_learning_diagnostic_models.R
├── 09_GSEA_input_preparation.R
└── 10_SHAP_model_interpretation.R

The scripts have been renamed for navigation. Their underlying code has not been refactored into a unified software package. The numbering indicates a conceptual order; some steps are separate branches and require manual preparation of intermediate files.

Not included: raw GEO datasets, platform files, patient-level clinical data, example datasets, generated figures, trained model files, or a validated automated workflow runner.

Conceptual Execution Order

Run 01_GEO_preprocessing.R after preparing the GEO expression table and platform annotation. Verify probe-to-gene mapping for the specific platform. Prepare the study-specific expression matrices and sample labels, including Sample Type Matrix.csv and the high-variability-gene matrix used for WGCNA.

Run 02_differential_expression.R and 03_WGCNA.R as separate analytical branches, then select the intended DEG and module gene lists.

Run 04_gene_intersection.R using only the intended gene-list files. Harmonize headers, delimiters, and gene-symbol formatting before downstream use.

Run 05_GO_enrichment.R and/or 06_KEGG_enrichment.R after checking the input format expected by each script.

Run 07_LASSO_feature_selection.R with the expression matrix, Control.txt, OA.txt, and the selected intersection-gene list.

Prepare train.csv, test.csv, and refer.txt for model construction. The Type column must encode control as 0 and OA as 1; training and validation samples must be independent. Run 08_machine_learning_diagnostic_models.R only after defining and verifying all analysis parameters, preprocessing, and nested cross-validation settings.

Run 10_SHAP_model_interpretation.R with the matching frozen model, preprocessing object, feature set, and cohort matrices. SHAP values describe predictive contributions and do not establish biological causality.

Separate GSEA preparation branch
09_GSEA_input_preparation.R prepares .gct, .cls, and expression-value files for downstream GSEA software. It does not execute GSEA itself.

Important: The scripts do not share a single working directory and cannot be assumed to run successfully without local path configuration, intermediate-file selection, and format harmonization.

Software Environment and Dependencies

Language: R

R version: The exact version used for the original analysis has not been documented in this release.

Python: Not required by the supplied R scripts.

Operating system: Several scripts contain Windows-specific paths. Update paths for the local environment. Cross-platform execution has not been verified.

Package versions: Exact versions are not pinned. No renv.lock is included.

The scripts use packages from CRAN and Bioconductor, including limma, clusterProfiler, org.Hs.eg.db, enrichplot, ComplexHeatmap, WGCNA, ggplot2, glmnet, caret, pROC, xgboost, randomForestSRC, kernelshap, shapviz, and other script-specific dependencies. Install only the packages required by the script being used, and record the tested environment with:
R.version.string
BiocManager::version()
writeLines(capture.output(sessionInfo()), "sessionInfo.txt")

Data Availability
Raw GEO expression data and platform annotations must be downloaded separately from the original sources, including the NCBI Gene Expression Omnibus (GEO), subject to the original data-use and citation conditions. No patient-level data or prepared study datasets are included in this release.

The analysis code refers to the following GEO accessions: GSE51588, GSE114007, GSE63359, GSE55457, and GSE12021. Before publication, verify that this list exactly matches the manuscript and document each accession's tissue type, cohort role, and source link.

Important Limitations

This release has not been independently validated as a complete end-to-end workflow.

Several scripts contain hard-coded paths, dataset-specific identifiers, and sample-name conventions that require adaptation.

Some scripts require manual intermediate-file preparation and input-format harmonization.

The supplied LASSO code includes fitted-data performance outputs and should not be interpreted as unbiased external validation unless the analysis is rerun with an appropriate resampling design.

The machine-learning analysis requires explicit configuration of cross-validation and related parameters before execution.

SHAP analysis requires a matching frozen model and preprocessing object; SHAP contributions are predictive explanations, not causal or mechanistic evidence.

The software is intended for exploratory biomedical research and methodological transparency. It is not a validated clinical diagnostic tool.

Citation
After Zenodo archives this GitHub release, cite the version-specific Zenodo DOI. Do not use a placeholder DOI as a real citation.

[Software author(s)]. (2026). Osteoarthritis Transcriptomics: R Analysis Scripts (Version v1.0.0) [Computer software]. Zenodo. [Add the version-specific DOI after archival].

If the software accompanies a publication, cite the corresponding article separately.

License
Add a LICENSE file before describing the repository as open-source. Choose a license only after confirming that all included code is yours to license and that any adapted third-party code permits redistribution. (Recommendation: MIT License)

Maintainer and Contact

Name: Fan Feng

Affiliation: The Fourth Military Medical University

Email: fengfan19910324@126.com

Issues: Use the GitHub Issues page for code-related questions, if enabled.
