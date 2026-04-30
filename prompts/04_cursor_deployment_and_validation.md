Prepare the Shiny app in `~/tcga_survival/shiny_app/` for robust local use and future deployment.

# Goals

1. validate that the app works end-to-end locally
2. add lightweight validation and defensive checks
3. make deployment assumptions explicit
4. improve reproducibility and maintainability

# Scope

This is not the prompt to redesign the app. Focus on validation, runtime robustness, packaging assumptions, and deployment-readiness.

# Tasks

## 1. Runtime validation checks
Add lightweight checks for:
- required directories existing
- key summary files existing
- cache availability
- graceful fallback if cache absent
- cohort discovery logic
- figure path validity

Do not crash the app if a nonessential artifact is missing.
Instead:
- log/message clearly
- show informative UI empty state

## 2. Local run instructions
Ensure the README clearly documents:
- required R packages
- optional prep step:
  - `Rscript shiny_app/scripts/prepare_app_data.R`
- local launch command:
  - `shiny::runApp('~/tcga_survival/shiny_app')`
- expected data locations
- how omitted cohorts are handled

## 3. Deployment assumptions
Document deployment assumptions for future hosting:
- currently optimized for local filesystem paths under `~/tcga_survival` and `~/mnt/datapool/...`
- note what would need to change for:
  - Posit Connect
  - Shiny Server
  - containerized deployment

If reasonable, centralize path configuration in one place so it can later be parameterized.

## 4. Reproducibility / settings visibility
Where useful, make current settings visible to the user for clustering and key analyses:
- cohort
- modality
- clustering method
- number of clusters
- feature subset mode
- transform

This can be a simple “Current settings” box.

## 5. Download/output validation
If download handlers exist, verify they work and have sensible filenames.
Examples:
- benchmark table CSV
- cohort summary CSV
- cluster assignment CSV
- plot export if implemented

## 6. Known limitations documentation
Ensure documentation clearly states:
- `TCGA-LAML` excluded due to ingestion issue
- some molecular views may depend on available normal samples
- small cohorts / low event counts limit interpretation
- benchmark pages rely on precomputed model outputs

## 7. Lightweight testing mindset
Without building a huge test framework, add sanity helpers or checks for:
- benchmark summary load
- cohort metadata load
- one representative cohort loading on multiple pages
- one clustering-capable modality available

# Success criteria

A successful pass should ensure:
- app launches locally without brittle assumptions hidden in multiple files
- missing optional files do not catastrophically break the app
- docs are sufficient for another developer to run it
- future deployment constraints are clearly documented
- key settings and caveats are visible to users

# Deliverables

1. updated README
2. improved path/config handling if needed
3. validation helpers or defensive logic
4. summary of deployment/readiness improvements
5. list of any remaining blockers for production deployment
