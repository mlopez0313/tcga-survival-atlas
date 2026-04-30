Polish the UI/UX and information architecture of the R Shiny app in `~/tcga_survival/shiny_app/`.

The app should **not** feel like a generic dashboard. It should feel like a **computational oncology exploration environment** with a guided narrative.

# Core goal

Refine the app so it communicates a coherent scientific story:
1. cohort landscape
2. tumor vs normal biology
3. molecular structure / clustering
4. survival consequences
5. cross-cohort benchmark findings

# Design principles

- Avoid business-dashboard aesthetics
- Avoid cluttered KPI-card overload
- Avoid tiny unreadable multi-panel grids
- Prioritize figure-first scientific communication
- Use space, typography, and annotation deliberately
- Emphasize narrative flow and interpretability

# Navigation / information architecture

Ensure the major navigation feels intentional and scientific.

Recommended major sections:
- Overview
- Cohort Atlas
- Tumor–Normal Biology
- Molecular Structure
- Survival Modeling
- Cross-Cohort Benchmark
- Methods & Caveats

Do not rename these to generic labels like “Data”, “Plots”, or “Dashboard”.

# Layout requirements

For most major pages, use a three-part structure where appropriate:
- left rail: controls only
- center: primary evidence/plots
- right rail: interpretation, caveats, selected-object summary

The right rail should not be decorative. It should help the user understand what they are seeing.

# Overview page

Make the landing page feel editorial.
It should include:
- a strong title and subtitle
- concise framing text
- a few high-value summary cards only
- benchmark heatmap
- win-count plot
- short key-conclusions panel
- clear note that `TCGA-LAML` is omitted for ingestion reasons
- clear “where to go next” entry points

# Interpretation panels

Where appropriate, add compact interpretation text such as:
- “No universal best model is observed across cohorts.”
- “This cohort has very few events; interpret the survival separation cautiously.”
- “Tumor–normal separation is visually strong in this modality.”
- “This cluster solution aligns with stage but not clearly with survival.”

These can be dynamically generated from the selected data where feasible, or templated conservatively.

# Plot styling

Apply a consistent visual language:
- use a restrained palette
- keep categorical colors consistent across pages when possible
  - tumor vs normal
  - cluster
  - stage
  - model family
- use readable font sizes
- avoid unnecessary gridlines
- use informative subtitles/captions where helpful
- keep legends clean and not overlarge

# Scientific visual preferences

Prefer:
- heatmaps with usable annotations
- box/violin plots with clean jitter overlays
- Kaplan–Meier curves with readable legends and risk info if feasible
- PCA/UMAP plots with clear encodings
- horizontal bar plots for ranked features

Avoid:
- pie charts
- gratuitous 3D plots
- overloaded all-in-one pages

# Consistency requirements

Make sure these are consistent across the app:
- cohort naming
- model naming
- modality naming
- excluded vs incomplete vs missing labels
- placement of download buttons
- placement of caveat text

# Empty states and warnings

Improve all empty/error states.
Examples:
- no normal samples available for this cohort/modality
- insufficient samples for clustering
- figure not available
- omitted cohort
- missing cache / missing modality

These should be informative and calm, not raw errors.

# Performance-aware UX

- show loading indicators for heavier pages
- avoid blocking the entire app when one module loads
- if clustering takes time, communicate that clearly
- progressively reveal complex outputs when available

# Download/export UX

If downloads exist, make them discoverable but not dominant.
Use clear labels like:
- Download table (CSV)
- Download cluster assignments
- Save current plot

# Accessibility / readability

- maintain good contrast
- avoid tiny side text
- ensure color is not the only grouping cue where practical
- use direct labels or helpful legends

# CSS/theming expectations

Use `bslib` and `www/styles.css` to create a cohesive scientific theme.
Potential style direction:
- light background
- strong typography
- subtle card borders
- muted accent colors
- publication-style headings

Do not make it look like a startup KPI dashboard.

# Success criteria

A successful UI/UX pass should result in:
- clearer narrative flow across pages
- less clutter
- more readable figures
- obvious scientific interpretation on each page
- consistent styling and labels
- a distinctive, non-generic feel

# Deliverables

1. update layout and styling files as needed
2. improve module UIs where needed
3. add or improve interpretation/caveat panels
4. summarize major UI/UX changes made
5. note any remaining UX limitations or future improvements
