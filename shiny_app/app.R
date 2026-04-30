# app.R -----------------------------------------------------------------------
# Top-level Shiny app for the TCGA Survival Atlas.
# Run with:  shiny::runApp('~/tcga_survival/shiny_app')
# -----------------------------------------------------------------------------

source("global.R", local = FALSE)

ui <- function(request) {
  bslib::page_navbar(
    title  = tagList(tags$span(class = "brand-mark", "TCGA"),
                     tags$span(class = "brand-name", "Survival Atlas")),
    id     = "nav",
    theme  = APP_THEME,
    fillable = FALSE,
    navbar_options = bslib::navbar_options(bg = "#1f2d3d", theme = "dark"),
    header   = tagList(
      tags$head(
        tags$link(rel = "stylesheet", type = "text/css", href = "styles.css"),
        tags$script(src = "helpers.js")
      )
    ),
    nav_panel("Overview",          mod_overview_ui("ov"),
              icon = bsicons::bs_icon("compass")),
    nav_panel("Cohort atlas",      mod_cohort_atlas_ui("ca"),
              icon = bsicons::bs_icon("grid-3x3-gap")),
    nav_panel("Tumor vs normal",   mod_tumor_normal_ui("tn"),
              icon = bsicons::bs_icon("droplet-half")),
    nav_panel("Molecular structure", mod_molecular_structure_ui("ms"),
              icon = bsicons::bs_icon("diagram-3")),
    nav_panel("Survival modeling", mod_survival_modeling_ui("sm"),
              icon = bsicons::bs_icon("activity")),
    nav_panel("Cross-cohort",      mod_cross_cohort_ui("cc"),
              icon = bsicons::bs_icon("bar-chart-line")),
    nav_panel("Methods",           mod_methods_caveats_ui("mc"),
              icon = bsicons::bs_icon("book")),
    nav_spacer(),
    nav_item(tags$span(class = "nav-meta",
                       "v1 · cache:",
                       tags$code(basename(CACHE_DIR)))),
    nav_item(actionLink("._bookmark_btn",
                        label = tagList(bsicons::bs_icon("link-45deg"),
                                        " bookmark"),
                        class = "nav-meta",
                        title = "Copy a shareable URL of the current view"))
  )
}

server <- function(input, output, session) {
  mod_overview_server("ov", APP_DATA)
  mod_cohort_atlas_server("ca", APP_DATA)
  mod_tumor_normal_server("tn", APP_DATA)
  mod_molecular_structure_server("ms", APP_DATA)
  mod_survival_modeling_server("sm", APP_DATA)
  mod_cross_cohort_server("cc", APP_DATA)
  mod_methods_caveats_server("mc", APP_DATA)

  # Inputs we don't want serialized in the URL bookmark.
  setBookmarkExclude(c(
    "._bookmark_btn",
    "ca-meta_table_rows_selected", "ca-meta_table_state",
    "tn-tn_table_rows_selected",   "tn-tn_table_state",
    "cc-table_state",
    "sm-summary_table_state",
    "ms-markers_table_state"
  ))

  observeEvent(input$._bookmark_btn, {
    session$doBookmark()
  })

  # Used by the "Where to go next" cards on the Overview page to switch tabs.
  GOTO_LABELS <- c(
    ov = "Overview", ca = "Cohort atlas",
    tn = "Tumor vs normal", ms = "Molecular structure",
    sm = "Survival modeling", cc = "Cross-cohort", mc = "Methods"
  )
  observeEvent(input$._goto, {
    target <- input$._goto
    if (target %in% names(GOTO_LABELS)) {
      bslib::nav_select("nav", selected = GOTO_LABELS[[target]])
    }
  }, ignoreInit = TRUE)

  onBookmarked(function(url) {
    showModal(modalDialog(
      title = "Bookmark this view",
      tags$p("Copy the URL below — it encodes your current cohort, ",
             "modality, clustering, and other selections."),
      tags$input(type = "text", value = url, readonly = NA,
                 class = "form-control",
                 onclick = "this.select()"),
      easyClose = TRUE,
      footer    = modalButton("Close")
    ))
  })
}

shinyApp(ui, server, enableBookmarking = "url")
