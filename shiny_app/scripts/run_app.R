#!/usr/bin/env Rscript
# scripts/run_app.R -----------------------------------------------------------
# Tiny launcher for hosted environments. Reads host/port from env vars so the
# same script works in `Rscript run_app.R` on a laptop and inside a container.
#
#   PORT          (default: 3838)
#   HOST          (default: 0.0.0.0)
#   TCGA_SURVIVAL_APP_ROOT, TCGA_SURVIVAL_CACHE_DIR, TCGA_SURVIVAL_COHORT_ROOT
# -----------------------------------------------------------------------------

port <- as.integer(Sys.getenv("PORT", unset = "3838"))
host <- Sys.getenv("HOST", unset = "0.0.0.0")
app_root <- Sys.getenv("TCGA_SURVIVAL_APP_ROOT",
                        unset = file.path(Sys.getenv("HOME"),
                                          "tcga_survival", "shiny_app"))

if (!dir.exists(app_root))
  stop("App root does not exist: ", app_root)

setwd(app_root)
options(shiny.port = port, shiny.host = host)
cat(sprintf("[run_app] %s on http://%s:%d\n", app_root, host, port))
shiny::runApp(app_root, launch.browser = FALSE)
