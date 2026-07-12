# Paquetes mínimos para G-Formula (10.0 / 10.2) ----
# Evita dependencias del proyecto completo (chilemapas, dlnm, etc.).

gform_lib_path <- function() {
  lib <- Sys.getenv("R_LIBS_USER", unset = NA_character_)
  if (is.na(lib) || !nzchar(lib)) {
    lib <- path.expand("~/R/x86_64-pc-linux-gnu-library/4.3")
  }
  dir.create(lib, recursive = TRUE, showWarnings = FALSE)
  lib
}

install_load <- function(packages, lib = gform_lib_path()) {
  old_paths <- .libPaths()
  on.exit(.libPaths(old_paths), add = TRUE)
  .libPaths(c(lib, old_paths))

  for (pkg in packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      message("Instalando paquete: ", pkg)
      install.packages(
        pkg,
        lib = lib,
        repos = "https://cloud.r-project.org",
        dependencies = TRUE
      )
    }
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop("No se pudo cargar el paquete: ", pkg, call. = FALSE)
    }
    suppressPackageStartupMessages(
      library(pkg, character.only = TRUE, lib.loc = lib)
    )
  }
  invisible(NULL)
}

install_load(c(
  "knitr",
  "rio",
  "janitor",
  "openxlsx",
  "writexl",
  "dplyr",
  "tidyr",
  "stringr",
  "purrr",
  "magrittr",
  "rlang",
  "broom",
  "scales",
  "data.table",
  "survival",
  "future",
  "furrr",
  "future.apply",
  "parallel",
  "beepr",
  "tictoc"
))
