#' Create cds project
#'
#' @param path
#' @param ...
#'
#' @return Create a project environment to perform a cds analysis.
#' @importFrom usethis use_template create_project
#' @export
#'
#' @examples
create_cds_project <- function(path = ".",
                               rstudio = rstudioapi::isAvailable(),
                               open = TRUE){

  # create project
  path_norm <- normalizePath(path, mustWork = FALSE)
  # If calling this from a current project, reset it on exit
  old_proj <- get_proj()
  if (!is.null(old_proj) && path_norm != normalizePath(getwd())) {
    on.exit(usethis::proj_set(old_proj), add = TRUE)
  }

  create_proj(path = path, rstudio = rstudio)
  done("Creating new project")

  # create repositories and files
  all_files <- c("R/cds/function/","res/","data/")
  files_real_dir <- paste(path, all_files, sep = "/")


  lapply(files_real_dir, dir.create, recursive = TRUE, showWarnings = FALSE)

  use_data_prep_template()
  use_cds_template()

  if (open) open_project(path)

}

#' Get the path to the current project if it exists, otherwise return NULL
#' @noRd
get_proj <- function() {
  if (usethis:::is_package() || usethis:::possibly_in_proj(".")) {
    return(usethis::proj_get())
  }
  NULL
}

#' Create a project if one doesn't exist
#' @noRd
create_proj <- function(path = ".", rstudio) {
  if (!(usethis:::is_package(path) || usethis:::possibly_in_proj(path))) {
    usethis::create_project(path = path, open = FALSE, rstudio = rstudio)
  }
  usethis::proj_set(path, force = TRUE)
  if (rstudio) usethis::use_rstudio()
  invisible(TRUE)
}

#' Open a project if in RStudio
#' @noRd
open_project <- function(path) {
  if (rstudioapi::isAvailable() && interactive()) {
    rstudioapi::openProject(path, newSession = TRUE)
  } else if (normalizePath(path) != getwd()) {
    congrats("Your new project is created in ", path)
  }
  invisible(TRUE)
}

#' @import crayon
#' @export
done <- function(...) {
  cat(paste0(crayon::green(clisymbols::symbol$tick), " ", ...),
      "\n", sep = "")
}
