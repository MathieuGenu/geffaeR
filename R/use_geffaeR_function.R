#' @export
#' @importFrom usethis use_template
use_geffaeR_template <- function(template,
                                 save_as = NULL,
                                 data = list(),
                                 ignore = FALSE,
                                 open = FALSE) {
  usethis::use_template(template = template,
                        save_as = save_as,
                        data = data,
                        ignore = ignore,
                        open = open,
                        package = "geffaeR")
}

#' @export
use_data_prep_template <- function(){


  usethis::use_directory("R")


  # data prep part
  use_geffaeR_template(
    template = "00_setup.R",
    save_as = "R/00_setup.R",
  )

  use_geffaeR_template(
    template = "01_change_colnames.R",
    save_as = "R/01_change_colnames.R",
  )

  use_geffaeR_template(
    template = "02_prepare_effort_and_obs.R",
    save_as = "R/02_prepare_effort_and_obs.R",
  )

}

#' @export
use_cds_template <- function(){

  usethis::use_directory("R/cds/function")

  # cds part
  use_geffaeR_template(
    template = "function_summarize_summary_ds.R",
    save_as = "R/cds/function/function_summarize_summary_ds.R",
  )

  use_geffaeR_template(
    template = "01_adjust_plot_detection.R",
    save_as = "R/cds/01_adjust_plot_detection.R",
  )

  use_geffaeR_template(
    template = "02_get_detection_info_and_graph.R",
    save_as = "R/cds/02_get_detection_info_and_graph.R",
  )

  use_geffaeR_template(
    template = "03_get_density_graph.R",
    save_as = "R/cds/03_get_density_graph.R",
  )

}
