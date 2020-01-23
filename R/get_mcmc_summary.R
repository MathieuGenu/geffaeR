#' @export

get_mcmc_summary <- function(x, alpha = 0.2, dig = 4, median = FALSE) {

  if(!coda::is.mcmc(x)) {
    x <- coda::as.mcmc(x)
  }

  if(median) {
    summary <- round(c(median(x), mad(x), mad(x) / median(x), coda::HPDinterval(x, prob = 1 - alpha)), dig)
    summary <- data.frame(summary)
    row_name_one_time <- c("median","mad","mad/meadian",
                           paste("HPDlower at ",(alpha/2),sep=""),
                           paste("HPDupper at ",(1-alpha/2),sep=""))
    row.names(summary) <- rep(row_name_one_time)

  } else {

    summary <- round(c(mean(x), sd(x), sd(x) / mean(x), coda::HPDinterval(x, prob = 1 - alpha)), dig)
    summary <- data.frame(summary)
    row.names(summary) <- c("mean","sd","sd/mean",
                            paste("HPDlower at ",(alpha/2),sep=""),
                            paste("HPDupper at ",(1-alpha/2),sep=""))
  }

  return(summary)

}
