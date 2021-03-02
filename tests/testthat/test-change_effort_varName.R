context("change_effort_varName compliance")

# Data
library(geffaeR)
data(effort_example)
expected_col <- c("lon", "lat", "legId", "segId", "survey", "strateSec", "subRegion",
                  "strate", "seaState", "segLength", "subjective", "transect",
                  "left", "right", "CenterTime")


test_that("good output",{

  # all column needed present
  expect_equal(colnames(change_effort_varName(effort_example)), expected_col)

  # check if it builds strateSec with subRegion and strate
  effort_missing_strateSec <- effort_example[,-6]
  expect_equal(colnames(change_effort_varName(effort_missing_strateSec)), c(expected_col[-6], "strateSec"))
})




test_that("incorrect entries", {

  # not a data.frame
  mat <- matrix(c(rep(1,5),rep(3,5)), ncol = 2)
  expect_error(change_effort_varName(effort_base = mat),
               "is not a data frame")

  # missing "legId"
  effort_missing_leg <- effort_example[,-3]
  expect_error(change_effort_varName(effort_missing_leg),
               "can't find equivalent name")

  # # multiple segId
  # effort_example_multiple_segid <- effort_example
  # effort_example_multiple_segid$segid10km <- rep(c("0","1"), nrow(effort_example_multiple_segid)/2)
  # expect_error(change_effort_varName(effort_example_multiple_segid),
  #              "Can't determine a segId column because there are several segId columns")
})
