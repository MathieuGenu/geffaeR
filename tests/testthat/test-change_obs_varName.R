context("change_obs_varName compliance")

# Data
library(geffaeR)
data(observation_example)
expected_col <- c("lon", "lat", "podSize", "segId", "taxon", "group", "family",
                  "species", "perpDist", "decAngle", "survey", "strateSec", "subRegion",
                  "strate", "transect", "side", "legId", "observerId")


test_that("good output",{

  # all column needed present
  expect_equal(colnames(change_obs_varName(observation_example)), expected_col)

  # # check if it builds strateSec with subRegion and strate
  obs_missing_strateSec <- observation_example[,-12]
  expect_equal(colnames(change_obs_varName(obs_missing_strateSec)), c(expected_col[-c(12,18)], "strateSec","observerId"))
})



test_that("incorrect entries", {

  # not a data.frame
  mat <- matrix(c(rep(1,5),rep(3,5)), ncol = 2)
  expect_error(change_obs_varName(obs_base = mat),
               "is not a data frame")

  # missing "legId"
  obs_missing_legid <- observation_example[,-18]
  expect_error(change_obs_varName(obs_missing_legid),
               "can't find equivalent name")

  # # multiple segId
  # observation_example_multiple_segid <- observation_example
  # observation_example_multiple_segid$segid10km <- rep(c("0","1"), nrow(observation_example_multiple_segid)/2)
  # expect_error(change_obs_varName(observation_example_multiple_segid),
  #              "Can't determine a segId column because there are several segId columns")
})
