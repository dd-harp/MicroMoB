test_that("multiplication works", {

  library(jsonlite)

  par <- list(
    p = 5,
    tmax = 10,
    aqua_path = system.file("extdata", "aqua_BH.json", package = "MicroMoB"),
    aqua_model = "BH",
    adult_path = system.file("extdata", "mosquito_RM.json", package = "MicroMoB"),
    adult_model = "RM"
  )

  json_path <- tempfile(pattern = "config", fileext = ".json")
  write_json(x = par, path = json_path, digits = NA, pretty = TRUE)

  par_in <- get_config_mosquito_MicroMoB(path = json_path)
  expect_true(all.equal(par, par_in))

  unlink(x = json_path)

})
