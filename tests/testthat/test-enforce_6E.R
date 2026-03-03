test_that("enforce_6E works", {

  seq1 <- "LSPVEEEAH" # "Good" sequence

  seq2 <- "LSPVETEAH" # "Good" sequence with the 6th aa changed to not E

  expect_equal(
    getITC_2024_06_27_cv(seq1)[1],
    7.121, #Original score
    tolerance=0.001
  )

  expect_equal(
    getITC_2024_06_27_cv(seq2)[1],
    25.329,
    tolerance=0.001
  )

  expect_equal(
    getITC_2024_06_27_cv(seq2, enforce_6E = TRUE)[1],
    95.0 # Default score returned when 6th position is not E.
  )

})
