# do formulas get parsed correctly in rrma()?

library(testthat)

test_that("no intercept formula parses correctly", {
  d <- data.frame(ln_OR = rnorm(10)*1:10, x = 1:10, study_id =
    c(rep(1, 3), rep(2, 3), rep(3, 4)), var_eff = runif(10), 
    ind_1 = rep(gl(3, 1), 4)[-c(1, 2)])

  m_continous_ind <- rrma(ln_OR ~ x, d, study_id, var_eff, rho = 0.5)
  m_factor_ind <- rrma(ln_OR ~ ind_1, d, study_id, var_eff, rho = 0.5)
  m_no_intercept1 <- rrma(ln_OR ~ ind_1 -1, d, study_id, var_eff, rho = 0.5)
  m_no_intercept0 <- rrma(ln_OR ~ 0 + ind_1, d, study_id, var_eff, rho = 0.5)

  expect_true("intercept" %in% names(m_continous_ind$X_full))
  expect_true("x" %in% names(m_continous_ind$X_full))
  expect_true("intercept" %in% names(m_factor_ind$X_full))
  expect_true(!"intercept" %in% names(m_no_intercept1$X_full))
  expect_true(!"intercept" %in% names(m_no_intercept1$X_full))
  #expect_true(1 == 2) # blatantly false test
})

