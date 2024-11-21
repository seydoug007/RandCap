test_that("RandCapProd maintains stratification and block sizes but alters treatment sequence", {
  # Generate a randomization object with RandCapGen
  rand_gen <- RandCapGen(
    n = 12,
    id_prefix = "ID_",
    block_prefix = "Block_",
    block_sizes = c(4, 8),
    arms = c("Control", "Treatment"),
    strat_vars = list(Sex = c("M", "F"),
                      AgeGroup=c("Young","Old")),
    seed = 2024
  )

  # Shuffle the randomization object with RandCapProd
  rand_prod <- RandCapProd(rand_gen, seed = 5678)

  # Verify that the stratification variable values are consistent, ignoring order
  gen_strata <- unique(rand_gen$tables$full_dataset$stratum)
  prod_strata <- unique(rand_prod$tables$full_dataset$stratum)
  expect_false(
    identical(gen_strata, prod_strata),
    info = "Strata order has not changed between RandCapGen and RandCapProd."
  )

  # Verify that the number of participants per arm in each block is consistent
  counts_gen <- rand_gen$tables$full_dataset %>%
    dplyr::group_by(strat_Sex, block.id, treatment_arm) %>%
    dplyr::count()

  counts_prod <- rand_prod$tables$full_dataset %>%
    dplyr::group_by(strat_Sex, block.id, treatment_arm) %>%
    dplyr::count()

  expect_equal(
    counts_gen %>% dplyr::arrange(strat_Sex, block.id, treatment_arm),
    counts_prod %>% dplyr::arrange(strat_Sex, block.id, treatment_arm),
    info = "Participant counts per treatment arm and block are not consistent."
  )

  # Verify that treatment sequences differ in each stratum
  seq_gen <- rand_gen$tables$full_dataset %>%
    dplyr::group_by(strat_Sex) %>%
    dplyr::summarize(sequence = paste(treatment_arm, collapse = "-"))

  seq_prod <- rand_prod$tables$full_dataset %>%
    dplyr::group_by(strat_Sex) %>%
    dplyr::summarize(sequence = paste(treatment_arm, collapse = "-"))

  seq_diff <- seq_gen$sequence != seq_prod$sequence
  expect_true(
    all(seq_diff),
    info = "Treatment sequences in each stratum are not different between RandCapGen and RandCapProd."
  )
})


test_that("RandCapProd generates the same object after repeated processes", {
  # Generate the first randomization object with RandCapGen
  rand_gen_1 <- RandCapGen(
    n = 20,
    id_prefix = "ID_",
    block_prefix = "Block_",
    block_sizes = c(10, 10),
    arms = c("Control", "Treatment"),
    strat_vars = list(Sex = c("M", "F")),
    seed = 12345
  )

  # Process the first object through RandCapProd
  rand_prod_1 <- RandCapProd(rand_gen_1, seed = 54321)

  # Repeat the process: Generate the randomization object again
  rand_gen_2 <- RandCapGen(
    n = 20,
    id_prefix = "ID_",
    block_prefix = "Block_",
    block_sizes = c(10, 10),
    arms = c("Control", "Treatment"),
    strat_vars = list(Sex = c("M", "F")),
    seed = 12345
  )

  # Process the second object through RandCapProd
  rand_prod_2 <- RandCapProd(rand_gen_2, seed = 54321)

  # Verify that the resulting RandCapProd objects are identical
  expect_equal(
    rand_prod_1,
    rand_prod_2,
    info = "The RandCapProd output is not reproducible across repeated processes."
  )
})

