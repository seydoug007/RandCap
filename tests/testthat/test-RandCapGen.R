test_that("RandCapGen outputs have correct structure", {
  # Case 1: Verify output structure
  result <- RandCapGen(
    n = 10,
    id_prefix = "ID_",
    block_prefix = "Block_",
    block_sizes = c(4, 6),
    arms = c("Control", "Treatment"),
    seed = 123
  )

  # Verify that the result is a list
  expect_type(result, "list")

  # Verify that the list contains the expected elements
  expect_named(result, c("settings", "tables"))

  # Verify the structure of datasets in "tables"
  expect_named(result$tables, c("full_dataset", "simplified_dataset"))

  # Verify columns in the full dataset
  expect_true(all(c("id", "block.id", "block.size", "treatment_arm") %in% colnames(result$tables$full_dataset)))
})

test_that("RandCapGen respects block sizes", {
  # Case 2: Verify that block sizes are respected
  result <- RandCapGen(
    n = 12,
    id_prefix = "ID_",
    block_prefix = "Block_",
    block_sizes = c(2, 4, 6),
    arms = c("Control", "Treatment"),
    seed = 456
  )

  # Verify that generated blocks have the expected sizes
  block_sizes <- result$tables$full_dataset %>%
    dplyr::count(block.id) %>%
    dplyr::pull(n)

  expect_equal(block_sizes, c(2, 4, 6))
})

test_that("RandCapGen handles allocation ratios correctly", {
  # Case 3: Verify correct allocation ratios
  result <- RandCapGen(
    n = 12,
    id_prefix = "ID_",
    block_prefix = "Block_",
    block_sizes = c(6, 6),
    arms = c("Control", "Treatment1", "Treatment2"),
    ratio = c(1, 2, 3), # Ratio 1:2:3
    seed = 789
  )

  # Verify that ratios are respected in the full dataset
  treatment_counts <- result$tables$full_dataset %>%
    dplyr::count(treatment_arm) %>%
    dplyr::pull(n)

  expect_equal(treatment_counts, c(2, 4, 6)) # 1:2:3 ratio for 12 participants
})

test_that("RandCapGen works without stratification", {
  # Case 4: Generation without stratification
  result <- RandCapGen(
    n = 8,
    id_prefix = "ID_",
    block_prefix = "Block_",
    block_sizes = c(4, 4),
    arms = c("A", "B"),
    seed = 1010
  )

  # Verify that no stratification column is present
  expect_false("stratum" %in% colnames(result$tables$full_dataset))
})

test_that("RandCapGen handles one stratification variable correctly", {
  # Case 5: Verify stratification
  result <- RandCapGen(
    n = 12,
    id_prefix = "ID_",
    block_prefix = "Block_",
    block_sizes = c(6, 6),
    arms = c("Control", "Treatment"),
    strat_vars = list(Sex = c("M", "F")),
    seed = 2024
  )

  # Verify that stratification is correctly defined
  expect_true("strat_Sex" %in% colnames(result$tables$full_dataset))

  # Verify that stratification values are correct
  strat_values <- unique(result$tables$full_dataset$strat_Sex)
  expect_equal(sort(as.character(strat_values)), c("F", "M"))
})

test_that("RandCapGen correctly handles multiple stratification variables", {
  # Generate data with multiple stratification variables
  result <- RandCapGen(
    n = 16,
    id_prefix = "ID_",
    block_prefix = "Block_",
    block_sizes = c(8, 8),
    arms = c("Control", "Treatment"),
    strat_vars = list(
      Sex = c("M", "F"),
      AgeGroup = c("Young", "Old"),
      Comorbidity = c("No", "Yes"),
      Centre = c("A", "B", "C")
    ),
    seed = 12345
  )

  # Verify that stratification columns are present in the full dataset
  expect_true("strat_Sex" %in% colnames(result$tables$full_dataset))
  expect_true("strat_AgeGroup" %in% colnames(result$tables$full_dataset))
  expect_true("strat_Comorbidity" %in% colnames(result$tables$full_dataset))
  expect_true("redcap_data_access_group" %in% colnames(result$tables$full_dataset))

  # Verify that expected columns are present in the simplified dataset
  expected_columns <- c("strat_Sex", "strat_AgeGroup", "strat_Comorbidity", "redcap_data_access_group")
  actual_columns <- colnames(result$tables$simplified_dataset)
  expect_equal(
    actual_columns,
    expected_columns,
    info = "Columns are either missing or not in the correct order in the simplified dataset."
  )

  # Verify that treatment_arm is identical to redcap_data_access_group in the simplified dataset
  simplified_data <- result$tables$simplified_dataset
  expect_true(
    all(simplified_data$treatment_arm == simplified_data$redcap_data_access_group),
    info = "The treatment_arm column is not identical to redcap_data_access_group."
  )
})


test_that("RandCapGen ensures reproducibility", {
  # Define the parameters
  params <- list(
    n = 90,
    id_prefix = "ID_",
    block_prefix = "Block_",
    block_sizes = c(30, 30, 30),
    arms = c("Control", "Treatment"),
    ratio = c(2, 1),
    strat_vars = list(Sex = c("M", "F"),
                      AgeGroup=c("Young","Old")),
    seed = 12345
  )

  # Run the function twice with the same parameters
  result1 <- do.call(RandCapGen, params)
  result2 <- do.call(RandCapGen, params)


  # Test 1: Reproducibility of treatment sequence within strata
  treatment_sequence1 <- result1$tables$full_dataset %>%
    dplyr::group_by(strat_Sex,strat_AgeGroup) %>%
    dplyr::summarize(treatment_sequence = list(treatment_arm))

  treatment_sequence2 <- result2$tables$full_dataset %>%
    dplyr::group_by(strat_Sex,strat_AgeGroup) %>%
    dplyr::summarize(treatment_sequence = list(treatment_arm))

  expect_equal(
    treatment_sequence1$treatment_sequence,
    treatment_sequence2$treatment_sequence,
    info = "Treatment sequence within strata is not reproducible."
  )

  # Test 3: Reproducibility of patient IDs and block IDs
  expect_equal(
    result1$tables$full_dataset$id,
    result2$tables$full_dataset$id,
    info = "Patient IDs are not reproducible."
  )

  expect_equal(
    result1$tables$full_dataset$block.id,
    result2$tables$full_dataset$block.id,
    info = "Block IDs are not reproducible."
  )
})


