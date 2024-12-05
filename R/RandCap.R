## RANDCAPGEN---------
#' Generate Randomized Allocation List with Blocking and Stratification
#'
#' `RandCapGen` generates a randomization list that supports blocked and stratified randomization with customizable allocation ratios.
#' This function is ideal for designing randomized allocation schemes in clinical trials or research studies that require balanced participant allocation.
#'
#' @param n Integer. The total number of participants to be randomized.
#' @param id_prefix Character string. A prefix for participant IDs in the randomization list.
#' @param block_prefix Character string. A prefix for each block in the randomization list.
#' @param block_sizes Integer vector. Specifies the sizes of each block for blocked randomization.
#' @param arms Character vector. The names of treatment arms/groups.
#' @param ratio Integer vector. The allocation ratio for each treatment arm, e.g., `c(1, 2)` for a 1:2 ratio.
#' @param strat_vars Character vector. Names of stratification variables to ensure balanced randomization across subgroups.
#' @param strat_vars_prefix Character string. Prefix to add to stratification variable names in the output.
#' @param seed Integer. Seed for random number generation to ensure reproducibility.
#' @param uneq_beg Integer. The beginning of unequal allocation in blocks, used when blocks differ in allocation ratios.
#' @param uneq_mid Integer. The midpoint for unequal allocation within blocks.
#' @param uneq_min Integer. Minimum size of blocks when using unequal allocation.
#' @param uneq_maxit Integer. Maximum iterations to achieve the target allocation ratio.
#' @param project_acronym Character string. Acronym for the project, used in file names for easy identification.
#'
#' @return A list containing two data frames:
#' \describe{
#'   \item{randomization_list}{A complete randomization table with participant IDs, blocks, treatment arms, and stratification variables.}
#'   \item{redcap_import}{A simplified randomization list ready for import into REDCap.}
#' }
#'
#' @seealso
#' \code{\link{RandCapBalance}} for assessing randomization balance,
#' \code{\link{RandCapSettings}} for exporting randomization settings,
#' \code{\link{RandCapTable}},
#' and \code{\link{RandCapProd}} for exporting randomization tables for REDCap.
#'
#' @export
#'
#' @examples
#' # Basic example with 100 participants, two treatment arms, and a 1:1 allocation ratio
#' RandCapGen(
#'   n = 300,
#'   id_prefix = "ID_",
#'   block_prefix = "Block_",
#'   block_sizes = c(2,4,6),
#'   arms = c("A", "B"),
#'   ratio = c(1, 2),
#'   strat_vars = c(sex=c(1,2)),
#'   seed = 1501
#' )
#' @importFrom dplyr count %>%


RandCapGen <- function(n, id_prefix = NULL, block_prefix = NULL, block_sizes,
                       arms, ratio = NULL, strat_vars = NULL,
                       strat_vars_prefix = NULL, seed = NULL,
                       uneq_beg = FALSE, uneq_mid = FALSE, uneq_min = 1,
                       uneq_maxit = 1000, project_acronym = NULL) {
  # Set the initial seed for generating seeds for each stratum
  set.seed(seed)

  # Check that the number of arms is valid
  if (length(arms) < 2) {
    stop("Error: At least two arms are required for randomization.", call. = FALSE)
  }

  # Default ratio to 1 for each arm if not provided
  if (is.null(ratio)) {
    ratio <- rep(1, length(arms))
  }

  # Check for ratio consistency
  if (length(ratio) != length(arms)) {
    stop("Error: The length of the 'ratio' vector must match the number of arms.", call. = FALSE)
  }

  if (n %% sum(ratio) != 0) {
    stop("Error: 'n' must be a multiple of the sum of the 'ratio' vector.", call. = FALSE)
  }

  if (any(block_sizes %% sum(ratio) != 0)) {
    stop("Error: Each block size must be a multiple of the sum of the 'ratio' vector.", call. = FALSE)
  }

  # Automatically determine the number of levels based on the length of 'arms'
  num_levels <- length(arms)

  # Generate balanced levels if a ratio is provided
  if (!is.null(ratio)) {
    # Function to generate levels according to the specified ratio
    generate_levels_with_ratio <- function(arms, ratio) {
      # Check if the 'arms' vector contains duplicates
      if (any(duplicated(arms))) {
        stop("Error: The 'arms' vector contains duplicated values. Please ensure all arms are unique.", call. = FALSE)
      }
      # Ensure that 'arms' and 'ratio' are of the same length
      if (length(arms) != length(ratio)) {
        stop("Error: The length of 'arms' and 'ratio' must be the same.", call. = FALSE)
      }
      # Generate levels based on the ratio
      unlist(mapply(rep, arms, ratio))
    }
    # Use the generated levels with the specified ratio
    balanced_levels <- generate_levels_with_ratio(arms, ratio)
  } else {
    # If no ratio is provided, use the original arms
    balanced_levels <- arms
  }

  # Calculate block sizes by dividing each element in block_sizes by the number of balanced levels
  adjusted_block_sizes <- block_sizes / length(balanced_levels)

  # Define the function to shuffle blocks
  shuffle_blocks <- function(data) {
    # Split the data into blocks based on block.id
    blocks <- split(data, data$block.id)
    # Shuffle the blocks
    shuffled_blocks <- sample(blocks)
    # Recombine the shuffled blocks into a final data frame
    shuffled_data <- do.call(rbind, shuffled_blocks)
    return(shuffled_data)
  }

  # Check if stratification variables are provided
  if (is.null(strat_vars)) {
    # No stratification: Perform a single global block randomization
    rand_df <- blockrand::blockrand(
      n = n,
      num.levels = length(balanced_levels),
      id.prefix = id_prefix,
      block.prefix = block_prefix,
      block.sizes = adjusted_block_sizes,
      levels = balanced_levels,
      uneq.beg = uneq_beg,
      uneq.mid = uneq_mid,
      uneq.min = uneq_min,
      uneq.maxit = uneq_maxit
    )

    # Rename the treatment column to "treatment_arm"
    colnames(rand_df)[colnames(rand_df) == "treatment"] <- "treatment_arm"

    # Shuffle the blocks
    rand_df <- shuffle_blocks(rand_df)

    # Full dataset with only the relevant columns
    final_df <- rand_df[, c("id", "block.id", "block.size", "treatment_arm")]

    # Simplified dataset with only "treatment_arm"
    simplified_df <- final_df[, "treatment_arm", drop = FALSE]
  } else {

    # Stratified randomization with combinations of stratification variables

    # Generate all combinations of stratification variables
    strat_combinations <- expand.grid(strat_vars)

    # Generate a unique seed for each stratum
    stratum_seeds <- sample(1:10000, nrow(strat_combinations), replace = FALSE)

    # Initialize an empty list to store randomization data frames
    all_strata <- list()

    # Iterate through each combination of stratification variables
    for (i in 1:nrow(strat_combinations)) {
      # Set the seed for this specific stratum
      set.seed(stratum_seeds[i])

      # Get the current combination of stratification levels
      strat_comb <- strat_combinations[i, , drop = FALSE]

      # Create a custom label for the stratum based on the combination of stratification variables
      if (ncol(strat_combinations) == 1) {
        stratum_label <- paste0(strat_vars_prefix[[names(strat_comb)]], strat_comb[[1]])
      } else {
        stratum_label <- paste(
          mapply(function(var_name, value) {
            paste0(strat_vars_prefix[[var_name]], value)
          }, names(strat_comb), strat_comb),
          collapse = "|"
        )
      }

      # Define custom prefixes for block and ID
      custom_block_prefix <- paste0("S", i, "_", block_prefix)
      custom_id_prefix <- paste0("S", i, "_", id_prefix)

      # Perform block randomization for this specific stratum
      rand_df <- blockrand::blockrand(
        n = n,
        num.levels = length(balanced_levels),
        id.prefix = custom_id_prefix,
        block.prefix = custom_block_prefix,
        stratum = stratum_label,
        block.sizes = adjusted_block_sizes,
        levels = balanced_levels,
        uneq.beg = uneq_beg,
        uneq.mid = uneq_mid,
        uneq.min = uneq_min,
        uneq.maxit = uneq_maxit
      )

      # Rename the treatment column to "treatment_arm"
      colnames(rand_df)[colnames(rand_df) == "treatment"] <- "treatment_arm"

      # Add stratification variables to the data frame
      for (strat_var_name in names(strat_comb)) {
        new_var_name <- if (tolower(strat_var_name) %in% c("center", "centre")) {
          "redcap_data_access_group"
        } else {
          paste0("strat_", strat_var_name)
        }
        rand_df[[new_var_name]] <- strat_comb[[strat_var_name]]
      }

      # Shuffle the blocks within the randomized data frame
      rand_df_shuffled <- shuffle_blocks(rand_df)

      # Append the shuffled randomization data frame for this stratum to the list
      all_strata[[stratum_label]] <- rand_df_shuffled
    }

    # Combine all strata into one final data frame
    final_df <- dplyr::bind_rows(all_strata)

    # Create the simplified dataset with only 'treatment_arm' and stratification variables
    stratification_columns <- grep("strat_|redcap_data_access_group", colnames(final_df), value = TRUE)
    simplified_df <- final_df[, c("treatment_arm", stratification_columns), drop = FALSE]

  }

  # Convert all columns in simplified_df to factors
  simplified_df <- dplyr::mutate(simplified_df, dplyr::across(everything(), as.factor))

  rownames(final_df) <- NULL
  rownames(simplified_df) <- NULL

  # Create the output structure
  randomization_object <- list(
    settings = list(
      R_version = R.version.string,
      date = Sys.time(),
      initial_seed = seed,
      treatment_arms = arms,
      ratios = ratio,
      strat_vars = if (!is.null(strat_vars)) names(strat_vars) else NULL,
      strat_vars_prefix = strat_vars_prefix,
      strat_levels = strat_vars,
      block_sizes_input = block_sizes,
      type_of_list = "development",
      project_acronym = project_acronym
    ),
    tables = list(
      full_dataset = final_df,
      simplified_dataset = simplified_df
    )
  )

  # Return the structured object
  return(randomization_object)
}

## RANDCAPPROD---------
#' Generate Production Randomization
#'
#' This function finalizes the production randomization based on the development settings
#' in the `randomization_object`. It shuffles the blocks and generates production-ready randomization tables.
#'
#' @param randomization_object A list containing randomization settings and tables, created by `RandCapGen`.
#' @param seed An integer specifying the seed for reproducibility.
#' @return A modified `randomization_object` with updated tables and settings for production.
#' @export
#' @examples
#' # Assuming you have a randomization object from RandCapGen
#' RandCapProd(randomization_object = my_randomization_object, seed = 123)
#' @seealso \code{\link{RandCapGen}},
#' \code{\link{RandCapSettings}},
#' \code{\link{RandCapTable}},
#' and \code{\link{RandCapBalance}}


RandCapProd <- function(randomization_object, seed) {

  shuffle_blocks <- function(data, seed = NULL) {
    # Set the seed for reproducibility if provided
    if (!is.null(seed)) set.seed(seed)

    # Split the dataset into blocks based on block.id
    blocks <- split(data, data$block.id)

    # Shuffle the rows within each block
    shuffled_blocks <- lapply(blocks, function(block) block[sample(nrow(block)), ])

    # Shuffle the order of the blocks themselves
    shuffled_blocks <- sample(shuffled_blocks)

    # Recombine the shuffled blocks into a final dataset
    shuffled_data <- do.call(rbind, shuffled_blocks)

    return(shuffled_data)
  }

  full_dataset <- randomization_object$tables$full_dataset
  full_dataset_shuffled <- shuffle_blocks(full_dataset, seed = seed)
  rownames(full_dataset_shuffled) <- NULL

  # Check if the 'stratum' column exists before trying to select
  columns_to_remove <- c("id", "block.id", "block.size")
  if ("stratum" %in% colnames(full_dataset_shuffled)) {
    columns_to_remove <- c(columns_to_remove, "stratum")
  }

  simplified_dataset_shuffled <- full_dataset_shuffled %>%
    dplyr::select(-all_of(columns_to_remove))

  randomization_object$tables$full_dataset <- full_dataset_shuffled
  randomization_object$tables$simplified_dataset <- simplified_dataset_shuffled
  randomization_object$settings$initial_seed <- seed
  randomization_object$settings$type_of_list <- "production"

  return(randomization_object)
}



#RANDCAP BALANCE--------
#' Assess Randomization Balance
#'
#' `RandCapBalance` checks the balance across stratified groups after randomization, based on specified stratification variables.
#'
#' @param randomization_object A list. The output object from `RandCapGen`, containing randomization settings and tables.
#' @param output_path Character string. The file path to save the balance summary as a PDF.
#'
#' @return A data frame showing the balance across stratified groups for each treatment arm.
#'
#' @export
#'
#' @examples
#' # Assuming you have a randomization object from RandCapGen
#' RandCapBalance(randomization_object = my_randomization_object, output_path = "Balance_Summary.pdf")
#' @seealso \code{\link{RandCapGen}},
#' \code{\link{RandCapSettings}},
#' \code{\link{RandCapTable}},
#' and \code{\link{RandCapProd}}

#' @importFrom tidyr pivot_wider
#' @importFrom gridExtra tableGrob ttheme_minimal
#' @importFrom grid grid.newpage grid.layout grid.text viewport pushViewport popViewport gpar unit

RandCapBalance <- function(randomization_object,
                           output_path = "Randomization_Balance.pdf") {

  if(randomization_object$settings$type_of_list=="development"){
    output_path = "Randomization_Balance_dev.pdf"
  }else{output_path = "Randomization_Balance_prod.pdf"}

  # Extract the full dataset and calculate the number of strata
  full_dataset <- randomization_object$tables$full_dataset
  num_strata <- length(unique(full_dataset$stratum))



  # Handle the case where no strata are defined
  if (num_strata == 0) {
    full_dataset <- full_dataset %>% dplyr::mutate(stratum = "No stratum")
    colnames(full_dataset)[4] <- "treatment_arm"
  }else{
    colnames(full_dataset)[5] <- "treatment_arm"
  }

  # Create the balance table with dplyr and tidyr
  balance_df <- full_dataset %>%
    dplyr::count(stratum, treatment_arm) %>%
    tidyr::pivot_wider(names_from = treatment_arm, values_from = n, values_fill = 0) %>%
    dplyr::mutate(Total = rowSums(dplyr::select(., -stratum)))

  # Get randomization settings
  strat_vars <- randomization_object$settings$strat_vars
  strat_vars_text <- if (is.null(strat_vars) || length(strat_vars) == 0) {
    "There is no stratification variable."
  } else {
    paste(strat_vars, collapse = ", ")
  }

  ratios <- if (is.null(randomization_object$settings$ratio)) {
    "Equal size ratio"
  } else {
    paste(randomization_object$settings$ratio, collapse = ", ")
  }

  treatment_arms <- paste(randomization_object$settings$treatment_arms, collapse = ", ")

  # Start the PDF
  pdf(output_path, width = 8.5, height = 11)

  # Determine the layout for the first page based on the number of rows in the table
  grid::grid.newpage()

  if (nrow(balance_df) >= 20) {
    # Case with at least 20 rows in the table
    text_layout_height <- 0.3
    table_layout_height <- 0.7
    table_rows <- 20
  } else {
    # Case with fewer than 20 rows in the table
    text_layout_height <- 0.3
    table_layout_height <- min(1, 0.025 * nrow(balance_df) + 0.1) # 2.5% * nrow + offset for spacing
    remaining_height <- 1 - (text_layout_height + table_layout_height)
    table_rows <- nrow(balance_df)
  }

  # Display text in the upper section (30% of the page)
  grid::pushViewport(grid::viewport(layout = grid.layout(10, 1)))  # Divide into 10 for proportional layout
  grid::pushViewport(grid::viewport(layout.pos.row = 1:3, layout.pos.col = 1))  # Upper section for text

  # Display headers and values with blue and bold style for headers and black for values
  if(randomization_object$settings$type_of_list=="development"){
    grid::grid.text("Randomization Balance Summary for development list",
              x = 0.15, y = 0.8, just = "left",
              gp = gpar(fontsize = 18, fontface = "bold", col = "#2E86C1"))
  }else{
    grid::grid.text("Randomization Balance Summary for production list", x = 0.15,
              y = 0.8, just = "left",
              gp = gpar(fontsize = 18, fontface = "bold", col = "#2E86C1"))
  }

  # Display number of strata
  if (num_strata == 0) {
    grid::grid.text("Number of Strata:", x = 0.15, y = 0.55, just = "left",
              gp = gpar(fontsize = 12, fontface = "bold", col = "#2874A6"))
    grid::grid.text("The randomization is not stratified", x = 0.4, y = 0.55, just = "left",
              gp = gpar(fontsize = 12, col = "#34495E"))
  } else {
    grid::grid.text("Number of Strata:", x = 0.15, y = 0.55, just = "left",
              gp = gpar(fontsize = 12, fontface = "bold", col = "#2874A6"))
    grid::grid.text(num_strata, x = 0.4, y = 0.55, just = "left",
              gp = gpar(fontsize = 12, col = "#34495E"))
  }

  grid::grid.text("Stratification variables:", x = 0.15, y = 0.4, just = "left",
            gp = gpar(fontsize = 12, fontface = "bold", col = "#2874A6"))
  grid::grid.text(strat_vars_text, x = 0.4, y = 0.4, just = "left",
            gp = gpar(fontsize = 12, col = "#34495E"))

  grid::grid.text("Treatment Arms:", x = 0.15, y = 0.25, just = "left",
            gp = gpar(fontsize = 12, fontface = "bold", col = "#2874A6"))
  grid::grid.text(treatment_arms, x = 0.4, y = 0.25, just = "left",
            gp = gpar(fontsize = 12, col = "#34495E"))

  grid::grid.text("Ratios:", x = 0.15, y = 0.1, just = "left",
            gp = gpar(fontsize = 12, fontface = "bold", col = "#2874A6"))
  grid::grid.text(ratios, x = 0.4, y = 0.1, just = "left",
            gp = gpar(fontsize = 12, col = "#34495E"))

  grid::popViewport()

  # Display the table in the lower section
  table_grob1 <- gridExtra::tableGrob(balance_df[1:table_rows, ], rows = NULL)
  table_theme_page1 <- ttheme_minimal(
    core = list(
      fg_params = list(col = "black"),
      bg_params = list(
        fill = c(
          rep("#D1E8FF", table_rows),                             # Color for strata column
          rep("#EBF5FB", table_rows * (ncol(balance_df) - 2)),    # Color for treatment columns
          rep("#D1E8FF", table_rows)                              # Color for Total column
        ),
        col = NA
      )
    ),
    colhead = list(
      fg_params = list(col = "white", fontface = "bold"),
      bg_params = list(fill = "#2E86C1")
    )
  )

  grid::pushViewport(grid::viewport(y = 1 - text_layout_height,
                        height = unit(table_layout_height, "npc"),
                        just = "top"))
  grid::grid.draw(gridExtra::tableGrob(balance_df[1:table_rows, ], rows = NULL,
                      theme = table_theme_page1))
  grid::popViewport()

  # Add an empty layout below if there are fewer than 20 rows in the table
  if (nrow(balance_df) < 20) {
    grid::pushViewport(grid::viewport(y = text_layout_height, height =
                            unit(remaining_height, "npc"), just = "top"))
    grid::popViewport()
  }

  # Pages for remaining rows
  start_row <- table_rows + 1
  num_rows <- nrow(balance_df)
  rows_per_page <- 31

  while (start_row <= num_rows) {
    grid::grid.newpage()

    # Check if this is the last page
    is_last_page <- (start_row + rows_per_page - 1 >= num_rows)
    end_row <- if (is_last_page) num_rows else (start_row + rows_per_page - 1)
    table_page <- balance_df[start_row:end_row, ]

    # Colors for each page
    page_rows <- nrow(table_page)
    page_fill_colors <- c(
      rep("#D1E8FF", page_rows),                       # Color for strata column
      rep("#EBF5FB", page_rows * (ncol(balance_df) - 2)), # Color for treatment columns
      rep("#D1E8FF", page_rows)                        # Color for Total column
    )

    table_theme <- ttheme_minimal(
      core = list(
        fg_params = list(col = "black"),
        bg_params = list(fill = page_fill_colors, col = NA)
      ),
      colhead = list(
        fg_params = list(col = "white", fontface = "bold"),
        bg_params = list(fill = "#2E86C1")
      )
    )
    table_grob <- gridExtra::tableGrob(table_page, rows = NULL, theme = table_theme)

    if (is_last_page) {
      # Last page: Adjust the height based on remaining rows with a factor of 3%
      remaining_rows <- end_row - start_row + 1
      last_page_height <- min(0.9, 0.025 * remaining_rows)

      # Position table in the upper section of the last page
      grid::pushViewport(grid::viewport(y = 0.9, height = unit(last_page_height, "npc"), just = "top"))
    } else {
      # Intermediate pages with table at the top
      grid::pushViewport(grid::viewport(y = 0.87, height = unit(0.75, "npc"), just = "top"))
    }

    grid::grid.draw(table_grob)
    grid::popViewport()

    # Move to the next set of rows
    start_row <- end_row + 1
  }

  # Close the PDF
  dev.off()

  message("Randomization balance summary saved to ", output_path)
}


# RANDCAP SETTINGS FUNCTION-------
#' Export Randomization Settings
#'
#' `RandCapSettings` generates a PDF summary of the randomization settings used in a randomization object, including details about stratification, treatment arms, and ratios.
#'
#' @param rand_obj A list. The randomization object generated by `RandCapGen`.
#' @param output_path Character string. The file path to save the settings PDF.
#' @param display_block_sizes Logical. Whether to display block sizes in the output PDF.
#'
#' @return None. Outputs a PDF file summarizing randomization settings.
#'
#' @export
#'
#' @examples
#' # Example usage
#' RandCapSettings(rand_obj = my_randomization_object, output_path = "Randomization_Settings.pdf")
#' @seealso \code{\link{RandCapGen}}, \code{\link{RandCapBalance}}, \code{\link{RandCapTable}}


RandCapSettings <- function(rand_obj, output_path = "Randomization_settings.pdf",
                            display_block_sizes = FALSE) {

  if(rand_obj$settings$type_of_list=="development"){
    output_path = "Randomization_settings_dev.pdf"
  }else{ output_path = "Randomization_settings_prod.pdf"}

  # Start PDF file
  pdf(output_path, width = 8.5, height = 11, family = "serif")

  # Title
  grid::grid.newpage()
  if(rand_obj$settings$type_of_list=="development"){
    grid::grid.text("Randomization settings for development list", x = 0.1, y = 0.9, just = "left",
            gp = gpar(fontsize = 20, fontface = "bold", col = "#2E86C1"))
  }else{
    grid::grid.text("Randomization settings for production list", x = 0.1, y = 0.9, just = "left",
              gp = gpar(fontsize = 20, fontface = "bold", col = "#2E86C1"))
    }

  # Variables
  r_version <- rand_obj$settings$R_version
  date_randomization <- rand_obj$settings$date
  seed <- if (is.null(rand_obj$settings$initial_seed)
  ) "No seed was specified" else rand_obj$settings$initial_seed
  treatment_arms <- paste(rand_obj$settings$treatment_arms, collapse = ", ")
  ratios <- if (is.null(rand_obj$settings$ratio)
  ) "Equal size ratio" else paste(rand_obj$settings$ratio,
                                  collapse = ", ")
  strat_vars <- rand_obj$settings$strat_vars
  strat_levels <- rand_obj$settings$strat_levels
  block_sizes <- paste(rand_obj$settings$block_sizes_input, collapse = ", ")

  # Display function with alignment adjustments
  display_setting <- function(label, value, y_pos) {
    # Bold blue label
    grid::grid.text(label, x = 0.1, y = y_pos, just = "left",
              gp = gpar(fontsize = 12, fontface = "bold", col = "#2874A6"))
    # Black value with adjusted spacing
    grid::grid.text(value, x = 0.45, y = y_pos, just = "left",
              gp = gpar(fontsize = 12, col = "#34495E"))
  }

  # Display settings with alignment adjustments
  y_pos <- 0.85
  display_setting("R Version:", r_version, y_pos); y_pos <- y_pos - 0.03
  display_setting("Date of Randomization:", date_randomization,
                  y_pos); y_pos <- y_pos - 0.03
  display_setting("Seed:", seed, y_pos); y_pos <- y_pos - 0.03
  display_setting("Treatment Arms:", treatment_arms, y_pos);
  y_pos <- y_pos - 0.03
  display_setting("Ratios:", ratios, y_pos); y_pos <- y_pos - 0.03

  # Stratification variables and levels with line wrapping for Centers
  if (is.null(strat_vars) || length(strat_vars) == 0) {
    # No stratification variable
    grid::grid.text("Stratification:", x = 0.1, y = y_pos, just = "left",
              gp = gpar(fontsize = 12, fontface = "bold", col = "#2874A6"))
    grid::grid.text("There is no stratification variable.", x = 0.45,
              y = y_pos, just = "left",
              gp = gpar(fontsize = 12, col = "#34495E"))
  } else {
    # Display stratification variables
    grid::grid.text("Stratification Variables:", x = 0.1, y = y_pos, just = "left",
              gp = gpar(fontsize = 12, fontface = "bold", col = "#2874A6"))
    y_pos <- y_pos - 0.03

    for (i in seq_along(strat_vars)) {
      # Check for Center-related variables
      if (tolower(strat_vars[i]) %in% c("center", "centre",
                                        "redcap_data_access_group")) {
        # Wrap levels for "Centers"
        levels_text <- paste(strat_levels[[i]], collapse = ", ")
        wrapped_levels <- strwrap(levels_text, width = 80)

        # Print "Centers:" as the label and wrap the levels across lines
        grid::grid.text("- Centers:", x = 0.15, y = y_pos, just = "left",
                  gp = gpar(fontsize = 12))
        y_pos <- y_pos - 0.03

        for (line in wrapped_levels) {
          grid::grid.text(line, x = 0.2, y = y_pos, just = "left",
                    gp = gpar(fontsize = 12, col = "#34495E"))
          y_pos <- y_pos - 0.03
        }

      } else {
        # For other stratification variables, display without line wrapping
        grid::grid.text(paste0("− ", strat_vars[i], ": ", paste(strat_levels[[i]],
                                                          collapse = ", ")),
                  x = 0.15, y = y_pos, just = "left",
                  gp = gpar(fontsize = 12, col = "#34495E"))
        y_pos <- y_pos - 0.025
      }
    }
  }

  # Conditional display of block sizes
  y_pos <- y_pos - 0.03
  grid::grid.text("Block Size per Stratum:", x = 0.1, y = y_pos, just = "left",
            gp = gpar(fontsize = 12, fontface = "bold", col = "#2874A6"))

  if (display_block_sizes) {
    grid::grid.text(block_sizes, x = 0.45, y = y_pos, just = "left",
              gp = gpar(fontsize = 12, col = "#34495E"))
  } else {
    grid::grid.text("Not allowed to be displayed", x = 0.45, y = y_pos, just = "left",
              gp = gpar(fontsize = 12, col = "#E74C3C"))
  }

  # Close the PDF
  dev.off()

  message("Randomization settings saved to ", output_path)
}



#RANDCAP TABLE--------
#' Export Randomization Tables for REDCap
#'
#' `RandCapTable` saves the randomization tables from a `RandCapGen` object to CSV files, one for a complete randomization list and another simplified for REDCap import.
#'
#' @param randomization_object A list. The randomization object generated by `RandCapGen`.
#' @param save_for_REDCap Logical. If `TRUE`, saves the simplified table for REDCap as a CSV file.
#' @param save_random_table Logical. If `TRUE`, saves the complete randomization table as a CSV file.
#'
#' @return A list containing two data frames:
#' \describe{
#'   \item{simplified_table}{A simplified dataset ready for REDCap import.}
#'   \item{full_table}{The complete randomization dataset.}
#' }
#'
#' @export
#'
#' @examples
#' # Example usage to save both tables
#' RandCapTable(randomization_object = my_randomization_object, save_for_REDCap = TRUE, save_random_table = TRUE)
#' @seealso \code{\link{RandCapGen}},
#' \code{\link{RandCapBalance}},
#' \code{\link{RandCapSettings}},
#' and \code{\link{RandCapProd}}

RandCapTable <- function(randomization_object, save_for_REDCap = FALSE,
                         save_random_table = FALSE) {
  # Save the simplified table if save_for_REDCap is TRUE
  if (save_for_REDCap) {
    redcap_table_name <-if(
      randomization_object$settings$type_of_list=="development"){
      paste0(randomization_object$settings$project_acronym,
             "_REDCap_dev_table", "_", Sys.Date(), ".csv")
    }else{paste0(randomization_object$settings$project_acronym,
                 "_REDCap_prod_table", "_", Sys.Date(), ".csv")}
    write.csv(randomization_object$tables$simplified_dataset,
              redcap_table_name, row.names = FALSE)
  }

  # Save the full table if save_random_table is TRUE
  if (save_random_table) {
    randomization_table_name <-if(
      randomization_object$settings$type_of_list=="development"){
      paste0(
        randomization_object$settings$project_acronym,
        "_randomization_dev_table", "_", Sys.Date(), ".csv")
    }else{
      paste0(
        randomization_object$settings$project_acronym,
        "_randomization_prod_table", "_", Sys.Date(), ".csv")
    }
    write.csv(randomization_object$tables$full_dataset,
              randomization_table_name, row.names = FALSE)
  }

  # Return both tables in a list
  # return(list(
  #   simplified_table = randomization_object$tables$simplified_dataset,
  #   full_table = randomization_object$tables$full_dataset
  # ))
}




#
# Hello World
# Follow this Hello World exercise to learn GitHub's pull request workflow.
#
# In this article
# Introduction
# Step 1: Create a repository
# Step 2: Create a branch
# Step 3: Make and commit changes
# Step 4: Open a pull request
# Step 5: Merge your pull request
# Conclusion
# Next steps
# Further reading
# Introduction
# This tutorial teaches you GitHub essentials like repositories, branches, commits, and pull requests. You'll create your own Hello World repository and learn GitHub's pull request workflow, a popular way to create and review code.
#
# In this quickstart guide, you will:
#
# Create and use a repository.
# Start and manage a new branch.
# Make changes to a file and push them to GitHub as commits.
# Open and merge a pull request.
# Prerequisites
# You must have a GitHub account. For more information, see "Creating an account on GitHub."
#
# You don't need to know how to code, use the command line, or install Git (the version control software that GitHub is built on).
#
# Step 1: Create a repository
# The first thing we'll do is create a repository. You can think of a repository as a folder that contains related items, such as files, images, videos, or even other folders. A repository usually groups together items that belong to the same "project" or thing you're working on.
#
# Often, repositories include a README file, a file with information about your project. README files are written in Markdown, which is an easy-to-read, easy-to-write language for formatting plain text. We'll learn more about Markdown in the next tutorial, "Setting up your profile."
#
# GitHub lets you add a README file at the same time you create your new repository. GitHub also offers other common options such as a license file, but you do not have to select any of them now.
#
# Your hello-world repository can be a place where you store ideas, resources, or even share and discuss things with others.
#
# In the upper-right corner of any page, select , then click New repository.
#
# Screenshot of a GitHub dropdown menu showing options to create new items. The menu item "New repository" is outlined in dark orange.
# In the "Repository name" box, type hello-world.
#
# In the "Description" box, type a short description. For example, type "This repository is for practicing the GitHub Flow."
#
# Select whether your repository will be Public or Private.
#
# Select Add a README file.
#
# Click Create repository.
