#' Encode DNA Sequences into Nucleotide Features
#'
#' Converts DNA sequences into a data frame where each nucleotide
#' position is represented as a categorical factor.
#'
#' @param dna_strings A character vector containing DNA sequences
#' of equal length, using A, T, C and G.
#'
#' @return A data frame in which each column represents one
#' nucleotide position and contains categorical factor values.
dna_encoding <- function(dna_strings){
  dna_strings <- as.character(dna_strings)

  if (!is.character(dna_strings) || length(dna_strings) == 0L ||
      anyNA(dna_strings) || any(!grepl("^[ATCG]+$", dna_strings))) {
    stop("dna_strings must contain non-missing DNA sequences using only A, T, C and G.")
  }

  sequence_lengths <- nchar(dna_strings)
  if (length(unique(sequence_lengths)) != 1L) {
    stop("All DNA sequences must have the same length.")
  }

  nn <- sequence_lengths[1]
  seq_m <- do.call(rbind, strsplit(dna_strings, "", fixed = TRUE))
  colnames(seq_m) <- paste0("nt_pos", seq_len(nn))
  seq_df <- as.data.frame(seq_m, stringsAsFactors = FALSE)
  seq_df[] <- lapply(seq_df, factor, levels = c("A", "T", "C", "G"))
  return(seq_df)
}

#' Predict m6A Sites in Multiple Sequences
#'
#' Uses a fitted random forest model to predict the probability
#' and status of m6A sites for multiple input sequences.
#'
#' @param ml_fit A fitted random forest classification model.
#' @param feature_df A data frame containing gc_content, RNA_type,
#' RNA_region, exon_length, distance_to_junction,
#' evolutionary_conservation and DNA_5mer.
#' The length and distance features should already be log2-transformed.
#' @param positive_threshold A numeric probability threshold.
#' Scores strictly greater than this threshold are classified
#' as Positive. Default is 0.5.
#'
#' @return A data frame containing the original input features,
#' predicted_m6A_prob and predicted_m6A_status.
#'
#' @import randomForest
#' @importFrom stats predict
#' @export
#' @examples
#' model_file <- system.file("extdata", "rf_fit.rds",
#'                           package = "m6APrediction")
#' ml_fit <- readRDS(model_file)
#'
#' data_file <- system.file("extdata", "m6A_input_example.csv",
#'                          package = "m6APrediction")
#' example_data <- read.csv(data_file)
#'
#' result <- prediction_multiple(
#'   ml_fit = ml_fit,
#'   feature_df = example_data
#' )
#' head(result)
prediction_multiple <- function(ml_fit, feature_df, positive_threshold = 0.5){
  required_columns <- c("gc_content", "RNA_type", "RNA_region", "exon_length",
                        "distance_to_junction", "evolutionary_conservation", "DNA_5mer")
  allowed_rna_types <- c("mRNA", "lincRNA", "lncRNA", "pseudogene")
  allowed_rna_regions <- c("CDS", "intron", "3'UTR", "5'UTR")
  numeric_columns <- c("gc_content", "exon_length", "distance_to_junction",
                       "evolutionary_conservation")

  if (!is.data.frame(feature_df) || nrow(feature_df) == 0L) {
    stop("feature_df must be a non-empty data frame.")
  }
  missing_columns <- setdiff(required_columns, names(feature_df))
  if (length(missing_columns) > 0L) {
    stop("feature_df is missing required columns: ",
         paste(missing_columns, collapse = ", "))
  }
  if (!is.numeric(positive_threshold) || length(positive_threshold) != 1L ||
      is.na(positive_threshold) || !is.finite(positive_threshold) ||
      positive_threshold < 0 || positive_threshold > 1) {
    stop("positive_threshold must be one finite number between 0 and 1.")
  }

  for (column in numeric_columns) {
    values <- feature_df[[column]]
    if (!is.numeric(values) || anyNA(values) || any(!is.finite(values))) {
      stop(column, " must contain only finite numeric values without missing data.")
    }
  }
  if (any(feature_df$gc_content < 0 | feature_df$gc_content > 1)) {
    stop("gc_content values must be between 0 and 1.")
  }
  if (any(feature_df$evolutionary_conservation < 0 |
          feature_df$evolutionary_conservation > 1)) {
    stop("evolutionary_conservation values must be between 0 and 1.")
  }
  if (any(feature_df$exon_length < 0) || any(feature_df$distance_to_junction < 0)) {
    stop("exon_length and distance_to_junction must be non-negative.")
  }

  rna_type <- as.character(feature_df$RNA_type)
  rna_region <- as.character(feature_df$RNA_region)
  if (anyNA(rna_type) || any(!rna_type %in% allowed_rna_types)) {
    stop("RNA_type must be one of: ", paste(allowed_rna_types, collapse = ", "), ".")
  }
  if (anyNA(rna_region) || any(!rna_region %in% allowed_rna_regions)) {
    stop("RNA_region must be one of: ", paste(allowed_rna_regions, collapse = ", "), ".")
  }

  dna_5mer <- as.character(feature_df$DNA_5mer)
  if (anyNA(dna_5mer) || any(!grepl("^[ATCG]{5}$", dna_5mer)) ||
      any(substr(dna_5mer, 3, 3) != "A")) {
    stop("DNA_5mer must contain exactly five A, T, C or G bases, with A at the central position.")
  }

  feature_df$RNA_type <- factor(rna_type, levels = allowed_rna_types)
  feature_df$RNA_region <- factor(rna_region, levels = allowed_rna_regions)
  infered_prob <- predict(
    ml_fit,
    newdata = cbind(feature_df, dna_encoding(dna_5mer)),
    type = "prob"
  )
  if (!"Positive" %in% colnames(infered_prob)) {
    stop("The fitted model must provide a 'Positive' probability column.")
  }

  feature_df$predicted_m6A_prob <- as.numeric(infered_prob[, "Positive"])
  feature_df$predicted_m6A_status <- ifelse(
    feature_df$predicted_m6A_prob > positive_threshold,
    "Positive", "Negative"
  )
  return(feature_df)
}

#' Predict an Individual m6A Site
#'
#' Uses a fitted random forest model to predict the probability
#' and status of m6A for a single input sequence.
#'
#' @param ml_fit A fitted random forest classification model.
#' @param gc_content The GC content of the sequence.
#' @param RNA_type The RNA type, such as mRNA or lncRNA.
#' @param RNA_region The RNA region, such as CDS or 3'UTR.
#' @param exon_length The log2-transformed exon length.
#' @param distance_to_junction The log2-transformed distance
#' to the nearest junction.
#' @param evolutionary_conservation The evolutionary conservation score.
#' @param DNA_5mer A five-nucleotide DNA sequence with A
#' at the central position.
#' @param positive_threshold A numeric probability threshold.
#' Scores strictly greater than this threshold are classified
#' as Positive. Default is 0.5.
#'
#' @return A named list containing a numeric
#' `predicted_m6A_prob` and a character
#' `predicted_m6A_status`.
#'
#' @export
#' @examples
#' model_file <- system.file("extdata", "rf_fit.rds",
#'                           package = "m6APrediction")
#' ml_fit <- readRDS(model_file)
#'
#' data_file <- system.file("extdata", "m6A_input_example.csv",
#'                          package = "m6APrediction")
#' example_data <- read.csv(data_file)
#'
#' single_result <- prediction_single(
#'   ml_fit = ml_fit,
#'   gc_content = example_data$gc_content[1],
#'   RNA_type = as.character(example_data$RNA_type[1]),
#'   RNA_region = as.character(example_data$RNA_region[1]),
#'   exon_length = example_data$exon_length[1],
#'   distance_to_junction = example_data$distance_to_junction[1],
#'   evolutionary_conservation = example_data$evolutionary_conservation[1],
#'   DNA_5mer = as.character(example_data$DNA_5mer[1])
#' )
#' print(single_result)
prediction_single <- function(ml_fit, gc_content, RNA_type, RNA_region, exon_length,
                              distance_to_junction, evolutionary_conservation, DNA_5mer, positive_threshold = 0.5){
  single_values <- list(
    gc_content = gc_content,
    RNA_type = RNA_type,
    RNA_region = RNA_region,
    exon_length = exon_length,
    distance_to_junction = distance_to_junction,
    evolutionary_conservation = evolutionary_conservation,
    DNA_5mer = DNA_5mer
  )
  if (any(lengths(single_values) != 1L)) {
    stop("Each input to prediction_single must contain exactly one value.")
  }

  pred_outcome <- prediction_multiple(
    ml_fit = ml_fit,
    feature_df = data.frame(
      gc_content = gc_content,
      RNA_type = as.character(RNA_type),
      RNA_region = as.character(RNA_region),
      exon_length = exon_length,
      distance_to_junction = distance_to_junction,
      evolutionary_conservation = evolutionary_conservation,
      DNA_5mer = as.character(DNA_5mer),
      stringsAsFactors = FALSE
    ),
    positive_threshold = positive_threshold
  )

  return(list(
    predicted_m6A_prob = pred_outcome$predicted_m6A_prob[[1]],
    predicted_m6A_status = pred_outcome$predicted_m6A_status[[1]]
  ))
}
