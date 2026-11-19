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
  nn <- nchar(dna_strings[1])
  seq_m <- matrix(unlist(strsplit(dna_strings, "")), ncol = nn, byrow = TRUE)
  colnames(seq_m) <- paste0("nt_pos", 1:nn)
  seq_df <- as.data.frame(seq_m)
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
  stopifnot(all(c("gc_content", "RNA_type", "RNA_region", "exon_length", "distance_to_junction",
                  "evolutionary_conservation", "DNA_5mer") %in% colnames(feature_df))) # Check errors if incorrect column names of input data.frame
  feature_df$RNA_type <- factor(feature_df$RNA_type, levels = c("mRNA", "lincRNA", "lncRNA", "pseudogene"))
  feature_df$RNA_region <- factor(feature_df$RNA_region, levels = c("CDS", "intron", "3'UTR", "5'UTR"))
  infered_prob <- predict(ml_fit, newdata = cbind(feature_df, dna_encoding(feature_df$DNA_5mer)), type = "prob")
  feature_df$predicted_m6A_prob <- infered_prob[,"Positive"]
  feature_df$predicted_m6A_status <- "Negative"
  feature_df$predicted_m6A_status[feature_df$predicted_m6A_prob > positive_threshold] <- "Positive"
  return(feature_df) # Return a data.frame with predicted m6A prob and predicted m6A status
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
#' @return A named character vector containing predicted_m6A_prob
#' and predicted_m6A_status.
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
  pred_outcome <- prediction_multiple(ml_fit = ml_fit,
                                      data.frame(gc_content = gc_content,
                                                 RNA_type = factor(RNA_type, levels = c("mRNA", "lincRNA", "lncRNA", "pseudogene")),
                                                 RNA_region = factor(RNA_region, levels = c("CDS", "intron", "3'UTR", "5'UTR")),
                                                 exon_length = exon_length,
                                                 distance_to_junction = distance_to_junction,
                                                 evolutionary_conservation = evolutionary_conservation,
                                                 DNA_5mer = DNA_5mer),
                                      positive_threshold = positive_threshold)
  returned_vector <- c(pred_outcome$predicted_m6A_prob, pred_outcome$predicted_m6A_status)
  names(returned_vector) <- c("predicted_m6A_prob", "predicted_m6A_status")
  return(returned_vector) # Return a named vector with predicted m6A prob and predicted m6A status
}
