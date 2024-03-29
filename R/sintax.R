#' @name sintaxFilter
#' @title Filter SINTAX taxonomy
#'
#' @description Filtering the taxonomy based on SINTAX scores.
#'
#' @param taxonomy.tbl A data.frame with SINTAX taxonomy results
#' @param sintax.threshold The threshold for accepting a classification
#'
#' @details A \code{taxonomy.tbl} from SINTAX has a score-column for each of the
#' ranks, e.g. the \code{domain} rank has a \code{domain_score}. If a score is below
#' the \code{sintax.threshold} the corresponding taxon is set to \code{NA}.
#'
#' Finally, all score-columns are removed from the table.
#'
#' @return A table with pure taxonomy results after filtering. The score columns
#' have been removed.
#'
#' @author Lars Snipen.
#'
#' @importFrom dplyr %>% mutate if_else select
#' @importFrom tidyselect contains
#'
#' @export sintaxFilter
#'
sintaxFilter <- function(taxonomy.tbl, sintax.threshold = 0.5){
  taxonomy.tbl <- taxonomy.tbl %>%
    mutate(across(contains("_score"), function(x){x = if_else(is.na(x), 0, as.numeric(x))})) %>%
    mutate(across(!contains("_score"), as.character)) %>%
    mutate(domain = if_else(domain_score >= sintax.threshold, domain, NA_character_)) %>%
    mutate(phylum = if_else(phylum_score >= sintax.threshold, phylum, NA_character_)) %>%
    mutate(class = if_else(class_score >= sintax.threshold, class, NA_character_)) %>%
    mutate(order = if_else(order_score >= sintax.threshold, order, NA_character_)) %>%
    mutate(family = if_else(family_score >= sintax.threshold, family, NA_character_)) %>%
    mutate(genus = if_else(genus_score >= sintax.threshold, genus, NA_character_)) %>%
    mutate(species = if_else(species_score >= sintax.threshold, species, NA_character_)) %>%
    select(!contains("_score"))
  return(taxonomy.tbl)
}
