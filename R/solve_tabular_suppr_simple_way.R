#' Find all the solutions to guess the missing values in a suppressed tabular data
#'
#' @param tab matrix. Inner cells of a tabular data with NA to indicate the suppressed counts
#' @param row_totals integer vector. Row margins.
#' @param col_totals integer vector. Column margins.
#'
#' @returns 
#' List of two components: 
#' - `solutions` : all the possible solutions
#' - `positions` : description of the positions
#' @export
#'
#' @examples
#' row_totals <- c(80,80,40)
#' col_totals <- c(80,27,48,45)
#' 
#' tab <- matrix(c(NA, 7,  NA, 18,
#'                 NA,12,  NA, 22,
#'                 15, 8,  12,  5), nrow=3, byrow=TRUE)
#' 
#' res2 <- find_solutions_bruteforce(tab, row_totals, col_totals)
find_solutions_bruteforce <- function(tab, row_totals, col_totals) {
  
  positions <- which(is.na(tab), arr.ind = TRUE)
  n_unknown <- nrow(positions)
  
  explore <- function(current_tab, pos_idx) {
    # Cas terminal: toutes les positions NA sont remplies
    if (pos_idx > n_unknown) {
      if (all(rowSums(current_tab) == row_totals) &&
          all(colSums(current_tab) == col_totals)) {
        return(list(current_tab))
      } else {
        return(list())
      }
    }
    i <- positions[pos_idx, 1]
    j <- positions[pos_idx, 2]
    
    # Valeurs maximales restant à attribuer pour cette case
    max_val_row <- row_totals[i] - sum(current_tab[i, ], na.rm = TRUE)
    max_val_col <- col_totals[j] - sum(current_tab[, j], na.rm = TRUE)
    max_val <- min(max_val_row, max_val_col)
    # Pas de valeurs négatives ou nulles
    possible_vals <- if (max_val > 0) 1:max_val else integer(0)
    
    solutions <- list()
    for (val in possible_vals) {
      new_tab <- current_tab
      new_tab[i, j] <- val
      
      if (all(rowSums(new_tab, na.rm = TRUE) <= row_totals) &&
          all(colSums(new_tab, na.rm = TRUE) <= col_totals)) {
        sols <- explore(new_tab, pos_idx + 1)
        solutions <- c(solutions, sols)
      }
    }
    return(solutions)
  }
  
  results <- explore(tab, 1)
  list(solutions = results, positions = positions)
}


#' Compute the probability that a missing value is a sensitive count
#'
#' @param all_solutions result of the `find_solutions_general`
#' @param min_conf integer. lower bound of the sensitive values interval 
#' @param max_conf integer. upper bound of the sensitive values interval 
#'
#' @returns data.frame of the probability for each cell
#' @export
#'
#' @examples
#' row_totals <- c(29,45)
#' col_totals <- c(24,25,25)
#' tab <- matrix(c(NA, 15, NA, 
#'                 NA,10,  NA), 
#'               nrow=2, byrow=TRUE)
#' 
#' res <- find_solutions_bruteforce(tab, row_totals, col_totals)
#' post_res <- posterior_prob_confidential(
#'   res,
#'   max_conf = 10
#' )
posterior_prob_confidential <- function(all_solutions,
                                        min_conf=1, max_conf=10) {
  sols <- all_solutions$solutions
  positions <- all_solutions$positions
  
  # Are excluded the solutions with no confidential case
   is_sols_possible <- sols |> 
    purrr::map(
      \(sol){
        is_confidential_possible <- vector(mode="logical", length=nrow(positions))
        for(r in seq_len(nrow(positions))){
          i <- positions[r,"row"]
          j <- positions[r,"col"]
          is_confidential_possible[r] <- sol[i,j] <= max_conf & sol[i,j] >= min_conf
        }
        sum(is_confidential_possible) > 0
      } ) |> purrr::list_c()
  sols_possible_with_confidential_cell <- sols[is_sols_possible]
  
  nsol <- length(sols_possible_with_confidential_cell)
  
  if (nsol == 0) stop("Aucune solution valide !") else cat("Nb de solutions avec une case confidentielle au moins: ", nsol)
  
  prob_conf <- rep(0, nrow(positions))
  
  for (s in sols_possible_with_confidential_cell) {
    vals <- s[positions]
    prob_conf <- prob_conf + as.numeric(vals >= min_conf & vals <= max_conf)
  }
  
  prob_conf <- prob_conf / nsol
  return(data.frame(ligne=positions[,1],
                    colonne=positions[,2],
                    proba_confidentielle=prob_conf))
}

