#' build_group_vec
#'
#' generates vector of group names from sample names
#'
#' @param DATA Dataframe, with rows of features and columns of samples
#' @param GRP_PATTERNS Patterns to identify group columns
#' @return GRP_VEC Vector of group names
build_group_vec  <- function(
  			     DATA,
			     GRP_PATTERNS
  			     ){
                               # Initialize group vector with length equal to number of samples
                               GRP_VEC <- vector(
                                                 mode = "character",
                             		         length = length(
                             				         colnames(DATA)
                             				         )
                             		         )

                               # Add group names to appropriate positions in GRP_VEC
                               sapply(
			              X = GRP_PATTERNS,
				      FUN = function(PATTERN){
				                              PATTERN_COLS <- grep(
					                                           pattern = PATTERN,
					                                           x = colnames(DATA)
					                                           )
					                      # <<- to change global GRP_VEC variable
					                      GRP_VEC[PATTERN_COLS] <<- PATTERN
					                      }
			              )

                               return(GRP_VEC)
                               }
