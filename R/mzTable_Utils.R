#' build_group_vec
#'
#' Generates vector of group names from sample names. Sample length as colnames(DATA).
#'
#' @param DATA Dataframe, where rows are mass feaures and columns are samples
#' @param GRP_PATTERNS Patterns to identify group columns
#' @return GRP_VEC Vector of group names
#' @export
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
				      FUN = function(GRP_PATTERN){
				                                  GRP_COLS <- grep(
					                                           pattern = GRP_PATTERN,
					                                           x = colnames(DATA)
					                                           )
					                          # <<- to change global GRP_VEC variable
					                          GRP_VEC[GRP_COLS] <<- GRP_PATTERN
					                          }
			              )

                               return(GRP_VEC)
                               }
