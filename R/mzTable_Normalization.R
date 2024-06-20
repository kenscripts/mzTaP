#' get_groupRLA
#'
#' Get relative log abundances of mass features compared to mass feature median.
#'
#' @param DATA Dataframe, where rows are mass feaures and columns are samples
#' @param GRP_PATTERNS Patterns to identify sample groups
#' @return Dataframe with relative log abundances of mass features
#' @export
get_groupRLA <- function(
	 		 DATA,
			 GRP_PATTERNS
	                 ){
                           # Calculate the median for each blank column
                           DATA_RLA <- DATA

                           lapply(
				  X = GRP_PATTERNS,
				  FUN = function(PATTERN){
				                          # get group columns
				                          GRP_COLS <- grep(
			                                                   pattern = PATTERN, 
						                           colnames(DATA)
						                           )

				                          # get group data
						          GRP_DATA <- DATA %>%
						                      select(
						                             all_of(GRP_COLS)
						                             )

							  # get relative log abundances for each group for each feature
						          GRP_RLA <- apply(
						                           X = GRP_DATA,
						                           MAR = 1,
						                           FUN = function(X){
						                                             X_RLA <- X - median(X)

						                                             return(X_RLA)
						                                             }
						                           ) %>%
							             # apply puts samples in x-dim and features in y-dim
							             # need to transpose
							             t()

						          # add rla to output
						          DATA_RLA[, GRP_COLS] <<- GRP_RLA
				                          }
	                          )

                           return(DATA_RLA)
		           }


#' rowwise_pw_test
#'
#' Perform pairwise statistical test between groups.
#'
#' @param DATA Dataframe, where rows are mass feaures and columns are samples
#' @param GRP_VEC Vector of sample group names
#' @param TEST_NAME String of R function to perform pairwise test (e.g. t.test, kw.test)
#' @param PAIRED Boolean, if paired test then set to True (Default = False)
#' @return Dataframe with relative log abundances of mass features
rowwise_pw_test  <- function(
			     DATA,
			     GRP_VEC,
			     TEST_NAME,
			     PAIRED = False
			     ){
			       # initialize a vector to store p-values
			       TEST.PVALS <- numeric(
						     nrow(DATA)
						     )

			       # initialize test function
                               TEST_FN <- match.fun(TEST_NAME)

                               # get test p-values
                               TEST.PVALS <- apply(
					           X = DATA,
					           MAR = 1,
  					           FUN = function(x){
  					                             TEST_FN(
	                                                                     x ~ GRP_VEC,
	                                                                     paired = PAIRED
	                                                                     )$p.value
					                             }
					           )

                               return(TEST.PVALS)
                               }
