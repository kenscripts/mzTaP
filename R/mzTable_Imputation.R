#' groupMVs_2_zero
#'
#' Replaces MVs absent in all group samples with 0
#'
#' @param MV_DF Dataframe with missing values, where rows are mass features and columns are samples
#' @param GRP_PATTERNS Patterns to identify group columns
#' @return Dataframe with missing values replaced by zeros in group columns
#' @export
groupMVs_2_zero <- function(
		            MV_DF,
			    GRP_PATTERNS
			    ){
                              # Initialize output dataframe
                              IMP_DF <- MV_DF

                              # For each group, identify columns matching the group pattern
                              sapply(
                                     X = GRP_PATTERNS,
                                     FUN = function(GRP_PATTERN){
                                                                 # Get group columns
                                                                 GRP_COLS <- grep(
						                                  pattern = GRP_PATTERN,
						                                  x = colnames(MV_DF)
    						                                  )

                                                                 # Get group data
                                                                 GRP_DATA <- MV_DF[, GRP_COLS]

                                                                 # Identify rows where all values in group data are NA
                                                                 ALL_NA_ROWS <- apply(
						                                      X = GRP_DATA,
						                                      MAR = 1,
						                                      FUN = function(X) all(
									                                    is.na(X)
									                                    )
						                                      )

				                                 # Replace the group MVs with 0
                                                                 # <<- to change global MIN_DF variable
                                                                 IMP_DF[ALL_NA_ROWS, GRP_COLS] <<- 0
				                                 }
                                     )

                              return(IMP_DF)
                              }


#' replaceMVs_with_groupMin 
#'
#' Replaces MVs in group with the minimum value of group
#'
#' @param MV_DF Dataframe with missing values, where rows are mass features and columns are samples
#' @param GRP_PATTERNS Patterns to identify group columns
#' @return Dataframe with missing values replaced by group minimum
#' @export
replaceMVs_with_groupMin <- function(
		                     MV_DF,
			             GRP_PATTERNS
        			     ){
                                       # Initialize output dataframe
                                       MIN_DF <- MV_DF

                                       # For each group, replace MVs with group min
                                       sapply(
                                              X = GRP_PATTERNS,
                                              FUN = function(GRP_PATTERN){
                                                                          # Identify columns matching the group pattern
                                                                          GRP_COLS <- grep(
        				                                                   pattern = GRP_PATTERN,
        						                                   x = colnames(MV_DF)
        						                                   )
        
                                                                          # Get group data
                                                                          GRP_DATA <- MV_DF[, GRP_COLS]

                                                                          # Generate imputed rows
                                                                          MIN_ROWS <- apply(
        						                                    X = GRP_DATA,
        						                                    MAR = 1,
          						                                    FUN = function(X){
								                                              # find row min
								                                              X_MIN <- min(
        									                                           X,
         									                                           na.rm = TRUE
        									                                           )

								                                              # replace row MVs with min
								                                              X[is.na(X)] <- X_MIN

								                                              return(X)
							                                                      }
        						                                    ) %>%
						                                      # apply puts samples in x-dim and features in y-dim
						                                      # need to transpose
					                                              t()

                                                                          # replace values in output dataframe
									  # <<- to change global MIN_DF variable
                                                                          MIN_DF[, GRP_COLS] <<- MIN_ROWS
					                                  }
                                             )

                                       return(MIN_DF)
                                       }
