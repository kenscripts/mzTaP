#' get_mf_blankInfo
#'
#' Get information on mass features from blanks.
#' Uses BLANK_PATTERN to identify blank columns.
#' Unmatched columns are interpreted as sample columns.
#'
#' @param DATA Dataframe, where rows are mass feaures and columns are samples
#' @param BLANK_PATTERN Pattern to identify blank samples
#' @return Dataframe with mass feature information from blanks
#' @export
get_mf_blankInfo <- function(
		 	     DATA,
			     BLANK_PATTERN
			     ){
                               # Identify columns matching the blank pattern
                               BLANK_COLS <- grep(
      			                          pattern = BLANK_PATTERN,
			                          x = colnames(DATA)
						  )

                               # Identify sample columns
                               SAMPLE_COLS <- grep(
      					           pattern = BLANK_PATTERN,
  						   x = colnames(DATA),
						   invert = TRUE
						   )

                               # Get mf info in blanks
                               MF_BLANK_INFO <- apply(
				                      X = DATA,
					              MAR = 1,
					              FUN = function(X){
                                                                        # get max intensity for each feature
                                                                        MAX_INT <- max(
								                       X,
						                                       na.rm = TRUE
								                       )

                                                                        # calculate the median for blank columns
                                                                        BLANK_MED <- median(
								                            X[BLANK_COLS],
								                            na.rm = TRUE
      								                            )

                                                                        # calculate sample median
                                                                        SAMPLE_MED <- median(
						 		                             X[SAMPLE_COLS],
  						 		                             na.rm = TRUE
								                             )

								        # feature found in blanks
								        DETECTED_IN_BLANK <- any(
											         X[BLANK_COLS] > 0
											         )

								        # feature missing in samples
								        MISSING_IN_SAMPLES <- !any(
											           X[SAMPLE_COLS] > 0
												   )

								        # if both true then feature is blank-specific
								        BLANK_FEATURES <- MISSING_IN_SAMPLES && DETECTED_IN_BLANK

                                                                        # is blank median > sample median
                                                                        BLANK_GREATER <- BLANK_MED > SAMPLE_MED
  
                                                                        # experiment group median:blank median
                                                                        MED_RATIO <- round(
								                           x = SAMPLE_MED/BLANK_MED,
								                           digits = 2
								                           )

                                                                        # Return list with info
                                                                        return(
				                                               c(
                                                                                 mf_max_intensity = MAX_INT,
                                                                                 blank_median = BLANK_MED,
                                                                                 sample_median = SAMPLE_MED,
                                                                                 only_in_blank = BLANK_FEATURES,
                                                                                 is_blank_greater = BLANK_GREATER,
                                                                                 sampleMedian_to_blankMedian = MED_RATIO
                                                                                 )
                                                                               )
                                                                        }
					              ) %>%
			                        # apply switches rows & columns
			                        t()

                               # output dataframe
                               MF_BLANK_INFO_DF <- as.data.frame(MF_BLANK_INFO)
                               MF_BLANK_INFO_DF$only_in_blank <- as.logical(MF_BLANK_INFO_DF$only_in_blank)
                               MF_BLANK_INFO_DF$is_blank_greater <- as.logical(MF_BLANK_INFO_DF$is_blank_greater)

                               return(
                                      MF_BLANK_INFO_DF
                                      )
                               }


#' subtract_blank_median
#'
#' For each mass feature, subtract median value in blanks from each sample.
#' Uses BLANK_PATTERN to identify blank columns.
#' Uses GRP_PATTERN to identify sample columns.
#'
#' @param DATA Dataframe, where rows are mass feaures and columns are samples
#' @param BLANK_PATTERN Pattern to identify blank 
#' @return Dataframe with mass feature information from blanks
#' @export
subtract_blank_median <- function(
				  DATA,
				  BLANK_PATTERN,
				  GRP_PATTERN
				  ){
                                    # Identify columns matching the blank pattern
                                    BLANK_COLS <- grep(
						       pattern = BLANK_PATTERN,
						       x = colnames(DATA),
						       value = TRUE
						       )

                                    # Identify columns matching the sample group pattern
                                    GRP_COLS <- grep(
      				  		     pattern = GRP_PATTERN,
  						     x = colnames(DATA),
						     value = TRUE
						     )

                                    # Calculate the median for each blank column
                                    BLANK_MED <- apply(
						       X = DATA[,BLANK_COLS],
						       MAR = 1, 
						       FUN = function(X) median(
						                                X,
						                                na.rm = TRUE
						                                )
						       )

				    # Copy dataframe for output
				    DATA.DENOISE <- DATA

                                    # Subtract the median of the corresponding blank column from each exp group column
				    sapply(
					   X = GRP_COLS,
					   FUN = function(COL){
                                                               DATA.DENOISE[,COL] <<- DATA.DENOISE[,COL] - BLANK_MED
					                       }
					   )

				    # remove blank columns
                                    DATA.DENOISE <- DATA.DENOISE %>%
                                                    select(
							   !contains(BLANK_COLS)
							   )

				    # transform negative values from blank subtraction to 0
                                    DATA.DENOISE[DATA.DENOISE < 0] <- 0

                                    return(DATA.DENOISE)
                                    }


#' get_mf_noiseInfo
#'
#' Determine the prevalence of each mass feature in sample groups.
#' If mass feature is not conserved in any sample group (detected in all replicates of sample group) then it is considered noisy.
#' If mass feature is missing in any sample group replicates then it's considered group noisy.
#'
#' @param DATA Dataframe, where rows are mass feaures and columns are samples
#' @param GRP_PATTERNS Patterns to identify sample groups 
#' @return Dataframe with information on mass feature noisiness.
#' @export
#' @export
get_mf_noiseInfo <- function(
		             DATA,
			     GRP_PATTERNS
        		     ){
                               # get mf prevalence in sample groups
                               PREVALENCE_LIST <- lapply(
                                                         X = GRP_PATTERNS,
                                                         FUN = function(GRP_PATTERN){
                                                                                     # Identify columns matching the group pattern
                                                                                     GRP_COLS <- grep(
        				                                                              pattern = GRP_PATTERN,
        						                                              x = colnames(DATA)
        						                                              )
        
                                                                                     # Get group data
                                                                                     GRP_DATA <- DATA[, GRP_COLS]

                                                                                     # Get prevalence of mf in group
                                                                                     GRP_PREVALENCE <- apply(
        						                                                     X = GRP_DATA,
        						                                                     MAR = 1,
          						                                                     FUN = function(X){
								                                                               # number of group samples with feature
								                                                               SAMPLE_COUNT <- sum(X > 0)
 
 								                                                               # total number of group samples
 								                                                               TOTAL_SAMPLES <- length(colnames(GRP_DATA))
 
								                                                               # prevalence of feature in group
								                                                               MF_PREVALENCE <- SAMPLE_COUNT / TOTAL_SAMPLES

								                                                               return(MF_PREVALENCE)
							                                                                       }
          						                                                     )

                                                                                     return(GRP_PREVALENCE)
					                                             }
                                                         )

                               # add names
                               names(PREVALENCE_LIST) <- paste(
							       "group",
							       GRP_PATTERNS,
							       sep = "."
							       )

                               # convert to dataframe
                               PREVALENCE_DF <- as.data.frame(PREVALENCE_LIST)

			       # get mf noisiness
			       NOISY_DF <- apply(
				                 X = PREVALENCE_DF,
				                 MAR = 1,
				                 FUN = function(X){
			                                           # check if feature is absent
						                   IS_ABSENT <- all(
				                                                    X == 0
				                                                    )

		 	                                           # check to see if feature is conserved in any group
			                                           # group-specific features should be present in all replicates
				                                   IS_NOISY <- !any(
				                                                    X == 1.0
				                                                    )

			                                           # check to see if feature is noisy in any group (not absent and not conserved in any group)
				                                   IS_GRP_NOISY <- any(
										       X > 0 & X < 1.0
										       )
                                                                   return(
									  c(
									    IS_ABSENT,
									    IS_NOISY,
									    IS_GRP_NOISY
									    )
									  )
							           }
				                 ) %>%
			                   # apply puts samples in x-dim and features in y-dim
                                           # need to transpose
			                   t()
					       
                               # add colnames
                               colnames(NOISY_DF) <- c(
                                                       "is_absent",
                                                       "is_noisy",
                                                       "is_grpNoisy"
                                                       )

			       # combine output
			       # use cbind to maintain row order
                               MF_NOISE_DF <- cbind(
						    PREVALENCE_DF,
						    NOISY_DF
						    )

                               return(MF_NOISE_DF)
                               }
