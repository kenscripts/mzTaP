#' get_mf_blankInfo
#'
#' Get information on mass features from blanks. Uses BLANK_PATTERN to identify blank columns. Unmatched columns are interpreted as sample columns.
#'
#' @param MF_DF Dataframe, where rows are mass feaures and columns are samples
#' @param BLANK_PATTERN Pattern to identify blank samples
#' @return Dataframe with mass feature information from blanks
#' @export
get_mf_blankInfo <- function(
		 	     MF_DF,
			     BLANK_PATTERN
			     ){
                               # Identify columns matching the blank pattern
                               BLANK_COLS <- grep(
      			                          pattern = BLANK_PATTERN,
			                          x = colnames(MF_DF)
						  )

                               # Identify sample columns
                               SAMPLE_COLS <- grep(
      					           pattern = BLANK_PATTERN,
  						   x = colnames(MF_DF),
						   invert = TRUE
						   )

                               # Get mass feature info in blanks
                               MF_BLANK_INFO <- apply(
				                      X = MF_DF,
					              MAR = 1,
					              FUN = function(X){
                                                                        # get max intensity of each mass feature
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

								        # mass feature found in blanks
								        DETECTED_IN_BLANK <- any(
											         X[BLANK_COLS] > 0
											         )

								        # mass feature missing in samples
								        MISSING_IN_SAMPLES <- !any(
											           X[SAMPLE_COLS] > 0
												   )

								        # if both true then mass feature is blank-specific
								        BLANK_FEATURES <- MISSING_IN_SAMPLES && DETECTED_IN_BLANK

                                                                        # is blank median > sample median
                                                                        BLANK_GREATER <- BLANK_MED > SAMPLE_MED
  
                                                                        # sample median:blank median
                                                                        MED_RATIO <- round(
								                           x = SAMPLE_MED/BLANK_MED,
								                           digits = 2
								                           )

                                                                        # return list with info
                                                                        return(
				                                               c(
                                                                                 max_intensity = MAX_INT,
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
#' For each mass feature, subtract median value in blanks from each sample. Uses BLANK_PATTERN to identify blank columns. Uses GRP_PATTERN to identify sample columns. Blank median is rounded up to nearest integer
#'
#' @param DATA Dataframe, where rows are mass feaures and columns are samples
#' @param BLANK_PATTERN Pattern to identify blank 
#' @param GRP_PATTERN Pattern to identify experimental samples
#' @return Dataframe with mass feature intensities subtracted by blank median
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

                                    # Calculate the median intensity in blank columns for each feature
                                    BLANK_MED <- apply(
						       X = DATA[,BLANK_COLS],
						       MAR = 1, 
						       FUN = function(X) median(
						                                X,
						                                na.rm = TRUE
						                                ) %>%
						                         # round up to nearest integer
						                         ceiling()
						       )

				    # copy dataframe for output
				    DATA.DENOISE <- DATA

                                    # Subtract the median intensity in blanks from samples
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
#' Determine the prevalence of each mass feature in sample groups. If mass feature is not conserved in any sample group (detected in all replicates of sample group) then it is considered noisy. If mass feature is missing in any sample group replicates then it's considered group noisy.
#'
#' @param MF_DF Mass feature dataframe, where rows are mass feaures and columns are samples
#' @param GRP_PATTERNS Patterns to identify sample groups 
#' @return Dataframe with information on mass feature noisiness.
#' @export
get_mf_noiseInfo <- function(
		             MF_DF,
			     GRP_PATTERNS
        		     ){
                               # get mf prevalence in sample groups
                               PREVALENCE_LIST <- lapply(
                                                         X = GRP_PATTERNS,
                                                         FUN = function(GRP_PAT){
                                                                                 # Identify columns matching the sample group pattern
                                                                                 GRP_COLS <- grep(
        				                                                          pattern = GRP_PAT,
        						                                          x = colnames(MF_DF)
        						                                          )
        
                                                                                 # Get group data
                                                                                 GRP_DATA <- MF_DF[, GRP_COLS]

                                                                                 # Get prevalence of mf in group
                                                                                 GRP_PREVALENCE <- apply(
        						                                                 X = GRP_DATA,
        						                                                 MAR = 1,
          						                                                 FUN = function(X){
								                                                           # number of group samples with feature
								                                                           SAMPLE_COUNT <- sum(X > 0)
 
 								                                                           # total number of group samples
 								                                                           TOTAL_SAMPLES <- length(colnames(GRP_DATA))
 
								                                                           # prevalence of feature in sample group
								                                                           MF_PREVALENCE <- SAMPLE_COUNT / TOTAL_SAMPLES

								                                                           return(MF_PREVALENCE)
							                                                                   }
          						                                                 )

                                                                                 return(GRP_PREVALENCE)
					                                         }
                                                         )

                               # add sample group names
                               names(PREVALENCE_LIST) <- paste(
							       "group",
							       GRP_PATTERNS,
							       sep = "."
							       )

                               # convert to dataframe
                               PREVALENCE_DF <- as.data.frame(PREVALENCE_LIST)

			       # get mass feature noisiness
			       NOISY_DF <- apply(
				                 X = PREVALENCE_DF,
				                 MAR = 1,
				                 FUN = function(X){
			                                           # check if mass feature is absent
						                   IS_ABSENT <- all(
				                                                    X == 0
				                                                    )

		 	                                           # check to see if mass feature is conserved in any group
			                                           # group-specific mass features should be present in all replicates
				                                   IS_NOISY <- !any(
				                                                    X == 1.0
				                                                    )

			                                           # check to see if mass feature is noisy in any group (not absent and not conserved in any group)
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


#' get_noisyMFtable
#'
#' Return long dataframe of noisy mass features for each sample group.
#'
#' @param MF_DF Mass feature dataframe, where rows are mass feaures and columns are samples
#' @param GRP_PATTERNS Patterns to identify sample groups 
#' @return Long dataframe of noisy mass features in each sample group
#' @export
get_noisyMFtable <- function(
		             MF_DF,
			     GRP_PATTERNS
        		     ){
	                       # get mf noise info
                               MF_NOISE_INFO <- get_mf_noiseInfo(
                                                                 MF_DF = MF_DF,
                                                                 GRP_PATTERNS = GRP_PATTERNS
                                                                 )


	                       # get noisy mf for each group
                               NOISYMF_LONG <- lapply(
                                                      X = GRP_PATTERNS,
                                                      FUN = function(GRP_PAT) {
                                                                               # get group columns
                                                                               GRP_COLS <- grep(
                                                                                                pattern = GRP_PAT,
                                                                                                x = colnames(MF_DF)
                                                                                                )

                                                                               # get group data
                                                                               GRP_DATA <- MF_DF[,GRP_COLS]

                                                                               # get group mf noise
                                                                               MF_NOISE.GRP <- MF_NOISE_INFO %>%
                                                                                               select(
												      contains(GRP_PAT)
												      ) %>%
                                                                                               pull(.) %>%
                                                                                               {. > 0 & . < 1}

                                                                               # get grp noisy mf long table
                                                                               GRP_NOISYMF_DF <- GRP_DATA[MF_NOISE.GRP,] %>%
                                                                                                 rownames_to_column(var = "mf_idx") %>%
                                                                                                 pivot_longer(
                                                                                                              cols = contains("_R"),
                                                                                                              names_to = "sample_id",
                                                                                                              values_to = "intensity"
                                                                                                              ) %>%
                                                                                                 filter(intensity > 0)

                                                                               return(GRP_NOISYMF_DF)
                                                                               }
                                                    ) %>%
                                               bind_rows() %>%
                                               as.data.frame()

                               return(NOISYMF_LONG)
                               }
