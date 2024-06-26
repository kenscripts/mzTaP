% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mzTable_Utils.R
\name{build_group_vec}
\alias{build_group_vec}
\title{build_group_vec}
\usage{
build_group_vec(DATA, GRP_PATTERNS)
}
\arguments{
\item{DATA}{Dataframe, where rows are mass feaures and columns are samples}

\item{GRP_PATTERNS}{Patterns to identify group columns}
}
\value{
GRP_VEC Vector of group names
}
\description{
Generates vector of group names from sample names. Sample length as colnames(DATA).
}
