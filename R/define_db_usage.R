#' Define the database processing order
#' @description Database processing order is defined by a numeric vector with one value for each available database DB1 to DBx.
#'
#' Possible values:
#' \describe{
#'  \item{1}{most important}
#'  \item{2}{2nd most important}
#'  \item{x}{xth most important}
#'  \item{0}{excluded}
#' }
#' i.e. c(1,2,1,0) means DB1 and DB3 are processed first, then DB2, DB4 isn't used at all
#'
#' @param ARG_LIST list with user settings created by prepare_input_arguments.R
#' @param ranking numeric vector with database processing (=rank) order
#' @return updated ARG_LIST
#' @export
define_db_usage <- function(ARG_LIST, ranking=NA) {
  db.table <- list_available_dbs(ARG_LIST)
  if (is.null(db.table)) stop("No usable databases exist.")

  if (missing(ranking)) {
    catn("Define ranks for all available databases DB1 to DBx! Available ranks:")
    catn("\t1 = most important\n\t2 = 2nd most important\n\tx = xth most important\n\t0 = excluded")
    catn("Assignment follows the present DB order, i.e. for DB1 to DB4:\n rank order '1 2 1 0' means DB1 and DB3 share rank 1, DB2 has got rank 2, DB4 is not used")
    catn("Available databases:")
    print(db.table)
    DB_RANKS <- readline("Define ranks in correct order (i.e. '1 2 3' for 3 dbs) :")
    ranking <- as.numeric(unlist(strsplit(DB_RANKS, " ")))

  }
  if( !is.numeric(ranking) || length(ranking)!=nrow(db.table) || any(is.na(ranking)) || min(ranking)<0 ||
      max(ranking)<1       || max(ranking)>nrow(db.table))   stop("Wrong rank definition")

  db.table$rank <- ranking
  ARG_LIST$DB.TABLE <- db.table[ranking>0,]

  save(ARG_LIST, file=ARG_LIST$pasteIN("ARG_LIST.Rda"))

  return(ARG_LIST)
}
