#' Prepare all input data for the revised probeset assignment
#'
#' The function wraps functions that can also be directly used in given order:
#' prepare_input_arguments, create_project, prepare_chip, add_local_db, add_remote_db, add_dna, index_db, map_db, define_db_usage
#'
#' @param WORK_PATH Working directory
#' @param PROJ_NAME Project name
#' @param CHIP_PATH Chip path with one pgf and one mps file
#' @param SPECIES Species name in latin language, i.e. "Sus_scrofa"
#' @param NCPU Number of threads to use
#' @param STRICT_MAPPING Avoid SNP correction of probes
#' @param ASK_QUESTIONS Answer questions automatically with their defaults
#' @export
#' @import data.table fs
#' @importFrom Rbowtie bowtie_build bowtie
#' @importFrom curl curl_download curl
#' @importFrom utils read.delim write.csv
#' @importFrom R.utils gunzip

prepare_data <- function(WORK_PATH,
                         PROJ_NAME,
                         CHIP_PATH,
                         SPECIES,
                         NCPU,
                         STRICT_MAPPING=T,
                         ASK_QUESTIONS=T) {

  ### CREATE INTRODUCTION: CREATE DIRS or LIST EXISTING DIRS #####################
  ###

  # usethis::use_package("RcppArmadillo", "LinkingTo")
  # roxygenise("/disk1/R/my_packages/test1")
  # usethis::use_rcpp()

  #ARG_LIST <- prepare_data(WORK_PATH="/disk1/R/reviseArrayAssignment/", PROJ_NAME="Pig1234", SPECIES="Sus_scrofa", NCPU=10, STRICT_MAPPING=F, ASK_QUESTIONS=F, CHIP_PATH = "/disk1/R/reviseArrayAssignment/CHIP_SSC/MO_PorGene-1_1-st-v1_strip_libraryfile/")
  #load("/disk1/R/reviseArrayAssignment/Pig1234/in/ARG_LIST.Rda")

  cat( sep = "",
       "#########################################################\n",
       "Welcome! This guide prepares the revised array assignment\n",
       "#########################################################\n\n")

  ARG_LIST <- do.call(prepare_input_arguments, as.list(match.call())[-1])

  create_project(ARG_LIST)

  prepare_chip(ARG_LIST)

  while(ARG_LIST$ASK_QUESTIONS) {
    if (readline("Import a local RNA database? y/[n]: ") != "y") break
    db.name <- add_local_db(ARG_LIST)
    index_db(ARG_LIST, db.name)
    map_db(ARG_LIST, db.name)
  }

  if (((!ARG_LIST$ASK_QUESTIONS && catn("Add or update remote RNA databases!")) ||
     readline("Add or update remote RNA databases? [y]/n: ")!="n")) {
    for (type in c("refseq", "ensembl")) {
      db.name <- add_remote_db(ARG_LIST, db.type = type)
      index_db(ARG_LIST, db.name)
      map_db(ARG_LIST, db.name)
    }
  }

  db.table <- list_available_dbs(ARG_LIST)
  dna.file <- vcf.file <- ""
  if (!ARG_LIST$ASK_QUESTIONS && any(grepl("ensembl", db.table$db))) {
    dna.file <- ARG_LIST$pasteIN(grep("ensembl", db.table$db, value = T), "genome.fa")
    vcf.file <- sub("genome.fa$","snp.vcf",dna.file)
  }
  if ((!ARG_LIST$ASK_QUESTIONS && catn("Add or update the local DNA database!")) ||
       readline("Add or update the local DNA database? [y]/n: ")!="n" ) {
    if (add_dna(ARG_LIST, dna.file, vcf.file)) {
      index_db(ARG_LIST, "dna")
      map_db(ARG_LIST, "dna")
    }
  }

  if (!ARG_LIST$ASK_QUESTIONS) {
    db.table <- list_available_dbs(ARG_LIST)
    if (is.null(db.table) || nrow(db.table)!=3)  stop("Number of databases is wrong!")
    ranking <- rep(1, 3); ranking[db.table$db %in% "dna"] <- 2
    ARG_LIST <- define_db_usage(ARG_LIST, ranking)
  } else
    ARG_LIST <- define_db_usage(ARG_LIST)

cat( sep = "", "\n",
     "##################################################################\n",
     "Preparation successful! Now start the analysis with command 'run'!\n",
     "##################################################################\n")

  return(ARG_LIST)
}
