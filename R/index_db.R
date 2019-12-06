#' Add bowtie index to an existing database
#'
#' @param ARG_LIST list with user settings created by prepare_input_arguments.R
#' @param db.name Database name
#' @importFrom Rbowtie bowtie
#' @export
index_db <- function(ARG_LIST, db.name) {

  with(ARG_LIST, {

    if (!file.exists(pasteIN(db.name, "sequences.fa"))) stop("File 'sequences.fa' is missing. Cannot create an index!")
    if (any(grepl("*.sa", list.files(pasteIN(db.name, "bwt"))))) unlink(pasteIN(db.name, "bwt"), recursive = T)
    if (file.exists(pasteIN(db.name, "bwt", "sequences.4.ebwt")) &&
        as.POSIXct(file_info(pasteIN(db.name, "bwt", "sequences.4.ebwt"))$change_time) > as.POSIXct(file_info(pasteIN(db.name, "sequences.fa"))$change_time))
      return(TRUE)

    cat.step("... build bowtie index!")
    #if (NCPU>1 && any(grepl("--threads", bowtie_build("dummy", "dummy", force = TRUE, usage = TRUE, strict = FALSE))))
    if (is_linux())
      bowtie_build(references = pasteIN(db.name,"sequences.fa"), prefix = "sequences", outdir = pasteIN(db.name, "bwt"), force = T, threads=NCPU)
    else
      bowtie_build(references = pasteIN(db.name,"sequences.fa"), prefix = "sequences", outdir = pasteIN(db.name, "bwt"), force = T)
    #else
    #  bowtie_build("sequences.fa", prefix = "sequences", outdir = "bwt", force = T)

  })

  return(TRUE)
}
