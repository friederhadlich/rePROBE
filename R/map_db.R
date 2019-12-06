#' Map sequences in 'probes.fasta' against a custom database reference
#' @param ARG_LIST list with user settings created by prepare_input_arguments.R
#' @param db.name Database name
#' @import data.table fs
#' @export
map_db <- function(ARG_LIST, db.name) {
  with(ARG_LIST, {

    if (!file.exists(pasteIN(db.name, "bwt", "sequences.4.ebwt"))) stop("Bowtie index is missing. Cannot run bowtie mapping!")
    if (!file.exists(pasteIN("probes.fasta"))) stop("File 'probes.fasta' is missing. Cannot run bowtie mapping!")

    run.bowtie <- function(...)
      bowtie( sequences = pasteIN("probes.fasta"),
              outfile = pasteIN(db.name, "bowtie.sam"), f=T, y=T, a=T, S=T, best=T, strata=T, p=NCPU, "sam-nohead"=T, force=T, ...)

    cat.step("... bowtie mapping!")
    if(!dir_exists(pasteIN(db.name, "vcf")) || STRICT_MAPPING)
      run.bowtie(index = pasteIN(db.name, "bwt", "sequences"), v=0, norc=db.name!="dna")
    else run.bowtie(index = pasteIN(db.name, "bwt", "sequences"), norc=db.name!="dna")

  })
  return(TRUE)
}
