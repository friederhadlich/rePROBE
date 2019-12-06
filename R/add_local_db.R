#' Add a RNA database from a local source into this project
#'
#' @param ARG_LIST list with user settings created by prepare_input_arguments()
#' @param db.name Database name
#' @param rna.file File in gff3 or fasta format with RNA information
#' @param update Renew existing data
#' @import data.table fs
#' @importFrom Rbowtie bowtie_build bowtie
#' @export
add_local_db <- function(ARG_LIST, db.name = NA, rna.file = NA, update = F) {

  if (!dir.exists(ARG_LIST$folder.in))  stop("Wrong input for add_local_db!")

  with(ARG_LIST, {

    if (is.na(db.name)) db.name <- readline("\tDatabase name: ")
    if (!grepl("^[0-9a-zA-Z_.-]+$", db.name)) stop("Database name invalid!")
    if ( update && (dir.exists(pasteIN(db.name))) ) unlink(pasteIN(db.name))

    if (!file_exists(pasteIN(db.name, "sequences.fa"))) {

      catn("Add/update local database",db.name)
      if (is.na(rna.file)) tryCatch({
        cat.step("Choose either a sequence file (unzipped fasta) or an annotation file (unzipped gff3) without comment lines!")
        tmp <- file.choose()
      }, error = function(e) { tmp <<- 0 })
      if (!grepl("\\.(fa|fasta|gff3)$", tmp))  stop("Wrong file type! Only 'fa', 'fasta' or 'gff3' format accepted!")

      dir_create(pasteIN(db.name))

      # RNA
      if (grepl("\\.(fa|fasta)$", tmp)) {
        link_create(tmp, pasteIN(db.name, "sequences.fa"), is_linux())
      } else {
        # ANNO+DNA(+SNP)
        # anno
        link_create(tmp, pasteIN(db.name, "annotation.gff3"), is_linux())
        # dna
        cat.step("You need a genome file AND a annotation file to generate a RNA sequence file.")
        cat.step("Choose a reference genome file (unzipped fasta)!")
        tmp <- file.choose()
        if (!grepl("\\.(fa|fasta)$", tmp)) { dir_delete(pasteIN(db.name)); stop("Wrong file type!") }
        link_create(tmp, pasteIN(db.name, "genome.fa"), is_linux())
        X <- fread(pasteIN(db.name, "annotation.gff3"), sep="\t")[V3=="exon"]
        writeLines(get.fasta(X, pasteIN(db.name, "genome.fa")), pasteIN(db.name, "sequences.fa"))
        # snp
        cat.step("Choose optional a SNP file (unzipped vcf) with chromosome names matching the genome file!")
        tryCatch(tmp <- file.choose(), error = function(e) { tmp <<- 0 })
        if (grepl("\\.vcf$", tmp)) link_create(tmp, pasteIN(db.name, "snp.vcf"), is_linux())
        # if (!STRICT_MAPPING && file.exists(pasteIN(db.name, "snp.vcf")) && !dir_exists(pasteIN(db.name, "vcf"))) {
        #   catn("\t... prepare variant calling file")
        #   dir_create(pasteIN(db.name, "vcf"))
        #   prepare_vcf(pasteIN(db.name,"snp.vcf"), .Platform$file.sep)
        # }
        prepare.vcf(ARG_LIST, db.name)
      }
    }

    return(db.name)
  })

}
