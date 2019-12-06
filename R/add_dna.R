#' Define a dna reference genome to be used as additional database, or change an existing one.
#' If parameters are missing, a menu will guide the user choice. Otherwise, user input is used!
#'
#' @param ARG_LIST list with user settings created by prepare_input_arguments()
#' @param dna.file File in fasta format with DNA information
#' @param vcf.file File in vcf format with SNP information
#' @param update Renew existing data
#' @export
add_dna <- function(ARG_LIST, dna.file = "", vcf.file = "", update = F) {

  if (!dir.exists(ARG_LIST$folder.in))  stop("Wrong input for add_dna!")

  f.seq <- ARG_LIST$pasteIN("dna","sequences.fa")
  if ( file_exists(f.seq) )
    if ( (!missing(update) && !update) || ( is_linux() && ( link_path(f.seq) == dna.file ||
      cat.step(paste0("DNA database exists and links to ", link_path(f.seq), "!")) &&
      ( !ARG_LIST$ASK_QUESTIONS && cat.step("DNA database exists and isn't updated!") ||
        readline("... Reset index? [y]/n: ") == "n" ) ) ||
         readline("... DNA database exists. Reset index? y/[n]: ") != "y" ))
      return(TRUE)


  if (!missing(dna.file) && file.exists(dna.file) && grepl(".fa", dna.file)) {
    # update only if different dna.file
    cat.step("Update reference genome by user-defined input.")
    dir_create(ARG_LIST$pasteIN("dna"))
    unlink(ARG_LIST$pasteIN("dna",c("sequences.fa","snp.vcf")))
    link_create(dna.file, ARG_LIST$pasteIN("dna","sequences.fa"), is_linux())
    if (!missing(vcf.file) && file.exists(vcf.file)) {
      link_create(vcf.file, ARG_LIST$pasteIN("dna","snp.vcf"), is_linux())
      if (dir.exists(ARG_LIST$pasteIN(sub("[0-9a-zA-Z]+.vcf$","vcf",vcf.file))))
        link_create(ARG_LIST$pasteIN(sub("[0-9a-zA-Z]+.vcf$","vcf",vcf.file)), ARG_LIST$pasteIN("dna","vcf"), is_linux())
    }
    return(TRUE)
  }

  with(ARG_LIST, {

    dbs.in.files <- sub("\\/genome.fa$","",list.files(path = folder.in, recursive = T, pattern = "genome.fa"))
    cat.step("These reference genomes are available:")
    cat.step("\t0 : no DNA database")
    if (length(dbs.in.files)>0) cat(paste0("\t",1:length(dbs.in.files)," : ",dbs.in.files,"\n"))
    cat.step(paste0("\t",length(dbs.in.files)+1," : new DNA database"))
    tmp <- as.numeric(readline("... Choose a number: "))
    if (is.na(tmp) || (tmp >= 0) && dir_exists(pasteIN("dna"))) unlink(pasteIN("dna"), recursive = T)
    if (tmp == 0) return(FALSE)
    dir_create(pasteIN("dna"))
    if (tmp == length(dbs.in.files)+1) {
      cat.step("Choose a DNA sequence file (unzipped fasta)!")
      tmp <- file.choose()
      if (! grepl("\\.(fa|fasta)$", tmp)) stop("Wrong user file selected")
      link_create(tmp, f.seq, is_linux())
      cat.step("Choose optional a SNP file (unzipped vcf) with chromosome names matching the genome file!")
      tryCatch(tmp <- file.choose(), error = function(e) { tmp <<- 0 })
      if (grepl("\\.vcf$", tmp)) link_create(tmp, pasteIN("dna", "snp.vcf"), is_linux())
      prepare.vcf(ARG_LIST, "dna")
    } else if (tmp <= length(dbs.in.files) && tmp > 0) {
      link_create(pasteIN(dbs.in.files[tmp],"genome.fa"), f.seq, is_linux())
      if (dir.exists(pasteIN(dbs.in.files[tmp],"vcf"))) {
        if (is_linux()) {
          link_create(pasteIN(dbs.in.files[tmp],"vcf"), pasteIN("dna","vcf"), is_linux())
        } else  dir_copy(pasteIN(dbs.in.files[tmp],"vcf"), pasteIN("dna","vcf"))
      }
    }

    return(TRUE)
  })

}
