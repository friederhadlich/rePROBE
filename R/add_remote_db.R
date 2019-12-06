#' Add a RNA database from a remote source into this project
#'
#' @param ARG_LIST list with user settings created by prepare_input_arguments()
#' @param db.type Database type, either "refseq" or "ensembl" (default)
#' @param db.version Release number
#' @param update Renew existing data
#' @import data.table fs
#' @importFrom Rbowtie bowtie_build bowtie
#' @importFrom curl curl_download curl
#' @importFrom R.utils gunzip
#' @export
add_remote_db <- function(ARG_LIST, db.type = "ensembl", db.version="", update = F) {

  if (!dir.exists(ARG_LIST$folder.in))  stop("Wrong input for add_remote_db!")

  if (!db.type %in% c("refseq", "ensembl"))  stop("Wrong db.type.type! Choose either refseq and/or ensembl!")

  ################

  download.from.db <- function(type, out.file, in.files, in.path, chr_accessions=NULL) {
    if (file.exists(out.file)) return()
    cat.step(paste("... download", type, "file!"))
    out.path <- sub("[0-9a-zA-Z._-]+$","",out.file)
    for (f in in.files) {
      ff <- paste0(out.path, f)
      if (!file.exists(ff)) curl_download(paste0(in.path, f), ff)
      if (grepl(".gz$", ff)) {
        if (type == "variant calling" || (length(in.files)==1 && file_info(ff)$blocks>1e6)) {
          gunzip(ff)
          file.rename(sub("\\.gz$","",ff), out.file)
        } else {
          con <- gzfile(ff);   content <- readLines(con);   close(con)
          if (is.null(chr_accessions))
            cat(content, file = out.file, sep = "\n", append = T)
          else if (grepl("annotation.gff3$", out.file)){
            tmp2 <- fread(text=grep("^#", content, value = T, invert = T), header = F, sep = "\t")
            #tmp2 <- as.data.table(tmp2)
            pos <- match(tmp2$V1, chr_accessions[,2], nomatch = 0)
            tmp2$V1[pos>0] <- chr_accessions[,1][pos] #tmp2[pos>0, V1:=chr_accessions[,1][pos]]
            fwrite(tmp2, file=out.file, sep = "\t", col.names=F, quote=F, append = T)
          } else {
            tmp2 <- content
            pos <- grep("^>.*ref\\|[0-9._A-Z]+\\|.*", tmp2)
            tmp2[pos] <- sub(".*ref\\|([0-9._A-Z]+)\\|.*","\\1",tmp2[pos])
            if (tmp2[1] %in% chr_accessions[,2]) tmp2[1] <- chr_accessions[,1][tmp2[1]==chr_accessions[,2]]
            tmp2[pos] <- paste0(">", tmp2[pos])
            cat(tmp2, file=out.file, sep = "\n", append = T)
          }
        }
      }
      else cat(readLines(ff), file = out.file, sep = "\n", append = T)
      #      cat("f=",f,"\n")
    }
  }

  ###########

  with(ARG_LIST, {

    tmp <- grep(db.type, list.files(folder.in, pattern = "genome.fa", recursive = T), value=T)
    if (length(tmp)>0 && update) #(!(!ASK_QUESTIONS && catn("\t",db.type," files exist and are reused!", sep="")) &&
         #readline(paste0("\t",db.type," files exist! Reuse? [y]/n: ")) == "n"))
      dir_delete(pasteIN(sub(".genome.fa$","",tmp)))

    if (db.type == "ensembl") {
      DIR <- "ftp://ftp.ensembl.org/pub/"
      con <- curl(DIR);   content <- readLines(con);  close(con)
      DIR <- paste0(DIR, sub(".*release","release", rev(grep(paste0("release-",db.version), content, value=T))[1]), "/")
      SPC <- tolower(SPECIES)
      con <- curl(paste0(DIR,"fasta/"));   content <- readLines(con);  close(con)
      if (! any(grepl(paste0(" ",SPC,"$"), content )) ) stop(paste0("Species '",SPECIES,"' not found!", sep=""))
      RELEASE <- sub(".*release-([0-9]+).*","\\1",DIR)
      GENOME_DIR <- paste0(DIR, "fasta/",SPC,"/dna/")
      GENOME_FILE_PRE <- paste0(SPECIES, ".")
      GENOME_FILE_SUF <- ".dna.toplevel.fa.gz"
      ANNO_DIR <- paste0(DIR, "gff3/",SPC,"/")
      ANNO_FILE_PRE <- paste0(SPECIES, ".")
      ANNO_FILE_SUF <- paste0(".", RELEASE, ".gff3.gz")
      VCF_ID  <- SPC
      VCF_DIR <- paste0(DIR, "variation/vcf/", SPC, "/")
      con <- curl(GENOME_DIR);  tmp <- readLines(con);  close(con)
    } else {
      DIR <- "ftp://ftp.ncbi.nih.gov/genomes/"
      SPC <- SPECIES
      if (SPECIES=="Homo_sapiens" || SPECIES=="Mus_musculus") SPC <- sub("[a-z]+_","_",SPC)
      con <- curl(DIR);  content <- readLines(con);  close(con)
      if (! any(grepl(paste0(" ",SPC,"$"), content )) ) stop(paste0("Species '",SPECIES,"' not found!", sep=""))
      GENOME_DIR <- paste0(DIR, SPC, "/Assembled_chromosomes/seq/")
      con <- curl(GENOME_DIR);  tmp <- readLines(con);  close(con)
      GENOME_FILE_PRE <- sub(".+ ([A-Za-z]+_).*","\\1",tmp[1])
      GENOME_FILE_SUF <- gsub(".*_","_",grep(".fa.gz", tmp, fixed = T, value=T)) #"_chr1.fa.gz"
      ANNO_DIR <- paste0(DIR, SPC, "/GFF/")
      ANNO_FILE_PRE <- ""
      ANNO_FILE_SUF <- "_top_level.gff3.gz"
      VCF_DIR <- "ftp://ftp.ncbi.nih.gov/snp/organisms/archive/"
      if (SPECIES=="Homo_sapiens") VCF_DIR <- "ftp://ftp.ncbi.nih.gov/snp/organisms/"
      if (!STRICT_MAPPING) {
        curl_download(paste0(DIR,SPC,"/README_CURRENT_RELEASE"), pasteIN("README_CURRENT_RELEASE"))
        VCF_ID  <- paste(sub(".*\\t","",readLines(pasteIN("README_CURRENT_RELEASE"))[2:3]), collapse="_")
        VCF_DIR <- paste0(VCF_DIR, VCF_ID, "/VCF/")
        VCF_GENOME <- sub(".*:\t","",grep("ASSEMBLY ACCESSION:", readLines(pasteIN("README_CURRENT_RELEASE")), value=T))
        file_delete(pasteIN("README_CURRENT_RELEASE"))
      }
    }

    db.release <- sub(paste0(".*",GENOME_FILE_PRE,"(.*)",GENOME_FILE_SUF[1]),"\\1",grep(GENOME_FILE_SUF[1], tmp, value=T))[1]
    #for (db.release in sub(paste0(".*",GENOME_FILE_PRE,"(.*)",GENOME_FILE_SUF[1]),"\\1",grep(GENOME_FILE_SUF[1], tmp, value=T))) {
    db.name <- paste0(db.type,"_",db.release)
    if ((ASK_QUESTIONS || !cat.step(paste0("update/download ",db.name,"!"))) &&
        readline(paste0("... ",db.name,": update/download this database? [y]/n :")) == "n")
      return(NULL)
    dir_create(pasteIN(db.name))
    if (db.type == "refseq") {
      con <- curl(sub("seq/","",GENOME_DIR));   content <- readLines(con);   close(con)
      download.from.db("chr_accessions", pasteIN(db.name, "chr_accessions"),
                       sub(".* chr_","chr_",grep("chr_accessions", content, value=T)),
                       sub("seq/","",GENOME_DIR))
      chr_accessions <- read.delim(pasteIN(db.name, "chr_accessions"), stringsAsFactors = F)
    } else chr_accessions <- NULL
    # download genome file
    download.from.db("genome", pasteIN(db.name, "genome.fa"), paste0(GENOME_FILE_PRE, db.release, GENOME_FILE_SUF), GENOME_DIR, chr_accessions)
    # download annotation file
    download.from.db("annotation", pasteIN(db.name, "annotation.gff3"), paste0(ANNO_FILE_PRE, db.release, ANNO_FILE_SUF), ANNO_DIR, chr_accessions)
    # download SNPs
    if (!STRICT_MAPPING) {
      con <- curl(sub(paste0(VCF_ID,".*"),"",VCF_DIR)); content <- gsub(".* ","",readLines(con)); close(con)
      VCF_ID.new <- VCF_ID
      if ((VCF_ID %in% content) || length((VCF_ID.new=grep(paste0("_",sub(".*_", "", VCF_ID), "$"), content, value=T)))>0) {
        VCF_DIR <- sub(VCF_ID, VCF_ID.new, VCF_DIR)
        con <- curl(VCF_DIR);  content <- readLines(con); close(con)
        if (db.type == "ensembl") {
          if(! any(grepl(paste0(SPC,".vcf.gz$"), content)) ) message("\t... no vcf file available")
          download.from.db("variant calling", pasteIN(db.name, "snp.vcf"), paste0(SPC,".vcf.gz"), VCF_DIR)
        } else {
          if (! any(grepl(" 00-All.vcf.gz$", content)) | grepl("^alt_", db.release) ) message("\t... no vcf file available")
          download.from.db("variant calling", pasteIN(db.name, "snp.vcf"), "00-All.vcf.gz", VCF_DIR)
        }
        if (db.type=="ensembl" || (db.type=="refseq") && any(grepl(VCF_GENOME, readLines(pasteIN(db.name, "snp.vcf"),n=5))))
          prepare.vcf(ARG_LIST, db.name)
      }
    }

    if (!file.exists(pasteIN(db.name, "sequences.fa"))) {
      cat.step("... create sequences.fa!")
      f <- tempfile()
      writeLines(fread(pasteIN(db.name, "annotation.gff3"), sep="\n", col.names = "X")[!grep("^#",X)]$X,f)
      X <- fread(f, sep="\t")[V3=="exon"] # fread(pasteIN(db.name, "annotation.gff3"), sep="\t")[V3=="exon"]
      writeLines(get.fasta(X, pasteIN(db.name, "genome.fa")), pasteIN(db.name, "sequences.fa"))
      unlink(f)
    }

    return(db.name)
  })
  ###########

}
