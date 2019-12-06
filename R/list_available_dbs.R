#' List all "installed" databases
#' @param ARG_LIST list with user settings created by prepare_input_arguments.R
#'
#' @export
list_available_dbs <- function(ARG_LIST) {

  if (!dir.exists(ARG_LIST$folder.in))  stop(paste("Folder", ARG_LIST$folder.in, "does not exist!"))

  dbs.in.files <- c(list.files(ARG_LIST$folder.in, recursive = T, pattern = "(annotation.gff3|genome.fa|sequences.fa)"),
                    sub(paste0(ARG_LIST$folder.in,"/"), "", grep("vcf$",list.dirs(ARG_LIST$folder.in, recursive = T), value=T)))
  if (length(dbs.in.files)==0) { message("No database found."); return(NULL) }
  TMP <- which(is_link(ARG_LIST$pasteIN(dbs.in.files)))
  if (length(TMP)>0)  TMP <- TMP[!is_file(link_path(ARG_LIST$pasteIN(dbs.in.files[TMP])))]
  dbs.in.files[dbs.in.files == "dna/sequences.fa"] <- "dna/genome.fa"

  if (!is_linux() && length(TMP)>0) dbs.in.files <- dbs.in.files[-TMP]
#  dbs.in <- unique(sub("\\/.*", "", dbs.in.files))
  db.content.ref <- c("sequences.fa","annotation.gff3","genome.fa","vcf")
  names(db.content.ref) <- c("rna","anno","dna","snp")
  TMP <- data.table(db = sub("\\/.*", "", dbs.in.files),
                    type = factor(names(db.content.ref)[match(sub(".*\\/", "", dbs.in.files),db.content.ref)], levels = names(db.content.ref)))
  TMP2 <- as.data.frame(dcast(TMP, ... ~ type, value.var = "type"))
  TMP3 <- cbind(TMP2[1], rna=F, anno=F, dna=F, snp=F)
  TMP2[-1] <- !is.na(TMP2[-1])
  TMP3[,names(TMP3) %in% names(TMP2)] <- TMP2
##  TMP2[names(db.content.ref)] <- lapply(TMP2[-1], function(x) !is.na(x))#lapply(names(db.content.ref), function(x) !is.na(x))]
  # for (db in dbs.in) {
  #   cat("\t-", db)
  #   if (db == "dna") {
  #     #dna.ref.file <- paste(TMP, "in", db, "ref.txt", sep=.Platform$file.sep)
  #     #if (file.exists(dna.ref.file)) cat("->", readLines(dna.ref.file)[1], sep="")
  #     cat("->", link_path(pastef(folder.in, "dna", "genome.fa")), sep="")
  #   }
  #   cat(paste0(" (", paste(collapse="+", names(db.content.ref)[db.content.ref %in% sub(".*\\/", "", grep(db, dbs.in.files, value=T))]), ")\n"))
  # }
  if (dir.exists(ARG_LIST$pasteIN("dna","vcf"))) TMP3$snp[TMP3$db=="dna"] <- TRUE
  rownames(TMP3) <- paste0("DB",1:nrow(TMP3),":")
  return(TMP3)
}
