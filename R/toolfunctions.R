cat.step <- function(string) catn(" ... ", string, sep="");

pastef <- function(...) paste(..., sep=.Platform$file.sep)

is_linux <- function() Sys.info()[1] == "Linux"
catn <- function(...) { cat(...,"\n"); return(T) }

choose_directory = function(caption = 'Select data directory') {
  if (!is_linux() && requireNamespace("utils", quietly = TRUE)) { # windows
    gsub("\\\\", "/", utils::choose.dir(caption = caption))
  } else { # linux
    if (requireNamespace("tcltk", quietly = TRUE))
    tcltk::tk_choose.dir(caption = caption)
  }
}

get.fasta <- function(EXONS, genome.file="genome.fa") {
  names(EXONS)[c(1,4,5,7)] <- c("chr", "start", "end", "strand")
  EXONS$id <- sub(".*Parent=([0-9a-zA-Z:._]+);.*", "\\1", as.data.frame(EXONS)[,9])

  X <- split(EXONS, EXONS$chr)
  Y <- get_fasta(genome.file, X)
  Z <- do.call(rbind, lapply(names(Y), function(i) cbind(X[[i]], seq=Y[[i]])))
  pos <- Z$strand=="-"
  Z <- as.data.frame(Z)
  Z$seq[pos] <- chartr("ACGTacgt", "TGCAtgca", sapply(strsplit(Z$seq[pos],""), function(x) paste(rev(x), collapse="")))
  #Z[strand=="-", seq:=chartr("ACGTacgt", "TGCAtgca", sapply(strsplit(seq,""), function(x) paste(rev(x), collapse="")))]
  #Z$id <- sub(".*Parent=([0-9a-zA-Z._]+);.*", "\\1", Z[,9])
  #Z[, c("id","rank"):=list(sub(".*Parent=([0-9a-zA-Z._]+);.*","\\1",V9), as.numeric(sub(".*;rank=(.*);version.*","\\1",V9)))]
  if (grepl(".*;rank=.*;version.*", Z[1,9])) {
    Z$rank <- as.numeric(sub(".*;rank=(.*);version.*", "\\1", Z[,9]))
    Z <- Z[order(Z$id,Z$rank), c("chr", "start", "end", "strand", "seq", "id", "rank")]
  } else {
    Z <- Z[order(Z$id), c("chr", "start", "end", "strand", "seq", "id")]
  }
  Z <- as.data.table(Z)
  ZZ <- as.data.frame(Z[, paste(seq, collapse=""), by="id"])
  lines <- rep(paste0(">", sub("transcript:","",ZZ$id)), each=2)
  lines[seq(2, length(lines), by=2)] <- ZZ$V1
  #writeLines(lines, fasta.file)
  return(invisible(lines))
}

f.list.to.df <- function(my.list, d=1e4) {
  res <- data.frame(stringsAsFactors = FALSE)
  for (i in 1:ceiling(length(my.list)/d)) {
    cat(i,"of",ceiling(length(my.list)/d),"\n")
    res <- rbind(res, do.call(rbind, my.list[((i-1)*d+1):min(i*d, length(my.list))]))
  }
  return(res)
}
