db_annotate <- function(ARG_LIST, db.name) {

  cat.step("... create annotation")

  with(ARG_LIST, {

    if (DB.TABLE$anno[db.name==DB.TABLE$db]) {

      gff.table <- read.delim(pasteIN(db.name, "annotation.gff3"), header=FALSE, comment.char = "#", stringsAsFactors = FALSE,
                              col.names = c("seqname","source","feature", "start", "end", "score", "strand", "frame", "attribute"))
      gff.table$attribute <- gsub("%3B", ";", gsub("%2C",",",gff.table$attribute)) # Sonderzeichen durch Kommas bzw. Semikolon ersetzen
      gff.table.ok <- gff.table[gff.table$feature=="exon", c(1,4,5,7,9)]

      is.ensembl <- grepl("ensembl", db.name) #grepl("transcript:ENS[A-Z]{3}T", gff.table.ok$attribute[1])
 #     is.refseq  <- TRUE
      if (is.ensembl)  {
        gff.table.ok$genbank <- sub(".*transcript:(ENS[A-Z]+[0-9]+);.*","\\1",gff.table.ok$attribute)
        gff.table.ok$gff_id  <- gff.table.ok$genbank
      } else {
        gff.table.ok$genbank <- NA
        pos <- grep("Genbank:", gff.table.ok$attribute)
        gff.table.ok$genbank[pos] <- sub(".*Genbank:([A-Z_0-9.]+)[;,].*", "\\1", gff.table.ok$attribute[pos])
        gff.table.ok$gff_id  <- sub(".*;Parent=([a-z0-9]+);.*", "\\1", gff.table.ok$attribute)
      }

      rna.exon.list <- split(gff.table.ok[,2:4], gff.table.ok$gff_id)
      rna.exon.list <- lapply(rna.exon.list, function(X) X[order(X$start, decreasing = (X$strand[1]=="-")),1:2]) # bring exons in correct order!
      rna.exon.list2 <- lapply(rna.exon.list, function(x) cumsum(x$end-x$start+1))

      rna.anno <- gff.table.ok[!duplicated(gff.table.ok$gff_id),] # only first line of each rna
      if (is.ensembl) {
        gff.table.transcript <- gff.table[grep("transcript|[a-z_]+RNA$|pseudogene|gene_segment",gff.table$feature), 9]
        gff.table.transcript <- grep("transcript:", gff.table.transcript, value=T) # grep important for pseudogene
        transcript.ids <- sub(".*transcript:([A-Z0-9]+);.*", "\\1", gff.table.transcript)
        pos <- match(rna.anno$gff_id, transcript.ids)
        gff.table.transcript.ok <- gff.table.transcript[pos]
        rna.anno$genbank <- sub(".*gene:(ENS[A-Z0-9]+);.*", "\\1", gff.table.transcript.ok)
        rna.anno$symbol  <- sub(".*gene:([A-Z0-9]+);.*","\\1",sub(".*Name=([A-Z0-9-]+);.*","\\1", gff.table.transcript.ok))

        gff.table.gene <- gff.table[grep("^(RNA|.*gene)$", gff.table$feature), 9]
        pos <- match(rna.anno$genbank, sub(".*gene:([A-Z0-9]+);.*", "\\1", gff.table.gene))
        gff.table.gene.ok <- gff.table.gene[pos]
        rna.anno$name <- sub(".*gene:([A-Z0-9]+);.*","\\1",sub(".*;description=(.+);gene_id=.*", "\\1", gff.table.gene.ok))
        rna.anno$name <- sub(" +\\[.*","", rna.anno$name) # remove source information
      } else {
        rna.anno$symbol  <- sub(".*;gene=([/().A-Za-z0-9_-]+);?.*","\\1",rna.anno$attribute)
        rna.anno$name    <- sub(".*;product=([/().A-Za-z0-9_, %+'-]+);?.*","\\1",rna.anno$attribute)
      }
      rna.anno$variant <- NA
      if (is.ensembl) {
        pos <- grep("-[0-9]+$", rna.anno$symbol)
        rna.anno$variant[pos] <- gsub(".*-", "", rna.anno$symbol[pos])
        rna.anno$symbol[pos] <- sub("-[0-9]+$","",rna.anno$symbol[pos])
        rna.anno$type    <- sub(".*;biotype=([a-z_A-Z]+);.*", "\\1", gff.table.gene.ok)
      } else {
        pos <- grep("transcript variant", rna.anno$name)
        rna.anno$variant[pos] <- sub(".* transcript variant ([.A-Za-z0-9]+).*", "\\1", rna.anno$name[pos])
        rna.anno$type    <- sub(".*;gbkey=([a-z_A-Z]+);.*", "\\1", rna.anno$attribute)

        # needed for pseudogenes ... and maybe others!
        pos <- which(rna.anno$type=="exon")
        tmp <- grep("gene_biotype=", gff.table$attribute, value=T)
        tmp <- do.call(rbind, strsplit(sub("ID=([a-z0-9]+);.*;gene_biotype=([a-zA-Z_]+).*","\\1 \\2", tmp), " "))
        rna.anno$type[pos] <- tmp[match(rna.anno$gff_id[pos], tmp[,1]),2]
        rna.anno$name[pos] <- sub(".*;gene=([A-Z0-9]+);.*","\\1",rna.anno$name[pos])
      }
      rna.anno$n_exon <- NA
      pos <- match(rna.anno$gff_id, names(rna.exon.list))
      rna.anno$n_exon[!is.na(pos)] <- sapply(rna.exon.list, nrow)[pos[!is.na(pos)]]
      rna.anno$chr <- rna.anno$attribute <- NULL
      names(rna.anno)[1] <- "chr"
      rna.anno$short <- paste0(rna.anno$symbol, "(", rna.anno$type, ")")
      pos <- !is.na(rna.anno$variant)
      rna.anno$short[pos] <- paste0(sub(")","", rna.anno$short[pos]), " ", rna.anno$variant[pos], ")")
      # filter
      if (!is.ensembl) {
        tmp <- rna.anno$genbank[duplicated(rna.anno$genbank) & !is.na(rna.anno$genbank)]
        write.csv(rna.anno[!is.na(match(rna.anno$genbank, tmp)),], pasteOUT(db.name,"rna-genbank-duplicates.csv"), row.names = FALSE)
        write.csv(rna.anno[is.na(rna.anno$genbank),], pasteOUT(db.name,"rna-genbank-missing.csv"), row.names = FALSE)
      }

      # no exonic information ... only rna.anno information
      if (length((pos=which(is.na(match(rna.anno$gff_id, names(rna.exon.list))))))>0) {
        warning(length(pos)," sequences without gff information.")
        write.csv(rna.anno[pos,], pasteOUT(db.name,"no-exonic-info.csv"), row.names=FALSE)
      }
      # only exonic information ... no rna.anno information
      if (length((pos=which(is.na(match(names(rna.exon.list), rna.anno$gff_id)))))>0) {
        warning(length(pos)," gff definitions without fasta sequence.")
        writeLines(names(rna.exon.list)[pos],  pasteOUT(db.name,"no-fasta.txt"))
      }
      save(list=c("rna.anno", "rna.exon.list", "rna.exon.list2"), file=pasteTMP(db.name,"annotation.Rda"))

    } else {
      lines <- grep("^>",readLines(pasteIN(db.name, "sequences.fa")), value=T)
      rna.anno <- as_tibble(do.call(rbind, strsplit(sub("^>([.|0-9A-Za-z_-]+) (.*)$","\\1 ___ \\2",lines), " ___ ")))
      if (length(rna.anno)==1) rna.anno[2] <- NA
      names(rna.anno) <- c("id", "name")
      rna.anno$id <- sub("\\..*", "", rna.anno$id)
      rna.anno$name <- chartr("()","[]",gsub(",", ".", rna.anno$name))
      save(rna.anno, file=pasteTMP(db.name,"annotation.Rda"))
    }
  })
}
