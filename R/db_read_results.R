db_read_results <- function(ARG_LIST, DB_LIST, probes.unchecked, cluster) {

#  save(DB_LIST, probes.unchecked, file="db_read_results.Rda")
  if (DB_LIST$RNA_EXISTS) cat.step("... read mapping results and add genomic/exonic position") else
    cat.step("... read mapping results")

  with(ARG_LIST, with(DB_LIST, {

#    if (!is_linux()) {
      mapping.table.raw <- as.data.frame(fread(text = grep("XM:i:0",readLines(pasteIN(db.name,"bowtie.sam")), invert=T, value=T), header = F))
#    } else {
#      mapping.table.raw <- read.delim(pasteIN(db.name,"bowtie.sam"), header=FALSE, stringsAsFactors=FALSE)
#    }
    names(mapping.table.raw) <- c("qname", "flag", "rname", "pos", "mapq", "cigar", "rnext",
                                  "pnext", "tlen", "seq", "qual", "etc1", "etc2", "nmismatch")
    mapping.table.raw <- mapping.table.raw[mapping.table.raw$flag != 4,]
    save(mapping.table.raw, file=pasteTMP(db.name,"mapping.table.raw.Rda"))

    mapping.table <- mapping.table.raw[,c("qname","rname","pos","nmismatch","etc2","cigar","flag")] %>%
      as.tbl %>% filter(qname %in% probes.unchecked) %>%
      mutate(probeset  = sub("___[0-9]+","", qname), nmismatch = as.numeric(sub("NM:i:","", nmismatch)),
             rname     = sub("\\.[0-9]+$","",rname), cigar     = as.integer(sub("M","",cigar)),
             flag      = c("+","-")[(flag==16)+1]) %>% select(1,2,8,3:7)
      #select(qname, rname, probeset, everything())
    names(mapping.table) <- c("probe", "rna.id", "probeset", "rna.pos", "nmismatch", "rna.mismatch", "len", "strand")

    if (RNA_EXISTS && DNA_EXISTS && ANNO_EXISTS) {

#      cat.step("... add genomic position")

      f <- function(X) {
        id <- X$rna.id[1]
        pos <- which(DB_LIST$rna.anno$gff_id==id)
        POS <- unlist(lapply(1:nrow(DB_LIST$rna.exon.list[[id]]), A=DB_LIST$rna.exon.list[[id]], function(i,A) seq(A[i,1],A[i,2])))
        if (DB_LIST$rna.anno$strand[pos]=="-") POS <- rev(POS)
        X$chr.name  <- DB_LIST$rna.anno$chr[pos]
        X$rna.chr.exon  <- apply(cbind(vapply(X$rna.pos, function(x) which(x<=DB_LIST$rna.exon.list2[[id]])[1], numeric(1)),
                                       vapply(X$rna.pos+X$len-1, function(x) which(x<=DB_LIST$rna.exon.list2[[id]])[1], numeric(1))), 1,
                                 function(x) paste0("E",unique(x), collapse="-"))
        X$chr.pos.start <- POS[X$rna.pos]
        X$chr.pos.stop  <- POS[X$rna.pos+X$len-1]
        return(X)
      }

      if (NCPU == 1) {
        mapping.table %<>% group_by(rna.id) %>% do(f(.)) %>% ungroup
      } else {
        cluster %>%  cluster_assign_value('f', f) %>%
          cluster_assign_value('rna.exon.list', rna.exon.list) %>%
          cluster_assign_value('rna.exon.list2', rna.exon.list2) %>%
          cluster_assign_value('rna.anno', rna.anno)
        mapping.table %<>% partition(rna.id, cluster=cluster) %>% do(f(.)) %>% collect %>% ungroup
        cluster %>% cluster_ls() %>% unique %>% unlist %>% cluster_rm(cluster=cluster)
      }
    } else if (!RNA_EXISTS) { # dna
      mapping.table <- mapping.table %>%
        mutate(chr.name = rna.id,
               chr.pos.start = ifelse(strand=="+", rna.pos, rna.pos + len - 1),
               chr.pos.stop  = ifelse(strand=="+", rna.pos + len - 1, rna.pos))
    } else { # rna without anno
      mapping.table <- mapping.table %>% mutate(chr.name = NA, chr.pos.start = NA, chr.pos.stop = NA)
    }

    save(mapping.table, file=pasteTMP(db.name,"mapping.table.step1.Rda"))

  }))

}
