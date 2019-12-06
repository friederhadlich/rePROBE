db_merge_results <- function(ARG_LIST, DB_LIST, cluster) {

  cat.step("... merge identical mapping results of single probe")

  with(ARG_LIST, with(DB_LIST, {

    load(pasteTMP(db.name,"mapping.table.step1.Rda"))

    if (RNA_EXISTS && DNA_EXISTS && ANNO_EXISTS) {
      mapping.table$rna.name <- as.character(rna.anno$short[match(mapping.table$rna.id, rna.anno$gff_id)])
      qname_pos <- paste(mapping.table$probe, mapping.table$chr.name, mapping.table$chr.pos.start, mapping.table$chr.pos.stop)
      pos <- match(qname_pos, unique(qname_pos[duplicated(qname_pos)]))

      f <- function(X) {
        res <- X[1,]
        col <- c("rna.id", "rna.pos", "rna.chr.exon")
        res[1,col] <- sapply(X[col], function(x) paste(unique(x), collapse=","))
        #res[1,c(2,4,9)] <- sapply(X[c(2,4,9)], function(x) paste(unique(x), collapse=","))
        tmp <- unique(sub("\\(.*","",X$rna.name))
        pos <- grep(tmp[1],X$rna.name)
        res$rna.name <- paste0(tmp[1], "(", paste(gsub(".*\\(|\\)","",X$rna.name[pos]), collapse=","), ")")
        if (length(tmp)>1) {
          for (tm in tmp[-1]) {
            pos <- grep(tm,X$rna.name)
            res$rna.name <- paste0(res$rna.name, ",", tm, "(", paste(gsub(".*\\(|\\)","",X$rna.name[pos]), collapse=","), ")")
          }
          #    cat(res$probe, ":", tmp,"\n")
        }
        return(res)
      }

      mapping.table.probe <- mapping.table[!is.na(pos),]
      mapping.table.probe$tmp <- qname_pos[!is.na(pos)]
      if (NCPU == 1) {
        mapping.table.probe %<>% group_by(tmp) %>% do(f(.data)) %>% ungroup
      } else {
        cluster %>%  cluster_assign_value('f', f)
        #start <- proc.time()
        mapping.table.probe %<>% partition(tmp, cluster=cluster) %>% do(f(.data)) %>% collect %>% ungroup
        #time_elapsed_parallel <- proc.time() - start # End clock
        cluster %>% cluster_ls() %>% unique %>% unlist %>% cluster_rm(cluster=cluster)
      }
      mapping.table.probe %<>% select(-tmp) %>% rbind(mapping.table[is.na(pos),])
    } else if (!RNA_EXISTS) { # dna
      mapping.table.probe <- mapping.table %>% mutate(rna.name = NA)
    } else { # rna without anno
      pos <- match(mapping.table$rna.id, rna.anno$id)
      #    if (any(is.na(pos))) pos <- match(mapping.table$rna.id, sub("\\..*", "", rna.anno$id))
      mapping.table.probe <- mapping.table
      if (any(is.na(pos))) {
        mapping.table.probe$rna.name <- mapping.table$rna.id
      } else {
        mapping.table.probe$rna.name <- rna.anno$name[pos]
      }
    }

    save(mapping.table.probe, file=pasteTMP(db.name,"mapping.table.step2.Rda"))

  }))

}
