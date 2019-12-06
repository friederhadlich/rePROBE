db_create_probe_table <- function(ARG_LIST, DB_LIST, probes.unchecked, cluster) {

  cat.step("... create probe table")

  save(DB_LIST, probes.unchecked, file="db_create_probe_table.Rda")
  with(ARG_LIST, with(DB_LIST, {

    load(pasteTMP(db.name,"mapping.table.step4.Rda"))
    load(pasteTMP("probes.Rda"))

    probes <- probes %>% as.tbl %>% filter(id %in% probes.unchecked)
    f <- function(X)  data.frame(nmismatch=X$nmismatch[1], nresult=nrow(X), result.id=paste(X$rna.id,collapse=";"),
                                 result.name=paste(X$rna.name,collapse=";"), result.pos=paste(ff(X),collapse=";"),
                                 nSNP.min=min(X$nSNP), nSNP.max=max(X$nSNP), stringsAsFactors = FALSE)
    ff <- function(X) {
      if("rna.chr.exon" %in% names(X)) { # (RNA_EXISTS && ANNO_EXISTS && DNA_EXISTS)
        paste0(X$chr.name, ":", X$chr.pos.start,"-", X$chr.pos.stop, "(", X$rna.chr.exon,")")
      } else if (!is.na(X$chr.name[1])) { # (!RNA_EXISTS) {
        paste0(X$chr.name, ":", X$chr.pos.start,"-", X$chr.pos.stop)
      } else return(NA)
    }

    probe.table <- data.frame(probe.id=probes$id, probe.seq=probes$seq, nmismatch=NA, nresult=NA, result.id=NA, result.name=NA,
                              result.pos=NA, nSNP.min=NA, nSNP.max=NA, stringsAsFactors=FALSE)
    pos <- mapping.table.probe$probe %>% match(unique(.[duplicated(.)]))
    tmp <- mapping.table.probe[is.na(pos),]  # single hit probes
    probe.table[match(tmp$probe,probes$id),-(1:2)] <- tmp[c("nmismatch", "nmismatch", "rna.id", "rna.name", "chr.name", "nSNP", "nSNP")]
    probe.table$nresult[!is.na(probe.table$nresult)] <- 1
    probe.table$result.pos[match(tmp$probe,probes$id)] <- ff(tmp)
    if (NCPU == 1) {
      tmp <- mapping.table.probe[!is.na(pos),] %>% group_by(probe) %>% do(f(.))
    } else {
      cluster %>%
        cluster_assign_value('f', f) %>%
        cluster_assign_value('ff', ff) # %>% cluster_assign_value('RNA_EXISTS', RNA_EXISTS)
      tmp <- mapping.table.probe[!is.na(pos),] %>% partition(probe, cluster=cluster) %>% do(f(.)) %>% collect()
      cluster %>% cluster_ls() %>% unique %>% unlist %>% cluster_rm(cluster=cluster)
    }

    probe.table[match(tmp$probe,probes$id),-(1:2)] <- tmp[-1]
    probe.table <- probe.table %>% as.tbl
    if (!RNA_EXISTS) probe.table$result.id <- probe.table$result.name <- NA
    if (!DNA_EXISTS) probe.table$result.pos <- NA

    write.csv(probe.table, file=pasteOUT(db.name,"probe.table.csv"), row.names = FALSE)
    save(probe.table, file=pasteTMP(db.name,"probe.table.Rda"))

  }))

}
