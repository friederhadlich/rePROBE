db_filter_results <- function(ARG_LIST, DB_LIST) {

  cat.step(paste("...", c("relaxed","strict")[ARG_LIST$STRICT_MAPPING+1], "filtration of mapping results"))

  with(ARG_LIST, with(DB_LIST, {

    load(pasteTMP(db.name,"mapping.table.step3.Rda"))
    mapping.table.probe <- mapping.table.probe %>% select(-one_of(c("seq","probeset"))) %>% ungroup

    if (!RNA_EXISTS) {
      mapping.table.probe$rna.id[nchar(mapping.table.probe$rna.id)<=5] <- 1 # permit only single hit in full genome
    }
    if (STRICT_MAPPING) {
      ## - remove mismatching mapping results ( keep only 0-mismatch-results)
      mapping.table.probe <- mapping.table.probe %>% filter(nmismatch==0)
      ## - remove probes multimapping on 1 rna
      probes.duplicated <- mapping.table.probe %>% .[duplicated(.[1:2]),"probe"] %>% unique %>% unlist
      mapping.table.probe <- mapping.table.probe %>% filter(! probe %in% probes.duplicated)
    } else {
      ## - remove all but best-strata mapping results
      pos <- mapping.table.probe$probe %>% match(unique(.[duplicated(.)]))
      mapping.table.probe %<>% .[!is.na(pos),] %>% group_by(probe) %>%
        do((function(X) X[X$nmismatch==min(X$nmismatch),])(.)) %>% ungroup %>% rbind(mapping.table.probe[is.na(pos),])
      catn("")
    }
    write.csv(mapping.table.probe, file=pasteOUT(db.name,"mapping.table.csv"), row.names = FALSE)
    save(mapping.table.probe, file=pasteTMP(db.name,"mapping.table.step4.Rda"))

  }))

}
