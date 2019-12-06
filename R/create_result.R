create_result <- function(ARG_LIST, keep_unqualified_probesets, probes, sep1, sep3) {

  save(probes, sep1, sep3, file="create_result.Rda")
  cat.step("create final probeset table") # Create final probeset table ...

  PS.final <- probes %>% filter((.data$rank-0.1) == round(.data$rank-0.1)) %>% group_by(.data$ps.id) %>%
    summarise(rank.best=(function(x) {tmp=table(x); sort(names(tmp)[tmp==max(tmp)])[1]})(rank)) %>%
    mutate(id = paste(.data$ps.id, .data$rank.best))

  with(ARG_LIST, {

#    if (!all(file.exists(c(pasteTMP("probeset.Rda"), pasteOUT("probeset.csv"))))) {

      # add PS tables from rounds!
      probeset.tables <- bind_rows(
        lapply(unique(DB.TABLE$rank), function(r) {load(pasteTMP(paste0("probeset_rank",r,".Rda"))); PS %>% mutate(rank.best=r+0.1)})) %>%
        mutate(id = paste(.data$probe.set, .data$rank.best))
      # find winner PS
      PS <- merge(PS.final, probeset.tables, by = "id") %>% select(.data$probe.set:.data$res.probes) # winner PS funzt nicht
      # update n.probe, n.probe.sh, ...
      PS.final <- probeset.tables %>% group_by(.data$probe.set) %>%
        summarise(n.probe=as.integer(max(.data$n.probe)), n.probe.sh=sum(.data$n.probe.sh), n.probe.mh=sum(.data$n.probe.mh)) %>% # compute n.probe, n.probe.sh, n.probe.mh
        mutate(n.probe.nh=.data$n.probe-.data$n.probe.sh-.data$n.probe.mh) # add n.probe.nh
      PS.final <- merge(PS.final, PS[-(2:5)], by="probe.set", all.x=T)  %>% as.tbl
      # update gene.alternative
      tmp <- probeset.tables %>% select(.data$probe.set:.data$res.probes) %>%
        bind_rows(PS) %>% group_by(.data$probe.set) %>%
        do((function(X) {Y=X[X$res.mapping==max(X$res.mapping),]; setdiff(Y,Y[duplicated(Y),])})(.data)) %>% # find alternatives
        summarise(alternatives=paste(c(.data$res.name,.data$res.alternative)[!is.na(c(.data$res.name,.data$res.alternative))], collapse=sep3)) %>%
        filter(.data$alternatives!="" & .data$alternatives!="NA")
      pos <- match(tmp$probe.set, PS.final$probe.set)
      PS.final$res.alternative[pos] <- sub(paste0("NA",sep3),"",paste(PS.final$res.alternative[pos], tmp$alternatives, sep=sep3), fixed = T)
      PS.final$res.symbol <- sapply(strsplit(PS.final$res.name, sep1, fixed = T), function(x) sub("\\(.*","",c(grep("^$|^LOC|^ENS[A-Z]{4}[0-9]+",x, value=TRUE, invert=TRUE), x)[1]))

      for (db.name in DB.TABLE$db) {
        if (!file.exists(pasteIN(db.name, "annotation.gff3"))) next
        id <- paste0("res.id.",db.name)
        PS.final[[id]] <- NA
        pos <- grep(db.name, PS.final$res.db)
        pos2 <- sapply(strsplit(PS.final$res.db[pos], split = sep1, fixed=T), function(x) which(x==db.name))
        rna.anno <- NULL
        load(pasteTMP(db.name,"annotation.Rda"))
        PS.final[pos, id] <- rna.anno$genbank[match(mapply(FUN=function(x,y) x[y], strsplit(PS.final$res.name[pos], split = sep1, fixed=T), pos2),rna.anno$short)]
      }

      if (file.exists(pasteTMP("probesubsetting.csv"))) {
        pss.table <- read.csv(pasteTMP("probesubsetting.csv"), stringsAsFactors = F)
        tmp <- sapply(split(pss.table$pss.name, pss.table$ps.id), paste, collapse=",")
        PS.final$probe.subsets <- tmp[match(PS.final$probe.set, names(tmp))]
      } else
        PS.final$probe.subsets <- NA
      write.csv(PS.final, pasteOUT("probeset.csv"), row.names = FALSE, na = "")
      save(PS.final, file=pasteTMP("probeset.Rda"))
      catn("")

#    } else load(pasteTMP("probeset.Rda"))

    #####################################################################################

#    if (!file.exists(pasteOUT("probes.csv"))) {
      probe.table.final <- data.frame()
      for (rank in sort(unique(DB.TABLE$rank))) {
        probe.table <- read.csv(pasteOUT(paste0("probes_rank",rank,".csv")), stringsAsFactors = F)
        probe.table.final <- rbind(probe.table.final, probe.table %>% filter(probe.id %in% probes$id[rank==floor(abs(probes$rank))]))
      }
      probe.table.final <- rbind(probe.table.final, probe.table %>% filter(probe.id %in% probes$id[0==probes$rank]))
      write.csv(probe.table.final, pasteOUT("probes.csv"), row.names = F)
#    }

    #####################################################################################

    load(pasteTMP("init.Rda"))

    if (is.na(CHIP_TYPE)) stop("Missing chip type!")
    cat.step(paste0("create file ", CHIP_TYPE,".Rda"))

    # used probes
    pos <- which(!is.na(PS.final$res.probes))
    tmp <- strsplit(PS.final$res.probes[pos],","); names(tmp) <- PS.final$probe.set[pos]
    tmp <- data.frame(ps = rep(names(tmp), times=sapply(tmp, length)), pm = unlist(tmp), stringsAsFactors = F)
    tmp$mm <- probe.table$mm[match(tmp$pm, probe.table$probe_id)]
    # control probes
    if(keep_unqualified_probesets)
      tmp2 <- probe.table[!probe.table$probeset_id %in% PS.final$probe.set[pos], c("probe_id", "probeset_id", "mm")]
    else
      tmp2 <- probe.table[!probe.table$probeset_id %in% PS.final$probe.set, c("probe_id", "probeset_id", "mm")]
    names(tmp2) <- c("pm", "ps", "mm")
    tmp12 <- rbind(tmp, tmp2[c(2,1,3)])
    tmp12 <- lapply(split(tmp12[-1], tmp12[[1]]), data.matrix)
    D <- list2env(tmp12)
    # tmp2 <- strsplit(PS.final$res.probes[-pos],",");  names(tmp2) <- PS.final$probe.set[-pos]
    # tmp2 <- data.frame(ps = rep(names(tmp2), times=sapply(tmp2, length)), p = unlist(tmp2))
    # # removed probes
    # #tmp[["removed_probes"]] <- probe.table %>% as.tbl %>% filter(! .data$probe_id %in% unlist(tmp)) %>%
    # #  filter(.data$probeset_id %in% PS.final$probe.set) %>% select(probe_id) %>% unlist %>% as.character
    # tmp <- lapply(tmp, function(x) as.matrix(data.frame(pm=as.numeric(x), mm=NA)))
    # # control probes
    # if (keep_unqualified_probesets)
    #   tmp2 <- probe.table[!probe.table$probeset_id %in% PS.final$probe.set[!is.na(PS.final$res.probes)], c("probe_id", "type", "probeset_id", "mm")]
    # tmp2$probe_id <- as.numeric(tmp2$probe_id)
    # tmp2$type <- as.character(tmp2$type)
    # tmp2 <- lapply(split(tmp2[1:2], tmp2$probeset_id), function(X) {
    #   pos.pm <- grep("^pm",X$type);
    #   if (length(pos.pm)>0) RES <- as.matrix(data.frame(pm=X$probe_id[pos.pm], mm=NA))
    #   pos.mm <- setdiff(1:nrow(X),pos.pm)
    #   pos.mm <- pos.mm[!is.na(match(pos.mm-1, pos.pm))]
    #   if (length(pos.mm)>0) RES[match(pos.mm-1, pos.pm),2] <- X$probe_id[pos.mm]
    #   else return(NULL)
    #   ## TODO: check output format ... maybe use AffyGenePDInfo
    #   RES
    # })
    # table(sapply(tmp2, is.null))

    # D <- list2env(c(tmp, tmp2))
    # create ps.type table to remove controls later on ...
    ps.type <- data.frame(id = names(D),
                          type = as.character(probeset.table$pss.type[match(names(D), probeset.table$ps.id)]),
                          stringsAsFactors = F)
    #ps.type$type[ps.type$id=="removed_probes"] <- "control->removed_probes"
    ps.type$category <- c("main","control")[grepl("^control->|^normgene->|^[rR]eporter|^rRNA$|^Introns$", ps.type$type)+1]
    eval(parse(text=paste0("`",CHIP_TYPE,"` <- D; save('",CHIP_TYPE,"', ps.type, file=pasteOUT('",CHIP_TYPE,".Rda'))")))
    eval(parse(text=paste0("assign(\"",CHIP_TYPE,"\", D, envir = globalenv())")))
  })

}
