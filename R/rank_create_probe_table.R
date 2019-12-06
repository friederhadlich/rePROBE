rank_create_probe_table <- function(ARG_LIST, rank, probe.tables, sep1, cluster) {

  cat.step("load rank results")

  db.names <- ARG_LIST$DB.TABLE$db[ARG_LIST$DB.TABLE$rank==rank]

  # check genomic positions
  P <- as.tbl(cbind(probe.tables[[1]][1:2], lapply(probe.tables,"[", -(1:3))))
  tmp <- do.call(cbind, lapply(probe.tables, "[", 3)); tmp[is.na(tmp)] <- 10
  if (ARG_LIST$STRICT_MAPPING) {
    P <- P %>% filter(apply(tmp, 1, min)==0) # remove probes without perfect hit
  } else
    P <- P %>% filter(apply(tmp, 1, min)<10) # remove probes without hit

  names(P) <- sub("\\.result\\.", ".", names(P))
  for (db.name in db.names) {
    P[[paste0(db.name,".exon")]] <- sub("[-._A-Z:0-9]+\\(", "", gsub(";[-._A-Z:0-9]+\\(", ";", gsub("\\)","",P[[paste0(db.name,".pos")]])))
    P[grep(":",fixed = T, P[[paste0(db.name,".exon")]]),paste0(db.name,".exon")] <- NA # dna only
    P[[paste0(db.name,".pos")]]  <- sub("\\(.*\\)", "", gsub("\\([-0-9E,]+\\);",";",P[[paste0(db.name,".pos")]]))
  }

  pos.pos <- grep("\\.pos", names(P))
  pos.name <- grep("\\.name", names(P))
  pos.exon <- grep("\\.exon", names(P))
  P <- as.data.frame(P) # otherwise error in access of uninitialized column
  P[paste(rep(c("exon","name","pos"),each=2),c("single","multi"),sep=".")] <- "" # generates warnings using dplyr!

  if (length(db.names)>1) {

    pos.multidbs <- which(apply(P[grep("\\.nresult", names(P))]==1,1,sum, na.rm=T)>1) # single hit for multiple databases
    p <- pos.multidbs[apply(as.matrix(P[pos.multidbs, pos.pos]), 1, function(x) length(unique(x[!is.na(x)]))==1)] # single hit for multiple databases & single genomic position
    # p <- pos.multidbs[do.call("==", P[pos.multidbs, pos.pos])] # single hit for multiple databases & single genomic position
    P$pos.multi[p]  <- apply(as.matrix(P[p, pos.pos]), 1, function(x) unique(x[!is.na(x)])) ## P[p, pos.pos[1]]
    P$name.multi[p] <- apply(as.matrix(P[p, pos.name]), 1, function(x) {x[is.na(x)] <- ""; paste(x, collapse=sep1)}) ## do.call(function(...) paste(..., sep=sep), P[p, pos.name])
    P$exon.multi[p] <- apply(as.matrix(P[p, pos.exon]), 1, function(x) if(length(x[!is.na(x)])>0) {x[is.na(x)] <- ""; paste(x, collapse=sep1)} else NA) ## do.call(function(...) paste(..., sep=sep), P[p, pos.exon])
    P[p,c(pos.pos, pos.name, pos.exon)] <- NA

    f.create.single.string <- function(pos.list, p.pos.list.dup) {
      tmp <- unlist(lapply(1:length(pos.list), function(i) paste(pos.list[[i]][-p.pos.list.dup[[i]]], collapse=";")))
      tmp[tmp==""] <- NA
      return(tmp)
    }

    f = function(X) {
      L <- lapply(list(pos.list=pos.pos), function(i) strsplit(unlist(X[i]), split=";"))
      pos <- c(unlist(L$pos.list)[!is.na(unlist(L$pos.list))])
      if (length(pos)==1 || length(unique(sub("\\.pos.*","",names(pos))))==1) {
        X$pos.multi <- paste(pos,collapse=";")
        tmp <- unlist(X[pos.name])[!is.na(X[pos.name])]
        X$name.multi <- ifelse(length(tmp)>0, paste(tmp, sep=sep1), NA)
        tmp <- unlist(X[pos.exon])[!is.na(X[pos.exon])]
        X$exon.multi <- ifelse(length(tmp)>0, paste(tmp, sep=sep1), NA)
        X[c(pos.pos,pos.name,pos.exon)] <- NA
        return(X)
      }
      if (length(db.names)!=2) return(X)
      if (!any(duplicated(pos))) return(X)
      L <- c(L, lapply(list(name.list=pos.name, exon.list=pos.exon), function(i) strsplit(unlist(X[i]), split=";")))
      pos.dup <- pos[duplicated(pos)]
      p.pos.list.dup <- lapply(L$pos.list, function(x) which(x %in% pos.dup))
      X$pos.multi  <- paste(pos.dup, collapse=";")
      X[c("name.multi","exon.multi")] <-
        lapply(L[-1], function(x) paste(sapply(split(unlist(x), match(pos, pos.dup)), paste, collapse=sep1), collapse=";"))
      X[1,c(pos.pos,pos.name,pos.exon)]  <- unlist(lapply(L, f.create.single.string, p.pos.list.dup))
      return(X)
    }

    # finde gleiche Positionen bei multihits und fasse sie zusammen
    cat.step("merge identical genomic positions")
    q <- setdiff(which(apply(P[grep("\\.nresult", names(P))]>=1,1,sum, na.rm=T)>1),p)

    if (ARG_LIST$NCPU == 1) {
      P[q,] <- P[q,] %>% group_by(.data$probe.id) %>% do(f(.data)) %>% ungroup %>% as.data.frame
    } else {
      probe.id <- NULL; rm("probe.id")
      cluster %>%
        cluster_assign(f = f)  %>%
        cluster_assign(db.names = db.names) %>%
        cluster_assign(f.create.single.string = f.create.single.string) %>%
        cluster_assign(pos.pos = pos.pos)  %>%
        cluster_assign(sep1 = sep1)  %>%
        cluster_assign(pos.name = pos.name)  %>%
        cluster_assign(pos.exon = pos.exon)
      P[q,] <- P[q,] %>% partition(probe.id, cluster=cluster) %>% do(f(.data)) %>% collect() %>% as.data.frame
      cluster %>% cluster_call(ls()) %>% unique %>% unlist %>% cluster_rm(cluster=cluster)
    }

    # fasse ungleiche Positionen zusammen
    cat.step("merge different genomic positions")
    for (x in c("pos","name","exon")) # removed "name"
      P[[paste0(x,".single")]]  <- sub("\\|\\|NA$",sep1,sub("NA\\|\\|",sep1, gsub(paste0(sep1,"NA",sep1), paste0(sep1, sep1), fixed=T,
                                                                                do.call(function(...) paste(...,sep=sep1), P[grep(paste0(x,"$"), names(P))]))))
    pos.nSNP <- grep("\\.nSNP", names(P))
    P$nSNP.min <- apply(P[grep("\\.nSNP.min", names(P))], 1, function(x) ifelse(length(x[!is.na(x)])>0, min(x,na.rm=T), NA))
    P$nSNP.max <- apply(P[grep("\\.nSNP.max", names(P))], 1, function(x) ifelse(length(x[!is.na(x)])>0, max(x,na.rm=T), NA))
    #P[P$pos.single==sep, grep("single$",names(P))] <- NA

    cat.step("remove temporary columns")
    P <- P[-c(sapply(db.names, grep, names(P)))]

  } else { # single RNA database

    probe.table <- P[c(1,2,9,11,5,13,6,15,7,8)]
    names(probe.table) <- sub(paste0("^",db.names,"\\."), "",names(P)[c(1,2,10:15,7,8)])
    P <- probe.table
  }

  probe.table <- probe.tables[[1]][1:2]
  probe.table[grep("probe", names(P), invert=TRUE)] <- NA
  probe.table[match(P$probe.id,probe.table$probe.id),] <- P
  names(probe.table) <- names(P)
  probe.table[probe.table==sep1] <- ""

  write.csv(probe.table, ARG_LIST$pasteOUT(paste0("probes_rank",rank,".csv")), row.names = FALSE)
  #save(probe.table, file=ARG_LIST$pasteTMP(paste0("probes_rank",rank,".Rda")))

  return(probe.table)
}
