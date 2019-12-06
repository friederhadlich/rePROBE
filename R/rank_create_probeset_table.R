f.zip.probe.list <- function(X, db.names, sep1, sep2, sep3, RNA.ANNO, MT, db.str) {
  nrow.X <- nrow(X)
  ps.name <- X$probeset.id[1]
  # init variables
  rna.count <- 0; rna.long <- rna.winner <- chr.pos <- chr.name <- ps.nSNP <- rna.variants.str <- exon.str <- NA
  rna.best.alternative <- c(); X.best <- NULL
  is.db.type <- rep(FALSE, length(db.names)); names(is.db.type) <- db.names
  # remove no hits
  pos <- (is.na(X$pos.multi) & is.na(X$pos.single)) | (is.na(X$name.multi) & is.na(X$name.single))
  nrow.Xnh <- sum(pos)
  X <- X[!pos,] # remove no hits
  # remove multi hits
  pos <- grepl(";",X$name.single) | grepl(";",X$name.multi) | (!is.na(X$name.multi)&!is.na(X$name.single)) |
    grepl(";",X$pos.single) | grepl(";",X$pos.multi)
  X.mh <- X[pos,]; X <- X[!pos,]
  cat(ps.name, "\n")
  if (nrow(X.mh)>0) {
    # check for containing >1 known chromosome (having no dot in name, i.e. 1..28,X,Y,MT)
    chr.list <- lapply(strsplit(paste(X.mh$pos.multi,X.mh$pos.single,sep=sep2),sep2),
                       function(X) {x=sub(":.*","",unlist(strsplit(X, split=sep1, fixed=T))); return(x[x!="" & x!="NA"])})
    pos <- which((sapply(chr.list, function(x) sum(!grepl("\\.", unique(x))))<2)) # multiple known chromosomes for single probe
    if (length(pos)>0) {
      tmp <- sapply(chr.list[pos], function(x) !any(duplicated(x)))
      #pos <- pos[tmp] # single chromosome with no multiple hits
      X <- rbind(X, X.mh[pos[tmp],])
      X.mh <- X.mh[pos[!tmp],]
    }
  }

  # recreate gene-symbol-variant names
  f <- function(S) {
    tmp <- sapply(strsplit(S, ")"), function(L) paste(lapply(strsplit(L, split = "\\("), function(Y)
    {paste0(Y[1], "(", gsub(",", paste0(");",sub("^[|;=]+","",Y[1]),"("), Y[-1]), ")")}), collapse=""))
    tmp[is.na(S)] <- NA
    tmp <- sub("()","",tmp, fixed = T)
    return(tmp)
  }
  X$name.symbol.multi  <- f(X$name.multi) #gsub("\\([A-Za-z 0-9,_]+\\)","",X$name.multi)
  X$name.symbol.single <- f(X$name.single) #gsub("\\([A-Za-z 0-9,_]+\\)","",X$name.single)

  # find best choice single hit (includes valid multihits)
  if (nrow(X)>0) {
    SEP <- gsub("|","\\|",sep1, fixed = T)
    X[c(7:8,12:13)] <- lapply(X[c(7:8,12:13)], function(x) sub(paste0(SEP,"$"), paste0(sep1,sep1), x))
    SYM <- as.data.frame(do.call(rbind, strsplit(c(X$name.symbol.single,X$name.symbol.multi), split=sep1, fixed=T)))
    if (ncol(SYM)==1)  SYM <- do.call(cbind, list(SYM, lapply(2:length(db.names), function(x) SYM[1])))#NA)))
    names(SYM) <- db.names
    SYM$probe.nr <- rep(1:(nrow(SYM)/2), times=2)
 #   cat(ps.name,"\n")
    ff <- function(x, probe.nr) {
      y <- unlist(lapply(split(x,probe.nr), function(z) unique(unlist(strsplit(as.character(z[!is.na(z)]), paste0(sep2,"|",sep3))))));
      sort(table(y[y!=""]),decreasing = TRUE)
    }
    rna.count.list <- lapply(SYM[-ncol(SYM)], ff, probe.nr=SYM$probe.nr)

    rna.count <- sort(unlist(rna.count.list),decreasing = TRUE)
    if (length(rna.count)==0) return(  data.frame(
      probe.set     = ps.name,
      n.probe       = nrow.X,
      n.probe.sh    = nrow.X - nrow(X.mh) - nrow.Xnh,  # single hit
      n.probe.mh    = nrow(X.mh),                      # multi hit
      n.probe.nh    = nrow.Xnh,                        # no hit
      res.db       = NA,
      res.name     = NA,      res.longname = NA,      res.mapping  = 0,
      res.nSNP     = NA,      res.chr.name = NA,      res.chr.pos  = NA,
      res.chr.exon = NA,      res.variant  = NA,      res.alternative = NA,
      res.probes   = NA,      stringsAsFactors = F
    )  )
    rna.count <- rna.count[rna.count==max(rna.count)]
    rna.best  <- sub(db.str, "", names(rna.count))
    rna.winner <- c(grep("^$|^LOC|^ENS[A-Z]{4}[0-9]+",rna.best, value=TRUE, invert=TRUE), rna.best)[1]
    rna.best.alternative <- setdiff(rna.best,rna.winner)

    db.type <- sub(paste0(db.str,".*"),"\\1",names(rna.count))[rna.winner==sub(db.str, "", names(rna.count))]
    pos <- unique(SYM$probe.nr[grep(rna.winner, SYM[[db.type[1]]], fixed = T)])
    if (length(pos)>max(rna.count)) {
      pos <- unique(SYM$probe.nr[rna.winner==SYM[[db.type[1]]]])
      pos <- pos[!is.na(pos)]
    }
    X.best <- X[pos,]
    rna.count <- nrow(X.best) # update rna.count
    if (length(db.names)>1 && (length(rna.best.alternative)>0)) { # || length(db.type)>length(rna.winner))) { # find winner for other dbs
      rna.count.list <- lapply(SYM[SYM$probe.nr %in% pos,-ncol(SYM),drop=F], ff, probe.nr=SYM$probe.nr[SYM$probe.nr %in% pos])
      rna.count.list <- lapply(rna.count.list, function(x) x[x==rna.count])
      rna.count.list <- rna.count.list[sapply(rna.count.list, length)>0]
      rna.count.list <- lapply(rna.count.list, function(x) {
        pos <- grepl("^$|^LOC|^ENS[A-Z]{4}[0-9]+", names(x))
        c(x[!pos], x[pos])
      })
      tmp <- names(sapply(rna.count.list, "[", 1))
      rna.winner <- sub(db.str, "", tmp)
      db.type <- sub(paste0(db.str,".*"),"\\1",tmp)
    }
    is.db.type[db.type] <- TRUE
    rna.best  <- rna.winner
    if (length(rna.best) < length(db.type)) {
       names(rna.best) <- setdiff(db.type, "dna")
    } else  names(rna.best) <- db.type

    if (length(rna.winner)>1 || length(rna.best.alternative)>0) {
      rna.variants.str <- paste(sapply(rna.count.list, function(x) paste(names(x)[-1], collapse=",")), collapse=sep1)
      rna.variants.str[gsub(sep1,"",rna.variants.str,fixed = T)==""] <- NA
      rna.best.alternative <- setdiff(rna.best.alternative, sub(db.str, "", names(unlist(rna.count.list))))
    }

    POS  <- as.data.frame(do.call(rbind, strsplit(c(X.best$pos.single,X.best$pos.multi), split=SEP)), stringsAsFactors = F)
    if ("dna" %in% db.names && !any(is.na(unlist(POS))))  is.db.type["dna"] <- TRUE
    if (length(POS)==0) { chr.list <- lapply(db.type, function(db) NA)
    } else {
      if (ncol(POS)==1) POS <- POS[rep(1, length(db.names))]
      names(POS) <- db.names
      chr.list <- lapply(db.type, function(db) {
        tmp.list <- split(POS[[db]], SYM$probe.nr[1:nrow(POS)])
        if (any(sapply(tmp.list, function(x) all(is.na(x))))) return(NA) # >0 probes have no genomic position
        tmp <- POS[[db]]
        tmp <- tmp[!is.na(tmp)]
        chr.pos <- as.numeric(unlist(strsplit(sub("^.*:","",tmp), "-")))
        if (any(is.na(chr.pos))) return(NA) # >0 probes have no genomic position
        chr.pos <- chr.pos[!is.na(chr.pos)]
        chr.names <- rep(sub(":.*","",tmp), each=2)
        chr.names <- chr.names[which(chr.names!="")]
        tmp <- sapply(split(chr.pos, chr.names), function(x)
          ifelse(chr.pos[2]>chr.pos[1], paste(min(x),max(x), sep="-"), paste(max(x),min(x), sep="-")))
        tmp2 <- paste(tmp, collapse=";")
        names(tmp2) <- paste(names(tmp), collapse=";")
        return(tmp2)
      })
    }
    names(chr.list) <- db.type

    X.best[c('exon.single','exon.multi')] <- lapply(X.best[c('exon.single','exon.multi')], sub, pattern="\\|\\|$", replacement=paste0(sep1,sep1))
    EXON  <- as.data.frame(do.call(rbind, strsplit(c(X.best$exon.single,X.best$exon.multi), split=SEP)), stringsAsFactors = F)
    if (any(grepl("E[0-9]",EXON))) {

      if (ncol(EXON)==1) EXON <- do.call(cbind, list(EXON, lapply(2:length(db.names), function(x) EXON[1])))
      names(EXON) <- db.names
      exon.str <- paste(sapply(db.type, function(db) {
        exons <- EXON[[db]]
        if (any(grepl(",",exons))) {
          tmp <- lapply(strsplit(c(X.best$name.symbol.single, X.best$name.symbol.multi), "\\||;|,"), function(x) which(x==rna.best[db]))
          exons <- unlist(lapply(1:nrow(X.best), function(i,n) {
            y=unlist(strsplit(exons[i],","))[unlist(tmp[c(i,n+i)])]}, n=nrow(X.best)))
          exons[is.na(exons)] <- sub(",.*","",EXON[[db]][is.na(exons)])
          exons <- exons[!is.na(exons)]
        }
        exon.range = range(as.numeric(unlist(strsplit(gsub("E","",exons), "-"))))
        ifelse(nrow(X)==0, NA, ifelse(diff(exon.range)==0, paste0("E",exon.range[1]), paste(paste0("E",exon.range),collapse="-")))
      }), collapse=sep1)

    }

    rna.long.list <- lapply(1:length(db.type), function(i) {
      db <- db.type[i]
      if (! "symbol" %in% names(RNA.ANNO[[db]])) return(NA)
      pos <- which(RNA.ANNO[[db]]$short == rna.best[i])
      paste(unique(RNA.ANNO[[db]]$name[pos]), collapse=";")
    })
    chr.name <- ifelse(is.na(unlist(chr.list))[1], NA, paste(unique(sapply(chr.list, names)),  collapse=sep1))
    chr.pos  <- ifelse(is.na(unlist(chr.list))[1], NA, paste(unique(unlist(chr.list)),  collapse=sep1))
    rna.long <- ifelse(is.na(unlist(rna.long.list))[1], NA, paste(unique(unlist(rna.long.list)), collapse=sep1))

    # update ps.nSNP if some hits did not belong to best choice
    ps.nSNP <- sum(X.best$nSNP.min)
    if (!is.na(ps.nSNP) && ps.nSNP<sum(X.best$nSNP.max)) {
      # sort chromosome names (first knowns)
      tmp <- unique(sapply(chr.list, names));      tmp <- c(tmp[!grepl("\\.",tmp)],tmp)[1]
      pos <- match(as.numeric(sub(".*___","",X.best$probe.id)), MT[[db.type[1]]]$probe.id, nomatch=0)
      ps.nSNP <- sum(MT[[db.type[1]]]$nSNP[pos][MT[[db.type[1]]]$chr.name[pos]==tmp])
    }

  }
  data.frame(
    probe.set     = ps.name,
    n.probe       = nrow.X,
    n.probe.sh    = nrow.X - nrow(X.mh) - nrow.Xnh,  # single hit
    n.probe.mh    = nrow(X.mh),                      # multi hit
    n.probe.nh    = nrow.Xnh,                        # no hit
    res.db       = ifelse(any(is.db.type), paste(names(is.db.type)[is.db.type], collapse=sep1), NA),
    res.name     = paste(rna.winner, collapse=sep1),
    res.longname = rna.long,
    res.mapping  = rna.count[1],
    res.nSNP     = ps.nSNP,
    res.chr.name = chr.name,
    res.chr.pos  = chr.pos,
    res.chr.exon = exon.str,
    res.variant  = rna.variants.str,
    res.alternative = ifelse(length(rna.best.alternative)==0, NA, paste(rna.best.alternative, collapse=",")),
    res.probes   = ifelse(nrow(X)==0, NA, paste(sort(as.numeric(sub(".*___","",X.best$probe.id))),collapse=",")),
    stringsAsFactors = F
  )
}

rank_create_probeset_table <- function(ARG_LIST, rank, probe.table, probe.tables, sep1, sep2, sep3, cluster) {

  ## create final probeset.table
  cat.step("create probeset table")

  db.names <- ARG_LIST$DB.TABLE$db[ARG_LIST$DB.TABLE$rank==rank]

  P <- probe.table
  P$name.multi[P$name.multi==""] <- NA; P$name.single[P$name.single==""] <- NA

  RNA.ANNO <- list()
  RNA.ANNO.EXON <- list()
  rna.exon.list <- NULL
  for (db.name in db.names) {
    if (!file.exists(ARG_LIST$pasteTMP(db.name,"annotation.Rda")))  next
    load(ARG_LIST$pasteTMP(db.name,"annotation.Rda"))
    if ("symbol" %in% names(rna.anno)) {
      rna.anno$symbol <- sub("\\(.*","",rna.anno$symbol)
      RNA.ANNO.EXON[[db.name]] <- rna.exon.list
    }
    RNA.ANNO[[db.name]] <- rna.anno
  }

  # important for correction of gene.nSNP in f.zip.probe.list
  MT <- list()
  for (db.name in db.names) {
    mapping.table.probe <- NULL
    load(ARG_LIST$pasteTMP(db.name, "mapping.table.step4.Rda"))
    MT[[db.name]] <- mapping.table.probe[,c("probe","chr.name","nSNP")]
    MT[[db.name]]$probe.id <- as.numeric(sub(".*___","",MT[[db.name]]$probe))
  }


  ## find dominating gene for each probeset
  P$probeset.id <- sub("\\..*","", sub("___.*","",P$probe.id))
  db.str <- paste0("(",paste(db.names, collapse="|"),")\\.")

  P$name.single <- as.character(P$name.single)
  # options(warn=2)
  #ids=15385015; XX=P[P$probeset.id==as.character(ids),]; X=XX; f.zip.probe.list(XX, db.names, sep1, sep2, sep3, RNA.ANNO, MT, db.str)
  if (ARG_LIST$NCPU == 1) {
    PS <- P %>% group_by(.data$probeset.id) %>%
      do(f.zip.probe.list(.data, db.names, sep1, sep2, sep3, RNA.ANNO, MT, db.str)) %>% collect()
  } else {
    #  cluster <- create_cluster(10)
    probeset.id <- NULL; rm("probeset.id")
    cluster %>%
      cluster_assign(f.zip.probe.list = f.zip.probe.list)  %>%
      cluster_assign(db.names = db.names) %>%
      cluster_assign(db.str = db.str)  %>%
      cluster_assign(RNA.ANNO = RNA.ANNO)  %>%
      cluster_assign(RNA.ANNO.EXON = RNA.ANNO.EXON) %>%
      cluster_assign(MT = MT) %>%
      cluster_assign(sep1 = sep1) %>%
      cluster_assign(sep2 = sep2) %>%
      cluster_assign(sep3 = sep3)
    PS <- P %>% partition(probeset.id, cluster=cluster) %>%
      do(f.zip.probe.list(.data, db.names, sep1, sep2, sep3, RNA.ANNO, MT, db.str)) %>% collect()
    cluster %>% cluster_call(ls()) %>% unique %>% unlist %>% cluster_rm(cluster=cluster)
  }
  PS %<>% ungroup %>% select(-1) %>% arrange(.data$probe.set)

  if ("dna" %in% db.names) { # check "dna"-only probe hits separately ...

    pos <- grepl("dna", PS$res.db) & is.na(PS$res.chr.name);
    PS$res.db[pos] <- sub("||dna","",sub("dna||","",PS$res.db[pos], fixed = T), fixed=T)

    f <- function(X) {
      #    cat(X$probeset.id[1],"\n")
      POS <- lapply(split(X$result.pos, X$probe.id), function(x) unlist(strsplit(x, split = ";")))
      POS <- POS[!is.na(X$nresult)]
      chr.names.list <- lapply(lapply(POS, sub, pattern=":.*", replacement=""), unique) #as.numeric(unlist(strsplit(sub("^.*:","",tmp), "-")))
      tmp <- table(unlist(chr.names.list))
      chr.names <- names(tmp)[tmp==tmp[1]]
      tmp <- as.data.frame(do.call(rbind, strsplit(unlist(POS), ":")), stringsAsFactors = F)
      tmp$V3 <- rep(names(POS), sapply(POS, length))
      tmp <- tmp[tmp[[1]] %in% chr.names,]
      #tmp <- split(as.character(tmp[[2]]), tmp[1])
      tmp <- do.call(rbind, lapply(split(tmp, tmp$V1), function(tmpp) {
        pos <- match(X$probe.id, tmpp$V3)
        y <- as.numeric(unlist(strsplit(sub("^.*:","",tmpp$V2), "-")))
        data.frame(
          probe.set     = X$probeset.id[1],
          n.probe       = nrow(X),
          n.probe.sh    = sum(X$nresult==1),  # single hit
          n.probe.mh    = sum(X$nresult>1),   # multi hit
          n.probe.nh    = sum(X$nresult==0),  # no hit
          res.db       = "dna",
          res.name     = NA,        res.longname = NA,
          res.mapping  = sum(X$nresult==1), #nrow(tmpp),
          res.nSNP     = sum(MT$dna$nSNP[MT$dna$probe %in% X$probe.id]),#paste(unique(sum(X$nSNP.min[pos]), sum(X$nSNP.max[pos])), collapse="-"),
          res.chr.name = tmpp$V1[1],
          res.chr.pos  = ifelse(y[1]<y[2], paste0(min(y),"-",max(y)), paste0(max(y),"-",min(y))),
          res.chr.exon = NA,        res.variant = NA,        res.alternative = NA,
          res.probes   = paste(sort(as.numeric(sub(".*___", "", names(table(tmpp$V3))[table(tmpp$V3)==1]))), collapse=","),
          stringsAsFactors = F)
      }))
      do.call(rbind, lapply(split(tmp, tmp$res.probes), function(X) {
        Y=X[1,]; Y$res.chr.name=paste(X$res.chr.name, collapse=";"); Y$res.chr.pos=paste(X$res.chr.pos, collapse=";"); Y
      }))
    }

    P.dna <- probe.tables$dna %>% mutate(probeset.id=sub("\\..*","", sub("___.*","", .data$probe.id))) %>% as.tbl %>% filter(!is.na(.data$nresult))
    if (ARG_LIST$NCPU == 1) {
      PS.dna <- P.dna %>% group_by(.data$probeset.id) %>% do(f(.data)) %>% collect()
    } else {
      cluster %>% cluster_assign(f = f) %>% cluster_assign(MT = MT)
      PS.dna <- P.dna %>% partition(probeset.id, cluster=cluster) %>% do(f(.data)) %>% collect()
      cluster %>% cluster_call(ls()) %>% unique %>% unlist %>% cluster_rm(cluster=cluster)
    }
    PS.dna <- PS.dna %>% ungroup %>% select(-.data$probeset.id)
    # pos <- which(PS.dna$n.probe.sh==0)
    # if (length(pos)>0) {
    #   PS.dna$res.db[pos] <- NA
    #   PS.dna$res.mapping[pos] <- 0
    # }

    #  X <- rbind(PS.nok, PS.dna[-1]) %>% group_by(probe.set) %>% filter(probe.set=="15180001")
    f <- function(X) {
      if (nrow(X)==1 || all(is.na(X$res.db))) return(X[1,,drop=F])
      n.probe <- X$n.probe[1] # update n.probe if needed
      X <- X[!is.na(X$res.db),]
      if (nrow(X)==1) return(X[1,,drop=F])
      pos.winner <- which(X$res.mapping == max(X$res.mapping))[1] # keeps rna result being winner!
      return(X[pos.winner,,drop=F])
    }
    PS.ok <- PS %>% filter(grepl("dna", .data$res.db) | .data$res.mapping==.data$n.probe)
    PS.nok <- PS %>% filter(! .data$probe.set %in% PS.ok$probe.set)
    PS.all <- PS.dna %>% filter(.data$probe.set %in% PS.nok$probe.set) %>% bind_rows(PS.nok) %>%
      group_by(.data$probe.set) %>% do(f(.data)) %>% collapse
    PS <- PS.all %>% bind_rows(PS.ok) %>% arrange(.data$probe.set)

  }

  save(PS, file=ARG_LIST$pasteTMP(paste0("probeset_rank",rank,".Rda")))
  write.csv(PS,file = ARG_LIST$pasteOUT(paste0("probeset_rank",rank,".csv")), row.names = FALSE)

  return(PS)
}
