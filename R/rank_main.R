rank_main <- function(ARG_LIST, rank, sep1, sep2, sep3, cluster) {


  db.names <- ARG_LIST$DB.TABLE$db[ARG_LIST$DB.TABLE$rank==rank]
  probe.tables <- list()
  for (db.name in db.names) {
    load(ARG_LIST$pasteTMP(db.name, "probe.table.Rda"))
    probe.table <- probe.table[order(probe.table$probe.id),]
    probe.tables[[db.name]] <- probe.table
  }
  probe.table <- rank_create_probe_table(ARG_LIST, rank, probe.tables, sep1, cluster)

  #############################################################################

  PS <- rank_create_probeset_table(ARG_LIST, rank, probe.table, probe.tables, sep1, sep2, sep3, cluster)

  #################################################################

  if (rank>1) load(ARG_LIST$pasteTMP(paste0("probes_rank",rank-1,".Rda"))) else load(ARG_LIST$pasteTMP("probes.Rda"))

  # create fasta file of probes having NM>0
  cat.step("create rank probes Rda file")
  # remove probes that cannot influence probeset assignment
  tmp <- probe.tables %>% lapply("[",3) %>% as.data.frame %>%
    apply(1, function(x) {x=x[!is.na(x)]; ifelse(length(x)==0, F, ifelse(ARG_LIST$STRICT_MAPPING,min(x)==0,T))})

  # add rank id for mapped probes
  probes$rank[probes$id %in% probe.tables[[1]]$probe.id[tmp]] <- rank
  # add rank id for usable probes
  rankx <- rank + 0.1
  probes$rank[sub(".*___","", probes$id) %in% unlist(strsplit(PS$res.probes, ","))] <- rankx

  PS0 <- probes %>% filter(rank==0) %>% group_by(.data$ps.id) %>% summarise(n0=sum(rank==0)) # count probesets unmapped probes
  PSx <- probes %>% filter(rank==rankx) %>% group_by(.data$ps.id) %>% summarise(nx=n()) # count probesets winner probes
  PSa <- probes %>% group_by(.data$ps.id) %>% summarise(n=n())                          # count probesets probes
  PS0 <- merge(merge(PSa,PS0,by="ps.id", all.x=T),PSx,by="ps.id", all.x=T)        # merge count tables
  PS0[is.na(PS0)] <- 0

  probes$rank[probes$ps.id %in% PS0$ps.id[PS0$n<=PS0$nx*2 | PS0$n0<=PS0$nx] & probes$rank == 0] <- -rank
  save(probes,  file = ARG_LIST$pasteTMP(paste0("probes_rank",rank,".Rda")))

  return(probes)
}
