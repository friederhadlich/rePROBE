#' Show facts about the finished revision
#' @param ARG_LIST list with user settings created by prepare_input_arguments() or prepare_data()
#' @export
show_report <- function(ARG_LIST) {

  if (missing(ARG_LIST))  stop("Input is missing!")
  if (!all(file.exists(ARG_LIST$pasteTMP(c("init.Rda", paste0("probes_rank", max(ARG_LIST$DB.TABLE$rank),".Rda"), "probeset.Rda")))))
    stop("Required output files from function \"run\" are missing.")

  probes <- probeset.table <- probe.table <- CHIP_TYPE <- PS.final <- NULL
  load(ARG_LIST$pasteTMP("init.Rda"))
  ######## TODO: load correct probes.table

  #probe.table <- read.csv(ARG_LIST$pasteOUT("probes.csv"), stringsAsFactors = F)
  load(ARG_LIST$pasteTMP(paste0("probes_rank", max(ARG_LIST$DB.TABLE$rank),".Rda")))
  probes_no_ctrl <- probes %>% filter(.data$rank!=-9.9)
  probesets_no_ctrl <- probeset.table %>% filter(.data$pss.type == "main")

  catn("")
  catn(paste("Processed chip:", CHIP_TYPE))
  catn("")

  catn("Original array ...")
  cat.step(paste("investigated probes     :", nrow(probes_no_ctrl)))
  cat.step(paste("investigated probe sets :", length(unique(probesets_no_ctrl$ps.id))))
  cat.step(paste("control probes          :", sum(probes$rank==-9.9)))
  cat.step(paste("main probe duplicates   :", sum(duplicated(probes$seq[probes$rank!=-9.9]))))
  catn("")

  load(ARG_LIST$pasteTMP("probeset.Rda"))
  pos <- probes$rank>0 & (floor(probes$rank) != probes$rank)
  # PS <- PS.final %>% filter(!.data$probe.set %in% probes$ps.id[pos]) #.data$n.probe > .data$n.probe.nh)
  # catn("Disqualified for revised array ...")
  # cat.step(paste("probes due to multimapping     : ", sum(PS.final$n.probe.mh)))
  # cat.step(paste("probe sets due to multimapping : ", sum(PS$n.probe.mh>0)))
  # cat.step(paste("probes due to no mapping       : ", sum(probes$rank == 0)))
  # cat.step(paste("probe sets due to no mapping   : ", sum(PS$n.probe == PS$n.probe.nh)))
  # cat.step(paste("probes due to found set winner : ", sum(probes$rank<0 & probes$rank==floor(probes$rank))))
  # cat.step(paste("probe sets 0% ok               : ", sum(is.na(PS.final$res.mapping))))
  # catn("")

  catn("Qualified for revised array ... ")
  PS <- PS.final %>% filter(.data$probe.set %in% probes$ps.id[pos]) #.data$n.probe > .data$n.probe.nh)
  cat.step(paste("probes                : ", sum(pos)))
  cat.step(paste("probe sets            : ", nrow(PS)))
  cat.step(paste("probe sets >= 50% ok  : ", sum(PS.final$res.mapping*2 >= PS.final$n.probe, na.rm=T)))
  cat.step(paste("probe sets   100% ok  : ", sum(PS.final$res.mapping == PS.final$n.probe, na.rm=T)))
  cat.step(paste("unique genes          : ", length(unique(PS$res.symbol))))

  # rank_counts <- table(probes$rank)
  # cat.step(paste("Never map : ", rank_counts["0"]))
  # cat.step(paste("Rank 1 out no map : ", rank_counts["-1"]))
  # cat.step(paste("Rank 2 out no map : ", rank_counts["-2"]))
  # cat.step(paste("Rank 1 out map    : ", rank_counts["1"]))
  # cat.step(paste("Rank 2 out map    : ", rank_counts["2"]))
  # cat.step(paste("Rank 1 in+out     : ", sum(rank_counts[c("1.1","1","-1")])))
  # cat.step(paste("Rank 1 in         : ", rank_counts["1.1"]))
  # cat.step(paste("Rank 2 in    : ", rank_counts["2.1"]))
  catn("")

  catn("Number of qualified probes and probe sets ...")
  print(PS.final %>% group_by(.data$res.db) %>% summarise(n.probeset=n(), n.probe=sum(.data$res.mapping)) %>% as.data.frame)
  cat("\n")

  # probe.tables <- lapply(ARG_LIST$DB.TABLE$db, function(x) {
  #   load(ARG_LIST$pasteTMP(x, "probe.table.Rda"));
  #   X <- data.frame(id = probe.table$probe.id, finite = is.finite(probe.table$nresult), stringsAsFactors = F)
  #   names(X)[2] <- x; return(X)})
  # names(probe.tables) <- ARG_LIST$DB.TABLE$db
  # for (rnk in sort(unique(ARG_LIST$DB.TABLE$rank))) {
  #   catn(paste0("Number of mapped probes rank ",rnk , ":"))
  #   pos <- which(ARG_LIST$DB.TABLE$rank == rnk)
  #   switch(as.character(length(pos)),
  #          '1' = print(table(probe.tables[[pos]][-1])),
  #          '2' = print(table(merge(probe.tables[[pos[1]]], probe.tables[[pos[2]]], by="id")[-1])),
  #          print("Too many databases used for this rank!"))
  # }

}
