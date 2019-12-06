#' Run the revised probeset assignment
#'
#' @param ARG_LIST List with user settings created by prepare_input_arguments() or prepare()
#' @param keep_unqualified_probesets Keep unqualified probesets for better background correction and later analysis
#' @param sep1 Separator for different database results at identical genomic position
#' @param sep2 Separator for different results
#' @param sep3 Separator for different gene.alternative results

#' @import multidplyr
#' @importFrom dplyr arrange as.tbl bind_rows collapse collect count do filter group_by mutate n order_by rename select summarise tbl ungroup %>% one_of
#' @importFrom magrittr %<>%
#' @importFrom rlang .data
#' @export
#'
run <- function(ARG_LIST, keep_unqualified_probesets=F, sep1="||", sep2=";", sep3=",") {

  cat( sep = "",
       "#########################################################\n",
       "Welcome! You are running the revised array assignment \n",
       "#########################################################\n\n")

  if (missing(ARG_LIST)) stop("Input is missing!")
  if (!"DB.TABLE" %in% ls(ARG_LIST)) stop("Input is incomplete!")

  with(ARG_LIST, {

    X <- DB.TABLE[1:5]
    Y <- list_available_dbs(ARG_LIST)

    if (!identical(X, Y[match(X$db, Y$db),])) stop("Input is corrupt!")

    if (length(list.files(folder.tmp))>3 && readline("Use existing files from previous run? [y]/n: ") == "n") {
      unlink(c(folder.tmp, folder.out, pasteIN("probes.fasta")), recursive = T)
      dir_create(folder.tmp)
      dir_create(folder.out)
      prepare_chip(ARG_LIST)
    }


    # https://github.com/tidyverse/dplyr/issues/2396 # info to warning message of Rstudio
    if (NCPU > 1) {
      if (!requireNamespace("multidplyr", quietly = T))
        stop("Install package multidplyr by running 'devtools::install_github(\"hadley/multidplyr\")' and rerun, or use NCPU=1 in ARG_LIST!")
      cluster <- new_cluster(NCPU) #create_cluster(NCPU) # important for multidlyr
    }
    load(pasteTMP("probes.Rda"))

    for (rank in sort(unique(DB.TABLE$rank))) {

      probes.unchecked <- probes$id[probes$rank == 0]
      catn("Rank", rank, "- analyse",length(probes.unchecked),"probes ...")
      if (length(probes.unchecked)==0) {
        cat.step("done!")
        break
      }
      db.names <- DB.TABLE$db[DB.TABLE$rank==rank]
      for (db.name in db.names) {
        cat.step(paste("process database", db.name, "..."))
        dir.create(pasteTMP(db.name), F)
        dir.create(pasteOUT(db.name), F)
        if (!file.exists(pasteTMP(db.name,"probe.table.Rda")))
          db_main(ARG_LIST, db.name, probes.unchecked, cluster)
        cat.step("... done!")
      }
      cat.step("done!")

      catn("Rank", rank, "- collapse database results on probeset level and filter ...")
      if ( !(file.exists(pasteTMP(paste0("probes_rank",rank,".Rda",sep="")))) ) {
#        catn("Run ...")
        probes <- rank_main(ARG_LIST, rank, sep1, sep2, sep3, cluster)
      } else
        load(pasteTMP(paste0("probes_rank",rank,".Rda")))
      cat.step(paste("number of mapping probes   : ", sum(round(probes$rank) == rank)))
      cat.step(paste("number of qualified probes : ", sum(probes$rank-0.1 == rank)))
      cat.step(paste("number of removed probes   : ", sum(probes$rank == -rank)))
      cat.step("done!")

    }

    #####################################################################################
    catn("Summarize ... ")
    save(keep_unqualified_probesets, probes, sep1, sep3, file="HG.Rda")
    create_result(ARG_LIST, keep_unqualified_probesets, probes, sep1, sep3)
    cat.step("done!\n")

    pos <- probes$rank>0 & (floor(probes$rank) != probes$rank)
    catn(paste("Final number of qualified probes    : ", sum(pos)))
    catn(paste("Final number of qualified probesets : ", length(unique(probes$ps.id[pos]))))

  })

  cat( sep = "", "\n",
       "##################################################################\n",
       "Revised array assignment finished! Now use package oligo or affy!\n",
       "##################################################################\n")

}
