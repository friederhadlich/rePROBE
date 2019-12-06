#' Prepare chip files and create file "probes.fasta"
#' @param ARG_LIST list with user settings created by prepare_input_arguments.R
#'
#' @export
prepare_chip <- function(ARG_LIST) {

  if (!dir_exists(ARG_LIST$folder.in) || !any("CHIP_PATH" == ls(ARG_LIST)))  stop("Wrong input for prepare_chip!")

  with(ARG_LIST, {

    if ( !(all(file.exists(pasteIN("chip",c("chip.pgf", "chip.mps")))) || all(file.exists(pasteIN("chip",c("chip.fa", "chip.cdf"))))) ||
         (is_linux() && is_link(pasteIN("chip","chip.pgf")) &&
         (sub(".[0-9a-zA-Z_-]+.pgf","",link_path(pasteIN("chip","chip.pgf"))) != CHIP_PATH)) )
      unlink(pasteIN("chip"))

    is_pgf <- length(list.files(CHIP_PATH, pattern = "*.pgf", full.names = T))==1
    if (!dir.exists(pasteIN("chip"))) {
      dir_create(pasteIN("chip"))
      if (is_pgf) {
        link_create(list.files(CHIP_PATH, pattern = "*.pgf", full.names = T), pasteIN("chip","chip.pgf"), is_linux())
        link_create(list.files(CHIP_PATH, pattern = "*.mps", full.names = T), pasteIN("chip","chip.mps"), is_linux())
      } else {
        link_create(list.files(CHIP_PATH, pattern = "*.cdf", full.names = T), pasteIN("chip","chip.cdf"), is_linux())
        link_create(list.files(CHIP_PATH, pattern = "*fasta$", full.names = T), pasteIN("chip","chip.fa"), is_linux())
      }
    }

    if (!file.exists(pasteIN("probes.fasta")) ||
        as.POSIXct(file_info(pasteIN("probes.fasta"))$change_time) <
        as.POSIXct(file_info(list.files(path = pasteIN("chip"), pattern = "chip", full.names = T)[1])$change_time)) {

      catn("Prepare chip files ...")

      if (is_pgf) {

        cat.step("read pgf + mps")

        lines <- readLines(pasteIN("chip", "chip.pgf"))
        lines.header <- grep(".*header.*=", lines, value=TRUE)
        CHIP_TYPE <- ifelse(grepl("chip_type=",lines[1]), sub("^.*chip_type=","",lines[1]), NA)
        if (is.na(CHIP_TYPE) && any(grep("chip_type=", lines[1:5])))
          CHIP_TYPE <- sub(".*chip_type=", "", grep("chip_type=", lines[1:5], value=T))
        lines <- grep("^#", lines, value=TRUE, invert=TRUE)
        lines.probeset.table <- grep("^\t", lines, value=TRUE, invert=TRUE)
        lines.probe.table <- sub("^\t\t", "", grep("^\t\t.*[mp]m:", value=TRUE, lines))

        mps.table <- read.delim(pasteIN("chip", "chip.mps"), comment.char = "#", stringsAsFactors = F)

        f.read.delim <- function(lines) as.data.frame(do.call(rbind, strsplit(lines,"\t")))

        # create probeset table
        cat.step("process pgf + mps")
        probeset.table <- f.read.delim(lines.probeset.table)
        if (ncol(probeset.table)==2) probeset.table$V3 <- probeset.table$V1
        names(probeset.table) <- strsplit(sub(".*=","",lines.header[1]), "\t")[[1]]
        tmp <- strsplit(mps.table$probeset_list, " ")
        tmp2 <- data.frame(probeset_id=rep(mps.table$probeset_id, sapply(tmp, length)), probesubset_id=unlist(tmp), stringsAsFactors = F)
        probeset.table$probeset_name2 <- tmp2$probeset_id[match(probeset.table$probeset_id, tmp2$probesubset_id)]
        names(probeset.table) <- c("pss.id","pss.type","pss.name", "ps.id")
        pos <- is.na(probeset.table$ps.id)
        probeset.table$ps.id[pos] <- probeset.table$pss.id[pos]
        probeset.table.filtered <- probeset.table[!grepl("^control->|^normgene->|^[rR]eporter|^rRNA$|^Introns$",
                                                         probeset.table$pss.type),]
        write.csv(probeset.table, pasteTMP("probesubsetting.csv"), row.names = F)

        # create probe table
        probe.table <- f.read.delim(lines.probe.table)
        names(probe.table) <- strsplit(sub(".*=\t\t","",lines.header[3]), "\t")[[1]]
        probe.table$probesubset_id <- as.integer(cut(grep("^\t\t.*[mp]m:", lines), c(grep("^\t", invert=TRUE, lines), length(lines)+1)))
        probe.table$probesubset_name <- probeset.table$pss.name[probe.table$probesubset_id]
        probe.table$probeset_id <- probeset.table$ps.id[probe.table$probesubset_id]
        probe.table$probesubset_id <- probe.table$probeset_id # comment this line if you want to process PSS instead of PS
        probe.table$name <- paste0(probe.table$probesubset_id, "___", probe.table$probe_id)
        pos <- grep("^mm",probe.table$type)
        probe.table$mm <- NA
        if (length(pos)>0)  {
          pos <- pos[substr(probe.table$probe_sequence[pos-1], 1, 10) == substr(probe.table$probe_sequence[pos], 1, 10)]
          probe.table$mm[pos-1] <- probe.table$probe_id[pos]
          probe.table <- probe.table[-pos,]
        }
        probe.table.filtered <- probe.table[probe.table$probesubset_id %in% probeset.table.filtered$ps.id,]

        # duplicate probes
        pos <- match(probe.table.filtered$probe_sequence, probe.table.filtered$probe_sequence[duplicated(probe.table.filtered$probe_sequence)])
        pos.dup <- pos[!is.na(pos)]
        probe.table.filtered.dup <- probe.table.filtered[!is.na(pos),]
        tmp <- sapply(split(probe.table.filtered.dup[,c(6,8)], pos.dup), function(x) paste(c(as.character(x[1,1]), as.character(x[,2])), collapse=","))
        writeLines(tmp, pasteOUT("duplicate_probes.csv"))

        # create probe fasta
        cat.step("create fasta")
        lines <- rep(paste0(">",probe.table.filtered$name), each=2)
        lines[seq(2, length(lines), by=2)] <- chartr("ACGT","TGCA",as.character(probe.table.filtered$probe_sequence))
        writeLines(lines, pasteIN("probes.fasta"))

      } else {

        cat.step("read fasta + cdf")

        lines <- readLines(pasteIN("chip", "chip.cdf"))
        CHIP_TYPE <- sub("Name=", "", grep("Name=",lines[1:10], value=T))
        lines.control <- lines[1:(grep("\\[Unit", lines)[1]-1)]
        lines.unit <- lines[grep("\\[Unit", lines)[1]:length(lines)]
        probe.table <- read.delim(textConnection(grep("^Cell[0-9]",lines.unit, value=T)), header=F, stringsAsFactors = F)
        probe.table$str <- paste(probe.table$V5, sub(".*=","",probe.table$V1), probe.table$V2, sep=":")

        # create probe fasta
        cat.step("create fasta")
        lines.fasta <-readLines(pasteIN("chip","chip.fa"))
        pos <- grep("^>", lines.fasta)
        lines.fasta[pos] <- sub("; .*", "", sub(paste0(".*",CHIP_TYPE,":"),"",lines.fasta[pos]))
        pos2 <- match(lines.fasta[pos], probe.table$str)
        pos3 <- match(probe.table$str, lines.fasta[pos])
        lines.fasta[pos] <- paste0(">", probe.table$V5[pos2], "___", probe.table$V12[pos2])
        writeLines(lines.fasta, pasteIN("probes.fasta"))

        # create probe.table
        names(probe.table)[c(12,7)] <- c("probe_id", "interrogation_position")
        probe.table <- cbind(probe.table, type="pm:st", gc_count=NA, probe_length=25, probe_sequence=lines.fasta[pos+1][pos3],
                             probesubset_id=probe.table$V5, probesubset_name=probe.table$V5, probeset_id=probe.table$V5,
                             name=sub("^>","",lines.fasta[pos][pos3]), stringsAsFactors=F)
        probe.table$type[is.na(probe.table$probe_sequence)] <- "mm:st"
        probe.table <- probe.table[order(probe.table$V5, probe.table$V6, probe.table$type),]
        pos <- which(probe.table$type=="mm:st")
        tmp <- probe.table[pos,]
        probe.table <- probe.table[pos+1,]
        probe.table$mm <- NA
        probe.table$mm[match(paste(tmp$probesubset_id, tmp$V11), paste(probe.table$probesubset_id, probe.table$V11))] <- tmp$probe_id
        ## change of 13th base ... not needed!
        # names(probe.table)[c(12,7)] <- c("probe_id", "interrogation_position")
        # probe.table <- cbind(probe.table, type="mm:st", gc_count=NA, probe_length=25, probe_sequence=lines.fasta[pos+1][pos3],
        #                      probesubset_id=probe.table$V5, probesubset_name=probe.table$V5, probeset_id=probe.table$V5,
        #                      name=sub("^>","",lines.fasta[pos][pos3]), stringsAsFactors=F)
        # probe.table$type[is.na(probe.table$probe_sequence)] <- "pm:st"
        # probe.table <- probe.table[order(probe.table$V5, probe.table$V6, probe.table$type),]
        # pos <- which(probe.table$type=="mm:st")
        # stopifnot(length(pos) != nrow(probe.table)/2)
        # probe.table$probe_sequence[pos+1] <-
        #   paste0(substr(probe.table$probe_sequence[pos],1,12),probe.table$V10[pos+1], substr(probe.table$probe_sequence[pos],14,25))
        # tmp <- probe.table[pos,]
        # probe.table <- probe.table[pos+1,]
        # probe.table$mm <- NA
        # probe.table$mm[match(paste(tmp$probesubset_id, tmp$V11), paste(probe.table$probesubset_id, probe.table$V11))] <- tmp$probe_id
        # lines.fasta <- rep(paste0(">", probe.table$V5, "___", probe.table$probe_id), each=2)
        # lines.fasta[seq(2, length(lines.fasta), 2)] <- probe.table$probe_sequence
        # writeLines(lines.fasta, pasteIN("probes.fasta"))


        probe.table <- probe.table[c("probe_id", "type", "gc_count", "probe_length", "interrogation_position", "probe_sequence",
                                     "probesubset_id", "probesubset_name", "probeset_id", "name", "mm")]
        probe.table.filtered <- probe.table

        probe.table.control <- read.delim(textConnection(grep("^Cell[0-9]", lines.control, value=T)), header=F, stringsAsFactors = F)

        control.ps.names <- gsub("(\\[|\\])","",grep("\\[", lines.control[-1:-5], value=T))
        control.ps.n <- diff(c(grep("^Cell1=", probe.table.control$V1), nrow(probe.table.control)+1))
        tmp <- rep(control.ps.names, control.ps.n)
        probe.table.control <- data.frame(probe_id=probe.table.control$V6, type="control", gc_count=NA, probe_length=25,
                                          interrogation_position=NA, probe_sequence=NA, probesubset_id=tmp,
                                          probesubset_name=tmp, probeset_id=tmp,
                                          name=paste(tmp, probe.table.control$V6, sep="___"), mm=NA, stringsAsFactors=F)
        probe.table <- rbind(probe.table.control, probe.table)

        # create probeset.table
        probeset.table <- data.frame(pss.id = unique(probe.table$probeset_id), stringsAsFactors = F)
        probeset.table$pss.type <- "main"
        probeset.table$pss.type[probeset.table$pss.id %in% control.ps.names] <- "control"
        probeset.table$pss.name <- probeset.table$pss.id
        probeset.table$ps.id <- probeset.table$pss.id

      }

      # save probes
      probes <- data.frame(id=paste0(probe.table$probesubset_id,"___",probe.table$probe_id),
                           ps.id=probe.table$probesubset_id,
                           seq=chartr("ACGT","TGCA",as.character(probe.table$probe_sequence)),
                           rank=0, stringsAsFactors = FALSE)
      probes$rank[! probes$id %in% probe.table.filtered$name] <- -9.9
      probes <- probes[order(probes$id),]
      save(probe.table, probeset.table, CHIP_TYPE, file=pasteTMP("init.Rda"))
      save(probes, CHIP_TYPE, file = pasteTMP("probes.Rda"))
    }

    return(TRUE)
  })

}
