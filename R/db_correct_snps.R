db_correct_snps <- function(ARG_LIST, DB_LIST, cluster) {

  cat.step("... correct mismatching mapping results using SNPs")

  with(ARG_LIST, with(DB_LIST, {

    load(pasteTMP(db.name,"mapping.table.step2.Rda"))
    load(pasteTMP(db.name,"mapping.table.raw.Rda"))
    mapping.table.probe$seq <- mapping.table.raw$seq[match(mapping.table.probe$probe, mapping.table.raw$qname)]

    if ((!STRICT_MAPPING && VCF_EXISTS)) {  # SNP check ...

      ff <- function(XX, vcf) {
        vcf.ff <- vcf %>% filter(pos>=min(XX$chr.pos.start) & pos<=max(XX$chr.pos.stop))
        if (nrow(vcf.ff)==0) return(XX)
        # create exonic POS
        gff.id <- sub(".*,","",XX$rna.id[1])
        if (exists("rna.anno")) {
          R <- rna.anno[rna.anno$gff_id==gff.id,,drop=TRUE]
          POS <- sort(unlist(apply(rna.exon.list[[gff.id]], 1, function(x) x[1]:x[2])), decreasing = R$strand=="-")
        } else { # dna (or missing annotation)
          POS <- sort(unique(c(unlist(apply(XX, 1, function(x) seq(x[10],x[11]))))))
          if (XX$strand[1]=="-") POS <- rev(POS)
        }

        as.data.frame(do.call(rbind, lapply(1:nrow(XX), function(i) {

          # count SNPs
          MM <- XX[i,]
          chr.start <- min(MM$chr.pos.start, MM$chr.pos.stop)
          chr.stop  <- max(MM$chr.pos.start, MM$chr.pos.stop)
          A <- vcf.ff %>% filter(pos>=chr.start & pos<=chr.stop)
          if (nrow(A)==0)    return(MM)

          if (chr.stop-chr.start > MM$len) # multiple exons (estimate SNP number)
            A <- A[(chr.start+MM$len-2)>=A$pos | (chr.stop-MM$len-2)<=A$pos,]
          if (nrow(A)==0)  return(MM)

          if (gff.id %in% unlist(strsplit(MM$rna.id, ","))) { # gleiches transkript ...
            POS <- POS[POS>=chr.start & POS<=chr.stop]
          }
          else { # unterschiedliche transkripte ...
            if (chr.stop-chr.start >= MM$len) { # multiple exon - calculate SNP number
              R <- rna.anno[rna.anno$gff_id==sub(".*,","",MM$rna.id),,drop=TRUE]
              if (R$strand == '-') {
                exon <- which(rna.exon.list[[R$gff_id]]$start<=chr.stop)[1]
                POS <- c(chr.stop : rna.exon.list[[R$gff_id]]$start[exon], rna.exon.list[[R$gff_id]]$end[exon+1]:chr.start)
              } else {
                exon <- which(rna.exon.list[[R$gff_id]]$end>=chr.start)[1]
                POS <- c(chr.start : rna.exon.list[[R$gff_id]]$end[exon], rna.exon.list[[R$gff_id]]$start[exon+1]:chr.stop)
              }
            }
          }
          A <- A[!is.na(match(A$pos,POS)),]
          if (nrow(A)==0)    return(MM)
          MM$nSNP <- length(unique(A$pos))

          # correct mismatch
          if(MM$nmismatch==0) return(MM)
          x <- sub(".*:","",MM$rna.mismatch)
          X <- data.frame(pos=cumsum(as.numeric(unlist(strsplit(x, split ="[A-Z]")))[1:MM$nmismatch]+1)) #,
          X$found <- unlist(strsplit(MM$seq, ""))[X$pos]
          strand <- c("+","-")[(chr.start>chr.stop)+1]
          if (strand == '-') {
            X$found <- chartr(old = "TGCA", new = "ACGT", X$found)
          }
          if (chr.stop-chr.start < MM$len) { # mapping at chromosome border
            if (strand == '-') {  X$chr.pos <- chr.stop - X$pos + 1 }
            else  X$chr.pos <- chr.start + X$pos - 1
          } else { X$chr.pos <- POS[X$pos] }

          B <- merge(X, A, by.x="chr.pos", by.y="pos")
          if (nrow(B)==0)  return(MM) # no mismatch position in SNP catalog
          if ((nfound=sum(apply(B[,-(1:2)],1,function(x) length(unique(x))==2)))>0) {
            MM$nmismatch <- MM$nmismatch - nfound
          }
          return(MM)
        })))
      }

      f <- function(X, db.name) {
        X$nSNP <- 0
        chr <- X$chr.name[1]
        vcf.file <- pasteIN(db.name, "vcf", paste0(chr,".vcf.cut"))
        if (!file.exists(vcf.file))    return(X)

        cat.step(paste0("... chromosome =", chr, " :"))
        vcf <- read.delim(vcf.file, comment.char = "#", header=FALSE, col.names = c("pos","ref","alt"), stringsAsFactors = FALSE)
        M <- X %>% group_by(probeset, strand)

        return(M %>% do(ff(.,vcf)) %>% ungroup)
      }

      if (NCPU == 1) { # serial
        mapping.table.probe %<>% group_by(chr.name) %>% do(f(., db.name))
      } else {    # parallel
        cluster %>%
          cluster_assign_value('pasteIN', pasteIN)  %>%
          cluster_assign_value('pastef', pastef) %>%
          cluster_assign_value('folder.in', folder.in)  %>%
          cluster_assign_value('f', f)  %>%
          cluster_assign_value('ff', ff)  %>%
          cluster_assign_value('db.name', db.name) %>%
          cluster_library('dplyr')
        if (ANNO_EXISTS) cluster %>%
          cluster_assign_value('rna.anno', rna.anno) %>%
          cluster_assign_value('rna.exon.list', rna.exon.list)
        mapping.table.probe %<>% partition(chr.name, cluster=cluster) %>% do(f(.data, db.name)) %>% collect()
        cluster %>% cluster_ls() %>% unique %>% unlist %>% cluster_rm(cluster=cluster)
      }

      cat("Found SNPs -> counts of mapping results:"); print(table(mapping.table.probe$nSNP))

    } else   # no SNP check ...
      mapping.table.probe$nSNP <- NA

    save(mapping.table.probe, file=pasteTMP(db.name,"mapping.table.step3.Rda"))

  }))
}
