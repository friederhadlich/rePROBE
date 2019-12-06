db_main <- function(ARG_LIST, db.name, probes.unchecked, cluster) {

  if (file.exists(ARG_LIST$pasteTMP(db.name,"probe.table.Rda")))  return(TRUE)

  ARGS <- ARG_LIST$DB.TABLE %>% filter(.data$db == db.name)
  ANNO_EXISTS <- ARGS$anno
  DNA_EXISTS  <- ARGS$dna #file.exists(pasteIN(db.name, "genome.fa"))
  VCF_EXISTS  <- dir.exists(ARG_LIST$pasteIN(db.name, "vcf"))
  RNA_EXISTS  <- ARGS$rna #file.exists(pasteIN(db.name, "sequences.fa"))
  DB_LIST <- list(db.name = db.name,
                  ANNO_EXISTS = ANNO_EXISTS, DNA_EXISTS = DNA_EXISTS,
                  VCF_EXISTS  = VCF_EXISTS, RNA_EXISTS  = RNA_EXISTS)

  # load/create annotation
  if (RNA_EXISTS) {
    rna.anno <- rna.exon.list <- rna.exon.list2 <- NULL
    if (!file.exists(ARG_LIST$pasteTMP(db.name,"annotation.Rda")))
      db_annotate(ARG_LIST, db.name)
    load(ARG_LIST$pasteTMP(db.name,"annotation.Rda"))
    DB_LIST$rna.anno <- rna.anno
    DB_LIST$rna.exon.list <- rna.exon.list
    DB_LIST$rna.exon.list2 <- rna.exon.list2
  }

  ## STEP 1 -----------------------------------------------------------------------------------------
  ## - read rna mapping results from file
  ## - add genomic and exonic position

  save(DB_LIST, probes.unchecked, file="db_main.Rda")
  if (!file.exists(ARG_LIST$pasteTMP(db.name,"mapping.table.step1.Rda")))
    db_read_results(ARG_LIST, DB_LIST, probes.unchecked, cluster)

  ## STEP 2 -----------------------------------------------------------------------------------------
  ## - merge genomically identical mapping results of single probes (= transcript variants)

  if (!file.exists(ARG_LIST$pasteTMP(db.name,"mapping.table.step2.Rda")))
    db_merge_results(ARG_LIST, DB_LIST, cluster)

  ## STEP 3 -----------------------------------------------------------------------------------------
  ## - count SNPs and correct

  if (!file.exists(ARG_LIST$pasteTMP(db.name,"mapping.table.step3.Rda")))
    db_correct_snps(ARG_LIST, DB_LIST, cluster)

  ## STEP 4 -----------------------------------------------------------------------------------------
  ## - remove mismatching mapping results with better strata alternatives (less mismatches)
  ## - remove probes multimapping on 1 rna
  ## - merge genomically different mapping results of single probes

  if (!file.exists(ARG_LIST$pasteTMP(db.name,"mapping.table.step4.Rda")))
    db_filter_results(ARG_LIST, DB_LIST)

  ## -----------------------------------------------------------------------------------------------
  ## - create output probe table with nSNP.min and nSNP.max

  if (!file.exists(ARG_LIST$pasteTMP(db.name,"probe.table.Rda")))
    db_create_probe_table(ARG_LIST, DB_LIST, probes.unchecked, cluster)

}
