prepare.vcf <- function(ARG_LIST, db.name) {
  if (!ARG_LIST$STRICT_MAPPING && file.exists(ARG_LIST$pasteIN(db.name, "snp.vcf")) &&
      !dir_exists(ARG_LIST$pasteIN(db.name, "vcf"))) {
    catn("\t... prepare variant calling file")
    dir_create(ARG_LIST$pasteIN(db.name, "vcf"))
    prepare_vcf(ARG_LIST$pasteIN(db.name,"snp.vcf"), ARG_LIST$pasteIN(db.name, "vcf"), .Platform$file.sep)
  }
}
