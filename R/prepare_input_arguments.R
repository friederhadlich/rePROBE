#' Prepare all input arguments for the revised probeset assignment
#'
#' @param WORK_PATH Working directory
#' @param PROJ_NAME Project name
#' @param CHIP_PATH Chip path with one pgf and one mps file
#' @param SPECIES Species name in latin language, i.e. "Sus_scrofa"
#' @param NCPU Number of threads to use
#' @param STRICT_MAPPING Avoid SNP correction of probes
#' @param ASK_QUESTIONS Answer questions automatically with their defaults
#' @useDynLib rePROBE
#' @export
#' @import data.table fs
#' @importFrom utils str
prepare_input_arguments <- function(
  WORK_PATH,
  PROJ_NAME,
  CHIP_PATH,
  SPECIES,
  NCPU,
  STRICT_MAPPING,
  ASK_QUESTIONS) {

  if (missing(WORK_PATH) || !dir_exists(WORK_PATH)) {
    WORK_PATH <- choose_directory("Select an existing working directory!")
    if (!dir.exists(WORK_PATH)) stop("Working directory does not exist!")
  }

  if (missing(PROJ_NAME) || nchar(PROJ_NAME)==0) {
    PROJ_NAME <- readline(prompt = "Define a project name: ")
    if (nchar(PROJ_NAME)==0) stop("Project name is empty!")
  }

  if (file.exists(pastef(WORK_PATH, PROJ_NAME, "in", "ARG_LIST.Rda"))) {
    TMP <- pastef(WORK_PATH, PROJ_NAME, "in", "ARG_LIST.Rda")
    catn("Load missing arguments from existing configuration file!")
    ARG_LIST <- list()
    load(TMP)
    str(ARG_LIST[!grepl("folder|paste|DB.TABLE",names(ARG_LIST))], give.head = F, no.list=T, comp.str = " ... ")
    if (readline("Configuration ok? [y]/n: ") != "n")   return(ARG_LIST)
    rm(TMP, ARG_LIST)
  }
  if (missing(CHIP_PATH) || length(list.files(CHIP_PATH, pattern = "*.pgf"))!=1 ||
         length(list.files(CHIP_PATH, pattern = "*.mps"))!=1) {
    catn("Define a microarray definition directory with pgf and mps file, or a directory with probes fasta file and cdf file!")
    CHIP_PATH <- choose_directory("Select the array definition directory!")
    if (!(length(list.files(CHIP_PATH, pattern = "*.pgf"))==1 && length(list.files(CHIP_PATH, pattern = "*.mps"))==1)  &&
        !(length(list.files(CHIP_PATH, pattern = "*fasta"))==1 && length(list.files(CHIP_PATH, pattern = "*.cdf"))==1))
    stop("Array definition files don't exist!")
  }

  if (missing(SPECIES) || !grepl("_",SPECIES)) {
    SPECIES <- readline("Define your SPECIES with a latin name (i.e. Sus_scrofa): ")
    if (!grepl("_",SPECIES)) stop("SPECIES definition isn't correct!")
  }

  if (missing(NCPU) || is.na(NCPU) || NCPU<0) {
    NCPU <- as.numeric(readline("Number of threads to use:"))
    if (is.na(NCPU) || NCPU<0) stop("Wrong number of threads!")
  }

  if (missing(STRICT_MAPPING) || !is.logical(STRICT_MAPPING) || length(STRICT_MAPPING)!=1)
    STRICT_MAPPING <- readline("Apply strict mapping without SNP correction? y/[n]: ") == "y"

  if (missing(ASK_QUESTIONS) || !is.logical(ASK_QUESTIONS) || length(ASK_QUESTIONS)!=1)
    ASK_QUESTIONS <- readline("Answer questions automatically by a default? y/[n]: ") != "y"

  folder.in  <- pastef(WORK_PATH, PROJ_NAME, "in")
  folder.out <- pastef(WORK_PATH, PROJ_NAME, "out")
  folder.tmp <- pastef(WORK_PATH, PROJ_NAME, "tmp")
  pasteIN  <- function(...) pastef(folder.in, ...)
  pasteTMP <- function(...) pastef(folder.tmp, ...)
  pasteOUT <- function(...) pastef(folder.out, ...)

  return(as.list(environment()))
}
