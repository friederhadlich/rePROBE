#' Create a project directory
#' @param ARG_LIST list with user settings created by prepare_input_arguments.R
#'
#' @export
create_project <- function(ARG_LIST) {

  if (!any("PROJ_NAME" == ls(ARG_LIST)))  stop("Wrong input for create_project!")

  TMP <- pastef(ARG_LIST$WORK_PATH, ARG_LIST$PROJ_NAME)
  if (dir.exists(TMP)) {
    catn("Project already exists!")
    catn("Available databases:")
    print(list_available_dbs(ARG_LIST))
    catn("")
  } else {
    catn("Create project directory!")
    dir_create(TMP)
  }

  with(ARG_LIST, {
    if (!dir.exists(folder.in))  {
      dir_create(folder.in)
      if (dir.exists(folder.out)) dir_delete(folder.out);
      if (dir.exists(folder.tmp)) dir_delete(folder.tmp);
    }
    dir_create(c(folder.out, folder.tmp))
    save(ARG_LIST, file=pasteIN("ARG_LIST.Rda"))
  })

  return(TRUE)
}
