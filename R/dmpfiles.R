#' @name read_names_dmp
#' @title Read and write the names.dmp
#' @aliases read_names_dmp write_names_dmp
#' 
#' @description Reads and writes the file names.dmp of the NCBI Taxonomy
#' 
#' @usage read_names_dmp(filename)
#' write_names_dmp(names.dmp, filename)
#' 
#' @param filename name of file to be read or written to (text).
#' @param names.dmp a names.dmp table (see details below).
#' 
#' @details The file pair \code{names.dmp} and \code{nodes.dmp} describe a taxonomy tree. 
#' The \code{read_names_dmp} reads a file formatted as the \code{names.dmp} file from the NCBI Taxonomy
#' database (https://www.ncbi.nlm.nih.gov/taxonomy/). This is represented as a \code{tibble} in R.
#' 
#' The \code{write_names_dmp} will write a table with the proper columns (see below) to a file, adding 
#' the separators of the NCBI format.
#' 
#' @return The \code{read_names_dmp} returns a \code{tibble} with the columns: \code{tax_id} (integers), \code{name_txt} (text),
#' \code{unique_name} (text) and \code{name_class} (text). 
#'
#' 
#' @author Lars Snipen.
#' 
#' @seealso \code{\link{read_nodes_dmp}}.
#' 
#' @examples
#' 
#' @importFrom tibble as_tibble
#' @importFrom stringr str_split
#' @importFrom dplyr %>% select mutate
#' 
#' @export read_names_dmp
#' @export write_names_dmp
#' 
read_names_dmp <- function(filename){
  filename <- normalizePath(filename)
  if(file.exists(filename)){
    readLines(filename) %>%
      str_split(pattern = "\t", simplify = T) %>%
      as_tibble(.name_repair = "universal") %>%
      select(c("tax_id" = 1, "name_txt" = 3, "unique_name" = 5, "name_class" = 7)) %>%
      mutate(tax_id = as.integer(tax_id)) %>% 
      return()
  } else {
    stop("Cannot find ", filename, ", please correct path and/or file name")
  }
}
write_names_dmp <- function(names.dmp, filename){
  paste(as.character(names.dmp$tax_id),
        as.character(names.dmp$name_txt),
        as.character(names.dmp$unique_name),
        as.character(names.dmp$name_class), sep = "\t|\t") %>% 
    paste0("\t|") %>% 
    writeLines(con = filename)
  return(filename)
}



#' @name read_nodes_dmp
#' @title Read and write the nodes.dmp
#' @aliases read_nodes_dmp write_nodes_dmp
#' 
#' @description Reads and writes the file nodes.dmp of the NCBI Taxonomy
#' 
#' @usage read_nodes_dmp(filename)
#' write_nodes_dmp(names.dmp, filename)
#' 
#' @param filename name of file to be read or written to.
#' @param nodes.dmp a nodes.dmp table (see details below).
#' 
#' @details The file pair \code{names.dmp} and \code{nodes.dmp} describe a taxonomy tree.
#' The \code{read_nodes_dmp} reads a file formatted as the \code{nodes.dmp} file from the NCBI Taxonomy
#' database (https://www.ncbi.nlm.nih.gov/taxonomy/). This is represented as a \code{tibble} in R.
#' 
#' The \code{write_nodes_dmp} will write a table with the proper columns (see below) to a file, adding 
#' the separators of the NCBI format.
#' 
#' The \code{nodes.dmp} table downloaded from NCBI will contain many columns, but only the first 3 of them
#' are relevant for parsing the taxonomy tree. Only these first three columns are read and used by these functions,
#' additional columns are ignored. 
#' 
#' @return The \code{read_nodes_dmp} returns a \code{tibble} with the columns: \code{tax_id} (integers),
#' \code{parent_tax_id} (integers) and \code{rank} (text). 
#'
#' 
#' @author Lars Snipen.
#' 
#' @seealso \code{\link{read_nodes_dmp}}.
#' 
#' @examples
#' 
#' @importFrom tibble as_tibble
#' @importFrom stringr str_split
#' @importFrom dplyr %>% select mutate
#' 
#' @export read_nodes_dmp
#' @export write_nodes_dmp
#' 
read_nodes_dmp <- function(filename){
  filename <- normalizePath(filename)
  if(file.exists(filename)){
    readLines(filename) %>%
      str_split(pattern = "\t", simplify = T) %>%
      as_tibble(.name_repair = "universal") %>% 
      select(c("tax_id" = 1, "parent_tax_id" = 3, "rank" = 5)) %>% 
      mutate(tax_id = as.integer(tax_id), parent_tax_id = as.integer(parent_tax_id)) %>% 
      return()
  } else {
    stop("Cannot find ", filename, ", please correct path and/or file name")
  }
}
write_nodes_dmp <- function(nodes.dmp, filename){
  paste(as.character(nodes.dmp$tax_id),
        as.character(nodes.dmp$parent_tax_id),
        as.character(nodes.dmp$rank), 
        rep("-", nrow(nodes.dmp)), sep = "\t|\t") %>% 
  paste0("\t|") %>% 
    writeLines(con = filename)
  return(filename)
}
