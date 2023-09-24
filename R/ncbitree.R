#' @name subset_clade
#' @title Subset the taxonomy tree from the root
#' 
#' @description Retrieves a clade from the taxonomy tree.
#' 
#' @param root.tax.id the tax_id of the clade root (integer).
#' @param nodes.dmp a nodes.dmp table.
#' 
#' @details This function retrieves a clade from the taxonomy tree listed in the
#' table \code{nodes.dmp} (see \code{\link{read_nodes_dmp}}). The clade starts at the single taxon specified in
#' \code{root_tax_id} and contains all descending branches.
#' 
#' If you want to subset the taxonomy tree based on leaf nodes, see \code{\link{subset_tree}}.
#' 
#' @return A subset of the table \code{nodes.dmp} containing the clade that descends from \code{root.tax.id}.
#' 
#' @author Lars Snipen.
#' 
#' @seealso \code{\link{read_nodes_dmp}}, \code{\link{subset_tree}}.
#' 
#' @examples
#' 
#' @importFrom dplyr %>% slice
#' 
#' @export subset_clade
#' 
subset_clade <- function(root.tax.id, nodes.dmp){
  idx <- which(nodes.dmp$tax_id == root.tax.id)
  if(length(idx) == 0) stop("The root.tax.id is not found in the supplied nodes.dmp")
  tax.id <- root.tax.id
  while(length(idx) > 0){
    idx <- unique(which(nodes.dmp$parent_tax_id %in% nodes.dmp$tax_id[idx] & nodes.dmp$parent_tax_id != nodes.dmp$tax_id))
    tax.id <- c(tax.id, nodes.dmp$tax_id[idx])
  }
  nodes.dmp %>% 
    filter(tax_id %in% unique(tax.id)) %>% 
    return()
}


#' @name subset_tree
#' @title Subset the taxonomy tree from leaves
#' 
#' @description Retrieves a sub-tree from the taxonomy tree.
#' 
#' @usage subset_tree(leaf_tax_id, nodes.dmp)
#' 
#' @param leaf.tax.id a vector of tax_id (integers) of the branch leaf/leaves.
#' @param nodes.dmp a nodes.dmp table.
#' 
#' @details This function retrieves a sub-tree by starting at the vector of specified
#' \code{leaf_tax_id} and collect all branches ending at these nodes. Thus, it will in general
#' not be a clade. If you want to subset an entire clade, see \code{\link{subset_clade}}.
#' 
#' @return A subset of the table \code{nodes.dmp} containing the tax_id's of the branches
#' ending at \code{leaf_tax_id}.
#' 
#' @author Lars Snipen.
#' 
#' @seealso \code{\link{read_nodes_dmp}}, \code{\link{subset_clade}}.
#' 
#' @examples
#' 
#' @importFrom dplyr %>% filter
#' 
#' @export subset_tree
#' 
subset_tree <- function(leaf.tax.id, nodes.dmp){
  idx <- unique(match(leaf.tax.id, nodes.dmp$tax_id))
  if(sum(is.na(idx)) > 0) stop("Not all leaf.tax.id are found in nodes.dmp!")
  tax.id <- unique(nodes.dmp$tax_id[idx])
  while(sum(nodes.dmp$tax_id[idx] != nodes.dmp$parent_tax_id[idx]) > 0){
    idx <- unique(match(nodes.dmp$parent_tax_id[idx], nodes.dmp$tax_id))
    if(sum(is.na(idx)) > 0) stop("Inconsistent nodes.dmp, some parent_tax_id are not among the tax_id!")
    tax.id <- c(nodes.dmp$tax_id[idx], tax.id)
  }
  nodes.dmp %>% 
    filter(tax_id %in% tax.id) %>% 
    return()
}




#' @name branch_retrieve
#' @title Retrieves branches
#' 
#' @description Retrieves branches in the taxonomy tree.
#' 
#' @usage branch_retrieve(leaf_tax_id, nodes.dmp)
#' 
#' @param leaf.tax.id a vector of tax_id (integers) of the branch leaf/leaves.
#' @param nodes.dmp a nodes.dmp table.
#' 
#' @details This function retrieves the branch from the taxonomy tree that ends at the \code{leaf.tax.id}
#' and starts at the root of tree defined in the table \code{nodes.dmp} (see \code{\link{read_nodes_dmp}}).
#' 
#' NB! Only the root node in nodes.dmp must have its own tax_id as parent_tax_id! This is the criterion for
#' ending a branch.
#' 
#' Multiple leaves may be given as argument, resulting in multiple branches retrieved.
#' 
#' @return A branch list, which is simply a list containing named vectors of tax_id integers, from the root to
#' the leaf of each branch. The name of each element is its rank. You may replace the tax_id integers by the 
#' taxon names using \code{\link{branch_taxid2name}}.
#' 
#' @author Lars Snipen.
#' 
#' @seealso \code{\link{branch_taxid2name}}, \code{\link{subset_clade}}.
#' 
#' @examples
#' 
#' @export branch_retrieve
#' 
branch_retrieve <- function(leaf.tax.id, nodes.dmp){
  idx <- match(leaf.tax.id, nodes.dmp$tax_id)
  idd <- which(is.na(idx))
  if(length(idd) > 0){
    stop("The tax_id ", paste(leaf.tax.id[idd], collapse = ","), " does not exist in nodes.dmp table\n")
  }
  cat("Pruning nodes.dmp")
  all.tax.id <- unique(leaf.tax.id)
  n.old <- length(all.tax.id)
  all.tax.id <- unique(c(all.tax.id, nodes.dmp$parent_tax_id[idx]))
  while(length(all.tax.id) > n.old){
    cat(".")
    n.old <- length(all.tax.id)
    idx <- match(all.tax.id, nodes.dmp$tax_id)
    all.tax.id <- unique(c(all.tax.id, nodes.dmp$parent_tax_id[idx]))
  }
  cat("keep", length(idx), "rows\n")
  nodes.dmp <- nodes.dmp[idx,]
  branch.lst <- lapply(leaf.tax.id, get_branch, nodes.dmp)
  return(branch.lst)
}
get_branch <- function(leaf.tax.id, nodes.dmp){   # a single leaf
  idx <- match(leaf.tax.id, nodes.dmp$tax_id)
  tax.id <- nodes.dmp$tax_id[idx]
  ranx <- nodes.dmp$rank[idx]
  while(nodes.dmp$tax_id[idx] != nodes.dmp$parent_tax_id[idx]){
    idx <- match(nodes.dmp$parent_tax_id[idx], nodes.dmp$tax_id)
    tax.id <- c(nodes.dmp$tax_id[idx], tax.id)
    ranx <- c(nodes.dmp$rank[idx], ranx)
  }
  names(tax.id) <- ranx
  return(tax.id)
}



#' @name branch_taxid2name
#' @aliases branch_taxid2name, branch_name2taxid
#' @title Replace tax_id with name
#' 
#' @description Converts vectors of tax_id to vectors of names, or vice versa.
#' 
#' @usage branch_taxid2name(branch.lst, names.dmp)
#' 
#' @param branch.lst a list of branches.
#' @param names.dmp a names.dmp table.
#' 
#' @details This function is used to convert the tax_id's to name_txt's in a branch-list, or vice versa. See 
#' \code{\link{branch_retrieve}} for more about branch-lists.
#' 
#' \code{branch_taxid2name}: The \code{names.dmp$name_txt} may contain many names for each tax_id, and the
#' argument \code{name.class} is used to select the type of name to use from the column \code{names.dmp$name_class}.
#' 
#' \code{branch_name2taxid}: The matching of names is exact (using \code{\link{match}}), which means only the first
#' occurrence of the name_txt is used. You are responsible for using a \code{names.dmp} with unique texts in
#' \code{names.dmp$name_txt}.
#' 
#' See \code{\link{read_names_dmp}} for more about \code{names.dmp} tables.
#' 
#' @return A list containing named vectors of either texts (\code{branch_taxid2name}) or integers 
#' (\code{branch_name2taxid}), from the root to the leaf of each branch.
#' The name of each element is its rank.
#' 
#' @author Lars Snipen.
#' 
#' @seealso 
#' 
#' @examples
#' 
#' @importFrom dplyr %>% filter
#' 
#' @export branch_taxid2name
#' @export branch_name2taxid
#' 
branch_taxid2name <- function(branch.lst, names.dmp, name.class = "scientific name"){
  all.unique.taxid <- unique(unlist(branch.lst))
  names.dmp <- names.dmp %>% 
    filter(name_class == name.class) %>% 
    filter(tax_id %in% all.unique.taxid)
  branch.lst <- lapply(branch.lst, get_taxon_name, names.dmp)
  return(branch.lst)
}
branch_name2taxid <- function(branch.lst, names.dmp){
  all.unique.names <- unique(unlist(branch.lst))
  names.dmp <- filter(names.dmp, name_txt %in% all.unique.names)
  branch.lst <- lapply(branch.lst, get_tax_id, names.dmp)
  return(branch.lst)
}

get_taxon_name <- function(branch, names.dmp){
  idx <- match(branch, names.dmp$tax_id)
  if(sum(is.na(idx)) > 0) stop("tax_id's of some branch is not in names.dmp$tax_id")
  new.branch <- names.dmp$name_txt[idx]
  names(new.branch) <- names(branch)
  return(new.branch)
}
get_tax_id <- function(branch, names.dmp){
  idx <- match(branch, names.dmp$name_txt)
  if(sum(is.na(idx)) > 0) stop("name_txt's of some branch is not in names.dmp$name_txt")
  new.branch <- names.dmp$tax_id[idx]
  names(new.branch) <- names(branch)
  return(new.branch)
}



#' @name branch_prune
#' @title Filtering ranks in branches
#' 
#' @description Keeps only specified ranks in specified ordering in all branches.
#' 
#' @param branch.lst a list of branches.
#' @param ranks texts specifying the ranks to keep and their ordering.
#' 
#' @details Branches in the taxonomy tree may have many different ranks or levels. This function is used
#' to prune them all to have the same set of ranks. Ranks may be missing in some branches.This function
#' will then fill in NA in these cells, ensuring all branches have the exact same ranks in the exact same
#' ordering. This is convenient for turning the list into a table, using \code{\link{as_tibble}}.
#' 
#' @return A new branch-list with pruned branches all containing the exact same ranks.
#' 
#' @author Lars Snipen.
#' 
#' @seealso 
#' 
#' @examples
#' 
#' @export branch_prune
#' 
branch_prune <- function(branch.lst, ranks = c("superkingdom", "phylum", "class", "order", "family", "genus", "species")){
  new.branch.lst <- lapply(branch.lst, function(x){
    idx <- match(ranks, names(x))
    y <- x[idx]
    names(y) <- ranks
    return(y)
  })
  return(new.branch.lst)
}



#' @name branch_list2table
#' @title Branch list na dtable conversion
#' @aliases branch_list2table branch_table2list
#' 
#' @description Turns a list of branches into a table with, or vice versa.
#' 
#' @param branch.lst a list of branches.
#' @param ranks texts specifying the ranks to keep and their ordering.
#' @param branch.tbl a table of branches.
#' 
#' @details Instead of having branches in a list, it is convenient in R to have data in tables. The
#' \code{branch_list2table} converts a list of branches into a table, with one row for each branch
#' and a column for each rank. Since branches may contain different ranks, they must first be pruned
#' to have the exact same ranks in the exact same order. This is done by \code{\link{branch_prune}}
#' inside this function. The ranks specified in \code{ranks} are the columns in the table returned here.
#' 
#' The function \code{branch_table2list} converts back again, i.e. from a table to a list.
#' 
#' @return The \code{branch_list2table} returns a tibble with one row for each branch and a column for
#' each rank. The cells contain the tax_id's or names in the \code{branch.lst}.
#' 
#' The \code{branch_table2list} returns a list where each element contains a named vector of tax_id or 
#' texts, and the name of each element is its rank in the taxonomy.
#' 
#' @author Lars Snipen.
#' 
#' @seealso 
#' 
#' @examples
#' 
#' @importFrom tibble as_tibble
#' @importFrom dplyr %>% 
#' 
#' @export branch_list2table
#' @export branch_table2list
#' 
branch_list2table <- function(branch.lst,
                              ranks = c("superkingdom", "phylum", "class", "order", "family", "genus", "species")){
  branch_prune(branch.lst, ranks = ranks) %>% 
    sapply(function(x){x}, simplify = T) %>% 
    t() %>% 
    as_tibble() -> branch.tbl
  return(branch.tbl)
}
branch_table2list <- function(branch.tbl){
  M <- as.matrix(branch.tbl)
  branch.lst <- lapply(1:nrow(M), function(i){M[i,]})
  return(branch.lst)
}



#' @name branch_list2qiime
#' @title Branch list to QIIME format
#' 
#' @description Creates QIIME formatted texts from a branch list.
#' 
#' @param branch.lst a list of branches.
#' 
#' @details A QIIME formatted branch is a text with the format 
#' 
#' \code{"d__domain;p__phylum;c__class;o__order;f__family;g__genus;s__species"}
#' 
#' where the words domain, phylum, etc are replaced by taxon names. This function converts
#' a list of N branches to N such texts.
#' 
#' NOTE 1: This is only meaningful if the \code{branch.lst} contains names instead of tax_id's, see
#' \code{\link{branch_taxid2name}}.
#' 
#' NOTE 2: The ranks are fixed to the ones listed above, where domain is the same as superkingdom in
#' the NCBI taxonomy.
#' 
#' @return The vector of texts, one for each element in \code{branch.lst}.
#' 
#' @author Lars Snipen.
#' 
#' @seealso 
#' 
#' @examples
#' 
#' @importFrom dplyr %>% 
#' 
#' @export branch_list2qiime
#' 
branch_list2qiime <- function(branch.lst){
  branch.lst %>% 
    branch_prune() %>% 
    branch_list2table() %>% 
    as.matrix() -> branch.mat
  prefix <- c("d__", "p__", "c__", "o__", "f__", "g__", "s__")
  txt <- apply(branch.mat, 1, function(x){
    x[is.na(x)] <- ""
    paste(paste0(prefix, x), collapse = ";")
  })
  return(txt)
}

