#' @name kraken2_read_table
#' @title Reads raw \code{kraken2} output
#' 
#' @description Reads the raw output table from \code{kraken2}.
#' 
#' @usage kraken2_read_table(krk.file, filter = TRUE)
#' 
#' @param krk.file name of file with raw \code{kraken2} output.
#' @param filter. returns only classified reads if TRUE.
#' 
#' @details The raw output from \code{kraken2} is a huge table with one row for each
#' read classified.
#' 
#' Note that in most cases only the report, produced by running \code{kraken2} with
#' the \code{--report} otpion, is what you want. Then use the function 
#' \code{\link{kraken2_read_report}}.
#' 
#' @return A data.frame with the kraken2 results for each read. It has the columns:
#' \itemize{
#'   \item status. C for classified reads, U for unclassified (if \code{filter = TRUE} only C).
#'   \item read_id. Short text identifying each read.
#'   \item name. Name of taxon, if classified.
#'   \item length. Length of read.
#'   \item tax_count. The k-mer counts for the read.
#' }
#'
#' @author Lars Snipen.
#' 
#' @seealso \code{\link{kraken2_read_report}}.
#' 
#' @examples
#' 
#' @importFrom data.table fread
#' @importFrom dplyr %>% filter
#' 
#' @export kraken2_read_table
#' 
kraken2_read_table <- function(krk.file, filter = TRUE){
  tbl <- fread(file = normalizePath(krk.file),
               sep = "\t",
               data.table = F,
               quote = "",
               col.names = c("status", "read_id", "name", "length", "tax_count"))
  if(filter) tbl %>% filter(status == "C") -> tbl
  return(tbl)
}


#' @name kraken2_read_report
#' @title Reads \code{kraken2} reports
#' 
#' @description Reads the report output from \code{kraken2}.
#' 
#' @usage kraken2_read_report(report.file)
#' 
#' @param report.file name of file with \code{kraken2} report output.
#' 
#' @details The report output from kraken2 is a table with one row for each
#' taxon detected by kraken2. Run \code{kraken2} using the \code{--report}
#' option to get this convenient summary table.
#' 
#' Note that the software \code{bracken} may also output a report in this format.
#' 
#' @return A data.frame with the \code{kraken2} results for each taxon. It has the columns:
#' \itemize{
#'   \item percent. The percentage of the reads classified to this taxon.
#'   \item clade_count. The total number of reads classified to the clade originating from this taxon.
#'   \item tax_count. The number of reads assigned directly to this taxon.
#'   \item rank. The rank of this taxon in the taxonomy.
#'   \item tax_id. The taxonomy identifier of this taxon in the taxonomy.
#'   \item name. The name of this taxon in the taxonomy.
#' }
#'
#' @author Lars Snipen.
#' 
#' @seealso \code{\link{kraken2_read_table}}, \code{\link{bracken_read_report}}.
#' 
#' @examples
#' 
#' @importFrom data.table fread
#' 
#' @export kraken2_read_report
#' 
kraken2_read_report <- function(report.file){
  fread(file = normalizePath(report.file),
        sep = "\t",
        data.table = F,
        col.names = c("percent", "clade_count", "tax_count", "rank", "tax_id", "name")) -> tbl
  return(tbl)
}


#' @name bracken_read_report
#' @title Reads bracken reports
#' 
#' @description Reads the report output from \code{bracken}.
#' 
#' @usage bracken_read_report(report.file)
#' 
#' @param report.file name of file with bracken report output.
#' 
#' @details The report output from bracken is a table with one row for each
#' taxon with assigned reads. This is the file produced by \code{bracken} using
#' the \code{-o} option.
#' 
#' Note that by the \code{-w} option \code{bracken} will output a report in
#' \code{kraken2} format, see \code{\link{kraken2_read_report}}.
#' 
#' @return A data.frame with the \code{bracken} results for each taxon. It has the columns:
#' \itemize{
#'   \item name. The name of this taxon in the taxonomy.
#'   \item tax_id. The taxonomy identifier of this taxon in the taxonomy.
#'   \item rank. The rank of this taxon in the taxonomy.
#'   \item tax_count_krk. The number of reads classified to this taxon by \code{kraken2}.
#'   \item added_brk. The additional reads assigned by \code{bracken}.
#'   \item tax_count. The final number of reads assigned to this taxon.
#'   \item fraction. The fraction of the total number of reads assigned to this taxon.
#' }
#'
#' @author Lars Snipen.
#' 
#' @seealso \code{\link{kraken2_read_report}}.
#' 
#' @examples
#' 
#' @importFrom data.table fread
#' 
#' @export bracken_read_report
#' 
bracken_read_report <- function(report.file){
  fread(file = normalizePath(report.file),
        sep = "\t",
        data.table = F,
        col.names = c("name", "tax_id", "rank", "tax_count_krk", "added_brk", "tax_count", "fraction")) -> tbl
  return(tbl)
}





reassign_counts <- function(report.tbl, target.rank = "S", verbose = TRUE){
  rank.order <- c("R", "D", "K", "P", "C", "O", "F", "G", "S")
  if(verbose) cat("reassign_counts:\n")
  report.tbl %>% 
    filter(rank != "U") -> tbl
  urank <- sort(unique(tbl$rank))
  if(!(target.rank %in% urank)){
    tibble(percent = NULL,
           clade_count = NULL,
           tax_count = NULL,
           rank = NULL,
           tax_id = NULL,
           name = NULL,
           new_count = NULL) -> tbl
    return(tbl)
  }
  levs <- unlist(sapply(rank.order, function(r){urank[str_detect(urank, r)]}))
  target.int <- which(levs == target.rank)
  tbl %>% 
    mutate(rank = factor(rank, levels = levs)) %>% 
    mutate(rank_int = as.integer(rank)) %>% 
    filter(rank_int <= target.int) -> tbl
  levs <- levs[1:target.int]
  N <- length(levs)
  bss <- c(0, which(diff(tbl$rank_int) <= 0), nrow(tbl))
  if(verbose) cat("   kraken2 report contains", N, "taxonomy ranks and", length(bss)-1, "branches...\n")
  tax.mat <- cnt.mat <- matrix(NA, nrow = length(bss)-1, ncol = N)
  colnames(tax.mat) <- colnames(cnt.mat) <- levs
  rownames(tax.mat) <- rownames(cnt.mat) <- tbl$tax_id[bss[-1]]
  tax <- cnt <- rep(NA, N)
  for(i in 2:length(bss)){
    rr <- (bss[i-1] + 1):bss[i]
    iii <- tbl$rank_int[rr]
    if(min(iii) > 1){
      tax <- c(tax[1:(min(iii) - 1)], rep(NA, N - min(iii) + 1))
      cnt <- c(cnt[1:(min(iii) - 1)], rep(NA, N - min(iii) + 1))
    }
    tax[iii] <- tbl$tax_id[rr]
    cnt[iii] <- tbl$clade_count[rr]
    tax.mat[i-1,] <- tax
    cnt.mat[i-1,] <- cnt
  }
  if(verbose) cat("   distributing clade counts down the taxonomy...\n")
  for(i in 1:(N-1)){
    uparent <- unique(tax.mat[,i])
    uparent <- uparent[which(!is.na(uparent))]
    for(j in 1:length(uparent)){
      p.rr <- which(tax.mat[,i] == uparent[j])
      c.cc <- rep(NA, length(p.rr))
      cc <- i + 1
      nidx <- which(is.na(c.cc))
      while(sum(nidx) > 0 & cc <= ncol(cnt.mat)){
        idx <- which(!is.na(tax.mat[p.rr,cc]))
        c.cc[idx] <- cc
        nidx <- which(is.na(c.cc))
        cc <- cc + 1
      }
      M <- matrix(c(p.rr, c.cc), ncol = 2, byrow = F) # p.rr are the rows and c.cc the columns of the children
      children <- factor(tax.mat[M], levels = unique(tax.mat[M]))
      u.rr <- which(!duplicated(children))
      MM <- matrix(c(p.rr[u.rr], c.cc[u.rr]), ncol = 2, byrow = F)
      up <- cnt.mat[MM] / sum(cnt.mat[MM])
      prob <- up[as.integer(children)]
      d <- cnt.mat[p.rr[1], i] - sum(cnt.mat[MM])
      cnt.mat[M] <- cnt.mat[M] + d * prob
    }
  }
  cnt.mat[,N,drop = F] %>% 
    as_tibble(rownames = "tax_id") %>% 
    magrittr::set_colnames(c("tax_id", "value")) %>% 
    mutate(tax_id = as.integer(tax_id)) -> cnt.tbl
  report.tbl %>% 
    filter(rank == target.rank) %>% 
    left_join(cnt.tbl, by = "tax_id") %>% 
    mutate(percent = 100 * value / sum(value)) %>% 
    select(percent, clade_count, tax_count, rank, tax_id, name, new_count = value) -> tbl
  if(exists("read_id", where = report.tbl)){
    tbl %>% mutate(read_id = report.tbl$read_id[1]) -> tbl
  }
  return(tbl)
}


kraken2_taxonomy <- function(inspect.file){
  rank.order <- c("R", "D", "K", "P", "C", "O", "F", "G", "S")
  read_kraken2_report(inspect.file) %>% 
    filter(rank != "U") -> insp.tbl
  urank <- sort(unique(insp.tbl$rank))
  levs <- unlist(sapply(rank.order, function(r){urank[str_detect(urank, r)]}))
  insp.tbl %>% 
    mutate(rank = factor(rank, levels = levs)) %>% 
    mutate(rank_int = as.integer(rank)) %>% 
    mutate(rank = as.character(rank)) %>% 
    mutate(percent = 0, tax_count = 0, clade_count = 0) -> insp.tbl
  bss <- c(0, which(diff(insp.tbl$rank_int) <= 0), nrow(insp.tbl))
  branch.mat <- matrix(NA, nrow = length(bss)-1, ncol = length(levs))
  colnames(branch.mat) <- levs
  rownames(branch.mat) <- insp.tbl$tax_id[bss[-1]]
  tax <- rep(NA, length(levs))
  for(i in 2:length(bss)){
    rr <- (bss[i-1] + 1):bss[i]
    iii <- insp.tbl$rank_int[rr]
    if(min(iii) > 1){
      tax <- c(tax[1:(min(iii) - 1)], rep(NA, length(levs) - min(iii) + 1))
    }
    tax[iii] <- insp.tbl$tax_id[rr]
    branch.mat[i-1,] <- tax
  }
  attr(insp.tbl, "branch.mat") <- branch.mat
  return(insp.tbl)
}


krk2kmer <- function(krk.tbl){
  krk.tbl %>% 
    select(read_id, tax_count) %>% 
    mutate(tax_count = str_remove(tax_count, "\\|:\\| ")) %>% 
    mutate(tax_id = str_remove_all(tax_count, ":[0-9]+")) %>% 
    mutate(count = str_remove_all(tax_count, "[0-9]+:")) %>% 
    select(-tax_count) %>% 
    mutate(tax_id = str_split(tax_id, pattern = " ")) %>% 
    mutate(count = str_split(count, pattern = " ")) %>% 
    unnest(c(tax_id, count)) %>% 
    mutate(tax_id = as.integer(tax_id), count = as.integer(count)) %>% 
    group_by(read_id, tax_id) %>% 
    summarise(tax_count = sum(count)) %>% 
    ungroup() -> kmer.tbl
  return(kmer.tbl)
}


kmer_report <- function(kmer.tbl, taxonomy.tbl){
  kmer.tbl %>% 
    filter(tax_id == 0) -> unclas.tbl
  kmer.tbl %>% 
    filter(tax_id != 0) -> kmr.tbl
  
  utxm <- numeric(0)
  branch.mat <- attr(taxonomy.tbl, "branch.mat")
  for(i in nrow(kmr.tbl):1){
    if(!(kmr.tbl$tax_id[i] %in% utxm)){
      M <- which(branch.mat == kmr.tbl$tax_id[i], arr.ind = T)
      txm <- unique(as.integer(branch.mat[M[,1], 1:M[1,2]]))
      utxm <- c(utxm, txm)
    }
  }
  utxm <- unique(utxm[!is.na(utxm)])
  taxonomy.tbl %>% 
    filter(tax_id %in% utxm) -> rep.tbl
  idx <- match(rep.tbl$tax_id, kmr.tbl$tax_id)
  idd <- which(!is.na(idx))
  rep.tbl$tax_count[idd] <- kmr.tbl$tax_count[idx[idd]]
  bss <- c(0, which(diff(rep.tbl$rank_int) <= 0), nrow(rep.tbl))
  
  t_count <- rep.tbl$tax_count
  for(i in length(bss):2){
    M <- which(branch.mat == rep.tbl$tax_id[bss[i]], arr.ind = T)
    bm <- branch.mat[M[1,1], 1:M[1,2]]
    bm <- bm[!is.na(bm)]
    rr <- match(bm, rep.tbl$tax_id)
    cc <- rev(cumsum(rev(t_count[rr])))
    rep.tbl$clade_count[rr] <- rep.tbl$clade_count[rr] + cc
    t_count[rr] <- 0
  }
  
  N.kmers <- rep.tbl$clade_count[1]
  if(nrow(unclas.tbl) > 0){
    N.kmers <- N.kmers + unclas.tbl$tax_count
    tibble(percent = 0,
           clade_count = unclas.tbl$tax_count,
           tax_count = unclas.tbl$tax_count,
           rank = "U",
           tax_id = unclas.tbl$tax_id,
           name = "unclassified",
           rank_int = 0) %>% 
      bind_rows(rep.tbl) -> rep.tbl
  }
  rep.tbl %>% 
    mutate(percent = clade_count / N.kmers) %>% 
    select(-rank_int) %>% 
    mutate(read_id = kmr.tbl$read_id[1]) -> rep.tbl
  return(rep.tbl)
}





# taxonomyTable <- function(names.dmp.file, nodes.dmp.file){
#   require(data.table)
#   fread(file = names.dmp.file, sep = "\t", header = F, data.table = F, quote = "") %>% 
#     filter(V7 == "scientific name") %>% 
#     select(tax_id = V1, name_txt = V3) -> names.tbl
#   fread(file = nodes.dmp.file, sep = "\t", header = F, data.table = F, quote = "") %>% 
#     select(tax_id = V1, parent_tax_id = V3, rank = V5) %>% 
#     left_join(names.tbl, by = "tax_id") -> tax.tbl
#   return(tax.tbl)
# }

# branchFinder <- function(tax.id, taxonomy.tbl){
#   idx <- which(tax.tbl$tax_id == tax.id)
#   if(length(idx) == 0) stop("tax_id ", tax.id, " is not in taxonomy.tbl")
#   branch <- tax.tbl$tax_id[idx]
#   while(tax.tbl$tax_id[idx] != 1){
#     idx <- which(tax.tbl$tax_id == tax.tbl$parent_tax_id[idx])
#     branch <- c(tax.tbl$tax_id[idx], branch)
#   }
#   return(branch)
# }

# rankRename <- function(tax.vec){
#   from <- c("superkingdom", "phylum", "class", "order", "family", "genus", "species")
#   to   <- c("D",            "P",      "C",     "O",     "F",      "G",     "S")
#   nam <- names(tax.vec)
#   nam[1] <- "R"
#   if(length(nam) > 1){
#     rnk <- "R"
#     num <- 1
#     for(i in 2:length(nam)){
#       idx <- match(nam[i], from)
#       if(is.na(idx)){
#         nam[i] <- str_c(rnk, num)
#         num <- num + 1
#       } else {
#         nam[i] <- rnk <- to[idx]
#         num <- 1
#       }
#     }
#   }
#   names(tax.vec) <- nam
#   return(tax.vec)
# }

# taxonomyGraph <- function(tax.tbl){
#   require(igraph)
#   M <- matrix(c(as.character(tax.tbl$parent_tax_id), as.character(tax.tbl$tax_id)), ncol = 2, byrow = F)
#   graph.edgelist(M, directed = T) %>% 
#     set_vertex_attr(name = "rank", index = M[,2], value = tax.tbl$rank) %>% 
#     set_vertex_attr(name = "name_txt", index = M[,2], value = tax.tbl$name_txt) -> graf
#   return(graf)
# }






# kraken2_run <- function(in.files, dbase, krk.file = "krk.txt",
#                     report.file = NULL, threads = 1, confidence = 0.0){
#   report.txt <- ifelse(is.null(report.file), "", str_c("--report ", report.file))
#   in.files <- normalizePath(in.files)
#   if(length(in.files) == 1){
#     paired.txt <- ""
#     query.txt <- in.files
#   } else {
#     paired.txt <- "--paired"
#     query.txt <- str_c(in.files, collapse = " ")
#   }
#   cmd <- paste("kraken2",
#                "--threads", threads,
#                "--use-names",
#                "--confidence", confidence,
#                "--db", dbase,
#                paired.txt,
#                query.txt,
#                report.txt,
#                ">", krk.file)
#   system(cmd)
# }


# kraken2_inspect <- function(dbase, report.file, threads = 1){
#   cmd <- paste("kraken2-inspect",
#                "--threads", threads,
#                "--db", dbase,
#                "--report-zero-counts",
#                ">", report.file)
#   system(cmd)
# }
