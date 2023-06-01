data(small.16S)
if (system("which makeblastdb")){
  dbase <- blastDbase16S("test", small.16S$Sequence, word(small.16S$Header, 2, 2))
  reads <- str_sub(small.16S$Sequence, 100, 550)
  tbl <- blastClassify16S(reads, dbase) |>
    bind_cols(small.16S) 
  
  
  test_that("BLAST classify is working", {
    expect_equal(tbl$Taxon.hat, gsub(".* (.*)", "\\1", tbl$Header))
})
} else{
  warning("Not testing BLAST; makeblastdb not found in path")
}