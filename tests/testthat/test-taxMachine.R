test_that("taxMachine detects E. coli", {
  expect_equal(taxMachine(sequence="ACTGGGCGTAAAGCGCACGCAGGCGGTTTGTTAAGTCAGATGTGAAATCCCCGGGCTCAACCTGG", model.on.disk="tax.model")$Genus, 
               "Escherichia|Shigella")
})
