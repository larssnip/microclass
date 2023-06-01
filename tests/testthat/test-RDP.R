data(small.16S)
genus <- sapply(strsplit(small.16S$Header,split=" "),function(x){x[2]})

rdp <- rdpTrain(small.16S$Sequence[seq(1,71,2)], genus[seq(1,71,2)])  # training step
predicted <- rdpClassify(small.16S$Sequence[seq(2,71,2)], rdp)        # classification step
cat( "Number of errors:", sum(predicted != genus[seq(2,71,2)]) )

test_that("RDP train/classify is working", {
  expect_equal(predicted, genus[seq(2,71,2)])
})