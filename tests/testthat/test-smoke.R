test_that("parse_imputeinfofile reads modernized data.table correctly", {
  # SETUP: Create a fake tiny file so we don't depend on the cluster
  fake_file <- tempfile()
  write.table(data.frame(
    chrom = c(1, 2), legend = "A", map = "B", hap = "C",
    start = 1, end = 100, is_par = c(0, 1)
  ), fake_file, row.names = FALSE, col.names = FALSE, sep = "\t")

  # EXECUTE: Call your function
  result <- parse_imputeinfofile(fake_file, is.male = TRUE)

  # ASSERT: Check basic facts
  expect_s3_class(result, "data.table")
  expect_equal(nrow(result), 1) # Should be 1 because is.male filters is_par == 1
})
