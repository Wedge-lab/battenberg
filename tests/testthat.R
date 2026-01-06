library(testthat)
library(Battenberg) # Load your library

# This line tells R to look into the tests/testthat/ folder
# and run every file that starts with "test-"
test_check("Battenberg")
