library(readr)
library(devtools)
# just read the character string
test.data <- read_file("data-raw/wufc46h_reference.txt")
usethis::use_data(test.data)