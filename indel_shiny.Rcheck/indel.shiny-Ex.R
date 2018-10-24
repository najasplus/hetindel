pkgname <- "indel.shiny"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('indel.shiny')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("indel.shiny-package")
### * indel.shiny-package

flush(stderr()); flush(stdout())

### Name: indel.shiny-package
### Title: Analyze heterozygous and homozygous indels
### Aliases: indel.shiny-package indel.shiny
### Keywords: package

### ** Examples

## Not run: 
##D start_indel_shiny()
## End(Not run)



cleanEx()
nameEx("start_indel_shiny")
### * start_indel_shiny

flush(stderr()); flush(stdout())

### Name: start_indel_shiny
### Title: Start Shiny App
### Aliases: start_indel_shiny
### Keywords: environment

### ** Examples

## Not run: 
##D start_indel_shiny()
## End(Not run)


### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
