pkgname <- "seqGen"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('seqGen')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("generateRandBioString")
### * generateRandBioString

flush(stderr()); flush(stdout())

### Name: generateRandBioString
### Title: Generate Random Bio-String
### Aliases: generateRandBioString
### Keywords: datagen

### ** Examples

biostring <- generateRandBioString(150)



cleanEx()
nameEx("generateSeq")
### * generateSeq

flush(stderr()); flush(stdout())

### Name: generateSeq
### Title: Generate DNA sequence
### Aliases: generateSeq
### Keywords: datagen

### ** Examples

seq <- generateSeq(150)



cleanEx()
nameEx("introduceShift")
### * introduceShift

flush(stderr()); flush(stdout())

### Name: introduceShift
### Title: Shift sequence
### Aliases: introduceShift
### Keywords: datagen

### ** Examples

sequence <- character(150)
sequence[] <- "A"
shiftedSeq <- introduceShift(sequence)
print(shiftedSeq)
shiftedSeq <- introduceShift(sequence,FALSE)
print(shiftedSeq)
shiftedSeq <- introduceShift(sequence,occurrence="b")
print(shiftedSeq)




cleanEx()
nameEx("iupacNotation")
### * iupacNotation

flush(stderr()); flush(stdout())

### Name: iupacNotation
### Title: IUPAC Notation
### Aliases: iupacNotation
### Keywords: datagen

### ** Examples

chars <- c("A","G")
iupac <- iupacNotation(chars)
print(iupac)



cleanEx()
nameEx("randomBase")
### * randomBase

flush(stderr()); flush(stdout())

### Name: randomBase
### Title: Randomized Base
### Aliases: randomBase
### Keywords: datagen

### ** Examples

base <- randomBase()
print(base)



cleanEx()
nameEx("seqGen-package")
### * seqGen-package

flush(stderr()); flush(stdout())

### Name: seqGen-package
### Title: Generates random DNA-String with introduced shifts.
### Aliases: seqGen-package seqGen
### Keywords: package

### ** Examples

generateRandBioString(150)



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
