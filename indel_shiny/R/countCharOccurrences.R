countCharOccurrences <-
function(char, s) {
    s2 <- gsub(char,"",s)
    return (nchar(s) - nchar(s2))
}
