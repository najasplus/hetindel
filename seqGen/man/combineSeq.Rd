\name{combineSeq}
\alias{combineSeq}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Combine Sequences}
\description{
Given two sequences as vectors, combine those to IUPAC notation.
}
\usage{
combineSeq(seq1, seq2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{seq1}{
A DNA sequence as a character vector.
}
  \item{seq2}{
A DNA sequence as a character vector.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
The combined sequence in IUPAC notation.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Robert Deibel}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
a <- c("A","T","T","C","G")
b <- c("A","C","A","G","C")
comb <- combineSeq(a,b)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ datagen }% use one of  RShowDoc("KEYWORDS")