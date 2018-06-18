\name{generateRandBioString}
\alias{generateRandBioString}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Generate Random Bio-String
}
\description{
Generate a random DNA string with introduced shifts at a given place in IUPAC notation.
}
\usage{
generateRandBioString(length, dist = runif(1, 0, 4), exact = TRUE, occurrence = "m")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{length}{
Length of the string (before shift). One of \code{c(150,300,700)}.
}
  \item{dist}{
The distribution used to generate the sequence.
}
  \item{exact}{
Whether the shift is to be introduced exactly or with some randomisation.
}
  \item{occurrence}{
Occurrence of the shift. One of \code{c("m","b","q1","q3")} indicating the middle, beginning, first quater and third quater of the sequence.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A string in IUPAC notation.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
biostring <- generateRandBioString(150)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ datagen }% use one of  RShowDoc("KEYWORDS")
