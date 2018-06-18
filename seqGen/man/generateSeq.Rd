\name{generateSeq}
\alias{generateSeq}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Generate DNA sequence 
}
\description{
Generate a DNA sequence of a specified length.
}
\usage{
generateSeq(length, dist = runif(1, 0, 4))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{length}{
  One of \code{c(150,300,700)}. Specifies the length of the generated sequence.
}
  \item{dist}{
The distribution by which the sequence is generated.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A randomly generated DNA sequence .
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
seq <- generateSeq(150)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ datagen }% use one of  RShowDoc("KEYWORDS")
