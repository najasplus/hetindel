\name{iupacNotation}
\alias{iupacNotation}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{IUPAC Notation}
\description{
Returns the IUPAC notation of a vector of two characters. 
}
\usage{
iupacNotation(chars)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{chars}{A vector of two characters of "A", "C", "G" or "T".}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A character representing the IUPAC notation of one or two characters,
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Robert Deibel
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
chars <- c("A","G")
iupac <- iupacNotation(chars)
print(iupac)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ datagen }% use one of  RShowDoc("KEYWORDS")
