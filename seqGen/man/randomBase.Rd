\name{randomBase}
\alias{randomBase}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Randomized Base}
\description{
Returns a random Base with specified distribution. If not further stated the uniform distribution is used.
}
\usage{
randomBase(dist = runif(1, 0, 4))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dist}{The distribution the base is derived from. Lower bound should be 0 and upper bound should always be 4.}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A character of either "A", "C", "G" or "T".
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
base <- randomBase()
print(base)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ datagen }% use one of  RShowDoc("KEYWORDS")

