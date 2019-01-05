\name{calculate.indel}
\alias{calculate.indel}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Calculates Indels for Base strings}
\description{
Calculates and prints the found indels for a given String of IUPAC bases
}
\usage{
calculate.indel(sequence, max.shift)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
	\item{sequence}{A String of IUPAC bases}
	\item{max.shift}{An integer stating the prefered maximum shift length. Is 15 by default}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Calculated indel Strings with marked indel positions
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
\dontrun{
calculate.indel("AAAGGGTTTTCCCGCT", 12)
}}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ manip }% use one of  RShowDoc("KEYWORDS")
