\name{introduceShift}
\alias{introduceShift}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Shift sequence}
\description{
Shifts a given sequence at specified occurrence.
}
\usage{
introduceShift(seqence, exact = TRUE, occurrence = "m")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{seqence}{
A vector of character values representing a DNA sequence.
}
  \item{exact}{
\code{logical} indicating whether the occurrence is to be at the exact place or if some randomisation should take place.
}
  \item{occurrence}{
One of \code{c("m","b","q1","q3")} indicating the middle beginning first quater and third quater of the sequence. Where the shift is to occur.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
The combined sequence in IUPAC notation with introduced shift.
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
sequence <- character(150)
sequence[] <- "A"
shiftedSeq <- introduceShift(sequence)
print(shiftedSeq)
shiftedSeq <- introduceShift(sequence,FALSE)
print(shiftedSeq)
shiftedSeq <- introduceShift(sequence,occurrence="b")
print(shiftedSeq)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ datagen}% use one of  RShowDoc("KEYWORDS")
