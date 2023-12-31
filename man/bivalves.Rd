% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{bivalves}
\alias{bivalves}
\title{Paleobiology Database occurrences of Pliocene fossil bivalves}
\format{
A data frame with 8095 rows and 9 variables:
\describe{
\item{genus}{Latin genus identification. Subgenera are not elevated.}
\item{paleolng, paleolat}{Coordinates of an occurrence, rotated to
its palaeogeographic location with the tectonic plate model of
\href{https://www.earthbyte.org/}{GPlates}}
\item{collection_no, reference_no}{Unique identifiers for the collection and published
reference containing the occurrence}
\item{environment}{One of 23 marine environment categories}
\item{max_ma, min_ma}{Bounds of the age estimate for an occurrence}
\item{accepted_name}{Original identification, including subgenus and
species epithet if applicable, according to the latest PBDB
accepted taxonomy at time of download}
}
}
\source{
\url{https://paleobiodb.org/}
}
\usage{
bivalves
}
\description{
A dataset containing the (palaeo)coordinates and genus identifications of
8,000 marine bivalves from the Pliocene (ca. 5.3-2.6 Ma). Records with
uncertain or unaccepted taxonomic names, non-marine palaeo-environments,
or missing coordinates are excluded from the original download
(24 June 2022).
}
\keyword{datasets}
