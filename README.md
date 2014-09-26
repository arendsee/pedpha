pedpha
======

Map protein intervals (e.g. domain positions) to their exon sources, calculate
exon phases.

exonstat
========

##input

One gff file which contains exon and CDS entries (it may contain other stuff as well)

##output

A tab-delimited file with the following columns:
 1. source
 1. transcript id
 1. exon number
 1. exon length
 1. exon 5' phase
 1. exon 3' phase
 1. exon strand
 1. exon start relative to gene
 1. absolute exon start

Phase is only defined if the exon boundary is inside a coding region.

exonphaser
==========

##input

 1. A tab-delimited file specifying the interval sequence id, interval id, start and stop postions (relative to the protein)
 1. A gff file containing both CDS and exon entries

##output

A tab-delimited file with the following columns:
 1. source 
 1. transcript id
 1. exon number
 1. exon length
 1. exon 5' phase
 1. exon 3' phase
 1. interval start relative to exon
 1. interval end relative to exon
 1. interval start relative to gene
 1. interval end relative to gene
 1. absolute start
 1. absolute end
