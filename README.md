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

formats
=======

GFF format is very standardized, so I've decided to error on the side of dying
if anything goes wrong. Currently pedpha works only with gff files that are
formatted exactly like the JGI phytozome/metazome files.

 1. There must be the following structures:
```
 gene -
      |- mRNA -
              | - exon
              | - CDS
              | - exon
              | - CDS
      |- mRNA -
              | - exon
              | - CDS
```
 1. col4 < col5 (i.e. start < stop) even on minus strand
 1. exons are ordered 5' to 3', biological order. I.e. minus and plus strands are in opposite order
 1. element identifiers are extracted from /ID=([^;]+)/ patterns in the 9th column

Example (TABs converted to spaces for ease of viewing):
```
 Chr1  phytozomev10  gene             5928  8737  .  -  .  ID=AT1G01020.TAIR10;Name=AT1G01020
 Chr1  phytozomev10  mRNA             5928  8737  .  -  .  ID=AT1G01020.1.TAIR10;Name=AT1G01020.1;pacid=19655142;longest=1;Parent=AT1G01020.TAIR10
 Chr1  phytozomev10  exon             8571  8737  .  -  .  ID=AT1G01020.1.TAIR10.exon.1;Parent=AT1G01020.1.TAIR10;pacid=19655142
 Chr1  phytozomev10  CDS              8571  8666  .  -  0  ID=AT1G01020.1.TAIR10.CDS.1;Parent=AT1G01020.1.TAIR10;pacid=19655142
 Chr1  phytozomev10  exon             8417  8464  .  -  .  ID=AT1G01020.1.TAIR10.exon.2;Parent=AT1G01020.1.TAIR10;pacid=19655142
 Chr1  phytozomev10  CDS              8417  8464  .  -  0  ID=AT1G01020.1.TAIR10.CDS.2;Parent=AT1G01020.1.TAIR10;pacid=19655142
 Chr1  phytozomev10  exon             8236  8325  .  -  .  ID=AT1G01020.1.TAIR10.exon.3;Parent=AT1G01020.1.TAIR10;pacid=19655142
 Chr1  phytozomev10  CDS              8236  8325  .  -  0  ID=AT1G01020.1.TAIR10.CDS.3;Parent=AT1G01020.1.TAIR10;pacid=19655142
 Chr1  phytozomev10  exon             7942  7987  .  -  .  ID=AT1G01020.1.TAIR10.exon.4;Parent=AT1G01020.1.TAIR10;pacid=19655142
 Chr1  phytozomev10  CDS              7942  7987  .  -  0  ID=AT1G01020.1.TAIR10.CDS.4;Parent=AT1G01020.1.TAIR10;pacid=19655142
 Chr1  phytozomev10  exon             7762  7835  .  -  .  ID=AT1G01020.1.TAIR10.exon.5;Parent=AT1G01020.1.TAIR10;pacid=19655142
 Chr1  phytozomev10  mRNA             6790  8737  .  -  .  ID=AT1G01020.2.TAIR10;Name=AT1G01020.2;pacid=19655143;longest=0;Parent=AT1G01020.TAIR10
 Chr1  phytozomev10  exon             8571  8737  .  -  .  ID=AT1G01020.2.TAIR10.exon.1;Parent=AT1G01020.2.TAIR10;pacid=19655143
 Chr1  phytozomev10  CDS              8571  8666  .  -  0  ID=AT1G01020.2.TAIR10.CDS.1;Parent=AT1G01020.2.TAIR10;pacid=19655143
 Chr1  phytozomev10  exon             8417  8464  .  -  .  ID=AT1G01020.2.TAIR10.exon.2;Parent=AT1G01020.2.TAIR10;pacid=19655143
 Chr1  phytozomev10  CDS              8417  8464  .  -  0  ID=AT1G01020.2.TAIR10.CDS.2;Parent=AT1G01020.2.TAIR10;pacid=19655143
 Chr1  phytozomev10  exon             8236  8325  .  -  .  ID=AT1G01020.2.TAIR10.exon.3;Parent=AT1G01020.2.TAIR10;pacid=19655143
```
