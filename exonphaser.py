#!/usr/bin/env python3

import lib.gffreader as reader
import sys
import argparse

def parse(argv=None):
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-i', '--intervals',
        help="File containing protein intervals (mRNA_ident, interval_ident, start, stop)",
        metavar="INTER",
        type=argparse.FileType('r')
    )
    parser.add_argument(
        '-g', '--gff',
        help="GFF file formatted according to JGI standards (like phytozome)",
        metavar="GFF",
        type=argparse.FileType('r')
    )
    parser.add_argument(
        '-d', '--delimiter',
        help="INTER file delimiter (defaults to whitespace)"
    )
    parser.add_argument(
        '-t', '--type',
        help="Interval type",
        choices=['prot', 'nucl']
    )

    args = parser.parse_args(argv)

    if not args.type:
        sys.exit("Please provide a type (e.g. -t prot)")

    return(args)

class Intervals:
    def __init__(self, data, delimiter=None):
        self.intervals = self._read_data(data, delimiter)

    def _read_data(self, data, delimiter):
        out = dict()
        for line in data:
            row = line.split(delimiter)
            try:
                out[row[0]] = row[1:4]
                out[row[0]][1:3] = [int(s) for s in out[row[0]][1:3]]
            except IndexError:
                sys.exit("Each interval line must have 4 columns")
        return(out)

    def get_bounds(self, ident):
        try:
            return(self.intervals[ident][1:3])
        except KeyError:
            return(None)

    def get_ident(self, ident):
        try:
            return(self.intervals[ident][0])
        except KeyError:
            return(None)

def map_protein(exon, bounds):
    bounds = (bounds[0] * 3, bounds[1] * 3)
    return([1,1])


def map_nucleotide(exon, bounds):
    return([2,2])

if __name__ == '__main__':
    args = parse()

    if args.type == "prot":
        mapper = map_protein
    else:
        mapper = map_nucleotide

    template = "%s %s %s %s %d %d %d %d %d %d"
    inter = Intervals(args.intervals, args.delimiter)
    for gene in reader.gff_reader(args.gff):
        for mrna in gene.mRNAs:
            bounds = inter.get_bounds(mrna.ident)
            if not bounds:
                continue
            for exon in mrna.exons:
                exon_bounds = mapper(exon, bounds)
                if exon_bounds:
                    print(template %
                        tuple([
                            inter.get_ident(mrna.ident),
                            mrna.ident,
                            exon.num,
                            gene.strand] +
                            bounds +
                            exon.bounds +
                            exon_bounds
                        ))
