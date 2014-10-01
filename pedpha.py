#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import lib.gffreader as reader
import sys
import argparse
import collections

__version__ = "1.2.0"

class Intervals:
    def __init__(self, data, delimiter=None):
        self.intervals = self._read_data(data, delimiter)

    def _read_data(self, data, delimiter):
        out = collections.defaultdict(list)
        for line in data:
            row = line.split(delimiter)
            try:
                seqid, domid = row[0:2]
                start, stop = (int(s) for s in row[2:4])
                bounds = to_dna_interval(sorted([start, stop]))
                value = (domid, bounds)
                out[seqid].append(value)
                if start < 1 or stop < 1:
                    raise ValueError
            except IndexError:
                sys.exit("Each interval line must have 4 columns")
            except ValueError:
                sys.exit("Interval coordinants must be integers greater than 0")
        return(out)

    def get_bounds(self, ident):
        try:
            for val in self.intervals[ident]:
                yield val
        except KeyError:
            yield None

def parse(argv=None):
    parser = argparse.ArgumentParser(prog='pedpha')
    parser.add_argument(
        '-v', '--version',
        action='version',
        version='%(prog)s {}'.format(__version__))
    parser.add_argument(
        '-i', '--intervals',
        help="""File containing protein intervals
             (mRNA_ident, interval_ident, start, stop)
             where columns are delimited by DEL (whitespace by default).""",
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
        help="INTER file delimiter (defaults to whitespace)",
        metavar="DEL"
    )

    args = parser.parse_args(argv)

    return(args)


def to_dna_interval(x):
    '''
    x is an interval (x1, x2) on a protein
    Where x1 and x2 are integers greater than 0 and x2 >= x1
    '''
    x = [(i - 1)*3 + 1 for i in x]
    x[1] += 2
    return(x)

def get_overlap(x, y, minus=False):
    '''
    @param x: indices for an interval numbered relative to y
    @param y: a positive integer interval (greater than 0)
    @returns: overlap
    '''
    if minus:
        a, b = [y[1] - xi + 1 for xi in x]


        if a < y[0]:
            return (None, None)

        if b <= y[0]:
            b = y[0]

        return((b,a))

    else:
        a,b = [xi + y[0] - 1 for xi in x]

        if a > y[1]:
            return (None, None)

        if b >= y[1]:
            b = y[1]

        return((a,b))

def phaser(gff, intervals, delimiter=None):
    inter = Intervals(intervals, delimiter)
    for gene in reader.gff_reader(gff):
        for mrna in gene.mRNAs:
            mrna.calculate_phases()
            for domid, bounds in inter.get_bounds(mrna.ident):
                if not bounds:
                    continue
                total, cds_length = 0, 0
                for exon in mrna.exons:
                    if not exon.CDS:
                        continue

                    a,b = get_overlap(bounds, exon.CDS.bounds, bool(gene.strand == "-"))

                    if a and b:
                        if gene.strand == "+":
                            ca, cb = sorted(i - exon.CDS.bounds[0] + total + 1 for i in (a,b))
                        else:
                            ca, cb = sorted(exon.CDS.bounds[1] - i + total + 1 for i in (a,b))

                        yield (
                            domid,
                            mrna.ident,
                            exon.num,
                            gene.strand,
                            exon.bounds[0],
                            exon.bounds[1],
                            a, b,
                            ca, cb,
                            '%s-%s' % exon.phase
                            )

                    cds_length = exon.CDS.bounds[1] - exon.CDS.bounds[0] + 1
                    total += cds_length
                    bounds[0] = 1 if (a and b) else bounds[0] - cds_length
                    bounds[1] -= cds_length

                    if bounds[1] < 1:
                        break


if __name__ == '__main__':
    args = parse()

    gff = args.gff if args.gff else sys.stdin

    if not args.intervals:
        for gene in reader.gff_reader(gff):
            for line in gene.tostr():
                print(line)
    else:
        for row in phaser(gff, args.intervals):
            print("%s %s %s %s %d %d %d %d %d %d %s" % row)
