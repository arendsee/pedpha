#!/usr/bin/env python3

import lib.gffreader as reader
import sys
import argparse

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
            return(sorted(3*x for x in self.intervals[ident][1:3]))
        except KeyError:
            return(None)

    def get_ident(self, ident):
        try:
            return(self.intervals[ident][0])
        except KeyError:
            return(None)


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

    args = parser.parse_args(argv)

    return(args)

def get_overlap(x, y):
    '''
    (x1, x2) and (y1, y2) define two intervals on the positive integer line
    where x2 >= x1 and y2 >= y1
    The y interval is shifted rightwards by a constant
    z is a normalized vector, i.e. z = (z1 + y1, z2 + y2)
    the output (a,b) is the overlap between z and y
    '''

    z = [xi + y[0] - 1 for xi in x]

    a = z[0]
    if a > y[1]:
        return (None, None)

    # Set end of exon match to max
    if z[1] >= y[1]:
        b = y[1]
    else:
        b = z[1]

    return((a,b))

def phaser(gff, intervals, delimiter=None):
    inter = Intervals(intervals, delimiter)
    for gene in reader.gff_reader(gff):
        for mrna in gene.mRNAs:
            bounds = inter.get_bounds(mrna.ident)
            if not bounds:
                continue
            for exon in mrna.exons:
                if not exon.CDS:
                    continue

                cds_length = exon.CDS.bounds[1] - exon.CDS.bounds[0] + 1

                a,b = get_overlap(bounds, exon.CDS.bounds)

                if a and b:
                    yield tuple([inter.get_ident(mrna.ident), mrna.ident, exon.num, gene.strand] + exon.bounds + [a, b])

                bounds = [x - cds_length for x in bounds]

                if b and b <= exon.bounds[1]:
                    break


if __name__ == '__main__':
    args = parse()

    for row in phaser(args.gff, args.intervals):
        print("%s %s %s %s %d %d %d %d" % row)
