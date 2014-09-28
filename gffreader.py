#!/usr/bin/env python3

'''
Assumptions about gff format:
1) Has types 'gene', 'mRNA', 'exon', and 'CDS'
2) These are all ordered in Phytozome style
3) Exons are ordered according to gene beginning not absolute
'''

import re

def parse_desc(desc):
    return re.sub('^ID=([^;]+).*', '\\1', desc).strip()

# TODO add format checking
def gff_reader(gfffile):
    for line in gfffile:
        d = line2gffdict(line)
        if not d:
            continue
        ident = parse_desc(d['desc'])
        if d['type'] == 'gene':
            try:
                yield g
            except NameError:
                pass
            g = Gene(ident, d['seqid'], d['start'], d['stop'], d['strand'])
        elif d['type'] == 'mRNA':
            mrna = mRNA(ident)
            try:
                g.add_mRNA(mrna)
            except NameError:
                sys.exit('No gene entry for this mRNA')
        elif d['type'] == 'exon':
            exon = Exon(ident, d['start'], d['stop'])
            g.mRNAs[-1].add_exon(exon)
        elif d['type'] == 'CDS':
            g.mRNAs[-1].exons[-1].CDS = CDS(ident, d['start'], d['stop'])
    try:
        yield g
    except NameError:
        pass

def line2gffdict(line):
    row = line.split('\t')
    if len(row) == 9:
        out = dict(zip(
            ['seqid', 'source', 'type', 'start', 'stop', 'score', 'strand', 'phase', 'desc'],
            row
        ))
    else:
        out = None
    return out

class Gene:
    def __init__(self, ident, seqid, start, stop, strand):
        self.ident = ident
        self.seqid = seqid
        self.start = start
        self.stop = stop
        self.strand = strand
        self.mRNAs = []

    def tostr(self):
        lines = []
        for mrna in self.mRNAs:
            lines += mrna.tostr(self)
        return(lines)

    def add_mRNA(self, mrna):
        mrna.tid = mrna.tid if mrna.tid else len(self.mRNAs) + 1
        self.mRNAs.append(mrna)

class mRNA:
    def __init__(self, ident, tid=None):
        self.ident = ident
        self.tid = tid
        self.exons = []

    def tostr(self, gene):
        self.calculate_phases()
        template = '{} {} {} {} {} {} {{}}'.format(gene.seqid,
                                                self.ident,
                                                self.tid,
                                                gene.start,
                                                gene.stop,
                                                gene.strand)
        lines = []
        for exon in self.exons:
            lines.append(template.format(exon.tostr()))
        return(lines)

    def add_exon(self, exon):
        exon.num = exon.num if exon.num else len(self.exons) + 1
        self.exons.append(exon)

    def calculate_phases(self):
        cds_length = 0
        for exon in self.exons:
            if exon.CDS.start != "." or exon.CDS.stop != ".":
                new_length = cds_length + abs(int(exon.CDS.stop) - int(exon.CDS.start)) + 1
                if exon.CDS.start == exon.start:
                    p5 = cds_length % 3
                else:
                    p5 = 0

                if exon.CDS.stop == exon.stop:
                    p3 = new_length % 3
                else:
                    p3 = 0
                exon.phase = (p5, p3)
                cds_length = new_length

class Exon:
    def __init__(self, ident, start, stop, num=None):
        self.num = num
        self.ident = ident
        self.start = start
        self.stop = stop
        self.CDS = CDS(None)
        self.phase = (".", ".")

    def tostr(self):
        out = " ".join(str(s) for s in (
            self.num,
            self.ident,
            self.start,
            self.stop,
            self.CDS.start,
            self.CDS.stop,
            self.phase[0], self.phase[1]))
        return out

class CDS:
    def __init__(self, ident, start=".", stop="."):
        self.ident = ident
        self.start = start
        self.stop = stop


if __name__ == '__main__':
    import sys
    for gene in gff_reader(sys.stdin):
        for line in gene.tostr():
            print(line)
