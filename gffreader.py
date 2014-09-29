#!/usr/bin/env python3

'''
Assumptions about gff format:
1) Has types 'gene', 'mRNA', 'exon', and 'CDS'
2) These are all ordered in Phytozome style
3) Exons are ordered according to gene beginning not absolute
'''

import re

def a_is_downstream_of_b(a, b, strand):
    if strand == "+":
        return(a >= b)
    else:
        return(a <= b)

def a_is_upstream_of_b(a, b, strand):
    if strand == "+":
        return(a <= b)
    else:
        return(a >= b)

def a_is_within_b(a, b):
    is_lowerbounded = bool(min(a) >= min(b) and min(a) <= max(b))
    is_upperbounded = bool(max(a) <= max(b))
    return(is_lowerbounded and is_upperbounded)



def parse_desc(desc):
    return re.sub('^ID=([^;]+).*', '\\1', desc).strip()

def gff_reader(gfffile):
    g = None
    for line in gfffile:
        d = line2gffdict(line)
        if not d:
            continue
        ident = parse_desc(d['desc'])

        # Assert the new element is properly placed within the gene
        if g and d['type'] in ('mRNA', 'CDS', 'exon'):
            assert d['strand'] == g.strand, \
                    "All gene elements must be on same strand:\n%s" % line
            assert a_is_within_b(d['bounds'], g.bounds), \
                   "All gene elements must be within gene boundaries gene:\n%s" % line

        if d['type'] == 'gene':
            if g:
                yield g
            g = Gene(ident, d['seqid'], d['bounds'], d['strand'])
        elif d['type'] == 'mRNA':
            mrna = mRNA(ident, d['bounds'])
            try:
                g.add_mRNA(mrna)
            except NameError:
                sys.exit('No gene entry for this mRNA')
        elif d['type'] == 'exon':
            exon = Exon(ident, d['bounds'])
            g.mRNAs[-1].add_exon(exon)
            try:
                assert a_is_downstream_of_b(exon.bounds[0], g.mRNAs[-1].exons[-2].bounds[1], g.strand), \
                        "Exons must be ordered:\n%s" % line
            except IndexError:
                pass
        elif d['type'] == 'CDS':
            g.mRNAs[-1].exons[-1].CDS = CDS(ident, d['bounds'])
            assert a_is_within_b(d['bounds'], g.mRNAs[-1].exons[-1].bounds), \
                   "CDS must be within exon"
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
        out['bounds'] = (int(out['start']), int(out['stop']))
    else:
        out = None
    return out

class Gene:
    def __init__(self, ident, seqid, bounds, strand):
        self.ident = ident
        self.seqid = seqid
        self.bounds = bounds
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
    def __init__(self, ident, bounds, tid=None):
        self.ident = ident
        self.bounds = bounds
        self.tid = tid
        self.exons = []

    def tostr(self, gene):
        self.calculate_phases()
        template = '{} {} {} {} {} {} {{}}'.format(gene.seqid,
                                                self.ident,
                                                self.tid,
                                                gene.bounds[0],
                                                gene.bounds[1],
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
            if exon.CDS:
                estart, estop = exon.bounds
                cstart, cstop = exon.CDS.bounds
                new_length = cds_length + abs(cstop - cstart) + 1
                if cstart == estart:
                    p5 = cds_length % 3
                else:
                    p5 = 0

                if cstop == estop:
                    p3 = new_length % 3
                else:
                    p3 = 0
                exon.phase = (p5, p3)
                cds_length = new_length

class Exon:
    def __init__(self, ident, bounds, num=None):
        self.num = num
        self.ident = ident
        self.bounds = bounds
        self.CDS = None
        self.phase = (".", ".")

    def tostr(self):
        if self.CDS:
            cstart, cstop = self.CDS.bounds
        else:
            cstart, cstop = ".", "."
        out = " ".join(str(s) for s in (
            self.num,
            self.ident,
            self.bounds[0],
            self.bounds[1],
            cstart,
            cstop,
            self.phase[0], self.phase[1]))
        return out

class CDS:
    def __init__(self, ident, bounds):
        self.ident = ident
        self.bounds = bounds


if __name__ == '__main__':
    import sys
    for gene in gff_reader(sys.stdin):
        for line in gene.tostr():
            print(line)
