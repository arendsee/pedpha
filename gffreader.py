#!/usr/bin/env python3

'''
Assumptions about gff format:
1) Has types 'gene', 'mRNA', 'exon', and 'CDS'
2) These are all ordered in Phytozome style
3) Exons are ordered according to gene beginning not absolute
'''

def gff_reader(gfffile):
    for line in gfffile:
        d = line2gffdict(line)
        if not d:
            continue
        if d['type'] == 'gene':
            try:
                yield g
            except NameError:
                pass
            g = Gene(d['seqid'], d['start'], d['stop'], d['strand'])
        elif d['type'] == 'mRNA':
            mrna = mRNA()
            try:
                g.add_mRNA(mrna)
            except NameError:
                sys.exit('No gene entry for this mRNA')
        elif d['type'] == 'exon':
            exon = Exon(d['start'], d['stop'])
            g.mRNAs[-1].add_exon(exon)
        elif d['type'] == 'CDS':
            g.mRNAs[-1].exons[-1].CDS = CDS(d['start'], d['stop'])
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
    def __init__(self, seqid, start, stop, strand):
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
    def __init__(self, tid=None):
        self.tid = tid
        self.exons = []

    def tostr(self, gene):
        template = '{} {} {} {} {} {{}}'.format(gene.seqid,
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


class Exon:
    def __init__(self, start, stop, num=None):
        self.num = num
        self.start = start
        self.stop = stop
        self.CDS = None

    def tostr(self):
        return("{} {} {}".format(self.num, self.start, self.stop))

class CDS:
    def __init__(self, start, stop):
        self.start = start
        self.stop = stop


if __name__ == '__main__':
    import sys
    for gene in gff_reader(sys.stdin):
        for line in gene.tostr():
            print(line)
