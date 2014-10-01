#!/usr/bin/env python3

'''
Assumes gff3 files follow the phytozome/metazom convetions. Any deviance from
this format should instantly kill the program. See README.
'''

import re
import sys

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

def a_overlaps_b(a, b):
    return(a_is_within_b([a[0]], b) or a_is_within_b([a[1]], b) or a_is_within_b(b, a))

def parse_desc(desc):
    try:
        return(re.match('^ID=([^;]+).*', desc).group(1))
    except AttributeError:
        return(None)

def line2gffdict(line):
    row = line.split('\t')
    if len(row) == 9:
        out = dict(zip(
            ['seqid', 'source', 'type', 'start', 'stop', 'score', 'strand', 'phase', 'desc'],
            row
        ))
        out['bounds'] = sorted([int(out['start']), int(out['stop'])])
    else:
        out = None
    return out

def gff_reader(gfffile, errout=sys.stderr):
    g = None
    valid = True
    fc = FormatChecker(errout)
    for line in gfffile:
        line = line.strip()
        d = line2gffdict(line)
        if not d:
            continue
        ident = parse_desc(d['desc'])

        # Skip the gene if any element is invalid (FormatChecker handles
        # warning messages)
        if valid and not fc.check_gene_element(g, d, ident, line):
            g = None
            valid = False

        if d['type'] == 'gene':
            # If the gene object is well formed, yield
            # Otherwise contine, writing warnings to STDERR
            if fc.check_gene(g):
                yield g
            g = Gene(ident, d['seqid'], d['bounds'], d['strand'])
            valid = True

        elif d['type'] == 'mRNA' and valid:
            mrna = mRNA(ident, d['bounds'], g.strand)
            g.add_mRNA(mrna)

        elif d['type'] == 'exon' and valid:
            exon = Exon(ident, d['bounds'])
            g.mRNAs[-1].add_exon(exon)

        elif d['type'] == 'CDS' and valid:
            g.mRNAs[-1].exons[-1].CDS = CDS(ident, d['bounds'])

    if fc.check_gene(g):
        yield g

def phase(ebounds, cbounds, offset, isplus):
    estart, estop = ebounds if isplus else reversed(ebounds)
    cstart, cstop = cbounds if isplus else reversed(cbounds)
    new_length = offset + abs(cstop - cstart) + 1
    if cstart == estart:
        p5 = offset % 3
    else:
        p5 = "."

    if cstop == estop:
        p3 = new_length % 3
    else:
        p3 = "."
    return(p5, p3)

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
    def __init__(self, ident, bounds, strand, tid=None):
        self.ident = ident
        self.bounds = bounds
        self.strand = strand
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
        offset = 0
        isplus = bool(self.strand == "+")
        for exon in self.exons:
            if exon.CDS:
                exon.phase = phase(exon.bounds, exon.CDS.bounds, offset, isplus)
                offset += abs(exon.CDS.bounds[1] - exon.CDS.bounds[0]) + 1

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

class FormatChecker:
    def __init__(self, errout=sys.stderr):
        self.errout = errout

    def _format_warning(self, msg, line=None):
        warnmsg = "GFF format error, skipping offending gene - %s"
        if line:
            warnmsg += "\n%s" % line
        print(warnmsg % msg, file=self.errout)


    def check_gene(self, g):
        valid = True
        if not g:
            return False
        for mrna in g.mRNAs:
            # ASSERT mRNA is within gene
            if not a_is_within_b(mrna.bounds, g.bounds):
                msg = "mRNA '%s' at (%d, %d) must be within gene bounds"
                self._format_warning(msg % tuple([mrna.ident] + mrna.bounds))
                valid = False

            priorexon = None
            for exon in mrna.exons:
                # ASSERT exon is within mRNA
                if not a_is_within_b(exon.bounds, mrna.bounds):
                    msg = "Exon '%s' at (%d, %d) must be within parent mRNA bounds (%d, %d)"
                    self._format_warning(msg % tuple([exon.ident] + exon.bounds + mrna.bounds))
                    valid = False
                    # ASSERT the exons are ordered (1 .. n), whether start
                    # positions are increasing or decreasing depends on the strand

                if priorexon:
                    if a_overlaps_b(exon.bounds, priorexon.bounds):
                        msg = "Exons '%s' and '%s' overlap"
                        self._format_warning(msg % (exon.ident, priorexon.ident))
                        valid = False
                    if not a_is_downstream_of_b(exon.bounds[0], priorexon.bounds[1], g.strand):
                        msg = "Exons '%s' and '%s' are out of order"
                        self._format_warning(msg % (exon.ident, priorexon.ident))
                        valid = False

                if exon.CDS:
                    if not a_is_within_b(exon.CDS.bounds, exon.bounds):
                        msg = "CDS %s at (%d, %d) must be within exon %s at (%d, %d)"
                        self._format_warning(msg % tuple([exon.CDS.ident] + exon.CDS.bounds + [exon.ident] + exon.bounds))
                        valid = False

                priorexon = exon

        return(valid)

    def check_gene_element(self, g, d, ident, line):
        valid = True
        # ASSERT the new element is properly placed within the gene
        if g and d['type'] in ('mRNA', 'CDS', 'exon'):
            if not ident:
                msg = "%s lacks identifier (ID=([^;]+))" % d['type']
                self._format_warning(msg, line)
                valid = False

            # ASSERT children are on same strand as parent gene
            if not d['strand'] == g.strand:
                self._format_warning("All gene elements must be on same strand", line)
                valid = False
            # ASSERT the seqid, i.e. the sequence to which the gff maps the
            # entry, is same in parent and child
            if not d['seqid'] == g.seqid:
                self._format_warning("mRNA is not from same sequence as expected parent", line)
                valid = False
            # ASSERT child is within the parent interval
            if not a_is_within_b(d['bounds'], g.bounds):
                self._format_warning("All gene elements must be within gene boundaries gene", line)
                valid = False

        if d['type'] in ('mRNA', 'CDS', 'exon') and not g:
            msg = "%s (%s at (%d, %d)) found outside of gene context"
            self._format_warning(msg % tuple([d['type'], d['seqid']] + d['bounds']))
            valid = False

        return(valid)

