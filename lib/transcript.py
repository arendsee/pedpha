class Transcript:
    def __init__(self, transid, transnum, bounds, strand):
        self.transid = transid
        self.transnum = transnum
        self.bounds = bounds
        self.strand = strand
        self.exons = []
        self._trans_length = 0
        self._cds_length = 0

    def add_exon(self, exon):
        exon.num = exon.num if exon.num else len(self.exons) + 1
        self._trans_length += exon.bounds.length()
        self.exons.append(exon)

    def add_cds(self, bounds):
        self.exons[-1].add_cds(bounds, self._cds_length, self.strand)
        self._cds_length += bounds.length()

class Exon:
    def __init__(self, ident, bounds, num=None):
        self.num = num
        self.ident = ident
        self.bounds = bounds
        self.CDS = None
        self.phase = (".", ".")

    def add_cds(ebounds, cbounds, cds_length, strand):
        self.CDS = cbounds
        self.phase = Exon.phase(ebounds, cbounds, cds_length, strand)

    @classmethod
    def phase(cls, ebounds, cbounds, offset, strand):
        estart, estop = ebounds if strand == "+" else reversed(ebounds)
        cstart, cstop = cbounds if strand == "+" else reversed(cbounds)
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
