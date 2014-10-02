class Interval:
    def __init__(self, a, b, strand=None):
        a, b = int(a), int(b)
        self.start = min(a, b)
        self.stop = max(a, b)
        self.strand = strand
