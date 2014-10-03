class Interval:
    def __init__(self, a, b, strand=None):
        a, b = int(a), int(b)
        self.a = min(a, b)
        self.b = max(a, b)
        self.strand = strand

    def to_dna_interval(self):
        '''
        x is an interval (x1, x2) on a protein
        Where x1 and x2 are integers greater than 0 and x2 >= x1
        '''
        a2, b2 = [(i - 1)*3 + 1 for i in (self.a, self.b)]
        b2 += 2
        newInterval = Interval(a=a2, b=b2, strand=self.strand)
        return(newInterval)

    def point_is_downstream(self, x):
        if self.strand == "+":
            return(x >= self.b)
        else:
            return(x <= self.a)

    def point_is_upstream(self, x):
        if strand == "+":
            return(x <= self.a)
        else:
            return(x >= self.b)

    def interval_is_within(self, x):
        return(x.a >= self.a and x.b <= self.b)

    def interval_overlaps(self, x):
        return(self.point_is_within(x.a) or self.point_is_within(x.b))

    def get_subinterval(self, x):
        if self.strand == "-":
            a, b = [self.b - xi + 1 for xi in x]


            if a < self.a:
                return (None, None)

            if b <= self.a:
                b = self.a

            return(Interval(b, a, self.strand))

        else:
            a, b = [xi + self.a - 1 for xi in x]

            if a > self.b:
                return (None, None)

            if b >= self.b:
                b = self.b

            return(Interval(b, a, self.strand))
