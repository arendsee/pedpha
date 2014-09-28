#!/usr/bin/env python3
import gffreader
import unittest

def readgff(gfflist):
    out = []
    for gene in gffreader.gff_reader(gfflist):
        for line in gene.tostr():
            out.append(line)
    return(out)

class Test_gffreader(unittest.TestCase):
    def setUp(self):
        self.good = [
            ['s1', '.', 'gene', '1', '1000', '.', '+', '.', 'ID=a'],
            ['s1', '.', 'mRNA', '1', '1000', '.', '+', '.', 'ID=a.1'],
            ['s1', '.', 'exon', '100', '109', '.', '+', '.', 'ID=a.1.exon.1'],
            ['s1', '.', 'exon', '110', '200', '.', '+', '.', 'ID=a.1.exon.2'],
            ['s1', '.', 'CDS', '150', '200', '.', '+', '.', 'ID=a.1.cds.1'],
            ['s1', '.', 'exon', '300', '400', '.', '+', '.', 'ID=a.1.exon.3'],
            ['s1', '.', 'CDS', '300', '400', '.', '+', '.', 'ID=a.1.cds.2'],
            ['s1', '.', 'exon', '600', '900', '.', '+', '.', 'ID=a.1.exon.4'],
            ['s1', '.', 'CDS', '600', '603', '.', '+', '.', 'ID=a.1.cds.2'],
            ['s1', '.', 'exon', '910', '930', '.', '+', '.', 'ID=a.1.exon.5']
        ]
        self.good = ['\t'.join(s) for s in self.good]
        self.good_output = [
            ['s1', 'a.1', '1', '1', '1000', '+', '1', 'a.1.exon.1', '100', '109', '.', '.', '.', '.'],
            ['s1', 'a.1', '1', '1', '1000', '+', '2', 'a.1.exon.2', '110', '200', '150', '200', '0', '0'],
            ['s1', 'a.1', '1', '1', '1000', '+', '3', 'a.1.exon.3', '300', '400', '300', '400', '0', '2'],
            ['s1', 'a.1', '1', '1', '1000', '+', '4', 'a.1.exon.4', '600', '900', '600', '603', '2', '0'],
            ['s1', 'a.1', '1', '1', '1000', '+', '5', 'a.1.exon.5', '910', '930', '.', '.', '.', '.']
        ]
        self.good_output = [' '.join(s) for s in self.good_output]

    def test_good(self):
        self.assertEqual(readgff(self.good), self.good_output)


if __name__ == '__main__':
    unittest.main()
