#!/usr/bin/env python3
import lib.gffreader as gffreader
import unittest
import os

def readgff(gfflist):
    # If the input is a list of lists, join them
    if gfflist and not hasattr(gfflist[0], "strip"):
        gfflist = ["\t".join([str(e) for e in x]) for x in gfflist]

    out = []
    with open(os.devnull, 'w') as errout:
        for gene in gffreader.gff_reader(gfflist, errout=errout):
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
        self.good_output = [
            ['s1', 'a.1', '1', '1', '1000', '+', '1', 'a.1.exon.1', '100', '109', '.', '.', '.', '.'],
            ['s1', 'a.1', '1', '1', '1000', '+', '2', 'a.1.exon.2', '110', '200', '150', '200', '.', '0'],
            ['s1', 'a.1', '1', '1', '1000', '+', '3', 'a.1.exon.3', '300', '400', '300', '400', '0', '2'],
            ['s1', 'a.1', '1', '1', '1000', '+', '4', 'a.1.exon.4', '600', '900', '600', '603', '2', '.'],
            ['s1', 'a.1', '1', '1', '1000', '+', '5', 'a.1.exon.5', '910', '930', '.', '.', '.', '.']
        ]
        self.good = ['\t'.join(s) for s in self.good]
        self.good_output = [' '.join(s) for s in self.good_output]

        self.minus = [
            ['s1', '.', 'gene', '1',   '1000', '.', '-', '.', 'ID=a'],
            ['s1', '.', 'mRNA', '1',   '1000', '.', '-', '.', 'ID=a.1'],
            ['s1', '.', 'exon', '800', '900',  '.', '-', '.', 'ID=a.1.exon.1'],
            ['s1', '.', 'exon', '600', '700',  '.', '-', '.', 'ID=a.1.exon.2'],
            ['s1', '.', 'CDS',  '600', '650',  '.', '-', '.', 'ID=a.1.cds.1'],
            ['s1', '.', 'exon', '400', '500',  '.', '-', '.', 'ID=a.1.exon.3'],
            ['s1', '.', 'CDS',  '400', '500',  '.', '-', '.', 'ID=a.1.cds.2'],
            ['s1', '.', 'exon', '200', '300',  '.', '-', '.', 'ID=a.1.exon.4'],
            ['s1', '.', 'CDS',  '297', '300',  '.', '-', '.', 'ID=a.1.cds.2'],
            ['s1', '.', 'exon', '10',  '150',  '.', '-', '.', 'ID=a.1.exon.5']
        ]
        self.minus_output = [
            ['s1', 'a.1', '1', '1', '1000', '-', '1', 'a.1.exon.1', '800', '900',   '.',   '.', '.', '.'],
            ['s1', 'a.1', '1', '1', '1000', '-', '2', 'a.1.exon.2', '600', '700', '600', '650', '.', '0'],
            ['s1', 'a.1', '1', '1', '1000', '-', '3', 'a.1.exon.3', '400', '500', '400', '500', '0', '2'],
            ['s1', 'a.1', '1', '1', '1000', '-', '4', 'a.1.exon.4', '200', '300', '297', '300', '2', '.'],
            ['s1', 'a.1', '1', '1', '1000', '-', '5', 'a.1.exon.5', '10',  '150',   '.',   '.', '.', '.']
        ]
        self.minus = ['\t'.join(s) for s in self.minus]
        self.minus_output = [' '.join(s) for s in self.minus_output]

    def test_good(self):
        self.assertEqual(readgff(self.good), self.good_output)

    def test_minus(self):
        self.assertEqual(readgff(self.minus), self.minus_output)

class Test_a_is_downstream_of_b(unittest.TestCase):
    def test_true_plus(self):
        self.assertTrue(gffreader.a_is_downstream_of_b(a=2, b=1, strand="+"))
        self.assertTrue(gffreader.a_is_downstream_of_b(a=2, b=2, strand="+"))
    def test_false_plus(self):
        self.assertFalse(gffreader.a_is_downstream_of_b(a=1, b=2, strand="+"))
    def test_true_minus(self):
        self.assertTrue(gffreader.a_is_downstream_of_b(a=99, b=100, strand="-"))
        self.assertTrue(gffreader.a_is_downstream_of_b(a=100, b=100, strand="-"))
    def test_false_minus(self):
        self.assertFalse(gffreader.a_is_downstream_of_b(a=100, b=99, strand="-"))
class Test_a_is_upstream_of_b(unittest.TestCase):
    def test_true_plus(self):
        self.assertTrue(gffreader.a_is_upstream_of_b(a=1, b=2, strand="+"))
        self.assertTrue(gffreader.a_is_upstream_of_b(a=2, b=2, strand="+"))
    def test_false_plus(self):
        self.assertFalse(gffreader.a_is_upstream_of_b(a=2, b=1, strand="+"))
    def test_true_minus(self):
        self.assertTrue(gffreader.a_is_upstream_of_b(a=100, b=99, strand="-"))
        self.assertTrue(gffreader.a_is_upstream_of_b(a=100, b=100, strand="-"))
    def test_false_minus(self):
        self.assertFalse(gffreader.a_is_upstream_of_b(a=99, b=100, strand="-"))
class Test_a_is_within_b(unittest.TestCase):
    def test_within(self):
        self.assertTrue(gffreader.a_is_within_b((2, 5), (1, 6)))
        self.assertTrue(gffreader.a_is_within_b((5, 2), (1, 6)))
    def test_same_interval(self):
        self.assertTrue(gffreader.a_is_within_b((1, 5), (1, 5)))
    def test_shared_boundary(self):
        self.assertTrue(gffreader.a_is_within_b((2, 5), (2, 6)))
        self.assertTrue(gffreader.a_is_within_b((2, 5), (1, 5)))
    def test_false(self):
        self.assertFalse(gffreader.a_is_within_b((2, 5), (3, 5)))
        self.assertFalse(gffreader.a_is_within_b((2, 5), (2, 4)))
        self.assertFalse(gffreader.a_is_within_b((2, 5), (3, 4)))
    def test_bound_order(self):
        self.assertTrue(gffreader.a_is_within_b((5, 2), (2, 6)))
        self.assertTrue(gffreader.a_is_within_b((2, 5), (6, 2)))

class Test_a_overlaps_b(unittest.TestCase):
    def test_true(self):
        self.assertTrue(gffreader.a_overlaps_b((1,15), (10,20)))
        self.assertTrue(gffreader.a_overlaps_b((15,30), (10,20)))
        self.assertTrue(gffreader.a_overlaps_b((1,30), (10,20)))
        self.assertTrue(gffreader.a_overlaps_b((11,19), (10,20)))
    def test_equal_border(self):
        self.assertTrue(gffreader.a_overlaps_b((10,20), (10,20)))
        self.assertTrue(gffreader.a_overlaps_b((1,10), (10,20)))
        self.assertTrue(gffreader.a_overlaps_b((20,30), (10,20)))
    def test_false(self):
        self.assertFalse(gffreader.a_overlaps_b((1,9), (10,20)))
        self.assertFalse(gffreader.a_overlaps_b((21,30), (10,20)))

class Test_format_checkint(unittest.TestCase):
    def test_misformatted_no_gene(self):
        test = [
            ['s1', '.', 'mRNA', '1', '1000', '.', '+', '.', 'ID=a.1'],
            ['s1', '.', 'exon', '100', '109', '.', '+', '.', 'ID=a.1.exon.1']
        ]
        self.assertEqual(readgff(test), [])

    def test_misformatted_mRNA_out(self):
        test = [
            ['s1', '.', 'gene', '1', '1000', '.', '+', '.', 'ID=a'],
            ['s1', '.', 'mRNA', '1', '1001', '.', '+', '.', 'ID=a.1']
        ]
        self.assertEqual(readgff(test), [])

    def test_misformatted_exon_out(self):
        test = [
            ['s1', '.', 'gene', '1', '1000', '.', '+', '.', 'ID=a'],
            ['s1', '.', 'mRNA', '1', '1000', '.', '+', '.', 'ID=a.1'],
            ['s1', '.', 'exon', '100', '1001', '.', '+', '.', 'ID=a.1.exon.1']
        ]
        self.assertEqual(readgff(test), [])

    def test_misformatted_cds_out(self):
        test = [
            ['s1', '.', 'gene', '1', '1000', '.', '+', '.', 'ID=a'],
            ['s1', '.', 'mRNA', '1', '1000', '.', '+', '.', 'ID=a.1'],
            ['s1', '.', 'exon', '100', '1001', '.', '+', '.', 'ID=a.1.exon.1'],
            ['s1', '.', 'CDS',  '30', '131',  '.', '+', '.', 'ID=a.1.cds.1']
        ]
        self.assertEqual(readgff(test), [])

    def test_misformatted_wrong_strand(self):
        test = [
            ['s1', '.', 'gene', '1', '1000', '.', '+', '.', 'ID=a'],
            ['s1', '.', 'mRNA', '1', '1000', '.', '-', '.', 'ID=a.1']
        ]
        self.assertEqual(readgff(test), [])

    def test_misformatted_no_strand(self):
        test = [
            ['s1', '.', 'gene', '1', '1000', '.', '.', '.', 'ID=a']
        ]
        self.assertEqual(readgff(test), [])

    def test_misformatted_no_mRNA(self):
        test = [
            ['s1', '.', 'gene', '1', '1000', '.', '+', '.', 'ID=a'],
            ['s1', '.', 'exon', '100', '1001', '.', '+', '.', 'ID=a.1.exon.1'],
        ]
        self.assertEqual(readgff(test), [])

    def test_misformatted_CDS_without_exon(self):
        test = [
            ['s1', '.', 'gene', '1', '1000', '.', '+', '.', 'ID=a'],
            ['s1', '.', 'mRNA', '1', '1000', '.', '+', '.', 'ID=a.1'],
            ['s1', '.', 'CDS',  '30', '131',  '.', '-', '.', 'ID=a.1.cds.1']
        ]
        self.assertEqual(readgff(test), [])

    def test_misformatted_exons_misordered_plus(self):
        test = [
            ['s1', '.', 'gene', '1', '1000', '.', '+', '.', 'ID=a'],
            ['s1', '.', 'mRNA', '1', '1000', '.', '+', '.', 'ID=a.1'],
            ['s1', '.', 'exon', '110', '200', '.', '+', '.', 'ID=a.1.exon.2'],
            ['s1', '.', 'exon', '100', '109', '.', '+', '.', 'ID=a.1.exon.1']
        ]
        self.assertEqual(readgff(test), [])

    def test_misformatted_overlapping_exons(self):
        test = [
            ['s1', '.', 'gene', '1', '1000', '.', '+', '.', 'ID=a'],
            ['s1', '.', 'mRNA', '1', '1000', '.', '+', '.', 'ID=a.1'],
            ['s1', '.', 'exon', '100', '200', '.', '+', '.', 'ID=a.1.exon.1'],
            ['s1', '.', 'exon', '150', '300', '.', '+', '.', 'ID=a.1.exon.2']
        ]
        self.assertEqual(readgff(test), [])

    def test_no_gene_ident(self):
        test = [
            ['s1', '.', 'gene', '1', '1000', '.', '+', '.', 'a']
        ]
        self.assertEqual(readgff(test), [])

    def test_no_mRNA_ident(self):
        test = [
            ['s1', '.', 'gene', '1', '1000', '.', '+', '.', 'ID=a'],
            ['s1', '.', 'mRNA', '1', '1000', '.', '+', '.', 'a.1']
        ]
        self.assertEqual(readgff(test), [])

    def test_no_exon_ident(self):
        test = [
            ['s1', '.', 'gene', '1', '1000', '.', '+', '.', 'ID=a'],
            ['s1', '.', 'mRNA', '1', '1000', '.', '+', '.', 'ID=a.1'],
            ['s1', '.', 'exon', '100', '200', '.', '+', '.', 'a.1.exon.1']
        ]
        self.assertEqual(readgff(test), [])

    def test_no_cds_ident(self):
        test = [
            ['s1', '.', 'gene', '1', '1000', '.', '+', '.', 'ID=a'],
            ['s1', '.', 'mRNA', '1', '1000', '.', '+', '.', 'ID=a.1'],
            ['s1', '.', 'exon', '150', '300', '.', '+', '.', 'ID=a.1.exon.2'],
            ['s1', '.', 'CDS',  '160', '162',  '.', '+', '.', 'a.1.cds.1']
        ]
        self.assertEqual(readgff(test), [])


if __name__ == '__main__':
    unittest.main()
