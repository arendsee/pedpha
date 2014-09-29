#!/usr/bin/env python3

'''
This is currently only tested for phytozome10 formatted *.gene_exons.gff3 files
'''

import lib.gffreader as reader
import sys

for gene in reader.gff_reader(sys.stdin):
    for line in gene.tostr():
        print(line)
