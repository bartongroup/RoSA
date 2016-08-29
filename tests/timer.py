#!/usr/bin/env python

import time
from summary_routines.parser import Parser
from parsing_routines.gff_gtf_tools import annotation

__author__ = 'kmourao'


def main():
    #data = "D:\SummaryData\TAIR10_GFF3_genes.gff"
    data = "C:\Users\kira\Documents\Dundee\TAIR10_GFF3_genes.gff"
    t0 = time.clock()
    p = Parser(data)
    t1 = time.clock()
    panda_time = t1-t0

    t0 = time.clock()
    a = annotation(data, filetype='gff3')
    t1 = time.clock()
    dict_time = t1-t0

    print("Panda time: {} Annotation time: {}".format(panda_time, dict_time))

if __name__ == "__main__":
    main()