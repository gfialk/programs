__author__ = 'gf861913'
'''
plot histogram of fragment length form a BAM file.
The program receives a file containing the histogram of the read length of a bam file
and plots the histogram.

'''

import sys
import os
import numpy as np
import matplotlib.pyplot as plt

HIST_FILE_ARG = 1


def main(args):
    hist_file = sys.argv[HIST_FILE_ARG]
    with open (hist_file) as hf:
        hist = hf.readlines()
        hist = map(lambda x: x.strip(',\n').strip(), hist)
        hist_dic = {}
        for i in hist:
            n_rds = i.split()[0]
            len_rds = i.split()[1]
            hist_dic[n_rds] = len_rds

        plt.scatter(hist_dic.keys(), hist_dic.values())
        plt.show()






    return 0

if __name__ == '__main__':
    import sys
    sys.exit(main(sys.argv))

