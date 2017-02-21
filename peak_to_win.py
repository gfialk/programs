#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# This program receives peak files (".narrowPeak" format) and creates a respective
# ".win" file which contains the fraction of window covered by peaks.
# PROGRAM ARGUMENTS
# This program receives 2 arguments:
# 1) directory containing peak files
# 2) file containing chromosome sizes (tab delem)
# 3) output directory to where the program wrights the '*win' files.
# For each chromosome, the program devides the chromosome to windows of size WIN_SIZE and
# calculates the fraction of window covered by a peak. (for example, a peak in chr13 that
# spans from 12300948 to 12301534 will add 0.052 to the window 12300000, and 0.534
# to the window 12301000 - considering windows of size 10000)
#
# The program outputs data into file where each line is the format of (tab delemetered):
# chromosome_number   widnow_start_position   window_end_position   fraction_of_window_covered_by_peaks 

import sys
import os

WIN_SIZE = 10000
PEAK_DIR_ARG = 1
CHRM_SIZE_ARG = 2
WIN_DIR_ARG = 3 

def main(args):
    peak_dir = sys.argv[PEAK_DIR_ARG]
    chrom_sizes = sys.argv[CHRM_SIZE_ARG]
    win_dir = sys.argv[WIN_DIR_ARG]
    with open(chrom_sizes, 'r') as f:
        chr_sizes = f.readlines()

        # dict of chromosomes and there sizes
        chroms = {}
        for chr in chr_sizes:
            chrom, size = chr.strip('\n').split('\t')
            chroms[int(chrom)] = int(size)

    # iterate through all files in directory
    
    for peak_file in os.listdir(peak_dir):
        if peak_file.endswith(".narrowPeak"):
            print (peak_file)

            # list of 25 dict's representing chromosomes (including mitochondrion). 
            chr_peaks = [{} for i in range(len(chroms))]

            for chr in sorted(chroms.keys()):
                for win in range(0,chroms[chr], WIN_SIZE):
                    chr_peaks[chr-1][win//WIN_SIZE] = 0


            with open(peak_dir + peak_file, 'r') as peak_f:
                for line in peak_f:
                    fields = line.split('\t')
                    chrom = int(fields[0].strip('chr').replace('X','23')\
                                .replace('Y','24').replace('M','25')) - 1
                    start = int(fields[1])
                    end = int(fields[2])

                    s_cord = start//WIN_SIZE
                    e_cord = end//WIN_SIZE

                    # peak within one window
                    if s_cord == e_cord:
                        chr_peaks[chrom][s_cord] += e_cord-s_cord

                    # peak spans over more than one window
                    else:
                        # first window
                        chr_peaks[chrom][s_cord] += WIN_SIZE-(start%WIN_SIZE)

                        n_wins = e_cord-s_cord
                        # spaning windows
                        for win in range(n_wins-1):
                            chr_peaks[chrom][s_cord+win+1] += WIN_SIZE

                        # last_window
                        chr_peaks[chrom][s_cord+n_wins] += e_cord%WIN_SIZE

                peak_f.close()

            output_file = win_dir + peak_file.strip('narrowPeak') + 'win'
            matrix_f = open(output_file, 'w')
            for chr_inx,chrm in enumerate(chr_peaks):
                for win in sorted(chrm.keys()):
                    c_inx = str(chr_inx +1)
                    s_crd = str(win*WIN_SIZE)
                    e_crd = str(win*WIN_SIZE + WIN_SIZE)
                    read_frctn = str(chrm[win]/WIN_SIZE)
                    line_elm = c_inx, s_crd, e_crd, read_frctn
                    line = '\t'.join(line_elm) + '\n'
                    matrix_f.write(line)

            
            matrix_f.close()
    #peak_f.close()
    
    return 0

if __name__ == '__main__':
    import sys
    sys.exit(main(sys.argv))

