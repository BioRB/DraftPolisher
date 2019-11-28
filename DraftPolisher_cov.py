#!/usr/bin/env python
# DraftPolisher: a tool for the polishing of draft sequences
#     Copyright (C) 2019 Rosario Nicola Brancaccio
#
#     This program is free software: you can redistribute it and/or modify
#     it under the terms of the GNU Affero General Public License as
#     published by the Free Software Foundation, either version 3 of the
#     License, or (at your option) any later version.
#
#     This program is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU Affero General Public License for more details.
#
#     You should have received a copy of the GNU Affero General Public License
#     along with this program.  If not, see <https://www.gnu.org/licenses/>
#
# first you need to add execute permissions to the tool file as follows:
# chmod u+x DRAFT_POLISHER.py
# perform muscle allignment as follows:
# input file for Muscle format:
# >QRY
# draft sequence
# >SBJ
# reference sequence
# use Muscle with "-clw" flag and put query sequence at first as follows:
# muscle -in input.fa -clw > outfile.fa
# the length of the sequences IDs in input.fa, has to be exactly 3 characters
# than use DRAFT_POLISHER as follows:
# python DRAFT_POLISHER.py (SPAdes)contigs.FILE muscle_out.FILE
# max gap allowed per each position, between draft and reference sequences is:
# 26(*2)bp in a max fragment length of 34(*2)bp
# meaning that, considering a "i" base position (in the draft or reference sequence):
# and considering a window of (34+i+34)bp = 69bp
# per each side of that window, the maximum number of gaps allowed is 26bp
import os
import re
from subprocess import Popen
import subprocess
from itertools import islice
import pandas as pd
import numpy as np
import argparse
from shlex import split

print("DRAFT_POLISHER v1.0.cov   by Rosario Nicola Brancaccio")
print("Start")
# MUSCLE alignment
parser = argparse.ArgumentParser(description='polish draft circular genomes')
parser.add_argument("--input1", "-q", help="draft sequence in Fasta format, use 'QRY' as ID", type=str)
parser.add_argument("--input2", "-s", help="reference genome in Fasta format, use 'SBJ' as ID", type=str)
parser.add_argument("--input3", "-f", help="sequences database to use as reference for the polishing of the draft sequence", type=str)
args = parser.parse_args()
os.system("cat {0} {1} > inseq.fa".format(args.input1, args.input2))
os.system("muscle -in inseq.fa -clw > muscle_out.fa")
os.system("sed -e '1,3d' < muscle_out.fa > m1.txt")
with open("m1.txt") as f:
    for line in islice(f, 0, None, 4):
        print(line[16:76].strip('\n'), file=open("output1.txt", "a+"))
# remove spaces and newlines
os.system("tr -d '\n\r' < output1.txt | sed 's/-/_/g' > ot1.txt")
# add newline in the end of file
command = '''
sed -i -e '$a\\ ' ot1.txt
'''
process = Popen(command, shell=True, stdout=subprocess.PIPE)
result = process.communicate()
with open('m1.txt') as f:
    for line in islice(f, 1, None, 4):
        print(line[16:76].strip('\n'), file=open("output2.txt", "a+"))
os.system("tr -d '\n\r' < output2.txt | sed 's/-/_/g' > ot2.txt")
command = '''
sed -i -e '$a\\ ' ot2.txt
'''
process = Popen(command, shell=True, stdout=subprocess.PIPE)
result = process.communicate()
with open('m1.txt') as f:
    for line in islice(f, 2, None, 4):
        print(line[16:76].strip('\n'), file=open("output3.txt", "a+"))
# remove newlines and substitute the spaces with m to indicate the mismatches
os.system("tr -d '\n\r' < output3.txt > ot3.txt| sed 's/ /m/g' ot3.txt > ot3p.txt ")
os.system("cat ot1.txt ot2.txt ot3p.txt > ff.txt")
# find position of mismatches and retrieve 18 bp
with open('ff.txt') as pr:
    lines = pr.readlines()
    hh = lines[4].rstrip("\n\r")
    pos = [m.start() for m in re.finditer('m', hh)]
    # processing of line 0 that correspond to the Query sequence, the draft genome to correct (query)
    ja = lines[0].rstrip("\n\r")
    ln = len(hh)
    for i in pos:
        if i == ln - 1:
            i1 = 0
            i2 = 1
            i3 = 2
            i4 = 3
            i5 = 4
            i6 = 5
            i7 = 6
            i8 = 7
            i9 = 8
            i10 = 9
            i11 = 10
            i12 = 11
            i13 = 12
            i14 = 13
            i15 = 14
            i16 = 15
            i17 = 16
            i18 = 17
            i19 = 18
            i20 = 19
            i21 = 20
            i22 = 21
            i23 = 22
            i24 = 23
            i25 = 24
            i26 = 25
            i27 = 26
            i28 = 27
            i29 = 28
            i30 = 29
            i31 = 30
            i32 = 31
            i33 = 32
            i34 = 33
        elif i == ln - 2:
            i1 = ln - 1
            i2 = 0
            i3 = 1
            i4 = 2
            i5 = 3
            i6 = 4
            i7 = 5
            i8 = 6
            i9 = 7
            i10 = 8
            i11 = 9
            i12 = 10
            i13 = 11
            i14 = 12
            i15 = 13
            i16 = 14
            i17 = 15
            i18 = 16
            i19 = 17
            i20 = 18
            i21 = 19
            i22 = 20
            i23 = 21
            i24 = 22
            i25 = 23
            i26 = 24
            i27 = 25
            i28 = 26
            i29 = 27
            i30 = 28
            i31 = 29
            i32 = 30
            i33 = 31
            i34 = 32
        elif i == ln - 3:
            i1 = ln - 2
            i2 = ln - 1
            i3 = 0
            i4 = 1
            i5 = 2
            i6 = 3
            i7 = 4
            i8 = 5
            i9 = 6
            i10 = 7
            i11 = 8
            i12 = 9
            i13 = 10
            i14 = 11
            i15 = 12
            i16 = 13
            i17 = 14
            i18 = 15
            i19 = 16
            i20 = 17
            i21 = 18
            i22 = 19
            i23 = 20
            i24 = 21
            i25 = 22
            i26 = 23
            i27 = 24
            i28 = 25
            i29 = 26
            i30 = 27
            i31 = 28
            i32 = 29
            i33 = 30
            i34 = 31
        elif i == ln - 4:
            i1 = ln - 3
            i2 = ln - 2
            i3 = ln - 1
            i4 = 0
            i5 = 1
            i6 = 2
            i7 = 3
            i8 = 4
            i9 = 5
            i10 = 6
            i11 = 7
            i12 = 8
            i13 = 9
            i14 = 10
            i15 = 11
            i16 = 12
            i17 = 13
            i18 = 14
            i19 = 15
            i20 = 16
            i21 = 17
            i22 = 18
            i23 = 19
            i24 = 20
            i25 = 21
            i26 = 22
            i27 = 23
            i28 = 24
            i29 = 25
            i30 = 26
            i31 = 27
            i32 = 28
            i33 = 29
            i34 = 30
        elif i == ln - 5:
            i1 = ln - 4
            i2 = ln - 3
            i3 = ln - 2
            i4 = ln - 1
            i5 = 0
            i6 = 1
            i7 = 2
            i8 = 3
            i9 = 4
            i10 = 5
            i11 = 6
            i12 = 7
            i13 = 8
            i14 = 9
            i15 = 10
            i16 = 11
            i17 = 12
            i18 = 13
            i19 = 14
            i20 = 15
            i21 = 16
            i22 = 17
            i23 = 18
            i24 = 19
            i25 = 20
            i26 = 21
            i27 = 22
            i28 = 23
            i29 = 24
            i30 = 25
            i31 = 26
            i32 = 27
            i33 = 28
            i34 = 29
        elif i == ln - 6:
            i1 = ln - 5
            i2 = ln - 4
            i3 = ln - 3
            i4 = ln - 2
            i5 = ln - 1
            i6 = 0
            i7 = 1
            i8 = 2
            i9 = 3
            i10 = 4
            i11 = 5
            i12 = 6
            i13 = 7
            i14 = 8
            i15 = 9
            i16 = 10
            i17 = 11
            i18 = 12
            i19 = 13
            i20 = 14
            i21 = 15
            i22 = 16
            i23 = 17
            i24 = 18
            i25 = 19
            i26 = 20
            i27 = 21
            i28 = 22
            i29 = 23
            i30 = 24
            i31 = 25
            i32 = 26
            i33 = 27
            i34 = 28
        elif i == ln - 7:
            i1 = ln - 6
            i2 = ln - 5
            i3 = ln - 4
            i4 = ln - 3
            i5 = ln - 2
            i6 = ln - 1
            i7 = 0
            i8 = 1
            i9 = 2
            i10 = 3
            i11 = 4
            i12 = 5
            i13 = 6
            i14 = 7
            i15 = 8
            i16 = 9
            i17 = 10
            i18 = 11
            i19 = 12
            i20 = 13
            i21 = 14
            i22 = 15
            i23 = 16
            i24 = 17
            i25 = 18
            i26 = 19
            i27 = 20
            i28 = 21
            i29 = 22
            i30 = 23
            i31 = 24
            i32 = 25
            i33 = 26
            i34 = 27
        elif i == ln - 8:
            i1 = ln - 7
            i2 = ln - 6
            i3 = ln - 5
            i4 = ln - 4
            i5 = ln - 3
            i6 = ln - 2
            i7 = ln - 1
            i8 = 0
            i9 = 1
            i10 = 2
            i11 = 3
            i12 = 4
            i13 = 5
            i14 = 6
            i15 = 7
            i16 = 8
            i17 = 9
            i18 = 10
            i19 = 11
            i20 = 12
            i21 = 13
            i22 = 14
            i23 = 15
            i24 = 16
            i25 = 17
            i26 = 18
            i27 = 19
            i28 = 20
            i29 = 21
            i30 = 22
            i31 = 23
            i32 = 24
            i33 = 25
            i34 = 26
        elif i == ln - 9:
            i1 = ln - 8
            i2 = ln - 7
            i3 = ln - 6
            i4 = ln - 5
            i5 = ln - 4
            i6 = ln - 3
            i7 = ln - 2
            i8 = ln - 1
            i9 = 0
            i10 = 1
            i11 = 2
            i12 = 3
            i13 = 4
            i14 = 5
            i15 = 6
            i16 = 7
            i17 = 8
            i18 = 9
            i19 = 10
            i20 = 11
            i21 = 12
            i22 = 13
            i23 = 14
            i24 = 15
            i25 = 16
            i26 = 17
            i27 = 18
            i28 = 19
            i29 = 20
            i30 = 21
            i31 = 22
            i32 = 23
            i33 = 24
            i34 = 25
        elif i == ln - 10:
            i1 = ln - 9
            i2 = ln - 8
            i3 = ln - 7
            i4 = ln - 6
            i5 = ln - 5
            i6 = ln - 4
            i7 = ln - 3
            i8 = ln - 2
            i9 = ln - 1
            i10 = 0
            i11 = 1
            i12 = 2
            i13 = 3
            i14 = 4
            i15 = 5
            i16 = 6
            i17 = 7
            i18 = 8
            i19 = 9
            i20 = 10
            i21 = 11
            i22 = 12
            i23 = 13
            i24 = 14
            i25 = 15
            i26 = 16
            i27 = 17
            i28 = 18
            i29 = 19
            i30 = 20
            i31 = 21
            i32 = 22
            i33 = 23
            i34 = 24
        elif i == ln - 11:
            i1 = ln - 10
            i2 = ln - 9
            i3 = ln - 8
            i4 = ln - 7
            i5 = ln - 6
            i6 = ln - 5
            i7 = ln - 4
            i8 = ln - 3
            i9 = ln - 2
            i10 = ln - 1
            i11 = 0
            i12 = 1
            i13 = 2
            i14 = 3
            i15 = 4
            i16 = 5
            i17 = 6
            i18 = 7
            i19 = 8
            i20 = 9
            i21 = 10
            i22 = 11
            i23 = 12
            i24 = 13
            i25 = 14
            i26 = 15
            i27 = 16
            i28 = 17
            i29 = 18
            i30 = 19
            i31 = 20
            i32 = 21
            i33 = 22
            i34 = 23
        elif i == ln - 12:
            i1 = ln - 11
            i2 = ln - 10
            i3 = ln - 9
            i4 = ln - 8
            i5 = ln - 7
            i6 = ln - 6
            i7 = ln - 5
            i8 = ln - 4
            i9 = ln - 3
            i10 = ln - 2
            i11 = ln - 1
            i12 = 0
            i13 = 1
            i14 = 2
            i15 = 3
            i16 = 4
            i17 = 5
            i18 = 6
            i19 = 7
            i20 = 8
            i21 = 9
            i22 = 10
            i23 = 11
            i24 = 12
            i25 = 13
            i26 = 14
            i27 = 15
            i28 = 16
            i29 = 17
            i30 = 18
            i31 = 19
            i32 = 20
            i33 = 21
            i34 = 22
        elif i == ln - 13:
            i1 = ln - 12
            i2 = ln - 11
            i3 = ln - 10
            i4 = ln - 9
            i5 = ln - 8
            i6 = ln - 7
            i7 = ln - 6
            i8 = ln - 5
            i9 = ln - 4
            i10 = ln - 3
            i11 = ln - 2
            i12 = ln - 1
            i13 = 0
            i14 = 1
            i15 = 2
            i16 = 3
            i17 = 4
            i18 = 5
            i19 = 6
            i20 = 7
            i21 = 8
            i22 = 9
            i23 = 10
            i24 = 11
            i25 = 12
            i26 = 13
            i27 = 14
            i28 = 15
            i29 = 16
            i30 = 17
            i31 = 18
            i32 = 19
            i33 = 20
            i34 = 21
        elif i == ln - 14:
            i1 = ln - 13
            i2 = ln - 12
            i3 = ln - 11
            i4 = ln - 10
            i5 = ln - 9
            i6 = ln - 8
            i7 = ln - 7
            i8 = ln - 6
            i9 = ln - 5
            i10 = ln - 4
            i11 = ln - 3
            i12 = ln - 2
            i13 = ln - 1
            i14 = 0
            i15 = 1
            i16 = 2
            i17 = 3
            i18 = 4
            i19 = 5
            i20 = 6
            i21 = 7
            i22 = 8
            i23 = 9
            i24 = 10
            i25 = 11
            i26 = 12
            i27 = 13
            i28 = 14
            i29 = 15
            i30 = 16
            i31 = 17
            i32 = 18
            i33 = 19
            i34 = 20
        elif i == ln - 15:
            i1 = ln - 14
            i2 = ln - 13
            i3 = ln - 12
            i4 = ln - 11
            i5 = ln - 10
            i6 = ln - 9
            i7 = ln - 8
            i8 = ln - 7
            i9 = ln - 6
            i10 = ln - 5
            i11 = ln - 4
            i12 = ln - 3
            i13 = ln - 2
            i14 = ln - 1
            i15 = 0
            i16 = 1
            i17 = 2
            i18 = 3
            i19 = 4
            i20 = 5
            i21 = 6
            i22 = 7
            i23 = 8
            i24 = 9
            i25 = 10
            i26 = 11
            i27 = 12
            i28 = 13
            i29 = 14
            i30 = 15
            i31 = 16
            i32 = 17
            i33 = 18
            i34 = 19
        elif i == ln - 16:
            i1 = ln - 15
            i2 = ln - 14
            i3 = ln - 13
            i4 = ln - 12
            i5 = ln - 11
            i6 = ln - 10
            i7 = ln - 9
            i8 = ln - 8
            i9 = ln - 7
            i10 = ln - 6
            i11 = ln - 5
            i12 = ln - 4
            i13 = ln - 3
            i14 = ln - 2
            i15 = ln - 1
            i16 = 0
            i17 = 1
            i18 = 2
            i19 = 3
            i20 = 4
            i21 = 5
            i22 = 6
            i23 = 7
            i24 = 8
            i25 = 9
            i26 = 10
            i27 = 11
            i28 = 12
            i29 = 13
            i30 = 14
            i31 = 15
            i32 = 16
            i33 = 17
            i34 = 18
        elif i == ln - 17:
            i1 = ln - 16
            i2 = ln - 15
            i3 = ln - 14
            i4 = ln - 13
            i5 = ln - 12
            i6 = ln - 11
            i7 = ln - 10
            i8 = ln - 9
            i9 = ln - 8
            i10 = ln - 7
            i11 = ln - 6
            i12 = ln - 5
            i13 = ln - 4
            i14 = ln - 3
            i15 = ln - 2
            i16 = ln - 1
            i17 = 0
            i18 = 1
            i19 = 2
            i20 = 3
            i21 = 4
            i22 = 5
            i23 = 6
            i24 = 7
            i25 = 8
            i26 = 9
            i27 = 10
            i28 = 11
            i29 = 12
            i30 = 13
            i31 = 14
            i32 = 15
            i33 = 16
            i34 = 17
        elif i == ln - 18:
            i1 = ln - 17
            i2 = ln - 16
            i3 = ln - 15
            i4 = ln - 14
            i5 = ln - 13
            i6 = ln - 12
            i7 = ln - 11
            i8 = ln - 10
            i9 = ln - 9
            i10 = ln - 8
            i11 = ln - 7
            i12 = ln - 6
            i13 = ln - 5
            i14 = ln - 4
            i15 = ln - 3
            i16 = ln - 2
            i17 = ln - 1
            i18 = 0
            i19 = 1
            i20 = 2
            i21 = 3
            i22 = 4
            i23 = 5
            i24 = 6
            i25 = 7
            i26 = 8
            i27 = 9
            i28 = 10
            i29 = 11
            i30 = 12
            i31 = 13
            i32 = 14
            i33 = 15
            i34 = 16
        elif i == ln - 19:
            i1 = ln - 18
            i2 = ln - 17
            i3 = ln - 16
            i4 = ln - 15
            i5 = ln - 14
            i6 = ln - 13
            i7 = ln - 12
            i8 = ln - 11
            i9 = ln - 10
            i10 = ln - 9
            i11 = ln - 8
            i12 = ln - 7
            i13 = ln - 6
            i14 = ln - 5
            i15 = ln - 4
            i16 = ln - 3
            i17 = ln - 2
            i18 = ln - 1
            i19 = 0
            i20 = 1
            i21 = 2
            i22 = 3
            i23 = 4
            i24 = 5
            i25 = 6
            i26 = 7
            i27 = 8
            i28 = 9
            i29 = 10
            i30 = 11
            i31 = 12
            i32 = 13
            i33 = 14
            i34 = 15
        elif i == ln - 20:
            i1 = ln - 19
            i2 = ln - 18
            i3 = ln - 17
            i4 = ln - 16
            i5 = ln - 15
            i6 = ln - 14
            i7 = ln - 13
            i8 = ln - 12
            i9 = ln - 11
            i10 = ln - 10
            i11 = ln - 9
            i12 = ln - 8
            i13 = ln - 7
            i14 = ln - 6
            i15 = ln - 5
            i16 = ln - 4
            i17 = ln - 3
            i18 = ln - 2
            i19 = ln - 1
            i20 = 0
            i21 = 1
            i22 = 2
            i23 = 3
            i24 = 4
            i25 = 5
            i26 = 6
            i27 = 7
            i28 = 8
            i29 = 9
            i30 = 10
            i31 = 11
            i32 = 12
            i33 = 13
            i34 = 14
        elif i == ln - 21:
            i1 = ln - 20
            i2 = ln - 19
            i3 = ln - 18
            i4 = ln - 17
            i5 = ln - 16
            i6 = ln - 15
            i7 = ln - 14
            i8 = ln - 13
            i9 = ln - 12
            i10 = ln - 11
            i11 = ln - 10
            i12 = ln - 9
            i13 = ln - 8
            i14 = ln - 7
            i15 = ln - 6
            i16 = ln - 5
            i17 = ln - 4
            i18 = ln - 3
            i19 = ln - 2
            i20 = ln - 1
            i21 = 0
            i22 = 1
            i23 = 2
            i24 = 3
            i25 = 4
            i26 = 5
            i27 = 6
            i28 = 7
            i29 = 8
            i30 = 9
            i31 = 10
            i32 = 11
            i33 = 12
            i34 = 13
        elif i == ln - 22:
            i1 = ln - 21
            i2 = ln - 20
            i3 = ln - 19
            i4 = ln - 18
            i5 = ln - 17
            i6 = ln - 16
            i7 = ln - 15
            i8 = ln - 14
            i9 = ln - 13
            i10 = ln - 12
            i11 = ln - 11
            i12 = ln - 10
            i13 = ln - 9
            i14 = ln - 8
            i15 = ln - 7
            i16 = ln - 6
            i17 = ln - 5
            i18 = ln - 4
            i19 = ln - 3
            i20 = ln - 2
            i21 = ln - 1
            i22 = 0
            i23 = 1
            i24 = 2
            i25 = 3
            i26 = 4
            i27 = 5
            i28 = 6
            i29 = 7
            i30 = 8
            i31 = 9
            i32 = 10
            i33 = 11
            i34 = 12
        elif i == ln - 23:
            i1 = ln - 22
            i2 = ln - 21
            i3 = ln - 20
            i4 = ln - 19
            i5 = ln - 18
            i6 = ln - 17
            i7 = ln - 16
            i8 = ln - 15
            i9 = ln - 14
            i10 = ln - 13
            i11 = ln - 12
            i12 = ln - 11
            i13 = ln - 10
            i14 = ln - 9
            i15 = ln - 8
            i16 = ln - 7
            i17 = ln - 6
            i18 = ln - 5
            i19 = ln - 4
            i20 = ln - 3
            i21 = ln - 2
            i22 = ln - 1
            i23 = 0
            i24 = 1
            i25 = 2
            i26 = 3
            i27 = 4
            i28 = 5
            i29 = 6
            i30 = 7
            i31 = 8
            i32 = 9
            i33 = 10
            i34 = 11
        elif i == ln - 24:
            i1 = ln - 23
            i2 = ln - 22
            i3 = ln - 21
            i4 = ln - 20
            i5 = ln - 19
            i6 = ln - 18
            i7 = ln - 17
            i8 = ln - 16
            i9 = ln - 15
            i10 = ln - 14
            i11 = ln - 13
            i12 = ln - 12
            i13 = ln - 11
            i14 = ln - 10
            i15 = ln - 9
            i16 = ln - 8
            i17 = ln - 7
            i18 = ln - 6
            i19 = ln - 5
            i20 = ln - 4
            i21 = ln - 3
            i22 = ln - 2
            i23 = ln - 1
            i24 = 0
            i25 = 1
            i26 = 2
            i27 = 3
            i28 = 4
            i29 = 5
            i30 = 6
            i31 = 7
            i32 = 8
            i33 = 9
            i34 = 10
        elif i == ln - 25:
            i1 = ln - 24
            i2 = ln - 23
            i3 = ln - 22
            i4 = ln - 21
            i5 = ln - 20
            i6 = ln - 19
            i7 = ln - 18
            i8 = ln - 17
            i9 = ln - 16
            i10 = ln - 15
            i11 = ln - 14
            i12 = ln - 13
            i13 = ln - 12
            i14 = ln - 11
            i15 = ln - 10
            i16 = ln - 9
            i17 = ln - 8
            i18 = ln - 7
            i19 = ln - 6
            i20 = ln - 5
            i21 = ln - 4
            i22 = ln - 3
            i23 = ln - 2
            i24 = ln - 1
            i25 = 0
            i26 = 1
            i27 = 2
            i28 = 3
            i29 = 4
            i30 = 5
            i31 = 6
            i32 = 7
            i33 = 8
            i34 = 9
        elif i == ln - 26:
            i1 = ln - 25
            i2 = ln - 24
            i3 = ln - 23
            i4 = ln - 22
            i5 = ln - 21
            i6 = ln - 20
            i7 = ln - 19
            i8 = ln - 18
            i9 = ln - 17
            i10 = ln - 16
            i11 = ln - 15
            i12 = ln - 14
            i13 = ln - 13
            i14 = ln - 12
            i15 = ln - 11
            i16 = ln - 10
            i17 = ln - 9
            i18 = ln - 8
            i19 = ln - 7
            i20 = ln - 6
            i21 = ln - 5
            i22 = ln - 4
            i23 = ln - 3
            i24 = ln - 2
            i25 = ln - 1
            i26 = 0
            i27 = 1
            i28 = 2
            i29 = 3
            i30 = 4
            i31 = 5
            i32 = 6
            i33 = 7
            i34 = 8
        elif i == ln - 27:
            i1 = ln - 26
            i2 = ln - 25
            i3 = ln - 24
            i4 = ln - 23
            i5 = ln - 22
            i6 = ln - 21
            i7 = ln - 20
            i8 = ln - 19
            i9 = ln - 18
            i10 = ln - 17
            i11 = ln - 16
            i12 = ln - 15
            i13 = ln - 14
            i14 = ln - 13
            i15 = ln - 12
            i16 = ln - 11
            i17 = ln - 10
            i18 = ln - 9
            i19 = ln - 8
            i20 = ln - 7
            i21 = ln - 6
            i22 = ln - 5
            i23 = ln - 4
            i24 = ln - 3
            i25 = ln - 2
            i26 = ln - 1
            i27 = 0
            i28 = 1
            i29 = 2
            i30 = 3
            i31 = 4
            i32 = 5
            i33 = 6
            i34 = 7
        elif i == ln - 28:
            i1 = ln - 27
            i2 = ln - 26
            i3 = ln - 25
            i4 = ln - 24
            i5 = ln - 23
            i6 = ln - 22
            i7 = ln - 21
            i8 = ln - 20
            i9 = ln - 19
            i10 = ln - 18
            i11 = ln - 17
            i12 = ln - 16
            i13 = ln - 15
            i14 = ln - 14
            i15 = ln - 13
            i16 = ln - 12
            i17 = ln - 11
            i18 = ln - 10
            i19 = ln - 9
            i20 = ln - 8
            i21 = ln - 7
            i22 = ln - 6
            i23 = ln - 5
            i24 = ln - 4
            i25 = ln - 3
            i26 = ln - 2
            i27 = ln - 1
            i28 = 0
            i29 = 1
            i30 = 2
            i31 = 3
            i32 = 4
            i33 = 5
            i34 = 6
        elif i == ln - 29:
            i1 = ln - 28
            i2 = ln - 27
            i3 = ln - 26
            i4 = ln - 25
            i5 = ln - 24
            i6 = ln - 23
            i7 = ln - 22
            i8 = ln - 21
            i9 = ln - 20
            i10 = ln - 19
            i11 = ln - 18
            i12 = ln - 17
            i13 = ln - 16
            i14 = ln - 15
            i15 = ln - 14
            i16 = ln - 13
            i17 = ln - 12
            i18 = ln - 11
            i19 = ln - 10
            i20 = ln - 9
            i21 = ln - 8
            i22 = ln - 7
            i23 = ln - 6
            i24 = ln - 5
            i25 = ln - 4
            i26 = ln - 3
            i27 = ln - 2
            i28 = ln - 1
            i29 = 0
            i30 = 1
            i31 = 2
            i32 = 3
            i33 = 4
            i34 = 5
        elif i == ln - 30:
            i1 = ln - 29
            i2 = ln - 28
            i3 = ln - 27
            i4 = ln - 26
            i5 = ln - 25
            i6 = ln - 24
            i7 = ln - 23
            i8 = ln - 22
            i9 = ln - 21
            i10 = ln - 20
            i11 = ln - 19
            i12 = ln - 18
            i13 = ln - 17
            i14 = ln - 16
            i15 = ln - 15
            i16 = ln - 14
            i17 = ln - 13
            i18 = ln - 12
            i19 = ln - 11
            i20 = ln - 10
            i21 = ln - 9
            i22 = ln - 8
            i23 = ln - 7
            i24 = ln - 6
            i25 = ln - 5
            i26 = ln - 4
            i27 = ln - 3
            i28 = ln - 2
            i29 = ln - 1
            i30 = 0
            i31 = 1
            i32 = 2
            i33 = 3
            i34 = 4
        elif i == ln - 31:
            i1 = ln - 30
            i2 = ln - 29
            i3 = ln - 28
            i4 = ln - 27
            i5 = ln - 26
            i6 = ln - 25
            i7 = ln - 24
            i8 = ln - 23
            i9 = ln - 22
            i10 = ln - 21
            i11 = ln - 20
            i12 = ln - 19
            i13 = ln - 18
            i14 = ln - 17
            i15 = ln - 16
            i16 = ln - 15
            i17 = ln - 14
            i18 = ln - 13
            i19 = ln - 12
            i20 = ln - 11
            i21 = ln - 10
            i22 = ln - 9
            i23 = ln - 8
            i24 = ln - 7
            i25 = ln - 6
            i26 = ln - 5
            i27 = ln - 4
            i28 = ln - 3
            i29 = ln - 2
            i30 = ln - 1
            i31 = 0
            i32 = 1
            i33 = 2
            i34 = 3
        elif i == ln - 32:
            i1 = ln - 31
            i2 = ln - 30
            i3 = ln - 29
            i4 = ln - 28
            i5 = ln - 27
            i6 = ln - 26
            i7 = ln - 25
            i8 = ln - 24
            i9 = ln - 23
            i10 = ln - 22
            i11 = ln - 21
            i12 = ln - 20
            i13 = ln - 19
            i14 = ln - 18
            i15 = ln - 17
            i16 = ln - 16
            i17 = ln - 15
            i18 = ln - 14
            i19 = ln - 13
            i20 = ln - 12
            i21 = ln - 11
            i22 = ln - 10
            i23 = ln - 9
            i24 = ln - 8
            i25 = ln - 7
            i26 = ln - 6
            i27 = ln - 5
            i28 = ln - 4
            i29 = ln - 3
            i30 = ln - 2
            i31 = ln - 1
            i32 = 0
            i33 = 1
            i34 = 2
        elif i == ln - 33:
            i1 = ln - 32
            i2 = ln - 31
            i3 = ln - 30
            i4 = ln - 29
            i5 = ln - 28
            i6 = ln - 27
            i7 = ln - 26
            i8 = ln - 25
            i9 = ln - 24
            i10 = ln - 23
            i11 = ln - 22
            i12 = ln - 21
            i13 = ln - 20
            i14 = ln - 19
            i15 = ln - 18
            i16 = ln - 17
            i17 = ln - 16
            i18 = ln - 15
            i19 = ln - 14
            i20 = ln - 13
            i21 = ln - 12
            i22 = ln - 11
            i23 = ln - 10
            i24 = ln - 9
            i25 = ln - 8
            i26 = ln - 7
            i27 = ln - 6
            i28 = ln - 5
            i29 = ln - 4
            i30 = ln - 3
            i31 = ln - 2
            i32 = ln - 1
            i33 = 0
            i34 = 1
        elif i == ln - 34:
            i1 = ln - 33
            i2 = ln - 32
            i3 = ln - 31
            i4 = ln - 30
            i5 = ln - 29
            i6 = ln - 28
            i7 = ln - 27
            i8 = ln - 26
            i9 = ln - 25
            i10 = ln - 24
            i11 = ln - 23
            i12 = ln - 22
            i13 = ln - 21
            i14 = ln - 20
            i15 = ln - 19
            i16 = ln - 18
            i17 = ln - 17
            i18 = ln - 16
            i19 = ln - 15
            i20 = ln - 14
            i21 = ln - 13
            i22 = ln - 12
            i23 = ln - 11
            i24 = ln - 10
            i25 = ln - 9
            i26 = ln - 8
            i27 = ln - 7
            i28 = ln - 6
            i29 = ln - 5
            i30 = ln - 4
            i31 = ln - 3
            i32 = ln - 2
            i33 = ln - 1
            i34 = 0
        if i + 34 > ln:
            if ja[i - 8] == "_" or ja[i - 7] == "_" or ja[i - 6] == "_" or ja[i - 5] == "_" or ja[i - 4] == "_" or ja[
                i - 3] == "_" or ja[i - 2] == "_" or ja[i - 1] == "_" or ja[i] == "_" or ja[i1] == "_" or ja[
                i2] == "_" or ja[i3] == "_" or ja[i4] == "_" or ja[i5] == "_" or ja[i6] == "_" or ja[
                i7] == "_" or ja[i8] == "_":
                print(ja[i - 34], ja[i - 33], ja[i - 32], ja[i - 31], ja[i - 30], ja[i - 29], ja[i - 28], ja[i - 27],
                      ja[i - 26], ja[i - 25], ja[i - 24], ja[i - 23], ja[i - 22], ja[i - 21], ja[i - 20], ja[i - 19],
                      ja[i - 18], ja[i - 17], ja[i - 16], ja[i - 15], ja[i - 14], ja[i - 13], ja[i - 12], ja[i - 11],
                      ja[i - 10], ja[i - 9], ja[i - 8], ja[i - 7], ja[i - 6], ja[i - 5], ja[i - 4], ja[i - 3],
                      ja[i - 2],
                      ja[i - 1], ja[i], ja[i1], ja[i2], ja[i3], ja[i4], ja[i5], ja[i6], ja[i7],
                      ja[i8], ja[i9], ja[i10], ja[i11], ja[i12], ja[i13], ja[i14], ja[i15],
                      ja[i16], ja[i17], ja[i18], ja[i19], ja[i20], ja[i21], ja[i22], ja[i23],
                      ja[i24], ja[i25], ja[i26], ja[i27], ja[i28], ja[i29], ja[i30], ja[i31],
                      ja[i32], ja[i33], ja[i34], file=open("output4.txt", "a+"), sep="")
                strix = open("output4.txt")
                striy = strix.readlines()
                string = striy[-1]
                m_pos = 35 - string.count("_", 0, 35)
                res = [sub.replace('_', '') for sub in string]
                x = "".join(res)
                print(x[m_pos - 8:m_pos + 9], file=open("oligo.txt", "a+"))
            elif ja[i - 8] != "_" or ja[i - 7] != "_" or ja[i - 6] != "_" or ja[i - 5] != "_" or ja[i - 4] != "_" or ja[
                i - 3] != "_" or ja[i - 2] != "_" or ja[i - 1] != "_" or ja[i] != "_" or ja[i1] != "_" or ja[
                i2] != "_" or ja[i3] != "_" or ja[i4] != "_" or ja[i5] != "_" or ja[i6] != "_" or ja[
                i7] != "_" or ja[i8] != "_":
                print(ja[i - 8], ja[i - 7], ja[i - 6], ja[i - 5], ja[i - 4], ja[i - 3], ja[i - 2], ja[i - 1], ja[i],
                      ja[i1], ja[i2], ja[i3], ja[i4], ja[i5], ja[i6], ja[i7], ja[i8],
                      file=open("oligo.txt", "a+"), sep="")

        if i + 34 <= ln:
            if ja[i - 8] == "_" or ja[i - 7] == "_" or ja[i - 6] == "_" or ja[i - 5] == "_" or ja[i - 4] == "_" or ja[
                i - 3] == "_" or ja[i - 2] == "_" or ja[i - 1] == "_" or ja[i] == "_" or ja[i + 1] == "_" or ja[
                i + 2] == "_" or ja[i + 3] == "_" or ja[i + 4] == "_" or ja[i + 5] == "_" or ja[i + 6] == "_" or ja[
                i + 7] == "_" or ja[i + 8] == "_":
                print(ja[i - 34], ja[i - 33], ja[i - 32], ja[i - 31], ja[i - 30], ja[i - 29], ja[i - 28], ja[i - 27],
                      ja[i - 26], ja[i - 25], ja[i - 24], ja[i - 23], ja[i - 22], ja[i - 21], ja[i - 20], ja[i - 19],
                      ja[i - 18], ja[i - 17], ja[i - 16], ja[i - 15], ja[i - 14], ja[i - 13], ja[i - 12], ja[i - 11],
                      ja[i - 10], ja[i - 9], ja[i - 8], ja[i - 7], ja[i - 6], ja[i - 5], ja[i - 4], ja[i - 3],
                      ja[i - 2],
                      ja[i - 1], ja[i], ja[i + 1], ja[i + 2], ja[i + 3], ja[i + 4], ja[i + 5], ja[i + 6], ja[i + 7],
                      ja[i + 8], ja[i + 9], ja[i + 10], ja[i + 11], ja[i + 12], ja[i + 13], ja[i + 14], ja[i + 15],
                      ja[i + 16], ja[i + 17], ja[i + 18], ja[i + 19], ja[i + 20], ja[i + 21], ja[i + 22], ja[i + 23],
                      ja[i + 24], ja[i + 25], ja[i + 26], ja[i + 27], ja[i + 28], ja[i + 29], ja[i + 30], ja[i + 31],
                      ja[i + 32], ja[i + 33], ja[i + 34], file=open("output4.txt", "a+"), sep="")
                strix = open("output4.txt")
                striy = strix.readlines()
                string = striy[-1]
                m_pos = 35 - string.count("_", 0, 35)
                res = [sub.replace('_', '') for sub in string]
                x = "".join(res)
                print(x[m_pos - 8:m_pos + 9], file=open("oligo.txt", "a+"))
            elif ja[i - 8] != "_" or ja[i - 7] != "_" or ja[i - 6] != "_" or ja[i - 5] != "_" or ja[i - 4] != "_" or ja[
                i - 3] != "_" or ja[i - 2] != "_" or ja[i - 1] != "_" or ja[i] != "_" or ja[i + 1] != "_" or ja[
                i + 2] != "_" or ja[i + 3] != "_" or ja[i + 4] != "_" or ja[i + 5] != "_" or ja[i + 6] != "_" or ja[
                i + 7] != "_" or ja[i + 8] != "_":
                print(ja[i - 8], ja[i - 7], ja[i - 6], ja[i - 5], ja[i - 4], ja[i - 3], ja[i - 2], ja[i - 1], ja[i],
                      ja[i + 1], ja[i + 2], ja[i + 3], ja[i + 4], ja[i + 5], ja[i + 6], ja[i + 7], ja[i + 8],
                      file=open("oligo.txt", "a+"), sep="")
# generate the reverse and complement of the oligos and save fw and rv in the same file
from Bio.Seq import Seq

with open('oligo.txt') as dew:
    kas = dew.readlines()
    for lina in kas:
        seq = Seq(lina)
        print(seq.reverse_complement(), file=open("oligorv1.txt", "a+"))
os.system("sed '/^$/d' oligorv1.txt > oligorv2.txt")
os.system("paste -d '|' oligo.txt oligorv2.txt > oligos_fw_rv.txt")
olix = open("oligos_fw_rv.txt")
oliy = olix.readlines()
linelist = [line.rstrip('\n') for line in open("oligos_fw_rv.txt")]
# spades = args.input3
spades = [line.rstrip('\n') for line in open(args.input3)]
for line in linelist:
    fieldsa = line.split('|')
    sum = 0
    for linew in spades:
        if linew[0] == '>':
            fields = linew.split('_')
        if (fieldsa[0] in linew) or (fieldsa[1] in linew) and fields:
            sum += float(fields[-1])  # or float
    print(sum, file=open("frequencies_query.txt", "a+"))

# os.system("grep -E '{0}' {1} | wc -l >> frequencies_query.txt".format(line, spades))
# subject sequence processing(subject)
with open('ff.txt') as pr:
    lines = pr.readlines()
    hh = lines[4].rstrip("\n\r")
    pos = [m.start() for m in re.finditer('m', hh)]
    ln = len(hh)
    ja = lines[2].rstrip("\n\r")
    # processing of line 2 that correspond to the subject sequence
    for i in pos:
        if i == ln - 1:
            i1 = 0
            i2 = 1
            i3 = 2
            i4 = 3
            i5 = 4
            i6 = 5
            i7 = 6
            i8 = 7
            i9 = 8
            i10 = 9
            i11 = 10
            i12 = 11
            i13 = 12
            i14 = 13
            i15 = 14
            i16 = 15
            i17 = 16
            i18 = 17
            i19 = 18
            i20 = 19
            i21 = 20
            i22 = 21
            i23 = 22
            i24 = 23
            i25 = 24
            i26 = 25
            i27 = 26
            i28 = 27
            i29 = 28
            i30 = 29
            i31 = 30
            i32 = 31
            i33 = 32
            i34 = 33
        elif i == ln - 2:
            i1 = ln - 1
            i2 = 0
            i3 = 1
            i4 = 2
            i5 = 3
            i6 = 4
            i7 = 5
            i8 = 6
            i9 = 7
            i10 = 8
            i11 = 9
            i12 = 10
            i13 = 11
            i14 = 12
            i15 = 13
            i16 = 14
            i17 = 15
            i18 = 16
            i19 = 17
            i20 = 18
            i21 = 19
            i22 = 20
            i23 = 21
            i24 = 22
            i25 = 23
            i26 = 24
            i27 = 25
            i28 = 26
            i29 = 27
            i30 = 28
            i31 = 29
            i32 = 30
            i33 = 31
            i34 = 32
        elif i == ln - 3:
            i1 = ln - 2
            i2 = ln - 1
            i3 = 0
            i4 = 1
            i5 = 2
            i6 = 3
            i7 = 4
            i8 = 5
            i9 = 6
            i10 = 7
            i11 = 8
            i12 = 9
            i13 = 10
            i14 = 11
            i15 = 12
            i16 = 13
            i17 = 14
            i18 = 15
            i19 = 16
            i20 = 17
            i21 = 18
            i22 = 19
            i23 = 20
            i24 = 21
            i25 = 22
            i26 = 23
            i27 = 24
            i28 = 25
            i29 = 26
            i30 = 27
            i31 = 28
            i32 = 29
            i33 = 30
            i34 = 31
        elif i == ln - 4:
            i1 = ln - 3
            i2 = ln - 2
            i3 = ln - 1
            i4 = 0
            i5 = 1
            i6 = 2
            i7 = 3
            i8 = 4
            i9 = 5
            i10 = 6
            i11 = 7
            i12 = 8
            i13 = 9
            i14 = 10
            i15 = 11
            i16 = 12
            i17 = 13
            i18 = 14
            i19 = 15
            i20 = 16
            i21 = 17
            i22 = 18
            i23 = 19
            i24 = 20
            i25 = 21
            i26 = 22
            i27 = 23
            i28 = 24
            i29 = 25
            i30 = 26
            i31 = 27
            i32 = 28
            i33 = 29
            i34 = 30
        elif i == ln - 5:
            i1 = ln - 4
            i2 = ln - 3
            i3 = ln - 2
            i4 = ln - 1
            i5 = 0
            i6 = 1
            i7 = 2
            i8 = 3
            i9 = 4
            i10 = 5
            i11 = 6
            i12 = 7
            i13 = 8
            i14 = 9
            i15 = 10
            i16 = 11
            i17 = 12
            i18 = 13
            i19 = 14
            i20 = 15
            i21 = 16
            i22 = 17
            i23 = 18
            i24 = 19
            i25 = 20
            i26 = 21
            i27 = 22
            i28 = 23
            i29 = 24
            i30 = 25
            i31 = 26
            i32 = 27
            i33 = 28
            i34 = 29
        elif i == ln - 6:
            i1 = ln - 5
            i2 = ln - 4
            i3 = ln - 3
            i4 = ln - 2
            i5 = ln - 1
            i6 = 0
            i7 = 1
            i8 = 2
            i9 = 3
            i10 = 4
            i11 = 5
            i12 = 6
            i13 = 7
            i14 = 8
            i15 = 9
            i16 = 10
            i17 = 11
            i18 = 12
            i19 = 13
            i20 = 14
            i21 = 15
            i22 = 16
            i23 = 17
            i24 = 18
            i25 = 19
            i26 = 20
            i27 = 21
            i28 = 22
            i29 = 23
            i30 = 24
            i31 = 25
            i32 = 26
            i33 = 27
            i34 = 28
        elif i == ln - 7:
            i1 = ln - 6
            i2 = ln - 5
            i3 = ln - 4
            i4 = ln - 3
            i5 = ln - 2
            i6 = ln - 1
            i7 = 0
            i8 = 1
            i9 = 2
            i10 = 3
            i11 = 4
            i12 = 5
            i13 = 6
            i14 = 7
            i15 = 8
            i16 = 9
            i17 = 10
            i18 = 11
            i19 = 12
            i20 = 13
            i21 = 14
            i22 = 15
            i23 = 16
            i24 = 17
            i25 = 18
            i26 = 19
            i27 = 20
            i28 = 21
            i29 = 22
            i30 = 23
            i31 = 24
            i32 = 25
            i33 = 26
            i34 = 27
        elif i == ln - 8:
            i1 = ln - 7
            i2 = ln - 6
            i3 = ln - 5
            i4 = ln - 4
            i5 = ln - 3
            i6 = ln - 2
            i7 = ln - 1
            i8 = 0
            i9 = 1
            i10 = 2
            i11 = 3
            i12 = 4
            i13 = 5
            i14 = 6
            i15 = 7
            i16 = 8
            i17 = 9
            i18 = 10
            i19 = 11
            i20 = 12
            i21 = 13
            i22 = 14
            i23 = 15
            i24 = 16
            i25 = 17
            i26 = 18
            i27 = 19
            i28 = 20
            i29 = 21
            i30 = 22
            i31 = 23
            i32 = 24
            i33 = 25
            i34 = 26
        elif i == ln - 9:
            i1 = ln - 8
            i2 = ln - 7
            i3 = ln - 6
            i4 = ln - 5
            i5 = ln - 4
            i6 = ln - 3
            i7 = ln - 2
            i8 = ln - 1
            i9 = 0
            i10 = 1
            i11 = 2
            i12 = 3
            i13 = 4
            i14 = 5
            i15 = 6
            i16 = 7
            i17 = 8
            i18 = 9
            i19 = 10
            i20 = 11
            i21 = 12
            i22 = 13
            i23 = 14
            i24 = 15
            i25 = 16
            i26 = 17
            i27 = 18
            i28 = 19
            i29 = 20
            i30 = 21
            i31 = 22
            i32 = 23
            i33 = 24
            i34 = 25
        elif i == ln - 10:
            i1 = ln - 9
            i2 = ln - 8
            i3 = ln - 7
            i4 = ln - 6
            i5 = ln - 5
            i6 = ln - 4
            i7 = ln - 3
            i8 = ln - 2
            i9 = ln - 1
            i10 = 0
            i11 = 1
            i12 = 2
            i13 = 3
            i14 = 4
            i15 = 5
            i16 = 6
            i17 = 7
            i18 = 8
            i19 = 9
            i20 = 10
            i21 = 11
            i22 = 12
            i23 = 13
            i24 = 14
            i25 = 15
            i26 = 16
            i27 = 17
            i28 = 18
            i29 = 19
            i30 = 20
            i31 = 21
            i32 = 22
            i33 = 23
            i34 = 24
        elif i == ln - 11:
            i1 = ln - 10
            i2 = ln - 9
            i3 = ln - 8
            i4 = ln - 7
            i5 = ln - 6
            i6 = ln - 5
            i7 = ln - 4
            i8 = ln - 3
            i9 = ln - 2
            i10 = ln - 1
            i11 = 0
            i12 = 1
            i13 = 2
            i14 = 3
            i15 = 4
            i16 = 5
            i17 = 6
            i18 = 7
            i19 = 8
            i20 = 9
            i21 = 10
            i22 = 11
            i23 = 12
            i24 = 13
            i25 = 14
            i26 = 15
            i27 = 16
            i28 = 17
            i29 = 18
            i30 = 19
            i31 = 20
            i32 = 21
            i33 = 22
            i34 = 23
        elif i == ln - 12:
            i1 = ln - 11
            i2 = ln - 10
            i3 = ln - 9
            i4 = ln - 8
            i5 = ln - 7
            i6 = ln - 6
            i7 = ln - 5
            i8 = ln - 4
            i9 = ln - 3
            i10 = ln - 2
            i11 = ln - 1
            i12 = 0
            i13 = 1
            i14 = 2
            i15 = 3
            i16 = 4
            i17 = 5
            i18 = 6
            i19 = 7
            i20 = 8
            i21 = 9
            i22 = 10
            i23 = 11
            i24 = 12
            i25 = 13
            i26 = 14
            i27 = 15
            i28 = 16
            i29 = 17
            i30 = 18
            i31 = 19
            i32 = 20
            i33 = 21
            i34 = 22
        elif i == ln - 13:
            i1 = ln - 12
            i2 = ln - 11
            i3 = ln - 10
            i4 = ln - 9
            i5 = ln - 8
            i6 = ln - 7
            i7 = ln - 6
            i8 = ln - 5
            i9 = ln - 4
            i10 = ln - 3
            i11 = ln - 2
            i12 = ln - 1
            i13 = 0
            i14 = 1
            i15 = 2
            i16 = 3
            i17 = 4
            i18 = 5
            i19 = 6
            i20 = 7
            i21 = 8
            i22 = 9
            i23 = 10
            i24 = 11
            i25 = 12
            i26 = 13
            i27 = 14
            i28 = 15
            i29 = 16
            i30 = 17
            i31 = 18
            i32 = 19
            i33 = 20
            i34 = 21
        elif i == ln - 14:
            i1 = ln - 13
            i2 = ln - 12
            i3 = ln - 11
            i4 = ln - 10
            i5 = ln - 9
            i6 = ln - 8
            i7 = ln - 7
            i8 = ln - 6
            i9 = ln - 5
            i10 = ln - 4
            i11 = ln - 3
            i12 = ln - 2
            i13 = ln - 1
            i14 = 0
            i15 = 1
            i16 = 2
            i17 = 3
            i18 = 4
            i19 = 5
            i20 = 6
            i21 = 7
            i22 = 8
            i23 = 9
            i24 = 10
            i25 = 11
            i26 = 12
            i27 = 13
            i28 = 14
            i29 = 15
            i30 = 16
            i31 = 17
            i32 = 18
            i33 = 19
            i34 = 20
        elif i == ln - 15:
            i1 = ln - 14
            i2 = ln - 13
            i3 = ln - 12
            i4 = ln - 11
            i5 = ln - 10
            i6 = ln - 9
            i7 = ln - 8
            i8 = ln - 7
            i9 = ln - 6
            i10 = ln - 5
            i11 = ln - 4
            i12 = ln - 3
            i13 = ln - 2
            i14 = ln - 1
            i15 = 0
            i16 = 1
            i17 = 2
            i18 = 3
            i19 = 4
            i20 = 5
            i21 = 6
            i22 = 7
            i23 = 8
            i24 = 9
            i25 = 10
            i26 = 11
            i27 = 12
            i28 = 13
            i29 = 14
            i30 = 15
            i31 = 16
            i32 = 17
            i33 = 18
            i34 = 19
        elif i == ln - 16:
            i1 = ln - 15
            i2 = ln - 14
            i3 = ln - 13
            i4 = ln - 12
            i5 = ln - 11
            i6 = ln - 10
            i7 = ln - 9
            i8 = ln - 8
            i9 = ln - 7
            i10 = ln - 6
            i11 = ln - 5
            i12 = ln - 4
            i13 = ln - 3
            i14 = ln - 2
            i15 = ln - 1
            i16 = 0
            i17 = 1
            i18 = 2
            i19 = 3
            i20 = 4
            i21 = 5
            i22 = 6
            i23 = 7
            i24 = 8
            i25 = 9
            i26 = 10
            i27 = 11
            i28 = 12
            i29 = 13
            i30 = 14
            i31 = 15
            i32 = 16
            i33 = 17
            i34 = 18
        elif i == ln - 17:
            i1 = ln - 16
            i2 = ln - 15
            i3 = ln - 14
            i4 = ln - 13
            i5 = ln - 12
            i6 = ln - 11
            i7 = ln - 10
            i8 = ln - 9
            i9 = ln - 8
            i10 = ln - 7
            i11 = ln - 6
            i12 = ln - 5
            i13 = ln - 4
            i14 = ln - 3
            i15 = ln - 2
            i16 = ln - 1
            i17 = 0
            i18 = 1
            i19 = 2
            i20 = 3
            i21 = 4
            i22 = 5
            i23 = 6
            i24 = 7
            i25 = 8
            i26 = 9
            i27 = 10
            i28 = 11
            i29 = 12
            i30 = 13
            i31 = 14
            i32 = 15
            i33 = 16
            i34 = 17
        elif i == ln - 18:
            i1 = ln - 17
            i2 = ln - 16
            i3 = ln - 15
            i4 = ln - 14
            i5 = ln - 13
            i6 = ln - 12
            i7 = ln - 11
            i8 = ln - 10
            i9 = ln - 9
            i10 = ln - 8
            i11 = ln - 7
            i12 = ln - 6
            i13 = ln - 5
            i14 = ln - 4
            i15 = ln - 3
            i16 = ln - 2
            i17 = ln - 1
            i18 = 0
            i19 = 1
            i20 = 2
            i21 = 3
            i22 = 4
            i23 = 5
            i24 = 6
            i25 = 7
            i26 = 8
            i27 = 9
            i28 = 10
            i29 = 11
            i30 = 12
            i31 = 13
            i32 = 14
            i33 = 15
            i34 = 16
        elif i == ln - 19:
            i1 = ln - 18
            i2 = ln - 17
            i3 = ln - 16
            i4 = ln - 15
            i5 = ln - 14
            i6 = ln - 13
            i7 = ln - 12
            i8 = ln - 11
            i9 = ln - 10
            i10 = ln - 9
            i11 = ln - 8
            i12 = ln - 7
            i13 = ln - 6
            i14 = ln - 5
            i15 = ln - 4
            i16 = ln - 3
            i17 = ln - 2
            i18 = ln - 1
            i19 = 0
            i20 = 1
            i21 = 2
            i22 = 3
            i23 = 4
            i24 = 5
            i25 = 6
            i26 = 7
            i27 = 8
            i28 = 9
            i29 = 10
            i30 = 11
            i31 = 12
            i32 = 13
            i33 = 14
            i34 = 15
        elif i == ln - 20:
            i1 = ln - 19
            i2 = ln - 18
            i3 = ln - 17
            i4 = ln - 16
            i5 = ln - 15
            i6 = ln - 14
            i7 = ln - 13
            i8 = ln - 12
            i9 = ln - 11
            i10 = ln - 10
            i11 = ln - 9
            i12 = ln - 8
            i13 = ln - 7
            i14 = ln - 6
            i15 = ln - 5
            i16 = ln - 4
            i17 = ln - 3
            i18 = ln - 2
            i19 = ln - 1
            i20 = 0
            i21 = 1
            i22 = 2
            i23 = 3
            i24 = 4
            i25 = 5
            i26 = 6
            i27 = 7
            i28 = 8
            i29 = 9
            i30 = 10
            i31 = 11
            i32 = 12
            i33 = 13
            i34 = 14
        elif i == ln - 21:
            i1 = ln - 20
            i2 = ln - 19
            i3 = ln - 18
            i4 = ln - 17
            i5 = ln - 16
            i6 = ln - 15
            i7 = ln - 14
            i8 = ln - 13
            i9 = ln - 12
            i10 = ln - 11
            i11 = ln - 10
            i12 = ln - 9
            i13 = ln - 8
            i14 = ln - 7
            i15 = ln - 6
            i16 = ln - 5
            i17 = ln - 4
            i18 = ln - 3
            i19 = ln - 2
            i20 = ln - 1
            i21 = 0
            i22 = 1
            i23 = 2
            i24 = 3
            i25 = 4
            i26 = 5
            i27 = 6
            i28 = 7
            i29 = 8
            i30 = 9
            i31 = 10
            i32 = 11
            i33 = 12
            i34 = 13
        elif i == ln - 22:
            i1 = ln - 21
            i2 = ln - 20
            i3 = ln - 19
            i4 = ln - 18
            i5 = ln - 17
            i6 = ln - 16
            i7 = ln - 15
            i8 = ln - 14
            i9 = ln - 13
            i10 = ln - 12
            i11 = ln - 11
            i12 = ln - 10
            i13 = ln - 9
            i14 = ln - 8
            i15 = ln - 7
            i16 = ln - 6
            i17 = ln - 5
            i18 = ln - 4
            i19 = ln - 3
            i20 = ln - 2
            i21 = ln - 1
            i22 = 0
            i23 = 1
            i24 = 2
            i25 = 3
            i26 = 4
            i27 = 5
            i28 = 6
            i29 = 7
            i30 = 8
            i31 = 9
            i32 = 10
            i33 = 11
            i34 = 12
        elif i == ln - 23:
            i1 = ln - 22
            i2 = ln - 21
            i3 = ln - 20
            i4 = ln - 19
            i5 = ln - 18
            i6 = ln - 17
            i7 = ln - 16
            i8 = ln - 15
            i9 = ln - 14
            i10 = ln - 13
            i11 = ln - 12
            i12 = ln - 11
            i13 = ln - 10
            i14 = ln - 9
            i15 = ln - 8
            i16 = ln - 7
            i17 = ln - 6
            i18 = ln - 5
            i19 = ln - 4
            i20 = ln - 3
            i21 = ln - 2
            i22 = ln - 1
            i23 = 0
            i24 = 1
            i25 = 2
            i26 = 3
            i27 = 4
            i28 = 5
            i29 = 6
            i30 = 7
            i31 = 8
            i32 = 9
            i33 = 10
            i34 = 11
        elif i == ln - 24:
            i1 = ln - 23
            i2 = ln - 22
            i3 = ln - 21
            i4 = ln - 20
            i5 = ln - 19
            i6 = ln - 18
            i7 = ln - 17
            i8 = ln - 16
            i9 = ln - 15
            i10 = ln - 14
            i11 = ln - 13
            i12 = ln - 12
            i13 = ln - 11
            i14 = ln - 10
            i15 = ln - 9
            i16 = ln - 8
            i17 = ln - 7
            i18 = ln - 6
            i19 = ln - 5
            i20 = ln - 4
            i21 = ln - 3
            i22 = ln - 2
            i23 = ln - 1
            i24 = 0
            i25 = 1
            i26 = 2
            i27 = 3
            i28 = 4
            i29 = 5
            i30 = 6
            i31 = 7
            i32 = 8
            i33 = 9
            i34 = 10
        elif i == ln - 25:
            i1 = ln - 24
            i2 = ln - 23
            i3 = ln - 22
            i4 = ln - 21
            i5 = ln - 20
            i6 = ln - 19
            i7 = ln - 18
            i8 = ln - 17
            i9 = ln - 16
            i10 = ln - 15
            i11 = ln - 14
            i12 = ln - 13
            i13 = ln - 12
            i14 = ln - 11
            i15 = ln - 10
            i16 = ln - 9
            i17 = ln - 8
            i18 = ln - 7
            i19 = ln - 6
            i20 = ln - 5
            i21 = ln - 4
            i22 = ln - 3
            i23 = ln - 2
            i24 = ln - 1
            i25 = 0
            i26 = 1
            i27 = 2
            i28 = 3
            i29 = 4
            i30 = 5
            i31 = 6
            i32 = 7
            i33 = 8
            i34 = 9
        elif i == ln - 26:
            i1 = ln - 25
            i2 = ln - 24
            i3 = ln - 23
            i4 = ln - 22
            i5 = ln - 21
            i6 = ln - 20
            i7 = ln - 19
            i8 = ln - 18
            i9 = ln - 17
            i10 = ln - 16
            i11 = ln - 15
            i12 = ln - 14
            i13 = ln - 13
            i14 = ln - 12
            i15 = ln - 11
            i16 = ln - 10
            i17 = ln - 9
            i18 = ln - 8
            i19 = ln - 7
            i20 = ln - 6
            i21 = ln - 5
            i22 = ln - 4
            i23 = ln - 3
            i24 = ln - 2
            i25 = ln - 1
            i26 = 0
            i27 = 1
            i28 = 2
            i29 = 3
            i30 = 4
            i31 = 5
            i32 = 6
            i33 = 7
            i34 = 8
        elif i == ln - 27:
            i1 = ln - 26
            i2 = ln - 25
            i3 = ln - 24
            i4 = ln - 23
            i5 = ln - 22
            i6 = ln - 21
            i7 = ln - 20
            i8 = ln - 19
            i9 = ln - 18
            i10 = ln - 17
            i11 = ln - 16
            i12 = ln - 15
            i13 = ln - 14
            i14 = ln - 13
            i15 = ln - 12
            i16 = ln - 11
            i17 = ln - 10
            i18 = ln - 9
            i19 = ln - 8
            i20 = ln - 7
            i21 = ln - 6
            i22 = ln - 5
            i23 = ln - 4
            i24 = ln - 3
            i25 = ln - 2
            i26 = ln - 1
            i27 = 0
            i28 = 1
            i29 = 2
            i30 = 3
            i31 = 4
            i32 = 5
            i33 = 6
            i34 = 7
        elif i == ln - 28:
            i1 = ln - 27
            i2 = ln - 26
            i3 = ln - 25
            i4 = ln - 24
            i5 = ln - 23
            i6 = ln - 22
            i7 = ln - 21
            i8 = ln - 20
            i9 = ln - 19
            i10 = ln - 18
            i11 = ln - 17
            i12 = ln - 16
            i13 = ln - 15
            i14 = ln - 14
            i15 = ln - 13
            i16 = ln - 12
            i17 = ln - 11
            i18 = ln - 10
            i19 = ln - 9
            i20 = ln - 8
            i21 = ln - 7
            i22 = ln - 6
            i23 = ln - 5
            i24 = ln - 4
            i25 = ln - 3
            i26 = ln - 2
            i27 = ln - 1
            i28 = 0
            i29 = 1
            i30 = 2
            i31 = 3
            i32 = 4
            i33 = 5
            i34 = 6
        elif i == ln - 29:
            i1 = ln - 28
            i2 = ln - 27
            i3 = ln - 26
            i4 = ln - 25
            i5 = ln - 24
            i6 = ln - 23
            i7 = ln - 22
            i8 = ln - 21
            i9 = ln - 20
            i10 = ln - 19
            i11 = ln - 18
            i12 = ln - 17
            i13 = ln - 16
            i14 = ln - 15
            i15 = ln - 14
            i16 = ln - 13
            i17 = ln - 12
            i18 = ln - 11
            i19 = ln - 10
            i20 = ln - 9
            i21 = ln - 8
            i22 = ln - 7
            i23 = ln - 6
            i24 = ln - 5
            i25 = ln - 4
            i26 = ln - 3
            i27 = ln - 2
            i28 = ln - 1
            i29 = 0
            i30 = 1
            i31 = 2
            i32 = 3
            i33 = 4
            i34 = 5
        elif i == ln - 30:
            i1 = ln - 29
            i2 = ln - 28
            i3 = ln - 27
            i4 = ln - 26
            i5 = ln - 25
            i6 = ln - 24
            i7 = ln - 23
            i8 = ln - 22
            i9 = ln - 21
            i10 = ln - 20
            i11 = ln - 19
            i12 = ln - 18
            i13 = ln - 17
            i14 = ln - 16
            i15 = ln - 15
            i16 = ln - 14
            i17 = ln - 13
            i18 = ln - 12
            i19 = ln - 11
            i20 = ln - 10
            i21 = ln - 9
            i22 = ln - 8
            i23 = ln - 7
            i24 = ln - 6
            i25 = ln - 5
            i26 = ln - 4
            i27 = ln - 3
            i28 = ln - 2
            i29 = ln - 1
            i30 = 0
            i31 = 1
            i32 = 2
            i33 = 3
            i34 = 4
        elif i == ln - 31:
            i1 = ln - 30
            i2 = ln - 29
            i3 = ln - 28
            i4 = ln - 27
            i5 = ln - 26
            i6 = ln - 25
            i7 = ln - 24
            i8 = ln - 23
            i9 = ln - 22
            i10 = ln - 21
            i11 = ln - 20
            i12 = ln - 19
            i13 = ln - 18
            i14 = ln - 17
            i15 = ln - 16
            i16 = ln - 15
            i17 = ln - 14
            i18 = ln - 13
            i19 = ln - 12
            i20 = ln - 11
            i21 = ln - 10
            i22 = ln - 9
            i23 = ln - 8
            i24 = ln - 7
            i25 = ln - 6
            i26 = ln - 5
            i27 = ln - 4
            i28 = ln - 3
            i29 = ln - 2
            i30 = ln - 1
            i31 = 0
            i32 = 1
            i33 = 2
            i34 = 3
        elif i == ln - 32:
            i1 = ln - 31
            i2 = ln - 30
            i3 = ln - 29
            i4 = ln - 28
            i5 = ln - 27
            i6 = ln - 26
            i7 = ln - 25
            i8 = ln - 24
            i9 = ln - 23
            i10 = ln - 22
            i11 = ln - 21
            i12 = ln - 20
            i13 = ln - 19
            i14 = ln - 18
            i15 = ln - 17
            i16 = ln - 16
            i17 = ln - 15
            i18 = ln - 14
            i19 = ln - 13
            i20 = ln - 12
            i21 = ln - 11
            i22 = ln - 10
            i23 = ln - 9
            i24 = ln - 8
            i25 = ln - 7
            i26 = ln - 6
            i27 = ln - 5
            i28 = ln - 4
            i29 = ln - 3
            i30 = ln - 2
            i31 = ln - 1
            i32 = 0
            i33 = 1
            i34 = 2
        elif i == ln - 33:
            i1 = ln - 32
            i2 = ln - 31
            i3 = ln - 30
            i4 = ln - 29
            i5 = ln - 28
            i6 = ln - 27
            i7 = ln - 26
            i8 = ln - 25
            i9 = ln - 24
            i10 = ln - 23
            i11 = ln - 22
            i12 = ln - 21
            i13 = ln - 20
            i14 = ln - 19
            i15 = ln - 18
            i16 = ln - 17
            i17 = ln - 16
            i18 = ln - 15
            i19 = ln - 14
            i20 = ln - 13
            i21 = ln - 12
            i22 = ln - 11
            i23 = ln - 10
            i24 = ln - 9
            i25 = ln - 8
            i26 = ln - 7
            i27 = ln - 6
            i28 = ln - 5
            i29 = ln - 4
            i30 = ln - 3
            i31 = ln - 2
            i32 = ln - 1
            i33 = 0
            i34 = 1
        elif i == ln - 34:
            i1 = ln - 33
            i2 = ln - 32
            i3 = ln - 31
            i4 = ln - 30
            i5 = ln - 29
            i6 = ln - 28
            i7 = ln - 27
            i8 = ln - 26
            i9 = ln - 25
            i10 = ln - 24
            i11 = ln - 23
            i12 = ln - 22
            i13 = ln - 21
            i14 = ln - 20
            i15 = ln - 19
            i16 = ln - 18
            i17 = ln - 17
            i18 = ln - 16
            i19 = ln - 15
            i20 = ln - 14
            i21 = ln - 13
            i22 = ln - 12
            i23 = ln - 11
            i24 = ln - 10
            i25 = ln - 9
            i26 = ln - 8
            i27 = ln - 7
            i28 = ln - 6
            i29 = ln - 5
            i30 = ln - 4
            i31 = ln - 3
            i32 = ln - 2
            i33 = ln - 1
            i34 = 0
        if i + 34 > ln:
            if ja[i - 8] == "_" or ja[i - 7] == "_" or ja[i - 6] == "_" or ja[i - 5] == "_" or ja[i - 4] == "_" or ja[
                i - 3] == "_" or ja[i - 2] == "_" or ja[i - 1] == "_" or ja[i] == "_" or ja[i1] == "_" or ja[
                i2] == "_" or ja[i3] == "_" or ja[i4] == "_" or ja[i5] == "_" or ja[i6] == "_" or ja[
                i7] == "_" or ja[i8] == "_":
                print(ja[i - 34], ja[i - 33], ja[i - 32], ja[i - 31], ja[i - 30], ja[i - 29], ja[i - 28], ja[i - 27],
                      ja[i - 26], ja[i - 25], ja[i - 24], ja[i - 23], ja[i - 22], ja[i - 21], ja[i - 20], ja[i - 19],
                      ja[i - 18], ja[i - 17], ja[i - 16], ja[i - 15], ja[i - 14], ja[i - 13], ja[i - 12], ja[i - 11],
                      ja[i - 10], ja[i - 9], ja[i - 8], ja[i - 7], ja[i - 6], ja[i - 5], ja[i - 4], ja[i - 3],
                      ja[i - 2],
                      ja[i - 1], ja[i], ja[i1], ja[i2], ja[i3], ja[i4], ja[i5], ja[i6], ja[i7],
                      ja[i8], ja[i9], ja[i10], ja[i11], ja[i12], ja[i13], ja[i14], ja[i15],
                      ja[i16], ja[i17], ja[i18], ja[i19], ja[i20], ja[i21], ja[i22], ja[i23],
                      ja[i24], ja[i25], ja[i26], ja[i27], ja[i28], ja[i29], ja[i30], ja[i31],
                      ja[i32], ja[i33], ja[i34], file=open("output4b.txt", "a+"), sep="")
                strix = open("output4b.txt")
                striy = strix.readlines()
                string = striy[-1]
                m_pos = 35 - string.count("_", 0, 35)
                res = [sub.replace('_', '') for sub in string]
                x = "".join(res)
                print(x[m_pos - 8:m_pos + 9], file=open("oligob.txt", "a+"))
            elif ja[i - 8] != "_" or ja[i - 7] != "_" or ja[i - 6] != "_" or ja[i - 5] != "_" or ja[i - 4] != "_" or ja[
                i - 3] != "_" or ja[i - 2] != "_" or ja[i - 1] != "_" or ja[i] != "_" or ja[i1] != "_" or ja[
                i2] != "_" or ja[i3] != "_" or ja[i4] != "_" or ja[i5] != "_" or ja[i6] != "_" or ja[
                i7] != "_" or ja[i8] != "_":
                print(ja[i - 8], ja[i - 7], ja[i - 6], ja[i - 5], ja[i - 4], ja[i - 3], ja[i - 2], ja[i - 1], ja[i],
                      ja[i1], ja[i2], ja[i3], ja[i4], ja[i5], ja[i6], ja[i7], ja[i8],
                      file=open("oligob.txt", "a+"), sep="")
        if i + 34 <= ln:
            if ja[i - 8] == "_" or ja[i - 7] == "_" or ja[i - 6] == "_" or ja[i - 5] == "_" or ja[i - 4] == "_" or ja[
                i - 3] == "_" or ja[i - 2] == "_" or ja[i - 1] == "_" or ja[i] == "_" or ja[i + 1] == "_" or \
                    ja[i + 2] == "_" or ja[i + 3] == "_" or ja[i + 4] == "_" or ja[i + 5] == "_" or ja[i + 6] == "_" \
                    or ja[i + 7] == "_" or ja[i + 8] == "_":
                print(ja[i - 34], ja[i - 33], ja[i - 32], ja[i - 31], ja[i - 30], ja[i - 29], ja[i - 28], ja[i - 27],
                      ja[i - 26], ja[i - 25], ja[i - 24], ja[i - 23], ja[i - 22], ja[i - 21], ja[i - 20], ja[i - 19],
                      ja[i - 18], ja[i - 17], ja[i - 16], ja[i - 15], ja[i - 14], ja[i - 13], ja[i - 12], ja[i - 11],
                      ja[i - 10], ja[i - 9], ja[i - 8], ja[i - 7], ja[i - 6], ja[i - 5], ja[i - 4], ja[i - 3],
                      ja[i - 2],
                      ja[i - 1], ja[i], ja[i + 1], ja[i + 2], ja[i + 3], ja[i + 4], ja[i + 5], ja[i + 6], ja[i + 7],
                      ja[i + 8], ja[i + 9], ja[i + 10], ja[i + 11], ja[i + 12], ja[i + 13], ja[i + 14], ja[i + 15],
                      ja[i + 16], ja[i + 17], ja[i + 18], ja[i + 19], ja[i + 20], ja[i + 21], ja[i + 22], ja[i + 23],
                      ja[i + 24], ja[i + 25], ja[i + 26], ja[i + 27], ja[i + 28], ja[i + 29], ja[i + 30], ja[i + 31],
                      ja[i + 32], ja[i + 33], ja[i + 34], file=open("output4b.txt", "a+"), sep="")
                strix = open("output4b.txt")
                striy = strix.readlines()
                string = striy[-1]
                m_pos = 35 - string.count("_", 0, 35)
                res = [sub.replace('_', '') for sub in string]
                x = "".join(res)
                print(x[m_pos - 8:m_pos + 9], file=open("oligob.txt", "a+"))
            elif ja[i - 8] != "_" or ja[i - 7] != "_" or ja[i - 6] != "_" or ja[i - 5] != "_" or ja[i - 4] != "_" or ja[
                i - 3] != "_" or ja[i - 2] != "_" or ja[i - 1] != "_" or ja[i] != "_" or ja[i + 1] != "_" or ja[
                i + 2] != "_" or ja[i + 3] != "_" or ja[i + 4] != "_" or ja[i + 5] != "_" or ja[i + 6] != "_" or ja[
                i + 7] \
                    != "_" or ja[i + 8] != "_":
                print(ja[i - 8], ja[i - 7], ja[i - 6], ja[i - 5], ja[i - 4], ja[i - 3], ja[i - 2], ja[i - 1], ja[i],
                      ja[i + 1], ja[i + 2], ja[i + 3], ja[i + 4], ja[i + 5], ja[i + 6], ja[i + 7], ja[i + 8],
                      file=open("oligob.txt", "a+"), sep="")
# generate the reverse and complement of the oligos and save fw and rv in the same file
from Bio.Seq import Seq

with open('oligob.txt') as dew:
    kas = dew.readlines()
    for lina in kas:
        seq = Seq(lina)
        print(seq.reverse_complement(), file=open("oligorv1b.txt", "a+"))
os.system("sed '/^$/d' oligorv1b.txt > oligorv2b.txt")
os.system("paste -d '|' oligob.txt oligorv2b.txt > oligos_fw_rvb.txt")
olix = open("oligos_fw_rvb.txt")
oliy = olix.readlines()
linelist = [line.rstrip('\n') for line in open("oligos_fw_rvb.txt")]
# spades = args.input3
spades = [line.rstrip('\n') for line in open(args.input3)]
for line in linelist:
    fieldsa = line.split('|')
    sum = 0
    for linew in spades:
        if linew[0] == '>':
            fields = linew.split('_')
        if (fieldsa[0] in linew) or (fieldsa[1] in linew) and fields:
            sum += float(fields[-1])  # or float
    print(sum, file=open("frequencies_subject.txt", "a+"))
#    os.system("grep -E '{0}' {1} | wc -l >> frequencies_subject.txt".format(line, spades))
# combine frequencies in 1 file and produce 1 column of the consensus
data1 = pd.read_csv('frequencies_query.txt', names=["Sequence"])
data2 = pd.read_csv('frequencies_subject.txt', names=["Sequence"])
data1["1"] = np.where(data1['Sequence'] < data2['Sequence'], 's',
                      np.where(data1['Sequence'] > data2['Sequence'], 'q', 'e'))
ghe = data1["1"]
ghe.to_csv("cons.txt", header=False, index=False)
# create a file with the correct consensus in a line
write_file = open("ff.txt", "r+")
line1 = list(write_file.readlines())
line = line1[4]
val = open("cons.txt")
vals = val.readlines()
i = -1
for j in line:
    if j == "*":
        print("*", file=open("fino.txt", "a+"))
    elif j == "m":
        print(vals[i + 1], file=open("fino.txt", "a+"))
        i += 1
os.system("tr -d '\n\r' < fino.txt >> finox.txt")
# really last
# convert string in list of characters
pr = open('ff.txt', "r+")
lines = pr.readlines()
query = lines[0].rstrip("\n\r")
subject = lines[2].rstrip("\n\r")
cons = open('finox.txt', 'r+').readline()
qu = list(query)
su = list(subject)
co = list(cons)


# generate final contig
def f(jx, k, l):
    if jx == "*":
        return l
    elif jx == "q":
        return l
    elif jx == "s":
        return k
    elif jx == "e":
        return l


for (a, b, c) in zip(co, su, qu):
    print(f(a, b, c), file=open("last_cons.txt", "a+"))

with open(r"last_consensus.txt", "w+") as f:
    f.write(">consensus\n")
os.system('tr -d "\n\r" < last_cons.txt >> last_consensus.txt')
os.system("sed -e 's/_//g' last_consensus.txt > consensus.fa")
print("Finish")
print("the polished sequence is in the file: consensus.fa")
answer = input(
    "Do you want to keep process files (default=False)? If yes, all the files inside the working dir will be destroied exept the output file [yes/no]: ")
if answer == "yes":
    print("all the process files will be kept!")
elif answer == "no":
    print("all the process files will be discarded exept the consensus.fa file!")
    os.system("find . -type f -not -name 'consensus.fa' -delete")
else:
    print("Please enter yes or no.")
