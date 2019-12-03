#!/usr/bin/env python
# DraftPolisher: a tool for the polishing of draft sequences
#     Copyright (C) 2019 Rosario Nicola Brancaccio
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
# chmod u+x DraftPolisher_circular.py
# # input files format:
# >QRY
# draft sequence
# >SBJ
# reference sequence
# # the length of the sequences IDs in input.fa, has to be exactly 3 characters
# use the tool as follows:
# python DraftPolisher_circular.py -q query.fa -s subject.fa -f sequences.fa
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

print("DraftPolisher_lin_cov v1.0 by Rosario Nicola Brancaccio")
print("Start")
# run muscle
parser = argparse.ArgumentParser(description='polish draft circular genomes')
parser.add_argument("--input1", "-q", help="query sequence file", type=str)
parser.add_argument("--input2", "-s", help="subject sequence file", type=str)
parser.add_argument("--input3", "-f", help="SPAdes contigs file (or any sequences file in FASTA format)", type=str)
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
        if i == 0:
            if ja[i] == "_" or ja[i + 1] == "_" or ja[i + 2] == "_" or ja[i + 3] == "_" or ja[i + 4] == "_" or \
                    ja[i + 5] == "_" or ja[i + 6] == "_" or ja[i + 7] == "_" or ja[i + 8] == "_" or ja[i + 9] == "_" or \
                    ja[i + 10] == "_" or ja[i + 11] == "_" or ja[i + 12] == "_" or ja[i + 13] == "_" or ja[
                i + 14] == "_" or ja[i + 15] == "_" or ja[i + 16] == "_":
                print(ja[i], ja[i + 1], ja[i + 2], ja[i + 3], ja[i + 4], ja[i + 5], ja[i + 6],
                      ja[i + 7],
                      ja[i + 8], ja[i + 9], ja[i + 10], ja[i + 11], ja[i + 12], ja[i + 13], ja[i + 14],
                      ja[i + 15],
                      ja[i + 16], ja[i + 17], ja[i + 18], ja[i + 19], ja[i + 20], ja[i + 21], ja[i + 22],
                      ja[i + 23],
                      ja[i + 24], ja[i + 25], ja[i + 26], ja[i + 27], ja[i + 28], ja[i + 29], ja[i + 30], ja[i + 31],
                      ja[i + 32],
                      ja[i + 33], ja[i + 34], ja[i + 35], ja[i + 36], ja[i + 37], ja[i + 38], ja[i + 39], ja[i + 40],
                      ja[i + 41],
                      ja[i + 42], ja[i + 43], ja[i + 44], ja[i + 45], ja[i + 46], ja[i + 47], ja[i + 48], ja[i + 49],
                      ja[i + 50], ja[i + 51], ja[i + 52], ja[i + 53], ja[i + 54], ja[i + 55], ja[i + 56], ja[i + 57],
                      ja[i + 58], ja[i + 59], ja[i + 60], ja[i + 61], ja[i + 62], ja[i + 63], ja[i + 64], ja[i + 65],
                      ja[i + 66], ja[i + 67], ja[i + 68], file=open("output4.txt", "a+"), sep="")
                strix = open("output4.txt")
                striy = strix.readlines()
                string = striy[-1]
                res = [sub.replace('_', '') for sub in string]
                x = "".join(res)
                print(x[:17], file=open("oligo.txt", "a+"))
            elif ja[i] != "_" or ja[i + 1] != "_" or ja[i + 2] != "_" or ja[i + 3] != "_" or ja[i + 4] != "_" or \
                    ja[i + 5] != "_" or ja[i + 6] != "_" or ja[i + 7] != "_" or ja[i + 9] != "_" or ja[i + 10] != "_" or \
                    ja[i + 11] != "_" or ja[i + 12] != "_" or ja[i + 13] != "_" or ja[i + 14] != "_" or ja[
                i + 15] != "_" or ja[
                i + 16] != "_":
                print(ja[i], ja[i + 1], ja[i + 2], ja[i + 3], ja[i + 4], ja[i + 5], ja[i + 6], ja[i + 7], ja[i + 8],
                      ja[i + 9], ja[i + 10], ja[i + 11], ja[i + 12], ja[i + 13], ja[i + 14], ja[i + 15], ja[i + 16],
                      file=open("oligo.txt", "a+"), sep="")
        if i == 1:
            if ja[i - 1] == "_" or ja[i] == "_" or ja[i + 1] == "_" or ja[i + 2] == "_" or ja[i + 3] == "_" or ja[
                i + 4] == "_" or \
                    ja[i + 5] == "_" or ja[i + 6] == "_" or ja[i + 7] == "_" or ja[i + 8] == "_" or ja[
                i + 9] == "_" or ja[i + 10] == "_" or ja[i + 11] == "_" or ja[i + 12] == "_" or ja[i + 13] == "_" \
                    or ja[i + 14] == "_" or ja[i + 15] == "_":
                print(ja[i - 1], ja[i], ja[i + 1], ja[i + 2], ja[i + 3], ja[i + 4], ja[i + 5], ja[i + 6],
                      ja[i + 7],
                      ja[i + 8], ja[i + 9], ja[i + 10], ja[i + 11], ja[i + 12], ja[i + 13], ja[i + 14],
                      ja[i + 15],
                      ja[i + 16], ja[i + 17], ja[i + 18], ja[i + 19], ja[i + 20], ja[i + 21], ja[i + 22],
                      ja[i + 23],
                      ja[i + 24], ja[i + 25], ja[i + 26], ja[i + 27], ja[i + 28], ja[i + 29], ja[i + 30],
                      ja[i + 31],
                      ja[i + 32],
                      ja[i + 33], ja[i + 34], ja[i + 35], ja[i + 36], ja[i + 37], ja[i + 38], ja[i + 39],
                      ja[i + 40], ja[i + 41],
                      ja[i + 42], ja[i + 43], ja[i + 44], ja[i + 45], ja[i + 46], ja[i + 47], ja[i + 48],
                      ja[i + 49],
                      ja[i + 50], ja[i + 51], ja[i + 52], ja[i + 53], ja[i + 54], ja[i + 55], ja[i + 56],
                      ja[i + 57],
                      ja[i + 58], ja[i + 59], ja[i + 60], ja[i + 61], ja[i + 62], ja[i + 63], ja[i + 64],
                      ja[i + 65],
                      ja[i + 66], ja[i + 67], file=open("output4.txt", "a+"), sep="")
                strix = open("output4.txt")
                striy = strix.readlines()
                string = striy[-1]
                res = [sub.replace('_', '') for sub in string]
                x = "".join(res)
                print(x[:17], file=open("oligo.txt", "a+"))
            else:
                print(ja[i - 1], ja[i], ja[i + 1], ja[i + 2], ja[i + 3], ja[i + 4], ja[i + 5], ja[i + 6], ja[i + 7],
                      ja[i + 8], ja[i + 9], ja[i + 10], ja[i + 11], ja[i + 12], ja[i + 13], ja[i + 14], ja[i + 15],
                      file=open("oligo.txt", "a+"), sep="")
        if i == 2:
            if ja[i - 2] == "_" or ja[i - 1] == "_" or ja[i] == "_" or ja[i + 1] == "_" or ja[i + 2] == "_" or ja[
                i + 3] == "_" or ja[i + 4] == "_" or \
                    ja[i + 5] == "_" or ja[i + 6] == "_" or ja[i + 7] == "_" or ja[i + 8] == "_" or ja[
                i + 9] == "_" or ja[i + 10] == "_" or ja[i + 11] == "_" or ja[i + 12] == "_" or ja[i + 13] == "_" \
                    or ja[i + 14] == "_":
                print(ja[i - 2], ja[i - 1], ja[i], ja[i + 1], ja[i + 2], ja[i + 3], ja[i + 4], ja[i + 5], ja[i + 6],
                      ja[i + 7],
                      ja[i + 8], ja[i + 9], ja[i + 10], ja[i + 11], ja[i + 12], ja[i + 13], ja[i + 14],
                      ja[i + 15],
                      ja[i + 16], ja[i + 17], ja[i + 18], ja[i + 19], ja[i + 20], ja[i + 21], ja[i + 22],
                      ja[i + 23],
                      ja[i + 24], ja[i + 25], ja[i + 26], ja[i + 27], ja[i + 28], ja[i + 29], ja[i + 30],
                      ja[i + 31],
                      ja[i + 32],
                      ja[i + 33], ja[i + 34], ja[i + 35], ja[i + 36], ja[i + 37], ja[i + 38], ja[i + 39],
                      ja[i + 40], ja[i + 41],
                      ja[i + 42], ja[i + 43], ja[i + 44], ja[i + 45], ja[i + 46], ja[i + 47], ja[i + 48],
                      ja[i + 49],
                      ja[i + 50], ja[i + 51], ja[i + 52], ja[i + 53], ja[i + 54], ja[i + 55], ja[i + 56],
                      ja[i + 57],
                      ja[i + 58], ja[i + 59], ja[i + 60], ja[i + 61], ja[i + 62], ja[i + 63], ja[i + 64],
                      ja[i + 65],
                      ja[i + 66], file=open("output4.txt", "a+"), sep="")
                strix = open("output4.txt")
                striy = strix.readlines()
                string = striy[-1]
                res = [sub.replace('_', '') for sub in string]
                x = "".join(res)
                print(x[:17], file=open("oligo.txt", "a+"))
            else:
                print(ja[i - 2], ja[i - 1], ja[i], ja[i + 1], ja[i + 2], ja[i + 3], ja[i + 4], ja[i + 5], ja[i + 6],
                      ja[i + 7], ja[i + 8], ja[i + 9], ja[i + 10], ja[i + 11], ja[i + 12], ja[i + 13],
                      ja[i + 14], file=open("oligo.txt", "a+"), sep="")
        if i == 3:
            if ja[i - 3] == "_" or ja[i - 2] == "_" or ja[i - 1] == "_" or ja[i] == "_" or ja[i + 1] == "_" or ja[
                i + 2] == "_" or ja[i + 3] == "_" or ja[i + 4] == "_" or \
                    ja[i + 5] == "_" or ja[i + 6] == "_" or ja[i + 7] == "_" or ja[i + 8] == "_" or ja[i + 9] == "_" \
                    or ja[i + 10] == "_" or ja[i + 11] == "_" or ja[i + 12] == "_" or ja[i + 13] == "_":
                print(ja[i - 3], ja[i - 2], ja[i - 1], ja[i], ja[i + 1], ja[i + 2], ja[i + 3], ja[i + 4], ja[i + 5],
                      ja[i + 6],
                      ja[i + 7],
                      ja[i + 8], ja[i + 9], ja[i + 10], ja[i + 11], ja[i + 12], ja[i + 13], ja[i + 14],
                      ja[i + 15],
                      ja[i + 16], ja[i + 17], ja[i + 18], ja[i + 19], ja[i + 20], ja[i + 21], ja[i + 22],
                      ja[i + 23],
                      ja[i + 24], ja[i + 25], ja[i + 26], ja[i + 27], ja[i + 28], ja[i + 29], ja[i + 30],
                      ja[i + 31],
                      ja[i + 32],
                      ja[i + 33], ja[i + 34], ja[i + 35], ja[i + 36], ja[i + 37], ja[i + 38], ja[i + 39],
                      ja[i + 40], ja[i + 41],
                      ja[i + 42], ja[i + 43], ja[i + 44], ja[i + 45], ja[i + 46], ja[i + 47], ja[i + 48],
                      ja[i + 49],
                      ja[i + 50], ja[i + 51], ja[i + 52], ja[i + 53], ja[i + 54], ja[i + 55], ja[i + 56],
                      ja[i + 57],
                      ja[i + 58], ja[i + 59], ja[i + 60], ja[i + 61], ja[i + 62], ja[i + 63], ja[i + 64],
                      ja[i + 65], file=open("output4.txt", "a+"), sep="")
                strix = open("output4.txt")
                striy = strix.readlines()
                string = striy[-1]
                res = [sub.replace('_', '') for sub in string]
                x = "".join(res)
                print(x[:17], file=open("oligo.txt", "a+"))
            else:
                print(ja[i - 3], ja[i - 2], ja[i - 1], ja[i], ja[i + 1], ja[i + 2], ja[i + 3], ja[i + 4], ja[i + 5],
                      ja[i + 6],
                      ja[i + 7], ja[i + 8], ja[i + 9], ja[i + 10], ja[i + 11], ja[i + 12], ja[i + 13],
                      file=open("oligo.txt", "a+"), sep="")
        if i == 4:
            if ja[i - 4] == "_" or ja[i - 3] == "_" or ja[i - 2] == "_" or ja[i - 1] == "_" or ja[i] == "_" or ja[
                i + 1] == "_" or ja[i + 2] == "_" or ja[i + 3] == "_" or ja[i + 4] == "_" or \
                    ja[i + 5] == "_" or ja[i + 6] == "_" or ja[i + 7] == "_" or ja[i + 8] == "_" or ja[i + 9] == "_" \
                    or ja[i + 10] == "_" or ja[i + 11] == "_" or ja[i + 12] == "_":
                print(ja[i - 4], ja[i - 3], ja[i - 2], ja[i - 1], ja[i], ja[i + 1], ja[i + 2], ja[i + 3], ja[i + 4],
                      ja[i + 5], ja[i + 6],
                      ja[i + 7],
                      ja[i + 8], ja[i + 9], ja[i + 10], ja[i + 11], ja[i + 12], ja[i + 13], ja[i + 14],
                      ja[i + 15],
                      ja[i + 16], ja[i + 17], ja[i + 18], ja[i + 19], ja[i + 20], ja[i + 21], ja[i + 22],
                      ja[i + 23],
                      ja[i + 24], ja[i + 25], ja[i + 26], ja[i + 27], ja[i + 28], ja[i + 29], ja[i + 30],
                      ja[i + 31],
                      ja[i + 32],
                      ja[i + 33], ja[i + 34], ja[i + 35], ja[i + 36], ja[i + 37], ja[i + 38], ja[i + 39],
                      ja[i + 40], ja[i + 41],
                      ja[i + 42], ja[i + 43], ja[i + 44], ja[i + 45], ja[i + 46], ja[i + 47], ja[i + 48],
                      ja[i + 49],
                      ja[i + 50], ja[i + 51], ja[i + 52], ja[i + 53], ja[i + 54], ja[i + 55], ja[i + 56],
                      ja[i + 57],
                      ja[i + 58], ja[i + 59], ja[i + 60], ja[i + 61], ja[i + 62], ja[i + 63], ja[i + 64],
                      file=open("output4.txt", "a+"),
                      sep="")
                strix = open("output4.txt")
                striy = strix.readlines()
                string = striy[-1]
                res = [sub.replace('_', '') for sub in string]
                x = "".join(res)
                print(x[:17], file=open("oligo.txt", "a+"))
            else:
                print(ja[i - 4], ja[i - 3], ja[i - 2], ja[i - 1], ja[i], ja[i + 1], ja[i + 2], ja[i + 3], ja[i + 4],
                      ja[i + 5], ja[i + 6],
                      ja[i + 7], ja[i + 8], ja[i + 9], ja[i + 10], ja[i + 11], ja[i + 12], file=open("oligo.txt", "a+"),
                      sep="")
        if i == 5:
            if ja[i - 5] == "_" or ja[i - 4] == "_" or ja[i - 3] == "_" or ja[i - 2] == "_" or ja[i - 1] == "_" or ja[
                i] == "_" or ja[i + 1] == "_" or ja[i + 2] == "_" or ja[i + 3] == "_" or ja[i + 4] == "_" or \
                    ja[i + 5] == "_" or ja[i + 6] == "_" or ja[i + 7] == "_" or ja[i + 8] == "_" or ja[i + 9] == "_" \
                    or ja[i + 10] == "_" or ja[i + 11] == "_":
                print(ja[i - 5], ja[i - 4], ja[i - 3], ja[i - 2], ja[i - 1], ja[i], ja[i + 1], ja[i + 2], ja[i + 3],
                      ja[i + 4], ja[i + 5], ja[i + 6],
                      ja[i + 7], ja[i + 8], ja[i + 9], ja[i + 10], ja[i + 11], ja[i + 12], ja[i + 13], ja[i + 14],
                      ja[i + 15], ja[i + 16], ja[i + 17], ja[i + 18], ja[i + 19], ja[i + 20], ja[i + 21], ja[i + 22],
                      ja[i + 23], ja[i + 24], ja[i + 25], ja[i + 26], ja[i + 27], ja[i + 28], ja[i + 29], ja[i + 30],
                      ja[i + 31], ja[i + 32], ja[i + 33], ja[i + 34], ja[i + 35], ja[i + 36], ja[i + 37],
                      ja[i + 38], ja[i + 39],
                      ja[i + 40], ja[i + 41], ja[i + 42], ja[i + 43], ja[i + 44], ja[i + 45], ja[i + 46], ja[i + 47],
                      ja[i + 48],
                      ja[i + 49], ja[i + 50], ja[i + 51], ja[i + 52], ja[i + 53], ja[i + 54], ja[i + 55], ja[i + 56],
                      ja[i + 57], ja[i + 58], ja[i + 59], ja[i + 60], ja[i + 61], ja[i + 62], ja[i + 63],
                      file=open("output4.txt", "a+"), sep="")
                strix = open("output4.txt")
                striy = strix.readlines()
                string = striy[-1]
                res = [sub.replace('_', '') for sub in string]
                x = "".join(res)
                print(x[:17], file=open("oligo.txt", "a+"))
            else:
                print(ja[i - 5], ja[i - 4], ja[i - 3], ja[i - 2], ja[i - 1], ja[i], ja[i + 1], ja[i + 2], ja[i + 3],
                      ja[i + 4],
                      ja[i + 5], ja[i + 6],
                      ja[i + 7], ja[i + 8], ja[i + 9], ja[i + 10], ja[i + 11], file=open("oligo.txt", "a+"),
                      sep="")
        if i == 6:
            if ja[i - 6] == "_" or ja[i - 5] == "_" or ja[i - 4] == "_" or ja[i - 3] == "_" or ja[i - 2] == "_" or \
                    ja[i - 1] == "_" or ja[i] == "_" or ja[i + 1] == "_" or ja[i + 2] == "_" or ja[i + 3] == "_" or \
                    ja[i + 4] == "_" or ja[i + 5] == "_" or ja[i + 6] == "_" or ja[i + 7] == "_" or ja[i + 8] == "_" or \
                    ja[i + 9] == "_" or ja[i + 10] == "_":
                print(ja[i - 6], ja[i - 5], ja[i - 4], ja[i - 3], ja[i - 2], ja[i - 1], ja[i], ja[i + 1], ja[i + 2],
                      ja[i + 3],
                      ja[i + 4], ja[i + 5], ja[i + 6],
                      ja[i + 7], ja[i + 8], ja[i + 9], ja[i + 10], ja[i + 11], ja[i + 12], ja[i + 13], ja[i + 14],
                      ja[i + 15], ja[i + 16], ja[i + 17], ja[i + 18], ja[i + 19], ja[i + 20], ja[i + 21], ja[i + 22],
                      ja[i + 23], ja[i + 24], ja[i + 25], ja[i + 26], ja[i + 27], ja[i + 28], ja[i + 29], ja[i + 30],
                      ja[i + 31], ja[i + 32], ja[i + 33], ja[i + 34], ja[i + 35], ja[i + 36], ja[i + 37],
                      ja[i + 38], ja[i + 39],
                      ja[i + 40], ja[i + 41], ja[i + 42], ja[i + 43], ja[i + 44], ja[i + 45], ja[i + 46], ja[i + 47],
                      ja[i + 48],
                      ja[i + 49], ja[i + 50], ja[i + 51], ja[i + 52], ja[i + 53], ja[i + 54], ja[i + 55], ja[i + 56],
                      ja[i + 57], ja[i + 58], ja[i + 59], ja[i + 60], ja[i + 61], ja[i + 62],
                      file=open("output4.txt", "a+"), sep="")
                strix = open("output4.txt")
                striy = strix.readlines()
                string = striy[-1]
                res = [sub.replace('_', '') for sub in string]
                x = "".join(res)
                print(x[:17], file=open("oligo.txt", "a+"))
            else:
                print(ja[i - 6], ja[i - 5], ja[i - 4], ja[i - 3], ja[i - 2], ja[i - 1], ja[i], ja[i + 1], ja[i + 2],
                      ja[i + 3], ja[i + 4], ja[i + 5], ja[i + 6],
                      ja[i + 7], ja[i + 8], ja[i + 9], ja[i + 10], file=open("oligo.txt", "a+"), sep="")
        if i == 7:
            if ja[i - 7] == "_" or ja[i - 6] == "_" or ja[i - 5] == "_" or ja[i - 4] == "_" or ja[i - 3] == "_" or ja[
                i - 2] == "_" or \
                    ja[i - 1] == "_" or ja[i] == "_" or ja[i + 1] == "_" or ja[i + 2] == "_" or ja[i + 3] == "_" or \
                    ja[i + 4] == "_" or ja[i + 5] == "_" or ja[i + 6] == "_" or ja[i + 7] == "_" or ja[i + 8] == "_" or \
                    ja[i + 9] == "_":
                print(ja[i - 7], ja[i - 6], ja[i - 5], ja[i - 4], ja[i - 3], ja[i - 2], ja[i - 1], ja[i], ja[i + 1],
                      ja[i + 2], ja[i + 3],
                      ja[i + 4], ja[i + 5], ja[i + 6],
                      ja[i + 7], ja[i + 8], ja[i + 9], ja[i + 10], ja[i + 11], ja[i + 12], ja[i + 13], ja[i + 14],
                      ja[i + 15], ja[i + 16], ja[i + 17], ja[i + 18], ja[i + 19], ja[i + 20], ja[i + 21], ja[i + 22],
                      ja[i + 23], ja[i + 24], ja[i + 25], ja[i + 26], ja[i + 27], ja[i + 28], ja[i + 29], ja[i + 30],
                      ja[i + 31], ja[i + 32], ja[i + 33], ja[i + 34], ja[i + 35], ja[i + 36], ja[i + 37],
                      ja[i + 38], ja[i + 39],
                      ja[i + 40], ja[i + 41], ja[i + 42], ja[i + 43], ja[i + 44], ja[i + 45], ja[i + 46], ja[i + 47],
                      ja[i + 48],
                      ja[i + 49], ja[i + 50], ja[i + 51], ja[i + 52], ja[i + 53], ja[i + 54], ja[i + 55], ja[i + 56],
                      ja[i + 57], ja[i + 58], ja[i + 59], ja[i + 60], ja[i + 61],
                      file=open("output4.txt", "a+"), sep="")
                strix = open("output4.txt")
                striy = strix.readlines()
                string = striy[-1]
                res = [sub.replace('_', '') for sub in string]
                x = "".join(res)
                print(x[:17], file=open("oligo.txt", "a+"))
            else:
                print(ja[i - 7], ja[i - 6], ja[i - 5], ja[i - 4], ja[i - 3], ja[i - 2], ja[i - 1], ja[i], ja[i + 1],
                      ja[i + 2], ja[i + 3], ja[i + 4], ja[i + 5], ja[i + 6],
                      ja[i + 7], ja[i + 8], ja[i + 9], file=open("oligo.txt", "a+"), sep="")
        if i == 8:
            if ja[i - 8] == "_" or ja[i - 7] == "_" or ja[i - 6] == "_" or ja[i - 5] == "_" or ja[i - 4] == "_" or ja[
                i - 3] == "_" or ja[i - 2] == "_" or \
                    ja[i - 1] == "_" or ja[i] == "_" or ja[i + 1] == "_" or ja[i + 2] == "_" or ja[i + 3] == "_" or \
                    ja[i + 4] == "_" or ja[i + 5] == "_" or ja[i + 6] == "_" or ja[i + 7] == "_" or ja[i + 8] == "_":
                print(ja[i - 8], ja[i - 7], ja[i - 6], ja[i - 5], ja[i - 4], ja[i - 3], ja[i - 2], ja[i - 1], ja[i],
                      ja[i + 1], ja[i + 2], ja[i + 3],
                      ja[i + 4], ja[i + 5], ja[i + 6],
                      ja[i + 7], ja[i + 8], ja[i + 9], ja[i + 10], ja[i + 11], ja[i + 12], ja[i + 13], ja[i + 14],
                      ja[i + 15], ja[i + 16], ja[i + 17], ja[i + 18], ja[i + 19], ja[i + 20], ja[i + 21], ja[i + 22],
                      ja[i + 23], ja[i + 24], ja[i + 25], ja[i + 26], ja[i + 27], ja[i + 28], ja[i + 29], ja[i + 30],
                      ja[i + 31], ja[i + 32], ja[i + 33], ja[i + 34], ja[i + 35], ja[i + 36], ja[i + 37],
                      ja[i + 38], ja[i + 39],
                      ja[i + 40], ja[i + 41], ja[i + 42], ja[i + 43], ja[i + 44], ja[i + 45], ja[i + 46], ja[i + 47],
                      ja[i + 48],
                      ja[i + 49], ja[i + 50], ja[i + 51], ja[i + 52], ja[i + 53], ja[i + 54], ja[i + 55], ja[i + 56],
                      ja[i + 57], ja[i + 58], ja[i + 59], ja[i + 60],
                      file=open("output4.txt", "a+"), sep="")
                strix = open("output4.txt")
                striy = strix.readlines()
                string = striy[-1]
                res = [sub.replace('_', '') for sub in string]
                x = "".join(res)
                print(x[:17], file=open("oligo.txt", "a+"))
            else:
                print(ja[i - 8], ja[i - 7], ja[i - 6], ja[i - 5], ja[i - 4], ja[i - 3], ja[i - 2], ja[i - 1], ja[i],
                      ja[i + 1],
                      ja[i + 2], ja[i + 3], ja[i + 4], ja[i + 5], ja[i + 6],
                      ja[i + 7], ja[i + 8], file=open("oligo.txt", "a+"), sep="")
        if i == 9:
            if ja[i - 9] == "_" or ja[i - 8] == "_" or ja[i - 7] == "_" or ja[i - 6] == "_" or ja[i - 5] == "_" or ja[
                i - 4] == "_" or ja[i - 3] == "_" or ja[i - 2] == "_" or \
                    ja[i - 1] == "_" or ja[i] == "_" or ja[i + 1] == "_" or ja[i + 2] == "_" or ja[i + 3] == "_" or \
                    ja[i + 4] == "_" or ja[i + 5] == "_" or ja[i + 6] == "_" or ja[i + 7] == "_":
                print(ja[i - 9], ja[i - 8], ja[i - 7], ja[i - 6], ja[i - 5], ja[i - 4], ja[i - 3], ja[i - 2], ja[i - 1],
                      ja[i], ja[i + 1], ja[i + 2], ja[i + 3],
                      ja[i + 4], ja[i + 5], ja[i + 6],
                      ja[i + 7], ja[i + 8], ja[i + 9], ja[i + 10], ja[i + 11], ja[i + 12], ja[i + 13], ja[i + 14],
                      ja[i + 15], ja[i + 16], ja[i + 17], ja[i + 18], ja[i + 19], ja[i + 20], ja[i + 21], ja[i + 22],
                      ja[i + 23], ja[i + 24], ja[i + 25], ja[i + 26], ja[i + 27], ja[i + 28], ja[i + 29], ja[i + 30],
                      ja[i + 31], ja[i + 32], ja[i + 33], ja[i + 34], ja[i + 35], ja[i + 36], ja[i + 37],
                      ja[i + 38], ja[i + 39],
                      ja[i + 40], ja[i + 41], ja[i + 42], ja[i + 43], ja[i + 44], ja[i + 45], ja[i + 46], ja[i + 47],
                      ja[i + 48],
                      ja[i + 49], ja[i + 50], ja[i + 51], ja[i + 52], ja[i + 53], ja[i + 54], ja[i + 55], ja[i + 56],
                      ja[i + 57], ja[i + 58], ja[i + 59],
                      file=open("output4.txt", "a+"), sep="")
                strix = open("output4.txt")
                striy = strix.readlines()
                string = striy[-1]
                res = [sub.replace('_', '') for sub in string]
                x = "".join(res)
                print(x[:17], file=open("oligo.txt", "a+"))
            else:
                print(ja[i - 9], ja[i - 8], ja[i - 7], ja[i - 6], ja[i - 5], ja[i - 4], ja[i - 3], ja[i - 2], ja[i - 1],
                      ja[i], ja[i + 1],
                      ja[i + 2], ja[i + 3], ja[i + 4], ja[i + 5], ja[i + 6],
                      ja[i + 7], file=open("oligo.txt", "a+"), sep="")
        if i == 10:
            if ja[i - 10] == "_" or ja[i - 9] == "_" or ja[i - 8] == "_" or ja[i - 7] == "_" or ja[i - 6] == "_" or ja[
                i - 5] == "_" or ja[i - 4] == "_" or ja[i - 3] == "_" or ja[i - 2] == "_" or \
                    ja[i - 1] == "_" or ja[i] == "_" or ja[i + 1] == "_" or ja[i + 2] == "_" or ja[i + 3] == "_" or \
                    ja[i + 4] == "_" or ja[i + 5] == "_" or ja[i + 6] == "_":
                print(ja[i - 10], ja[i - 9], ja[i - 8], ja[i - 7], ja[i - 6], ja[i - 5], ja[i - 4], ja[i - 3],
                      ja[i - 2], ja[i - 1], ja[i], ja[i + 1], ja[i + 2], ja[i + 3],
                      ja[i + 4], ja[i + 5], ja[i + 6],
                      ja[i + 7], ja[i + 8], ja[i + 9], ja[i + 10], ja[i + 11], ja[i + 12], ja[i + 13], ja[i + 14],
                      ja[i + 15], ja[i + 16], ja[i + 17], ja[i + 18], ja[i + 19], ja[i + 20], ja[i + 21], ja[i + 22],
                      ja[i + 23], ja[i + 24], ja[i + 25], ja[i + 26], ja[i + 27], ja[i + 28], ja[i + 29], ja[i + 30],
                      ja[i + 31], ja[i + 32], ja[i + 33], ja[i + 34], ja[i + 35], ja[i + 36], ja[i + 37],
                      ja[i + 38], ja[i + 39],
                      ja[i + 40], ja[i + 41], ja[i + 42], ja[i + 43], ja[i + 44], ja[i + 45], ja[i + 46], ja[i + 47],
                      ja[i + 48],
                      ja[i + 49], ja[i + 50], ja[i + 51], ja[i + 52], ja[i + 53], ja[i + 54], ja[i + 55], ja[i + 56],
                      ja[i + 57], ja[i + 58],
                      file=open("output4.txt", "a+"), sep="")
                strix = open("output4.txt")
                striy = strix.readlines()
                string = striy[-1]
                res = [sub.replace('_', '') for sub in string]
                x = "".join(res)
                print(x[:17], file=open("oligo.txt", "a+"))
            else:
                print(ja[i - 10], ja[i - 9], ja[i - 8], ja[i - 7], ja[i - 6], ja[i - 5], ja[i - 4], ja[i - 3],
                      ja[i - 2], ja[i - 1], ja[i], ja[i + 1], ja[i + 2], ja[i + 3], ja[i + 4], ja[i + 5], ja[i + 6],
                      file=open("oligo.txt", "a+"), sep="")
        if i == 11:
            if ja[i - 11] == "_" or ja[i - 10] == "_" or ja[i - 9] == "_" or ja[i - 8] == "_" or ja[i - 7] == "_" or ja[
                i - 6] == "_" or ja[i - 5] == "_" or ja[i - 4] == "_" or ja[i - 3] == "_" or ja[i - 2] == "_" or \
                    ja[i - 1] == "_" or ja[i] == "_" or ja[i + 1] == "_" or ja[i + 2] == "_" or ja[i + 3] == "_" or \
                    ja[i + 4] == "_" or ja[i + 5] == "_":
                print(ja[i - 11], ja[i - 10], ja[i - 9], ja[i - 8], ja[i - 7], ja[i - 6], ja[i - 5], ja[i - 4],
                      ja[i - 3], ja[i - 2], ja[i - 1], ja[i], ja[i + 1], ja[i + 2], ja[i + 3],
                      ja[i + 4], ja[i + 5], ja[i + 6],
                      ja[i + 7], ja[i + 8], ja[i + 9], ja[i + 10], ja[i + 11], ja[i + 12], ja[i + 13], ja[i + 14],
                      ja[i + 15], ja[i + 16], ja[i + 17], ja[i + 18], ja[i + 19], ja[i + 20], ja[i + 21], ja[i + 22],
                      ja[i + 23], ja[i + 24], ja[i + 25], ja[i + 26], ja[i + 27], ja[i + 28], ja[i + 29], ja[i + 30],
                      ja[i + 31], ja[i + 32], ja[i + 33], ja[i + 34], ja[i + 35], ja[i + 36], ja[i + 37],
                      ja[i + 38], ja[i + 39],
                      ja[i + 40], ja[i + 41], ja[i + 42], ja[i + 43], ja[i + 44], ja[i + 45], ja[i + 46], ja[i + 47],
                      ja[i + 48],
                      ja[i + 49], ja[i + 50], ja[i + 51], ja[i + 52], ja[i + 53], ja[i + 54], ja[i + 55], ja[i + 56],
                      ja[i + 57],
                      file=open("output4.txt", "a+"), sep="")
                strix = open("output4.txt")
                striy = strix.readlines()
                string = striy[-1]
                res = [sub.replace('_', '') for sub in string]
                x = "".join(res)
                print(x[:17], file=open("oligo.txt", "a+"))
            else:
                print(ja[i - 11], ja[i - 10], ja[i - 9], ja[i - 8], ja[i - 7], ja[i - 6], ja[i - 5], ja[i - 4],
                      ja[i - 3],
                      ja[i - 2], ja[i - 1], ja[i], ja[i + 1], ja[i + 2], ja[i + 3], ja[i + 4], ja[i + 5],
                      file=open("oligo.txt", "a+"), sep="")
        if i == 12:
            if ja[i - 12] == "_" or ja[i - 11] == "_" or ja[i - 10] == "_" or ja[i - 9] == "_" or ja[i - 8] == "_" or \
                    ja[i - 7] == "_" or ja[i - 6] == "_" or ja[i - 5] == "_" or ja[i - 4] == "_" or ja[i - 3] == "_" or \
                    ja[i - 2] == "_" or \
                    ja[i - 1] == "_" or ja[i] == "_" or ja[i + 1] == "_" or ja[i + 2] == "_" or ja[i + 3] == "_" or \
                    ja[i + 4] == "_":
                print(ja[i - 12], ja[i - 11], ja[i - 10], ja[i - 9], ja[i - 8], ja[i - 7], ja[i - 6], ja[i - 5],
                      ja[i - 4], ja[i - 3], ja[i - 2], ja[i - 1], ja[i], ja[i + 1], ja[i + 2], ja[i + 3],
                      ja[i + 4], ja[i + 5], ja[i + 6],
                      ja[i + 7], ja[i + 8], ja[i + 9], ja[i + 10], ja[i + 11], ja[i + 12], ja[i + 13], ja[i + 14],
                      ja[i + 15], ja[i + 16], ja[i + 17], ja[i + 18], ja[i + 19], ja[i + 20], ja[i + 21], ja[i + 22],
                      ja[i + 23], ja[i + 24], ja[i + 25], ja[i + 26], ja[i + 27], ja[i + 28], ja[i + 29], ja[i + 30],
                      ja[i + 31], ja[i + 32], ja[i + 33], ja[i + 34], ja[i + 35], ja[i + 36], ja[i + 37],
                      ja[i + 38], ja[i + 39],
                      ja[i + 40], ja[i + 41], ja[i + 42], ja[i + 43], ja[i + 44], ja[i + 45], ja[i + 46], ja[i + 47],
                      ja[i + 48],
                      ja[i + 49], ja[i + 50], ja[i + 51], ja[i + 52], ja[i + 53], ja[i + 54], ja[i + 55], ja[i + 56],
                      file=open("output4.txt", "a+"), sep="")
                strix = open("output4.txt")
                striy = strix.readlines()
                string = striy[-1]
                res = [sub.replace('_', '') for sub in string]
                x = "".join(res)
                print(x[:17], file=open("oligo.txt", "a+"))
            else:
                print(ja[i - 12], ja[i - 11], ja[i - 10], ja[i - 9], ja[i - 8], ja[i - 7], ja[i - 6], ja[i - 5],
                      ja[i - 4], ja[i - 3],
                      ja[i - 2], ja[i - 1], ja[i], ja[i + 1], ja[i + 2], ja[i + 3], ja[i + 4],
                      file=open("oligo.txt", "a+"), sep="")
        if i == 13:
            if ja[i - 13] == "_" or ja[i - 12] == "_" or ja[i - 11] == "_" or ja[i - 10] == "_" or ja[i - 9] == "_" or \
                    ja[i - 8] == "_" or ja[i - 7] == "_" or ja[i - 6] == "_" or ja[i - 5] == "_" or ja[i - 4] == "_" or \
                    ja[i - 3] == "_" or ja[i - 2] == "_" or \
                    ja[i - 1] == "_" or ja[i] == "_" or ja[i + 1] == "_" or ja[i + 2] == "_" or ja[i + 3] == "_":
                print(ja[i - 13], ja[i - 12], ja[i - 11], ja[i - 10], ja[i - 9], ja[i - 8], ja[i - 7], ja[i - 6],
                      ja[i - 5], ja[i - 4], ja[i - 3], ja[i - 2], ja[i - 1], ja[i], ja[i + 1], ja[i + 2], ja[i + 3],
                      ja[i + 4], ja[i + 5], ja[i + 6],
                      ja[i + 7], ja[i + 8], ja[i + 9], ja[i + 10], ja[i + 11], ja[i + 12], ja[i + 13], ja[i + 14],
                      ja[i + 15], ja[i + 16], ja[i + 17], ja[i + 18], ja[i + 19], ja[i + 20], ja[i + 21], ja[i + 22],
                      ja[i + 23], ja[i + 24], ja[i + 25], ja[i + 26], ja[i + 27], ja[i + 28], ja[i + 29], ja[i + 30],
                      ja[i + 31], ja[i + 32], ja[i + 33], ja[i + 34], ja[i + 35], ja[i + 36], ja[i + 37],
                      ja[i + 38], ja[i + 39],
                      ja[i + 40], ja[i + 41], ja[i + 42], ja[i + 43], ja[i + 44], ja[i + 45], ja[i + 46], ja[i + 47],
                      ja[i + 48],
                      ja[i + 49], ja[i + 50], ja[i + 51], ja[i + 52], ja[i + 53], ja[i + 54], ja[i + 55],
                      file=open("output4.txt", "a+"), sep="")
                strix = open("output4.txt")
                striy = strix.readlines()
                string = striy[-1]
                res = [sub.replace('_', '') for sub in string]
                x = "".join(res)
                print(x[:17], file=open("oligo.txt", "a+"))
            else:
                print(ja[i - 13], ja[i - 12], ja[i - 11], ja[i - 10], ja[i - 9], ja[i - 8], ja[i - 7], ja[i - 6],
                      ja[i - 5], ja[i - 4], ja[i - 3],
                      ja[i - 2], ja[i - 1], ja[i], ja[i + 1], ja[i + 2], ja[i + 3],
                      file=open("oligo.txt", "a+"), sep="")
        if i == 14:
            if ja[i - 14] == "_" or ja[i - 13] == "_" or ja[i - 12] == "_" or ja[i - 11] == "_" or ja[i - 10] == "_" or \
                    ja[i - 9] == "_" or ja[i - 8] == "_" or ja[i - 7] == "_" or ja[i - 6] == "_" or ja[i - 5] == "_" or \
                    ja[i - 4] == "_" or ja[i - 3] == "_" or ja[i - 2] == "_" or \
                    ja[i - 1] == "_" or ja[i] == "_" or ja[i + 1] == "_" or ja[i + 2] == "_":
                print(ja[i - 14], ja[i - 13], ja[i - 12], ja[i - 11], ja[i - 10], ja[i - 9], ja[i - 8], ja[i - 7],
                      ja[i - 6], ja[i - 5], ja[i - 4], ja[i - 3], ja[i - 2], ja[i - 1], ja[i], ja[i + 1], ja[i + 2],
                      ja[i + 3],
                      ja[i + 4], ja[i + 5], ja[i + 6],
                      ja[i + 7], ja[i + 8], ja[i + 9], ja[i + 10], ja[i + 11], ja[i + 12], ja[i + 13], ja[i + 14],
                      ja[i + 15], ja[i + 16], ja[i + 17], ja[i + 18], ja[i + 19], ja[i + 20], ja[i + 21], ja[i + 22],
                      ja[i + 23], ja[i + 24], ja[i + 25], ja[i + 26], ja[i + 27], ja[i + 28], ja[i + 29], ja[i + 30],
                      ja[i + 31], ja[i + 32], ja[i + 33], ja[i + 34], ja[i + 35], ja[i + 36], ja[i + 37],
                      ja[i + 38], ja[i + 39],
                      ja[i + 40], ja[i + 41], ja[i + 42], ja[i + 43], ja[i + 44], ja[i + 45], ja[i + 46], ja[i + 47],
                      ja[i + 48],
                      ja[i + 49], ja[i + 50], ja[i + 51], ja[i + 52], ja[i + 53], ja[i + 54],
                      file=open("output4.txt", "a+"), sep="")
                strix = open("output4.txt")
                striy = strix.readlines()
                string = striy[-1]
                res = [sub.replace('_', '') for sub in string]
                x = "".join(res)
                print(x[:17], file=open("oligo.txt", "a+"))
            else:
                print(ja[i - 14], ja[i - 13], ja[i - 12], ja[i - 11], ja[i - 10], ja[i - 9], ja[i - 8], ja[i - 7],
                      ja[i - 6], ja[i - 5], ja[i - 4], ja[i - 3],
                      ja[i - 2], ja[i - 1], ja[i], ja[i + 1], ja[i + 2],
                      file=open("oligo.txt", "a+"), sep="")
        if i == 15:
            if ja[i - 15] == "_" or ja[i - 14] == "_" or ja[i - 13] == "_" or ja[i - 12] == "_" or ja[i - 11] == "_" or \
                    ja[i - 10] == "_" or ja[i - 9] == "_" or ja[i - 8] == "_" or ja[i - 7] == "_" or ja[i - 6] == "_" or \
                    ja[i - 5] == "_" or ja[i - 4] == "_" or ja[i - 3] == "_" or ja[i - 2] == "_" or \
                    ja[i - 1] == "_" or ja[i] == "_" or ja[i + 1] == "_":
                print(ja[i - 15], ja[i - 14], ja[i - 13], ja[i - 12], ja[i - 11], ja[i - 10], ja[i - 9], ja[i - 8],
                      ja[i - 7], ja[i - 6], ja[i - 5], ja[i - 4], ja[i - 3], ja[i - 2], ja[i - 1], ja[i], ja[i + 1],
                      ja[i + 2], ja[i + 3],
                      ja[i + 4], ja[i + 5], ja[i + 6],
                      ja[i + 7], ja[i + 8], ja[i + 9], ja[i + 10], ja[i + 11], ja[i + 12], ja[i + 13], ja[i + 14],
                      ja[i + 15], ja[i + 16], ja[i + 17], ja[i + 18], ja[i + 19], ja[i + 20], ja[i + 21], ja[i + 22],
                      ja[i + 23], ja[i + 24], ja[i + 25], ja[i + 26], ja[i + 27], ja[i + 28], ja[i + 29], ja[i + 30],
                      ja[i + 31], ja[i + 32], ja[i + 33], ja[i + 34], ja[i + 35], ja[i + 36], ja[i + 37],
                      ja[i + 38], ja[i + 39],
                      ja[i + 40], ja[i + 41], ja[i + 42], ja[i + 43], ja[i + 44], ja[i + 45], ja[i + 46], ja[i + 47],
                      ja[i + 48],
                      ja[i + 49], ja[i + 50], ja[i + 51], ja[i + 52], ja[i + 53], file=open("output4.txt", "a+"),
                      sep="")
                strix = open("output4.txt")
                striy = strix.readlines()
                string = striy[-1]
                res = [sub.replace('_', '') for sub in string]
                x = "".join(res)
                print(x[:17], file=open("oligo.txt", "a+"))
            else:
                print(ja[i - 15], ja[i - 14], ja[i - 13], ja[i - 12], ja[i - 11], ja[i - 10], ja[i - 9], ja[i - 8],
                      ja[i - 7], ja[i - 6], ja[i - 5], ja[i - 4], ja[i - 3],
                      ja[i - 2], ja[i - 1], ja[i], ja[i + 1],
                      file=open("oligo.txt", "a+"), sep="")
        if i == 16:
            if ja[i - 16] == "_" or ja[i - 15] == "_" or ja[i - 14] == "_" or ja[i - 13] == "_" or ja[i - 12] == "_" or \
                    ja[i - 11] == "_" or ja[i - 10] == "_" or ja[i - 9] == "_" or ja[i - 8] == "_" or \
                    ja[i - 7] == "_" or ja[i - 6] == "_" or ja[i - 5] == "_" or ja[i - 4] == "_" or ja[i - 3] == "_" or \
                    ja[i - 2] == "_" or \
                    ja[i - 1] == "_" or ja[i] == "_":
                print(ja[i - 16], ja[i - 15], ja[i - 14], ja[i - 13], ja[i - 12], ja[i - 11], ja[i - 10], ja[i - 9],
                      ja[i - 8], ja[i - 7], ja[i - 6], ja[i - 5], ja[i - 4], ja[i - 3], ja[i - 2], ja[i - 1], ja[i],
                      ja[i + 1], ja[i + 2], ja[i + 3],
                      ja[i + 4], ja[i + 5], ja[i + 6],
                      ja[i + 7], ja[i + 8], ja[i + 9], ja[i + 10], ja[i + 11], ja[i + 12], ja[i + 13], ja[i + 14],
                      ja[i + 15], ja[i + 16], ja[i + 17], ja[i + 18], ja[i + 19], ja[i + 20], ja[i + 21], ja[i + 22],
                      ja[i + 23], ja[i + 24], ja[i + 25], ja[i + 26], ja[i + 27], ja[i + 28], ja[i + 29], ja[i + 30],
                      ja[i + 31], ja[i + 32], ja[i + 33], ja[i + 34], ja[i + 35], ja[i + 36], ja[i + 37],
                      ja[i + 38], ja[i + 39],
                      ja[i + 40], ja[i + 41], ja[i + 42], ja[i + 43], ja[i + 44], ja[i + 45], ja[i + 46], ja[i + 47],
                      ja[i + 48],
                      ja[i + 49], ja[i + 50], ja[i + 51], ja[i + 52], file=open("output4.txt", "a+"), sep="")
                strix = open("output4.txt")
                striy = strix.readlines()
                string = striy[-1]
                res = [sub.replace('_', '') for sub in string]
                x = "".join(res)
                print(x[:17], file=open("oligo.txt", "a+"))
            else:
                print(ja[i - 16], ja[i - 15], ja[i - 14], ja[i - 13], ja[i - 12], ja[i - 11], ja[i - 10], ja[i - 9],
                      ja[i - 8], ja[i - 7], ja[i - 6], ja[i - 5], ja[i - 4], ja[i - 3],
                      ja[i - 2], ja[i - 1], ja[i], file=open("oligo.txt", "a+"), sep="")
        if i == 17:
            if ja[i] == "_" or ja[i + 1] == "_" or ja[i + 2] == "_" or ja[i + 3] == "_" or ja[i + 4] == "_" or ja[i + 5] \
                    == "_" or ja[i + 6] == "_" or ja[i + 7] == "_" or ja[i + 8] == "_" or ja[i + 9] == "_" or ja[i + 10] \
                    == "_" or ja[i + 11] == "_" or ja[i + 12] == "_" or ja[i + 13] == "_" or ja[i + 14] == "_" or \
                    ja[i + 15] == "_" or ja[i + 16] == "_":
                print(ja[i], ja[i + 1], ja[i + 2], ja[i + 3], ja[i + 4], ja[i + 5], ja[i + 6], ja[i + 7], ja[i + 8],
                      ja[i + 9], ja[i + 10], ja[i + 11], ja[i + 12], ja[i + 13], ja[i + 14], ja[i + 15], ja[i + 16],
                      ja[i + 17], ja[i + 18], ja[i + 19], ja[i + 20], ja[i + 21], ja[i + 22], ja[i + 23], ja[i + 24],
                      ja[i + 25], ja[i + 26], ja[i + 27], ja[i + 28], ja[i + 29], ja[i + 30],
                      ja[i + 31], ja[i + 32], ja[i + 33], ja[i + 34], ja[i + 35], ja[i + 36], ja[i + 37], ja[i + 38],
                      ja[i + 39], ja[i + 40], ja[i + 41], ja[i + 42], ja[i + 43], ja[i + 44], ja[i + 45], ja[i + 46],
                      ja[i + 47], ja[i + 48], ja[i + 49], ja[i + 50], ja[i + 51], ja[i + 52], ja[i + 53],
                      ja[i + 54], ja[i + 55], ja[i + 56], ja[i + 57], ja[i + 58], ja[i + 59], ja[i + 60], ja[i + 61],
                      ja[i + 62], ja[i + 63], ja[i + 64], ja[i + 65], ja[i + 66], ja[i + 67], ja[i + 68],
                      file=open("output4.txt", "a+"), sep="")
                strix = open("output4.txt")
                striy = strix.readlines()
                string = striy[-1]
                res = [sub.replace('_', '') for sub in string]
                x = "".join(res)
                print(x[:17], file=open("oligo.txt", "a+"))
            else:
                print(ja[i], ja[i + 1], ja[i + 2], ja[i + 3], ja[i + 4], ja[i + 5], ja[i + 6], ja[i + 7], ja[i + 8],
                      ja[i + 9], ja[i + 10], ja[i + 11], ja[i + 12], ja[i + 13], ja[i + 14], ja[i + 15], ja[i + 16],
                      file=open("oligo.txt", "a+"), sep="")
        if i == 18:
            if ja[i - 1] == "_" or ja[i] == "_" or ja[i + 1] == "_" or ja[i + 2] == "_" or ja[i + 3] == "_" or \
                    ja[i + 4] == "_" or ja[i + 5] \
                    == "_" or ja[i + 6] == "_" or ja[i + 7] == "_" or ja[i + 8] == "_" or ja[i + 9] == "_" or \
                    ja[i + 10] == "_" or ja[i + 11] == "_" or ja[i + 12] == "_" or ja[i + 13] == "_" or \
                    ja[i + 14] == "_" or ja[i + 15] == "_":
                print(ja[i - 1], ja[i], ja[i + 1], ja[i + 2], ja[i + 3], ja[i + 4], ja[i + 5], ja[i + 6], ja[i + 7],
                      ja[i + 8],
                      ja[i + 9], ja[i + 10], ja[i + 11], ja[i + 12], ja[i + 13], ja[i + 14], ja[i + 15], ja[i + 16],
                      ja[i + 17], ja[i + 18], ja[i + 19], ja[i + 20], ja[i + 21], ja[i + 22], ja[i + 23], ja[i + 24],
                      ja[i + 25], ja[i + 26], ja[i + 27], ja[i + 28], ja[i + 29], ja[i + 30],
                      ja[i + 31], ja[i + 32], ja[i + 33], ja[i + 34], ja[i + 35], ja[i + 36], ja[i + 37], ja[i + 38],
                      ja[i + 39], ja[i + 40], ja[i + 41], ja[i + 42], ja[i + 43], ja[i + 44], ja[i + 45], ja[i + 46],
                      ja[i + 47], ja[i + 48], ja[i + 49], ja[i + 50], ja[i + 51], ja[i + 52], ja[i + 53],
                      ja[i + 54], ja[i + 55], ja[i + 56], ja[i + 57], ja[i + 58], ja[i + 59], ja[i + 60], ja[i + 61],
                      ja[i + 62], ja[i + 63], ja[i + 64], ja[i + 65], ja[i + 66], ja[i + 67],
                      file=open("output4.txt", "a+"), sep="")
                strix = open("output4.txt")
                striy = strix.readlines()
                string = striy[-1]
                res = [sub.replace('_', '') for sub in string]
                x = "".join(res)
                print(x[:17], file=open("oligo.txt", "a+"))
            else:
                print(ja[i - 1], ja[i], ja[i + 1], ja[i + 2], ja[i + 3], ja[i + 4], ja[i + 5], ja[i + 6], ja[i + 7],
                      ja[i + 8],
                      ja[i + 9], ja[i + 10], ja[i + 11], ja[i + 12], ja[i + 13], ja[i + 14], ja[i + 15],
                      file=open("oligo.txt", "a+"), sep="")
        if i == 19:
            if ja[i - 2] == "_" or ja[i - 1] == "_" or ja[i] == "_" or ja[i + 1] == "_" or ja[i + 2] == "_" or \
                    ja[i + 3] == "_" or ja[i + 4] == "_" or ja[i + 5] == "_" or ja[i + 6] == "_" or ja[i + 7] == "_" or \
                    ja[i + 8] == "_" or ja[i + 9] == "_" or ja[i + 10] \
                    == "_" or ja[i + 11] == "_" or ja[i + 12] == "_" or ja[i + 13] == "_" or ja[i + 14] == "_":
                print(ja[i - 2], ja[i - 1], ja[i], ja[i + 1], ja[i + 2], ja[i + 3], ja[i + 4], ja[i + 5], ja[i + 6],
                      ja[i + 7], ja[i + 8],
                      ja[i + 9], ja[i + 10], ja[i + 11], ja[i + 12], ja[i + 13], ja[i + 14], ja[i + 15], ja[i + 16],
                      ja[i + 17], ja[i + 18], ja[i + 19], ja[i + 20], ja[i + 21], ja[i + 22], ja[i + 23], ja[i + 24],
                      ja[i + 25], ja[i + 26], ja[i + 27], ja[i + 28], ja[i + 29], ja[i + 30],
                      ja[i + 31], ja[i + 32], ja[i + 33], ja[i + 34], ja[i + 35], ja[i + 36], ja[i + 37], ja[i + 38],
                      ja[i + 39], ja[i + 40], ja[i + 41], ja[i + 42], ja[i + 43], ja[i + 44], ja[i + 45], ja[i + 46],
                      ja[i + 47], ja[i + 48], ja[i + 49], ja[i + 50], ja[i + 51], ja[i + 52], ja[i + 53],
                      ja[i + 54], ja[i + 55], ja[i + 56], ja[i + 57], ja[i + 58], ja[i + 59], ja[i + 60], ja[i + 61],
                      ja[i + 62], ja[i + 63], ja[i + 64], ja[i + 65], ja[i + 66],
                      file=open("output4.txt", "a+"), sep="")
                strix = open("output4.txt")
                striy = strix.readlines()
                string = striy[-1]
                res = [sub.replace('_', '') for sub in string]
                x = "".join(res)
                print(x[:17], file=open("oligo.txt", "a+"))
            else:
                print(ja[i - 2], ja[i - 1], ja[i], ja[i + 1], ja[i + 2], ja[i + 3], ja[i + 4], ja[i + 5], ja[i + 6],
                      ja[i + 7], ja[i + 8],
                      ja[i + 9], ja[i + 10], ja[i + 11], ja[i + 12], ja[i + 13], ja[i + 14],
                      file=open("oligo.txt", "a+"), sep="")
        if i == 20:
            if ja[i - 3] == "_" or ja[i - 2] == "_" or ja[i - 1] == "_" or ja[i] == "_" or ja[i + 1] == "_" or \
                    ja[i + 2] == "_" or ja[i + 3] == "_" or ja[i + 4] == "_" or ja[i + 5] \
                    == "_" or ja[i + 6] == "_" or ja[i + 7] == "_" or ja[i + 8] == "_" or ja[i + 9] == "_" or ja[i + 10] \
                    == "_" or ja[i + 11] == "_" or ja[i + 12] == "_" or ja[i + 13] == "_":
                print(ja[i - 3], ja[i - 2], ja[i - 1], ja[i], ja[i + 1], ja[i + 2], ja[i + 3], ja[i + 4], ja[i + 5],
                      ja[i + 6], ja[i + 7], ja[i + 8],
                      ja[i + 9], ja[i + 10], ja[i + 11], ja[i + 12], ja[i + 13], ja[i + 14], ja[i + 15], ja[i + 16],
                      ja[i + 17], ja[i + 18], ja[i + 19], ja[i + 20], ja[i + 21], ja[i + 22], ja[i + 23], ja[i + 24],
                      ja[i + 25], ja[i + 26], ja[i + 27], ja[i + 28], ja[i + 29], ja[i + 30],
                      ja[i + 31], ja[i + 32], ja[i + 33], ja[i + 34], ja[i + 35], ja[i + 36], ja[i + 37], ja[i + 38],
                      ja[i + 39], ja[i + 40], ja[i + 41], ja[i + 42], ja[i + 43], ja[i + 44], ja[i + 45], ja[i + 46],
                      ja[i + 47], ja[i + 48], ja[i + 49], ja[i + 50], ja[i + 51], ja[i + 52], ja[i + 53],
                      ja[i + 54], ja[i + 55], ja[i + 56], ja[i + 57], ja[i + 58], ja[i + 59], ja[i + 60], ja[i + 61],
                      ja[i + 62], ja[i + 63], ja[i + 64], ja[i + 65],
                      file=open("output4.txt", "a+"), sep="")
                strix = open("output4.txt")
                striy = strix.readlines()
                string = striy[-1]
                res = [sub.replace('_', '') for sub in string]
                x = "".join(res)
                print(x[:17], file=open("oligo.txt", "a+"))
            else:
                print(ja[i - 3], ja[i - 2], ja[i - 1], ja[i], ja[i + 1], ja[i + 2], ja[i + 3], ja[i + 4], ja[i + 5],
                      ja[i + 6], ja[i + 7], ja[i + 8],
                      ja[i + 9], ja[i + 10], ja[i + 11], ja[i + 12], ja[i + 13],
                      file=open("oligo.txt", "a+"), sep="")
        if i == 21:
            if ja[i - 4] == "_" or ja[i - 3] == "_" or ja[i - 2] == "_" or ja[i - 1] == "_" or ja[i] == "_" or \
                    ja[i + 1] == "_" or ja[i + 2] == "_" or ja[i + 3] == "_" or ja[i + 4] == "_" or ja[i + 5] \
                    == "_" or ja[i + 6] == "_" or ja[i + 7] == "_" or ja[i + 8] == "_" or ja[i + 9] == "_" or \
                    ja[i + 10] == "_" or ja[i + 11] == "_" or ja[i + 12] == "_":
                print(ja[i - 4], ja[i - 3], ja[i - 2], ja[i - 1], ja[i], ja[i + 1], ja[i + 2], ja[i + 3], ja[i + 4],
                      ja[i + 5], ja[i + 6], ja[i + 7], ja[i + 8],
                      ja[i + 9], ja[i + 10], ja[i + 11], ja[i + 12], ja[i + 13], ja[i + 14], ja[i + 15], ja[i + 16],
                      ja[i + 17], ja[i + 18], ja[i + 19], ja[i + 20], ja[i + 21], ja[i + 22], ja[i + 23], ja[i + 24],
                      ja[i + 25], ja[i + 26], ja[i + 27], ja[i + 28], ja[i + 29], ja[i + 30],
                      ja[i + 31], ja[i + 32], ja[i + 33], ja[i + 34], ja[i + 35], ja[i + 36], ja[i + 37], ja[i + 38],
                      ja[i + 39], ja[i + 40], ja[i + 41], ja[i + 42], ja[i + 43], ja[i + 44], ja[i + 45], ja[i + 46],
                      ja[i + 47], ja[i + 48], ja[i + 49], ja[i + 50], ja[i + 51], ja[i + 52], ja[i + 53],
                      ja[i + 54], ja[i + 55], ja[i + 56], ja[i + 57], ja[i + 58], ja[i + 59], ja[i + 60], ja[i + 61],
                      ja[i + 62], ja[i + 63], ja[i + 64],
                      file=open("output4.txt", "a+"), sep="")
                strix = open("output4.txt")
                striy = strix.readlines()
                string = striy[-1]
                res = [sub.replace('_', '') for sub in string]
                x = "".join(res)
                print(x[:17], file=open("oligo.txt", "a+"))
            else:
                print(ja[i - 4], ja[i - 3], ja[i - 2], ja[i - 1], ja[i], ja[i + 1], ja[i + 2], ja[i + 3], ja[i + 4],
                      ja[i + 5], ja[i + 6], ja[i + 7], ja[i + 8],
                      ja[i + 9], ja[i + 10], ja[i + 11], ja[i + 12],
                      file=open("oligo.txt", "a+"), sep="")
        if i == 22:
            if ja[i - 5] == "_" or ja[i - 4] == "_" or ja[i - 3] == "_" or ja[i - 2] == "_" or ja[i - 1] == "_" or \
                    ja[i] == "_" or ja[i + 1] == "_" or ja[i + 2] == "_" or ja[i + 3] == "_" or ja[i + 4] == "_" or \
                    ja[i + 5] == "_" or ja[i + 6] == "_" or ja[i + 7] == "_" or ja[i + 8] == "_" or ja[i + 9] == "_" or \
                    ja[i + 10] == "_" or ja[i + 11] == "_":
                print(ja[i - 5], ja[i - 4], ja[i - 3], ja[i - 2], ja[i - 1], ja[i], ja[i + 1], ja[i + 2], ja[i + 3],
                      ja[i + 4], ja[i + 5], ja[i + 6], ja[i + 7], ja[i + 8],
                      ja[i + 9], ja[i + 10], ja[i + 11], ja[i + 12], ja[i + 13], ja[i + 14], ja[i + 15], ja[i + 16],
                      ja[i + 17], ja[i + 18], ja[i + 19], ja[i + 20], ja[i + 21], ja[i + 22], ja[i + 23], ja[i + 24],
                      ja[i + 25], ja[i + 26], ja[i + 27], ja[i + 28], ja[i + 29], ja[i + 30],
                      ja[i + 31], ja[i + 32], ja[i + 33], ja[i + 34], ja[i + 35], ja[i + 36], ja[i + 37], ja[i + 38],
                      ja[i + 39], ja[i + 40], ja[i + 41], ja[i + 42], ja[i + 43], ja[i + 44], ja[i + 45], ja[i + 46],
                      ja[i + 47], ja[i + 48], ja[i + 49], ja[i + 50], ja[i + 51], ja[i + 52], ja[i + 53],
                      ja[i + 54], ja[i + 55], ja[i + 56], ja[i + 57], ja[i + 58], ja[i + 59], ja[i + 60], ja[i + 61],
                      ja[i + 62], ja[i + 63],
                      file=open("output4.txt", "a+"), sep="")
                strix = open("output4.txt")
                striy = strix.readlines()
                string = striy[-1]
                res = [sub.replace('_', '') for sub in string]
                x = "".join(res)
                print(x[:17], file=open("oligo.txt", "a+"))
            else:
                print(ja[i - 5], ja[i - 4], ja[i - 3], ja[i - 2], ja[i - 1], ja[i], ja[i + 1], ja[i + 2], ja[i + 3],
                      ja[i + 4], ja[i + 5], ja[i + 6], ja[i + 7], ja[i + 8], ja[i + 9], ja[i + 10], ja[i + 11],
                      file=open("oligo.txt", "a+"), sep="")
        if i == 23:
            if ja[i - 6] == "_" or ja[i - 5] == "_" or ja[i - 4] == "_" or ja[i - 3] == "_" or ja[i - 2] == "_" or \
                    ja[i - 1] == "_" or ja[i] == "_" or ja[i + 1] == "_" or ja[i + 2] == "_" or ja[i + 3] == "_" or \
                    ja[i + 4] == "_" or ja[i + 5] \
                    == "_" or ja[i + 6] == "_" or ja[i + 7] == "_" or ja[i + 8] == "_" or ja[i + 9] == "_" or \
                    ja[i + 10] == "_":
                print(ja[i - 6], ja[i - 5], ja[i - 4], ja[i - 3], ja[i - 2], ja[i - 1], ja[i], ja[i + 1], ja[i + 2],
                      ja[i + 3], ja[i + 4], ja[i + 5], ja[i + 6], ja[i + 7], ja[i + 8],
                      ja[i + 9], ja[i + 10], ja[i + 11], ja[i + 12], ja[i + 13], ja[i + 14], ja[i + 15], ja[i + 16],
                      ja[i + 17], ja[i + 18], ja[i + 19], ja[i + 20], ja[i + 21], ja[i + 22], ja[i + 23], ja[i + 24],
                      ja[i + 25], ja[i + 26], ja[i + 27], ja[i + 28], ja[i + 29], ja[i + 30],
                      ja[i + 31], ja[i + 32], ja[i + 33], ja[i + 34], ja[i + 35], ja[i + 36], ja[i + 37], ja[i + 38],
                      ja[i + 39], ja[i + 40], ja[i + 41], ja[i + 42], ja[i + 43], ja[i + 44], ja[i + 45], ja[i + 46],
                      ja[i + 47], ja[i + 48], ja[i + 49], ja[i + 50], ja[i + 51], ja[i + 52], ja[i + 53],
                      ja[i + 54], ja[i + 55], ja[i + 56], ja[i + 57], ja[i + 58], ja[i + 59], ja[i + 60], ja[i + 61],
                      ja[i + 62],
                      file=open("output4.txt", "a+"), sep="")
                strix = open("output4.txt")
                striy = strix.readlines()
                string = striy[-1]
                res = [sub.replace('_', '') for sub in string]
                x = "".join(res)
                print(x[:17], file=open("oligo.txt", "a+"))
            else:
                print(ja[i - 6], ja[i - 5], ja[i - 4], ja[i - 3], ja[i - 2], ja[i - 1], ja[i], ja[i + 1], ja[i + 2],
                      ja[i + 3], ja[i + 4], ja[i + 5], ja[i + 6], ja[i + 7], ja[i + 8],
                      ja[i + 9], ja[i + 10],
                      file=open("oligo.txt", "a+"), sep="")
        if i == 24:
            if ja[i - 7] == "_" or ja[i - 6] == "_" or ja[i - 5] == "_" or ja[i - 4] == "_" or ja[i - 3] == "_" or \
                    ja[i - 2] == "_" or ja[i - 1] == "_" or ja[i] == "_" or ja[i + 1] == "_" or ja[i + 2] == "_" or \
                    ja[i + 3] == "_" or ja[i + 4] == "_" or ja[i + 5] \
                    == "_" or ja[i + 6] == "_" or ja[i + 7] == "_" or ja[i + 8] == "_" or ja[i + 9] == "_":
                print(ja[i - 7], ja[i - 6], ja[i - 5], ja[i - 4], ja[i - 3], ja[i - 2], ja[i - 1], ja[i], ja[i + 1],
                      ja[i + 2], ja[i + 3], ja[i + 4], ja[i + 5], ja[i + 6], ja[i + 7], ja[i + 8],
                      ja[i + 9], ja[i + 10], ja[i + 11], ja[i + 12], ja[i + 13], ja[i + 14], ja[i + 15], ja[i + 16],
                      ja[i + 17], ja[i + 18], ja[i + 19], ja[i + 20], ja[i + 21], ja[i + 22], ja[i + 23], ja[i + 24],
                      ja[i + 25], ja[i + 26], ja[i + 27], ja[i + 28], ja[i + 29], ja[i + 30],
                      ja[i + 31], ja[i + 32], ja[i + 33], ja[i + 34], ja[i + 35], ja[i + 36], ja[i + 37], ja[i + 38],
                      ja[i + 39], ja[i + 40], ja[i + 41], ja[i + 42], ja[i + 43], ja[i + 44], ja[i + 45], ja[i + 46],
                      ja[i + 47], ja[i + 48], ja[i + 49], ja[i + 50], ja[i + 51], ja[i + 52], ja[i + 53],
                      ja[i + 54], ja[i + 55], ja[i + 56], ja[i + 57], ja[i + 58], ja[i + 59], ja[i + 60], ja[i + 61],
                      file=open("output4.txt", "a+"), sep="")
                strix = open("output4.txt")
                striy = strix.readlines()
                string = striy[-1]
                res = [sub.replace('_', '') for sub in string]
                x = "".join(res)
                print(x[:17], file=open("oligo.txt", "a+"))
            else:
                print(ja[i - 7], ja[i - 6], ja[i - 5], ja[i - 4], ja[i - 3], ja[i - 2], ja[i - 1], ja[i], ja[i + 1],
                      ja[i + 2], ja[i + 3], ja[i + 4], ja[i + 5], ja[i + 6], ja[i + 7], ja[i + 8],
                      ja[i + 9], file=open("oligo.txt", "a+"), sep="")
        if i == 25:
            if ja[i - 8] == "_" or ja[i - 7] == "_" or ja[i - 6] == "_" or ja[i - 5] == "_" or ja[i - 4] == "_" or \
                    ja[i - 3] == "_" or ja[i - 2] == "_" or ja[i - 1] == "_" or ja[i] == "_" or ja[i + 1] == "_" or \
                    ja[i + 2] == "_" or ja[i + 3] == "_" or ja[i + 4] == "_" or ja[i + 5] \
                    == "_" or ja[i + 6] == "_" or ja[i + 7] == "_" or ja[i + 8] == "_":
                print(ja[i - 8], ja[i - 7], ja[i - 6], ja[i - 5], ja[i - 4], ja[i - 3], ja[i - 2], ja[i - 1], ja[i],
                      ja[i + 1], ja[i + 2], ja[i + 3], ja[i + 4], ja[i + 5], ja[i + 6], ja[i + 7], ja[i + 8],
                      ja[i + 9], ja[i + 10], ja[i + 11], ja[i + 12], ja[i + 13], ja[i + 14], ja[i + 15], ja[i + 16],
                      ja[i + 17], ja[i + 18], ja[i + 19], ja[i + 20], ja[i + 21], ja[i + 22], ja[i + 23], ja[i + 24],
                      ja[i + 25], ja[i + 26], ja[i + 27], ja[i + 28], ja[i + 29], ja[i + 30],
                      ja[i + 31], ja[i + 32], ja[i + 33], ja[i + 34], ja[i + 35], ja[i + 36], ja[i + 37], ja[i + 38],
                      ja[i + 39], ja[i + 40], ja[i + 41], ja[i + 42], ja[i + 43], ja[i + 44], ja[i + 45], ja[i + 46],
                      ja[i + 47], ja[i + 48], ja[i + 49], ja[i + 50], ja[i + 51], ja[i + 52], ja[i + 53],
                      ja[i + 54], ja[i + 55], ja[i + 56], ja[i + 57], ja[i + 58], ja[i + 59], ja[i + 60],
                      file=open("output4.txt", "a+"), sep="")
                strix = open("output4.txt")
                striy = strix.readlines()
                string = striy[-1]
                res = [sub.replace('_', '') for sub in string]
                x = "".join(res)
                print(x[:17], file=open("oligo.txt", "a+"))
            else:
                print(ja[i - 8], ja[i - 7], ja[i - 6], ja[i - 5], ja[i - 4], ja[i - 3], ja[i - 2], ja[i - 1], ja[i],
                      ja[i + 1], ja[i + 2], ja[i + 3], ja[i + 4], ja[i + 5], ja[i + 6], ja[i + 7], ja[i + 8],
                      file=open("oligo.txt", "a+"), sep="")

        if i == 26:
            if ja[i - 9] == "_" or ja[i - 8] == "_" or ja[i - 7] == "_" or ja[i - 6] == "_" or ja[i - 5] == "_" or \
                    ja[i - 4] == "_" or ja[i - 3] == "_" or ja[i - 2] == "_" or ja[i - 1] == "_" or ja[i] == "_" or \
                    ja[i + 1] == "_" or ja[i + 2] == "_" or ja[i + 3] == "_" or ja[i + 4] == "_" or ja[i + 5] \
                    == "_" or ja[i + 6] == "_" or ja[i + 7] == "_":
                print(ja[i - 9], ja[i - 8], ja[i - 7], ja[i - 6], ja[i - 5], ja[i - 4], ja[i - 3], ja[i - 2], ja[i - 1],
                      ja[i], ja[i + 1], ja[i + 2], ja[i + 3], ja[i + 4], ja[i + 5], ja[i + 6], ja[i + 7], ja[i + 8],
                      ja[i + 9], ja[i + 10], ja[i + 11], ja[i + 12], ja[i + 13], ja[i + 14], ja[i + 15], ja[i + 16],
                      ja[i + 17], ja[i + 18], ja[i + 19], ja[i + 20], ja[i + 21], ja[i + 22], ja[i + 23], ja[i + 24],
                      ja[i + 25], ja[i + 26], ja[i + 27], ja[i + 28], ja[i + 29], ja[i + 30],
                      ja[i + 31], ja[i + 32], ja[i + 33], ja[i + 34], ja[i + 35], ja[i + 36], ja[i + 37], ja[i + 38],
                      ja[i + 39], ja[i + 40], ja[i + 41], ja[i + 42], ja[i + 43], ja[i + 44], ja[i + 45], ja[i + 46],
                      ja[i + 47], ja[i + 48], ja[i + 49], ja[i + 50], ja[i + 51], ja[i + 52], ja[i + 53],
                      ja[i + 54], ja[i + 55], ja[i + 56], ja[i + 57], ja[i + 58], ja[i + 59],
                      file=open("output4.txt", "a+"), sep="")
                strix = open("output4.txt")
                striy = strix.readlines()
                string = striy[-1]
                res = [sub.replace('_', '') for sub in string]
                x = "".join(res)
                print(x[:17], file=open("oligo.txt", "a+"))
            else:
                print(ja[i - 9], ja[i - 8], ja[i - 7], ja[i - 6], ja[i - 5], ja[i - 4], ja[i - 3], ja[i - 2], ja[i - 1],
                      ja[i], ja[i + 1], ja[i + 2], ja[i + 3], ja[i + 4], ja[i + 5], ja[i + 6], ja[i + 7],
                      file=open("oligo.txt", "a+"), sep="")
        if i == 27:
            if ja[i - 10] == "_" or ja[i - 9] == "_" or ja[i - 8] == "_" or ja[i - 7] == "_" or ja[i - 6] == "_" or \
                    ja[i - 5] == "_" or ja[i - 4] == "_" or ja[i - 3] == "_" or ja[i - 2] == "_" or ja[i - 1] == "_" or \
                    ja[i] == "_" or ja[i + 1] == "_" or ja[i + 2] == "_" or ja[i + 3] == "_" or ja[i + 4] == "_" or \
                    ja[i + 5] == "_" or ja[i + 6] == "_":
                print(ja[i - 10], ja[i - 9], ja[i - 8], ja[i - 7], ja[i - 6], ja[i - 5], ja[i - 4], ja[i - 3],
                      ja[i - 2], ja[i - 1],
                      ja[i], ja[i + 1], ja[i + 2], ja[i + 3], ja[i + 4], ja[i + 5], ja[i + 6], ja[i + 7], ja[i + 8],
                      ja[i + 9], ja[i + 10], ja[i + 11], ja[i + 12], ja[i + 13], ja[i + 14], ja[i + 15], ja[i + 16],
                      ja[i + 17], ja[i + 18], ja[i + 19], ja[i + 20], ja[i + 21], ja[i + 22], ja[i + 23], ja[i + 24],
                      ja[i + 25], ja[i + 26], ja[i + 27], ja[i + 28], ja[i + 29], ja[i + 30],
                      ja[i + 31], ja[i + 32], ja[i + 33], ja[i + 34], ja[i + 35], ja[i + 36], ja[i + 37], ja[i + 38],
                      ja[i + 39], ja[i + 40], ja[i + 41], ja[i + 42], ja[i + 43], ja[i + 44], ja[i + 45], ja[i + 46],
                      ja[i + 47], ja[i + 48], ja[i + 49], ja[i + 50], ja[i + 51], ja[i + 52], ja[i + 53],
                      ja[i + 54], ja[i + 55], ja[i + 56], ja[i + 57], ja[i + 58],
                      file=open("output4.txt", "a+"), sep="")
                strix = open("output4.txt")
                striy = strix.readlines()
                string = striy[-1]
                res = [sub.replace('_', '') for sub in string]
                x = "".join(res)
                print(x[:17], file=open("oligo.txt", "a+"))
            else:
                print(ja[i - 10], ja[i - 9], ja[i - 8], ja[i - 7], ja[i - 6], ja[i - 5], ja[i - 4], ja[i - 3],
                      ja[i - 2], ja[i - 1],
                      ja[i], ja[i + 1], ja[i + 2], ja[i + 3], ja[i + 4], ja[i + 5], ja[i + 6],
                      file=open("oligo.txt", "a+"), sep="")
        if i == 28:
            if ja[i - 11] == "_" or ja[i - 10] == "_" or ja[i - 9] == "_" or ja[i - 8] == "_" or ja[i - 7] == "_" or \
                    ja[i - 6] == "_" or \
                    ja[i - 5] == "_" or ja[i - 4] == "_" or ja[i - 3] == "_" or ja[i - 2] == "_" or ja[i - 1] == "_" or \
                    ja[i] == "_" or ja[i + 1] == "_" or ja[i + 2] == "_" or ja[i + 3] == "_" or ja[i + 4] == "_" or \
                    ja[i + 5] == "_":
                print(ja[i - 11], ja[i - 9], ja[i - 8], ja[i - 7], ja[i - 6], ja[i - 5], ja[i - 4], ja[i - 3],
                      ja[i - 2], ja[i - 1],
                      ja[i], ja[i + 1], ja[i + 2], ja[i + 3], ja[i + 4], ja[i + 5], ja[i + 6], ja[i + 7], ja[i + 8],
                      ja[i + 9], ja[i + 10], ja[i + 11], ja[i + 12], ja[i + 13], ja[i + 14], ja[i + 15], ja[i + 16],
                      ja[i + 17], ja[i + 18], ja[i + 19], ja[i + 20], ja[i + 21], ja[i + 22], ja[i + 23], ja[i + 24],
                      ja[i + 25], ja[i + 26], ja[i + 27], ja[i + 28], ja[i + 29], ja[i + 30],
                      ja[i + 31], ja[i + 32], ja[i + 33], ja[i + 34], ja[i + 35], ja[i + 36], ja[i + 37], ja[i + 38],
                      ja[i + 39], ja[i + 40], ja[i + 41], ja[i + 42], ja[i + 43], ja[i + 44], ja[i + 45], ja[i + 46],
                      ja[i + 47], ja[i + 48], ja[i + 49], ja[i + 50], ja[i + 51], ja[i + 52], ja[i + 53],
                      ja[i + 54], ja[i + 55], ja[i + 56], ja[i + 57],
                      file=open("output4.txt", "a+"), sep="")
                strix = open("output4.txt")
                striy = strix.readlines()
                string = striy[-1]
                res = [sub.replace('_', '') for sub in string]
                x = "".join(res)
                print(x[:17], file=open("oligo.txt", "a+"))
            else:
                print(ja[i - 11], ja[i - 10], ja[i - 9], ja[i - 8], ja[i - 7], ja[i - 6], ja[i - 5], ja[i - 4],
                      ja[i - 3],
                      ja[i - 2], ja[i - 1],
                      ja[i], ja[i + 1], ja[i + 2], ja[i + 3], ja[i + 4], ja[i + 5],
                      file=open("oligo.txt", "a+"), sep="")
        if i == 29:
            if ja[i - 12] == "_" or ja[i - 11] == "_" or ja[i - 10] == "_" or ja[i - 9] == "_" or ja[i - 8] == "_" or \
                    ja[i - 7] == "_" or ja[i - 6] == "_" or \
                    ja[i - 5] == "_" or ja[i - 4] == "_" or ja[i - 3] == "_" or ja[i - 2] == "_" or ja[i - 1] == "_" or \
                    ja[i] == "_" or ja[i + 1] == "_" or ja[i + 2] == "_" or ja[i + 3] == "_" or ja[i + 4] == "_":
                print(ja[i - 12], ja[i - 11], ja[i - 9], ja[i - 8], ja[i - 7], ja[i - 6], ja[i - 5], ja[i - 4],
                      ja[i - 3], ja[i - 2], ja[i - 1],
                      ja[i], ja[i + 1], ja[i + 2], ja[i + 3], ja[i + 4], ja[i + 5], ja[i + 6], ja[i + 7], ja[i + 8],
                      ja[i + 9], ja[i + 10], ja[i + 11], ja[i + 12], ja[i + 13], ja[i + 14], ja[i + 15], ja[i + 16],
                      ja[i + 17], ja[i + 18], ja[i + 19], ja[i + 20], ja[i + 21], ja[i + 22], ja[i + 23], ja[i + 24],
                      ja[i + 25], ja[i + 26], ja[i + 27], ja[i + 28], ja[i + 29], ja[i + 30],
                      ja[i + 31], ja[i + 32], ja[i + 33], ja[i + 34], ja[i + 35], ja[i + 36], ja[i + 37], ja[i + 38],
                      ja[i + 39], ja[i + 40], ja[i + 41], ja[i + 42], ja[i + 43], ja[i + 44], ja[i + 45], ja[i + 46],
                      ja[i + 47], ja[i + 48], ja[i + 49], ja[i + 50], ja[i + 51], ja[i + 52], ja[i + 53],
                      ja[i + 54], ja[i + 55], ja[i + 56],
                      file=open("output4.txt", "a+"), sep="")
                strix = open("output4.txt")
                striy = strix.readlines()
                string = striy[-1]
                res = [sub.replace('_', '') for sub in string]
                x = "".join(res)
                print(x[:17], file=open("oligo.txt", "a+"))
            else:
                print(ja[i - 12], ja[i - 11], ja[i - 10], ja[i - 9], ja[i - 8], ja[i - 7], ja[i - 6], ja[i - 5],
                      ja[i - 4], ja[i - 3],
                      ja[i - 2], ja[i - 1],
                      ja[i], ja[i + 1], ja[i + 2], ja[i + 3], ja[i + 4],
                      file=open("oligo.txt", "a+"), sep="")
        if i == 30:
            if ja[i - 13] == "_" or ja[i - 12] == "_" or ja[i - 11] == "_" or ja[i - 10] == "_" or ja[i - 9] == "_" or \
                    ja[i - 8] == "_" or ja[i - 7] == "_" or ja[i - 6] == "_" or \
                    ja[i - 5] == "_" or ja[i - 4] == "_" or ja[i - 3] == "_" or ja[i - 2] == "_" or ja[i - 1] == "_" or \
                    ja[i] == "_" or ja[i + 1] == "_" or ja[i + 2] == "_" or ja[i + 3] == "_":
                print(ja[i - 13], ja[i - 12], ja[i - 11], ja[i - 9], ja[i - 8], ja[i - 7], ja[i - 6], ja[i - 5],
                      ja[i - 4], ja[i - 3], ja[i - 2], ja[i - 1],
                      ja[i], ja[i + 1], ja[i + 2], ja[i + 3], ja[i + 4], ja[i + 5], ja[i + 6], ja[i + 7], ja[i + 8],
                      ja[i + 9], ja[i + 10], ja[i + 11], ja[i + 12], ja[i + 13], ja[i + 14], ja[i + 15], ja[i + 16],
                      ja[i + 17], ja[i + 18], ja[i + 19], ja[i + 20], ja[i + 21], ja[i + 22], ja[i + 23], ja[i + 24],
                      ja[i + 25], ja[i + 26], ja[i + 27], ja[i + 28], ja[i + 29], ja[i + 30],
                      ja[i + 31], ja[i + 32], ja[i + 33], ja[i + 34], ja[i + 35], ja[i + 36], ja[i + 37], ja[i + 38],
                      ja[i + 39], ja[i + 40], ja[i + 41], ja[i + 42], ja[i + 43], ja[i + 44], ja[i + 45], ja[i + 46],
                      ja[i + 47], ja[i + 48], ja[i + 49], ja[i + 50], ja[i + 51], ja[i + 52], ja[i + 53],
                      ja[i + 54], ja[i + 55],
                      file=open("output4.txt", "a+"), sep="")
                strix = open("output4.txt")
                striy = strix.readlines()
                string = striy[-1]
                res = [sub.replace('_', '') for sub in string]
                x = "".join(res)
                print(x[:17], file=open("oligo.txt", "a+"))
            else:
                print(ja[i - 13], ja[i - 12], ja[i - 11], ja[i - 10], ja[i - 9], ja[i - 8], ja[i - 7], ja[i - 6],
                      ja[i - 5], ja[i - 4], ja[i - 3],
                      ja[i - 2], ja[i - 1], ja[i], ja[i + 1], ja[i + 2], ja[i + 3],
                      file=open("oligo.txt", "a+"), sep="")
        if i == 31:
            if ja[i - 14] == "_" or ja[i - 13] == "_" or ja[i - 12] == "_" or ja[i - 11] == "_" or ja[i - 10] == "_" or \
                    ja[i - 9] == "_" or \
                    ja[i - 8] == "_" or ja[i - 7] == "_" or ja[i - 6] == "_" or \
                    ja[i - 5] == "_" or ja[i - 4] == "_" or ja[i - 3] == "_" or ja[i - 2] == "_" or ja[i - 1] == "_" or \
                    ja[i] == "_" or ja[i + 1] == "_" or ja[i + 2] == "_":
                print(ja[i - 14], ja[i - 13], ja[i - 12], ja[i - 11], ja[i - 9], ja[i - 8], ja[i - 7], ja[i - 6],
                      ja[i - 5],
                      ja[i - 4], ja[i - 3], ja[i - 2], ja[i - 1],
                      ja[i], ja[i + 1], ja[i + 2], ja[i + 3], ja[i + 4], ja[i + 5], ja[i + 6], ja[i + 7], ja[i + 8],
                      ja[i + 9], ja[i + 10], ja[i + 11], ja[i + 12], ja[i + 13], ja[i + 14], ja[i + 15], ja[i + 16],
                      ja[i + 17], ja[i + 18], ja[i + 19], ja[i + 20], ja[i + 21], ja[i + 22], ja[i + 23], ja[i + 24],
                      ja[i + 25], ja[i + 26], ja[i + 27], ja[i + 28], ja[i + 29], ja[i + 30],
                      ja[i + 31], ja[i + 32], ja[i + 33], ja[i + 34], ja[i + 35], ja[i + 36], ja[i + 37], ja[i + 38],
                      ja[i + 39], ja[i + 40], ja[i + 41], ja[i + 42], ja[i + 43], ja[i + 44], ja[i + 45], ja[i + 46],
                      ja[i + 47], ja[i + 48], ja[i + 49], ja[i + 50], ja[i + 51], ja[i + 52], ja[i + 53],
                      ja[i + 54], file=open("output4.txt", "a+"), sep="")
                strix = open("output4.txt")
                striy = strix.readlines()
                string = striy[-1]
                res = [sub.replace('_', '') for sub in string]
                x = "".join(res)
                print(x[:17], file=open("oligo.txt", "a+"))
            else:
                print(ja[i - 14], ja[i - 13], ja[i - 12], ja[i - 11], ja[i - 10], ja[i - 9], ja[i - 8], ja[i - 7],
                      ja[i - 6],
                      ja[i - 5], ja[i - 4], ja[i - 3],
                      ja[i - 2], ja[i - 1], ja[i], ja[i + 1], ja[i + 2], file=open("oligo.txt", "a+"), sep="")
        if i == 32:
            if ja[i - 15] == "_" or ja[i - 14] == "_" or ja[i - 13] == "_" or ja[i - 12] == "_" or ja[i - 11] == "_" or \
                    ja[i - 10] == "_" or ja[i - 9] == "_" or ja[i - 8] == "_" or ja[i - 7] == "_" or ja[i - 6] == "_" or \
                    ja[i - 5] == "_" or ja[i - 4] == "_" or ja[i - 3] == "_" or ja[i - 2] == "_" or ja[i - 1] == "_" or \
                    ja[i] == "_" or ja[i + 1] == "_":
                print(ja[i - 15], ja[i - 14], ja[i - 13], ja[i - 12], ja[i - 11], ja[i - 9], ja[i - 8], ja[i - 7],
                      ja[i - 6], ja[i - 5],
                      ja[i - 4], ja[i - 3], ja[i - 2], ja[i - 1],
                      ja[i], ja[i + 1], ja[i + 2], ja[i + 3], ja[i + 4], ja[i + 5], ja[i + 6], ja[i + 7], ja[i + 8],
                      ja[i + 9], ja[i + 10], ja[i + 11], ja[i + 12], ja[i + 13], ja[i + 14], ja[i + 15], ja[i + 16],
                      ja[i + 17], ja[i + 18], ja[i + 19], ja[i + 20], ja[i + 21], ja[i + 22], ja[i + 23], ja[i + 24],
                      ja[i + 25], ja[i + 26], ja[i + 27], ja[i + 28], ja[i + 29], ja[i + 30],
                      ja[i + 31], ja[i + 32], ja[i + 33], ja[i + 34], ja[i + 35], ja[i + 36], ja[i + 37], ja[i + 38],
                      ja[i + 39], ja[i + 40], ja[i + 41], ja[i + 42], ja[i + 43], ja[i + 44], ja[i + 45], ja[i + 46],
                      ja[i + 47], ja[i + 48], ja[i + 49], ja[i + 50], ja[i + 51], ja[i + 52], ja[i + 53],
                      file=open("output4.txt", "a+"), sep="")
                strix = open("output4.txt")
                striy = strix.readlines()
                string = striy[-1]
                res = [sub.replace('_', '') for sub in string]
                x = "".join(res)
                print(x[:17], file=open("oligo.txt", "a+"))
            else:
                print(ja[i - 15], ja[i - 14], ja[i - 13], ja[i - 12], ja[i - 11], ja[i - 10], ja[i - 9], ja[i - 8],
                      ja[i - 7], ja[i - 6], ja[i - 5], ja[i - 4], ja[i - 3],
                      ja[i - 2], ja[i - 1], ja[i], ja[i + 1], file=open("oligo.txt", "a+"), sep="")
        if i == 33:
            if ja[i - 16] == "_" or ja[i - 15] == "_" or ja[i - 14] == "_" or ja[i - 13] == "_" or ja[i - 12] == "_" or \
                    ja[i - 11] == "_" or \
                    ja[i - 10] == "_" or ja[i - 9] == "_" or ja[i - 8] == "_" or ja[i - 7] == "_" or ja[i - 6] == "_" or \
                    ja[i - 5] == "_" or ja[i - 4] == "_" or ja[i - 3] == "_" or ja[i - 2] == "_" or ja[i - 1] == "_" or \
                    ja[i] == "_":
                print(ja[i - 16], ja[i - 15], ja[i - 14], ja[i - 13], ja[i - 12], ja[i - 11], ja[i - 9], ja[i - 8],
                      ja[i - 7], ja[i - 6], ja[i - 5],
                      ja[i - 4], ja[i - 3], ja[i - 2], ja[i - 1],
                      ja[i], ja[i + 1], ja[i + 2], ja[i + 3], ja[i + 4], ja[i + 5], ja[i + 6], ja[i + 7], ja[i + 8],
                      ja[i + 9], ja[i + 10], ja[i + 11], ja[i + 12], ja[i + 13], ja[i + 14], ja[i + 15], ja[i + 16],
                      ja[i + 17], ja[i + 18], ja[i + 19], ja[i + 20], ja[i + 21], ja[i + 22], ja[i + 23], ja[i + 24],
                      ja[i + 25], ja[i + 26], ja[i + 27], ja[i + 28], ja[i + 29], ja[i + 30],
                      ja[i + 31], ja[i + 32], ja[i + 33], ja[i + 34], ja[i + 35], ja[i + 36], ja[i + 37], ja[i + 38],
                      ja[i + 39], ja[i + 40], ja[i + 41], ja[i + 42], ja[i + 43], ja[i + 44], ja[i + 45], ja[i + 46],
                      ja[i + 47], ja[i + 48], ja[i + 49], ja[i + 50], ja[i + 51], ja[i + 52],
                      file=open("output4.txt", "a+"), sep="")
                strix = open("output4.txt")
                striy = strix.readlines()
                string = striy[-1]
                res = [sub.replace('_', '') for sub in string]
                x = "".join(res)
                print(x[:17], file=open("oligo.txt", "a+"))
            else:
                print(ja[i - 16], ja[i - 15], ja[i - 14], ja[i - 13], ja[i - 12], ja[i - 11], ja[i - 10], ja[i - 9],
                      ja[i - 8],
                      ja[i - 7], ja[i - 6], ja[i - 5], ja[i - 4], ja[i - 3],
                      ja[i - 2], ja[i - 1], ja[i], file=open("oligo.txt", "a+"), sep="")
        if i == 34:
            if ja[i] == "_" or ja[i + 1] == "_" or ja[i + 2] == "_" or ja[i + 3] == "_" or ja[i + 4] == "_" or ja[i + 5] \
                    == "_" or ja[i + 6] == "_" or ja[i + 7] == "_" or ja[i + 8] == "_" or ja[i + 9] == "_" or ja[i + 10] \
                    == "_" or ja[i + 11] == "_" or ja[i + 12] == "_" or ja[i + 13] == "_" or ja[i + 14] == "_" or \
                    ja[i + 15] == "_" or ja[i + 16] == "_":
                print(ja[i], ja[i + 1], ja[i + 2], ja[i + 3], ja[i + 4], ja[i + 5], ja[i + 6], ja[i + 7], ja[i + 8],
                      ja[i + 9], ja[i + 10], ja[i + 11], ja[i + 12], ja[i + 13], ja[i + 14], ja[i + 15], ja[i + 16],
                      ja[i + 17], ja[i + 18], ja[i + 19], ja[i + 20], ja[i + 21], ja[i + 22], ja[i + 23], ja[i + 24],
                      ja[i + 25], ja[i + 26], ja[i + 27], ja[i + 28], ja[i + 29], ja[i + 30],
                      ja[i + 31], ja[i + 32], ja[i + 33], ja[i + 34], ja[i + 35], ja[i + 36], ja[i + 37], ja[i + 38],
                      ja[i + 39], ja[i + 40], ja[i + 41], ja[i + 42], ja[i + 43], ja[i + 44], ja[i + 45], ja[i + 46],
                      ja[i + 47], ja[i + 48], ja[i + 49], ja[i + 50], ja[i + 51], ja[i + 52], ja[i + 53],
                      ja[i + 54], ja[i + 55], ja[i + 56], ja[i + 57], ja[i + 58], ja[i + 59], ja[i + 60], ja[i + 61],
                      ja[i + 62], ja[i + 63], ja[i + 64], ja[i + 65], ja[i + 66], ja[i + 67], ja[i + 68],
                      file=open("output4.txt", "a+"), sep="")
                strix = open("output4.txt")
                striy = strix.readlines()
                string = striy[-1]
                res = [sub.replace('_', '') for sub in string]
                x = "".join(res)
                print(x[:17], file=open("oligo.txt", "a+"))
            else:
                print(ja[i], ja[i + 1], ja[i + 2], ja[i + 3], ja[i + 4], ja[i + 5], ja[i + 6], ja[i + 7], ja[i + 8],
                      ja[i + 9], ja[i + 10], ja[i + 11], ja[i + 12], ja[i + 13], ja[i + 14], ja[i + 15], ja[i + 16],
                      file=open("oligo.txt", "a+"), sep="")
        if i == ln:
            if ja[i - 16] == "_" or ja[i - 15] == "_" or ja[i - 14] == "_" or ja[i - 13] == "_" or ja[i - 12] == "_" or \
                    ja[i - 11] == "_" or ja[i - 10] == "_" or ja[i - 9] == "_" or ja[i - 8] == "_" or \
                    ja[i - 7] == "_" or ja[i - 6] == "_" or ja[i - 5] == "_" or \
                    ja[i - 4] == "_" or ja[i - 3] == "_" or ja[i - 2] == "_" or \
                    ja[i - 1] == "_" or ja[i] == "_":
                print(ja[i - 68], ja[i - 67], ja[i - 66], ja[i - 65], ja[i - 64], ja[i - 63], ja[i - 62], ja[i - 61],
                      ja[i - 60],
                      ja[i - 59], ja[i - 58], ja[i - 57], ja[i - 56], ja[i - 55], ja[i - 54], ja[i - 53], ja[i - 52],
                      ja[i - 51],
                      ja[i - 50], ja[i - 49], ja[i - 48], ja[i - 47], ja[i - 46], ja[i - 45], ja[i - 44], ja[i - 43],
                      ja[i - 42],
                      ja[i - 41], ja[i - 40], ja[i - 39], ja[i - 38], ja[i - 37], ja[i - 36], ja[i - 35], ja[i - 34],
                      ja[i - 33],
                      ja[i - 32], ja[i - 31], ja[i - 30], ja[i - 29], ja[i - 28], ja[i - 27], ja[i - 26], ja[i - 25],
                      ja[i - 24],
                      ja[i - 23], ja[i - 22], ja[i - 21], ja[i - 20], ja[i - 19], ja[i - 18], ja[i - 17], ja[i - 16],
                      ja[i - 15],
                      ja[i - 14], ja[i - 13], ja[i - 12], ja[i - 11], ja[i - 10], ja[i - 9], ja[i - 8], ja[i - 7],
                      ja[i - 6],
                      ja[i - 5], ja[i - 4], ja[i - 3], ja[i - 2], ja[i - 1], ja[i], file=open("output4.txt", "a+"),
                      sep="")
                strix = open("output4.txt")
                striy = strix.readlines()
                string = striy[-1]
                res = [sub.replace('_', '') for sub in string]
                x = "".join(res)
                print(x[-17:], file=open("oligo.txt", "a+"))
            else:
                print(ja[i - 16], ja[i - 15], ja[i - 14], ja[i - 13], ja[i - 11], ja[i - 10], ja[i - 9], ja[i - 8],
                      ja[i - 7], ja[i - 6], ja[i - 5], ja[i - 4], ja[i - 3], ja[i - 2], ja[i - 1], ja[i],
                      file=open("oligo.txt", "a+"), sep="")
        if i == (ln - 1):
            if ja[i - 15] == "_" or ja[i - 14] == "_" or ja[i - 13] == "_" or ja[i - 12] == "_" or \
                    ja[i - 11] == "_" or ja[i - 10] == "_" or ja[i - 9] == "_" or ja[i - 8] == "_" or \
                    ja[i - 7] == "_" or ja[i - 6] == "_" or ja[i - 5] == "_" or \
                    ja[i - 4] == "_" or ja[i - 3] == "_" or ja[i - 2] == "_" or ja[i - 1] == "_" or ja[i] == "_" or \
                    ja[i + 1] == "_":
                print(ja[i - 67], ja[i - 66], ja[i - 65], ja[i - 64], ja[i - 63], ja[i - 62], ja[i - 61],
                      ja[i - 60],
                      ja[i - 59], ja[i - 58], ja[i - 57], ja[i - 56], ja[i - 55], ja[i - 54], ja[i - 53], ja[i - 52],
                      ja[i - 51],
                      ja[i - 50], ja[i - 49], ja[i - 48], ja[i - 47], ja[i - 46], ja[i - 45], ja[i - 44], ja[i - 43],
                      ja[i - 42],
                      ja[i - 41], ja[i - 40], ja[i - 39], ja[i - 38], ja[i - 37], ja[i - 36], ja[i - 35], ja[i - 34],
                      ja[i - 33],
                      ja[i - 32], ja[i - 31], ja[i - 30], ja[i - 29], ja[i - 28], ja[i - 27], ja[i - 26], ja[i - 25],
                      ja[i - 24],
                      ja[i - 23], ja[i - 22], ja[i - 21], ja[i - 20], ja[i - 19], ja[i - 18], ja[i - 17], ja[i - 16],
                      ja[i - 15],
                      ja[i - 14], ja[i - 13], ja[i - 12], ja[i - 11], ja[i - 10], ja[i - 9], ja[i - 8], ja[i - 7],
                      ja[i - 6],
                      ja[i - 5], ja[i - 4], ja[i - 3], ja[i - 2], ja[i - 1], ja[i], ja[i + 1],
                      file=open("output4.txt", "a+"), sep="")
                strix = open("output4.txt")
                striy = strix.readlines()
                string = striy[-1]
                res = [sub.replace('_', '') for sub in string]
                x = "".join(res)
                print(x[-17:], file=open("oligo.txt", "a+"))
            else:
                print(ja[i - 15], ja[i - 14], ja[i - 13], ja[i - 11], ja[i - 10], ja[i - 9], ja[i - 8],
                      ja[i - 7], ja[i - 6], ja[i - 5], ja[i - 4], ja[i - 3], ja[i - 2], ja[i - 1], ja[i], ja[i + 1],
                      file=open("oligo.txt", "a+"), sep="")
        if i == (ln - 2):
            if ja[i - 14] == "_" or ja[i - 13] == "_" or ja[i - 12] == "_" or \
                    ja[i - 11] == "_" or ja[i - 10] == "_" or ja[i - 9] == "_" or ja[i - 8] == "_" or \
                    ja[i - 7] == "_" or ja[i - 6] == "_" or ja[i - 5] == "_" or \
                    ja[i - 4] == "_" or ja[i - 3] == "_" or ja[i - 2] == "_" or ja[i - 1] == "_" or ja[i] == "_" or \
                    ja[i + 1] == "_" or ja[i + 2] == "_":
                print(ja[i - 66], ja[i - 65], ja[i - 64], ja[i - 63], ja[i - 62], ja[i - 61],
                      ja[i - 60],
                      ja[i - 59], ja[i - 58], ja[i - 57], ja[i - 56], ja[i - 55], ja[i - 54], ja[i - 53], ja[i - 52],
                      ja[i - 51],
                      ja[i - 50], ja[i - 49], ja[i - 48], ja[i - 47], ja[i - 46], ja[i - 45], ja[i - 44], ja[i - 43],
                      ja[i - 42],
                      ja[i - 41], ja[i - 40], ja[i - 39], ja[i - 38], ja[i - 37], ja[i - 36], ja[i - 35], ja[i - 34],
                      ja[i - 33],
                      ja[i - 32], ja[i - 31], ja[i - 30], ja[i - 29], ja[i - 28], ja[i - 27], ja[i - 26], ja[i - 25],
                      ja[i - 24],
                      ja[i - 23], ja[i - 22], ja[i - 21], ja[i - 20], ja[i - 19], ja[i - 18], ja[i - 17], ja[i - 16],
                      ja[i - 15],
                      ja[i - 14], ja[i - 13], ja[i - 12], ja[i - 11], ja[i - 10], ja[i - 9], ja[i - 8], ja[i - 7],
                      ja[i - 6],
                      ja[i - 5], ja[i - 4], ja[i - 3], ja[i - 2], ja[i - 1], ja[i], ja[i + 1], ja[i + 2],
                      file=open("output4.txt", "a+"), sep="")
                strix = open("output4.txt")
                striy = strix.readlines()
                string = striy[-1]
                res = [sub.replace('_', '') for sub in string]
                x = "".join(res)
                print(x[-17:], file=open("oligo.txt", "a+"))
            else:
                print(ja[i - 14], ja[i - 13], ja[i - 11], ja[i - 10], ja[i - 9], ja[i - 8],
                      ja[i - 7], ja[i - 6], ja[i - 5], ja[i - 4], ja[i - 3], ja[i - 2], ja[i - 1], ja[i], ja[i + 1],
                      ja[i + 2], file=open("oligo.txt", "a+"), sep="")
        if i == (ln - 3):
            if ja[i - 13] == "_" or ja[i - 12] == "_" or \
                    ja[i - 11] == "_" or ja[i - 10] == "_" or ja[i - 9] == "_" or ja[i - 8] == "_" or \
                    ja[i - 7] == "_" or ja[i - 6] == "_" or ja[i - 5] == "_" or \
                    ja[i - 4] == "_" or ja[i - 3] == "_" or ja[i - 2] == "_" or ja[i - 1] == "_" or ja[i] == "_" or \
                    ja[i + 1] == "_" or ja[i + 2] == "_" or ja[i + 3] == "_":
                print(ja[i - 65], ja[i - 64], ja[i - 63], ja[i - 62], ja[i - 61],
                      ja[i - 60],
                      ja[i - 59], ja[i - 58], ja[i - 57], ja[i - 56], ja[i - 55], ja[i - 54], ja[i - 53], ja[i - 52],
                      ja[i - 51],
                      ja[i - 50], ja[i - 49], ja[i - 48], ja[i - 47], ja[i - 46], ja[i - 45], ja[i - 44], ja[i - 43],
                      ja[i - 42],
                      ja[i - 41], ja[i - 40], ja[i - 39], ja[i - 38], ja[i - 37], ja[i - 36], ja[i - 35], ja[i - 34],
                      ja[i - 33],
                      ja[i - 32], ja[i - 31], ja[i - 30], ja[i - 29], ja[i - 28], ja[i - 27], ja[i - 26], ja[i - 25],
                      ja[i - 24],
                      ja[i - 23], ja[i - 22], ja[i - 21], ja[i - 20], ja[i - 19], ja[i - 18], ja[i - 17], ja[i - 16],
                      ja[i - 15],
                      ja[i - 14], ja[i - 13], ja[i - 12], ja[i - 11], ja[i - 10], ja[i - 9], ja[i - 8], ja[i - 7],
                      ja[i - 6],
                      ja[i - 5], ja[i - 4], ja[i - 3], ja[i - 2], ja[i - 1], ja[i], ja[i + 1], ja[i + 2], ja[i + 3],
                      file=open("output4.txt", "a+"), sep="")
                strix = open("output4.txt")
                striy = strix.readlines()
                string = striy[-1]
                res = [sub.replace('_', '') for sub in string]
                x = "".join(res)
                print(x[-17:], file=open("oligo.txt", "a+"))
            else:
                print(ja[i - 13], ja[i - 11], ja[i - 10], ja[i - 9], ja[i - 8],
                      ja[i - 7], ja[i - 6], ja[i - 5], ja[i - 4], ja[i - 3], ja[i - 2], ja[i - 1], ja[i], ja[i + 1],
                      ja[i + 2], ja[i + 3], file=open("oligo.txt", "a+"), sep="")
        if i == (ln - 4):
            if ja[i - 12] == "_" or ja[i - 11] == "_" or ja[i - 10] == "_" or ja[i - 9] == "_" or ja[i - 8] == "_" or \
                    ja[i - 7] == "_" or ja[i - 6] == "_" or ja[i - 5] == "_" or \
                    ja[i - 4] == "_" or ja[i - 3] == "_" or ja[i - 2] == "_" or ja[i - 1] == "_" or ja[i] == "_" or \
                    ja[i + 1] == "_" or ja[i + 2] == "_" or ja[i + 3] == "_" or ja[i + 4] == "_":
                print(ja[i - 64], ja[i - 63], ja[i - 62], ja[i - 61],
                      ja[i - 60],
                      ja[i - 59], ja[i - 58], ja[i - 57], ja[i - 56], ja[i - 55], ja[i - 54], ja[i - 53], ja[i - 52],
                      ja[i - 51],
                      ja[i - 50], ja[i - 49], ja[i - 48], ja[i - 47], ja[i - 46], ja[i - 45], ja[i - 44], ja[i - 43],
                      ja[i - 42],
                      ja[i - 41], ja[i - 40], ja[i - 39], ja[i - 38], ja[i - 37], ja[i - 36], ja[i - 35], ja[i - 34],
                      ja[i - 33],
                      ja[i - 32], ja[i - 31], ja[i - 30], ja[i - 29], ja[i - 28], ja[i - 27], ja[i - 26], ja[i - 25],
                      ja[i - 24],
                      ja[i - 23], ja[i - 22], ja[i - 21], ja[i - 20], ja[i - 19], ja[i - 18], ja[i - 17], ja[i - 16],
                      ja[i - 15],
                      ja[i - 14], ja[i - 13], ja[i - 12], ja[i - 11], ja[i - 10], ja[i - 9], ja[i - 8], ja[i - 7],
                      ja[i - 6],
                      ja[i - 5], ja[i - 4], ja[i - 3], ja[i - 2], ja[i - 1], ja[i], ja[i + 1], ja[i + 2], ja[i + 3],
                      ja[i + 4], file=open("output4.txt", "a+"), sep="")
                strix = open("output4.txt")
                striy = strix.readlines()
                string = striy[-1]
                res = [sub.replace('_', '') for sub in string]
                x = "".join(res)
                print(x[-17:], file=open("oligo.txt", "a+"))
            else:
                print(ja[i - 12], ja[i - 11], ja[i - 10], ja[i - 9], ja[i - 8],
                      ja[i - 7], ja[i - 6], ja[i - 5], ja[i - 4], ja[i - 3], ja[i - 2], ja[i - 1], ja[i], ja[i + 1],
                      ja[i + 2], ja[i + 3], ja[i + 4], file=open("oligo.txt", "a+"), sep="")
        if i == (ln - 5):
            if ja[i - 11] == "_" or ja[i - 10] == "_" or ja[i - 9] == "_" or ja[i - 8] == "_" or \
                    ja[i - 7] == "_" or ja[i - 6] == "_" or ja[i - 5] == "_" or \
                    ja[i - 4] == "_" or ja[i - 3] == "_" or ja[i - 2] == "_" or ja[i - 1] == "_" or ja[i] == "_" or \
                    ja[i + 1] == "_" or ja[i + 2] == "_" or ja[i + 3] == "_" or ja[i + 4] == "_" or ja[i + 5] == "_":
                print(ja[i - 63], ja[i - 62], ja[i - 61],
                      ja[i - 60],
                      ja[i - 59], ja[i - 58], ja[i - 57], ja[i - 56], ja[i - 55], ja[i - 54], ja[i - 53], ja[i - 52],
                      ja[i - 51],
                      ja[i - 50], ja[i - 49], ja[i - 48], ja[i - 47], ja[i - 46], ja[i - 45], ja[i - 44], ja[i - 43],
                      ja[i - 42],
                      ja[i - 41], ja[i - 40], ja[i - 39], ja[i - 38], ja[i - 37], ja[i - 36], ja[i - 35], ja[i - 34],
                      ja[i - 33],
                      ja[i - 32], ja[i - 31], ja[i - 30], ja[i - 29], ja[i - 28], ja[i - 27], ja[i - 26], ja[i - 25],
                      ja[i - 24],
                      ja[i - 23], ja[i - 22], ja[i - 21], ja[i - 20], ja[i - 19], ja[i - 18], ja[i - 17], ja[i - 16],
                      ja[i - 15],
                      ja[i - 14], ja[i - 13], ja[i - 12], ja[i - 11], ja[i - 10], ja[i - 9], ja[i - 8], ja[i - 7],
                      ja[i - 6],
                      ja[i - 5], ja[i - 4], ja[i - 3], ja[i - 2], ja[i - 1], ja[i], ja[i + 1], ja[i + 2], ja[i + 3],
                      ja[i + 4], ja[i + 5], file=open("output4.txt", "a+"), sep="")
                strix = open("output4.txt")
                striy = strix.readlines()
                string = striy[-1]
                res = [sub.replace('_', '') for sub in string]
                x = "".join(res)
                print(x[-17:], file=open("oligo.txt", "a+"))
            else:
                print(ja[i - 11], ja[i - 10], ja[i - 9], ja[i - 8],
                      ja[i - 7], ja[i - 6], ja[i - 5], ja[i - 4], ja[i - 3], ja[i - 2], ja[i - 1], ja[i], ja[i + 1],
                      ja[i + 2], ja[i + 3], ja[i + 4], ja[i + 5], file=open("oligo.txt", "a+"), sep="")
        if i == (ln - 6):
            if ja[i - 10] == "_" or ja[i - 9] == "_" or ja[i - 8] == "_" or \
                    ja[i - 7] == "_" or ja[i - 6] == "_" or ja[i - 5] == "_" or \
                    ja[i - 4] == "_" or ja[i - 3] == "_" or ja[i - 2] == "_" or ja[i - 1] == "_" or ja[i] == "_" or \
                    ja[i + 1] == "_" or ja[i + 2] == "_" or ja[i + 3] == "_" or ja[i + 4] == "_" or ja[i + 5] == "_" or \
                    ja[i + 6] == "_":
                print(ja[i - 62], ja[i - 61],
                      ja[i - 60],
                      ja[i - 59], ja[i - 58], ja[i - 57], ja[i - 56], ja[i - 55], ja[i - 54], ja[i - 53], ja[i - 52],
                      ja[i - 51],
                      ja[i - 50], ja[i - 49], ja[i - 48], ja[i - 47], ja[i - 46], ja[i - 45], ja[i - 44], ja[i - 43],
                      ja[i - 42],
                      ja[i - 41], ja[i - 40], ja[i - 39], ja[i - 38], ja[i - 37], ja[i - 36], ja[i - 35], ja[i - 34],
                      ja[i - 33],
                      ja[i - 32], ja[i - 31], ja[i - 30], ja[i - 29], ja[i - 28], ja[i - 27], ja[i - 26], ja[i - 25],
                      ja[i - 24],
                      ja[i - 23], ja[i - 22], ja[i - 21], ja[i - 20], ja[i - 19], ja[i - 18], ja[i - 17], ja[i - 16],
                      ja[i - 15],
                      ja[i - 14], ja[i - 13], ja[i - 12], ja[i - 11], ja[i - 10], ja[i - 9], ja[i - 8], ja[i - 7],
                      ja[i - 6],
                      ja[i - 5], ja[i - 4], ja[i - 3], ja[i - 2], ja[i - 1], ja[i], ja[i + 1], ja[i + 2], ja[i + 3],
                      ja[i + 4], ja[i + 5], ja[i + 6], file=open("output4.txt", "a+"), sep="")
                strix = open("output4.txt")
                striy = strix.readlines()
                string = striy[-1]
                res = [sub.replace('_', '') for sub in string]
                x = "".join(res)
                print(x[-17:], file=open("oligo.txt", "a+"))
            else:
                print(ja[i - 10], ja[i - 9], ja[i - 8],
                      ja[i - 7], ja[i - 6], ja[i - 5], ja[i - 4], ja[i - 3], ja[i - 2], ja[i - 1], ja[i], ja[i + 1],
                      ja[i + 2], ja[i + 3], ja[i + 4], ja[i + 5], ja[i + 6], file=open("oligo.txt", "a+"), sep="")
        if i == (ln - 7):
            if ja[i - 9] == "_" or ja[i - 8] == "_" or \
                    ja[i - 7] == "_" or ja[i - 6] == "_" or ja[i - 5] == "_" or \
                    ja[i - 4] == "_" or ja[i - 3] == "_" or ja[i - 2] == "_" or ja[i - 1] == "_" or ja[i] == "_" or \
                    ja[i + 1] == "_" or ja[i + 2] == "_" or ja[i + 3] == "_" or ja[i + 4] == "_" or ja[i + 5] == "_" or \
                    ja[i + 6] == "_" or ja[i + 7] == "_":
                print(ja[i - 61],
                      ja[i - 60],
                      ja[i - 59], ja[i - 58], ja[i - 57], ja[i - 56], ja[i - 55], ja[i - 54], ja[i - 53], ja[i - 52],
                      ja[i - 51],
                      ja[i - 50], ja[i - 49], ja[i - 48], ja[i - 47], ja[i - 46], ja[i - 45], ja[i - 44], ja[i - 43],
                      ja[i - 42],
                      ja[i - 41], ja[i - 40], ja[i - 39], ja[i - 38], ja[i - 37], ja[i - 36], ja[i - 35], ja[i - 34],
                      ja[i - 33],
                      ja[i - 32], ja[i - 31], ja[i - 30], ja[i - 29], ja[i - 28], ja[i - 27], ja[i - 26], ja[i - 25],
                      ja[i - 24],
                      ja[i - 23], ja[i - 22], ja[i - 21], ja[i - 20], ja[i - 19], ja[i - 18], ja[i - 17], ja[i - 16],
                      ja[i - 15],
                      ja[i - 14], ja[i - 13], ja[i - 12], ja[i - 11], ja[i - 10], ja[i - 9], ja[i - 8], ja[i - 7],
                      ja[i - 6],
                      ja[i - 5], ja[i - 4], ja[i - 3], ja[i - 2], ja[i - 1], ja[i], ja[i + 1], ja[i + 2], ja[i + 3],
                      ja[i + 4], ja[i + 5], ja[i + 6], ja[i + 7], file=open("output4.txt", "a+"), sep="")
                strix = open("output4.txt")
                striy = strix.readlines()
                string = striy[-1]
                res = [sub.replace('_', '') for sub in string]
                x = "".join(res)
                print(x[-17:], file=open("oligo.txt", "a+"))
            else:
                print(ja[i - 9], ja[i - 8],
                      ja[i - 7], ja[i - 6], ja[i - 5], ja[i - 4], ja[i - 3], ja[i - 2], ja[i - 1], ja[i], ja[i + 1],
                      ja[i + 2], ja[i + 3], ja[i + 4], ja[i + 5], ja[i + 6], ja[i + 7], file=open("oligo.txt", "a+"),
                      sep="")
        if i == (ln - 8):
            if ja[i - 8] == "_" or \
                    ja[i - 7] == "_" or ja[i - 6] == "_" or ja[i - 5] == "_" or \
                    ja[i - 4] == "_" or ja[i - 3] == "_" or ja[i - 2] == "_" or ja[i - 1] == "_" or ja[i] == "_" or \
                    ja[i + 1] == "_" or ja[i + 2] == "_" or ja[i + 3] == "_" or ja[i + 4] == "_" or ja[i + 5] == "_" or \
                    ja[i + 6] == "_" or ja[i + 7] == "_" or ja[i + 8] == "_":
                print(ja[i - 60],
                      ja[i - 59], ja[i - 58], ja[i - 57], ja[i - 56], ja[i - 55], ja[i - 54], ja[i - 53], ja[i - 52],
                      ja[i - 51],
                      ja[i - 50], ja[i - 49], ja[i - 48], ja[i - 47], ja[i - 46], ja[i - 45], ja[i - 44], ja[i - 43],
                      ja[i - 42],
                      ja[i - 41], ja[i - 40], ja[i - 39], ja[i - 38], ja[i - 37], ja[i - 36], ja[i - 35], ja[i - 34],
                      ja[i - 33],
                      ja[i - 32], ja[i - 31], ja[i - 30], ja[i - 29], ja[i - 28], ja[i - 27], ja[i - 26], ja[i - 25],
                      ja[i - 24],
                      ja[i - 23], ja[i - 22], ja[i - 21], ja[i - 20], ja[i - 19], ja[i - 18], ja[i - 17], ja[i - 16],
                      ja[i - 15],
                      ja[i - 14], ja[i - 13], ja[i - 12], ja[i - 11], ja[i - 10], ja[i - 9], ja[i - 8], ja[i - 7],
                      ja[i - 6],
                      ja[i - 5], ja[i - 4], ja[i - 3], ja[i - 2], ja[i - 1], ja[i], ja[i + 1], ja[i + 2], ja[i + 3],
                      ja[i + 4], ja[i + 5], ja[i + 6], ja[i + 7], ja[i + 8], file=open("output4.txt", "a+"), sep="")
                strix = open("output4.txt")
                striy = strix.readlines()
                string = striy[-1]
                res = [sub.replace('_', '') for sub in string]
                x = "".join(res)
                print(x[-17:], file=open("oligo.txt", "a+"))
            else:
                print(ja[i - 8],
                      ja[i - 7], ja[i - 6], ja[i - 5], ja[i - 4], ja[i - 3], ja[i - 2], ja[i - 1], ja[i], ja[i + 1],
                      ja[i + 2], ja[i + 3], ja[i + 4], ja[i + 5], ja[i + 6], ja[i + 7], ja[i + 8],
                      file=open("oligo.txt", "a+"), sep="")
        if i == (ln - 9):
            if ja[i - 7] == "_" or ja[i - 6] == "_" or ja[i - 5] == "_" or \
                    ja[i - 4] == "_" or ja[i - 3] == "_" or ja[i - 2] == "_" or ja[i - 1] == "_" or ja[i] == "_" or \
                    ja[i + 1] == "_" or ja[i + 2] == "_" or ja[i + 3] == "_" or ja[i + 4] == "_" or ja[i + 5] == "_" or \
                    ja[i + 6] == "_" or ja[i + 7] == "_" or ja[i + 8] == "_" or ja[i + 9] == "_":
                print(ja[i - 59], ja[i - 58], ja[i - 57], ja[i - 56], ja[i - 55], ja[i - 54], ja[i - 53], ja[i - 52],
                      ja[i - 51],
                      ja[i - 50], ja[i - 49], ja[i - 48], ja[i - 47], ja[i - 46], ja[i - 45], ja[i - 44], ja[i - 43],
                      ja[i - 42],
                      ja[i - 41], ja[i - 40], ja[i - 39], ja[i - 38], ja[i - 37], ja[i - 36], ja[i - 35], ja[i - 34],
                      ja[i - 33],
                      ja[i - 32], ja[i - 31], ja[i - 30], ja[i - 29], ja[i - 28], ja[i - 27], ja[i - 26], ja[i - 25],
                      ja[i - 24],
                      ja[i - 23], ja[i - 22], ja[i - 21], ja[i - 20], ja[i - 19], ja[i - 18], ja[i - 17], ja[i - 16],
                      ja[i - 15],
                      ja[i - 14], ja[i - 13], ja[i - 12], ja[i - 11], ja[i - 10], ja[i - 9], ja[i - 8], ja[i - 7],
                      ja[i - 6],
                      ja[i - 5], ja[i - 4], ja[i - 3], ja[i - 2], ja[i - 1], ja[i], ja[i + 1], ja[i + 2], ja[i + 3],
                      ja[i + 4], ja[i + 5], ja[i + 6], ja[i + 7], ja[i + 8], ja[i + 9], file=open("output4.txt", "a+"),
                      sep="")
                strix = open("output4.txt")
                striy = strix.readlines()
                string = striy[-1]
                res = [sub.replace('_', '') for sub in string]
                x = "".join(res)
                print(x[-17:], file=open("oligo.txt", "a+"))
            else:
                print(ja[i - 7], ja[i - 6], ja[i - 5], ja[i - 4], ja[i - 3], ja[i - 2], ja[i - 1], ja[i], ja[i + 1],
                      ja[i + 2], ja[i + 3], ja[i + 4], ja[i + 5], ja[i + 6], ja[i + 7], ja[i + 8], ja[i + 9],
                      file=open("oligo.txt", "a+"), sep="")
        if i == (ln - 10):
            if ja[i - 6] == "_" or ja[i - 5] == "_" or \
                    ja[i - 4] == "_" or ja[i - 3] == "_" or ja[i - 2] == "_" or ja[i - 1] == "_" or ja[i] == "_" or \
                    ja[i + 1] == "_" or ja[i + 2] == "_" or ja[i + 3] == "_" or ja[i + 4] == "_" or ja[i + 5] == "_" or \
                    ja[i + 6] == "_" or ja[i + 7] == "_" or ja[i + 8] == "_" or ja[i + 9] == "_" or ja[i + 10] == "_":
                print(ja[i - 58], ja[i - 57], ja[i - 56], ja[i - 55], ja[i - 54], ja[i - 53], ja[i - 52],
                      ja[i - 51],
                      ja[i - 50], ja[i - 49], ja[i - 48], ja[i - 47], ja[i - 46], ja[i - 45], ja[i - 44], ja[i - 43],
                      ja[i - 42],
                      ja[i - 41], ja[i - 40], ja[i - 39], ja[i - 38], ja[i - 37], ja[i - 36], ja[i - 35], ja[i - 34],
                      ja[i - 33],
                      ja[i - 32], ja[i - 31], ja[i - 30], ja[i - 29], ja[i - 28], ja[i - 27], ja[i - 26], ja[i - 25],
                      ja[i - 24],
                      ja[i - 23], ja[i - 22], ja[i - 21], ja[i - 20], ja[i - 19], ja[i - 18], ja[i - 17], ja[i - 16],
                      ja[i - 15],
                      ja[i - 14], ja[i - 13], ja[i - 12], ja[i - 11], ja[i - 10], ja[i - 9], ja[i - 8], ja[i - 7],
                      ja[i - 6],
                      ja[i - 5], ja[i - 4], ja[i - 3], ja[i - 2], ja[i - 1], ja[i], ja[i + 1], ja[i + 2], ja[i + 3],
                      ja[i + 4], ja[i + 5], ja[i + 6], ja[i + 7], ja[i + 8], ja[i + 9], ja[i + 10],
                      file=open("output4.txt", "a+"), sep="")
                strix = open("output4.txt")
                striy = strix.readlines()
                string = striy[-1]
                res = [sub.replace('_', '') for sub in string]
                x = "".join(res)
                print(x[-17:], file=open("oligo.txt", "a+"))
            else:
                print(ja[i - 6], ja[i - 5], ja[i - 4], ja[i - 3], ja[i - 2], ja[i - 1], ja[i], ja[i + 1],
                      ja[i + 2], ja[i + 3], ja[i + 4], ja[i + 5], ja[i + 6], ja[i + 7], ja[i + 8], ja[i + 9],
                      ja[i + 10],
                      file=open("oligo.txt", "a+"), sep="")
        if i == (ln - 11):
            if ja[i - 5] == "_" or \
                    ja[i - 4] == "_" or ja[i - 3] == "_" or ja[i - 2] == "_" or ja[i - 1] == "_" or ja[i] == "_" or \
                    ja[i + 1] == "_" or ja[i + 2] == "_" or ja[i + 3] == "_" or ja[i + 4] == "_" or ja[i + 5] == "_" or \
                    ja[i + 6] == "_" or ja[i + 7] == "_" or ja[i + 8] == "_" or ja[i + 9] == "_" or ja[i + 10] == "_" or \
                    ja[i + 11] == "_":
                print(ja[i - 57], ja[i - 56], ja[i - 55], ja[i - 54], ja[i - 53], ja[i - 52],
                      ja[i - 51],
                      ja[i - 50], ja[i - 49], ja[i - 48], ja[i - 47], ja[i - 46], ja[i - 45], ja[i - 44], ja[i - 43],
                      ja[i - 42],
                      ja[i - 41], ja[i - 40], ja[i - 39], ja[i - 38], ja[i - 37], ja[i - 36], ja[i - 35], ja[i - 34],
                      ja[i - 33],
                      ja[i - 32], ja[i - 31], ja[i - 30], ja[i - 29], ja[i - 28], ja[i - 27], ja[i - 26], ja[i - 25],
                      ja[i - 24],
                      ja[i - 23], ja[i - 22], ja[i - 21], ja[i - 20], ja[i - 19], ja[i - 18], ja[i - 17], ja[i - 16],
                      ja[i - 15],
                      ja[i - 14], ja[i - 13], ja[i - 12], ja[i - 11], ja[i - 10], ja[i - 9], ja[i - 8], ja[i - 7],
                      ja[i - 6],
                      ja[i - 5], ja[i - 4], ja[i - 3], ja[i - 2], ja[i - 1], ja[i], ja[i + 1], ja[i + 2], ja[i + 3],
                      ja[i + 4], ja[i + 5], ja[i + 6], ja[i + 7], ja[i + 8], ja[i + 9], ja[i + 10], ja[i + 11],
                      file=open("output4.txt", "a+"), sep="")
                strix = open("output4.txt")
                striy = strix.readlines()
                string = striy[-1]
                res = [sub.replace('_', '') for sub in string]
                x = "".join(res)
                print(x[-17:], file=open("oligo.txt", "a+"))
            else:
                print(ja[i - 5], ja[i - 4], ja[i - 3], ja[i - 2], ja[i - 1], ja[i], ja[i + 1],
                      ja[i + 2], ja[i + 3], ja[i + 4], ja[i + 5], ja[i + 6], ja[i + 7], ja[i + 8], ja[i + 9],
                      ja[i + 10], ja[i + 11], file=open("oligo.txt", "a+"), sep="")
        if i == (ln - 12):
            if ja[i - 4] == "_" or ja[i - 3] == "_" or ja[i - 2] == "_" or ja[i - 1] == "_" or ja[i] == "_" or \
                    ja[i + 1] == "_" or ja[i + 2] == "_" or ja[i + 3] == "_" or ja[i + 4] == "_" or ja[i + 5] == "_" or \
                    ja[i + 6] == "_" or ja[i + 7] == "_" or ja[i + 8] == "_" or ja[i + 9] == "_" or ja[i + 10] == "_" or \
                    ja[i + 11] == "_" or ja[i + 12] == "_":
                print(ja[i - 56], ja[i - 55], ja[i - 54], ja[i - 53], ja[i - 52],
                      ja[i - 51],
                      ja[i - 50], ja[i - 49], ja[i - 48], ja[i - 47], ja[i - 46], ja[i - 45], ja[i - 44], ja[i - 43],
                      ja[i - 42],
                      ja[i - 41], ja[i - 40], ja[i - 39], ja[i - 38], ja[i - 37], ja[i - 36], ja[i - 35], ja[i - 34],
                      ja[i - 33],
                      ja[i - 32], ja[i - 31], ja[i - 30], ja[i - 29], ja[i - 28], ja[i - 27], ja[i - 26], ja[i - 25],
                      ja[i - 24],
                      ja[i - 23], ja[i - 22], ja[i - 21], ja[i - 20], ja[i - 19], ja[i - 18], ja[i - 17], ja[i - 16],
                      ja[i - 15],
                      ja[i - 14], ja[i - 13], ja[i - 12], ja[i - 11], ja[i - 10], ja[i - 9], ja[i - 8], ja[i - 7],
                      ja[i - 6],
                      ja[i - 5], ja[i - 4], ja[i - 3], ja[i - 2], ja[i - 1], ja[i], ja[i + 1], ja[i + 2], ja[i + 3],
                      ja[i + 4], ja[i + 5], ja[i + 6], ja[i + 7], ja[i + 8], ja[i + 9], ja[i + 10], ja[i + 11],
                      ja[i + 12],
                      file=open("output4.txt", "a+"), sep="")
                strix = open("output4.txt")
                striy = strix.readlines()
                string = striy[-1]
                res = [sub.replace('_', '') for sub in string]
                x = "".join(res)
                print(x[-17:], file=open("oligo.txt", "a+"))
            else:
                print(ja[i - 4], ja[i - 3], ja[i - 2], ja[i - 1], ja[i], ja[i + 1],
                      ja[i + 2], ja[i + 3], ja[i + 4], ja[i + 5], ja[i + 6], ja[i + 7], ja[i + 8], ja[i + 9],
                      ja[i + 10], ja[i + 11], ja[i + 12], file=open("oligo.txt", "a+"), sep="")
        if i == (ln - 13):
            if ja[i - 3] == "_" or ja[i - 2] == "_" or ja[i - 1] == "_" or ja[i] == "_" or \
                    ja[i + 1] == "_" or ja[i + 2] == "_" or ja[i + 3] == "_" or ja[i + 4] == "_" or ja[i + 5] == "_" or \
                    ja[i + 6] == "_" or ja[i + 7] == "_" or ja[i + 8] == "_" or ja[i + 9] == "_" or ja[i + 10] == "_" or \
                    ja[i + 11] == "_" or ja[i + 12] == "_" or ja[i + 13] == "_":
                print(ja[i - 55], ja[i - 54], ja[i - 53], ja[i - 52],
                      ja[i - 51],
                      ja[i - 50], ja[i - 49], ja[i - 48], ja[i - 47], ja[i - 46], ja[i - 45], ja[i - 44], ja[i - 43],
                      ja[i - 42],
                      ja[i - 41], ja[i - 40], ja[i - 39], ja[i - 38], ja[i - 37], ja[i - 36], ja[i - 35], ja[i - 34],
                      ja[i - 33],
                      ja[i - 32], ja[i - 31], ja[i - 30], ja[i - 29], ja[i - 28], ja[i - 27], ja[i - 26], ja[i - 25],
                      ja[i - 24],
                      ja[i - 23], ja[i - 22], ja[i - 21], ja[i - 20], ja[i - 19], ja[i - 18], ja[i - 17], ja[i - 16],
                      ja[i - 15],
                      ja[i - 14], ja[i - 13], ja[i - 12], ja[i - 11], ja[i - 10], ja[i - 9], ja[i - 8], ja[i - 7],
                      ja[i - 6],
                      ja[i - 5], ja[i - 4], ja[i - 3], ja[i - 2], ja[i - 1], ja[i], ja[i + 1], ja[i + 2], ja[i + 3],
                      ja[i + 4], ja[i + 5], ja[i + 6], ja[i + 7], ja[i + 8], ja[i + 9], ja[i + 10], ja[i + 11],
                      ja[i + 12],
                      ja[i + 13], file=open("output4.txt", "a+"), sep="")
                strix = open("output4.txt")
                striy = strix.readlines()
                string = striy[-1]
                res = [sub.replace('_', '') for sub in string]
                x = "".join(res)
                print(x[-17:], file=open("oligo.txt", "a+"))
            else:
                print(ja[i - 3], ja[i - 2], ja[i - 1], ja[i], ja[i + 1],
                      ja[i + 2], ja[i + 3], ja[i + 4], ja[i + 5], ja[i + 6], ja[i + 7], ja[i + 8], ja[i + 9],
                      ja[i + 10], ja[i + 11], ja[i + 12], ja[i + 13], file=open("oligo.txt", "a+"), sep="")
        if i == (ln - 14):
            if ja[i - 2] == "_" or ja[i - 1] == "_" or ja[i] == "_" or \
                    ja[i + 1] == "_" or ja[i + 2] == "_" or ja[i + 3] == "_" or ja[i + 4] == "_" or ja[i + 5] == "_" or \
                    ja[i + 6] == "_" or ja[i + 7] == "_" or ja[i + 8] == "_" or ja[i + 9] == "_" or ja[i + 10] == "_" or \
                    ja[i + 11] == "_" or ja[i + 12] == "_" or ja[i + 13] == "_" or ja[i + 14] == "_":
                print(ja[i - 54], ja[i - 53], ja[i - 52],
                      ja[i - 51],
                      ja[i - 50], ja[i - 49], ja[i - 48], ja[i - 47], ja[i - 46], ja[i - 45], ja[i - 44], ja[i - 43],
                      ja[i - 42],
                      ja[i - 41], ja[i - 40], ja[i - 39], ja[i - 38], ja[i - 37], ja[i - 36], ja[i - 35], ja[i - 34],
                      ja[i - 33],
                      ja[i - 32], ja[i - 31], ja[i - 30], ja[i - 29], ja[i - 28], ja[i - 27], ja[i - 26], ja[i - 25],
                      ja[i - 24],
                      ja[i - 23], ja[i - 22], ja[i - 21], ja[i - 20], ja[i - 19], ja[i - 18], ja[i - 17], ja[i - 16],
                      ja[i - 15],
                      ja[i - 14], ja[i - 13], ja[i - 12], ja[i - 11], ja[i - 10], ja[i - 9], ja[i - 8], ja[i - 7],
                      ja[i - 6],
                      ja[i - 5], ja[i - 4], ja[i - 3], ja[i - 2], ja[i - 1], ja[i], ja[i + 1], ja[i + 2], ja[i + 3],
                      ja[i + 4], ja[i + 5], ja[i + 6], ja[i + 7], ja[i + 8], ja[i + 9], ja[i + 10], ja[i + 11],
                      ja[i + 12],
                      ja[i + 13], ja[i + 14], file=open("output4.txt", "a+"), sep="")
                strix = open("output4.txt")
                striy = strix.readlines()
                string = striy[-1]
                res = [sub.replace('_', '') for sub in string]
                x = "".join(res)
                print(x[-17:], file=open("oligo.txt", "a+"))
            else:
                print(ja[i - 2], ja[i - 1], ja[i], ja[i + 1],
                      ja[i + 2], ja[i + 3], ja[i + 4], ja[i + 5], ja[i + 6], ja[i + 7], ja[i + 8], ja[i + 9],
                      ja[i + 10], ja[i + 11], ja[i + 12], ja[i + 13], ja[i + 14], file=open("oligo.txt", "a+"), sep="")
        if i == (ln - 15):
            if ja[i - 1] == "_" or ja[i] == "_" or \
                    ja[i + 1] == "_" or ja[i + 2] == "_" or ja[i + 3] == "_" or ja[i + 4] == "_" or ja[i + 5] == "_" or \
                    ja[i + 6] == "_" or ja[i + 7] == "_" or ja[i + 8] == "_" or ja[i + 9] == "_" or ja[i + 10] == "_" or \
                    ja[i + 11] == "_" or ja[i + 12] == "_" or ja[i + 13] == "_" or ja[i + 14] == "_" or ja[
                i + 15] == "_":
                print(ja[i - 53], ja[i - 52],
                      ja[i - 51],
                      ja[i - 50], ja[i - 49], ja[i - 48], ja[i - 47], ja[i - 46], ja[i - 45], ja[i - 44], ja[i - 43],
                      ja[i - 42],
                      ja[i - 41], ja[i - 40], ja[i - 39], ja[i - 38], ja[i - 37], ja[i - 36], ja[i - 35], ja[i - 34],
                      ja[i - 33],
                      ja[i - 32], ja[i - 31], ja[i - 30], ja[i - 29], ja[i - 28], ja[i - 27], ja[i - 26], ja[i - 25],
                      ja[i - 24],
                      ja[i - 23], ja[i - 22], ja[i - 21], ja[i - 20], ja[i - 19], ja[i - 18], ja[i - 17], ja[i - 16],
                      ja[i - 15],
                      ja[i - 14], ja[i - 13], ja[i - 12], ja[i - 11], ja[i - 10], ja[i - 9], ja[i - 8], ja[i - 7],
                      ja[i - 6],
                      ja[i - 5], ja[i - 4], ja[i - 3], ja[i - 2], ja[i - 1], ja[i], ja[i + 1], ja[i + 2], ja[i + 3],
                      ja[i + 4], ja[i + 5], ja[i + 6], ja[i + 7], ja[i + 8], ja[i + 9], ja[i + 10], ja[i + 11],
                      ja[i + 12],
                      ja[i + 13], ja[i + 14], ja[i + 15], file=open("output4.txt", "a+"), sep="")
                strix = open("output4.txt")
                striy = strix.readlines()
                string = striy[-1]
                res = [sub.replace('_', '') for sub in string]
                x = "".join(res)
                print(x[-17:], file=open("oligo.txt", "a+"))
            else:
                print(ja[i - 1], ja[i], ja[i + 1],
                      ja[i + 2], ja[i + 3], ja[i + 4], ja[i + 5], ja[i + 6], ja[i + 7], ja[i + 8], ja[i + 9],
                      ja[i + 10], ja[i + 11], ja[i + 12], ja[i + 13], ja[i + 14], ja[i + 15],
                      file=open("oligo.txt", "a+"), sep="")
        if i == (ln - 16):
            if ja[i] == "_" or \
                    ja[i + 1] == "_" or ja[i + 2] == "_" or ja[i + 3] == "_" or ja[i + 4] == "_" or ja[i + 5] == "_" or \
                    ja[i + 6] == "_" or ja[i + 7] == "_" or ja[i + 8] == "_" or ja[i + 9] == "_" or ja[i + 10] == "_" or \
                    ja[i + 11] == "_" or ja[i + 12] == "_" or ja[i + 13] == "_" or ja[i + 14] == "_" \
                    or ja[i + 15] == "_" or ja[i + 16] == "_":
                print(ja[i - 52],
                      ja[i - 51],
                      ja[i - 50], ja[i - 49], ja[i - 48], ja[i - 47], ja[i - 46], ja[i - 45], ja[i - 44], ja[i - 43],
                      ja[i - 42],
                      ja[i - 41], ja[i - 40], ja[i - 39], ja[i - 38], ja[i - 37], ja[i - 36], ja[i - 35], ja[i - 34],
                      ja[i - 33],
                      ja[i - 32], ja[i - 31], ja[i - 30], ja[i - 29], ja[i - 28], ja[i - 27], ja[i - 26], ja[i - 25],
                      ja[i - 24],
                      ja[i - 23], ja[i - 22], ja[i - 21], ja[i - 20], ja[i - 19], ja[i - 18], ja[i - 17], ja[i - 16],
                      ja[i - 15],
                      ja[i - 14], ja[i - 13], ja[i - 12], ja[i - 11], ja[i - 10], ja[i - 9], ja[i - 8], ja[i - 7],
                      ja[i - 6],
                      ja[i - 5], ja[i - 4], ja[i - 3], ja[i - 2], ja[i - 1], ja[i], ja[i + 1], ja[i + 2], ja[i + 3],
                      ja[i + 4], ja[i + 5], ja[i + 6], ja[i + 7], ja[i + 8], ja[i + 9], ja[i + 10], ja[i + 11],
                      ja[i + 12],
                      ja[i + 13], ja[i + 14], ja[i + 15], ja[i + 16], file=open("output4.txt", "a+"), sep="")
                strix = open("output4.txt")
                striy = strix.readlines()
                string = striy[-1]
                res = [sub.replace('_', '') for sub in string]
                x = "".join(res)
                print(x[-17:], file=open("oligo.txt", "a+"))
            else:
                print(ja[i], ja[i + 1],
                      ja[i + 2], ja[i + 3], ja[i + 4], ja[i + 5], ja[i + 6], ja[i + 7], ja[i + 8], ja[i + 9],
                      ja[i + 10], ja[i + 11], ja[i + 12], ja[i + 13], ja[i + 14], ja[i + 15], ja[i + 16],
                      file=open("oligo.txt", "a+"), sep="")
        if i == (ln - 17):
            if ja[i - 16] == "_" or ja[i - 15] == "_" or ja[i - 14] == "_" or ja[i - 13] == "_" or ja[i - 12] == "_" or \
                    ja[i - 11] == "_" or ja[i - 10] == "_" or ja[i - 9] == "_" or ja[i - 8] == "_" or \
                    ja[i - 7] == "_" or ja[i - 6] == "_" or ja[i - 5] == "_" or ja[i - 4] == "_" or ja[i - 3] == "_" or \
                    ja[i - 2] == "_" or ja[i - 1] == "_" or ja[i] == "_":
                print(ja[i - 68],
                      ja[i - 67],
                      ja[i - 66], ja[i - 65], ja[i - 64], ja[i - 63], ja[i - 62], ja[i - 61], ja[i - 60], ja[i - 59],
                      ja[i - 58],
                      ja[i - 57], ja[i - 56], ja[i - 55], ja[i - 54], ja[i - 53], ja[i - 52], ja[i - 51], ja[i - 50],
                      ja[i - 49],
                      ja[i - 48], ja[i - 47], ja[i - 46], ja[i - 45], ja[i - 44], ja[i - 43], ja[i - 42], ja[i - 41],
                      ja[i - 40],
                      ja[i - 39], ja[i - 38], ja[i - 37], ja[i - 36], ja[i - 35], ja[i - 34], ja[i - 33], ja[i - 32],
                      ja[i - 31],
                      ja[i - 30], ja[i - 29], ja[i - 28], ja[i - 27], ja[i - 26], ja[i - 25], ja[i - 24], ja[i - 23],
                      ja[i - 22],
                      ja[i - 21], ja[i - 20], ja[i - 19], ja[i - 18], ja[i - 17], ja[i - 16], ja[i - 15], ja[i - 14],
                      ja[i - 13],
                      ja[i - 12], ja[i - 11], ja[i - 10], ja[i - 9], ja[i - 8], ja[i - 7], ja[i - 6], ja[i - 5],
                      ja[i - 4],
                      ja[i - 3], ja[i - 2], ja[i - 1], ja[i], file=open("output4.txt", "a+"), sep="")
                strix = open("output4.txt")
                striy = strix.readlines()
                string = striy[-1]
                res = [sub.replace('_', '') for sub in string]
                x = "".join(res)
                print(x[-17:], file=open("oligo.txt", "a+"))
            else:
                print(ja[i - 16], ja[i - 15], ja[i - 14], ja[i - 13], ja[i - 12], ja[i - 11], ja[i - 10], ja[i - 9],
                      ja[i - 8], ja[i - 7], ja[i - 6], ja[i - 5], ja[i - 4], ja[i - 3], ja[i - 2], ja[i - 1], ja[i],
                      file=open("oligo.txt", "a+"), sep="")
        if i == (ln - 18):
            if ja[i - 15] == "_" or ja[i - 14] == "_" or ja[i - 13] == "_" or ja[i - 12] == "_" or \
                    ja[i - 11] == "_" or ja[i - 10] == "_" or ja[i - 9] == "_" or ja[i - 8] == "_" or \
                    ja[i - 7] == "_" or ja[i - 6] == "_" or ja[i - 5] == "_" or ja[i - 4] == "_" or ja[i - 3] == "_" or \
                    ja[i - 2] == "_" or ja[i - 1] == "_" or ja[i] == "_" or ja[i + 1] == "_":
                print(ja[i - 67],
                      ja[i - 66], ja[i - 65], ja[i - 64], ja[i - 63], ja[i - 62], ja[i - 61], ja[i - 60], ja[i - 59],
                      ja[i - 58],
                      ja[i - 57], ja[i - 56], ja[i - 55], ja[i - 54], ja[i - 53], ja[i - 52], ja[i - 51], ja[i - 50],
                      ja[i - 49],
                      ja[i - 48], ja[i - 47], ja[i - 46], ja[i - 45], ja[i - 44], ja[i - 43], ja[i - 42], ja[i - 41],
                      ja[i - 40],
                      ja[i - 39], ja[i - 38], ja[i - 37], ja[i - 36], ja[i - 35], ja[i - 34], ja[i - 33], ja[i - 32],
                      ja[i - 31],
                      ja[i - 30], ja[i - 29], ja[i - 28], ja[i - 27], ja[i - 26], ja[i - 25], ja[i - 24], ja[i - 23],
                      ja[i - 22],
                      ja[i - 21], ja[i - 20], ja[i - 19], ja[i - 18], ja[i - 17], ja[i - 16], ja[i - 15], ja[i - 14],
                      ja[i - 13],
                      ja[i - 12], ja[i - 11], ja[i - 10], ja[i - 9], ja[i - 8], ja[i - 7], ja[i - 6], ja[i - 5],
                      ja[i - 4],
                      ja[i - 3], ja[i - 2], ja[i - 1], ja[i], ja[i + 1], file=open("output4.txt", "a+"), sep="")
                strix = open("output4.txt")
                striy = strix.readlines()
                string = striy[-1]
                res = [sub.replace('_', '') for sub in string]
                x = "".join(res)
                print(x[-17:], file=open("oligo.txt", "a+"))
            else:
                print(ja[i - 15], ja[i - 14], ja[i - 13], ja[i - 12], ja[i - 11], ja[i - 10], ja[i - 9],
                      ja[i - 8], ja[i - 7], ja[i - 6], ja[i - 5], ja[i - 4], ja[i - 3], ja[i - 2], ja[i - 1], ja[i],
                      ja[i + 1], file=open("oligo.txt", "a+"), sep="")
        if i == (ln - 19):
            if ja[i - 14] == "_" or ja[i - 13] == "_" or ja[i - 12] == "_" or \
                    ja[i - 11] == "_" or ja[i - 10] == "_" or ja[i - 9] == "_" or ja[i - 8] == "_" or \
                    ja[i - 7] == "_" or ja[i - 6] == "_" or ja[i - 5] == "_" or ja[i - 4] == "_" or ja[i - 3] == "_" or \
                    ja[i - 2] == "_" or ja[i - 1] == "_" or ja[i] == "_" or ja[i + 1] == "_" or ja[i + 2] == "_":
                print(ja[i - 66], ja[i - 65], ja[i - 64], ja[i - 63], ja[i - 62], ja[i - 61], ja[i - 60], ja[i - 59],
                      ja[i - 58],
                      ja[i - 57], ja[i - 56], ja[i - 55], ja[i - 54], ja[i - 53], ja[i - 52], ja[i - 51], ja[i - 50],
                      ja[i - 49],
                      ja[i - 48], ja[i - 47], ja[i - 46], ja[i - 45], ja[i - 44], ja[i - 43], ja[i - 42], ja[i - 41],
                      ja[i - 40],
                      ja[i - 39], ja[i - 38], ja[i - 37], ja[i - 36], ja[i - 35], ja[i - 34], ja[i - 33], ja[i - 32],
                      ja[i - 31],
                      ja[i - 30], ja[i - 29], ja[i - 28], ja[i - 27], ja[i - 26], ja[i - 25], ja[i - 24], ja[i - 23],
                      ja[i - 22],
                      ja[i - 21], ja[i - 20], ja[i - 19], ja[i - 18], ja[i - 17], ja[i - 16], ja[i - 15], ja[i - 14],
                      ja[i - 13],
                      ja[i - 12], ja[i - 11], ja[i - 10], ja[i - 9], ja[i - 8], ja[i - 7], ja[i - 6], ja[i - 5],
                      ja[i - 4],
                      ja[i - 3], ja[i - 2], ja[i - 1], ja[i], ja[i + 1], ja[i + 2], file=open("output4.txt", "a+"),
                      sep="")
                strix = open("output4.txt")
                striy = strix.readlines()
                string = striy[-1]
                res = [sub.replace('_', '') for sub in string]
                x = "".join(res)
                print(x[-17:], file=open("oligo.txt", "a+"))
            else:
                print(ja[i - 14], ja[i - 13], ja[i - 12], ja[i - 11], ja[i - 10], ja[i - 9],
                      ja[i - 8], ja[i - 7], ja[i - 6], ja[i - 5], ja[i - 4], ja[i - 3], ja[i - 2], ja[i - 1], ja[i],
                      ja[i + 1], ja[i + 2], file=open("oligo.txt", "a+"), sep="")
        if i == (ln - 20):
            if ja[i - 13] == "_" or ja[i - 12] == "_" or \
                    ja[i - 11] == "_" or ja[i - 10] == "_" or ja[i - 9] == "_" or ja[i - 8] == "_" or \
                    ja[i - 7] == "_" or ja[i - 6] == "_" or ja[i - 5] == "_" or ja[i - 4] == "_" or ja[i - 3] == "_" or \
                    ja[i - 2] == "_" or ja[i - 1] == "_" or ja[i] == "_" or ja[i + 1] == "_" or ja[i + 2] == "_" or \
                    ja[i + 3] == "_":
                print(ja[i - 65], ja[i - 64], ja[i - 63], ja[i - 62], ja[i - 61], ja[i - 60], ja[i - 59],
                      ja[i - 58],
                      ja[i - 57], ja[i - 56], ja[i - 55], ja[i - 54], ja[i - 53], ja[i - 52], ja[i - 51], ja[i - 50],
                      ja[i - 49],
                      ja[i - 48], ja[i - 47], ja[i - 46], ja[i - 45], ja[i - 44], ja[i - 43], ja[i - 42], ja[i - 41],
                      ja[i - 40],
                      ja[i - 39], ja[i - 38], ja[i - 37], ja[i - 36], ja[i - 35], ja[i - 34], ja[i - 33], ja[i - 32],
                      ja[i - 31],
                      ja[i - 30], ja[i - 29], ja[i - 28], ja[i - 27], ja[i - 26], ja[i - 25], ja[i - 24], ja[i - 23],
                      ja[i - 22],
                      ja[i - 21], ja[i - 20], ja[i - 19], ja[i - 18], ja[i - 17], ja[i - 16], ja[i - 15], ja[i - 14],
                      ja[i - 13],
                      ja[i - 12], ja[i - 11], ja[i - 10], ja[i - 9], ja[i - 8], ja[i - 7], ja[i - 6], ja[i - 5],
                      ja[i - 4],
                      ja[i - 3], ja[i - 2], ja[i - 1], ja[i], ja[i + 1], ja[i + 2], ja[i + 3],
                      file=open("output4.txt", "a+"), sep="")
                strix = open("output4.txt")
                striy = strix.readlines()
                string = striy[-1]
                res = [sub.replace('_', '') for sub in string]
                x = "".join(res)
                print(x[-17:], file=open("oligo.txt", "a+"))
            else:
                print(ja[i - 13], ja[i - 12], ja[i - 11], ja[i - 10], ja[i - 9],
                      ja[i - 8], ja[i - 7], ja[i - 6], ja[i - 5], ja[i - 4], ja[i - 3], ja[i - 2], ja[i - 1], ja[i],
                      ja[i + 1], ja[i + 2], ja[i + 3], file=open("oligo.txt", "a+"), sep="")
        if i == (ln - 21):
            if ja[i - 12] == "_" or \
                    ja[i - 11] == "_" or ja[i - 10] == "_" or ja[i - 9] == "_" or ja[i - 8] == "_" or \
                    ja[i - 7] == "_" or ja[i - 6] == "_" or ja[i - 5] == "_" or ja[i - 4] == "_" or ja[i - 3] == "_" or \
                    ja[i - 2] == "_" or ja[i - 1] == "_" or ja[i] == "_" or ja[i + 1] == "_" or ja[i + 2] == "_" or \
                    ja[i + 3] == "_" or ja[i + 4] == "_":
                print(ja[i - 64], ja[i - 63], ja[i - 62], ja[i - 61], ja[i - 60], ja[i - 59],
                      ja[i - 58],
                      ja[i - 57], ja[i - 56], ja[i - 55], ja[i - 54], ja[i - 53], ja[i - 52], ja[i - 51], ja[i - 50],
                      ja[i - 49],
                      ja[i - 48], ja[i - 47], ja[i - 46], ja[i - 45], ja[i - 44], ja[i - 43], ja[i - 42], ja[i - 41],
                      ja[i - 40],
                      ja[i - 39], ja[i - 38], ja[i - 37], ja[i - 36], ja[i - 35], ja[i - 34], ja[i - 33], ja[i - 32],
                      ja[i - 31],
                      ja[i - 30], ja[i - 29], ja[i - 28], ja[i - 27], ja[i - 26], ja[i - 25], ja[i - 24], ja[i - 23],
                      ja[i - 22],
                      ja[i - 21], ja[i - 20], ja[i - 19], ja[i - 18], ja[i - 17], ja[i - 16], ja[i - 15], ja[i - 14],
                      ja[i - 13],
                      ja[i - 12], ja[i - 11], ja[i - 10], ja[i - 9], ja[i - 8], ja[i - 7], ja[i - 6], ja[i - 5],
                      ja[i - 4],
                      ja[i - 3], ja[i - 2], ja[i - 1], ja[i], ja[i + 1], ja[i + 2], ja[i + 3], ja[i + 4],
                      file=open("output4.txt", "a+"), sep="")
                strix = open("output4.txt")
                striy = strix.readlines()
                string = striy[-1]
                res = [sub.replace('_', '') for sub in string]
                x = "".join(res)
                print(x[-17:], file=open("oligo.txt", "a+"))
            else:
                print(ja[i - 12], ja[i - 11], ja[i - 10], ja[i - 9],
                      ja[i - 8], ja[i - 7], ja[i - 6], ja[i - 5], ja[i - 4], ja[i - 3], ja[i - 2], ja[i - 1], ja[i],
                      ja[i + 1], ja[i + 2], ja[i + 3], ja[i + 4], file=open("oligo.txt", "a+"), sep="")
        if i == (ln - 22):
            if ja[i - 11] == "_" or ja[i - 10] == "_" or ja[i - 9] == "_" or ja[i - 8] == "_" or \
                    ja[i - 7] == "_" or ja[i - 6] == "_" or ja[i - 5] == "_" or ja[i - 4] == "_" or ja[i - 3] == "_" or \
                    ja[i - 2] == "_" or ja[i - 1] == "_" or ja[i] == "_" or ja[i + 1] == "_" or ja[i + 2] == "_" or \
                    ja[i + 3] == "_" or ja[i + 4] == "_" or ja[i + 5] == "_":
                print(ja[i - 63], ja[i - 62], ja[i - 61], ja[i - 60], ja[i - 59],
                      ja[i - 58],
                      ja[i - 57], ja[i - 56], ja[i - 55], ja[i - 54], ja[i - 53], ja[i - 52], ja[i - 51], ja[i - 50],
                      ja[i - 49],
                      ja[i - 48], ja[i - 47], ja[i - 46], ja[i - 45], ja[i - 44], ja[i - 43], ja[i - 42], ja[i - 41],
                      ja[i - 40],
                      ja[i - 39], ja[i - 38], ja[i - 37], ja[i - 36], ja[i - 35], ja[i - 34], ja[i - 33], ja[i - 32],
                      ja[i - 31],
                      ja[i - 30], ja[i - 29], ja[i - 28], ja[i - 27], ja[i - 26], ja[i - 25], ja[i - 24], ja[i - 23],
                      ja[i - 22],
                      ja[i - 21], ja[i - 20], ja[i - 19], ja[i - 18], ja[i - 17], ja[i - 16], ja[i - 15], ja[i - 14],
                      ja[i - 13],
                      ja[i - 12], ja[i - 11], ja[i - 10], ja[i - 9], ja[i - 8], ja[i - 7], ja[i - 6], ja[i - 5],
                      ja[i - 4],
                      ja[i - 3], ja[i - 2], ja[i - 1], ja[i], ja[i + 1], ja[i + 2], ja[i + 3], ja[i + 4], ja[i + 5],
                      file=open("output4.txt", "a+"), sep="")
                strix = open("output4.txt")
                striy = strix.readlines()
                string = striy[-1]
                res = [sub.replace('_', '') for sub in string]
                x = "".join(res)
                print(x[-17:], file=open("oligo.txt", "a+"))
            else:
                print(ja[i - 11], ja[i - 10], ja[i - 9],
                      ja[i - 8], ja[i - 7], ja[i - 6], ja[i - 5], ja[i - 4], ja[i - 3], ja[i - 2], ja[i - 1], ja[i],
                      ja[i + 1], ja[i + 2], ja[i + 3], ja[i + 4], ja[i + 5], file=open("oligo.txt", "a+"), sep="")
        if i == (ln - 23):
            if ja[i - 10] == "_" or ja[i - 9] == "_" or ja[i - 8] == "_" or \
                    ja[i - 7] == "_" or ja[i - 6] == "_" or ja[i - 5] == "_" or ja[i - 4] == "_" or ja[i - 3] == "_" or \
                    ja[i - 2] == "_" or ja[i - 1] == "_" or ja[i] == "_" or ja[i + 1] == "_" or ja[i + 2] == "_" or \
                    ja[i + 3] == "_" or ja[i + 4] == "_" or ja[i + 5] == "_" or ja[i + 6] == "_":
                print(ja[i - 62], ja[i - 61], ja[i - 60], ja[i - 59],
                      ja[i - 58],
                      ja[i - 57], ja[i - 56], ja[i - 55], ja[i - 54], ja[i - 53], ja[i - 52], ja[i - 51], ja[i - 50],
                      ja[i - 49],
                      ja[i - 48], ja[i - 47], ja[i - 46], ja[i - 45], ja[i - 44], ja[i - 43], ja[i - 42], ja[i - 41],
                      ja[i - 40],
                      ja[i - 39], ja[i - 38], ja[i - 37], ja[i - 36], ja[i - 35], ja[i - 34], ja[i - 33], ja[i - 32],
                      ja[i - 31],
                      ja[i - 30], ja[i - 29], ja[i - 28], ja[i - 27], ja[i - 26], ja[i - 25], ja[i - 24], ja[i - 23],
                      ja[i - 22],
                      ja[i - 21], ja[i - 20], ja[i - 19], ja[i - 18], ja[i - 17], ja[i - 16], ja[i - 15], ja[i - 14],
                      ja[i - 13],
                      ja[i - 12], ja[i - 11], ja[i - 10], ja[i - 9], ja[i - 8], ja[i - 7], ja[i - 6], ja[i - 5],
                      ja[i - 4],
                      ja[i - 3], ja[i - 2], ja[i - 1], ja[i], ja[i + 1], ja[i + 2], ja[i + 3], ja[i + 4], ja[i + 5],
                      ja[i + 6],
                      file=open("output4.txt", "a+"), sep="")
                strix = open("output4.txt")
                striy = strix.readlines()
                string = striy[-1]
                res = [sub.replace('_', '') for sub in string]
                x = "".join(res)
                print(x[-17:], file=open("oligo.txt", "a+"))
            else:
                print(ja[i - 10], ja[i - 9],
                      ja[i - 8], ja[i - 7], ja[i - 6], ja[i - 5], ja[i - 4], ja[i - 3], ja[i - 2], ja[i - 1], ja[i],
                      ja[i + 1], ja[i + 2], ja[i + 3], ja[i + 4], ja[i + 5], ja[i + 6], file=open("oligo.txt", "a+"),
                      sep="")
        if i == (ln - 24):
            if ja[i - 9] == "_" or ja[i - 8] == "_" or \
                    ja[i - 7] == "_" or ja[i - 6] == "_" or ja[i - 5] == "_" or ja[i - 4] == "_" or ja[i - 3] == "_" or \
                    ja[i - 2] == "_" or ja[i - 1] == "_" or ja[i] == "_" or ja[i + 1] == "_" or ja[i + 2] == "_" or \
                    ja[i + 3] == "_" or ja[i + 4] == "_" or ja[i + 5] == "_" or ja[i + 6] == "_" or ja[i + 7] == "_":
                print(ja[i - 61], ja[i - 60], ja[i - 59],
                      ja[i - 58],
                      ja[i - 57], ja[i - 56], ja[i - 55], ja[i - 54], ja[i - 53], ja[i - 52], ja[i - 51], ja[i - 50],
                      ja[i - 49],
                      ja[i - 48], ja[i - 47], ja[i - 46], ja[i - 45], ja[i - 44], ja[i - 43], ja[i - 42], ja[i - 41],
                      ja[i - 40],
                      ja[i - 39], ja[i - 38], ja[i - 37], ja[i - 36], ja[i - 35], ja[i - 34], ja[i - 33], ja[i - 32],
                      ja[i - 31],
                      ja[i - 30], ja[i - 29], ja[i - 28], ja[i - 27], ja[i - 26], ja[i - 25], ja[i - 24], ja[i - 23],
                      ja[i - 22],
                      ja[i - 21], ja[i - 20], ja[i - 19], ja[i - 18], ja[i - 17], ja[i - 16], ja[i - 15], ja[i - 14],
                      ja[i - 13],
                      ja[i - 12], ja[i - 11], ja[i - 10], ja[i - 9], ja[i - 8], ja[i - 7], ja[i - 6], ja[i - 5],
                      ja[i - 4],
                      ja[i - 3], ja[i - 2], ja[i - 1], ja[i], ja[i + 1], ja[i + 2], ja[i + 3], ja[i + 4],
                      ja[i + 5], ja[i + 6], ja[i + 7], file=open("output4.txt", "a+"), sep="")
                strix = open("output4.txt")
                striy = strix.readlines()
                string = striy[-1]
                res = [sub.replace('_', '') for sub in string]
                x = "".join(res)
                print(x[-17:], file=open("oligo.txt", "a+"))
            else:
                print(ja[i - 9], ja[i - 8], ja[i - 7], ja[i - 6], ja[i - 5], ja[i - 4], ja[i - 3], ja[i - 2],
                      ja[i - 1], ja[i], ja[i + 1], ja[i + 2], ja[i + 3], ja[i + 4], ja[i + 5], ja[i + 6],
                      ja[i + 7], file=open("oligo.txt", "a+"), sep="")
        if i == (ln - 25):
            if ja[i - 8] == "_" or \
                    ja[i - 7] == "_" or ja[i - 6] == "_" or ja[i - 5] == "_" or ja[i - 4] == "_" or ja[i - 3] == "_" or \
                    ja[i - 2] == "_" or ja[i - 1] == "_" or ja[i] == "_" or ja[i + 1] == "_" or ja[i + 2] == "_" or \
                    ja[i + 3] == "_" or ja[i + 4] == "_" or ja[i + 5] == "_" or ja[i + 6] == "_" or ja[i + 7] == "_" or \
                    ja[i + 8] == "_":
                print(ja[i - 60], ja[i - 59],
                      ja[i - 58],
                      ja[i - 57], ja[i - 56], ja[i - 55], ja[i - 54], ja[i - 53], ja[i - 52], ja[i - 51], ja[i - 50],
                      ja[i - 49],
                      ja[i - 48], ja[i - 47], ja[i - 46], ja[i - 45], ja[i - 44], ja[i - 43], ja[i - 42], ja[i - 41],
                      ja[i - 40],
                      ja[i - 39], ja[i - 38], ja[i - 37], ja[i - 36], ja[i - 35], ja[i - 34], ja[i - 33], ja[i - 32],
                      ja[i - 31],
                      ja[i - 30], ja[i - 29], ja[i - 28], ja[i - 27], ja[i - 26], ja[i - 25], ja[i - 24], ja[i - 23],
                      ja[i - 22],
                      ja[i - 21], ja[i - 20], ja[i - 19], ja[i - 18], ja[i - 17], ja[i - 16], ja[i - 15], ja[i - 14],
                      ja[i - 13],
                      ja[i - 12], ja[i - 11], ja[i - 10], ja[i - 9], ja[i - 8], ja[i - 7], ja[i - 6], ja[i - 5],
                      ja[i - 4],
                      ja[i - 3], ja[i - 2], ja[i - 1], ja[i], ja[i + 1], ja[i + 2], ja[i + 3], ja[i + 4],
                      ja[i + 5], ja[i + 6], ja[i + 7], ja[i + 8], file=open("output4.txt", "a+"), sep="")
                strix = open("output4.txt")
                striy = strix.readlines()
                string = striy[-1]
                res = [sub.replace('_', '') for sub in string]
                x = "".join(res)
                print(x[-17:], file=open("oligo.txt", "a+"))
            else:
                print(ja[i - 8], ja[i - 7], ja[i - 6], ja[i - 5], ja[i - 4], ja[i - 3], ja[i - 2],
                      ja[i - 1], ja[i], ja[i + 1], ja[i + 2], ja[i + 3], ja[i + 4], ja[i + 5], ja[i + 6],
                      ja[i + 7], ja[i + 8], file=open("oligo.txt", "a+"), sep="")
        if i == (ln - 26):
            if ja[i - 7] == "_" or ja[i - 6] == "_" or ja[i - 5] == "_" or ja[i - 4] == "_" or ja[i - 3] == "_" or \
                    ja[i - 2] == "_" or ja[i - 1] == "_" or ja[i] == "_" or ja[i + 1] == "_" or ja[i + 2] == "_" or \
                    ja[i + 3] == "_" or ja[i + 4] == "_" or ja[i + 5] == "_" or ja[i + 6] == "_" or ja[i + 7] == "_" or \
                    ja[i + 8] == "_" or ja[i + 9] == "_":
                print(ja[i - 59],
                      ja[i - 58],
                      ja[i - 57], ja[i - 56], ja[i - 55], ja[i - 54], ja[i - 53], ja[i - 52], ja[i - 51], ja[i - 50],
                      ja[i - 49],
                      ja[i - 48], ja[i - 47], ja[i - 46], ja[i - 45], ja[i - 44], ja[i - 43], ja[i - 42], ja[i - 41],
                      ja[i - 40],
                      ja[i - 39], ja[i - 38], ja[i - 37], ja[i - 36], ja[i - 35], ja[i - 34], ja[i - 33], ja[i - 32],
                      ja[i - 31],
                      ja[i - 30], ja[i - 29], ja[i - 28], ja[i - 27], ja[i - 26], ja[i - 25], ja[i - 24], ja[i - 23],
                      ja[i - 22],
                      ja[i - 21], ja[i - 20], ja[i - 19], ja[i - 18], ja[i - 17], ja[i - 16], ja[i - 15], ja[i - 14],
                      ja[i - 13],
                      ja[i - 12], ja[i - 11], ja[i - 10], ja[i - 9], ja[i - 8], ja[i - 7], ja[i - 6], ja[i - 5],
                      ja[i - 4],
                      ja[i - 3], ja[i - 2], ja[i - 1], ja[i], ja[i + 1], ja[i + 2], ja[i + 3], ja[i + 4],
                      ja[i + 5], ja[i + 6], ja[i + 7], ja[i + 8], ja[i + 9], file=open("output4.txt", "a+"), sep="")
                strix = open("output4.txt")
                striy = strix.readlines()
                string = striy[-1]
                res = [sub.replace('_', '') for sub in string]
                x = "".join(res)
                print(x[-17:], file=open("oligo.txt", "a+"))
            else:
                print(ja[i - 7], ja[i - 6], ja[i - 5], ja[i - 4], ja[i - 3], ja[i - 2],
                      ja[i - 1], ja[i], ja[i + 1], ja[i + 2], ja[i + 3], ja[i + 4], ja[i + 5], ja[i + 6],
                      ja[i + 7], ja[i + 8], ja[i + 9], file=open("oligo.txt", "a+"), sep="")
        if i == (ln - 27):
            if ja[i - 6] == "_" or ja[i - 5] == "_" or ja[i - 4] == "_" or ja[i - 3] == "_" or \
                    ja[i - 2] == "_" or ja[i - 1] == "_" or ja[i] == "_" or ja[i + 1] == "_" or ja[i + 2] == "_" or \
                    ja[i + 3] == "_" or ja[i + 4] == "_" or ja[i + 5] == "_" or ja[i + 6] == "_" or ja[i + 7] == "_" or \
                    ja[i + 8] == "_" or ja[i + 9] == "_" or ja[i + 10] == "_":
                print(ja[i - 58],
                      ja[i - 57], ja[i - 56], ja[i - 55], ja[i - 54], ja[i - 53], ja[i - 52], ja[i - 51], ja[i - 50],
                      ja[i - 49],
                      ja[i - 48], ja[i - 47], ja[i - 46], ja[i - 45], ja[i - 44], ja[i - 43], ja[i - 42], ja[i - 41],
                      ja[i - 40],
                      ja[i - 39], ja[i - 38], ja[i - 37], ja[i - 36], ja[i - 35], ja[i - 34], ja[i - 33], ja[i - 32],
                      ja[i - 31],
                      ja[i - 30], ja[i - 29], ja[i - 28], ja[i - 27], ja[i - 26], ja[i - 25], ja[i - 24], ja[i - 23],
                      ja[i - 22],
                      ja[i - 21], ja[i - 20], ja[i - 19], ja[i - 18], ja[i - 17], ja[i - 16], ja[i - 15], ja[i - 14],
                      ja[i - 13],
                      ja[i - 12], ja[i - 11], ja[i - 10], ja[i - 9], ja[i - 8], ja[i - 7], ja[i - 6], ja[i - 5],
                      ja[i - 4],
                      ja[i - 3], ja[i - 2], ja[i - 1], ja[i], ja[i + 1], ja[i + 2], ja[i + 3], ja[i + 4],
                      ja[i + 5], ja[i + 6], ja[i + 7], ja[i + 8], ja[i + 9], ja[i + 10],
                      file=open("output4.txt", "a+"), sep="")
                strix = open("output4.txt")
                striy = strix.readlines()
                string = striy[-1]
                res = [sub.replace('_', '') for sub in string]
                x = "".join(res)
                print(x[-17:], file=open("oligo.txt", "a+"))
            else:
                print(ja[i - 6], ja[i - 5], ja[i - 4], ja[i - 3], ja[i - 2],
                      ja[i - 1], ja[i], ja[i + 1], ja[i + 2], ja[i + 3], ja[i + 4], ja[i + 5], ja[i + 6],
                      ja[i + 7], ja[i + 8], ja[i + 9], ja[i + 10], file=open("oligo.txt", "a+"), sep="")
        if i == (ln - 28):
            if ja[i - 5] == "_" or ja[i - 4] == "_" or ja[i - 3] == "_" or \
                    ja[i - 2] == "_" or ja[i - 1] == "_" or ja[i] == "_" or ja[i + 1] == "_" or ja[i + 2] == "_" or \
                    ja[i + 3] == "_" or ja[i + 4] == "_" or ja[i + 5] == "_" or ja[i + 6] == "_" or ja[i + 7] == "_" or \
                    ja[i + 8] == "_" or ja[i + 9] == "_" or ja[i + 10] == "_" or ja[i + 11] == "_":
                print(ja[i - 57], ja[i - 56], ja[i - 55], ja[i - 54], ja[i - 53], ja[i - 52], ja[i - 51], ja[i - 50],
                      ja[i - 49],
                      ja[i - 48], ja[i - 47], ja[i - 46], ja[i - 45], ja[i - 44], ja[i - 43], ja[i - 42], ja[i - 41],
                      ja[i - 40],
                      ja[i - 39], ja[i - 38], ja[i - 37], ja[i - 36], ja[i - 35], ja[i - 34], ja[i - 33], ja[i - 32],
                      ja[i - 31],
                      ja[i - 30], ja[i - 29], ja[i - 28], ja[i - 27], ja[i - 26], ja[i - 25], ja[i - 24], ja[i - 23],
                      ja[i - 22],
                      ja[i - 21], ja[i - 20], ja[i - 19], ja[i - 18], ja[i - 17], ja[i - 16], ja[i - 15], ja[i - 14],
                      ja[i - 13],
                      ja[i - 12], ja[i - 11], ja[i - 10], ja[i - 9], ja[i - 8], ja[i - 7], ja[i - 6], ja[i - 5],
                      ja[i - 4],
                      ja[i - 3], ja[i - 2], ja[i - 1], ja[i], ja[i + 1], ja[i + 2], ja[i + 3], ja[i + 4],
                      ja[i + 5], ja[i + 6], ja[i + 7], ja[i + 8], ja[i + 9], ja[i + 10], ja[i + 11],
                      file=open("output4.txt", "a+"), sep="")
                strix = open("output4.txt")
                striy = strix.readlines()
                string = striy[-1]
                res = [sub.replace('_', '') for sub in string]
                x = "".join(res)
                print(x[-17:], file=open("oligo.txt", "a+"))
            else:
                print(ja[i - 5], ja[i - 4], ja[i - 3], ja[i - 2],
                      ja[i - 1], ja[i], ja[i + 1], ja[i + 2], ja[i + 3], ja[i + 4], ja[i + 5], ja[i + 6],
                      ja[i + 7], ja[i + 8], ja[i + 9], ja[i + 10], ja[i + 11], file=open("oligo.txt", "a+"), sep="")
        if i == (ln - 29):
            if ja[i - 4] == "_" or ja[i - 3] == "_" or \
                    ja[i - 2] == "_" or ja[i - 1] == "_" or ja[i] == "_" or ja[i + 1] == "_" or ja[i + 2] == "_" or \
                    ja[i + 3] == "_" or ja[i + 4] == "_" or ja[i + 5] == "_" or ja[i + 6] == "_" or ja[i + 7] == "_" or \
                    ja[i + 8] == "_" or ja[i + 9] == "_" or ja[i + 10] == "_" or ja[i + 11] == "_" or ja[i + 12] == "_":
                print(ja[i - 56], ja[i - 55], ja[i - 54], ja[i - 53], ja[i - 52], ja[i - 51], ja[i - 50],
                      ja[i - 49],
                      ja[i - 48], ja[i - 47], ja[i - 46], ja[i - 45], ja[i - 44], ja[i - 43], ja[i - 42], ja[i - 41],
                      ja[i - 40],
                      ja[i - 39], ja[i - 38], ja[i - 37], ja[i - 36], ja[i - 35], ja[i - 34], ja[i - 33], ja[i - 32],
                      ja[i - 31],
                      ja[i - 30], ja[i - 29], ja[i - 28], ja[i - 27], ja[i - 26], ja[i - 25], ja[i - 24], ja[i - 23],
                      ja[i - 22],
                      ja[i - 21], ja[i - 20], ja[i - 19], ja[i - 18], ja[i - 17], ja[i - 16], ja[i - 15], ja[i - 14],
                      ja[i - 13],
                      ja[i - 12], ja[i - 11], ja[i - 10], ja[i - 9], ja[i - 8], ja[i - 7], ja[i - 6], ja[i - 5],
                      ja[i - 4],
                      ja[i - 3], ja[i - 2], ja[i - 1], ja[i], ja[i + 1], ja[i + 2], ja[i + 3], ja[i + 4],
                      ja[i + 5], ja[i + 6], ja[i + 7], ja[i + 8], ja[i + 9], ja[i + 10], ja[i + 11], ja[i + 12],
                      file=open("output4.txt", "a+"), sep="")
                strix = open("output4.txt")
                striy = strix.readlines()
                string = striy[-1]
                res = [sub.replace('_', '') for sub in string]
                x = "".join(res)
                print(x[-17:], file=open("oligo.txt", "a+"))
            else:
                print(ja[i - 4], ja[i - 3], ja[i - 2],
                      ja[i - 1], ja[i], ja[i + 1], ja[i + 2], ja[i + 3], ja[i + 4], ja[i + 5], ja[i + 6],
                      ja[i + 7], ja[i + 8], ja[i + 9], ja[i + 10], ja[i + 11], ja[i + 12],
                      file=open("oligo.txt", "a+"), sep="")
        if i == (ln - 30):
            if ja[i - 3] == "_" or \
                    ja[i - 2] == "_" or ja[i - 1] == "_" or ja[i] == "_" or ja[i + 1] == "_" or ja[i + 2] == "_" or \
                    ja[i + 3] == "_" or ja[i + 4] == "_" or ja[i + 5] == "_" or ja[i + 6] == "_" or ja[i + 7] == "_" or \
                    ja[i + 8] == "_" or ja[i + 9] == "_" or ja[i + 10] == "_" or ja[i + 11] == "_" or ja[i + 12] == "_" \
                    or ja[i + 13] == "_":
                print(ja[i - 55], ja[i - 54], ja[i - 53], ja[i - 52], ja[i - 51], ja[i - 50],
                      ja[i - 49],
                      ja[i - 48], ja[i - 47], ja[i - 46], ja[i - 45], ja[i - 44], ja[i - 43], ja[i - 42], ja[i - 41],
                      ja[i - 40],
                      ja[i - 39], ja[i - 38], ja[i - 37], ja[i - 36], ja[i - 35], ja[i - 34], ja[i - 33], ja[i - 32],
                      ja[i - 31],
                      ja[i - 30], ja[i - 29], ja[i - 28], ja[i - 27], ja[i - 26], ja[i - 25], ja[i - 24], ja[i - 23],
                      ja[i - 22],
                      ja[i - 21], ja[i - 20], ja[i - 19], ja[i - 18], ja[i - 17], ja[i - 16], ja[i - 15], ja[i - 14],
                      ja[i - 13],
                      ja[i - 12], ja[i - 11], ja[i - 10], ja[i - 9], ja[i - 8], ja[i - 7], ja[i - 6], ja[i - 5],
                      ja[i - 4],
                      ja[i - 3], ja[i - 2], ja[i - 1], ja[i], ja[i + 1], ja[i + 2], ja[i + 3], ja[i + 4],
                      ja[i + 5], ja[i + 6], ja[i + 7], ja[i + 8], ja[i + 9], ja[i + 10], ja[i + 11], ja[i + 12],
                      ja[i + 13], file=open("output4.txt", "a+"), sep="")
                strix = open("output4.txt")
                striy = strix.readlines()
                string = striy[-1]
                res = [sub.replace('_', '') for sub in string]
                x = "".join(res)
                print(x[-17:], file=open("oligo.txt", "a+"))
            else:
                print(ja[i - 3], ja[i - 2],
                      ja[i - 1], ja[i], ja[i + 1], ja[i + 2], ja[i + 3], ja[i + 4], ja[i + 5], ja[i + 6],
                      ja[i + 7], ja[i + 8], ja[i + 9], ja[i + 10], ja[i + 11], ja[i + 12], ja[i + 13],
                      file=open("oligo.txt", "a+"), sep="")
        if i == (ln - 31):
            if ja[i - 2] == "_" or ja[i - 1] == "_" or ja[i] == "_" or ja[i + 1] == "_" or ja[i + 2] == "_" or \
                    ja[i + 3] == "_" or ja[i + 4] == "_" or ja[i + 5] == "_" or ja[i + 6] == "_" or ja[i + 7] == "_" or \
                    ja[i + 8] == "_" or ja[i + 9] == "_" or ja[i + 10] == "_" or ja[i + 11] == "_" or ja[i + 12] == "_" \
                    or ja[i + 13] == "_" or ja[i + 14] == "_":
                print(ja[i - 54], ja[i - 53], ja[i - 52], ja[i - 51], ja[i - 50],
                      ja[i - 49],
                      ja[i - 48], ja[i - 47], ja[i - 46], ja[i - 45], ja[i - 44], ja[i - 43], ja[i - 42], ja[i - 41],
                      ja[i - 40],
                      ja[i - 39], ja[i - 38], ja[i - 37], ja[i - 36], ja[i - 35], ja[i - 34], ja[i - 33], ja[i - 32],
                      ja[i - 31],
                      ja[i - 30], ja[i - 29], ja[i - 28], ja[i - 27], ja[i - 26], ja[i - 25], ja[i - 24], ja[i - 23],
                      ja[i - 22],
                      ja[i - 21], ja[i - 20], ja[i - 19], ja[i - 18], ja[i - 17], ja[i - 16], ja[i - 15], ja[i - 14],
                      ja[i - 13],
                      ja[i - 12], ja[i - 11], ja[i - 10], ja[i - 9], ja[i - 8], ja[i - 7], ja[i - 6], ja[i - 5],
                      ja[i - 4],
                      ja[i - 3], ja[i - 2], ja[i - 1], ja[i], ja[i + 1], ja[i + 2], ja[i + 3], ja[i + 4],
                      ja[i + 5], ja[i + 6], ja[i + 7], ja[i + 8], ja[i + 9], ja[i + 10], ja[i + 11], ja[i + 12],
                      ja[i + 13], ja[i + 14], file=open("output4.txt", "a+"), sep="")
                strix = open("output4.txt")
                striy = strix.readlines()
                string = striy[-1]
                res = [sub.replace('_', '') for sub in string]
                x = "".join(res)
                print(x[-17:], file=open("oligo.txt", "a+"))
            else:
                print(ja[i - 2],
                      ja[i - 1], ja[i], ja[i + 1], ja[i + 2], ja[i + 3], ja[i + 4], ja[i + 5], ja[i + 6],
                      ja[i + 7], ja[i + 8], ja[i + 9], ja[i + 10], ja[i + 11], ja[i + 12], ja[i + 13], ja[i + 14],
                      file=open("oligo.txt", "a+"), sep="")
        if i == (ln - 32):
            if ja[i - 1] == "_" or ja[i] == "_" or ja[i + 1] == "_" or ja[i + 2] == "_" or \
                    ja[i + 3] == "_" or ja[i + 4] == "_" or ja[i + 5] == "_" or ja[i + 6] == "_" or ja[i + 7] == "_" or \
                    ja[i + 8] == "_" or ja[i + 9] == "_" or ja[i + 10] == "_" or ja[i + 11] == "_" or ja[i + 12] == "_" \
                    or ja[i + 13] == "_" or ja[i + 14] == "_" or ja[i + 15] == "_":
                print(ja[i - 53], ja[i - 52], ja[i - 51], ja[i - 50],
                      ja[i - 49],
                      ja[i - 48], ja[i - 47], ja[i - 46], ja[i - 45], ja[i - 44], ja[i - 43], ja[i - 42], ja[i - 41],
                      ja[i - 40],
                      ja[i - 39], ja[i - 38], ja[i - 37], ja[i - 36], ja[i - 35], ja[i - 34], ja[i - 33], ja[i - 32],
                      ja[i - 31],
                      ja[i - 30], ja[i - 29], ja[i - 28], ja[i - 27], ja[i - 26], ja[i - 25], ja[i - 24], ja[i - 23],
                      ja[i - 22],
                      ja[i - 21], ja[i - 20], ja[i - 19], ja[i - 18], ja[i - 17], ja[i - 16], ja[i - 15], ja[i - 14],
                      ja[i - 13],
                      ja[i - 12], ja[i - 11], ja[i - 10], ja[i - 9], ja[i - 8], ja[i - 7], ja[i - 6], ja[i - 5],
                      ja[i - 4],
                      ja[i - 3], ja[i - 2], ja[i - 1], ja[i], ja[i + 1], ja[i + 2], ja[i + 3], ja[i + 4],
                      ja[i + 5], ja[i + 6], ja[i + 7], ja[i + 8], ja[i + 9], ja[i + 10], ja[i + 11], ja[i + 12],
                      ja[i + 13], ja[i + 14], ja[i + 15], file=open("output4.txt", "a+"), sep="")
                strix = open("output4.txt")
                striy = strix.readlines()
                string = striy[-1]
                res = [sub.replace('_', '') for sub in string]
                x = "".join(res)
                print(x[-17:], file=open("oligo.txt", "a+"))
            else:
                print(ja[i - 1], ja[i], ja[i + 1], ja[i + 2], ja[i + 3], ja[i + 4], ja[i + 5], ja[i + 6],
                      ja[i + 7], ja[i + 8], ja[i + 9], ja[i + 10], ja[i + 11], ja[i + 12], ja[i + 13], ja[i + 14],
                      ja[i + 15], file=open("oligo.txt", "a+"), sep="")
        if i == (ln - 33):
            if ja[i] == "_" or ja[i + 1] == "_" or ja[i + 2] == "_" or \
                    ja[i + 3] == "_" or ja[i + 4] == "_" or ja[i + 5] == "_" or ja[i + 6] == "_" or ja[i + 7] == "_" or \
                    ja[i + 8] == "_" or ja[i + 9] == "_" or ja[i + 10] == "_" or ja[i + 11] == "_" or ja[i + 12] == "_" \
                    or ja[i + 13] == "_" or ja[i + 14] == "_" or ja[i + 15] == "_" or ja[i + 16] == "_":
                print(ja[i - 52], ja[i - 51], ja[i - 50],
                      ja[i - 49],
                      ja[i - 48], ja[i - 47], ja[i - 46], ja[i - 45], ja[i - 44], ja[i - 43], ja[i - 42], ja[i - 41],
                      ja[i - 40],
                      ja[i - 39], ja[i - 38], ja[i - 37], ja[i - 36], ja[i - 35], ja[i - 34], ja[i - 33], ja[i - 32],
                      ja[i - 31],
                      ja[i - 30], ja[i - 29], ja[i - 28], ja[i - 27], ja[i - 26], ja[i - 25], ja[i - 24], ja[i - 23],
                      ja[i - 22],
                      ja[i - 21], ja[i - 20], ja[i - 19], ja[i - 18], ja[i - 17], ja[i - 16], ja[i - 15], ja[i - 14],
                      ja[i - 13],
                      ja[i - 12], ja[i - 11], ja[i - 10], ja[i - 9], ja[i - 8], ja[i - 7], ja[i - 6], ja[i - 5],
                      ja[i - 4],
                      ja[i - 3], ja[i - 2], ja[i - 1], ja[i], ja[i + 1], ja[i + 2], ja[i + 3], ja[i + 4],
                      ja[i + 5], ja[i + 6], ja[i + 7], ja[i + 8], ja[i + 9], ja[i + 10], ja[i + 11], ja[i + 12],
                      ja[i + 13], ja[i + 14], ja[i + 15], ja[i + 16], file=open("output4.txt", "a+"), sep="")
                strix = open("output4.txt")
                striy = strix.readlines()
                string = striy[-1]
                res = [sub.replace('_', '') for sub in string]
                x = "".join(res)
                print(x[-17:], file=open("oligo.txt", "a+"))
            else:
                print(ja[i], ja[i + 1], ja[i + 2], ja[i + 3], ja[i + 4], ja[i + 5], ja[i + 6],
                      ja[i + 7], ja[i + 8], ja[i + 9], ja[i + 10], ja[i + 11], ja[i + 12], ja[i + 13], ja[i + 14],
                      ja[i + 15], ja[i + 16], file=open("oligo.txt", "a+"), sep="")
        if i == (ln - 34):
            if ja[i] == "_" or ja[i - 1] == "_" or ja[i - 2] == "_" or \
                    ja[i - 3] == "_" or ja[i - 4] == "_" or ja[i - 5] == "_" or ja[i - 6] == "_" or ja[i - 7] == "_" or \
                    ja[i - 8] == "_" or ja[i - 9] == "_" or ja[i - 10] == "_" or ja[i - 11] == "_" or ja[i - 12] == "_" \
                    or ja[i - 13] == "_" or ja[i - 14] == "_" or ja[i - 15] == "_" or ja[i - 16] == "_":
                print(ja[i - 68], ja[i - 67], ja[i - 66],
                      ja[i - 65],
                      ja[i - 64], ja[i - 63], ja[i - 62], ja[i - 61], ja[i - 60], ja[i - 59], ja[i - 58], ja[i - 57],
                      ja[i - 56],
                      ja[i - 55], ja[i - 54], ja[i - 53], ja[i - 52], ja[i - 51], ja[i - 50], ja[i - 49], ja[i - 48],
                      ja[i - 47],
                      ja[i - 46], ja[i - 45], ja[i - 44], ja[i - 43], ja[i - 42], ja[i - 41], ja[i - 40], ja[i - 39],
                      ja[i - 38],
                      ja[i - 37], ja[i - 36], ja[i - 35], ja[i - 34], ja[i - 33], ja[i - 32], ja[i - 31], ja[i - 30],
                      ja[i - 29],
                      ja[i - 28], ja[i - 27], ja[i - 26], ja[i - 25], ja[i - 24], ja[i - 23], ja[i - 22], ja[i - 21],
                      ja[i - 20],
                      ja[i - 19], ja[i - 18], ja[i - 17], ja[i - 16], ja[i - 15], ja[i - 14], ja[i - 13], ja[i - 12],
                      ja[i - 11], ja[i - 10], ja[i - 9], ja[i - 8], ja[i - 7], ja[i - 6], ja[i - 5], ja[i - 4],
                      ja[i - 3], ja[i - 2], ja[i - 1], ja[i], file=open("output4.txt", "a+"), sep="")
                strix = open("output4.txt")
                striy = strix.readlines()
                string = striy[-1]
                res = [sub.replace('_', '') for sub in string]
                x = "".join(res)
                print(x[-17:], file=open("oligo.txt", "a+"))
            else:
                print(ja[i - 16], ja[i - 15], ja[i - 14], ja[i - 13], ja[i - 12], ja[i - 11], ja[i - 10],
                      ja[i - 9], ja[i - 8], ja[i - 7], ja[i - 6], ja[i - 5], ja[i - 4], ja[i - 3], ja[i - 2],
                      ja[i - 1], ja[i], file=open("oligo.txt", "a+"), sep="")
        if i >= 35 or i <= (ln - 35):
            if ja[i - 8] == "_" or ja[i - 7] == "_" or ja[i - 6] == "_" or ja[i - 5] == "_" or ja[i - 4] == "_" or \
                    ja[i - 3] == "_" or ja[i - 2] == "_" or ja[i - 1] == "_" or ja[i] == "_" or ja[i + 1] == "_" or \
                    ja[i + 2] == "_" or ja[i + 3] == "_" or ja[i + 4] == "_" or ja[i + 5] == "_" or ja[i + 6] == "_" or \
                    ja[i + 7] == "_" or ja[i + 8] == "_":
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
            elif ja[i - 8] != "_" or ja[i - 7] != "_" or ja[i - 6] != "_" or ja[i - 5] != "_" or ja[i - 4] != "_" or \
                    ja[i - 3] != "_" or ja[i - 2] != "_" or ja[i - 1] != "_" or ja[i] != "_" or ja[i + 1] != "_" or \
                    ja[i + 2] != "_" or ja[i + 3] != "_" or ja[i + 4] != "_" or ja[i + 5] != "_" or ja[i + 6] != "_" or \
                    ja[i + 7] != "_" or ja[i + 8] != "_":
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
        if i == 0:
            if ja[i] == "_" or ja[i + 1] == "_" or ja[i + 2] == "_" or ja[i + 3] == "_" or ja[i + 4] == "_" or \
                    ja[i + 5] == "_" or ja[i + 6] == "_" or ja[i + 7] == "_" or ja[i + 8] == "_" or ja[i + 9] == "_" or \
                    ja[i + 10] == "_" or ja[i + 11] == "_" or ja[i + 12] == "_" or ja[i + 13] == "_" or \
                    ja[i + 14] == "_" or ja[i + 15] == "_" or ja[i + 16] == "_":
                print(ja[i], ja[i + 1], ja[i + 2], ja[i + 3], ja[i + 4], ja[i + 5], ja[i + 6],
                      ja[i + 7],
                      ja[i + 8], ja[i + 9], ja[i + 10], ja[i + 11], ja[i + 12], ja[i + 13], ja[i + 14],
                      ja[i + 15],
                      ja[i + 16], ja[i + 17], ja[i + 18], ja[i + 19], ja[i + 20], ja[i + 21], ja[i + 22],
                      ja[i + 23],
                      ja[i + 24], ja[i + 25], ja[i + 26], ja[i + 27], ja[i + 28], ja[i + 29], ja[i + 30], ja[i + 31],
                      ja[i + 32],
                      ja[i + 33], ja[i + 34], ja[i + 35], ja[i + 36], ja[i + 37], ja[i + 38], ja[i + 39], ja[i + 40],
                      ja[i + 41],
                      ja[i + 42], ja[i + 43], ja[i + 44], ja[i + 45], ja[i + 46], ja[i + 47], ja[i + 48], ja[i + 49],
                      ja[i + 50], ja[i + 51], ja[i + 52], ja[i + 53], ja[i + 54], ja[i + 55], ja[i + 56], ja[i + 57],
                      ja[i + 58], ja[i + 59], ja[i + 60], ja[i + 61], ja[i + 62], ja[i + 63], ja[i + 64], ja[i + 65],
                      ja[i + 66], ja[i + 67], ja[i + 68], file=open("output4b.txt", "a+"), sep="")
                strix = open("output4b.txt")
                striy = strix.readlines()
                string = striy[-1]
                res = [sub.replace('_', '') for sub in string]
                x = "".join(res)
                print(x[:17], file=open("oligob.txt", "a+"))
            elif ja[i] != "_" or ja[i + 1] != "_" or ja[i + 2] != "_" or ja[i + 3] != "_" or ja[i + 4] != "_" or \
                    ja[i + 5] != "_" or ja[i + 6] != "_" or ja[i + 7] != "_" or ja[i + 9] != "_" or ja[i + 10] != "_" or \
                    ja[i + 11] != "_" or ja[i + 12] != "_" or ja[i + 13] != "_" or ja[i + 14] != "_" or \
                    ja[i + 15] != "_" or ja[i + 16] != "_":
                print(ja[i], ja[i + 1], ja[i + 2], ja[i + 3], ja[i + 4], ja[i + 5], ja[i + 6], ja[i + 7], ja[i + 8],
                      ja[i + 9], ja[i + 10], ja[i + 11], ja[i + 12], ja[i + 13], ja[i + 14], ja[i + 15], ja[i + 16],
                      file=open("oligob.txt", "a+"), sep="")
        if i == 1:
            if ja[i - 1] == "_" or ja[i] == "_" or ja[i + 1] == "_" or ja[i + 2] == "_" or ja[i + 3] == "_" or \
                    ja[i + 4] == "_" or ja[i + 5] == "_" or ja[i + 6] == "_" or ja[i + 7] == "_" or ja[i + 8] == "_" or \
                    ja[i + 9] == "_" or ja[i + 10] == "_" or ja[i + 11] == "_" or ja[i + 12] == "_" or \
                    ja[i + 13] == "_" or ja[i + 14] == "_" or ja[i + 15] == "_":
                print(ja[i - 1], ja[i], ja[i + 1], ja[i + 2], ja[i + 3], ja[i + 4], ja[i + 5], ja[i + 6],
                      ja[i + 7],
                      ja[i + 8], ja[i + 9], ja[i + 10], ja[i + 11], ja[i + 12], ja[i + 13], ja[i + 14],
                      ja[i + 15],
                      ja[i + 16], ja[i + 17], ja[i + 18], ja[i + 19], ja[i + 20], ja[i + 21], ja[i + 22],
                      ja[i + 23],
                      ja[i + 24], ja[i + 25], ja[i + 26], ja[i + 27], ja[i + 28], ja[i + 29], ja[i + 30],
                      ja[i + 31],
                      ja[i + 32],
                      ja[i + 33], ja[i + 34], ja[i + 35], ja[i + 36], ja[i + 37], ja[i + 38], ja[i + 39],
                      ja[i + 40], ja[i + 41],
                      ja[i + 42], ja[i + 43], ja[i + 44], ja[i + 45], ja[i + 46], ja[i + 47], ja[i + 48],
                      ja[i + 49],
                      ja[i + 50], ja[i + 51], ja[i + 52], ja[i + 53], ja[i + 54], ja[i + 55], ja[i + 56],
                      ja[i + 57],
                      ja[i + 58], ja[i + 59], ja[i + 60], ja[i + 61], ja[i + 62], ja[i + 63], ja[i + 64],
                      ja[i + 65],
                      ja[i + 66], ja[i + 67], file=open("output4b.txt", "a+"), sep="")
                strix = open("output4b.txt")
                striy = strix.readlines()
                string = striy[-1]
                res = [sub.replace('_', '') for sub in string]
                x = "".join(res)
                print(x[:17], file=open("oligob.txt", "a+"))
            else:
                print(ja[i - 1], ja[i], ja[i + 1], ja[i + 2], ja[i + 3], ja[i + 4], ja[i + 5], ja[i + 6], ja[i + 7],
                      ja[i + 8], ja[i + 9], ja[i + 10], ja[i + 11], ja[i + 12], ja[i + 13], ja[i + 14], ja[i + 15],
                      file=open("oligob.txt", "a+"), sep="")
        if i == 2:
            if ja[i - 2] == "_" or ja[i - 1] == "_" or ja[i] == "_" or ja[i + 1] == "_" or ja[i + 2] == "_" or \
                    ja[i + 3] == "_" or ja[i + 4] == "_" or ja[i + 5] == "_" or ja[i + 6] == "_" or ja[i + 7] == "_" or \
                    ja[i + 8] == "_" or ja[i + 9] == "_" or ja[i + 10] == "_" or ja[i + 11] == "_" or \
                    ja[i + 12] == "_" or ja[i + 13] == "_" or ja[i + 14] == "_":
                print(ja[i - 2], ja[i - 1], ja[i], ja[i + 1], ja[i + 2], ja[i + 3], ja[i + 4], ja[i + 5], ja[i + 6],
                      ja[i + 7],
                      ja[i + 8], ja[i + 9], ja[i + 10], ja[i + 11], ja[i + 12], ja[i + 13], ja[i + 14],
                      ja[i + 15],
                      ja[i + 16], ja[i + 17], ja[i + 18], ja[i + 19], ja[i + 20], ja[i + 21], ja[i + 22],
                      ja[i + 23],
                      ja[i + 24], ja[i + 25], ja[i + 26], ja[i + 27], ja[i + 28], ja[i + 29], ja[i + 30],
                      ja[i + 31],
                      ja[i + 32],
                      ja[i + 33], ja[i + 34], ja[i + 35], ja[i + 36], ja[i + 37], ja[i + 38], ja[i + 39],
                      ja[i + 40], ja[i + 41],
                      ja[i + 42], ja[i + 43], ja[i + 44], ja[i + 45], ja[i + 46], ja[i + 47], ja[i + 48],
                      ja[i + 49],
                      ja[i + 50], ja[i + 51], ja[i + 52], ja[i + 53], ja[i + 54], ja[i + 55], ja[i + 56],
                      ja[i + 57],
                      ja[i + 58], ja[i + 59], ja[i + 60], ja[i + 61], ja[i + 62], ja[i + 63], ja[i + 64],
                      ja[i + 65],
                      ja[i + 66], file=open("output4b.txt", "a+"), sep="")
                strix = open("output4b.txt")
                striy = strix.readlines()
                string = striy[-1]
                res = [sub.replace('_', '') for sub in string]
                x = "".join(res)
                print(x[:17], file=open("oligob.txt", "a+"))
            else:
                print(ja[i - 2], ja[i - 1], ja[i], ja[i + 1], ja[i + 2], ja[i + 3], ja[i + 4], ja[i + 5], ja[i + 6],
                      ja[i + 7], ja[i + 8], ja[i + 9], ja[i + 10], ja[i + 11], ja[i + 12], ja[i + 13],
                      ja[i + 14], file=open("oligob.txt", "a+"), sep="")
        if i == 3:
            if ja[i - 3] == "_" or ja[i - 2] == "_" or ja[i - 1] == "_" or ja[i] == "_" or ja[i + 1] == "_" or \
                    ja[i + 2] == "_" or ja[i + 3] == "_" or ja[i + 4] == "_" or ja[i + 5] == "_" or ja[i + 6] == "_" or \
                    ja[i + 7] == "_" or ja[i + 8] == "_" or ja[i + 9] == "_" or ja[i + 10] == "_" or \
                    ja[i + 11] == "_" or ja[i + 12] == "_" or ja[i + 13] == "_":
                print(ja[i - 3], ja[i - 2], ja[i - 1], ja[i], ja[i + 1], ja[i + 2], ja[i + 3], ja[i + 4], ja[i + 5],
                      ja[i + 6],
                      ja[i + 7],
                      ja[i + 8], ja[i + 9], ja[i + 10], ja[i + 11], ja[i + 12], ja[i + 13], ja[i + 14],
                      ja[i + 15],
                      ja[i + 16], ja[i + 17], ja[i + 18], ja[i + 19], ja[i + 20], ja[i + 21], ja[i + 22],
                      ja[i + 23],
                      ja[i + 24], ja[i + 25], ja[i + 26], ja[i + 27], ja[i + 28], ja[i + 29], ja[i + 30],
                      ja[i + 31],
                      ja[i + 32],
                      ja[i + 33], ja[i + 34], ja[i + 35], ja[i + 36], ja[i + 37], ja[i + 38], ja[i + 39],
                      ja[i + 40], ja[i + 41],
                      ja[i + 42], ja[i + 43], ja[i + 44], ja[i + 45], ja[i + 46], ja[i + 47], ja[i + 48],
                      ja[i + 49],
                      ja[i + 50], ja[i + 51], ja[i + 52], ja[i + 53], ja[i + 54], ja[i + 55], ja[i + 56],
                      ja[i + 57],
                      ja[i + 58], ja[i + 59], ja[i + 60], ja[i + 61], ja[i + 62], ja[i + 63], ja[i + 64],
                      ja[i + 65], file=open("output4b.txt", "a+"), sep="")
                strix = open("output4b.txt")
                striy = strix.readlines()
                string = striy[-1]
                res = [sub.replace('_', '') for sub in string]
                x = "".join(res)
                print(x[:17], file=open("oligob.txt", "a+"))
            else:
                print(ja[i - 3], ja[i - 2], ja[i - 1], ja[i], ja[i + 1], ja[i + 2], ja[i + 3], ja[i + 4], ja[i + 5],
                      ja[i + 6],
                      ja[i + 7], ja[i + 8], ja[i + 9], ja[i + 10], ja[i + 11], ja[i + 12], ja[i + 13],
                      file=open("oligob.txt", "a+"), sep="")
        if i == 4:
            if ja[i - 4] == "_" or ja[i - 3] == "_" or ja[i - 2] == "_" or ja[i - 1] == "_" or ja[i] == "_" or \
                    ja[i + 1] == "_" or ja[i + 2] == "_" or ja[i + 3] == "_" or ja[i + 4] == "_" or \
                    ja[i + 5] == "_" or ja[i + 6] == "_" or ja[i + 7] == "_" or ja[i + 8] == "_" or \
                    ja[i + 9] == "_" or ja[i + 10] == "_" or ja[i + 11] == "_" or ja[i + 12] == "_":
                print(ja[i - 4], ja[i - 3], ja[i - 2], ja[i - 1], ja[i], ja[i + 1], ja[i + 2], ja[i + 3], ja[i + 4],
                      ja[i + 5], ja[i + 6],
                      ja[i + 7],
                      ja[i + 8], ja[i + 9], ja[i + 10], ja[i + 11], ja[i + 12], ja[i + 13], ja[i + 14],
                      ja[i + 15],
                      ja[i + 16], ja[i + 17], ja[i + 18], ja[i + 19], ja[i + 20], ja[i + 21], ja[i + 22],
                      ja[i + 23],
                      ja[i + 24], ja[i + 25], ja[i + 26], ja[i + 27], ja[i + 28], ja[i + 29], ja[i + 30],
                      ja[i + 31],
                      ja[i + 32],
                      ja[i + 33], ja[i + 34], ja[i + 35], ja[i + 36], ja[i + 37], ja[i + 38], ja[i + 39],
                      ja[i + 40], ja[i + 41],
                      ja[i + 42], ja[i + 43], ja[i + 44], ja[i + 45], ja[i + 46], ja[i + 47], ja[i + 48],
                      ja[i + 49],
                      ja[i + 50], ja[i + 51], ja[i + 52], ja[i + 53], ja[i + 54], ja[i + 55], ja[i + 56],
                      ja[i + 57],
                      ja[i + 58], ja[i + 59], ja[i + 60], ja[i + 61], ja[i + 62], ja[i + 63], ja[i + 64],
                      file=open("output4b.txt", "a+"),
                      sep="")
                strix = open("output4b.txt")
                striy = strix.readlines()
                string = striy[-1]
                res = [sub.replace('_', '') for sub in string]
                x = "".join(res)
                print(x[:17], file=open("oligob.txt", "a+"))
            else:
                print(ja[i - 4], ja[i - 3], ja[i - 2], ja[i - 1], ja[i], ja[i + 1], ja[i + 2], ja[i + 3], ja[i + 4],
                      ja[i + 5], ja[i + 6],
                      ja[i + 7], ja[i + 8], ja[i + 9], ja[i + 10], ja[i + 11], ja[i + 12],
                      file=open("oligob.txt", "a+"),
                      sep="")
        if i == 5:
            if ja[i - 5] == "_" or ja[i - 4] == "_" or ja[i - 3] == "_" or ja[i - 2] == "_" or ja[i - 1] == "_" or ja[
                i] == "_" or ja[i + 1] == "_" or ja[i + 2] == "_" or ja[i + 3] == "_" or ja[i + 4] == "_" or \
                    ja[i + 5] == "_" or ja[i + 6] == "_" or ja[i + 7] == "_" or ja[i + 8] == "_" or ja[i + 9] == "_" \
                    or ja[i + 10] == "_" or ja[i + 11] == "_":
                print(ja[i - 5], ja[i - 4], ja[i - 3], ja[i - 2], ja[i - 1], ja[i], ja[i + 1], ja[i + 2], ja[i + 3],
                      ja[i + 4], ja[i + 5], ja[i + 6],
                      ja[i + 7], ja[i + 8], ja[i + 9], ja[i + 10], ja[i + 11], ja[i + 12], ja[i + 13], ja[i + 14],
                      ja[i + 15], ja[i + 16], ja[i + 17], ja[i + 18], ja[i + 19], ja[i + 20], ja[i + 21], ja[i + 22],
                      ja[i + 23], ja[i + 24], ja[i + 25], ja[i + 26], ja[i + 27], ja[i + 28], ja[i + 29], ja[i + 30],
                      ja[i + 31], ja[i + 32], ja[i + 33], ja[i + 34], ja[i + 35], ja[i + 36], ja[i + 37],
                      ja[i + 38], ja[i + 39],
                      ja[i + 40], ja[i + 41], ja[i + 42], ja[i + 43], ja[i + 44], ja[i + 45], ja[i + 46], ja[i + 47],
                      ja[i + 48],
                      ja[i + 49], ja[i + 50], ja[i + 51], ja[i + 52], ja[i + 53], ja[i + 54], ja[i + 55], ja[i + 56],
                      ja[i + 57], ja[i + 58], ja[i + 59], ja[i + 60], ja[i + 61], ja[i + 62], ja[i + 63],
                      file=open("output4b.txt", "a+"), sep="")
                strix = open("output4b.txt")
                striy = strix.readlines()
                string = striy[-1]
                res = [sub.replace('_', '') for sub in string]
                x = "".join(res)
                print(x[:17], file=open("oligob.txt", "a+"))
            else:
                print(ja[i - 5], ja[i - 4], ja[i - 3], ja[i - 2], ja[i - 1], ja[i], ja[i + 1], ja[i + 2], ja[i + 3],
                      ja[i + 4],
                      ja[i + 5], ja[i + 6],
                      ja[i + 7], ja[i + 8], ja[i + 9], ja[i + 10], ja[i + 11], file=open("oligob.txt", "a+"),
                      sep="")
        if i == 6:
            if ja[i - 6] == "_" or ja[i - 5] == "_" or ja[i - 4] == "_" or ja[i - 3] == "_" or ja[i - 2] == "_" or \
                    ja[i - 1] == "_" or ja[i] == "_" or ja[i + 1] == "_" or ja[i + 2] == "_" or ja[i + 3] == "_" or \
                    ja[i + 4] == "_" or ja[i + 5] == "_" or ja[i + 6] == "_" or ja[i + 7] == "_" or ja[i + 8] == "_" or \
                    ja[i + 9] == "_" or ja[i + 10] == "_":
                print(ja[i - 6], ja[i - 5], ja[i - 4], ja[i - 3], ja[i - 2], ja[i - 1], ja[i], ja[i + 1], ja[i + 2],
                      ja[i + 3],
                      ja[i + 4], ja[i + 5], ja[i + 6],
                      ja[i + 7], ja[i + 8], ja[i + 9], ja[i + 10], ja[i + 11], ja[i + 12], ja[i + 13], ja[i + 14],
                      ja[i + 15], ja[i + 16], ja[i + 17], ja[i + 18], ja[i + 19], ja[i + 20], ja[i + 21], ja[i + 22],
                      ja[i + 23], ja[i + 24], ja[i + 25], ja[i + 26], ja[i + 27], ja[i + 28], ja[i + 29], ja[i + 30],
                      ja[i + 31], ja[i + 32], ja[i + 33], ja[i + 34], ja[i + 35], ja[i + 36], ja[i + 37],
                      ja[i + 38], ja[i + 39],
                      ja[i + 40], ja[i + 41], ja[i + 42], ja[i + 43], ja[i + 44], ja[i + 45], ja[i + 46], ja[i + 47],
                      ja[i + 48],
                      ja[i + 49], ja[i + 50], ja[i + 51], ja[i + 52], ja[i + 53], ja[i + 54], ja[i + 55], ja[i + 56],
                      ja[i + 57], ja[i + 58], ja[i + 59], ja[i + 60], ja[i + 61], ja[i + 62],
                      file=open("output4b.txt", "a+"), sep="")
                strix = open("output4b.txt")
                striy = strix.readlines()
                string = striy[-1]
                res = [sub.replace('_', '') for sub in string]
                x = "".join(res)
                print(x[:17], file=open("oligob.txt", "a+"))
            else:
                print(ja[i - 6], ja[i - 5], ja[i - 4], ja[i - 3], ja[i - 2], ja[i - 1], ja[i], ja[i + 1], ja[i + 2],
                      ja[i + 3], ja[i + 4], ja[i + 5], ja[i + 6],
                      ja[i + 7], ja[i + 8], ja[i + 9], ja[i + 10], file=open("oligob.txt", "a+"), sep="")
        if i == 7:
            if ja[i - 7] == "_" or ja[i - 6] == "_" or ja[i - 5] == "_" or ja[i - 4] == "_" or ja[i - 3] == "_" or ja[
                i - 2] == "_" or \
                    ja[i - 1] == "_" or ja[i] == "_" or ja[i + 1] == "_" or ja[i + 2] == "_" or ja[i + 3] == "_" or \
                    ja[i + 4] == "_" or ja[i + 5] == "_" or ja[i + 6] == "_" or ja[i + 7] == "_" or ja[i + 8] == "_" or \
                    ja[i + 9] == "_":
                print(ja[i - 7], ja[i - 6], ja[i - 5], ja[i - 4], ja[i - 3], ja[i - 2], ja[i - 1], ja[i], ja[i + 1],
                      ja[i + 2], ja[i + 3],
                      ja[i + 4], ja[i + 5], ja[i + 6],
                      ja[i + 7], ja[i + 8], ja[i + 9], ja[i + 10], ja[i + 11], ja[i + 12], ja[i + 13], ja[i + 14],
                      ja[i + 15], ja[i + 16], ja[i + 17], ja[i + 18], ja[i + 19], ja[i + 20], ja[i + 21], ja[i + 22],
                      ja[i + 23], ja[i + 24], ja[i + 25], ja[i + 26], ja[i + 27], ja[i + 28], ja[i + 29], ja[i + 30],
                      ja[i + 31], ja[i + 32], ja[i + 33], ja[i + 34], ja[i + 35], ja[i + 36], ja[i + 37],
                      ja[i + 38], ja[i + 39],
                      ja[i + 40], ja[i + 41], ja[i + 42], ja[i + 43], ja[i + 44], ja[i + 45], ja[i + 46], ja[i + 47],
                      ja[i + 48],
                      ja[i + 49], ja[i + 50], ja[i + 51], ja[i + 52], ja[i + 53], ja[i + 54], ja[i + 55], ja[i + 56],
                      ja[i + 57], ja[i + 58], ja[i + 59], ja[i + 60], ja[i + 61],
                      file=open("output4b.txt", "a+"), sep="")
                strix = open("output4b.txt")
                striy = strix.readlines()
                string = striy[-1]
                res = [sub.replace('_', '') for sub in string]
                x = "".join(res)
                print(x[:17], file=open("oligob.txt", "a+"))
            else:
                print(ja[i - 7], ja[i - 6], ja[i - 5], ja[i - 4], ja[i - 3], ja[i - 2], ja[i - 1], ja[i], ja[i + 1],
                      ja[i + 2], ja[i + 3], ja[i + 4], ja[i + 5], ja[i + 6],
                      ja[i + 7], ja[i + 8], ja[i + 9], file=open("oligob.txt", "a+"), sep="")
        if i == 8:
            if ja[i - 8] == "_" or ja[i - 7] == "_" or ja[i - 6] == "_" or ja[i - 5] == "_" or ja[i - 4] == "_" or ja[
                i - 3] == "_" or ja[i - 2] == "_" or \
                    ja[i - 1] == "_" or ja[i] == "_" or ja[i + 1] == "_" or ja[i + 2] == "_" or ja[i + 3] == "_" or \
                    ja[i + 4] == "_" or ja[i + 5] == "_" or ja[i + 6] == "_" or ja[i + 7] == "_" or ja[i + 8] == "_":
                print(ja[i - 8], ja[i - 7], ja[i - 6], ja[i - 5], ja[i - 4], ja[i - 3], ja[i - 2], ja[i - 1], ja[i],
                      ja[i + 1], ja[i + 2], ja[i + 3],
                      ja[i + 4], ja[i + 5], ja[i + 6],
                      ja[i + 7], ja[i + 8], ja[i + 9], ja[i + 10], ja[i + 11], ja[i + 12], ja[i + 13], ja[i + 14],
                      ja[i + 15], ja[i + 16], ja[i + 17], ja[i + 18], ja[i + 19], ja[i + 20], ja[i + 21], ja[i + 22],
                      ja[i + 23], ja[i + 24], ja[i + 25], ja[i + 26], ja[i + 27], ja[i + 28], ja[i + 29], ja[i + 30],
                      ja[i + 31], ja[i + 32], ja[i + 33], ja[i + 34], ja[i + 35], ja[i + 36], ja[i + 37],
                      ja[i + 38], ja[i + 39],
                      ja[i + 40], ja[i + 41], ja[i + 42], ja[i + 43], ja[i + 44], ja[i + 45], ja[i + 46], ja[i + 47],
                      ja[i + 48],
                      ja[i + 49], ja[i + 50], ja[i + 51], ja[i + 52], ja[i + 53], ja[i + 54], ja[i + 55], ja[i + 56],
                      ja[i + 57], ja[i + 58], ja[i + 59], ja[i + 60],
                      file=open("output4b.txt", "a+"), sep="")
                strix = open("output4b.txt")
                striy = strix.readlines()
                string = striy[-1]
                res = [sub.replace('_', '') for sub in string]
                x = "".join(res)
                print(x[:17], file=open("oligob.txt", "a+"))
            else:
                print(ja[i - 8], ja[i - 7], ja[i - 6], ja[i - 5], ja[i - 4], ja[i - 3], ja[i - 2], ja[i - 1], ja[i],
                      ja[i + 1],
                      ja[i + 2], ja[i + 3], ja[i + 4], ja[i + 5], ja[i + 6],
                      ja[i + 7], ja[i + 8], file=open("oligob.txt", "a+"), sep="")
        if i == 9:
            if ja[i - 9] == "_" or ja[i - 8] == "_" or ja[i - 7] == "_" or ja[i - 6] == "_" or ja[i - 5] == "_" or ja[
                i - 4] == "_" or ja[i - 3] == "_" or ja[i - 2] == "_" or \
                    ja[i - 1] == "_" or ja[i] == "_" or ja[i + 1] == "_" or ja[i + 2] == "_" or ja[i + 3] == "_" or \
                    ja[i + 4] == "_" or ja[i + 5] == "_" or ja[i + 6] == "_" or ja[i + 7] == "_":
                print(ja[i - 9], ja[i - 8], ja[i - 7], ja[i - 6], ja[i - 5], ja[i - 4], ja[i - 3], ja[i - 2], ja[i - 1],
                      ja[i], ja[i + 1], ja[i + 2], ja[i + 3],
                      ja[i + 4], ja[i + 5], ja[i + 6],
                      ja[i + 7], ja[i + 8], ja[i + 9], ja[i + 10], ja[i + 11], ja[i + 12], ja[i + 13], ja[i + 14],
                      ja[i + 15], ja[i + 16], ja[i + 17], ja[i + 18], ja[i + 19], ja[i + 20], ja[i + 21], ja[i + 22],
                      ja[i + 23], ja[i + 24], ja[i + 25], ja[i + 26], ja[i + 27], ja[i + 28], ja[i + 29], ja[i + 30],
                      ja[i + 31], ja[i + 32], ja[i + 33], ja[i + 34], ja[i + 35], ja[i + 36], ja[i + 37],
                      ja[i + 38], ja[i + 39],
                      ja[i + 40], ja[i + 41], ja[i + 42], ja[i + 43], ja[i + 44], ja[i + 45], ja[i + 46], ja[i + 47],
                      ja[i + 48],
                      ja[i + 49], ja[i + 50], ja[i + 51], ja[i + 52], ja[i + 53], ja[i + 54], ja[i + 55], ja[i + 56],
                      ja[i + 57], ja[i + 58], ja[i + 59],
                      file=open("output4b.txt", "a+"), sep="")
                strix = open("output4b.txt")
                striy = strix.readlines()
                string = striy[-1]
                res = [sub.replace('_', '') for sub in string]
                x = "".join(res)
                print(x[:17], file=open("oligob.txt", "a+"))
            else:
                print(ja[i - 9], ja[i - 8], ja[i - 7], ja[i - 6], ja[i - 5], ja[i - 4], ja[i - 3], ja[i - 2], ja[i - 1],
                      ja[i], ja[i + 1],
                      ja[i + 2], ja[i + 3], ja[i + 4], ja[i + 5], ja[i + 6],
                      ja[i + 7], file=open("oligob.txt", "a+"), sep="")
        if i == 10:
            if ja[i - 10] == "_" or ja[i - 9] == "_" or ja[i - 8] == "_" or ja[i - 7] == "_" or ja[i - 6] == "_" or ja[
                i - 5] == "_" or ja[i - 4] == "_" or ja[i - 3] == "_" or ja[i - 2] == "_" or \
                    ja[i - 1] == "_" or ja[i] == "_" or ja[i + 1] == "_" or ja[i + 2] == "_" or ja[i + 3] == "_" or \
                    ja[i + 4] == "_" or ja[i + 5] == "_" or ja[i + 6] == "_":
                print(ja[i - 10], ja[i - 9], ja[i - 8], ja[i - 7], ja[i - 6], ja[i - 5], ja[i - 4], ja[i - 3],
                      ja[i - 2],
                      ja[i - 1], ja[i], ja[i + 1], ja[i + 2], ja[i + 3],
                      ja[i + 4], ja[i + 5], ja[i + 6],
                      ja[i + 7], ja[i + 8], ja[i + 9], ja[i + 10], ja[i + 11], ja[i + 12], ja[i + 13], ja[i + 14],
                      ja[i + 15], ja[i + 16], ja[i + 17], ja[i + 18], ja[i + 19], ja[i + 20], ja[i + 21], ja[i + 22],
                      ja[i + 23], ja[i + 24], ja[i + 25], ja[i + 26], ja[i + 27], ja[i + 28], ja[i + 29], ja[i + 30],
                      ja[i + 31], ja[i + 32], ja[i + 33], ja[i + 34], ja[i + 35], ja[i + 36], ja[i + 37],
                      ja[i + 38], ja[i + 39],
                      ja[i + 40], ja[i + 41], ja[i + 42], ja[i + 43], ja[i + 44], ja[i + 45], ja[i + 46], ja[i + 47],
                      ja[i + 48],
                      ja[i + 49], ja[i + 50], ja[i + 51], ja[i + 52], ja[i + 53], ja[i + 54], ja[i + 55], ja[i + 56],
                      ja[i + 57], ja[i + 58],
                      file=open("output4b.txt", "a+"), sep="")
                strix = open("output4b.txt")
                striy = strix.readlines()
                string = striy[-1]
                res = [sub.replace('_', '') for sub in string]
                x = "".join(res)
                print(x[:17], file=open("oligob.txt", "a+"))
            else:
                print(ja[i - 10], ja[i - 9], ja[i - 8], ja[i - 7], ja[i - 6], ja[i - 5], ja[i - 4], ja[i - 3],
                      ja[i - 2], ja[i - 1], ja[i], ja[i + 1], ja[i + 2], ja[i + 3], ja[i + 4], ja[i + 5], ja[i + 6],
                      file=open("oligob.txt", "a+"), sep="")
        if i == 11:
            if ja[i - 11] == "_" or ja[i - 10] == "_" or ja[i - 9] == "_" or ja[i - 8] == "_" or ja[i - 7] == "_" or ja[
                i - 6] == "_" or ja[i - 5] == "_" or ja[i - 4] == "_" or ja[i - 3] == "_" or ja[i - 2] == "_" or \
                    ja[i - 1] == "_" or ja[i] == "_" or ja[i + 1] == "_" or ja[i + 2] == "_" or ja[i + 3] == "_" or \
                    ja[i + 4] == "_" or ja[i + 5] == "_":
                print(ja[i - 11], ja[i - 10], ja[i - 9], ja[i - 8], ja[i - 7], ja[i - 6], ja[i - 5], ja[i - 4],
                      ja[i - 3],
                      ja[i - 2], ja[i - 1], ja[i], ja[i + 1], ja[i + 2], ja[i + 3],
                      ja[i + 4], ja[i + 5], ja[i + 6],
                      ja[i + 7], ja[i + 8], ja[i + 9], ja[i + 10], ja[i + 11], ja[i + 12], ja[i + 13], ja[i + 14],
                      ja[i + 15], ja[i + 16], ja[i + 17], ja[i + 18], ja[i + 19], ja[i + 20], ja[i + 21], ja[i + 22],
                      ja[i + 23], ja[i + 24], ja[i + 25], ja[i + 26], ja[i + 27], ja[i + 28], ja[i + 29], ja[i + 30],
                      ja[i + 31], ja[i + 32], ja[i + 33], ja[i + 34], ja[i + 35], ja[i + 36], ja[i + 37],
                      ja[i + 38], ja[i + 39],
                      ja[i + 40], ja[i + 41], ja[i + 42], ja[i + 43], ja[i + 44], ja[i + 45], ja[i + 46], ja[i + 47],
                      ja[i + 48],
                      ja[i + 49], ja[i + 50], ja[i + 51], ja[i + 52], ja[i + 53], ja[i + 54], ja[i + 55], ja[i + 56],
                      ja[i + 57],
                      file=open("output4b.txt", "a+"), sep="")
                strix = open("output4b.txt")
                striy = strix.readlines()
                string = striy[-1]
                res = [sub.replace('_', '') for sub in string]
                x = "".join(res)
                print(x[:17], file=open("oligob.txt", "a+"))
            else:
                print(ja[i - 11], ja[i - 10], ja[i - 9], ja[i - 8], ja[i - 7], ja[i - 6], ja[i - 5], ja[i - 4],
                      ja[i - 3],
                      ja[i - 2], ja[i - 1], ja[i], ja[i + 1], ja[i + 2], ja[i + 3], ja[i + 4], ja[i + 5],
                      file=open("oligob.txt", "a+"), sep="")
        if i == 12:
            if ja[i - 12] == "_" or ja[i - 11] == "_" or ja[i - 10] == "_" or ja[i - 9] == "_" or ja[i - 8] == "_" or \
                    ja[
                        i - 7] == "_" or ja[i - 6] == "_" or ja[i - 5] == "_" or ja[i - 4] == "_" or ja[i - 3] == "_" or \
                    ja[
                        i - 2] == "_" or \
                    ja[i - 1] == "_" or ja[i] == "_" or ja[i + 1] == "_" or ja[i + 2] == "_" or ja[i + 3] == "_" or \
                    ja[i + 4] == "_":
                print(ja[i - 12], ja[i - 11], ja[i - 10], ja[i - 9], ja[i - 8], ja[i - 7], ja[i - 6], ja[i - 5],
                      ja[i - 4],
                      ja[i - 3], ja[i - 2], ja[i - 1], ja[i], ja[i + 1], ja[i + 2], ja[i + 3],
                      ja[i + 4], ja[i + 5], ja[i + 6],
                      ja[i + 7], ja[i + 8], ja[i + 9], ja[i + 10], ja[i + 11], ja[i + 12], ja[i + 13], ja[i + 14],
                      ja[i + 15], ja[i + 16], ja[i + 17], ja[i + 18], ja[i + 19], ja[i + 20], ja[i + 21], ja[i + 22],
                      ja[i + 23], ja[i + 24], ja[i + 25], ja[i + 26], ja[i + 27], ja[i + 28], ja[i + 29], ja[i + 30],
                      ja[i + 31], ja[i + 32], ja[i + 33], ja[i + 34], ja[i + 35], ja[i + 36], ja[i + 37],
                      ja[i + 38], ja[i + 39],
                      ja[i + 40], ja[i + 41], ja[i + 42], ja[i + 43], ja[i + 44], ja[i + 45], ja[i + 46], ja[i + 47],
                      ja[i + 48],
                      ja[i + 49], ja[i + 50], ja[i + 51], ja[i + 52], ja[i + 53], ja[i + 54], ja[i + 55], ja[i + 56],
                      file=open("output4b.txt", "a+"), sep="")
                strix = open("output4b.txt")
                striy = strix.readlines()
                string = striy[-1]
                res = [sub.replace('_', '') for sub in string]
                x = "".join(res)
                print(x[:17], file=open("oligob.txt", "a+"))
            else:
                print(ja[i - 12], ja[i - 11], ja[i - 10], ja[i - 9], ja[i - 8], ja[i - 7], ja[i - 6], ja[i - 5],
                      ja[i - 4],
                      ja[i - 3],
                      ja[i - 2], ja[i - 1], ja[i], ja[i + 1], ja[i + 2], ja[i + 3], ja[i + 4],
                      file=open("oligob.txt", "a+"), sep="")
        if i == 13:
            if ja[i - 13] == "_" or ja[i - 12] == "_" or ja[i - 11] == "_" or ja[i - 10] == "_" or ja[i - 9] == "_" or \
                    ja[
                        i - 8] == "_" or ja[i - 7] == "_" or ja[i - 6] == "_" or ja[i - 5] == "_" or ja[i - 4] == "_" or \
                    ja[
                        i - 3] == "_" or ja[i - 2] == "_" or \
                    ja[i - 1] == "_" or ja[i] == "_" or ja[i + 1] == "_" or ja[i + 2] == "_" or ja[i + 3] == "_":
                print(ja[i - 13], ja[i - 12], ja[i - 11], ja[i - 10], ja[i - 9], ja[i - 8], ja[i - 7], ja[i - 6],
                      ja[i - 5],
                      ja[i - 4], ja[i - 3], ja[i - 2], ja[i - 1], ja[i], ja[i + 1], ja[i + 2], ja[i + 3],
                      ja[i + 4], ja[i + 5], ja[i + 6],
                      ja[i + 7], ja[i + 8], ja[i + 9], ja[i + 10], ja[i + 11], ja[i + 12], ja[i + 13], ja[i + 14],
                      ja[i + 15], ja[i + 16], ja[i + 17], ja[i + 18], ja[i + 19], ja[i + 20], ja[i + 21], ja[i + 22],
                      ja[i + 23], ja[i + 24], ja[i + 25], ja[i + 26], ja[i + 27], ja[i + 28], ja[i + 29], ja[i + 30],
                      ja[i + 31], ja[i + 32], ja[i + 33], ja[i + 34], ja[i + 35], ja[i + 36], ja[i + 37],
                      ja[i + 38], ja[i + 39],
                      ja[i + 40], ja[i + 41], ja[i + 42], ja[i + 43], ja[i + 44], ja[i + 45], ja[i + 46], ja[i + 47],
                      ja[i + 48],
                      ja[i + 49], ja[i + 50], ja[i + 51], ja[i + 52], ja[i + 53], ja[i + 54], ja[i + 55],
                      file=open("output4b.txt", "a+"), sep="")
                strix = open("output4b.txt")
                striy = strix.readlines()
                string = striy[-1]
                res = [sub.replace('_', '') for sub in string]
                x = "".join(res)
                print(x[:17], file=open("oligob.txt", "a+"))
            else:
                print(ja[i - 13], ja[i - 12], ja[i - 11], ja[i - 10], ja[i - 9], ja[i - 8], ja[i - 7], ja[i - 6],
                      ja[i - 5],
                      ja[i - 4], ja[i - 3],
                      ja[i - 2], ja[i - 1], ja[i], ja[i + 1], ja[i + 2], ja[i + 3],
                      file=open("oligob.txt", "a+"), sep="")
        if i == 14:
            if ja[i - 14] == "_" or ja[i - 13] == "_" or ja[i - 12] == "_" or ja[i - 11] == "_" or ja[i - 10] == "_" or \
                    ja[
                        i - 9] == "_" or ja[i - 8] == "_" or ja[i - 7] == "_" or ja[i - 6] == "_" or ja[i - 5] == "_" or \
                    ja[
                        i - 4] == "_" or ja[i - 3] == "_" or ja[i - 2] == "_" or \
                    ja[i - 1] == "_" or ja[i] == "_" or ja[i + 1] == "_" or ja[i + 2] == "_":
                print(ja[i - 14], ja[i - 13], ja[i - 12], ja[i - 11], ja[i - 10], ja[i - 9], ja[i - 8], ja[i - 7],
                      ja[i - 6], ja[i - 5], ja[i - 4], ja[i - 3], ja[i - 2], ja[i - 1], ja[i], ja[i + 1], ja[i + 2],
                      ja[i + 3],
                      ja[i + 4], ja[i + 5], ja[i + 6],
                      ja[i + 7], ja[i + 8], ja[i + 9], ja[i + 10], ja[i + 11], ja[i + 12], ja[i + 13], ja[i + 14],
                      ja[i + 15], ja[i + 16], ja[i + 17], ja[i + 18], ja[i + 19], ja[i + 20], ja[i + 21], ja[i + 22],
                      ja[i + 23], ja[i + 24], ja[i + 25], ja[i + 26], ja[i + 27], ja[i + 28], ja[i + 29], ja[i + 30],
                      ja[i + 31], ja[i + 32], ja[i + 33], ja[i + 34], ja[i + 35], ja[i + 36], ja[i + 37],
                      ja[i + 38], ja[i + 39],
                      ja[i + 40], ja[i + 41], ja[i + 42], ja[i + 43], ja[i + 44], ja[i + 45], ja[i + 46], ja[i + 47],
                      ja[i + 48],
                      ja[i + 49], ja[i + 50], ja[i + 51], ja[i + 52], ja[i + 53], ja[i + 54],
                      file=open("output4b.txt", "a+"), sep="")
                strix = open("output4b.txt")
                striy = strix.readlines()
                string = striy[-1]
                res = [sub.replace('_', '') for sub in string]
                x = "".join(res)
                print(x[:17], file=open("oligob.txt", "a+"))
            else:
                print(ja[i - 14], ja[i - 13], ja[i - 12], ja[i - 11], ja[i - 10], ja[i - 9], ja[i - 8], ja[i - 7],
                      ja[i - 6], ja[i - 5], ja[i - 4], ja[i - 3],
                      ja[i - 2], ja[i - 1], ja[i], ja[i + 1], ja[i + 2],
                      file=open("oligob.txt", "a+"), sep="")
        if i == 15:
            if ja[i - 15] == "_" or ja[i - 14] == "_" or ja[i - 13] == "_" or ja[i - 12] == "_" or ja[i - 11] == "_" or \
                    ja[i - 10] == "_" or ja[i - 9] == "_" or ja[i - 8] == "_" or ja[i - 7] == "_" or ja[i - 6] == "_" or \
                    ja[i - 5] == "_" or ja[i - 4] == "_" or ja[i - 3] == "_" or ja[i - 2] == "_" or \
                    ja[i - 1] == "_" or ja[i] == "_" or ja[i + 1] == "_":
                print(ja[i - 15], ja[i - 14], ja[i - 13], ja[i - 12], ja[i - 11], ja[i - 10], ja[i - 9], ja[i - 8],
                      ja[i - 7], ja[i - 6], ja[i - 5], ja[i - 4], ja[i - 3], ja[i - 2], ja[i - 1], ja[i], ja[i + 1],
                      ja[i + 2], ja[i + 3],
                      ja[i + 4], ja[i + 5], ja[i + 6],
                      ja[i + 7], ja[i + 8], ja[i + 9], ja[i + 10], ja[i + 11], ja[i + 12], ja[i + 13], ja[i + 14],
                      ja[i + 15], ja[i + 16], ja[i + 17], ja[i + 18], ja[i + 19], ja[i + 20], ja[i + 21], ja[i + 22],
                      ja[i + 23], ja[i + 24], ja[i + 25], ja[i + 26], ja[i + 27], ja[i + 28], ja[i + 29], ja[i + 30],
                      ja[i + 31], ja[i + 32], ja[i + 33], ja[i + 34], ja[i + 35], ja[i + 36], ja[i + 37],
                      ja[i + 38], ja[i + 39],
                      ja[i + 40], ja[i + 41], ja[i + 42], ja[i + 43], ja[i + 44], ja[i + 45], ja[i + 46], ja[i + 47],
                      ja[i + 48],
                      ja[i + 49], ja[i + 50], ja[i + 51], ja[i + 52], ja[i + 53], file=open("output4b.txt", "a+"),
                      sep="")
                strix = open("output4b.txt")
                striy = strix.readlines()
                string = striy[-1]
                res = [sub.replace('_', '') for sub in string]
                x = "".join(res)
                print(x[:17], file=open("oligob.txt", "a+"))
            else:
                print(ja[i - 15], ja[i - 14], ja[i - 13], ja[i - 12], ja[i - 11], ja[i - 10], ja[i - 9], ja[i - 8],
                      ja[i - 7], ja[i - 6], ja[i - 5], ja[i - 4], ja[i - 3],
                      ja[i - 2], ja[i - 1], ja[i], ja[i + 1],
                      file=open("oligob.txt", "a+"), sep="")
        if i == 16:
            if ja[i - 16] == "_" or ja[i - 15] == "_" or ja[i - 14] == "_" or ja[i - 13] == "_" or ja[i - 12] == "_" or \
                    ja[i - 11] == "_" or ja[i - 10] == "_" or ja[i - 9] == "_" or ja[i - 8] == "_" or \
                    ja[i - 7] == "_" or ja[i - 6] == "_" or ja[i - 5] == "_" or ja[i - 4] == "_" or ja[i - 3] == "_" or \
                    ja[i - 2] == "_" or \
                    ja[i - 1] == "_" or ja[i] == "_":
                print(ja[i - 16], ja[i - 15], ja[i - 14], ja[i - 13], ja[i - 12], ja[i - 11], ja[i - 10], ja[i - 9],
                      ja[i - 8], ja[i - 7], ja[i - 6], ja[i - 5], ja[i - 4], ja[i - 3], ja[i - 2], ja[i - 1], ja[i],
                      ja[i + 1], ja[i + 2], ja[i + 3],
                      ja[i + 4], ja[i + 5], ja[i + 6],
                      ja[i + 7], ja[i + 8], ja[i + 9], ja[i + 10], ja[i + 11], ja[i + 12], ja[i + 13], ja[i + 14],
                      ja[i + 15], ja[i + 16], ja[i + 17], ja[i + 18], ja[i + 19], ja[i + 20], ja[i + 21], ja[i + 22],
                      ja[i + 23], ja[i + 24], ja[i + 25], ja[i + 26], ja[i + 27], ja[i + 28], ja[i + 29], ja[i + 30],
                      ja[i + 31], ja[i + 32], ja[i + 33], ja[i + 34], ja[i + 35], ja[i + 36], ja[i + 37],
                      ja[i + 38], ja[i + 39],
                      ja[i + 40], ja[i + 41], ja[i + 42], ja[i + 43], ja[i + 44], ja[i + 45], ja[i + 46], ja[i + 47],
                      ja[i + 48],
                      ja[i + 49], ja[i + 50], ja[i + 51], ja[i + 52], file=open("output4b.txt", "a+"), sep="")
                strix = open("output4b.txt")
                striy = strix.readlines()
                string = striy[-1]
                res = [sub.replace('_', '') for sub in string]
                x = "".join(res)
                print(x[:17], file=open("oligob.txt", "a+"))
            else:
                print(ja[i - 16], ja[i - 15], ja[i - 14], ja[i - 13], ja[i - 12], ja[i - 11], ja[i - 10], ja[i - 9],
                      ja[i - 8], ja[i - 7], ja[i - 6], ja[i - 5], ja[i - 4], ja[i - 3],
                      ja[i - 2], ja[i - 1], ja[i], file=open("oligob.txt", "a+"), sep="")
        if i == 17:
            if ja[i] == "_" or ja[i + 1] == "_" or ja[i + 2] == "_" or ja[i + 3] == "_" or ja[i + 4] == "_" or ja[i + 5] \
                    == "_" or ja[i + 6] == "_" or ja[i + 7] == "_" or ja[i + 8] == "_" or ja[i + 9] == "_" or ja[i + 10] \
                    == "_" or ja[i + 11] == "_" or ja[i + 12] == "_" or ja[i + 13] == "_" or ja[i + 14] == "_" or \
                    ja[i + 15] == "_" or ja[i + 16] == "_":
                print(ja[i], ja[i + 1], ja[i + 2], ja[i + 3], ja[i + 4], ja[i + 5], ja[i + 6], ja[i + 7], ja[i + 8],
                      ja[i + 9], ja[i + 10], ja[i + 11], ja[i + 12], ja[i + 13], ja[i + 14], ja[i + 15], ja[i + 16],
                      ja[i + 17], ja[i + 18], ja[i + 19], ja[i + 20], ja[i + 21], ja[i + 22], ja[i + 23], ja[i + 24],
                      ja[i + 25], ja[i + 26], ja[i + 27], ja[i + 28], ja[i + 29], ja[i + 30],
                      ja[i + 31], ja[i + 32], ja[i + 33], ja[i + 34], ja[i + 35], ja[i + 36], ja[i + 37], ja[i + 38],
                      ja[i + 39], ja[i + 40], ja[i + 41], ja[i + 42], ja[i + 43], ja[i + 44], ja[i + 45], ja[i + 46],
                      ja[i + 47], ja[i + 48], ja[i + 49], ja[i + 50], ja[i + 51], ja[i + 52], ja[i + 53],
                      ja[i + 54], ja[i + 55], ja[i + 56], ja[i + 57], ja[i + 58], ja[i + 59], ja[i + 60], ja[i + 61],
                      ja[i + 62], ja[i + 63], ja[i + 64], ja[i + 65], ja[i + 66], ja[i + 67], ja[i + 68],
                      file=open("output4b.txt", "a+"), sep="")
                strix = open("output4b.txt")
                striy = strix.readlines()
                string = striy[-1]
                res = [sub.replace('_', '') for sub in string]
                x = "".join(res)
                print(x[:17], file=open("oligob.txt", "a+"))
            else:
                print(ja[i], ja[i + 1], ja[i + 2], ja[i + 3], ja[i + 4], ja[i + 5], ja[i + 6], ja[i + 7], ja[i + 8],
                      ja[i + 9], ja[i + 10], ja[i + 11], ja[i + 12], ja[i + 13], ja[i + 14], ja[i + 15], ja[i + 16],
                      file=open("oligob.txt", "a+"), sep="")
        if i == 18:
            if ja[i - 1] == "_" or ja[i] == "_" or ja[i + 1] == "_" or ja[i + 2] == "_" or ja[i + 3] == "_" or \
                    ja[i + 4] == "_" or ja[i + 5] \
                    == "_" or ja[i + 6] == "_" or ja[i + 7] == "_" or ja[i + 8] == "_" or ja[i + 9] == "_" or \
                    ja[i + 10] == "_" or ja[i + 11] == "_" or ja[i + 12] == "_" or ja[i + 13] == "_" or \
                    ja[i + 14] == "_" or ja[i + 15] == "_":
                print(ja[i - 1], ja[i], ja[i + 1], ja[i + 2], ja[i + 3], ja[i + 4], ja[i + 5], ja[i + 6], ja[i + 7],
                      ja[i + 8],
                      ja[i + 9], ja[i + 10], ja[i + 11], ja[i + 12], ja[i + 13], ja[i + 14], ja[i + 15], ja[i + 16],
                      ja[i + 17], ja[i + 18], ja[i + 19], ja[i + 20], ja[i + 21], ja[i + 22], ja[i + 23], ja[i + 24],
                      ja[i + 25], ja[i + 26], ja[i + 27], ja[i + 28], ja[i + 29], ja[i + 30],
                      ja[i + 31], ja[i + 32], ja[i + 33], ja[i + 34], ja[i + 35], ja[i + 36], ja[i + 37], ja[i + 38],
                      ja[i + 39], ja[i + 40], ja[i + 41], ja[i + 42], ja[i + 43], ja[i + 44], ja[i + 45], ja[i + 46],
                      ja[i + 47], ja[i + 48], ja[i + 49], ja[i + 50], ja[i + 51], ja[i + 52], ja[i + 53],
                      ja[i + 54], ja[i + 55], ja[i + 56], ja[i + 57], ja[i + 58], ja[i + 59], ja[i + 60], ja[i + 61],
                      ja[i + 62], ja[i + 63], ja[i + 64], ja[i + 65], ja[i + 66], ja[i + 67],
                      file=open("output4b.txt", "a+"), sep="")
                strix = open("output4b.txt")
                striy = strix.readlines()
                string = striy[-1]
                res = [sub.replace('_', '') for sub in string]
                x = "".join(res)
                print(x[:17], file=open("oligob.txt", "a+"))
            else:
                print(ja[i - 1], ja[i], ja[i + 1], ja[i + 2], ja[i + 3], ja[i + 4], ja[i + 5], ja[i + 6], ja[i + 7],
                      ja[i + 8],
                      ja[i + 9], ja[i + 10], ja[i + 11], ja[i + 12], ja[i + 13], ja[i + 14], ja[i + 15],
                      file=open("oligob.txt", "a+"), sep="")
        if i == 19:
            if ja[i - 2] == "_" or ja[i - 1] == "_" or ja[i] == "_" or ja[i + 1] == "_" or ja[i + 2] == "_" or \
                    ja[i + 3] == "_" or ja[i + 4] == "_" or ja[i + 5] == "_" or ja[i + 6] == "_" or ja[i + 7] == "_" or \
                    ja[i + 8] == "_" or ja[i + 9] == "_" or ja[i + 10] \
                    == "_" or ja[i + 11] == "_" or ja[i + 12] == "_" or ja[i + 13] == "_" or ja[i + 14] == "_":
                print(ja[i - 2], ja[i - 1], ja[i], ja[i + 1], ja[i + 2], ja[i + 3], ja[i + 4], ja[i + 5], ja[i + 6],
                      ja[i + 7], ja[i + 8],
                      ja[i + 9], ja[i + 10], ja[i + 11], ja[i + 12], ja[i + 13], ja[i + 14], ja[i + 15], ja[i + 16],
                      ja[i + 17], ja[i + 18], ja[i + 19], ja[i + 20], ja[i + 21], ja[i + 22], ja[i + 23], ja[i + 24],
                      ja[i + 25], ja[i + 26], ja[i + 27], ja[i + 28], ja[i + 29], ja[i + 30],
                      ja[i + 31], ja[i + 32], ja[i + 33], ja[i + 34], ja[i + 35], ja[i + 36], ja[i + 37], ja[i + 38],
                      ja[i + 39], ja[i + 40], ja[i + 41], ja[i + 42], ja[i + 43], ja[i + 44], ja[i + 45], ja[i + 46],
                      ja[i + 47], ja[i + 48], ja[i + 49], ja[i + 50], ja[i + 51], ja[i + 52], ja[i + 53],
                      ja[i + 54], ja[i + 55], ja[i + 56], ja[i + 57], ja[i + 58], ja[i + 59], ja[i + 60], ja[i + 61],
                      ja[i + 62], ja[i + 63], ja[i + 64], ja[i + 65], ja[i + 66],
                      file=open("output4b.txt", "a+"), sep="")
                strix = open("output4b.txt")
                striy = strix.readlines()
                string = striy[-1]
                res = [sub.replace('_', '') for sub in string]
                x = "".join(res)
                print(x[:17], file=open("oligob.txt", "a+"))
            else:
                print(ja[i - 2], ja[i - 1], ja[i], ja[i + 1], ja[i + 2], ja[i + 3], ja[i + 4], ja[i + 5], ja[i + 6],
                      ja[i + 7], ja[i + 8],
                      ja[i + 9], ja[i + 10], ja[i + 11], ja[i + 12], ja[i + 13], ja[i + 14],
                      file=open("oligob.txt", "a+"), sep="")
        if i == 20:
            if ja[i - 3] == "_" or ja[i - 2] == "_" or ja[i - 1] == "_" or ja[i] == "_" or ja[i + 1] == "_" or \
                    ja[i + 2] == "_" or ja[i + 3] == "_" or ja[i + 4] == "_" or ja[i + 5] \
                    == "_" or ja[i + 6] == "_" or ja[i + 7] == "_" or ja[i + 8] == "_" or ja[i + 9] == "_" or ja[i + 10] \
                    == "_" or ja[i + 11] == "_" or ja[i + 12] == "_" or ja[i + 13] == "_":
                print(ja[i - 3], ja[i - 2], ja[i - 1], ja[i], ja[i + 1], ja[i + 2], ja[i + 3], ja[i + 4], ja[i + 5],
                      ja[i + 6], ja[i + 7], ja[i + 8],
                      ja[i + 9], ja[i + 10], ja[i + 11], ja[i + 12], ja[i + 13], ja[i + 14], ja[i + 15], ja[i + 16],
                      ja[i + 17], ja[i + 18], ja[i + 19], ja[i + 20], ja[i + 21], ja[i + 22], ja[i + 23], ja[i + 24],
                      ja[i + 25], ja[i + 26], ja[i + 27], ja[i + 28], ja[i + 29], ja[i + 30],
                      ja[i + 31], ja[i + 32], ja[i + 33], ja[i + 34], ja[i + 35], ja[i + 36], ja[i + 37], ja[i + 38],
                      ja[i + 39], ja[i + 40], ja[i + 41], ja[i + 42], ja[i + 43], ja[i + 44], ja[i + 45], ja[i + 46],
                      ja[i + 47], ja[i + 48], ja[i + 49], ja[i + 50], ja[i + 51], ja[i + 52], ja[i + 53],
                      ja[i + 54], ja[i + 55], ja[i + 56], ja[i + 57], ja[i + 58], ja[i + 59], ja[i + 60], ja[i + 61],
                      ja[i + 62], ja[i + 63], ja[i + 64], ja[i + 65],
                      file=open("output4b.txt", "a+"), sep="")
                strix = open("output4b.txt")
                striy = strix.readlines()
                string = striy[-1]
                res = [sub.replace('_', '') for sub in string]
                x = "".join(res)
                print(x[:17], file=open("oligob.txt", "a+"))
            else:
                print(ja[i - 3], ja[i - 2], ja[i - 1], ja[i], ja[i + 1], ja[i + 2], ja[i + 3], ja[i + 4], ja[i + 5],
                      ja[i + 6], ja[i + 7], ja[i + 8],
                      ja[i + 9], ja[i + 10], ja[i + 11], ja[i + 12], ja[i + 13],
                      file=open("oligob.txt", "a+"), sep="")
        if i == 21:
            if ja[i - 4] == "_" or ja[i - 3] == "_" or ja[i - 2] == "_" or ja[i - 1] == "_" or ja[i] == "_" or \
                    ja[i + 1] == "_" or ja[i + 2] == "_" or ja[i + 3] == "_" or ja[i + 4] == "_" or ja[i + 5] \
                    == "_" or ja[i + 6] == "_" or ja[i + 7] == "_" or ja[i + 8] == "_" or ja[i + 9] == "_" or \
                    ja[i + 10] == "_" or ja[i + 11] == "_" or ja[i + 12] == "_":
                print(ja[i - 4], ja[i - 3], ja[i - 2], ja[i - 1], ja[i], ja[i + 1], ja[i + 2], ja[i + 3], ja[i + 4],
                      ja[i + 5], ja[i + 6], ja[i + 7], ja[i + 8],
                      ja[i + 9], ja[i + 10], ja[i + 11], ja[i + 12], ja[i + 13], ja[i + 14], ja[i + 15], ja[i + 16],
                      ja[i + 17], ja[i + 18], ja[i + 19], ja[i + 20], ja[i + 21], ja[i + 22], ja[i + 23], ja[i + 24],
                      ja[i + 25], ja[i + 26], ja[i + 27], ja[i + 28], ja[i + 29], ja[i + 30],
                      ja[i + 31], ja[i + 32], ja[i + 33], ja[i + 34], ja[i + 35], ja[i + 36], ja[i + 37], ja[i + 38],
                      ja[i + 39], ja[i + 40], ja[i + 41], ja[i + 42], ja[i + 43], ja[i + 44], ja[i + 45], ja[i + 46],
                      ja[i + 47], ja[i + 48], ja[i + 49], ja[i + 50], ja[i + 51], ja[i + 52], ja[i + 53],
                      ja[i + 54], ja[i + 55], ja[i + 56], ja[i + 57], ja[i + 58], ja[i + 59], ja[i + 60], ja[i + 61],
                      ja[i + 62], ja[i + 63], ja[i + 64],
                      file=open("output4b.txt", "a+"), sep="")
                strix = open("output4b.txt")
                striy = strix.readlines()
                string = striy[-1]
                res = [sub.replace('_', '') for sub in string]
                x = "".join(res)
                print(x[:17], file=open("oligob.txt", "a+"))
            else:
                print(ja[i - 4], ja[i - 3], ja[i - 2], ja[i - 1], ja[i], ja[i + 1], ja[i + 2], ja[i + 3], ja[i + 4],
                      ja[i + 5], ja[i + 6], ja[i + 7], ja[i + 8],
                      ja[i + 9], ja[i + 10], ja[i + 11], ja[i + 12],
                      file=open("oligob.txt", "a+"), sep="")
        if i == 22:
            if ja[i - 5] == "_" or ja[i - 4] == "_" or ja[i - 3] == "_" or ja[i - 2] == "_" or ja[i - 1] == "_" or \
                    ja[i] == "_" or ja[i + 1] == "_" or ja[i + 2] == "_" or ja[i + 3] == "_" or ja[i + 4] == "_" or \
                    ja[i + 5] == "_" or ja[i + 6] == "_" or ja[i + 7] == "_" or ja[i + 8] == "_" or ja[i + 9] == "_" or \
                    ja[i + 10] == "_" or ja[i + 11] == "_":
                print(ja[i - 5], ja[i - 4], ja[i - 3], ja[i - 2], ja[i - 1], ja[i], ja[i + 1], ja[i + 2], ja[i + 3],
                      ja[i + 4], ja[i + 5], ja[i + 6], ja[i + 7], ja[i + 8],
                      ja[i + 9], ja[i + 10], ja[i + 11], ja[i + 12], ja[i + 13], ja[i + 14], ja[i + 15], ja[i + 16],
                      ja[i + 17], ja[i + 18], ja[i + 19], ja[i + 20], ja[i + 21], ja[i + 22], ja[i + 23], ja[i + 24],
                      ja[i + 25], ja[i + 26], ja[i + 27], ja[i + 28], ja[i + 29], ja[i + 30],
                      ja[i + 31], ja[i + 32], ja[i + 33], ja[i + 34], ja[i + 35], ja[i + 36], ja[i + 37], ja[i + 38],
                      ja[i + 39], ja[i + 40], ja[i + 41], ja[i + 42], ja[i + 43], ja[i + 44], ja[i + 45], ja[i + 46],
                      ja[i + 47], ja[i + 48], ja[i + 49], ja[i + 50], ja[i + 51], ja[i + 52], ja[i + 53],
                      ja[i + 54], ja[i + 55], ja[i + 56], ja[i + 57], ja[i + 58], ja[i + 59], ja[i + 60], ja[i + 61],
                      ja[i + 62], ja[i + 63],
                      file=open("output4b.txt", "a+"), sep="")
                strix = open("output4b.txt")
                striy = strix.readlines()
                string = striy[-1]
                res = [sub.replace('_', '') for sub in string]
                x = "".join(res)
                print(x[:17], file=open("oligob.txt", "a+"))
            else:
                print(ja[i - 5], ja[i - 4], ja[i - 3], ja[i - 2], ja[i - 1], ja[i], ja[i + 1], ja[i + 2], ja[i + 3],
                      ja[i + 4], ja[i + 5], ja[i + 6], ja[i + 7], ja[i + 8], ja[i + 9], ja[i + 10], ja[i + 11],
                      file=open("oligob.txt", "a+"), sep="")
        if i == 23:
            if ja[i - 6] == "_" or ja[i - 5] == "_" or ja[i - 4] == "_" or ja[i - 3] == "_" or ja[i - 2] == "_" or \
                    ja[i - 1] == "_" or ja[i] == "_" or ja[i + 1] == "_" or ja[i + 2] == "_" or ja[i + 3] == "_" or \
                    ja[i + 4] == "_" or ja[i + 5] \
                    == "_" or ja[i + 6] == "_" or ja[i + 7] == "_" or ja[i + 8] == "_" or ja[i + 9] == "_" or \
                    ja[i + 10] == "_":
                print(ja[i - 6], ja[i - 5], ja[i - 4], ja[i - 3], ja[i - 2], ja[i - 1], ja[i], ja[i + 1], ja[i + 2],
                      ja[i + 3], ja[i + 4], ja[i + 5], ja[i + 6], ja[i + 7], ja[i + 8],
                      ja[i + 9], ja[i + 10], ja[i + 11], ja[i + 12], ja[i + 13], ja[i + 14], ja[i + 15], ja[i + 16],
                      ja[i + 17], ja[i + 18], ja[i + 19], ja[i + 20], ja[i + 21], ja[i + 22], ja[i + 23], ja[i + 24],
                      ja[i + 25], ja[i + 26], ja[i + 27], ja[i + 28], ja[i + 29], ja[i + 30],
                      ja[i + 31], ja[i + 32], ja[i + 33], ja[i + 34], ja[i + 35], ja[i + 36], ja[i + 37], ja[i + 38],
                      ja[i + 39], ja[i + 40], ja[i + 41], ja[i + 42], ja[i + 43], ja[i + 44], ja[i + 45], ja[i + 46],
                      ja[i + 47], ja[i + 48], ja[i + 49], ja[i + 50], ja[i + 51], ja[i + 52], ja[i + 53],
                      ja[i + 54], ja[i + 55], ja[i + 56], ja[i + 57], ja[i + 58], ja[i + 59], ja[i + 60], ja[i + 61],
                      ja[i + 62],
                      file=open("output4b.txt", "a+"), sep="")
                strix = open("output4b.txt")
                striy = strix.readlines()
                string = striy[-1]
                res = [sub.replace('_', '') for sub in string]
                x = "".join(res)
                print(x[:17], file=open("oligob.txt", "a+"))
            else:
                print(ja[i - 6], ja[i - 5], ja[i - 4], ja[i - 3], ja[i - 2], ja[i - 1], ja[i], ja[i + 1], ja[i + 2],
                      ja[i + 3], ja[i + 4], ja[i + 5], ja[i + 6], ja[i + 7], ja[i + 8],
                      ja[i + 9], ja[i + 10],
                      file=open("oligob.txt", "a+"), sep="")
        if i == 24:
            if ja[i - 7] == "_" or ja[i - 6] == "_" or ja[i - 5] == "_" or ja[i - 4] == "_" or ja[i - 3] == "_" or \
                    ja[i - 2] == "_" or ja[i - 1] == "_" or ja[i] == "_" or ja[i + 1] == "_" or ja[i + 2] == "_" or \
                    ja[i + 3] == "_" or ja[i + 4] == "_" or ja[i + 5] \
                    == "_" or ja[i + 6] == "_" or ja[i + 7] == "_" or ja[i + 8] == "_" or ja[i + 9] == "_":
                print(ja[i - 7], ja[i - 6], ja[i - 5], ja[i - 4], ja[i - 3], ja[i - 2], ja[i - 1], ja[i], ja[i + 1],
                      ja[i + 2], ja[i + 3], ja[i + 4], ja[i + 5], ja[i + 6], ja[i + 7], ja[i + 8],
                      ja[i + 9], ja[i + 10], ja[i + 11], ja[i + 12], ja[i + 13], ja[i + 14], ja[i + 15], ja[i + 16],
                      ja[i + 17], ja[i + 18], ja[i + 19], ja[i + 20], ja[i + 21], ja[i + 22], ja[i + 23], ja[i + 24],
                      ja[i + 25], ja[i + 26], ja[i + 27], ja[i + 28], ja[i + 29], ja[i + 30],
                      ja[i + 31], ja[i + 32], ja[i + 33], ja[i + 34], ja[i + 35], ja[i + 36], ja[i + 37], ja[i + 38],
                      ja[i + 39], ja[i + 40], ja[i + 41], ja[i + 42], ja[i + 43], ja[i + 44], ja[i + 45], ja[i + 46],
                      ja[i + 47], ja[i + 48], ja[i + 49], ja[i + 50], ja[i + 51], ja[i + 52], ja[i + 53],
                      ja[i + 54], ja[i + 55], ja[i + 56], ja[i + 57], ja[i + 58], ja[i + 59], ja[i + 60], ja[i + 61],
                      file=open("output4b.txt", "a+"), sep="")
                strix = open("output4b.txt")
                striy = strix.readlines()
                string = striy[-1]
                res = [sub.replace('_', '') for sub in string]
                x = "".join(res)
                print(x[:17], file=open("oligob.txt", "a+"))
            else:
                print(ja[i - 7], ja[i - 6], ja[i - 5], ja[i - 4], ja[i - 3], ja[i - 2], ja[i - 1], ja[i], ja[i + 1],
                      ja[i + 2], ja[i + 3], ja[i + 4], ja[i + 5], ja[i + 6], ja[i + 7], ja[i + 8],
                      ja[i + 9], file=open("oligob.txt", "a+"), sep="")
        if i == 25:
            if ja[i - 8] == "_" or ja[i - 7] == "_" or ja[i - 6] == "_" or ja[i - 5] == "_" or ja[i - 4] == "_" or \
                    ja[i - 3] == "_" or ja[i - 2] == "_" or ja[i - 1] == "_" or ja[i] == "_" or ja[i + 1] == "_" or \
                    ja[i + 2] == "_" or ja[i + 3] == "_" or ja[i + 4] == "_" or ja[i + 5] \
                    == "_" or ja[i + 6] == "_" or ja[i + 7] == "_" or ja[i + 8] == "_":
                print(ja[i - 8], ja[i - 7], ja[i - 6], ja[i - 5], ja[i - 4], ja[i - 3], ja[i - 2], ja[i - 1], ja[i],
                      ja[i + 1], ja[i + 2], ja[i + 3], ja[i + 4], ja[i + 5], ja[i + 6], ja[i + 7], ja[i + 8],
                      ja[i + 9], ja[i + 10], ja[i + 11], ja[i + 12], ja[i + 13], ja[i + 14], ja[i + 15], ja[i + 16],
                      ja[i + 17], ja[i + 18], ja[i + 19], ja[i + 20], ja[i + 21], ja[i + 22], ja[i + 23], ja[i + 24],
                      ja[i + 25], ja[i + 26], ja[i + 27], ja[i + 28], ja[i + 29], ja[i + 30],
                      ja[i + 31], ja[i + 32], ja[i + 33], ja[i + 34], ja[i + 35], ja[i + 36], ja[i + 37], ja[i + 38],
                      ja[i + 39], ja[i + 40], ja[i + 41], ja[i + 42], ja[i + 43], ja[i + 44], ja[i + 45], ja[i + 46],
                      ja[i + 47], ja[i + 48], ja[i + 49], ja[i + 50], ja[i + 51], ja[i + 52], ja[i + 53],
                      ja[i + 54], ja[i + 55], ja[i + 56], ja[i + 57], ja[i + 58], ja[i + 59], ja[i + 60],
                      file=open("output4b.txt", "a+"), sep="")
                strix = open("output4b.txt")
                striy = strix.readlines()
                string = striy[-1]
                res = [sub.replace('_', '') for sub in string]
                x = "".join(res)
                print(x[:17], file=open("oligob.txt", "a+"))
            else:
                print(ja[i - 8], ja[i - 7], ja[i - 6], ja[i - 5], ja[i - 4], ja[i - 3], ja[i - 2], ja[i - 1], ja[i],
                      ja[i + 1], ja[i + 2], ja[i + 3], ja[i + 4], ja[i + 5], ja[i + 6], ja[i + 7], ja[i + 8],
                      file=open("oligob.txt", "a+"), sep="")

        if i == 26:
            if ja[i - 9] == "_" or ja[i - 8] == "_" or ja[i - 7] == "_" or ja[i - 6] == "_" or ja[i - 5] == "_" or \
                    ja[i - 4] == "_" or ja[i - 3] == "_" or ja[i - 2] == "_" or ja[i - 1] == "_" or ja[i] == "_" or \
                    ja[i + 1] == "_" or ja[i + 2] == "_" or ja[i + 3] == "_" or ja[i + 4] == "_" or ja[i + 5] \
                    == "_" or ja[i + 6] == "_" or ja[i + 7] == "_":
                print(ja[i - 9], ja[i - 8], ja[i - 7], ja[i - 6], ja[i - 5], ja[i - 4], ja[i - 3], ja[i - 2], ja[i - 1],
                      ja[i], ja[i + 1], ja[i + 2], ja[i + 3], ja[i + 4], ja[i + 5], ja[i + 6], ja[i + 7], ja[i + 8],
                      ja[i + 9], ja[i + 10], ja[i + 11], ja[i + 12], ja[i + 13], ja[i + 14], ja[i + 15], ja[i + 16],
                      ja[i + 17], ja[i + 18], ja[i + 19], ja[i + 20], ja[i + 21], ja[i + 22], ja[i + 23], ja[i + 24],
                      ja[i + 25], ja[i + 26], ja[i + 27], ja[i + 28], ja[i + 29], ja[i + 30],
                      ja[i + 31], ja[i + 32], ja[i + 33], ja[i + 34], ja[i + 35], ja[i + 36], ja[i + 37], ja[i + 38],
                      ja[i + 39], ja[i + 40], ja[i + 41], ja[i + 42], ja[i + 43], ja[i + 44], ja[i + 45], ja[i + 46],
                      ja[i + 47], ja[i + 48], ja[i + 49], ja[i + 50], ja[i + 51], ja[i + 52], ja[i + 53],
                      ja[i + 54], ja[i + 55], ja[i + 56], ja[i + 57], ja[i + 58], ja[i + 59],
                      file=open("output4b.txt", "a+"), sep="")
                strix = open("output4b.txt")
                striy = strix.readlines()
                string = striy[-1]
                res = [sub.replace('_', '') for sub in string]
                x = "".join(res)
                print(x[:17], file=open("oligob.txt", "a+"))
            else:
                print(ja[i - 9], ja[i - 8], ja[i - 7], ja[i - 6], ja[i - 5], ja[i - 4], ja[i - 3], ja[i - 2], ja[i - 1],
                      ja[i], ja[i + 1], ja[i + 2], ja[i + 3], ja[i + 4], ja[i + 5], ja[i + 6], ja[i + 7],
                      file=open("oligob.txt", "a+"), sep="")
        if i == 27:
            if ja[i - 10] == "_" or ja[i - 9] == "_" or ja[i - 8] == "_" or ja[i - 7] == "_" or ja[i - 6] == "_" or \
                    ja[i - 5] == "_" or ja[i - 4] == "_" or ja[i - 3] == "_" or ja[i - 2] == "_" or ja[i - 1] == "_" or \
                    ja[i] == "_" or ja[i + 1] == "_" or ja[i + 2] == "_" or ja[i + 3] == "_" or ja[i + 4] == "_" or \
                    ja[i + 5] == "_" or ja[i + 6] == "_":
                print(ja[i - 10], ja[i - 9], ja[i - 8], ja[i - 7], ja[i - 6], ja[i - 5], ja[i - 4], ja[i - 3],
                      ja[i - 2], ja[i - 1],
                      ja[i], ja[i + 1], ja[i + 2], ja[i + 3], ja[i + 4], ja[i + 5], ja[i + 6], ja[i + 7], ja[i + 8],
                      ja[i + 9], ja[i + 10], ja[i + 11], ja[i + 12], ja[i + 13], ja[i + 14], ja[i + 15], ja[i + 16],
                      ja[i + 17], ja[i + 18], ja[i + 19], ja[i + 20], ja[i + 21], ja[i + 22], ja[i + 23], ja[i + 24],
                      ja[i + 25], ja[i + 26], ja[i + 27], ja[i + 28], ja[i + 29], ja[i + 30],
                      ja[i + 31], ja[i + 32], ja[i + 33], ja[i + 34], ja[i + 35], ja[i + 36], ja[i + 37], ja[i + 38],
                      ja[i + 39], ja[i + 40], ja[i + 41], ja[i + 42], ja[i + 43], ja[i + 44], ja[i + 45], ja[i + 46],
                      ja[i + 47], ja[i + 48], ja[i + 49], ja[i + 50], ja[i + 51], ja[i + 52], ja[i + 53],
                      ja[i + 54], ja[i + 55], ja[i + 56], ja[i + 57], ja[i + 58],
                      file=open("output4b.txt", "a+"), sep="")
                strix = open("output4b.txt")
                striy = strix.readlines()
                string = striy[-1]
                res = [sub.replace('_', '') for sub in string]
                x = "".join(res)
                print(x[:17], file=open("oligob.txt", "a+"))
            else:
                print(ja[i - 10], ja[i - 9], ja[i - 8], ja[i - 7], ja[i - 6], ja[i - 5], ja[i - 4], ja[i - 3],
                      ja[i - 2], ja[i - 1],
                      ja[i], ja[i + 1], ja[i + 2], ja[i + 3], ja[i + 4], ja[i + 5], ja[i + 6],
                      file=open("oligob.txt", "a+"), sep="")
        if i == 28:
            if ja[i - 11] == "_" or ja[i - 10] == "_" or ja[i - 9] == "_" or ja[i - 8] == "_" or ja[i - 7] == "_" or \
                    ja[i - 6] == "_" or \
                    ja[i - 5] == "_" or ja[i - 4] == "_" or ja[i - 3] == "_" or ja[i - 2] == "_" or ja[i - 1] == "_" or \
                    ja[i] == "_" or ja[i + 1] == "_" or ja[i + 2] == "_" or ja[i + 3] == "_" or ja[i + 4] == "_" or \
                    ja[i + 5] == "_":
                print(ja[i - 11], ja[i - 9], ja[i - 8], ja[i - 7], ja[i - 6], ja[i - 5], ja[i - 4], ja[i - 3],
                      ja[i - 2], ja[i - 1],
                      ja[i], ja[i + 1], ja[i + 2], ja[i + 3], ja[i + 4], ja[i + 5], ja[i + 6], ja[i + 7], ja[i + 8],
                      ja[i + 9], ja[i + 10], ja[i + 11], ja[i + 12], ja[i + 13], ja[i + 14], ja[i + 15], ja[i + 16],
                      ja[i + 17], ja[i + 18], ja[i + 19], ja[i + 20], ja[i + 21], ja[i + 22], ja[i + 23], ja[i + 24],
                      ja[i + 25], ja[i + 26], ja[i + 27], ja[i + 28], ja[i + 29], ja[i + 30],
                      ja[i + 31], ja[i + 32], ja[i + 33], ja[i + 34], ja[i + 35], ja[i + 36], ja[i + 37], ja[i + 38],
                      ja[i + 39], ja[i + 40], ja[i + 41], ja[i + 42], ja[i + 43], ja[i + 44], ja[i + 45], ja[i + 46],
                      ja[i + 47], ja[i + 48], ja[i + 49], ja[i + 50], ja[i + 51], ja[i + 52], ja[i + 53],
                      ja[i + 54], ja[i + 55], ja[i + 56], ja[i + 57],
                      file=open("output4b.txt", "a+"), sep="")
                strix = open("output4b.txt")
                striy = strix.readlines()
                string = striy[-1]
                res = [sub.replace('_', '') for sub in string]
                x = "".join(res)
                print(x[:17], file=open("oligob.txt", "a+"))
            else:
                print(ja[i - 11], ja[i - 10], ja[i - 9], ja[i - 8], ja[i - 7], ja[i - 6], ja[i - 5], ja[i - 4],
                      ja[i - 3],
                      ja[i - 2], ja[i - 1],
                      ja[i], ja[i + 1], ja[i + 2], ja[i + 3], ja[i + 4], ja[i + 5],
                      file=open("oligob.txt", "a+"), sep="")
        if i == 29:
            if ja[i - 12] == "_" or ja[i - 11] == "_" or ja[i - 10] == "_" or ja[i - 9] == "_" or ja[i - 8] == "_" or \
                    ja[i - 7] == "_" or ja[i - 6] == "_" or \
                    ja[i - 5] == "_" or ja[i - 4] == "_" or ja[i - 3] == "_" or ja[i - 2] == "_" or ja[i - 1] == "_" or \
                    ja[i] == "_" or ja[i + 1] == "_" or ja[i + 2] == "_" or ja[i + 3] == "_" or ja[i + 4] == "_":
                print(ja[i - 12], ja[i - 11], ja[i - 9], ja[i - 8], ja[i - 7], ja[i - 6], ja[i - 5], ja[i - 4],
                      ja[i - 3], ja[i - 2], ja[i - 1],
                      ja[i], ja[i + 1], ja[i + 2], ja[i + 3], ja[i + 4], ja[i + 5], ja[i + 6], ja[i + 7], ja[i + 8],
                      ja[i + 9], ja[i + 10], ja[i + 11], ja[i + 12], ja[i + 13], ja[i + 14], ja[i + 15], ja[i + 16],
                      ja[i + 17], ja[i + 18], ja[i + 19], ja[i + 20], ja[i + 21], ja[i + 22], ja[i + 23], ja[i + 24],
                      ja[i + 25], ja[i + 26], ja[i + 27], ja[i + 28], ja[i + 29], ja[i + 30],
                      ja[i + 31], ja[i + 32], ja[i + 33], ja[i + 34], ja[i + 35], ja[i + 36], ja[i + 37], ja[i + 38],
                      ja[i + 39], ja[i + 40], ja[i + 41], ja[i + 42], ja[i + 43], ja[i + 44], ja[i + 45], ja[i + 46],
                      ja[i + 47], ja[i + 48], ja[i + 49], ja[i + 50], ja[i + 51], ja[i + 52], ja[i + 53],
                      ja[i + 54], ja[i + 55], ja[i + 56],
                      file=open("output4b.txt", "a+"), sep="")
                strix = open("output4b.txt")
                striy = strix.readlines()
                string = striy[-1]
                res = [sub.replace('_', '') for sub in string]
                x = "".join(res)
                print(x[:17], file=open("oligob.txt", "a+"))
            else:
                print(ja[i - 12], ja[i - 11], ja[i - 10], ja[i - 9], ja[i - 8], ja[i - 7], ja[i - 6], ja[i - 5],
                      ja[i - 4],
                      ja[i - 3],
                      ja[i - 2], ja[i - 1],
                      ja[i], ja[i + 1], ja[i + 2], ja[i + 3], ja[i + 4],
                      file=open("oligob.txt", "a+"), sep="")
        if i == 30:
            if ja[i - 13] == "_" or ja[i - 12] == "_" or ja[i - 11] == "_" or ja[i - 10] == "_" or ja[i - 9] == "_" or \
                    ja[i - 8] == "_" or ja[i - 7] == "_" or ja[i - 6] == "_" or \
                    ja[i - 5] == "_" or ja[i - 4] == "_" or ja[i - 3] == "_" or ja[i - 2] == "_" or ja[i - 1] == "_" or \
                    ja[i] == "_" or ja[i + 1] == "_" or ja[i + 2] == "_" or ja[i + 3] == "_":
                print(ja[i - 13], ja[i - 12], ja[i - 11], ja[i - 9], ja[i - 8], ja[i - 7], ja[i - 6], ja[i - 5],
                      ja[i - 4], ja[i - 3], ja[i - 2], ja[i - 1],
                      ja[i], ja[i + 1], ja[i + 2], ja[i + 3], ja[i + 4], ja[i + 5], ja[i + 6], ja[i + 7], ja[i + 8],
                      ja[i + 9], ja[i + 10], ja[i + 11], ja[i + 12], ja[i + 13], ja[i + 14], ja[i + 15], ja[i + 16],
                      ja[i + 17], ja[i + 18], ja[i + 19], ja[i + 20], ja[i + 21], ja[i + 22], ja[i + 23], ja[i + 24],
                      ja[i + 25], ja[i + 26], ja[i + 27], ja[i + 28], ja[i + 29], ja[i + 30],
                      ja[i + 31], ja[i + 32], ja[i + 33], ja[i + 34], ja[i + 35], ja[i + 36], ja[i + 37], ja[i + 38],
                      ja[i + 39], ja[i + 40], ja[i + 41], ja[i + 42], ja[i + 43], ja[i + 44], ja[i + 45], ja[i + 46],
                      ja[i + 47], ja[i + 48], ja[i + 49], ja[i + 50], ja[i + 51], ja[i + 52], ja[i + 53],
                      ja[i + 54], ja[i + 55],
                      file=open("output4b.txt", "a+"), sep="")
                strix = open("output4b.txt")
                striy = strix.readlines()
                string = striy[-1]
                res = [sub.replace('_', '') for sub in string]
                x = "".join(res)
                print(x[:17], file=open("oligob.txt", "a+"))
            else:
                print(ja[i - 13], ja[i - 12], ja[i - 11], ja[i - 10], ja[i - 9], ja[i - 8], ja[i - 7], ja[i - 6],
                      ja[i - 5], ja[i - 4], ja[i - 3],
                      ja[i - 2], ja[i - 1], ja[i], ja[i + 1], ja[i + 2], ja[i + 3],
                      file=open("oligob.txt", "a+"), sep="")
        if i == 31:
            if ja[i - 14] == "_" or ja[i - 13] == "_" or ja[i - 12] == "_" or ja[i - 11] == "_" or ja[i - 10] == "_" or \
                    ja[
                        i - 9] == "_" or \
                    ja[i - 8] == "_" or ja[i - 7] == "_" or ja[i - 6] == "_" or \
                    ja[i - 5] == "_" or ja[i - 4] == "_" or ja[i - 3] == "_" or ja[i - 2] == "_" or ja[i - 1] == "_" or \
                    ja[i] == "_" or ja[i + 1] == "_" or ja[i + 2] == "_":
                print(ja[i - 14], ja[i - 13], ja[i - 12], ja[i - 11], ja[i - 9], ja[i - 8], ja[i - 7], ja[i - 6],
                      ja[i - 5],
                      ja[i - 4], ja[i - 3], ja[i - 2], ja[i - 1],
                      ja[i], ja[i + 1], ja[i + 2], ja[i + 3], ja[i + 4], ja[i + 5], ja[i + 6], ja[i + 7], ja[i + 8],
                      ja[i + 9], ja[i + 10], ja[i + 11], ja[i + 12], ja[i + 13], ja[i + 14], ja[i + 15], ja[i + 16],
                      ja[i + 17], ja[i + 18], ja[i + 19], ja[i + 20], ja[i + 21], ja[i + 22], ja[i + 23], ja[i + 24],
                      ja[i + 25], ja[i + 26], ja[i + 27], ja[i + 28], ja[i + 29], ja[i + 30],
                      ja[i + 31], ja[i + 32], ja[i + 33], ja[i + 34], ja[i + 35], ja[i + 36], ja[i + 37], ja[i + 38],
                      ja[i + 39], ja[i + 40], ja[i + 41], ja[i + 42], ja[i + 43], ja[i + 44], ja[i + 45], ja[i + 46],
                      ja[i + 47], ja[i + 48], ja[i + 49], ja[i + 50], ja[i + 51], ja[i + 52], ja[i + 53],
                      ja[i + 54], file=open("output4b.txt", "a+"), sep="")
                strix = open("output4b.txt")
                striy = strix.readlines()
                string = striy[-1]
                res = [sub.replace('_', '') for sub in string]
                x = "".join(res)
                print(x[:17], file=open("oligob.txt", "a+"))
            else:
                print(ja[i - 14], ja[i - 13], ja[i - 12], ja[i - 11], ja[i - 10], ja[i - 9], ja[i - 8], ja[i - 7],
                      ja[i - 6],
                      ja[i - 5], ja[i - 4], ja[i - 3],
                      ja[i - 2], ja[i - 1], ja[i], ja[i + 1], ja[i + 2], file=open("oligob.txt", "a+"), sep="")
        if i == 32:
            if ja[i - 15] == "_" or ja[i - 14] == "_" or ja[i - 13] == "_" or ja[i - 12] == "_" or ja[i - 11] == "_" or \
                    ja[i - 10] == "_" or ja[i - 9] == "_" or ja[i - 8] == "_" or ja[i - 7] == "_" or ja[i - 6] == "_" or \
                    ja[i - 5] == "_" or ja[i - 4] == "_" or ja[i - 3] == "_" or ja[i - 2] == "_" or ja[i - 1] == "_" or \
                    ja[i] == "_" or ja[i + 1] == "_":
                print(ja[i - 15], ja[i - 14], ja[i - 13], ja[i - 12], ja[i - 11], ja[i - 9], ja[i - 8], ja[i - 7],
                      ja[i - 6], ja[i - 5],
                      ja[i - 4], ja[i - 3], ja[i - 2], ja[i - 1],
                      ja[i], ja[i + 1], ja[i + 2], ja[i + 3], ja[i + 4], ja[i + 5], ja[i + 6], ja[i + 7], ja[i + 8],
                      ja[i + 9], ja[i + 10], ja[i + 11], ja[i + 12], ja[i + 13], ja[i + 14], ja[i + 15], ja[i + 16],
                      ja[i + 17], ja[i + 18], ja[i + 19], ja[i + 20], ja[i + 21], ja[i + 22], ja[i + 23], ja[i + 24],
                      ja[i + 25], ja[i + 26], ja[i + 27], ja[i + 28], ja[i + 29], ja[i + 30],
                      ja[i + 31], ja[i + 32], ja[i + 33], ja[i + 34], ja[i + 35], ja[i + 36], ja[i + 37], ja[i + 38],
                      ja[i + 39], ja[i + 40], ja[i + 41], ja[i + 42], ja[i + 43], ja[i + 44], ja[i + 45], ja[i + 46],
                      ja[i + 47], ja[i + 48], ja[i + 49], ja[i + 50], ja[i + 51], ja[i + 52], ja[i + 53],
                      file=open("output4b.txt", "a+"), sep="")
                strix = open("output4b.txt")
                striy = strix.readlines()
                string = striy[-1]
                res = [sub.replace('_', '') for sub in string]
                x = "".join(res)
                print(x[:17], file=open("oligob.txt", "a+"))
            else:
                print(ja[i - 15], ja[i - 14], ja[i - 13], ja[i - 12], ja[i - 11], ja[i - 10], ja[i - 9], ja[i - 8],
                      ja[i - 7], ja[i - 6], ja[i - 5], ja[i - 4], ja[i - 3],
                      ja[i - 2], ja[i - 1], ja[i], ja[i + 1], file=open("oligob.txt", "a+"), sep="")
        if i == 33:
            if ja[i - 16] == "_" or ja[i - 15] == "_" or ja[i - 14] == "_" or ja[i - 13] == "_" or ja[i - 12] == "_" or \
                    ja[
                        i - 11] == "_" or \
                    ja[i - 10] == "_" or ja[i - 9] == "_" or ja[i - 8] == "_" or ja[i - 7] == "_" or ja[i - 6] == "_" or \
                    ja[i - 5] == "_" or ja[i - 4] == "_" or ja[i - 3] == "_" or ja[i - 2] == "_" or ja[i - 1] == "_" or \
                    ja[i] == "_":
                print(ja[i - 16], ja[i - 15], ja[i - 14], ja[i - 13], ja[i - 12], ja[i - 11], ja[i - 9], ja[i - 8],
                      ja[i - 7], ja[i - 6], ja[i - 5],
                      ja[i - 4], ja[i - 3], ja[i - 2], ja[i - 1],
                      ja[i], ja[i + 1], ja[i + 2], ja[i + 3], ja[i + 4], ja[i + 5], ja[i + 6], ja[i + 7], ja[i + 8],
                      ja[i + 9], ja[i + 10], ja[i + 11], ja[i + 12], ja[i + 13], ja[i + 14], ja[i + 15], ja[i + 16],
                      ja[i + 17], ja[i + 18], ja[i + 19], ja[i + 20], ja[i + 21], ja[i + 22], ja[i + 23], ja[i + 24],
                      ja[i + 25], ja[i + 26], ja[i + 27], ja[i + 28], ja[i + 29], ja[i + 30],
                      ja[i + 31], ja[i + 32], ja[i + 33], ja[i + 34], ja[i + 35], ja[i + 36], ja[i + 37], ja[i + 38],
                      ja[i + 39], ja[i + 40], ja[i + 41], ja[i + 42], ja[i + 43], ja[i + 44], ja[i + 45], ja[i + 46],
                      ja[i + 47], ja[i + 48], ja[i + 49], ja[i + 50], ja[i + 51], ja[i + 52],
                      file=open("output4b.txt", "a+"), sep="")
                strix = open("output4b.txt")
                striy = strix.readlines()
                string = striy[-1]
                res = [sub.replace('_', '') for sub in string]
                x = "".join(res)
                print(x[:17], file=open("oligob.txt", "a+"))
            else:
                print(ja[i - 16], ja[i - 15], ja[i - 14], ja[i - 13], ja[i - 12], ja[i - 11], ja[i - 10], ja[i - 9],
                      ja[i - 8],
                      ja[i - 7], ja[i - 6], ja[i - 5], ja[i - 4], ja[i - 3],
                      ja[i - 2], ja[i - 1], ja[i], file=open("oligob.txt", "a+"), sep="")
        if i == 34:
            if ja[i] == "_" or ja[i + 1] == "_" or ja[i + 2] == "_" or ja[i + 3] == "_" or ja[i + 4] == "_" or ja[i + 5] \
                    == "_" or ja[i + 6] == "_" or ja[i + 7] == "_" or ja[i + 8] == "_" or ja[i + 9] == "_" or ja[i + 10] \
                    == "_" or ja[i + 11] == "_" or ja[i + 12] == "_" or ja[i + 13] == "_" or ja[i + 14] == "_" or \
                    ja[i + 15] == "_" or ja[i + 16] == "_":
                print(ja[i], ja[i + 1], ja[i + 2], ja[i + 3], ja[i + 4], ja[i + 5], ja[i + 6], ja[i + 7], ja[i + 8],
                      ja[i + 9], ja[i + 10], ja[i + 11], ja[i + 12], ja[i + 13], ja[i + 14], ja[i + 15], ja[i + 16],
                      ja[i + 17], ja[i + 18], ja[i + 19], ja[i + 20], ja[i + 21], ja[i + 22], ja[i + 23], ja[i + 24],
                      ja[i + 25], ja[i + 26], ja[i + 27], ja[i + 28], ja[i + 29], ja[i + 30],
                      ja[i + 31], ja[i + 32], ja[i + 33], ja[i + 34], ja[i + 35], ja[i + 36], ja[i + 37], ja[i + 38],
                      ja[i + 39], ja[i + 40], ja[i + 41], ja[i + 42], ja[i + 43], ja[i + 44], ja[i + 45], ja[i + 46],
                      ja[i + 47], ja[i + 48], ja[i + 49], ja[i + 50], ja[i + 51], ja[i + 52], ja[i + 53],
                      ja[i + 54], ja[i + 55], ja[i + 56], ja[i + 57], ja[i + 58], ja[i + 59], ja[i + 60], ja[i + 61],
                      ja[i + 62], ja[i + 63], ja[i + 64], ja[i + 65], ja[i + 66], ja[i + 67], ja[i + 68],
                      file=open("output4b.txt", "a+"), sep="")
                strix = open("output4b.txt")
                striy = strix.readlines()
                string = striy[-1]
                res = [sub.replace('_', '') for sub in string]
                x = "".join(res)
                print(x[:17], file=open("oligob.txt", "a+"))
            else:
                print(ja[i], ja[i + 1], ja[i + 2], ja[i + 3], ja[i + 4], ja[i + 5], ja[i + 6], ja[i + 7], ja[i + 8],
                      ja[i + 9], ja[i + 10], ja[i + 11], ja[i + 12], ja[i + 13], ja[i + 14], ja[i + 15], ja[i + 16],
                      file=open("oligob.txt", "a+"), sep="")
        if i == ln:
            if ja[i - 16] == "_" or ja[i - 15] == "_" or ja[i - 14] == "_" or ja[i - 13] == "_" or ja[i - 12] == "_" or \
                    ja[i - 11] == "_" or ja[i - 10] == "_" or ja[i - 9] == "_" or ja[i - 8] == "_" or \
                    ja[i - 7] == "_" or ja[i - 6] == "_" or ja[i - 5] == "_" or \
                    ja[i - 4] == "_" or ja[i - 3] == "_" or ja[i - 2] == "_" or \
                    ja[i - 1] == "_" or ja[i] == "_":
                print(ja[i - 68], ja[i - 67], ja[i - 66], ja[i - 65], ja[i - 64], ja[i - 63], ja[i - 62], ja[i - 61],
                      ja[i - 60],
                      ja[i - 59], ja[i - 58], ja[i - 57], ja[i - 56], ja[i - 55], ja[i - 54], ja[i - 53], ja[i - 52],
                      ja[i - 51],
                      ja[i - 50], ja[i - 49], ja[i - 48], ja[i - 47], ja[i - 46], ja[i - 45], ja[i - 44], ja[i - 43],
                      ja[i - 42],
                      ja[i - 41], ja[i - 40], ja[i - 39], ja[i - 38], ja[i - 37], ja[i - 36], ja[i - 35], ja[i - 34],
                      ja[i - 33],
                      ja[i - 32], ja[i - 31], ja[i - 30], ja[i - 29], ja[i - 28], ja[i - 27], ja[i - 26], ja[i - 25],
                      ja[i - 24],
                      ja[i - 23], ja[i - 22], ja[i - 21], ja[i - 20], ja[i - 19], ja[i - 18], ja[i - 17], ja[i - 16],
                      ja[i - 15],
                      ja[i - 14], ja[i - 13], ja[i - 12], ja[i - 11], ja[i - 10], ja[i - 9], ja[i - 8], ja[i - 7],
                      ja[i - 6],
                      ja[i - 5], ja[i - 4], ja[i - 3], ja[i - 2], ja[i - 1], ja[i], file=open("output4b.txt", "a+"),
                      sep="")
                strix = open("output4b.txt")
                striy = strix.readlines()
                string = striy[-1]
                res = [sub.replace('_', '') for sub in string]
                x = "".join(res)
                print(x[-17:], file=open("oligob.txt", "a+"))
            else:
                print(ja[i - 16], ja[i - 15], ja[i - 14], ja[i - 13], ja[i - 11], ja[i - 10], ja[i - 9], ja[i - 8],
                      ja[i - 7], ja[i - 6], ja[i - 5], ja[i - 4], ja[i - 3], ja[i - 2], ja[i - 1], ja[i],
                      file=open("oligob.txt", "a+"), sep="")
        if i == (ln - 1):
            if ja[i - 15] == "_" or ja[i - 14] == "_" or ja[i - 13] == "_" or ja[i - 12] == "_" or \
                    ja[i - 11] == "_" or ja[i - 10] == "_" or ja[i - 9] == "_" or ja[i - 8] == "_" or \
                    ja[i - 7] == "_" or ja[i - 6] == "_" or ja[i - 5] == "_" or \
                    ja[i - 4] == "_" or ja[i - 3] == "_" or ja[i - 2] == "_" or ja[i - 1] == "_" or ja[i] == "_" or \
                    ja[i + 1] == "_":
                print(ja[i - 67], ja[i - 66], ja[i - 65], ja[i - 64], ja[i - 63], ja[i - 62], ja[i - 61],
                      ja[i - 60],
                      ja[i - 59], ja[i - 58], ja[i - 57], ja[i - 56], ja[i - 55], ja[i - 54], ja[i - 53], ja[i - 52],
                      ja[i - 51],
                      ja[i - 50], ja[i - 49], ja[i - 48], ja[i - 47], ja[i - 46], ja[i - 45], ja[i - 44], ja[i - 43],
                      ja[i - 42],
                      ja[i - 41], ja[i - 40], ja[i - 39], ja[i - 38], ja[i - 37], ja[i - 36], ja[i - 35], ja[i - 34],
                      ja[i - 33],
                      ja[i - 32], ja[i - 31], ja[i - 30], ja[i - 29], ja[i - 28], ja[i - 27], ja[i - 26], ja[i - 25],
                      ja[i - 24],
                      ja[i - 23], ja[i - 22], ja[i - 21], ja[i - 20], ja[i - 19], ja[i - 18], ja[i - 17], ja[i - 16],
                      ja[i - 15],
                      ja[i - 14], ja[i - 13], ja[i - 12], ja[i - 11], ja[i - 10], ja[i - 9], ja[i - 8], ja[i - 7],
                      ja[i - 6],
                      ja[i - 5], ja[i - 4], ja[i - 3], ja[i - 2], ja[i - 1], ja[i], ja[i + 1],
                      file=open("output4b.txt", "a+"), sep="")
                strix = open("output4b.txt")
                striy = strix.readlines()
                string = striy[-1]
                res = [sub.replace('_', '') for sub in string]
                x = "".join(res)
                print(x[-17:], file=open("oligob.txt", "a+"))
            else:
                print(ja[i - 15], ja[i - 14], ja[i - 13], ja[i - 11], ja[i - 10], ja[i - 9], ja[i - 8],
                      ja[i - 7], ja[i - 6], ja[i - 5], ja[i - 4], ja[i - 3], ja[i - 2], ja[i - 1], ja[i], ja[i + 1],
                      file=open("oligob.txt", "a+"), sep="")
        if i == (ln - 2):
            if ja[i - 14] == "_" or ja[i - 13] == "_" or ja[i - 12] == "_" or \
                    ja[i - 11] == "_" or ja[i - 10] == "_" or ja[i - 9] == "_" or ja[i - 8] == "_" or \
                    ja[i - 7] == "_" or ja[i - 6] == "_" or ja[i - 5] == "_" or \
                    ja[i - 4] == "_" or ja[i - 3] == "_" or ja[i - 2] == "_" or ja[i - 1] == "_" or ja[i] == "_" or \
                    ja[i + 1] == "_" or ja[i + 2] == "_":
                print(ja[i - 66], ja[i - 65], ja[i - 64], ja[i - 63], ja[i - 62], ja[i - 61],
                      ja[i - 60],
                      ja[i - 59], ja[i - 58], ja[i - 57], ja[i - 56], ja[i - 55], ja[i - 54], ja[i - 53], ja[i - 52],
                      ja[i - 51],
                      ja[i - 50], ja[i - 49], ja[i - 48], ja[i - 47], ja[i - 46], ja[i - 45], ja[i - 44], ja[i - 43],
                      ja[i - 42],
                      ja[i - 41], ja[i - 40], ja[i - 39], ja[i - 38], ja[i - 37], ja[i - 36], ja[i - 35], ja[i - 34],
                      ja[i - 33],
                      ja[i - 32], ja[i - 31], ja[i - 30], ja[i - 29], ja[i - 28], ja[i - 27], ja[i - 26], ja[i - 25],
                      ja[i - 24],
                      ja[i - 23], ja[i - 22], ja[i - 21], ja[i - 20], ja[i - 19], ja[i - 18], ja[i - 17], ja[i - 16],
                      ja[i - 15],
                      ja[i - 14], ja[i - 13], ja[i - 12], ja[i - 11], ja[i - 10], ja[i - 9], ja[i - 8], ja[i - 7],
                      ja[i - 6],
                      ja[i - 5], ja[i - 4], ja[i - 3], ja[i - 2], ja[i - 1], ja[i], ja[i + 1], ja[i + 2],
                      file=open("output4b.txt", "a+"), sep="")
                strix = open("output4b.txt")
                striy = strix.readlines()
                string = striy[-1]
                res = [sub.replace('_', '') for sub in string]
                x = "".join(res)
                print(x[-17:], file=open("oligob.txt", "a+"))
            else:
                print(ja[i - 14], ja[i - 13], ja[i - 11], ja[i - 10], ja[i - 9], ja[i - 8],
                      ja[i - 7], ja[i - 6], ja[i - 5], ja[i - 4], ja[i - 3], ja[i - 2], ja[i - 1], ja[i], ja[i + 1],
                      ja[i + 2], file=open("oligob.txt", "a+"), sep="")
        if i == (ln - 3):
            if ja[i - 13] == "_" or ja[i - 12] == "_" or \
                    ja[i - 11] == "_" or ja[i - 10] == "_" or ja[i - 9] == "_" or ja[i - 8] == "_" or \
                    ja[i - 7] == "_" or ja[i - 6] == "_" or ja[i - 5] == "_" or \
                    ja[i - 4] == "_" or ja[i - 3] == "_" or ja[i - 2] == "_" or ja[i - 1] == "_" or ja[i] == "_" or \
                    ja[i + 1] == "_" or ja[i + 2] == "_" or ja[i + 3] == "_":
                print(ja[i - 65], ja[i - 64], ja[i - 63], ja[i - 62], ja[i - 61],
                      ja[i - 60],
                      ja[i - 59], ja[i - 58], ja[i - 57], ja[i - 56], ja[i - 55], ja[i - 54], ja[i - 53], ja[i - 52],
                      ja[i - 51],
                      ja[i - 50], ja[i - 49], ja[i - 48], ja[i - 47], ja[i - 46], ja[i - 45], ja[i - 44], ja[i - 43],
                      ja[i - 42],
                      ja[i - 41], ja[i - 40], ja[i - 39], ja[i - 38], ja[i - 37], ja[i - 36], ja[i - 35], ja[i - 34],
                      ja[i - 33],
                      ja[i - 32], ja[i - 31], ja[i - 30], ja[i - 29], ja[i - 28], ja[i - 27], ja[i - 26], ja[i - 25],
                      ja[i - 24],
                      ja[i - 23], ja[i - 22], ja[i - 21], ja[i - 20], ja[i - 19], ja[i - 18], ja[i - 17], ja[i - 16],
                      ja[i - 15],
                      ja[i - 14], ja[i - 13], ja[i - 12], ja[i - 11], ja[i - 10], ja[i - 9], ja[i - 8], ja[i - 7],
                      ja[i - 6],
                      ja[i - 5], ja[i - 4], ja[i - 3], ja[i - 2], ja[i - 1], ja[i], ja[i + 1], ja[i + 2], ja[i + 3],
                      file=open("output4b.txt", "a+"), sep="")
                strix = open("output4b.txt")
                striy = strix.readlines()
                string = striy[-1]
                res = [sub.replace('_', '') for sub in string]
                x = "".join(res)
                print(x[-17:], file=open("oligob.txt", "a+"))
            else:
                print(ja[i - 13], ja[i - 11], ja[i - 10], ja[i - 9], ja[i - 8],
                      ja[i - 7], ja[i - 6], ja[i - 5], ja[i - 4], ja[i - 3], ja[i - 2], ja[i - 1], ja[i], ja[i + 1],
                      ja[i + 2], ja[i + 3], file=open("oligob.txt", "a+"), sep="")
        if i == (ln - 4):
            if ja[i - 12] == "_" or ja[i - 11] == "_" or ja[i - 10] == "_" or ja[i - 9] == "_" or ja[i - 8] == "_" or \
                    ja[i - 7] == "_" or ja[i - 6] == "_" or ja[i - 5] == "_" or \
                    ja[i - 4] == "_" or ja[i - 3] == "_" or ja[i - 2] == "_" or ja[i - 1] == "_" or ja[i] == "_" or \
                    ja[i + 1] == "_" or ja[i + 2] == "_" or ja[i + 3] == "_" or ja[i + 4] == "_":
                print(ja[i - 64], ja[i - 63], ja[i - 62], ja[i - 61],
                      ja[i - 60],
                      ja[i - 59], ja[i - 58], ja[i - 57], ja[i - 56], ja[i - 55], ja[i - 54], ja[i - 53], ja[i - 52],
                      ja[i - 51],
                      ja[i - 50], ja[i - 49], ja[i - 48], ja[i - 47], ja[i - 46], ja[i - 45], ja[i - 44], ja[i - 43],
                      ja[i - 42],
                      ja[i - 41], ja[i - 40], ja[i - 39], ja[i - 38], ja[i - 37], ja[i - 36], ja[i - 35], ja[i - 34],
                      ja[i - 33],
                      ja[i - 32], ja[i - 31], ja[i - 30], ja[i - 29], ja[i - 28], ja[i - 27], ja[i - 26], ja[i - 25],
                      ja[i - 24],
                      ja[i - 23], ja[i - 22], ja[i - 21], ja[i - 20], ja[i - 19], ja[i - 18], ja[i - 17], ja[i - 16],
                      ja[i - 15],
                      ja[i - 14], ja[i - 13], ja[i - 12], ja[i - 11], ja[i - 10], ja[i - 9], ja[i - 8], ja[i - 7],
                      ja[i - 6],
                      ja[i - 5], ja[i - 4], ja[i - 3], ja[i - 2], ja[i - 1], ja[i], ja[i + 1], ja[i + 2], ja[i + 3],
                      ja[i + 4], file=open("output4b.txt", "a+"), sep="")
                strix = open("output4b.txt")
                striy = strix.readlines()
                string = striy[-1]
                res = [sub.replace('_', '') for sub in string]
                x = "".join(res)
                print(x[-17:], file=open("oligob.txt", "a+"))
            else:
                print(ja[i - 12], ja[i - 11], ja[i - 10], ja[i - 9], ja[i - 8],
                      ja[i - 7], ja[i - 6], ja[i - 5], ja[i - 4], ja[i - 3], ja[i - 2], ja[i - 1], ja[i], ja[i + 1],
                      ja[i + 2], ja[i + 3], ja[i + 4], file=open("oligob.txt", "a+"), sep="")
        if i == (ln - 5):
            if ja[i - 11] == "_" or ja[i - 10] == "_" or ja[i - 9] == "_" or ja[i - 8] == "_" or \
                    ja[i - 7] == "_" or ja[i - 6] == "_" or ja[i - 5] == "_" or \
                    ja[i - 4] == "_" or ja[i - 3] == "_" or ja[i - 2] == "_" or ja[i - 1] == "_" or ja[i] == "_" or \
                    ja[i + 1] == "_" or ja[i + 2] == "_" or ja[i + 3] == "_" or ja[i + 4] == "_" or ja[i + 5] == "_":
                print(ja[i - 63], ja[i - 62], ja[i - 61],
                      ja[i - 60],
                      ja[i - 59], ja[i - 58], ja[i - 57], ja[i - 56], ja[i - 55], ja[i - 54], ja[i - 53], ja[i - 52],
                      ja[i - 51],
                      ja[i - 50], ja[i - 49], ja[i - 48], ja[i - 47], ja[i - 46], ja[i - 45], ja[i - 44], ja[i - 43],
                      ja[i - 42],
                      ja[i - 41], ja[i - 40], ja[i - 39], ja[i - 38], ja[i - 37], ja[i - 36], ja[i - 35], ja[i - 34],
                      ja[i - 33],
                      ja[i - 32], ja[i - 31], ja[i - 30], ja[i - 29], ja[i - 28], ja[i - 27], ja[i - 26], ja[i - 25],
                      ja[i - 24],
                      ja[i - 23], ja[i - 22], ja[i - 21], ja[i - 20], ja[i - 19], ja[i - 18], ja[i - 17], ja[i - 16],
                      ja[i - 15],
                      ja[i - 14], ja[i - 13], ja[i - 12], ja[i - 11], ja[i - 10], ja[i - 9], ja[i - 8], ja[i - 7],
                      ja[i - 6],
                      ja[i - 5], ja[i - 4], ja[i - 3], ja[i - 2], ja[i - 1], ja[i], ja[i + 1], ja[i + 2], ja[i + 3],
                      ja[i + 4], ja[i + 5], file=open("output4b.txt", "a+"), sep="")
                strix = open("output4b.txt")
                striy = strix.readlines()
                string = striy[-1]
                res = [sub.replace('_', '') for sub in string]
                x = "".join(res)
                print(x[-17:], file=open("oligob.txt", "a+"))
            else:
                print(ja[i - 11], ja[i - 10], ja[i - 9], ja[i - 8],
                      ja[i - 7], ja[i - 6], ja[i - 5], ja[i - 4], ja[i - 3], ja[i - 2], ja[i - 1], ja[i], ja[i + 1],
                      ja[i + 2], ja[i + 3], ja[i + 4], ja[i + 5], file=open("oligob.txt", "a+"), sep="")
        if i == (ln - 6):
            if ja[i - 10] == "_" or ja[i - 9] == "_" or ja[i - 8] == "_" or \
                    ja[i - 7] == "_" or ja[i - 6] == "_" or ja[i - 5] == "_" or \
                    ja[i - 4] == "_" or ja[i - 3] == "_" or ja[i - 2] == "_" or ja[i - 1] == "_" or ja[i] == "_" or \
                    ja[i + 1] == "_" or ja[i + 2] == "_" or ja[i + 3] == "_" or ja[i + 4] == "_" or ja[i + 5] == "_" or \
                    ja[i + 6] == "_":
                print(ja[i - 62], ja[i - 61],
                      ja[i - 60],
                      ja[i - 59], ja[i - 58], ja[i - 57], ja[i - 56], ja[i - 55], ja[i - 54], ja[i - 53], ja[i - 52],
                      ja[i - 51],
                      ja[i - 50], ja[i - 49], ja[i - 48], ja[i - 47], ja[i - 46], ja[i - 45], ja[i - 44], ja[i - 43],
                      ja[i - 42],
                      ja[i - 41], ja[i - 40], ja[i - 39], ja[i - 38], ja[i - 37], ja[i - 36], ja[i - 35], ja[i - 34],
                      ja[i - 33],
                      ja[i - 32], ja[i - 31], ja[i - 30], ja[i - 29], ja[i - 28], ja[i - 27], ja[i - 26], ja[i - 25],
                      ja[i - 24],
                      ja[i - 23], ja[i - 22], ja[i - 21], ja[i - 20], ja[i - 19], ja[i - 18], ja[i - 17], ja[i - 16],
                      ja[i - 15],
                      ja[i - 14], ja[i - 13], ja[i - 12], ja[i - 11], ja[i - 10], ja[i - 9], ja[i - 8], ja[i - 7],
                      ja[i - 6],
                      ja[i - 5], ja[i - 4], ja[i - 3], ja[i - 2], ja[i - 1], ja[i], ja[i + 1], ja[i + 2], ja[i + 3],
                      ja[i + 4], ja[i + 5], ja[i + 6], file=open("output4b.txt", "a+"), sep="")
                strix = open("output4b.txt")
                striy = strix.readlines()
                string = striy[-1]
                res = [sub.replace('_', '') for sub in string]
                x = "".join(res)
                print(x[-17:], file=open("oligob.txt", "a+"))
            else:
                print(ja[i - 10], ja[i - 9], ja[i - 8],
                      ja[i - 7], ja[i - 6], ja[i - 5], ja[i - 4], ja[i - 3], ja[i - 2], ja[i - 1], ja[i], ja[i + 1],
                      ja[i + 2], ja[i + 3], ja[i + 4], ja[i + 5], ja[i + 6], file=open("oligob.txt", "a+"), sep="")
        if i == (ln - 7):
            if ja[i - 9] == "_" or ja[i - 8] == "_" or \
                    ja[i - 7] == "_" or ja[i - 6] == "_" or ja[i - 5] == "_" or \
                    ja[i - 4] == "_" or ja[i - 3] == "_" or ja[i - 2] == "_" or ja[i - 1] == "_" or ja[i] == "_" or \
                    ja[i + 1] == "_" or ja[i + 2] == "_" or ja[i + 3] == "_" or ja[i + 4] == "_" or ja[i + 5] == "_" or \
                    ja[i + 6] == "_" or ja[i + 7] == "_":
                print(ja[i - 61],
                      ja[i - 60],
                      ja[i - 59], ja[i - 58], ja[i - 57], ja[i - 56], ja[i - 55], ja[i - 54], ja[i - 53], ja[i - 52],
                      ja[i - 51],
                      ja[i - 50], ja[i - 49], ja[i - 48], ja[i - 47], ja[i - 46], ja[i - 45], ja[i - 44], ja[i - 43],
                      ja[i - 42],
                      ja[i - 41], ja[i - 40], ja[i - 39], ja[i - 38], ja[i - 37], ja[i - 36], ja[i - 35], ja[i - 34],
                      ja[i - 33],
                      ja[i - 32], ja[i - 31], ja[i - 30], ja[i - 29], ja[i - 28], ja[i - 27], ja[i - 26], ja[i - 25],
                      ja[i - 24],
                      ja[i - 23], ja[i - 22], ja[i - 21], ja[i - 20], ja[i - 19], ja[i - 18], ja[i - 17], ja[i - 16],
                      ja[i - 15],
                      ja[i - 14], ja[i - 13], ja[i - 12], ja[i - 11], ja[i - 10], ja[i - 9], ja[i - 8], ja[i - 7],
                      ja[i - 6],
                      ja[i - 5], ja[i - 4], ja[i - 3], ja[i - 2], ja[i - 1], ja[i], ja[i + 1], ja[i + 2], ja[i + 3],
                      ja[i + 4], ja[i + 5], ja[i + 6], ja[i + 7], file=open("output4b.txt", "a+"), sep="")
                strix = open("output4b.txt")
                striy = strix.readlines()
                string = striy[-1]
                res = [sub.replace('_', '') for sub in string]
                x = "".join(res)
                print(x[-17:], file=open("oligob.txt", "a+"))
            else:
                print(ja[i - 9], ja[i - 8],
                      ja[i - 7], ja[i - 6], ja[i - 5], ja[i - 4], ja[i - 3], ja[i - 2], ja[i - 1], ja[i], ja[i + 1],
                      ja[i + 2], ja[i + 3], ja[i + 4], ja[i + 5], ja[i + 6], ja[i + 7], file=open("oligob.txt", "a+"),
                      sep="")
        if i == (ln - 8):
            if ja[i - 8] == "_" or \
                    ja[i - 7] == "_" or ja[i - 6] == "_" or ja[i - 5] == "_" or \
                    ja[i - 4] == "_" or ja[i - 3] == "_" or ja[i - 2] == "_" or ja[i - 1] == "_" or ja[i] == "_" or \
                    ja[i + 1] == "_" or ja[i + 2] == "_" or ja[i + 3] == "_" or ja[i + 4] == "_" or ja[i + 5] == "_" or \
                    ja[i + 6] == "_" or ja[i + 7] == "_" or ja[i + 8] == "_":
                print(ja[i - 60],
                      ja[i - 59], ja[i - 58], ja[i - 57], ja[i - 56], ja[i - 55], ja[i - 54], ja[i - 53], ja[i - 52],
                      ja[i - 51],
                      ja[i - 50], ja[i - 49], ja[i - 48], ja[i - 47], ja[i - 46], ja[i - 45], ja[i - 44], ja[i - 43],
                      ja[i - 42],
                      ja[i - 41], ja[i - 40], ja[i - 39], ja[i - 38], ja[i - 37], ja[i - 36], ja[i - 35], ja[i - 34],
                      ja[i - 33],
                      ja[i - 32], ja[i - 31], ja[i - 30], ja[i - 29], ja[i - 28], ja[i - 27], ja[i - 26], ja[i - 25],
                      ja[i - 24],
                      ja[i - 23], ja[i - 22], ja[i - 21], ja[i - 20], ja[i - 19], ja[i - 18], ja[i - 17], ja[i - 16],
                      ja[i - 15],
                      ja[i - 14], ja[i - 13], ja[i - 12], ja[i - 11], ja[i - 10], ja[i - 9], ja[i - 8], ja[i - 7],
                      ja[i - 6],
                      ja[i - 5], ja[i - 4], ja[i - 3], ja[i - 2], ja[i - 1], ja[i], ja[i + 1], ja[i + 2], ja[i + 3],
                      ja[i + 4], ja[i + 5], ja[i + 6], ja[i + 7], ja[i + 8], file=open("output4b.txt", "a+"), sep="")
                strix = open("output4b.txt")
                striy = strix.readlines()
                string = striy[-1]
                res = [sub.replace('_', '') for sub in string]
                x = "".join(res)
                print(x[-17:], file=open("oligob.txt", "a+"))
            else:
                print(ja[i - 8],
                      ja[i - 7], ja[i - 6], ja[i - 5], ja[i - 4], ja[i - 3], ja[i - 2], ja[i - 1], ja[i], ja[i + 1],
                      ja[i + 2], ja[i + 3], ja[i + 4], ja[i + 5], ja[i + 6], ja[i + 7], ja[i + 8],
                      file=open("oligob.txt", "a+"), sep="")
        if i == (ln - 9):
            if ja[i - 7] == "_" or ja[i - 6] == "_" or ja[i - 5] == "_" or \
                    ja[i - 4] == "_" or ja[i - 3] == "_" or ja[i - 2] == "_" or ja[i - 1] == "_" or ja[i] == "_" or \
                    ja[i + 1] == "_" or ja[i + 2] == "_" or ja[i + 3] == "_" or ja[i + 4] == "_" or ja[i + 5] == "_" or \
                    ja[i + 6] == "_" or ja[i + 7] == "_" or ja[i + 8] == "_" or ja[i + 9] == "_":
                print(ja[i - 59], ja[i - 58], ja[i - 57], ja[i - 56], ja[i - 55], ja[i - 54], ja[i - 53], ja[i - 52],
                      ja[i - 51],
                      ja[i - 50], ja[i - 49], ja[i - 48], ja[i - 47], ja[i - 46], ja[i - 45], ja[i - 44], ja[i - 43],
                      ja[i - 42],
                      ja[i - 41], ja[i - 40], ja[i - 39], ja[i - 38], ja[i - 37], ja[i - 36], ja[i - 35], ja[i - 34],
                      ja[i - 33],
                      ja[i - 32], ja[i - 31], ja[i - 30], ja[i - 29], ja[i - 28], ja[i - 27], ja[i - 26], ja[i - 25],
                      ja[i - 24],
                      ja[i - 23], ja[i - 22], ja[i - 21], ja[i - 20], ja[i - 19], ja[i - 18], ja[i - 17], ja[i - 16],
                      ja[i - 15],
                      ja[i - 14], ja[i - 13], ja[i - 12], ja[i - 11], ja[i - 10], ja[i - 9], ja[i - 8], ja[i - 7],
                      ja[i - 6],
                      ja[i - 5], ja[i - 4], ja[i - 3], ja[i - 2], ja[i - 1], ja[i], ja[i + 1], ja[i + 2], ja[i + 3],
                      ja[i + 4], ja[i + 5], ja[i + 6], ja[i + 7], ja[i + 8], ja[i + 9], file=open("output4b.txt", "a+"),
                      sep="")
                strix = open("output4b.txt")
                striy = strix.readlines()
                string = striy[-1]
                res = [sub.replace('_', '') for sub in string]
                x = "".join(res)
                print(x[-17:], file=open("oligob.txt", "a+"))
            else:
                print(ja[i - 7], ja[i - 6], ja[i - 5], ja[i - 4], ja[i - 3], ja[i - 2], ja[i - 1], ja[i], ja[i + 1],
                      ja[i + 2], ja[i + 3], ja[i + 4], ja[i + 5], ja[i + 6], ja[i + 7], ja[i + 8], ja[i + 9],
                      file=open("oligob.txt", "a+"), sep="")
        if i == (ln - 10):
            if ja[i - 6] == "_" or ja[i - 5] == "_" or \
                    ja[i - 4] == "_" or ja[i - 3] == "_" or ja[i - 2] == "_" or ja[i - 1] == "_" or ja[i] == "_" or \
                    ja[i + 1] == "_" or ja[i + 2] == "_" or ja[i + 3] == "_" or ja[i + 4] == "_" or ja[i + 5] == "_" or \
                    ja[i + 6] == "_" or ja[i + 7] == "_" or ja[i + 8] == "_" or ja[i + 9] == "_" or ja[i + 10] == "_":
                print(ja[i - 58], ja[i - 57], ja[i - 56], ja[i - 55], ja[i - 54], ja[i - 53], ja[i - 52],
                      ja[i - 51],
                      ja[i - 50], ja[i - 49], ja[i - 48], ja[i - 47], ja[i - 46], ja[i - 45], ja[i - 44], ja[i - 43],
                      ja[i - 42],
                      ja[i - 41], ja[i - 40], ja[i - 39], ja[i - 38], ja[i - 37], ja[i - 36], ja[i - 35], ja[i - 34],
                      ja[i - 33],
                      ja[i - 32], ja[i - 31], ja[i - 30], ja[i - 29], ja[i - 28], ja[i - 27], ja[i - 26], ja[i - 25],
                      ja[i - 24],
                      ja[i - 23], ja[i - 22], ja[i - 21], ja[i - 20], ja[i - 19], ja[i - 18], ja[i - 17], ja[i - 16],
                      ja[i - 15],
                      ja[i - 14], ja[i - 13], ja[i - 12], ja[i - 11], ja[i - 10], ja[i - 9], ja[i - 8], ja[i - 7],
                      ja[i - 6],
                      ja[i - 5], ja[i - 4], ja[i - 3], ja[i - 2], ja[i - 1], ja[i], ja[i + 1], ja[i + 2], ja[i + 3],
                      ja[i + 4], ja[i + 5], ja[i + 6], ja[i + 7], ja[i + 8], ja[i + 9], ja[i + 10],
                      file=open("output4b.txt", "a+"), sep="")
                strix = open("output4b.txt")
                striy = strix.readlines()
                string = striy[-1]
                res = [sub.replace('_', '') for sub in string]
                x = "".join(res)
                print(x[-17:], file=open("oligob.txt", "a+"))
            else:
                print(ja[i - 6], ja[i - 5], ja[i - 4], ja[i - 3], ja[i - 2], ja[i - 1], ja[i], ja[i + 1],
                      ja[i + 2], ja[i + 3], ja[i + 4], ja[i + 5], ja[i + 6], ja[i + 7], ja[i + 8], ja[i + 9],
                      ja[i + 10],
                      file=open("oligob.txt", "a+"), sep="")
        if i == (ln - 11):
            if ja[i - 5] == "_" or \
                    ja[i - 4] == "_" or ja[i - 3] == "_" or ja[i - 2] == "_" or ja[i - 1] == "_" or ja[i] == "_" or \
                    ja[i + 1] == "_" or ja[i + 2] == "_" or ja[i + 3] == "_" or ja[i + 4] == "_" or ja[i + 5] == "_" or \
                    ja[i + 6] == "_" or ja[i + 7] == "_" or ja[i + 8] == "_" or ja[i + 9] == "_" or ja[i + 10] == "_" or \
                    ja[i + 11] == "_":
                print(ja[i - 57], ja[i - 56], ja[i - 55], ja[i - 54], ja[i - 53], ja[i - 52],
                      ja[i - 51],
                      ja[i - 50], ja[i - 49], ja[i - 48], ja[i - 47], ja[i - 46], ja[i - 45], ja[i - 44], ja[i - 43],
                      ja[i - 42],
                      ja[i - 41], ja[i - 40], ja[i - 39], ja[i - 38], ja[i - 37], ja[i - 36], ja[i - 35], ja[i - 34],
                      ja[i - 33],
                      ja[i - 32], ja[i - 31], ja[i - 30], ja[i - 29], ja[i - 28], ja[i - 27], ja[i - 26], ja[i - 25],
                      ja[i - 24],
                      ja[i - 23], ja[i - 22], ja[i - 21], ja[i - 20], ja[i - 19], ja[i - 18], ja[i - 17], ja[i - 16],
                      ja[i - 15],
                      ja[i - 14], ja[i - 13], ja[i - 12], ja[i - 11], ja[i - 10], ja[i - 9], ja[i - 8], ja[i - 7],
                      ja[i - 6],
                      ja[i - 5], ja[i - 4], ja[i - 3], ja[i - 2], ja[i - 1], ja[i], ja[i + 1], ja[i + 2], ja[i + 3],
                      ja[i + 4], ja[i + 5], ja[i + 6], ja[i + 7], ja[i + 8], ja[i + 9], ja[i + 10], ja[i + 11],
                      file=open("output4b.txt", "a+"), sep="")
                strix = open("output4b.txt")
                striy = strix.readlines()
                string = striy[-1]
                res = [sub.replace('_', '') for sub in string]
                x = "".join(res)
                print(x[-17:], file=open("oligob.txt", "a+"))
            else:
                print(ja[i - 5], ja[i - 4], ja[i - 3], ja[i - 2], ja[i - 1], ja[i], ja[i + 1],
                      ja[i + 2], ja[i + 3], ja[i + 4], ja[i + 5], ja[i + 6], ja[i + 7], ja[i + 8], ja[i + 9],
                      ja[i + 10], ja[i + 11], file=open("oligob.txt", "a+"), sep="")
        if i == (ln - 12):
            if ja[i - 4] == "_" or ja[i - 3] == "_" or ja[i - 2] == "_" or ja[i - 1] == "_" or ja[i] == "_" or \
                    ja[i + 1] == "_" or ja[i + 2] == "_" or ja[i + 3] == "_" or ja[i + 4] == "_" or ja[i + 5] == "_" or \
                    ja[i + 6] == "_" or ja[i + 7] == "_" or ja[i + 8] == "_" or ja[i + 9] == "_" or ja[i + 10] == "_" or \
                    ja[i + 11] == "_" or ja[i + 12] == "_":
                print(ja[i - 56], ja[i - 55], ja[i - 54], ja[i - 53], ja[i - 52],
                      ja[i - 51],
                      ja[i - 50], ja[i - 49], ja[i - 48], ja[i - 47], ja[i - 46], ja[i - 45], ja[i - 44], ja[i - 43],
                      ja[i - 42],
                      ja[i - 41], ja[i - 40], ja[i - 39], ja[i - 38], ja[i - 37], ja[i - 36], ja[i - 35], ja[i - 34],
                      ja[i - 33],
                      ja[i - 32], ja[i - 31], ja[i - 30], ja[i - 29], ja[i - 28], ja[i - 27], ja[i - 26], ja[i - 25],
                      ja[i - 24],
                      ja[i - 23], ja[i - 22], ja[i - 21], ja[i - 20], ja[i - 19], ja[i - 18], ja[i - 17], ja[i - 16],
                      ja[i - 15],
                      ja[i - 14], ja[i - 13], ja[i - 12], ja[i - 11], ja[i - 10], ja[i - 9], ja[i - 8], ja[i - 7],
                      ja[i - 6],
                      ja[i - 5], ja[i - 4], ja[i - 3], ja[i - 2], ja[i - 1], ja[i], ja[i + 1], ja[i + 2], ja[i + 3],
                      ja[i + 4], ja[i + 5], ja[i + 6], ja[i + 7], ja[i + 8], ja[i + 9], ja[i + 10], ja[i + 11],
                      ja[i + 12],
                      file=open("output4b.txt", "a+"), sep="")
                strix = open("output4b.txt")
                striy = strix.readlines()
                string = striy[-1]
                res = [sub.replace('_', '') for sub in string]
                x = "".join(res)
                print(x[-17:], file=open("oligob.txt", "a+"))
            else:
                print(ja[i - 4], ja[i - 3], ja[i - 2], ja[i - 1], ja[i], ja[i + 1],
                      ja[i + 2], ja[i + 3], ja[i + 4], ja[i + 5], ja[i + 6], ja[i + 7], ja[i + 8], ja[i + 9],
                      ja[i + 10], ja[i + 11], ja[i + 12], file=open("oligob.txt", "a+"), sep="")
        if i == (ln - 13):
            if ja[i - 3] == "_" or ja[i - 2] == "_" or ja[i - 1] == "_" or ja[i] == "_" or \
                    ja[i + 1] == "_" or ja[i + 2] == "_" or ja[i + 3] == "_" or ja[i + 4] == "_" or ja[i + 5] == "_" or \
                    ja[i + 6] == "_" or ja[i + 7] == "_" or ja[i + 8] == "_" or ja[i + 9] == "_" or ja[i + 10] == "_" or \
                    ja[i + 11] == "_" or ja[i + 12] == "_" or ja[i + 13] == "_":
                print(ja[i - 55], ja[i - 54], ja[i - 53], ja[i - 52],
                      ja[i - 51],
                      ja[i - 50], ja[i - 49], ja[i - 48], ja[i - 47], ja[i - 46], ja[i - 45], ja[i - 44], ja[i - 43],
                      ja[i - 42],
                      ja[i - 41], ja[i - 40], ja[i - 39], ja[i - 38], ja[i - 37], ja[i - 36], ja[i - 35], ja[i - 34],
                      ja[i - 33],
                      ja[i - 32], ja[i - 31], ja[i - 30], ja[i - 29], ja[i - 28], ja[i - 27], ja[i - 26], ja[i - 25],
                      ja[i - 24],
                      ja[i - 23], ja[i - 22], ja[i - 21], ja[i - 20], ja[i - 19], ja[i - 18], ja[i - 17], ja[i - 16],
                      ja[i - 15],
                      ja[i - 14], ja[i - 13], ja[i - 12], ja[i - 11], ja[i - 10], ja[i - 9], ja[i - 8], ja[i - 7],
                      ja[i - 6],
                      ja[i - 5], ja[i - 4], ja[i - 3], ja[i - 2], ja[i - 1], ja[i], ja[i + 1], ja[i + 2], ja[i + 3],
                      ja[i + 4], ja[i + 5], ja[i + 6], ja[i + 7], ja[i + 8], ja[i + 9], ja[i + 10], ja[i + 11],
                      ja[i + 12],
                      ja[i + 13], file=open("output4b.txt", "a+"), sep="")
                strix = open("output4b.txt")
                striy = strix.readlines()
                string = striy[-1]
                res = [sub.replace('_', '') for sub in string]
                x = "".join(res)
                print(x[-17:], file=open("oligob.txt", "a+"))
            else:
                print(ja[i - 3], ja[i - 2], ja[i - 1], ja[i], ja[i + 1],
                      ja[i + 2], ja[i + 3], ja[i + 4], ja[i + 5], ja[i + 6], ja[i + 7], ja[i + 8], ja[i + 9],
                      ja[i + 10], ja[i + 11], ja[i + 12], ja[i + 13], file=open("oligob.txt", "a+"), sep="")
        if i == (ln - 14):
            if ja[i - 2] == "_" or ja[i - 1] == "_" or ja[i] == "_" or \
                    ja[i + 1] == "_" or ja[i + 2] == "_" or ja[i + 3] == "_" or ja[i + 4] == "_" or ja[i + 5] == "_" or \
                    ja[i + 6] == "_" or ja[i + 7] == "_" or ja[i + 8] == "_" or ja[i + 9] == "_" or ja[i + 10] == "_" or \
                    ja[i + 11] == "_" or ja[i + 12] == "_" or ja[i + 13] == "_" or ja[i + 14] == "_":
                print(ja[i - 54], ja[i - 53], ja[i - 52],
                      ja[i - 51],
                      ja[i - 50], ja[i - 49], ja[i - 48], ja[i - 47], ja[i - 46], ja[i - 45], ja[i - 44], ja[i - 43],
                      ja[i - 42],
                      ja[i - 41], ja[i - 40], ja[i - 39], ja[i - 38], ja[i - 37], ja[i - 36], ja[i - 35], ja[i - 34],
                      ja[i - 33],
                      ja[i - 32], ja[i - 31], ja[i - 30], ja[i - 29], ja[i - 28], ja[i - 27], ja[i - 26], ja[i - 25],
                      ja[i - 24],
                      ja[i - 23], ja[i - 22], ja[i - 21], ja[i - 20], ja[i - 19], ja[i - 18], ja[i - 17], ja[i - 16],
                      ja[i - 15],
                      ja[i - 14], ja[i - 13], ja[i - 12], ja[i - 11], ja[i - 10], ja[i - 9], ja[i - 8], ja[i - 7],
                      ja[i - 6],
                      ja[i - 5], ja[i - 4], ja[i - 3], ja[i - 2], ja[i - 1], ja[i], ja[i + 1], ja[i + 2], ja[i + 3],
                      ja[i + 4], ja[i + 5], ja[i + 6], ja[i + 7], ja[i + 8], ja[i + 9], ja[i + 10], ja[i + 11],
                      ja[i + 12],
                      ja[i + 13], ja[i + 14], file=open("output4b.txt", "a+"), sep="")
                strix = open("output4b.txt")
                striy = strix.readlines()
                string = striy[-1]
                res = [sub.replace('_', '') for sub in string]
                x = "".join(res)
                print(x[-17:], file=open("oligob.txt", "a+"))
            else:
                print(ja[i - 2], ja[i - 1], ja[i], ja[i + 1],
                      ja[i + 2], ja[i + 3], ja[i + 4], ja[i + 5], ja[i + 6], ja[i + 7], ja[i + 8], ja[i + 9],
                      ja[i + 10], ja[i + 11], ja[i + 12], ja[i + 13], ja[i + 14], file=open("oligob.txt", "a+"), sep="")
        if i == (ln - 15):
            if ja[i - 1] == "_" or ja[i] == "_" or \
                    ja[i + 1] == "_" or ja[i + 2] == "_" or ja[i + 3] == "_" or ja[i + 4] == "_" or ja[i + 5] == "_" or \
                    ja[i + 6] == "_" or ja[i + 7] == "_" or ja[i + 8] == "_" or ja[i + 9] == "_" or ja[i + 10] == "_" or \
                    ja[i + 11] == "_" or ja[i + 12] == "_" or ja[i + 13] == "_" or ja[i + 14] == "_" or \
                    ja[i + 15] == "_":
                print(ja[i - 53], ja[i - 52],
                      ja[i - 51],
                      ja[i - 50], ja[i - 49], ja[i - 48], ja[i - 47], ja[i - 46], ja[i - 45], ja[i - 44], ja[i - 43],
                      ja[i - 42],
                      ja[i - 41], ja[i - 40], ja[i - 39], ja[i - 38], ja[i - 37], ja[i - 36], ja[i - 35], ja[i - 34],
                      ja[i - 33],
                      ja[i - 32], ja[i - 31], ja[i - 30], ja[i - 29], ja[i - 28], ja[i - 27], ja[i - 26], ja[i - 25],
                      ja[i - 24],
                      ja[i - 23], ja[i - 22], ja[i - 21], ja[i - 20], ja[i - 19], ja[i - 18], ja[i - 17], ja[i - 16],
                      ja[i - 15],
                      ja[i - 14], ja[i - 13], ja[i - 12], ja[i - 11], ja[i - 10], ja[i - 9], ja[i - 8], ja[i - 7],
                      ja[i - 6],
                      ja[i - 5], ja[i - 4], ja[i - 3], ja[i - 2], ja[i - 1], ja[i], ja[i + 1], ja[i + 2], ja[i + 3],
                      ja[i + 4], ja[i + 5], ja[i + 6], ja[i + 7], ja[i + 8], ja[i + 9], ja[i + 10], ja[i + 11],
                      ja[i + 12],
                      ja[i + 13], ja[i + 14], ja[i + 15], file=open("output4b.txt", "a+"), sep="")
                strix = open("output4b.txt")
                striy = strix.readlines()
                string = striy[-1]
                res = [sub.replace('_', '') for sub in string]
                x = "".join(res)
                print(x[-17:], file=open("oligob.txt", "a+"))
            else:
                print(ja[i - 1], ja[i], ja[i + 1],
                      ja[i + 2], ja[i + 3], ja[i + 4], ja[i + 5], ja[i + 6], ja[i + 7], ja[i + 8], ja[i + 9],
                      ja[i + 10], ja[i + 11], ja[i + 12], ja[i + 13], ja[i + 14], ja[i + 15],
                      file=open("oligob.txt", "a+"), sep="")
        if i == (ln - 16):
            if ja[i] == "_" or \
                    ja[i + 1] == "_" or ja[i + 2] == "_" or ja[i + 3] == "_" or ja[i + 4] == "_" or ja[i + 5] == "_" or \
                    ja[i + 6] == "_" or ja[i + 7] == "_" or ja[i + 8] == "_" or ja[i + 9] == "_" or ja[i + 10] == "_" or \
                    ja[i + 11] == "_" or ja[i + 12] == "_" or ja[i + 13] == "_" or ja[i + 14] == "_" \
                    or ja[i + 15] == "_" or ja[i + 16] == "_":
                print(ja[i - 52],
                      ja[i - 51],
                      ja[i - 50], ja[i - 49], ja[i - 48], ja[i - 47], ja[i - 46], ja[i - 45], ja[i - 44], ja[i - 43],
                      ja[i - 42],
                      ja[i - 41], ja[i - 40], ja[i - 39], ja[i - 38], ja[i - 37], ja[i - 36], ja[i - 35], ja[i - 34],
                      ja[i - 33],
                      ja[i - 32], ja[i - 31], ja[i - 30], ja[i - 29], ja[i - 28], ja[i - 27], ja[i - 26], ja[i - 25],
                      ja[i - 24],
                      ja[i - 23], ja[i - 22], ja[i - 21], ja[i - 20], ja[i - 19], ja[i - 18], ja[i - 17], ja[i - 16],
                      ja[i - 15],
                      ja[i - 14], ja[i - 13], ja[i - 12], ja[i - 11], ja[i - 10], ja[i - 9], ja[i - 8], ja[i - 7],
                      ja[i - 6],
                      ja[i - 5], ja[i - 4], ja[i - 3], ja[i - 2], ja[i - 1], ja[i], ja[i + 1], ja[i + 2], ja[i + 3],
                      ja[i + 4], ja[i + 5], ja[i + 6], ja[i + 7], ja[i + 8], ja[i + 9], ja[i + 10], ja[i + 11],
                      ja[i + 12],
                      ja[i + 13], ja[i + 14], ja[i + 15], ja[i + 16], file=open("output4b.txt", "a+"), sep="")
                strix = open("output4b.txt")
                striy = strix.readlines()
                string = striy[-1]
                res = [sub.replace('_', '') for sub in string]
                x = "".join(res)
                print(x[-17:], file=open("oligob.txt", "a+"))
            else:
                print(ja[i], ja[i + 1],
                      ja[i + 2], ja[i + 3], ja[i + 4], ja[i + 5], ja[i + 6], ja[i + 7], ja[i + 8], ja[i + 9],
                      ja[i + 10], ja[i + 11], ja[i + 12], ja[i + 13], ja[i + 14], ja[i + 15], ja[i + 16],
                      file=open("oligob.txt", "a+"), sep="")
        if i == (ln - 17):
            if ja[i - 16] == "_" or ja[i - 15] == "_" or ja[i - 14] == "_" or ja[i - 13] == "_" or ja[i - 12] == "_" or \
                    ja[i - 11] == "_" or ja[i - 10] == "_" or ja[i - 9] == "_" or ja[i - 8] == "_" or \
                    ja[i - 7] == "_" or ja[i - 6] == "_" or ja[i - 5] == "_" or ja[i - 4] == "_" or ja[i - 3] == "_" or \
                    ja[i - 2] == "_" or ja[i - 1] == "_" or ja[i] == "_":
                print(ja[i - 68],
                      ja[i - 67],
                      ja[i - 66], ja[i - 65], ja[i - 64], ja[i - 63], ja[i - 62], ja[i - 61], ja[i - 60], ja[i - 59],
                      ja[i - 58],
                      ja[i - 57], ja[i - 56], ja[i - 55], ja[i - 54], ja[i - 53], ja[i - 52], ja[i - 51], ja[i - 50],
                      ja[i - 49],
                      ja[i - 48], ja[i - 47], ja[i - 46], ja[i - 45], ja[i - 44], ja[i - 43], ja[i - 42], ja[i - 41],
                      ja[i - 40],
                      ja[i - 39], ja[i - 38], ja[i - 37], ja[i - 36], ja[i - 35], ja[i - 34], ja[i - 33], ja[i - 32],
                      ja[i - 31],
                      ja[i - 30], ja[i - 29], ja[i - 28], ja[i - 27], ja[i - 26], ja[i - 25], ja[i - 24], ja[i - 23],
                      ja[i - 22],
                      ja[i - 21], ja[i - 20], ja[i - 19], ja[i - 18], ja[i - 17], ja[i - 16], ja[i - 15], ja[i - 14],
                      ja[i - 13],
                      ja[i - 12], ja[i - 11], ja[i - 10], ja[i - 9], ja[i - 8], ja[i - 7], ja[i - 6], ja[i - 5],
                      ja[i - 4],
                      ja[i - 3], ja[i - 2], ja[i - 1], ja[i], file=open("output4b.txt", "a+"), sep="")
                strix = open("output4b.txt")
                striy = strix.readlines()
                string = striy[-1]
                res = [sub.replace('_', '') for sub in string]
                x = "".join(res)
                print(x[-17:], file=open("oligob.txt", "a+"))
            else:
                print(ja[i - 16], ja[i - 15], ja[i - 14], ja[i - 13], ja[i - 12], ja[i - 11], ja[i - 10], ja[i - 9],
                      ja[i - 8], ja[i - 7], ja[i - 6], ja[i - 5], ja[i - 4], ja[i - 3], ja[i - 2], ja[i - 1], ja[i],
                      file=open("oligob.txt", "a+"), sep="")
        if i == (ln - 18):
            if ja[i - 15] == "_" or ja[i - 14] == "_" or ja[i - 13] == "_" or ja[i - 12] == "_" or \
                    ja[i - 11] == "_" or ja[i - 10] == "_" or ja[i - 9] == "_" or ja[i - 8] == "_" or \
                    ja[i - 7] == "_" or ja[i - 6] == "_" or ja[i - 5] == "_" or ja[i - 4] == "_" or ja[i - 3] == "_" or \
                    ja[i - 2] == "_" or ja[i - 1] == "_" or ja[i] == "_" or ja[i + 1] == "_":
                print(ja[i - 67],
                      ja[i - 66], ja[i - 65], ja[i - 64], ja[i - 63], ja[i - 62], ja[i - 61], ja[i - 60], ja[i - 59],
                      ja[i - 58],
                      ja[i - 57], ja[i - 56], ja[i - 55], ja[i - 54], ja[i - 53], ja[i - 52], ja[i - 51], ja[i - 50],
                      ja[i - 49],
                      ja[i - 48], ja[i - 47], ja[i - 46], ja[i - 45], ja[i - 44], ja[i - 43], ja[i - 42], ja[i - 41],
                      ja[i - 40],
                      ja[i - 39], ja[i - 38], ja[i - 37], ja[i - 36], ja[i - 35], ja[i - 34], ja[i - 33], ja[i - 32],
                      ja[i - 31],
                      ja[i - 30], ja[i - 29], ja[i - 28], ja[i - 27], ja[i - 26], ja[i - 25], ja[i - 24], ja[i - 23],
                      ja[i - 22],
                      ja[i - 21], ja[i - 20], ja[i - 19], ja[i - 18], ja[i - 17], ja[i - 16], ja[i - 15], ja[i - 14],
                      ja[i - 13],
                      ja[i - 12], ja[i - 11], ja[i - 10], ja[i - 9], ja[i - 8], ja[i - 7], ja[i - 6], ja[i - 5],
                      ja[i - 4],
                      ja[i - 3], ja[i - 2], ja[i - 1], ja[i], ja[i + 1], file=open("output4b.txt", "a+"), sep="")
                strix = open("output4b.txt")
                striy = strix.readlines()
                string = striy[-1]
                res = [sub.replace('_', '') for sub in string]
                x = "".join(res)
                print(x[-17:], file=open("oligob.txt", "a+"))
            else:
                print(ja[i - 15], ja[i - 14], ja[i - 13], ja[i - 12], ja[i - 11], ja[i - 10], ja[i - 9],
                      ja[i - 8], ja[i - 7], ja[i - 6], ja[i - 5], ja[i - 4], ja[i - 3], ja[i - 2], ja[i - 1], ja[i],
                      ja[i + 1], file=open("oligob.txt", "a+"), sep="")
        if i == (ln - 19):
            if ja[i - 14] == "_" or ja[i - 13] == "_" or ja[i - 12] == "_" or \
                    ja[i - 11] == "_" or ja[i - 10] == "_" or ja[i - 9] == "_" or ja[i - 8] == "_" or \
                    ja[i - 7] == "_" or ja[i - 6] == "_" or ja[i - 5] == "_" or ja[i - 4] == "_" or ja[i - 3] == "_" or \
                    ja[i - 2] == "_" or ja[i - 1] == "_" or ja[i] == "_" or ja[i + 1] == "_" or ja[i + 2] == "_":
                print(ja[i - 66], ja[i - 65], ja[i - 64], ja[i - 63], ja[i - 62], ja[i - 61], ja[i - 60], ja[i - 59],
                      ja[i - 58],
                      ja[i - 57], ja[i - 56], ja[i - 55], ja[i - 54], ja[i - 53], ja[i - 52], ja[i - 51], ja[i - 50],
                      ja[i - 49],
                      ja[i - 48], ja[i - 47], ja[i - 46], ja[i - 45], ja[i - 44], ja[i - 43], ja[i - 42], ja[i - 41],
                      ja[i - 40],
                      ja[i - 39], ja[i - 38], ja[i - 37], ja[i - 36], ja[i - 35], ja[i - 34], ja[i - 33], ja[i - 32],
                      ja[i - 31],
                      ja[i - 30], ja[i - 29], ja[i - 28], ja[i - 27], ja[i - 26], ja[i - 25], ja[i - 24], ja[i - 23],
                      ja[i - 22],
                      ja[i - 21], ja[i - 20], ja[i - 19], ja[i - 18], ja[i - 17], ja[i - 16], ja[i - 15], ja[i - 14],
                      ja[i - 13],
                      ja[i - 12], ja[i - 11], ja[i - 10], ja[i - 9], ja[i - 8], ja[i - 7], ja[i - 6], ja[i - 5],
                      ja[i - 4],
                      ja[i - 3], ja[i - 2], ja[i - 1], ja[i], ja[i + 1], ja[i + 2], file=open("output4b.txt", "a+"),
                      sep="")
                strix = open("output4b.txt")
                striy = strix.readlines()
                string = striy[-1]
                res = [sub.replace('_', '') for sub in string]
                x = "".join(res)
                print(x[-17:], file=open("oligob.txt", "a+"))
            else:
                print(ja[i - 14], ja[i - 13], ja[i - 12], ja[i - 11], ja[i - 10], ja[i - 9],
                      ja[i - 8], ja[i - 7], ja[i - 6], ja[i - 5], ja[i - 4], ja[i - 3], ja[i - 2], ja[i - 1], ja[i],
                      ja[i + 1], ja[i + 2], file=open("oligob.txt", "a+"), sep="")
        if i == (ln - 20):
            if ja[i - 13] == "_" or ja[i - 12] == "_" or \
                    ja[i - 11] == "_" or ja[i - 10] == "_" or ja[i - 9] == "_" or ja[i - 8] == "_" or \
                    ja[i - 7] == "_" or ja[i - 6] == "_" or ja[i - 5] == "_" or ja[i - 4] == "_" or ja[i - 3] == "_" or \
                    ja[i - 2] == "_" or ja[i - 1] == "_" or ja[i] == "_" or ja[i + 1] == "_" or ja[i + 2] == "_" or \
                    ja[i + 3] == "_":
                print(ja[i - 65], ja[i - 64], ja[i - 63], ja[i - 62], ja[i - 61], ja[i - 60], ja[i - 59],
                      ja[i - 58],
                      ja[i - 57], ja[i - 56], ja[i - 55], ja[i - 54], ja[i - 53], ja[i - 52], ja[i - 51], ja[i - 50],
                      ja[i - 49],
                      ja[i - 48], ja[i - 47], ja[i - 46], ja[i - 45], ja[i - 44], ja[i - 43], ja[i - 42], ja[i - 41],
                      ja[i - 40],
                      ja[i - 39], ja[i - 38], ja[i - 37], ja[i - 36], ja[i - 35], ja[i - 34], ja[i - 33], ja[i - 32],
                      ja[i - 31],
                      ja[i - 30], ja[i - 29], ja[i - 28], ja[i - 27], ja[i - 26], ja[i - 25], ja[i - 24], ja[i - 23],
                      ja[i - 22],
                      ja[i - 21], ja[i - 20], ja[i - 19], ja[i - 18], ja[i - 17], ja[i - 16], ja[i - 15], ja[i - 14],
                      ja[i - 13],
                      ja[i - 12], ja[i - 11], ja[i - 10], ja[i - 9], ja[i - 8], ja[i - 7], ja[i - 6], ja[i - 5],
                      ja[i - 4],
                      ja[i - 3], ja[i - 2], ja[i - 1], ja[i], ja[i + 1], ja[i + 2], ja[i + 3],
                      file=open("output4b.txt", "a+"), sep="")
                strix = open("output4b.txt")
                striy = strix.readlines()
                string = striy[-1]
                res = [sub.replace('_', '') for sub in string]
                x = "".join(res)
                print(x[-17:], file=open("oligob.txt", "a+"))
            else:
                print(ja[i - 13], ja[i - 12], ja[i - 11], ja[i - 10], ja[i - 9],
                      ja[i - 8], ja[i - 7], ja[i - 6], ja[i - 5], ja[i - 4], ja[i - 3], ja[i - 2], ja[i - 1], ja[i],
                      ja[i + 1], ja[i + 2], ja[i + 3], file=open("oligob.txt", "a+"), sep="")
        if i == (ln - 21):
            if ja[i - 12] == "_" or \
                    ja[i - 11] == "_" or ja[i - 10] == "_" or ja[i - 9] == "_" or ja[i - 8] == "_" or \
                    ja[i - 7] == "_" or ja[i - 6] == "_" or ja[i - 5] == "_" or ja[i - 4] == "_" or ja[i - 3] == "_" or \
                    ja[i - 2] == "_" or ja[i - 1] == "_" or ja[i] == "_" or ja[i + 1] == "_" or ja[i + 2] == "_" or \
                    ja[i + 3] == "_" or ja[i + 4] == "_":
                print(ja[i - 64], ja[i - 63], ja[i - 62], ja[i - 61], ja[i - 60], ja[i - 59],
                      ja[i - 58],
                      ja[i - 57], ja[i - 56], ja[i - 55], ja[i - 54], ja[i - 53], ja[i - 52], ja[i - 51], ja[i - 50],
                      ja[i - 49],
                      ja[i - 48], ja[i - 47], ja[i - 46], ja[i - 45], ja[i - 44], ja[i - 43], ja[i - 42], ja[i - 41],
                      ja[i - 40],
                      ja[i - 39], ja[i - 38], ja[i - 37], ja[i - 36], ja[i - 35], ja[i - 34], ja[i - 33], ja[i - 32],
                      ja[i - 31],
                      ja[i - 30], ja[i - 29], ja[i - 28], ja[i - 27], ja[i - 26], ja[i - 25], ja[i - 24], ja[i - 23],
                      ja[i - 22],
                      ja[i - 21], ja[i - 20], ja[i - 19], ja[i - 18], ja[i - 17], ja[i - 16], ja[i - 15], ja[i - 14],
                      ja[i - 13],
                      ja[i - 12], ja[i - 11], ja[i - 10], ja[i - 9], ja[i - 8], ja[i - 7], ja[i - 6], ja[i - 5],
                      ja[i - 4],
                      ja[i - 3], ja[i - 2], ja[i - 1], ja[i], ja[i + 1], ja[i + 2], ja[i + 3], ja[i + 4],
                      file=open("output4b.txt", "a+"), sep="")
                strix = open("output4b.txt")
                striy = strix.readlines()
                string = striy[-1]
                res = [sub.replace('_', '') for sub in string]
                x = "".join(res)
                print(x[-17:], file=open("oligob.txt", "a+"))
            else:
                print(ja[i - 12], ja[i - 11], ja[i - 10], ja[i - 9],
                      ja[i - 8], ja[i - 7], ja[i - 6], ja[i - 5], ja[i - 4], ja[i - 3], ja[i - 2], ja[i - 1], ja[i],
                      ja[i + 1], ja[i + 2], ja[i + 3], ja[i + 4], file=open("oligob.txt", "a+"), sep="")
        if i == (ln - 22):
            if ja[i - 11] == "_" or ja[i - 10] == "_" or ja[i - 9] == "_" or ja[i - 8] == "_" or \
                    ja[i - 7] == "_" or ja[i - 6] == "_" or ja[i - 5] == "_" or ja[i - 4] == "_" or ja[i - 3] == "_" or \
                    ja[i - 2] == "_" or ja[i - 1] == "_" or ja[i] == "_" or ja[i + 1] == "_" or ja[i + 2] == "_" or \
                    ja[i + 3] == "_" or ja[i + 4] == "_" or ja[i + 5] == "_":
                print(ja[i - 63], ja[i - 62], ja[i - 61], ja[i - 60], ja[i - 59],
                      ja[i - 58],
                      ja[i - 57], ja[i - 56], ja[i - 55], ja[i - 54], ja[i - 53], ja[i - 52], ja[i - 51], ja[i - 50],
                      ja[i - 49],
                      ja[i - 48], ja[i - 47], ja[i - 46], ja[i - 45], ja[i - 44], ja[i - 43], ja[i - 42], ja[i - 41],
                      ja[i - 40],
                      ja[i - 39], ja[i - 38], ja[i - 37], ja[i - 36], ja[i - 35], ja[i - 34], ja[i - 33], ja[i - 32],
                      ja[i - 31],
                      ja[i - 30], ja[i - 29], ja[i - 28], ja[i - 27], ja[i - 26], ja[i - 25], ja[i - 24], ja[i - 23],
                      ja[i - 22],
                      ja[i - 21], ja[i - 20], ja[i - 19], ja[i - 18], ja[i - 17], ja[i - 16], ja[i - 15], ja[i - 14],
                      ja[i - 13],
                      ja[i - 12], ja[i - 11], ja[i - 10], ja[i - 9], ja[i - 8], ja[i - 7], ja[i - 6], ja[i - 5],
                      ja[i - 4],
                      ja[i - 3], ja[i - 2], ja[i - 1], ja[i], ja[i + 1], ja[i + 2], ja[i + 3], ja[i + 4], ja[i + 5],
                      file=open("output4b.txt", "a+"), sep="")
                strix = open("output4b.txt")
                striy = strix.readlines()
                string = striy[-1]
                res = [sub.replace('_', '') for sub in string]
                x = "".join(res)
                print(x[-17:], file=open("oligob.txt", "a+"))
            else:
                print(ja[i - 11], ja[i - 10], ja[i - 9],
                      ja[i - 8], ja[i - 7], ja[i - 6], ja[i - 5], ja[i - 4], ja[i - 3], ja[i - 2], ja[i - 1], ja[i],
                      ja[i + 1], ja[i + 2], ja[i + 3], ja[i + 4], ja[i + 5], file=open("oligob.txt", "a+"), sep="")
        if i == (ln - 23):
            if ja[i - 10] == "_" or ja[i - 9] == "_" or ja[i - 8] == "_" or \
                    ja[i - 7] == "_" or ja[i - 6] == "_" or ja[i - 5] == "_" or ja[i - 4] == "_" or ja[i - 3] == "_" or \
                    ja[i - 2] == "_" or ja[i - 1] == "_" or ja[i] == "_" or ja[i + 1] == "_" or ja[i + 2] == "_" or \
                    ja[i + 3] == "_" or ja[i + 4] == "_" or ja[i + 5] == "_" or ja[i + 6] == "_":
                print(ja[i - 62], ja[i - 61], ja[i - 60], ja[i - 59],
                      ja[i - 58],
                      ja[i - 57], ja[i - 56], ja[i - 55], ja[i - 54], ja[i - 53], ja[i - 52], ja[i - 51], ja[i - 50],
                      ja[i - 49],
                      ja[i - 48], ja[i - 47], ja[i - 46], ja[i - 45], ja[i - 44], ja[i - 43], ja[i - 42], ja[i - 41],
                      ja[i - 40],
                      ja[i - 39], ja[i - 38], ja[i - 37], ja[i - 36], ja[i - 35], ja[i - 34], ja[i - 33], ja[i - 32],
                      ja[i - 31],
                      ja[i - 30], ja[i - 29], ja[i - 28], ja[i - 27], ja[i - 26], ja[i - 25], ja[i - 24], ja[i - 23],
                      ja[i - 22],
                      ja[i - 21], ja[i - 20], ja[i - 19], ja[i - 18], ja[i - 17], ja[i - 16], ja[i - 15], ja[i - 14],
                      ja[i - 13],
                      ja[i - 12], ja[i - 11], ja[i - 10], ja[i - 9], ja[i - 8], ja[i - 7], ja[i - 6], ja[i - 5],
                      ja[i - 4],
                      ja[i - 3], ja[i - 2], ja[i - 1], ja[i], ja[i + 1], ja[i + 2], ja[i + 3], ja[i + 4], ja[i + 5],
                      ja[i + 6],
                      file=open("output4b.txt", "a+"), sep="")
                strix = open("output4b.txt")
                striy = strix.readlines()
                string = striy[-1]
                res = [sub.replace('_', '') for sub in string]
                x = "".join(res)
                print(x[-17:], file=open("oligob.txt", "a+"))
            else:
                print(ja[i - 10], ja[i - 9],
                      ja[i - 8], ja[i - 7], ja[i - 6], ja[i - 5], ja[i - 4], ja[i - 3], ja[i - 2], ja[i - 1], ja[i],
                      ja[i + 1], ja[i + 2], ja[i + 3], ja[i + 4], ja[i + 5], ja[i + 6], file=open("oligob.txt", "a+"),
                      sep="")
        if i == (ln - 24):
            if ja[i - 9] == "_" or ja[i - 8] == "_" or \
                    ja[i - 7] == "_" or ja[i - 6] == "_" or ja[i - 5] == "_" or ja[i - 4] == "_" or ja[i - 3] == "_" or \
                    ja[i - 2] == "_" or ja[i - 1] == "_" or ja[i] == "_" or ja[i + 1] == "_" or ja[i + 2] == "_" or \
                    ja[i + 3] == "_" or ja[i + 4] == "_" or ja[i + 5] == "_" or ja[i + 6] == "_" or ja[i + 7] == "_":
                print(ja[i - 61], ja[i - 60], ja[i - 59],
                      ja[i - 58],
                      ja[i - 57], ja[i - 56], ja[i - 55], ja[i - 54], ja[i - 53], ja[i - 52], ja[i - 51], ja[i - 50],
                      ja[i - 49],
                      ja[i - 48], ja[i - 47], ja[i - 46], ja[i - 45], ja[i - 44], ja[i - 43], ja[i - 42], ja[i - 41],
                      ja[i - 40],
                      ja[i - 39], ja[i - 38], ja[i - 37], ja[i - 36], ja[i - 35], ja[i - 34], ja[i - 33], ja[i - 32],
                      ja[i - 31],
                      ja[i - 30], ja[i - 29], ja[i - 28], ja[i - 27], ja[i - 26], ja[i - 25], ja[i - 24], ja[i - 23],
                      ja[i - 22],
                      ja[i - 21], ja[i - 20], ja[i - 19], ja[i - 18], ja[i - 17], ja[i - 16], ja[i - 15], ja[i - 14],
                      ja[i - 13],
                      ja[i - 12], ja[i - 11], ja[i - 10], ja[i - 9], ja[i - 8], ja[i - 7], ja[i - 6], ja[i - 5],
                      ja[i - 4],
                      ja[i - 3], ja[i - 2], ja[i - 1], ja[i], ja[i + 1], ja[i + 2], ja[i + 3], ja[i + 4],
                      ja[i + 5], ja[i + 6], ja[i + 7], file=open("output4b.txt", "a+"), sep="")
                strix = open("output4b.txt")
                striy = strix.readlines()
                string = striy[-1]
                res = [sub.replace('_', '') for sub in string]
                x = "".join(res)
                print(x[-17:], file=open("oligob.txt", "a+"))
            else:
                print(ja[i - 9], ja[i - 8], ja[i - 7], ja[i - 6], ja[i - 5], ja[i - 4], ja[i - 3], ja[i - 2],
                      ja[i - 1], ja[i], ja[i + 1], ja[i + 2], ja[i + 3], ja[i + 4], ja[i + 5], ja[i + 6],
                      ja[i + 7], file=open("oligob.txt", "a+"), sep="")
        if i == (ln - 25):
            if ja[i - 8] == "_" or \
                    ja[i - 7] == "_" or ja[i - 6] == "_" or ja[i - 5] == "_" or ja[i - 4] == "_" or ja[i - 3] == "_" or \
                    ja[i - 2] == "_" or ja[i - 1] == "_" or ja[i] == "_" or ja[i + 1] == "_" or ja[i + 2] == "_" or \
                    ja[i + 3] == "_" or ja[i + 4] == "_" or ja[i + 5] == "_" or ja[i + 6] == "_" or ja[i + 7] == "_" or \
                    ja[i + 8] == "_":
                print(ja[i - 60], ja[i - 59],
                      ja[i - 58],
                      ja[i - 57], ja[i - 56], ja[i - 55], ja[i - 54], ja[i - 53], ja[i - 52], ja[i - 51], ja[i - 50],
                      ja[i - 49],
                      ja[i - 48], ja[i - 47], ja[i - 46], ja[i - 45], ja[i - 44], ja[i - 43], ja[i - 42], ja[i - 41],
                      ja[i - 40],
                      ja[i - 39], ja[i - 38], ja[i - 37], ja[i - 36], ja[i - 35], ja[i - 34], ja[i - 33], ja[i - 32],
                      ja[i - 31],
                      ja[i - 30], ja[i - 29], ja[i - 28], ja[i - 27], ja[i - 26], ja[i - 25], ja[i - 24], ja[i - 23],
                      ja[i - 22],
                      ja[i - 21], ja[i - 20], ja[i - 19], ja[i - 18], ja[i - 17], ja[i - 16], ja[i - 15], ja[i - 14],
                      ja[i - 13],
                      ja[i - 12], ja[i - 11], ja[i - 10], ja[i - 9], ja[i - 8], ja[i - 7], ja[i - 6], ja[i - 5],
                      ja[i - 4],
                      ja[i - 3], ja[i - 2], ja[i - 1], ja[i], ja[i + 1], ja[i + 2], ja[i + 3], ja[i + 4],
                      ja[i + 5], ja[i + 6], ja[i + 7], ja[i + 8], file=open("output4b.txt", "a+"), sep="")
                strix = open("output4b.txt")
                striy = strix.readlines()
                string = striy[-1]
                res = [sub.replace('_', '') for sub in string]
                x = "".join(res)
                print(x[-17:], file=open("oligob.txt", "a+"))
            else:
                print(ja[i - 8], ja[i - 7], ja[i - 6], ja[i - 5], ja[i - 4], ja[i - 3], ja[i - 2],
                      ja[i - 1], ja[i], ja[i + 1], ja[i + 2], ja[i + 3], ja[i + 4], ja[i + 5], ja[i + 6],
                      ja[i + 7], ja[i + 8], file=open("oligob.txt", "a+"), sep="")
        if i == (ln - 26):
            if ja[i - 7] == "_" or ja[i - 6] == "_" or ja[i - 5] == "_" or ja[i - 4] == "_" or ja[i - 3] == "_" or \
                    ja[i - 2] == "_" or ja[i - 1] == "_" or ja[i] == "_" or ja[i + 1] == "_" or ja[i + 2] == "_" or \
                    ja[i + 3] == "_" or ja[i + 4] == "_" or ja[i + 5] == "_" or ja[i + 6] == "_" or ja[i + 7] == "_" or \
                    ja[i + 8] == "_" or ja[i + 9] == "_":
                print(ja[i - 59],
                      ja[i - 58],
                      ja[i - 57], ja[i - 56], ja[i - 55], ja[i - 54], ja[i - 53], ja[i - 52], ja[i - 51], ja[i - 50],
                      ja[i - 49],
                      ja[i - 48], ja[i - 47], ja[i - 46], ja[i - 45], ja[i - 44], ja[i - 43], ja[i - 42], ja[i - 41],
                      ja[i - 40],
                      ja[i - 39], ja[i - 38], ja[i - 37], ja[i - 36], ja[i - 35], ja[i - 34], ja[i - 33], ja[i - 32],
                      ja[i - 31],
                      ja[i - 30], ja[i - 29], ja[i - 28], ja[i - 27], ja[i - 26], ja[i - 25], ja[i - 24], ja[i - 23],
                      ja[i - 22],
                      ja[i - 21], ja[i - 20], ja[i - 19], ja[i - 18], ja[i - 17], ja[i - 16], ja[i - 15], ja[i - 14],
                      ja[i - 13],
                      ja[i - 12], ja[i - 11], ja[i - 10], ja[i - 9], ja[i - 8], ja[i - 7], ja[i - 6], ja[i - 5],
                      ja[i - 4],
                      ja[i - 3], ja[i - 2], ja[i - 1], ja[i], ja[i + 1], ja[i + 2], ja[i + 3], ja[i + 4],
                      ja[i + 5], ja[i + 6], ja[i + 7], ja[i + 8], ja[i + 9], file=open("output4b.txt", "a+"), sep="")
                strix = open("output4b.txt")
                striy = strix.readlines()
                string = striy[-1]
                res = [sub.replace('_', '') for sub in string]
                x = "".join(res)
                print(x[-17:], file=open("oligob.txt", "a+"))
            else:
                print(ja[i - 7], ja[i - 6], ja[i - 5], ja[i - 4], ja[i - 3], ja[i - 2],
                      ja[i - 1], ja[i], ja[i + 1], ja[i + 2], ja[i + 3], ja[i + 4], ja[i + 5], ja[i + 6],
                      ja[i + 7], ja[i + 8], ja[i + 9], file=open("oligob.txt", "a+"), sep="")
        if i == (ln - 27):
            if ja[i - 6] == "_" or ja[i - 5] == "_" or ja[i - 4] == "_" or ja[i - 3] == "_" or \
                    ja[i - 2] == "_" or ja[i - 1] == "_" or ja[i] == "_" or ja[i + 1] == "_" or ja[i + 2] == "_" or \
                    ja[i + 3] == "_" or ja[i + 4] == "_" or ja[i + 5] == "_" or ja[i + 6] == "_" or ja[i + 7] == "_" or \
                    ja[i + 8] == "_" or ja[i + 9] == "_" or ja[i + 10] == "_":
                print(ja[i - 58],
                      ja[i - 57], ja[i - 56], ja[i - 55], ja[i - 54], ja[i - 53], ja[i - 52], ja[i - 51], ja[i - 50],
                      ja[i - 49],
                      ja[i - 48], ja[i - 47], ja[i - 46], ja[i - 45], ja[i - 44], ja[i - 43], ja[i - 42], ja[i - 41],
                      ja[i - 40],
                      ja[i - 39], ja[i - 38], ja[i - 37], ja[i - 36], ja[i - 35], ja[i - 34], ja[i - 33], ja[i - 32],
                      ja[i - 31],
                      ja[i - 30], ja[i - 29], ja[i - 28], ja[i - 27], ja[i - 26], ja[i - 25], ja[i - 24], ja[i - 23],
                      ja[i - 22],
                      ja[i - 21], ja[i - 20], ja[i - 19], ja[i - 18], ja[i - 17], ja[i - 16], ja[i - 15], ja[i - 14],
                      ja[i - 13],
                      ja[i - 12], ja[i - 11], ja[i - 10], ja[i - 9], ja[i - 8], ja[i - 7], ja[i - 6], ja[i - 5],
                      ja[i - 4],
                      ja[i - 3], ja[i - 2], ja[i - 1], ja[i], ja[i + 1], ja[i + 2], ja[i + 3], ja[i + 4],
                      ja[i + 5], ja[i + 6], ja[i + 7], ja[i + 8], ja[i + 9], ja[i + 10],
                      file=open("output4b.txt", "a+"), sep="")
                strix = open("output4b.txt")
                striy = strix.readlines()
                string = striy[-1]
                res = [sub.replace('_', '') for sub in string]
                x = "".join(res)
                print(x[-17:], file=open("oligob.txt", "a+"))
            else:
                print(ja[i - 6], ja[i - 5], ja[i - 4], ja[i - 3], ja[i - 2],
                      ja[i - 1], ja[i], ja[i + 1], ja[i + 2], ja[i + 3], ja[i + 4], ja[i + 5], ja[i + 6],
                      ja[i + 7], ja[i + 8], ja[i + 9], ja[i + 10], file=open("oligob.txt", "a+"), sep="")
        if i == (ln - 28):
            if ja[i - 5] == "_" or ja[i - 4] == "_" or ja[i - 3] == "_" or \
                    ja[i - 2] == "_" or ja[i - 1] == "_" or ja[i] == "_" or ja[i + 1] == "_" or ja[i + 2] == "_" or \
                    ja[i + 3] == "_" or ja[i + 4] == "_" or ja[i + 5] == "_" or ja[i + 6] == "_" or ja[i + 7] == "_" or \
                    ja[i + 8] == "_" or ja[i + 9] == "_" or ja[i + 10] == "_" or ja[i + 11] == "_":
                print(ja[i - 57], ja[i - 56], ja[i - 55], ja[i - 54], ja[i - 53], ja[i - 52], ja[i - 51], ja[i - 50],
                      ja[i - 49],
                      ja[i - 48], ja[i - 47], ja[i - 46], ja[i - 45], ja[i - 44], ja[i - 43], ja[i - 42], ja[i - 41],
                      ja[i - 40],
                      ja[i - 39], ja[i - 38], ja[i - 37], ja[i - 36], ja[i - 35], ja[i - 34], ja[i - 33], ja[i - 32],
                      ja[i - 31],
                      ja[i - 30], ja[i - 29], ja[i - 28], ja[i - 27], ja[i - 26], ja[i - 25], ja[i - 24], ja[i - 23],
                      ja[i - 22],
                      ja[i - 21], ja[i - 20], ja[i - 19], ja[i - 18], ja[i - 17], ja[i - 16], ja[i - 15], ja[i - 14],
                      ja[i - 13],
                      ja[i - 12], ja[i - 11], ja[i - 10], ja[i - 9], ja[i - 8], ja[i - 7], ja[i - 6], ja[i - 5],
                      ja[i - 4],
                      ja[i - 3], ja[i - 2], ja[i - 1], ja[i], ja[i + 1], ja[i + 2], ja[i + 3], ja[i + 4],
                      ja[i + 5], ja[i + 6], ja[i + 7], ja[i + 8], ja[i + 9], ja[i + 10], ja[i + 11],
                      file=open("output4b.txt", "a+"), sep="")
                strix = open("output4b.txt")
                striy = strix.readlines()
                string = striy[-1]
                res = [sub.replace('_', '') for sub in string]
                x = "".join(res)
                print(x[-17:], file=open("oligob.txt", "a+"))
            else:
                print(ja[i - 5], ja[i - 4], ja[i - 3], ja[i - 2],
                      ja[i - 1], ja[i], ja[i + 1], ja[i + 2], ja[i + 3], ja[i + 4], ja[i + 5], ja[i + 6],
                      ja[i + 7], ja[i + 8], ja[i + 9], ja[i + 10], ja[i + 11], file=open("oligob.txt", "a+"), sep="")
        if i == (ln - 29):
            if ja[i - 4] == "_" or ja[i - 3] == "_" or \
                    ja[i - 2] == "_" or ja[i - 1] == "_" or ja[i] == "_" or ja[i + 1] == "_" or ja[i + 2] == "_" or \
                    ja[i + 3] == "_" or ja[i + 4] == "_" or ja[i + 5] == "_" or ja[i + 6] == "_" or ja[i + 7] == "_" or \
                    ja[i + 8] == "_" or ja[i + 9] == "_" or ja[i + 10] == "_" or ja[i + 11] == "_" or ja[i + 12] == "_":
                print(ja[i - 56], ja[i - 55], ja[i - 54], ja[i - 53], ja[i - 52], ja[i - 51], ja[i - 50],
                      ja[i - 49],
                      ja[i - 48], ja[i - 47], ja[i - 46], ja[i - 45], ja[i - 44], ja[i - 43], ja[i - 42], ja[i - 41],
                      ja[i - 40],
                      ja[i - 39], ja[i - 38], ja[i - 37], ja[i - 36], ja[i - 35], ja[i - 34], ja[i - 33], ja[i - 32],
                      ja[i - 31],
                      ja[i - 30], ja[i - 29], ja[i - 28], ja[i - 27], ja[i - 26], ja[i - 25], ja[i - 24], ja[i - 23],
                      ja[i - 22],
                      ja[i - 21], ja[i - 20], ja[i - 19], ja[i - 18], ja[i - 17], ja[i - 16], ja[i - 15], ja[i - 14],
                      ja[i - 13],
                      ja[i - 12], ja[i - 11], ja[i - 10], ja[i - 9], ja[i - 8], ja[i - 7], ja[i - 6], ja[i - 5],
                      ja[i - 4],
                      ja[i - 3], ja[i - 2], ja[i - 1], ja[i], ja[i + 1], ja[i + 2], ja[i + 3], ja[i + 4],
                      ja[i + 5], ja[i + 6], ja[i + 7], ja[i + 8], ja[i + 9], ja[i + 10], ja[i + 11], ja[i + 12],
                      file=open("output4b.txt", "a+"), sep="")
                strix = open("output4b.txt")
                striy = strix.readlines()
                string = striy[-1]
                res = [sub.replace('_', '') for sub in string]
                x = "".join(res)
                print(x[-17:], file=open("oligob.txt", "a+"))
            else:
                print(ja[i - 4], ja[i - 3], ja[i - 2],
                      ja[i - 1], ja[i], ja[i + 1], ja[i + 2], ja[i + 3], ja[i + 4], ja[i + 5], ja[i + 6],
                      ja[i + 7], ja[i + 8], ja[i + 9], ja[i + 10], ja[i + 11], ja[i + 12],
                      file=open("oligob.txt", "a+"), sep="")
        if i == (ln - 30):
            if ja[i - 3] == "_" or \
                    ja[i - 2] == "_" or ja[i - 1] == "_" or ja[i] == "_" or ja[i + 1] == "_" or ja[i + 2] == "_" or \
                    ja[i + 3] == "_" or ja[i + 4] == "_" or ja[i + 5] == "_" or ja[i + 6] == "_" or ja[i + 7] == "_" or \
                    ja[i + 8] == "_" or ja[i + 9] == "_" or ja[i + 10] == "_" or ja[i + 11] == "_" or ja[i + 12] == "_" \
                    or ja[i + 13] == "_":
                print(ja[i - 55], ja[i - 54], ja[i - 53], ja[i - 52], ja[i - 51], ja[i - 50],
                      ja[i - 49],
                      ja[i - 48], ja[i - 47], ja[i - 46], ja[i - 45], ja[i - 44], ja[i - 43], ja[i - 42], ja[i - 41],
                      ja[i - 40],
                      ja[i - 39], ja[i - 38], ja[i - 37], ja[i - 36], ja[i - 35], ja[i - 34], ja[i - 33], ja[i - 32],
                      ja[i - 31],
                      ja[i - 30], ja[i - 29], ja[i - 28], ja[i - 27], ja[i - 26], ja[i - 25], ja[i - 24], ja[i - 23],
                      ja[i - 22],
                      ja[i - 21], ja[i - 20], ja[i - 19], ja[i - 18], ja[i - 17], ja[i - 16], ja[i - 15], ja[i - 14],
                      ja[i - 13],
                      ja[i - 12], ja[i - 11], ja[i - 10], ja[i - 9], ja[i - 8], ja[i - 7], ja[i - 6], ja[i - 5],
                      ja[i - 4],
                      ja[i - 3], ja[i - 2], ja[i - 1], ja[i], ja[i + 1], ja[i + 2], ja[i + 3], ja[i + 4],
                      ja[i + 5], ja[i + 6], ja[i + 7], ja[i + 8], ja[i + 9], ja[i + 10], ja[i + 11], ja[i + 12],
                      ja[i + 13], file=open("output4b.txt", "a+"), sep="")
                strix = open("output4b.txt")
                striy = strix.readlines()
                string = striy[-1]
                res = [sub.replace('_', '') for sub in string]
                x = "".join(res)
                print(x[-17:], file=open("oligob.txt", "a+"))
            else:
                print(ja[i - 3], ja[i - 2],
                      ja[i - 1], ja[i], ja[i + 1], ja[i + 2], ja[i + 3], ja[i + 4], ja[i + 5], ja[i + 6],
                      ja[i + 7], ja[i + 8], ja[i + 9], ja[i + 10], ja[i + 11], ja[i + 12], ja[i + 13],
                      file=open("oligob.txt", "a+"), sep="")
        if i == (ln - 31):
            if ja[i - 2] == "_" or ja[i - 1] == "_" or ja[i] == "_" or ja[i + 1] == "_" or ja[i + 2] == "_" or \
                    ja[i + 3] == "_" or ja[i + 4] == "_" or ja[i + 5] == "_" or ja[i + 6] == "_" or ja[i + 7] == "_" or \
                    ja[i + 8] == "_" or ja[i + 9] == "_" or ja[i + 10] == "_" or ja[i + 11] == "_" or ja[i + 12] == "_" \
                    or ja[i + 13] == "_" or ja[i + 14] == "_":
                print(ja[i - 54], ja[i - 53], ja[i - 52], ja[i - 51], ja[i - 50],
                      ja[i - 49],
                      ja[i - 48], ja[i - 47], ja[i - 46], ja[i - 45], ja[i - 44], ja[i - 43], ja[i - 42], ja[i - 41],
                      ja[i - 40],
                      ja[i - 39], ja[i - 38], ja[i - 37], ja[i - 36], ja[i - 35], ja[i - 34], ja[i - 33], ja[i - 32],
                      ja[i - 31],
                      ja[i - 30], ja[i - 29], ja[i - 28], ja[i - 27], ja[i - 26], ja[i - 25], ja[i - 24], ja[i - 23],
                      ja[i - 22],
                      ja[i - 21], ja[i - 20], ja[i - 19], ja[i - 18], ja[i - 17], ja[i - 16], ja[i - 15], ja[i - 14],
                      ja[i - 13],
                      ja[i - 12], ja[i - 11], ja[i - 10], ja[i - 9], ja[i - 8], ja[i - 7], ja[i - 6], ja[i - 5],
                      ja[i - 4],
                      ja[i - 3], ja[i - 2], ja[i - 1], ja[i], ja[i + 1], ja[i + 2], ja[i + 3], ja[i + 4],
                      ja[i + 5], ja[i + 6], ja[i + 7], ja[i + 8], ja[i + 9], ja[i + 10], ja[i + 11], ja[i + 12],
                      ja[i + 13], ja[i + 14], file=open("output4b.txt", "a+"), sep="")
                strix = open("output4b.txt")
                striy = strix.readlines()
                string = striy[-1]
                res = [sub.replace('_', '') for sub in string]
                x = "".join(res)
                print(x[-17:], file=open("oligob.txt", "a+"))
            else:
                print(ja[i - 2],
                      ja[i - 1], ja[i], ja[i + 1], ja[i + 2], ja[i + 3], ja[i + 4], ja[i + 5], ja[i + 6],
                      ja[i + 7], ja[i + 8], ja[i + 9], ja[i + 10], ja[i + 11], ja[i + 12], ja[i + 13], ja[i + 14],
                      file=open("oligob.txt", "a+"), sep="")
        if i == (ln - 32):
            if ja[i - 1] == "_" or ja[i] == "_" or ja[i + 1] == "_" or ja[i + 2] == "_" or \
                    ja[i + 3] == "_" or ja[i + 4] == "_" or ja[i + 5] == "_" or ja[i + 6] == "_" or ja[i + 7] == "_" or \
                    ja[i + 8] == "_" or ja[i + 9] == "_" or ja[i + 10] == "_" or ja[i + 11] == "_" or ja[i + 12] == "_" \
                    or ja[i + 13] == "_" or ja[i + 14] == "_" or ja[i + 15] == "_":
                print(ja[i - 53], ja[i - 52], ja[i - 51], ja[i - 50],
                      ja[i - 49],
                      ja[i - 48], ja[i - 47], ja[i - 46], ja[i - 45], ja[i - 44], ja[i - 43], ja[i - 42], ja[i - 41],
                      ja[i - 40],
                      ja[i - 39], ja[i - 38], ja[i - 37], ja[i - 36], ja[i - 35], ja[i - 34], ja[i - 33], ja[i - 32],
                      ja[i - 31],
                      ja[i - 30], ja[i - 29], ja[i - 28], ja[i - 27], ja[i - 26], ja[i - 25], ja[i - 24], ja[i - 23],
                      ja[i - 22],
                      ja[i - 21], ja[i - 20], ja[i - 19], ja[i - 18], ja[i - 17], ja[i - 16], ja[i - 15], ja[i - 14],
                      ja[i - 13],
                      ja[i - 12], ja[i - 11], ja[i - 10], ja[i - 9], ja[i - 8], ja[i - 7], ja[i - 6], ja[i - 5],
                      ja[i - 4],
                      ja[i - 3], ja[i - 2], ja[i - 1], ja[i], ja[i + 1], ja[i + 2], ja[i + 3], ja[i + 4],
                      ja[i + 5], ja[i + 6], ja[i + 7], ja[i + 8], ja[i + 9], ja[i + 10], ja[i + 11], ja[i + 12],
                      ja[i + 13], ja[i + 14], ja[i + 15], file=open("output4b.txt", "a+"), sep="")
                strix = open("output4b.txt")
                striy = strix.readlines()
                string = striy[-1]
                res = [sub.replace('_', '') for sub in string]
                x = "".join(res)
                print(x[-17:], file=open("oligob.txt", "a+"))
            else:
                print(ja[i - 1], ja[i], ja[i + 1], ja[i + 2], ja[i + 3], ja[i + 4], ja[i + 5], ja[i + 6],
                      ja[i + 7], ja[i + 8], ja[i + 9], ja[i + 10], ja[i + 11], ja[i + 12], ja[i + 13], ja[i + 14],
                      ja[i + 15], file=open("oligob.txt", "a+"), sep="")
        if i == (ln - 33):
            if ja[i] == "_" or ja[i + 1] == "_" or ja[i + 2] == "_" or \
                    ja[i + 3] == "_" or ja[i + 4] == "_" or ja[i + 5] == "_" or ja[i + 6] == "_" or ja[i + 7] == "_" or \
                    ja[i + 8] == "_" or ja[i + 9] == "_" or ja[i + 10] == "_" or ja[i + 11] == "_" or ja[i + 12] == "_" \
                    or ja[i + 13] == "_" or ja[i + 14] == "_" or ja[i + 15] == "_" or ja[i + 16] == "_":
                print(ja[i - 52], ja[i - 51], ja[i - 50],
                      ja[i - 49],
                      ja[i - 48], ja[i - 47], ja[i - 46], ja[i - 45], ja[i - 44], ja[i - 43], ja[i - 42], ja[i - 41],
                      ja[i - 40],
                      ja[i - 39], ja[i - 38], ja[i - 37], ja[i - 36], ja[i - 35], ja[i - 34], ja[i - 33], ja[i - 32],
                      ja[i - 31],
                      ja[i - 30], ja[i - 29], ja[i - 28], ja[i - 27], ja[i - 26], ja[i - 25], ja[i - 24], ja[i - 23],
                      ja[i - 22],
                      ja[i - 21], ja[i - 20], ja[i - 19], ja[i - 18], ja[i - 17], ja[i - 16], ja[i - 15], ja[i - 14],
                      ja[i - 13],
                      ja[i - 12], ja[i - 11], ja[i - 10], ja[i - 9], ja[i - 8], ja[i - 7], ja[i - 6], ja[i - 5],
                      ja[i - 4],
                      ja[i - 3], ja[i - 2], ja[i - 1], ja[i], ja[i + 1], ja[i + 2], ja[i + 3], ja[i + 4],
                      ja[i + 5], ja[i + 6], ja[i + 7], ja[i + 8], ja[i + 9], ja[i + 10], ja[i + 11], ja[i + 12],
                      ja[i + 13], ja[i + 14], ja[i + 15], ja[i + 16], file=open("output4b.txt", "a+"), sep="")
                strix = open("output4b.txt")
                striy = strix.readlines()
                string = striy[-1]
                res = [sub.replace('_', '') for sub in string]
                x = "".join(res)
                print(x[-17:], file=open("oligob.txt", "a+"))
            else:
                print(ja[i], ja[i + 1], ja[i + 2], ja[i + 3], ja[i + 4], ja[i + 5], ja[i + 6],
                      ja[i + 7], ja[i + 8], ja[i + 9], ja[i + 10], ja[i + 11], ja[i + 12], ja[i + 13], ja[i + 14],
                      ja[i + 15], ja[i + 16], file=open("oligob.txt", "a+"), sep="")
        if i == (ln - 34):
            if ja[i] == "_" or ja[i - 1] == "_" or ja[i - 2] == "_" or \
                    ja[i - 3] == "_" or ja[i - 4] == "_" or ja[i - 5] == "_" or ja[i - 6] == "_" or ja[i - 7] == "_" or \
                    ja[i - 8] == "_" or ja[i - 9] == "_" or ja[i - 10] == "_" or ja[i - 11] == "_" or ja[i - 12] == "_" \
                    or ja[i - 13] == "_" or ja[i - 14] == "_" or ja[i - 15] == "_" or ja[i - 16] == "_":
                print(ja[i - 68], ja[i - 67], ja[i - 66],
                      ja[i - 65],
                      ja[i - 64], ja[i - 63], ja[i - 62], ja[i - 61], ja[i - 60], ja[i - 59], ja[i - 58], ja[i - 57],
                      ja[i - 56],
                      ja[i - 55], ja[i - 54], ja[i - 53], ja[i - 52], ja[i - 51], ja[i - 50], ja[i - 49], ja[i - 48],
                      ja[i - 47],
                      ja[i - 46], ja[i - 45], ja[i - 44], ja[i - 43], ja[i - 42], ja[i - 41], ja[i - 40], ja[i - 39],
                      ja[i - 38],
                      ja[i - 37], ja[i - 36], ja[i - 35], ja[i - 34], ja[i - 33], ja[i - 32], ja[i - 31], ja[i - 30],
                      ja[i - 29],
                      ja[i - 28], ja[i - 27], ja[i - 26], ja[i - 25], ja[i - 24], ja[i - 23], ja[i - 22], ja[i - 21],
                      ja[i - 20],
                      ja[i - 19], ja[i - 18], ja[i - 17], ja[i - 16], ja[i - 15], ja[i - 14], ja[i - 13], ja[i - 12],
                      ja[i - 11], ja[i - 10], ja[i - 9], ja[i - 8], ja[i - 7], ja[i - 6], ja[i - 5], ja[i - 4],
                      ja[i - 3], ja[i - 2], ja[i - 1], ja[i], file=open("output4b.txt", "a+"), sep="")
                strix = open("output4b.txt")
                striy = strix.readlines()
                string = striy[-1]
                res = [sub.replace('_', '') for sub in string]
                x = "".join(res)
                print(x[-17:], file=open("oligob.txt", "a+"))
            else:
                print(ja[i - 16], ja[i - 15], ja[i - 14], ja[i - 13], ja[i - 12], ja[i - 11], ja[i - 10],
                      ja[i - 9], ja[i - 8], ja[i - 7], ja[i - 6], ja[i - 5], ja[i - 4], ja[i - 3], ja[i - 2],
                      ja[i - 1], ja[i], file=open("oligob.txt", "a+"), sep="")
        if i >= 35 or i <= (ln - 35):
            if ja[i - 8] == "_" or ja[i - 7] == "_" or ja[i - 6] == "_" or ja[i - 5] == "_" or ja[i - 4] == "_" or \
                    ja[i - 3] == "_" or ja[i - 2] == "_" or ja[i - 1] == "_" or ja[i] == "_" or ja[i + 1] == "_" or \
                    ja[i + 2] == "_" or ja[i + 3] == "_" or ja[i + 4] == "_" or ja[i + 5] == "_" or ja[i + 6] == "_" or \
                    ja[i + 7] == "_" or ja[i + 8] == "_":
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
            elif ja[i - 8] != "_" or ja[i - 7] != "_" or ja[i - 6] != "_" or ja[i - 5] != "_" or ja[i - 4] != "_" or \
                    ja[i - 3] != "_" or ja[i - 2] != "_" or ja[i - 1] != "_" or ja[i] != "_" or ja[i + 1] != "_" or \
                    ja[i + 2] != "_" or ja[i + 3] != "_" or ja[i + 4] != "_" or ja[i + 5] != "_" or ja[i + 6] != "_" or \
                    ja[i + 7] != "_" or ja[i + 8] != "_":
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
    "Do you want to keep process files (default=False)? If no, all the files inside the working dir will be destroied exept the output file [yes/no]: ")
if answer == "yes":
    print("all the process files will be kept!")
elif answer == "no":
    print("all the process files will be discarded exept the consensus.fa file!")
    os.system("find . -type f -not -name 'consensus.fa' -delete")
else:
    print("Please enter yes or no.")a