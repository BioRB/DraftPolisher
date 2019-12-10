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
# input files format:
# >QRY
# draft sequence
# >SBJ
# reference sequence
# # the length of the sequences IDs in input.fa, has to be exactly 3 characters
# use the tool as follows:
# python DraftPolisher_cov.py -q query.fa -s subject.fa -f sequences.fa

import os
import re
from subprocess import Popen
import subprocess
from itertools import islice
import pandas as pd
import numpy as np
import argparse
from Bio.Seq import Seq


print("DraftPolisher v1.0 by Rosario Nicola Brancaccio")


# function for fragments generation
def get_slice(s, idx, n=8, ignored_chars='_'):
    if idx == 0:
        line = re.sub('[_]', '', s)
        return line[:(2*n+1)]
    if ((len(s[idx+1:]) - s[idx+1:].count('_')) >= n) and ((len(s[:idx-1]) - s[:idx-1].count('_')) >= n):
        if s[idx] in ignored_chars:
            #  adjust idx to first valid on right side
            idx = next((i for i, ch in enumerate(s[idx:], idx) if ch not in ignored_chars), None)
            if idx is None:
                return ''
        d = {i: ch for i, ch in enumerate(s) if ch not in ignored_chars}
        if idx in d:
            keys = [k for k in d.keys()]
            idx = keys.index(idx)
            return ''.join(d[k] for k in keys[max(0, idx-n):min(idx+n+1, len(s))])
    if ((len(s[idx+1:]) - s[idx+1:].count('_')) >= n) and ((len(s[:idx-1]) - s[:idx-1].count('_')) < n):
        if s[idx] in ignored_chars:
            #  adjust idx to first valid on right side:
            idx = next((i for i, ch in enumerate(s[idx:], idx) if ch not in ignored_chars), None)
            if idx is None:
                return ''
        d = {i: ch for i, ch in enumerate(s) if ch not in ignored_chars}
        if idx in d:
            keys = [k for k in d.keys()]
            idx = keys.index(idx)
            return ''.join(d[k] for k in keys[max(0, idx):min(idx+(2*n)+1, len(s))])
    if ((len(s[idx+1:]) - s[idx+1:].count('_')) < n) and ((len(s[:idx-1]) - s[:idx-1].count('_')) >= n):
        line = re.sub('[_]', '', s)
        return line[-(n+n+1):]


# run muscle
parser = argparse.ArgumentParser(description='polish draft genomes')
parser.add_argument("--query", "-q", help="query sequence", type=str)
parser.add_argument("--subject", "-s", help="reference sequence", type=str)
parser.add_argument("--sequences_database", "-f", help="sequencese database (any sequences file in FASTA format)", type=str)
parser.add_argument("--kmer_size", "-k", help="k-mer size", type=int)
args = parser.parse_args()
os.system("cat {0} {1} > inseq.fa".format(args.input1, args.input2))
myfile = open('inseq1.fa', 'w')
with open("inseq.fa") as inseq:
    lines = inseq.readlines()
    lines[0] = ">QRY\n"
    lines[2] = ">SBY\n"
    myfile.writelines(lines)
    myfile.close()
os.system("muscle -in inseq1.fa -clw > muscle_out.fa")
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
# checkpoint for start or end with missing base in one off the sequences
    for linea in lines:
        sima = "_"
        if linea[0] in sima:
            print("After the MUSCLE alignment a gap was identified in the BEGINNING of the alignment.")
            print("This is incompatible with this tool,"
                  " please remove this portion of the sequence and run again the tool.")
            exit()
        if linea[-1] in sima:
            print("After the MUSCLE alignment a gap was identified in the END of the alignment.")
            print("This is incompatible with this tool,"
                  " please remove this portion of the sequence and run again the tool.")
            exit()

# processing of line 0 that correspond to the Query sequence, the draft genome to correct (query)
    ja = lines[0].rstrip("\n\r")
    ln = len(hh)
    for i in pos:
        print(get_slice(ja, i, args.input4, '_'), file=open("oligo.txt", "a+"))
# generate the reverse and complement of the oligos and save fw and rv in the same file

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
    fields = None
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
        print(get_slice(ja, i, args.input4, '_'), file=open("oligob.txt", "a+"))
# generate the reverse and complement of the oligos and save fw and rv in the same file
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
    fields = None
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
    "Do you want to keep process files (default=False)?\n"
    "If no, all the files inside the working dir will be destroied exept the output file [yes/no]: ")
if answer == "yes":
    print("all the process files will be kept!")
elif answer == "no":
    print("all the process files will be discarded exept the consensus.fa file!")
    os.system("find . -type f -not -name 'consensus.fa' -delete")
else:
    print("Please enter yes or no.")
