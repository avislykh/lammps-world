# Copyright 2023 Alexander Vislykh
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import math
import pip
package_names=['biopython==1.68', "phylopandas", "pyasr", "cogent3"] #packages to install
pip.main(['install'] + package_names)
import phylopandas as pd
import dendropy as d
import pyasr
import cogent3
from cogent3 import make_aligned_seqs, make_unaligned_seqs
from cogent3.evolve.models import get_model
from cogent3 import load_unaligned_seqs, make_tree
from cogent3.align.progressive import TreeAlign
from cogent3 import make_unaligned_seqs
from cogent3 import load_aligned_seqs, make_tree
from cogent3.evolve.models import get_model
import random
import pip
import os
os.system("pip install rdkit")
import rdkit
from rdkit import Chem
from os import listdir
dname = 1
era = 1
moltext = ""
txtid = ""
naname = ""
m = ""
finalform = ""
worldfile_list = [""]
atom_list = [[1, '2', '0', '0', '0', '0', '4.003', '0', '0', '0']]

os.system("cd /workspace/sudoku")




def part_str():
    global worldfile, worldfile_list
    wpartl = worldfile.splitlines()
    j_nl = """
"""
    if len(wpartl) > 1000000:
        worldfile_list = [j_nl.join(wpartl[:1000000]), worldfile_list]
        worldfile = j_nl.join(wpartl[1000000:])

def part_list():
    global list_j, atom_list
    j_nl = """
"""
    if len(list_j) > 1000000:
        atom_list = [list_j[:1000000], atom_list]
        list_j = list_j[1000000:]





def syndata():
    global moltext, txtid
    typemap = """H
He
Li
Be
B
C
N
O
F
Ne
Na
Mg
Al
Si
P
S
Cl
Ar
K
Ca
Sc
Ti
V
Cr
Mn
Fe
Co
Ni
Cu
Zn
Ga
Ge
As
Se
Br
Kr
Rb
Sr
Y
Zr
Nb
Mo
Tc
Ru
Rh
Pd
Ag
Cd
In
Sn
Sb
Te
I
Xe
Cs
Ba
La
Ce
Pr
Nd
Pm
Sm
Eu
Gd
Tb
Dy
Ho
Er
Tm
Yb
Lu
Hf
Ta
W
Re
Os
Ir
Pt
Au
Hg
Tl
Pb
Bi
Po
At
Rn
Fr
Ra
Ac
Th
Pa
U
Np
Pu
Am
Cm
x
xx
X
Xx"""
    masstable = """mass 1 1.008
mass 2 4.003
mass 3 6.941
mass 4 9.012
mass 5 10.811
mass 6 12.011
mass 7 14.007
mass 8 15.999
mass 9 18.998
mass 10 20.180
mass 11 22.990
mass 12 24.305
mass 13 26.982
mass 14 28.086
mass 15 30.974
mass 16 32.065
mass 17 35.453
mass 18 39.948
mass 19 39.098
mass 20 40.078
mass 21 44.956
mass 22 47.867
mass 23 50.942
mass 24 51.996
mass 25 54.938
mass 26 55.845
mass 27 58.933
mass 28 58.693
mass 29 63.546
mass 30 65.390
mass 31 69.723
mass 32 72.640
mass 33 74.922
mass 34 78.960
mass 35 79.904
mass 36 83.800
mass 37 85.468
mass 38 87.620
mass 39 88.906
mass 40 91.224
mass 41 92.906
mass 42 95.940
mass 43 98.000
mass 44 101.070
mass 45 102.906
mass 46 106.420
mass 47 107.868
mass 48 112.411
mass 49 114.818
mass 50 118.710
mass 51 121.760
mass 52 127.600
mass 53 126.905
mass 54 131.293
mass 55 132.906
mass 56 137.327
mass 57 138.906
mass 58 140.116
mass 59 140.908
mass 60 144.240
mass 61 145.000
mass 62 150.360
mass 63 151.964
mass 64 157.250
mass 65 158.925
mass 66 162.500
mass 67 164.930
mass 68 167.259
mass 69 168.934
mass 70 173.040
mass 71 174.967
mass 72 178.490
mass 73 180.948
mass 74 183.840
mass 75 186.207
mass 76 190.230
mass 77 192.217
mass 78 195.078
mass 79 196.967
mass 80 200.590
mass 81 204.383
mass 82 207.200
mass 83 208.980
mass 84 209.000
mass 85 210.000
mass 86 222.000
mass 87 223.000
mass 88 226.000
mass 89 227.000
mass 90 232.038
mass 91 231.036
mass 92 238.029
mass 93 237.000
mass 94 244.000
mass 95 243.000
mass 96 247.000
mass 97 247.000
mass 98 251.000
mass 99 252.000
mass 100 257.000"""
    sec1 = """
Masses
"""
    sec2 = """

Charges
"""
    sec3 = """

Coords
"""
    sec4 = """

Bonds
"""
    sec5 = """

Types
"""
    na = 1
    nb = 1
    nline = """
"""
    moltext = moltext.split(nline)
    cutline = moltext.index("M  END")
    while len(moltext) > cutline:
        moltext.pop()
    if moltext[3].endswith("V3000"):
        del moltext[6]
        del moltext[5]
        del moltext[4]
        del moltext[3]
        del moltext[2]
        del moltext[1]
        del moltext[0]
        for u in moltext:
            elem = u.split(" ")
            for y in elem:
                if y == '':
                    elem.remove(y)
                if y == " ":
                    elem.remove(y)
            for y in elem:
                if y == '':
                    elem.remove(y)
                if y == " ":
                    elem.remove(y)
            if len(elem) == 8:
                tidx = typemap.splitlines()
                atype = tidx.index(elem[3]) + 1
                midx = masstable.splitlines()
                eidx = midx[atype - 1]
                ridx = eidx.split(" ")
                amass = ridx[2]
                sec1 = sec1 + nline + "  " + str(na) + " " + amass
                sec2 = sec2 + nline + "  " + str(na) + " " + "0"
                sec3 = sec3 + nline + "  " + str(na) + " " + elem[4] + " " + elem[5] + " " + elem[6]
                sec5 = sec5 + nline + "  " + str(na) + " " + str(atype)
                na = na + 1
            if len(elem) == 9:
                tidx = typemap.splitlines()
                atype = tidx.index(elem[3]) + 1
                midx = masstable.splitlines()
                eidx = midx[atype - 1]
                ridx = eidx.split(" ")
                amass = ridx[2]
                sec1 = sec1 + nline + "  " + str(na) + " " + amass
                if elem[8] == "CHG=-3":
                    sec2 = sec2 + nline + "  " + str(na) + " " + "-3"
                elif elem[8] == "CHG=-2":
                    sec2 = sec2 + nline + "  " + str(na) + " " + "-2"
                elif elem[8] == "CHG=-1":
                    sec2 = sec2 + nline + "  " + str(na) + " " + "-1"
                elif elem[8] == "CHG=1":
                    sec2 = sec2 + nline + "  " + str(na) + " " + "1"
                elif elem[8] == "CHG=2":
                    sec2 = sec2 + nline + "  " + str(na) + " " + "2"
                else:
                    sec2 = sec2 + nline + "  " + str(na) + " " + "3"
                sec3 = sec3 + nline + "  " + str(na) + " " + elem[4] + " " + elem[5] + " " + elem[6]
                sec5 = sec5 + nline + "  " + str(na) + " " + str(atype)
                na = na + 1
            if len(elem) == 6:
                sec4 = sec4 + nline + "  " + str(nb) + " " + elem[3] + " " + elem[4] + " " + elem[5]
                nb = nb + 1
    else:
        for u in moltext:
            elem = u.split(" ")
            for y in elem:
                if y == '':
                    elem.remove(y)
                if y == " ":
                    elem.remove(y)
            for y in elem:
                if y == '':
                    elem.remove(y)
                if y == " ":
                    elem.remove(y)
            if len(elem) == 7:
                sec4 = sec4 + nline + "  " + str(nb) + " " + elem[2] + " " + elem[0] + " " + elem[1]
                nb = nb + 1
            if len(elem) == 16:
                tidx = typemap.splitlines()
                atype = tidx.index(elem[3]) + 1
                midx = masstable.splitlines()
                eidx = midx[atype - 1]
                ridx = eidx.split(" ")
                amass = ridx[2]
                sec1 = sec1 + nline + "  " + str(na) + " " + amass
                if elem[5] == 5:
                    sec2 = sec2 + nline + "  " + str(na) + " " + "1"
                elif elem[5] == 3:
                    sec2 = sec2 + nline + "  " + str(na) + " " + "-1"
                if elem[5] == 4:
                    sec2 = sec2 + nline + "  " + str(na) + " " + "2"
                else:
                    sec2 = sec2 + nline + "  " + str(na) + " " + elem[5]
                sec3 = sec3 + nline + "  " + str(na) + " " + elem[0] + " " + elem[1] + " " + elem[2]
                sec5 = sec5 + nline + "  " + str(na) + " " + str(atype)
                na = na + 1
    header = ""
    txtid2 = txtid
    txtid2 = txtid2[:-1]
    txtid2 = txtid2[:-1]
    txtid2 = txtid2[:-1]
    txtid2 = txtid2 + "txt"
    lmpfile = open(txtid2, "w")
    header = """LAMMPS Molecule File

"""
    header = header + str(na - 1) + " atoms" + nline
    header = header + "0" + " bonds" + nline
    header = header + """0 angles
0 dihedrals
0 impropers

"""
    header = header + sec1 + sec2 + sec3 + "" + sec5
    lmpfile.write(header)
    lmpfile.close()


def convertformat():
    global moltext, txtid, naname
    txtid = naname
    with open(naname) as mdlf:
        moltext = mdlf.read()
    na = 0
    nb = 0
    nline = """
"""
    moltext = moltext.split(nline)
    for n in moltext:
        elemcor = n.split(" ")
        for p in elemcor:
            if p == '':
                elemcor.remove(p)
            if p == " ":
                elemcor.remove(p)
        for p in elemcor:
            if p == '':
                elemcor.remove(p)
            if p == " ":
                elemcor.remove(p)
        if len(elemcor) == 8 or len(elemcor) == 9:
            elemcor[4] = str(float(elemcor[4]) * 1.5)
            elemcor[5] = str(float(elemcor[5]) * 1.5)
            moltext[moltext.index(n)] = ' '.join([str(d) for d in elemcor])
    listall = moltext
    la = listall.index("M  V30 END ATOM")
    lb = listall.index("M  V30 END BOND")
    lats = listall[6:(la)]
    lbonds = listall[la:(lb)]
    cutline = moltext.index("M  END")
    while len(moltext) > cutline:
        moltext.pop()
    if moltext[3] == " 0 0  0  0  0  0  0  0  0  0999 V3000":
        del moltext[6]
        del moltext[5]
        del moltext[4]
        del moltext[3]
        del moltext[2]
        del moltext[1]
        del moltext[0]
        for i in moltext:
            elem = i.split(" ")
            for y in elem:
                if y == '':
                    elem.remove(y)
                if y == " ":
                    elem.remove(y)
            for y in elem:
                if y == '':
                    elem.remove(y)
                if y == " ":
                    elem.remove(y)
            if len(elem) == 8 or len(elem) == 9:
                na = na + 1
            if len(elem) == 6:
                nb = nb + 1
        for u in moltext:
            elem = u.split(" ")
            val = 0
            for y in elem:
                if y == '':
                    elem.remove(y)
                if y == " ":
                    elem.remove(y)
            for y in elem:
                if y == '':
                    elem.remove(y)
                if y == " ":
                    elem.remove(y)
            if len(elem) == 8 or len(elem) == 9:
                eutype = elem[3]
                for t in moltext:
                    bondx = t.split(" ")
                    for w in bondx:
                        if w == '':
                            bondx.remove(w)
                        if w == ' ':
                            bondx.remove(w)
                    for w in bondx:
                        if w == '':
                            bondx.remove(w)
                        if w == ' ':
                            bondx.remove(w)
                    if len(bondx) == 6:
                        if bondx[4] == elem[2] or bondx[5] == elem[2]:
                            val = val + int(bondx[3])
                if elem[3] == "C":
                    if val == 1:
                        na = na + 1
                        extraline = "M  V30 " + str(na) + " H " + elem[4] + " " + elem[5] + " 1.09 0"
                        lats.append(extraline)
                        la = la + 1
                        nb = nb + 1
                        extraline = "M  V30 " + str(nb) + " 1 " + elem[2] + " " + str(na)
                        lbonds.append(extraline)
                        na = na + 1
                        extraline = "M  V30 " + str(na) + " H " + elem[4] + " " + elem[5] + " -1.09 0"
                        lats.append(extraline)
                        la = la + 1
                        nb = nb + 1
                        extraline = "M  V30 " + str(nb) + " 1 " + elem[2] + " " + str(na)
                        lbonds.append(extraline)
                        na = na + 1
                        extraline = "M  V30 " + str(na) + " H " + elem[4] + " " + str(float(elem[5]) - float("1.09")) + " 0 0"
                        lats.append(extraline)
                        la = la + 1
                        nb = nb + 1
                        extraline = "M  V30 " + str(nb) + " 1 " + elem[2] + " " + str(na)
                        lbonds.append(extraline)
                    if val == 2:
                        na = na + 1
                        extraline = "M  V30 " + str(na) + " H " + elem[4] + " " + elem[5] + " 1.09 0"
                        lats.append(extraline)
                        la = la + 1
                        nb = nb + 1
                        extraline = "M  V30 " + str(nb) + " 1 " + elem[2] + " " + str(na)
                        lbonds.append(extraline)
                        na = na + 1
                        extraline = "M  V30 " + str(na) + " H " + elem[4] + " " + elem[5] + " -1.09 0"
                        lats.append(extraline)
                        la = la + 1
                        nb = nb + 1
                        extraline = "M  V30 " + str(nb) + " 1 " + elem[2] + " " + str(na)
                        lbonds.append(extraline)
                    if val == 3:
                        na = na + 1
                        extraline = "M  V30 " + str(na) + " H " + elem[4] + " " + elem[5] + " 1.09 0"
                        lats.append(extraline)
                        la = la + 1
                        nb = nb + 1
                        extraline = "M  V30 " + str(nb) + " 1 " + elem[2] + " " + str(na)
                        lbonds.append(extraline)
                if elem[3] == "N":
                    if val == 1:
                        na = na + 1
                        extraline = "M  V30 " + str(na) + " H " + elem[4] + " " + elem[5] + " 1.09 0"
                        lats.append(extraline)
                        la = la + 1
                        nb = nb + 1
                        extraline = "M  V30 " + str(nb) + " 1 " + elem[2] + " " + str(na)
                        lbonds.append(extraline)
                        na = na + 1
                        extraline = "M  V30 " + str(na) + " H " + elem[4] + " " + elem[5] + " -1.09 0"
                        lats.append(extraline)
                        la = la + 1
                        nb = nb + 1
                        extraline = "M  V30 " + str(nb) + " 1 " + elem[2] + " " + str(na)
                        lbonds.append(extraline)
                    if val == 2:
                        na = na + 1
                        extraline = "M  V30 " + str(na) + " H " + elem[4] + " " + elem[5] + " 1.09 0"
                        lats.append(extraline)
                        la = la + 1
                        nb = nb + 1
                        extraline = "M  V30 " + str(nb) + " 1 " + elem[2] + " " + str(na)
                        lbonds.append(extraline)
    lvit1 = ["x", "x", "x", " 0 0  0  0  0  0  0  0  0  0999 V3000", "M  V30 BEGIN CTAB", "M  V30 COUNTS " + str(na) +" " + str(nb) + " 0 0 1"]                    
    lvit2 = ["M  V30 END BOND", "M  V30 END CTAB", "M  END"]
    listall = lvit1 + lats + lbonds + lvit2
    instext = ""
    for y in listall:
        instext = instext + y + nline
    tarfile = open(txtid, "w")
    tarfile.write(instext)
    moltext = instext
    syndata()
    tarfile.close()






path = "/workspace/sudoku/mols"
dir_list = os.listdir(path)
for i in dir_list:
    with open(path + "/" + i) as mdlf:
        moltext = mdlf.read()
    txtid = i
    syndata()




import bs4
import requests

url = "https://ftp.ncbi.nlm.nih.gov/genomes/"
r = requests.get(url)
data = bs4.BeautifulSoup(r.text, "html.parser")
urllist = []
for l in data.find_all("a"):
    r = requests.get(url + l["href"])
    hyperlink = url + l["href"]
    if r.status_code == 200:
        urllist.append(hyperlink)
for count in range(10):
    listlength = len(urllist)
    listloc = 0
    for count in range(int(listlength)):
        url = urllist[listloc]
        r = requests.get(url)
        data = bs4.BeautifulSoup(r.text, "html.parser")
        for l in data.find_all("a"):
            r = requests.get(url + l["href"])
            hyperlink = url + l["href"]
            if r.status_code == 200:
                urllist.append(hyperlink)
        listloc = listloc + 1
rawgenomes = []
listlength = len(urllist)
listloc = 0
for count in range(int(listlength)):
    url = urllist[listloc]
    if url.endswith('.fna.gz'):
        r = requests.get(url)
        data = r.content
        rawgenomes.append(data)
    listloc = listloc + 1
import gzip
fnagenomes = []
listlength = len(rawgenomes)
listloc = 0
for count in range(int(listlength)):
    gztxt = rawgenomes[listloc]
    fnadata = gzip.decompress(gztxt)
    fnagenomes.append(fnadata)
    listloc = listloc + 1
genpure = []
genomenum = 1
genamount = len(fnagenomes)
for count2 in range(int(genamount)):
  fnabyte = fnagenomes[int(genomenum - 1)]
  fnatext = str(fnabyte, "UTF-8")
  lineloc = 1
  puregentext = ''
  textlength = len(fnatext)
  for count in range(int(textlength)):
    if fnatext[int(lineloc - 1)] == 'a' or fnatext[int(lineloc - 1)] == 'A':
      puregentext = str(puregentext) + 'A'
    elif fnatext[int(lineloc - 1)] == 'g' or fnatext[int(lineloc - 1)] == 'G':
      puregentext = str(puregentext) + 'G'
    elif fnatext[int(lineloc - 1)] == 'c' or fnatext[int(lineloc - 1)] == 'C':
      puregentext = str(puregentext) + 'C'
    elif fnatext[int(lineloc - 1)] == 't' or fnatext[int(lineloc - 1)] == 'T' or fnatext[int(lineloc - 1)] == 'u' or fnatext[int(lineloc - 1)] == 'U':
      puregentext = str(puregentext) + 'T'
    lineloc = lineloc + 1
  genpure.append(puregentext)
  genomenum = genomenum + 1



url = "https://ftp.ncbi.nlm.nih.gov/genomes/genbank/protozoa/Physarum_polycephalum/all_assembly_versions/GCA_000413255.3_Physarum_polycephalum-10.0/GCA_000413255.3_Physarum_polycephalum-10.0_genomic.fna.gz"
r = requests.get(url)
data = r.content
fnabyte = gzip.decompress(data)
fnatext = str(fnabyte, "UTF-8")
lineloc = 1
puregentext = ''
textlength = len(fnatext)
for count in range(int(textlength)):
   if fnatext[int(lineloc - 1)] == 'a' or fnatext[int(lineloc - 1)] == 'A':
    puregentext = str(puregentext) + 'A'
   elif fnatext[int(lineloc - 1)] == 'g' or fnatext[int(lineloc - 1)] == 'G':
     puregentext = str(puregentext) + 'G'
   elif fnatext[int(lineloc - 1)] == 'c' or fnatext[int(lineloc - 1)] == 'C':
     puregentext = str(puregentext) + 'C'
   elif fnatext[int(lineloc - 1)] == 't' or fnatext[int(lineloc - 1)] == 'T' or fnatext[int(lineloc - 1)] == 'u' or fnatext[int(lineloc - 1)] == 'U':
     puregentext = str(puregentext) + 'T'
   lineloc = lineloc + 1

slime = puregentext



fnatext = """1 tgacgtccaa gggaagaacc gggtacccgg aatagcacct agcggtctaa gacgcaccgc
       61 gagaaactac accgagctgg gtctctgaag taagattcag aacctagcaa tctcaagagc
      121 gctgcactct caaaacggct attaggctct ggacccccca tccagagacc ccaggggcac
      181 tcgggctgca ggccgagtgg taaaagctgt tattaacgtg gtccaagctc cacgagtgca
      241 ggttttacgg tctagcttct gcaacgaagc tggatccatc gacggcgtca cgtcgatgaa
      301 cccctcaaga cgccttgggg aaggcgtccc gggacccccg tcctgcggtt tccccagccc
      361 ccgctgtcgc cagcaacgac agtccccgcc atcccctcac ctgaggcgtg atggtttaaa
      421 ctgggccagg atcgtcaagc tcacgcgggt"""


lineloc = 1
puregentext = ''
textlength = len(fnatext)
for count in range(int(textlength)):
    if fnatext[int(lineloc - 1)] == 'a' or fnatext[int(lineloc - 1)] == 'A':
      puregentext = str(puregentext) + 'A'
    elif fnatext[int(lineloc - 1)] == 'g' or fnatext[int(lineloc - 1)] == 'G':
      puregentext = str(puregentext) + 'G'
    elif fnatext[int(lineloc - 1)] == 'c' or fnatext[int(lineloc - 1)] == 'C':
      puregentext = str(puregentext) + 'C'
    elif fnatext[int(lineloc - 1)] == 't' or fnatext[int(lineloc - 1)] == 'T' or fnatext[int(lineloc - 1)] == 'u' or fnatext[int(lineloc - 1)] == 'U':
      puregentext = str(puregentext) + 'T'
    lineloc = lineloc + 1
genpure.append(puregentext)


fnatext = """1 taccctctcc ttacccccct ccagaagcta ctaccatagc gcagcaactg aagagttggc
       61 aggaaaggat caaagatcct aaagtcacaa tctgttgttt gcccggttgc agctccgggg
      121 gttcggtcac ttgtctcgtc gggttttcac cccaggattt gttgtcctga gcgtggagcg
      181 gccgttatac aggaaacaaa gcctgtaccg gtttctcacc ccaagtccac aagtggccca
      241 cctcacggcc ctggaggctc ccttacggga gcacaggttt caccgctgcc agcacagtct
      301 atttataatt aatcgccact gtgcactgcg atgggaccag gacagggttc cccaaccctg
      361 ttttttttca caattcttcc ttggcaggaa agaacggcgt cctaatgccc aggggcttta
      421 ggagagctt"""


lineloc = 1
puregentext = ''
textlength = len(fnatext)
for count in range(int(textlength)):
    if fnatext[int(lineloc - 1)] == 'a' or fnatext[int(lineloc - 1)] == 'A':
      puregentext = str(puregentext) + 'A'
    elif fnatext[int(lineloc - 1)] == 'g' or fnatext[int(lineloc - 1)] == 'G':
      puregentext = str(puregentext) + 'G'
    elif fnatext[int(lineloc - 1)] == 'c' or fnatext[int(lineloc - 1)] == 'C':
      puregentext = str(puregentext) + 'C'
    elif fnatext[int(lineloc - 1)] == 't' or fnatext[int(lineloc - 1)] == 'T' or fnatext[int(lineloc - 1)] == 'u' or fnatext[int(lineloc - 1)] == 'U':
      puregentext = str(puregentext) + 'T'
    lineloc = lineloc + 1
genpure.append(puregentext)

fnatext = """1 atgccccccc cctccagaag ctactaccat agcgcagcaa ctgaagagtt ggcaggaaag
       61 gatcaaagat cctaaagtca caatctgttg tttgcccggt ttgcagctcc gggggttcgg
      121 tcacttgtct cgtcgggttt tcaccccagg atttgttgtc ctgagcgtgg agcggccgtt
      181 atacaggaaa caaagcctgt accggtttct caccccaagt ccacaagtgg cccacctcac
      241 ggccctggag gctcccttac gggagcacag gtttcaccgc tgccagcaca gtctatttat
      301 aattaatcgc cactgtgcac tgcgatggga ccaggacagg gttccccaac cctgtttttt
      361 cacaattctt ccttggcagg aaagaacggc gtccta"""


lineloc = 1
puregentext = ''
textlength = len(fnatext)
for count in range(int(textlength)):
    if fnatext[int(lineloc - 1)] == 'a' or fnatext[int(lineloc - 1)] == 'A':
      puregentext = str(puregentext) + 'A'
    elif fnatext[int(lineloc - 1)] == 'g' or fnatext[int(lineloc - 1)] == 'G':
      puregentext = str(puregentext) + 'G'
    elif fnatext[int(lineloc - 1)] == 'c' or fnatext[int(lineloc - 1)] == 'C':
      puregentext = str(puregentext) + 'C'
    elif fnatext[int(lineloc - 1)] == 't' or fnatext[int(lineloc - 1)] == 'T' or fnatext[int(lineloc - 1)] == 'u' or fnatext[int(lineloc - 1)] == 'U':
      puregentext = str(puregentext) + 'T'
    lineloc = lineloc + 1
genpure.append(puregentext)

fnatext = """1 tcacgggagc acaggtttca ccgctgccag cacagtctat ttataattaa tcgccactgt
       61 gcactgcgat gggaccagga cagggtcccc caaccctgtg ttttcacaat tcttccttgg
      121 caggaaagaa cggcgtccta atgccccccc cctccagaag ctactaccat agcgcagcaa
      181 ctgaagagtt ggcaggaaag gatcaaagat cctaaagtca caatctgttg tttgcccggt
      241 tgcagctccg ggggttcggt cacttgtctc gtcgggtttt cgccccagga tttgttgtcc
      301 tgagcgtgga gcggccgtta tacaggaaac aaagcctgta ccggtttctc accccaagtc
      361 cacaagtggc ccaa"""


lineloc = 1
puregentext = ''
textlength = len(fnatext)
for count in range(int(textlength)):
    if fnatext[int(lineloc - 1)] == 'a' or fnatext[int(lineloc - 1)] == 'A':
      puregentext = str(puregentext) + 'A'
    elif fnatext[int(lineloc - 1)] == 'g' or fnatext[int(lineloc - 1)] == 'G':
      puregentext = str(puregentext) + 'G'
    elif fnatext[int(lineloc - 1)] == 'c' or fnatext[int(lineloc - 1)] == 'C':
      puregentext = str(puregentext) + 'C'
    elif fnatext[int(lineloc - 1)] == 't' or fnatext[int(lineloc - 1)] == 'T' or fnatext[int(lineloc - 1)] == 'u' or fnatext[int(lineloc - 1)] == 'U':
      puregentext = str(puregentext) + 'T'
    lineloc = lineloc + 1
genpure.append(puregentext)

fnatext = """1 aacattagtg ttcgctggtg aaatctatac caatcctcta agtggctcga tctcccactc
       61 taaatatcga agaatttaga atgggagggt gactacctgt accttttggg tcaggtaaca
      121 ccccatctac acaagtcaac tatttggttg aattatttta cgggcgtacg tccatgaaat
      181 gattaagtta catatagact ggcccaagta ggtggacaac a"""


lineloc = 1
puregentext = ''
textlength = len(fnatext)
for count in range(int(textlength)):
    if fnatext[int(lineloc - 1)] == 'a' or fnatext[int(lineloc - 1)] == 'A':
      puregentext = str(puregentext) + 'A'
    elif fnatext[int(lineloc - 1)] == 'g' or fnatext[int(lineloc - 1)] == 'G':
      puregentext = str(puregentext) + 'G'
    elif fnatext[int(lineloc - 1)] == 'c' or fnatext[int(lineloc - 1)] == 'C':
      puregentext = str(puregentext) + 'C'
    elif fnatext[int(lineloc - 1)] == 't' or fnatext[int(lineloc - 1)] == 'T' or fnatext[int(lineloc - 1)] == 'u' or fnatext[int(lineloc - 1)] == 'U':
      puregentext = str(puregentext) + 'T'
    lineloc = lineloc + 1
genpure.append(puregentext)

fnatext = """1 tctcccactc taaatatcga agaatttaga ataggagggt gactacctgt accttttggg
       61 tcaggtaaca ccccatctac acaagtcaac tatttggttg aattatttta cgggcgtacg
      121 tccatgaaat gattaagtta catatagact ggctcga"""


lineloc = 1
puregentext = ''
textlength = len(fnatext)
for count in range(int(textlength)):
    if fnatext[int(lineloc - 1)] == 'a' or fnatext[int(lineloc - 1)] == 'A':
      puregentext = str(puregentext) + 'A'
    elif fnatext[int(lineloc - 1)] == 'g' or fnatext[int(lineloc - 1)] == 'G':
      puregentext = str(puregentext) + 'G'
    elif fnatext[int(lineloc - 1)] == 'c' or fnatext[int(lineloc - 1)] == 'C':
      puregentext = str(puregentext) + 'C'
    elif fnatext[int(lineloc - 1)] == 't' or fnatext[int(lineloc - 1)] == 'T' or fnatext[int(lineloc - 1)] == 'u' or fnatext[int(lineloc - 1)] == 'U':
      puregentext = str(puregentext) + 'T'
    lineloc = lineloc + 1
genpure.append(puregentext)










worldfile = """#LAMMPS Input file for world generation
# Lammps World
# Intialization
units real
dimension 3
boundary f f f
thermo_modify lost ignore
neigh_modify one 100000 page 10000000
atom_style full
read_data set.txt
pair_style hybrid reaxff NULL enobonds yes lj/cut/coul/cut 10.0
pair_coeff * * lj/cut/coul/cut 100.0 1.5 9.0 9.0
pair_coeff * * reaxff ffield.reax.FC H NULL NULL NULL NULL C N O F NULL NULL NULL NULL C N S Cl NULL NULL NULL Ni Ni Ni Ni Ni Ni Ni Ni Pt Ni NULL C N S Cl NULL NULL NULL Ni Ni Ni Ni Ni Ni Ni Pt Pt Ni NULL C N S Cl NULL NULL NULL Ni Ni Ni Ni Ni Ni Ni Ni Ni Ni Ni Ni Ni Ni Ni Ni Ni Ni Ni Ni Ni Pt Pt Ni NULL C N S Cl NULL NULL NULL NULL NULL NULL NULL NULL NULL NULL NULL NULL NULL NULL NULL
molecule 1 earth1.txt
molecule 2 earth2.txt
molecule 3 earth3.txt
molecule 4 earth4.txt
molecule 5 earth5.txt
molecule 6 earth6.txt
molecule 7 earth7.txt
molecule 8 rare1.txt
molecule 9 rare2.txt
molecule 11 water.txt
molecule 12 co2.txt
molecule 13 n2.txt
molecule 14 oxygen.txt
molecule 15 helium.txt
molecule 16 neon.txt
molecule 17 argon.txt
molecule 18 krypton.txt
molecule 19 xenon.txt
molecule 21 aluminium.txt
molecule 22 copper.txt
molecule 23 gold.txt
molecule 24 nickel.txt
molecule 25 palladium.txt
molecule 26 platinum.txt
molecule 27 silver.txt
molecule 30 ascacid.txt
molecule 31 biotin.txt
molecule 32 cobalamin.txt
molecule 33 ergocalciferol.txt
molecule 34 folacin.txt
molecule 35 heme.txt
molecule 36 moco.txt
molecule 37 niacin.txt
molecule 38 phylloquinone.txt
molecule 39 pntacid.txt
molecule 40 pyridoxine.txt
molecule 41 retinal.txt
molecule 42 riboflavin.txt
molecule 43 thiamine.txt
molecule 44 tocopherol.txt
molecule 50 calcpho.txt
molecule 51 cr2o3.txt
molecule 52 cu2.txt
molecule 53 iodide.txt
molecule 54 mg2.txt
molecule 55 mn2.txt
molecule 56 molybdate.txt
molecule 57 nacl.txt
molecule 58 nak.txt
molecule 59 zn.txt
molecule 60 adp.txt
molecule 61 amp.txt
molecule 62 atp.txt
molecule 63 ctp.txt
molecule 64 fad.txt
molecule 65 fadh2.txt
molecule 66 glucose.txt
molecule 67 gtp.txt
molecule 68 nad.txt
molecule 69 nadh.txt
molecule 70 utp.txt
molecule 100 formylmethionine.txt
molecule 101 selenocysteine.txt
molecule 102 ala.txt
molecule 103 arg.txt
molecule 104 asp.txt
molecule 105 asr.txt
molecule 106 cys.txt
molecule 107 glu.txt
molecule 108 gly.txt
molecule 109 gte.txt
molecule 110 his.txt
molecule 111 iso.txt
molecule 112 leu.txt
molecule 113 lys.txt
molecule 114 met.txt
molecule 115 phe.txt
molecule 116 pro.txt
molecule 117 ser.txt
molecule 118 tre.txt
molecule 119 try.txt
molecule 120 tyr.txt
molecule 121 val.txt
molecule 200 glyc1.txt
molecule 201 glyc2.txt
molecule 202 glyc3.txt
molecule 203 glyc4.txt
molecule 204 glyc5.txt
molecule 205 glyc6.txt
molecule 206 glyc7.txt
molecule 207 glyc8.txt
molecule 208 glyc9.txt
molecule 209 glyc10.txt
molecule 210 formyl.txt
molecule 211 metsyn.txt
molecule 212 synfory.txt
molecule 213 synthr.txt
molecule 214 fors.txt
molecule 215 sulfox.txt
molecule 216 alars.txt
molecule 217 argrs.txt
molecule 218 asnrs.txt
molecule 219 asprs.txt
molecule 220 cysrs.txt
molecule 221 glnrs.txt
molecule 222 glurs.txt
molecule 223 glyrs.txt
molecule 224 hisrs.txt
molecule 225 ilers.txt
molecule 226 leurs.txt
molecule 227 lysrs.txt
molecule 228 metrs.txt
molecule 229 phers1.txt
molecule 230 prors.txt
molecule 231 serrs.txt
molecule 232 thrrs.txt
molecule 233 trprs.txt
molecule 234 tyrrs.txt
molecule 235 valrs.txt
molecule 236 phers2.txt
molecule 237 r30s.txt
molecule 238 deox.txt
molecule 239 if1.txt
molecule 240 if2.txt
molecule 241 if3.txt
molecule 242 pol1.txt
molecule 243 pol2.txt
molecule 244 pol3.txt
molecule 245 pol4.txt
molecule 246 pol5.txt
molecule 247 ribos1.txt
molecule 248 rnapol.txt
molecule 249 r50s1.txt
molecule 250 r50s2.txt
molecule 251 eftu.txt
molecule 252 factor1.txt
molecule 253 factor2.txt
molecule 254 factor3.txt
molecule 255 factor4.txt
molecule 256 factor5.txt
molecule 257 factor6.txt
molecule 258 factor7.txt
molecule 259 factor8.txt
molecule 260 factor9.txt
molecule 261 factor9.txt
molecule 262 frr.txt
molecule 263 rl1.txt
molecule 264 rl10.txt
molecule 265 rl11.txt
molecule 266 rl16.txt
molecule 267 rl31.txt
molecule 268 rl36.txt
molecule 269 rl712.txt
region bedrock block 100 100000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000 100 10000000000000 100 100000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
lattice fcc 3.524
create_atoms 0 region bedrock mol 24 3
region bar block 100 100000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000 10000000000100 120000000000000 50000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000 50000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000110000000000000
lattice fcc 3.524
create_atoms 0 region bar mol 24 3
region ocean block 100 100000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000 10000000000100 120000000000000 50000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000110000000000100 100000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
lattice sc 3.1061591792
create_atoms 0 region ocean mol 11 3
region nit block 100 100000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000 120000000000100 186300000000000 100 100000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
lattice sc 34.181298585
create_atoms 0 region nit mol 13 3
region oxy block 100 100000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000 186300000000100 204300000000000 100 100000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
lattice sc 34.181298585
create_atoms 0 region oxy mol 14 3
region wat block 100 100000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000 204300000000100 204510000000000 100 100000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
lattice sc 34.181298585
create_atoms 0 region wat mol 11 3
region co block 100 100000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000 204510000000100 204544000000000 100 100000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
lattice sc 34.181298585
create_atoms 0 region co mol 12 3
region he block 100 100000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000 204544000000100 204545000000000 100 100000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
lattice sc 34.181298585
create_atoms 0 region he mol 15 3
region ne block 100 100000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000 204545000000100 204546000000000 100 100000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
lattice sc 34.181298585
create_atoms 0 region ne mol 16 3
region ar block 100 100000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000 204546000000100 204547000000000 100 100000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
lattice sc 34.181298585
create_atoms 0 region ar mol 17 3
region kr block 100 100000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000 204547000000100 204548000000000 100 100000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
lattice sc 34.181298585
create_atoms 0 region kr mol 18 3
region xe block 100 100000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000 204548000000100 204549000000000 100 100000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
lattice sc 34.181298585
create_atoms 0 region xe mol 19 3
region exo block 100 10000000000000 1000000000000000 100000000000000000 100 10000000000000
lattice sc 100000000000
create_atoms 0 region exo mol 19 3
"""




def makestrand():
    global finalform, dname, avatarx, avatary, avatarz, worldfile, m, naname, txtid
    lnovo = """
"""
    atomblock = """
M  V30 BEGIN ATOM"""
    bondblock = """M  V30 BEGIN BOND"""
    dnamdlfile = """x
x
x"""
    merj1 = 0
    merj2 = 0
    bpairmono = 0
    allatom = 0
    allbond = 0
    atommer = 0
    for count in range(int(len(finalform))):
      monomer = finalform[int(bpairmono)]
      if monomer == 'A':
        atomblock = ''.join([str(x) for x in [atomblock, lnovo, 'M  V30 ', allatom + 1, ' ', 'C', ' ', 0.9578 + 6 * bpairmono, ' ', -1.9677, ' ', 0, ' ', 0]])
        allatom = allatom + 1
        atomblock = ''.join([str(x) for x in [atomblock, lnovo, 'M  V30 ', allatom + 1, ' ', 'C', ' ', 1.8238 + 6 * bpairmono, ' ', -2.4677, ' ', 0, ' ', 0]])
        allatom = allatom + 1
        atomblock = ''.join([str(x) for x in [atomblock, lnovo, 'M  V30 ', allatom + 1, ' ', 'C', ' ', 1.8238 + 6 * bpairmono, ' ', -3.4677, ' ', 0, ' ', 0]])
        allatom = allatom + 1
        atomblock = ''.join([str(x) for x in [atomblock, lnovo, 'M  V30 ', allatom + 1, ' ', 'N', ' ', 0.9578 + 6 * bpairmono, ' ', -3.9677, ' ', 0, ' ', 0]])
        allatom = allatom + 1
        atomblock = ''.join([str(x) for x in [atomblock, lnovo, 'M  V30 ', allatom + 1, ' ', 'C', ' ', 0.0918 + 6 * bpairmono, ' ', -3.4677, ' ', 0, ' ', 0]])
        allatom = allatom + 1
        atomblock = ''.join([str(x) for x in [atomblock, lnovo, 'M  V30 ', allatom + 1, ' ', 'N', ' ', 0.0918 + 6 * bpairmono, ' ', -2.4677, ' ', 0, ' ', 0]])
        allatom = allatom + 1
        atomblock = ''.join([str(x) for x in [atomblock, lnovo, 'M  V30 ', allatom + 1, ' ', 'O', ' ', 0.9578 + 6 * bpairmono, ' ', -0.9677, ' ', 0, ' ', 0]])
        allatom = allatom + 1
        atomblock = ''.join([str(x) for x in [atomblock, lnovo, 'M  V30 ', allatom + 1, ' ', 'O', ' ', -0.7743 + 6 * bpairmono, ' ', -3.9677, ' ', 0, ' ', 0]])
        allatom = allatom + 1
        atomblock = ''.join([str(x) for x in [atomblock, lnovo, 'M  V30 ', allatom + 1, ' ', 'C', ' ', 2.6898 + 6 * bpairmono, ' ', -1.9677, ' ', 0, ' ', 0]])
        allatom = allatom + 1
        atomblock = ''.join([str(x) for x in [atomblock, lnovo, 'M  V30 ', allatom + 1, ' ', 'C', ' ', 0.9578 + 6 * bpairmono, ' ', -4.9677, ' ', 0, ' ', 0]])
        allatom = allatom + 1
        atomblock = ''.join([str(x) for x in [atomblock, lnovo, 'M  V30 ', allatom + 1, ' ', 'C', ' ', 1.7668 + 6 * bpairmono, ' ', -5.5555, ' ', 0, ' ', 0]])
        allatom = allatom + 1
        atomblock = ''.join([str(x) for x in [atomblock, lnovo, 'M  V30 ', allatom + 1, ' ', 'C', ' ', 1.4578 + 6 * bpairmono, ' ', -6.5066, ' ', 0, ' ', 0]])
        allatom = allatom + 1
        atomblock = ''.join([str(x) for x in [atomblock, lnovo, 'M  V30 ', allatom + 1, ' ', 'C', ' ', 0.4578 + 6 * bpairmono, ' ', -6.5066, ' ', 0, ' ', 0]])
        allatom = allatom + 1
        atomblock = ''.join([str(x) for x in [atomblock, lnovo, 'M  V30 ', allatom + 1, ' ', 'O', ' ', 0.1488 + 6 * bpairmono, ' ', -5.5555, ' ', 0, ' ', 0]])
        allatom = allatom + 1
        atomblock = ''.join([str(x) for x in [atomblock, lnovo, 'M  V30 ', allatom + 1, ' ', 'O', ' ', 2.4578 + 6 * bpairmono, ' ', -6.5066, ' ', 0, ' ', 0]])
        allatom = allatom + 1
        atomblock = ''.join([str(x) for x in [atomblock, lnovo, 'M  V30 ', allatom + 1, ' ', 'C', ' ', -0.5422 + 6 * bpairmono, ' ', -6.5066, ' ', 0, ' ', 0]])
        allatom = allatom + 1
        atomblock = ''.join([str(x) for x in [atomblock, lnovo, 'M  V30 ', allatom + 1, ' ', 'O', ' ', -1.5422 + 6 * bpairmono, ' ', -6.5066, ' ', 0, ' ', 0]])
        allatom = allatom + 1
        atomblock = ''.join([str(x) for x in [atomblock, lnovo, 'M  V30 ', allatom + 1, ' ', 'P', ' ', -2.5422 + 6 * bpairmono, ' ', -6.5066, ' ', 0, ' ', 0]])
        allatom = allatom + 1
        if bpairmono > 0:
          bondblock = ''.join([str(x) for x in [bondblock, lnovo, 'M  V30 ', 1 + allbond, ' ', 1, ' ', merj1, ' ', allatom]])
          allbond = allbond + 1
        merj1 = allatom - 3
        atomblock = ''.join([str(x) for x in [atomblock, lnovo, 'M  V30 ', allatom + 1, ' ', 'C', ' ', 1.5422 + 6 * bpairmono, ' ', 4.9677, ' ', 0, ' ', 0]])
        allatom = allatom + 1
        atomblock = ''.join([str(x) for x in [atomblock, lnovo, 'M  V30 ', allatom + 1, ' ', 'C', ' ', 0.7332 + 6 * bpairmono, ' ', 5.5555, ' ', 0, ' ', 0]])
        allatom = allatom + 1
        atomblock = ''.join([str(x) for x in [atomblock, lnovo, 'M  V30 ', allatom + 1, ' ', 'C', ' ', 1.0422 + 6 * bpairmono, ' ', 6.5066, ' ', 0, ' ', 0]])
        allatom = allatom + 1
        atomblock = ''.join([str(x) for x in [atomblock, lnovo, 'M  V30 ', allatom + 1, ' ', 'C', ' ', 2.0422 + 6 * bpairmono, ' ', 6.5066, ' ', 0, ' ', 0]])
        allatom = allatom + 1
        atomblock = ''.join([str(x) for x in [atomblock, lnovo, 'M  V30 ', allatom + 1, ' ', 'O', ' ', 2.3512 + 6 * bpairmono, ' ', 5.5555, ' ', 0, ' ', 0]])
        allatom = allatom + 1
        atomblock = ''.join([str(x) for x in [atomblock, lnovo, 'M  V30 ', allatom + 1, ' ', 'C', ' ', 3.0422 + 6 * bpairmono, ' ', 6.5066, ' ', 0, ' ', 0]])
        allatom = allatom + 1
        atomblock = ''.join([str(x) for x in [atomblock, lnovo, 'M  V30 ', allatom + 1, ' ', 'O', ' ', 0.0422 + 6 * bpairmono, ' ', 6.5066, ' ', 0, ' ', 0]])
        allatom = allatom + 1
        atomblock = ''.join([str(x) for x in [atomblock, lnovo, 'M  V30 ', allatom + 1, ' ', 'P', ' ', -0.9578 + 6 * bpairmono, ' ', 6.5066, ' ', 0, ' ', 0]])
        allatom = allatom + 1
        atomblock = ''.join([str(x) for x in [atomblock, lnovo, 'M  V30 ', allatom + 1, ' ', 'O', ' ', -1.9578 + 6 * bpairmono, ' ', 6.5066, ' ', 0, ' ', 0]])
        allatom = allatom + 1
        if bpairmono > 0:
          bondblock = ''.join([str(x) for x in [bondblock, lnovo, 'M  V30 ', 1 + allbond, ' ', 1, ' ', merj2, ' ', allatom]])
          allbond = allbond + 1
        merj2 = allatom - 3
        atomblock = ''.join([str(x) for x in [atomblock, lnovo, 'M  V30 ', allatom + 1, ' ', 'O', ' ', -0.4578 + 6 * bpairmono, ' ', 7.3726, ' ', 0, ' ', 0]])
        allatom = allatom + 1
        atomblock = ''.join([str(x) for x in [atomblock, lnovo, 'M  V30 ', allatom + 1, ' ', 'O', ' ', -1.4578 + 6 * bpairmono, ' ', 7.3726, ' ', 0, ' ', 0, ' ', 'CHG=-1']])
        allatom = allatom + 1
        atomblock = ''.join([str(x) for x in [atomblock, lnovo, 'M  V30 ', allatom + 1, ' ', 'N', ' ', 1.5422 + 6 * bpairmono, ' ', 3.9677, ' ', 0, ' ', 0]])
        allatom = allatom + 1
        atomblock = ''.join([str(x) for x in [atomblock, lnovo, 'M  V30 ', allatom + 1, ' ', 'C', ' ', 2.3512 + 6 * bpairmono, ' ', 3.3799, ' ', 0, ' ', 0]])
        allatom = allatom + 1
        atomblock = ''.join([str(x) for x in [atomblock, lnovo, 'M  V30 ', allatom + 1, ' ', 'N', ' ', 2.0422 + 6 * bpairmono, ' ', 2.4289, ' ', 0, ' ', 0]])
        allatom = allatom + 1
        atomblock = ''.join([str(x) for x in [atomblock, lnovo, 'M  V30 ', allatom + 1, ' ', 'C', ' ', 1.0422 + 6 * bpairmono, ' ', 2.4289, ' ', 0, ' ', 0]])
        allatom = allatom + 1
        atomblock = ''.join([str(x) for x in [atomblock, lnovo, 'M  V30 ', allatom + 1, ' ', 'C', ' ', 0.7332 + 6 * bpairmono, ' ', 3.3799, ' ', 0, ' ', 0]])
        allatom = allatom + 1
        atomblock = ''.join([str(x) for x in [atomblock, lnovo, 'M  V30 ', allatom + 1, ' ', 'C', ' ', 0.3731 + 6 * bpairmono, ' ', 1.6857, ' ', 0, ' ', 0]])
        allatom = allatom + 1
        atomblock = ''.join([str(x) for x in [atomblock, lnovo, 'M  V30 ', allatom + 1, ' ', 'N', ' ', -0.6051 + 6 * bpairmono, ' ', 1.8936, ' ', 0, ' ', 0]])
        allatom = allatom + 1
        atomblock = ''.join([str(x) for x in [atomblock, lnovo, 'M  V30 ', allatom + 1, ' ', 'C', ' ', -0.9141 + 6 * bpairmono, ' ', 2.8447, ' ', 0, ' ', 0]])
        allatom = allatom + 1
        atomblock = ''.join([str(x) for x in [atomblock, lnovo, 'M  V30 ', allatom + 1, ' ', 'N', ' ', -0.2450 + 6 * bpairmono, ' ', 3.5878, ' ', 0, ' ', 0]])
        allatom = allatom + 1
        atomblock = ''.join([str(x) for x in [atomblock, lnovo, 'M  V30 ', allatom + 1, ' ', 'N', ' ', 0.8731 + 6 * bpairmono, ' ', 0.8197, ' ', 0, ' ', 0]])
        allatom = allatom + 1
        atomblock = ''.join([str(x) for x in [atomblock, lnovo, 'M  V30 ', allatom + 1, ' ', 'O', ' ', -2.0422 + 6 * bpairmono, ' ', -7.3726, ' ', 0, ' ', 0, ' ', 'CHG=-1']])
        allatom = allatom + 1
        atomblock = ''.join([str(x) for x in [atomblock, lnovo, 'M  V30 ', allatom + 1, ' ', 'O', ' ', -3.0422 + 6 * bpairmono, ' ', -7.3726, ' ', 0, ' ', 0]])
        allatom = allatom + 1
        bondblock = ''.join([str(x) for x in [bondblock, lnovo, 'M  V30 ', 1 + allbond, ' ', 1, ' ', 1 + atommer, ' ', 2 + atommer]])
        allbond = allbond + 1
        bondblock = ''.join([str(x) for x in [bondblock, lnovo, 'M  V30 ', 1 + allbond, ' ', 2, ' ', 1 + atommer, ' ', 7 + atommer]])
        allbond = allbond + 1
        bondblock = ''.join([str(x) for x in [bondblock, lnovo, 'M  V30 ', 1 + allbond, ' ', 2, ' ', 2 + atommer, ' ', 3 + atommer]])
        allbond = allbond + 1
        bondblock = ''.join([str(x) for x in [bondblock, lnovo, 'M  V30 ', 1 + allbond, ' ', 1, ' ', 2 + atommer, ' ', 9 + atommer]])
        allbond = allbond + 1
        bondblock = ''.join([str(x) for x in [bondblock, lnovo, 'M  V30 ', 1 + allbond, ' ', 1, ' ', 3 + atommer, ' ', 4 + atommer]])
        allbond = allbond + 1
        bondblock = ''.join([str(x) for x in [bondblock, lnovo, 'M  V30 ', 1 + allbond, ' ', 1, ' ', 4 + atommer, ' ', 5 + atommer]])
        allbond = allbond + 1
        bondblock = ''.join([str(x) for x in [bondblock, lnovo, 'M  V30 ', 1 + allbond, ' ', 1, ' ', 5 + atommer, ' ', 6 + atommer]])
        allbond = allbond + 1
        bondblock = ''.join([str(x) for x in [bondblock, lnovo, 'M  V30 ', 1 + allbond, ' ', 2, ' ', 5 + atommer, ' ', 8 + atommer]])
        allbond = allbond + 1
        bondblock = ''.join([str(x) for x in [bondblock, lnovo, 'M  V30 ', 1 + allbond, ' ', 1, ' ', 6 + atommer, ' ', 1 + atommer]])
        allbond = allbond + 1
        bondblock = ''.join([str(x) for x in [bondblock, lnovo, 'M  V30 ', 1 + allbond, ' ', 1, ' ', 10 + atommer, ' ', 4 + atommer]])
        allbond = allbond + 1
        bondblock = ''.join([str(x) for x in [bondblock, lnovo, 'M  V30 ', 1 + allbond, ' ', 1, ' ', 10 + atommer, ' ', 11 + atommer]])
        allbond = allbond + 1
        bondblock = ''.join([str(x) for x in [bondblock, lnovo, 'M  V30 ', 1 + allbond, ' ', 1, ' ', 11 + atommer, ' ', 12 + atommer]])
        allbond = allbond + 1
        bondblock = ''.join([str(x) for x in [bondblock, lnovo, 'M  V30 ', 1 + allbond, ' ', 1, ' ', 12 + atommer, ' ', 13 + atommer]])
        allbond = allbond + 1
        bondblock = ''.join([str(x) for x in [bondblock, lnovo, 'M  V30 ', 1 + allbond, ' ', 1, ' ', 12 + atommer, ' ', 15 + atommer]])
        allbond = allbond + 1
        bondblock = ''.join([str(x) for x in [bondblock, lnovo, 'M  V30 ', 1 + allbond, ' ', 1, ' ', 13 + atommer, ' ', 14 + atommer]])
        allbond = allbond + 1
        bondblock = ''.join([str(x) for x in [bondblock, lnovo, 'M  V30 ', 1 + allbond, ' ', 1, ' ', 13 + atommer, ' ', 16 + atommer]])
        allbond = allbond + 1
        bondblock = ''.join([str(x) for x in [bondblock, lnovo, 'M  V30 ', 1 + allbond, ' ', 1, ' ', 14 + atommer, ' ', 10 + atommer]])
        allbond = allbond + 1
        bondblock = ''.join([str(x) for x in [bondblock, lnovo, 'M  V30 ', 1 + allbond, ' ', 1, ' ', 16 + atommer, ' ', 17 + atommer]])
        allbond = allbond + 1
        bondblock = ''.join([str(x) for x in [bondblock, lnovo, 'M  V30 ', 1 + allbond, ' ', 1, ' ', 17 + atommer, ' ', 18 + atommer]])
        allbond = allbond + 1
        bondblock = ''.join([str(x) for x in [bondblock, lnovo, 'M  V30 ', 1 + allbond, ' ', 1, ' ', 18 + atommer, ' ', 40 + atommer]])
        allbond = allbond + 1
        bondblock = ''.join([str(x) for x in [bondblock, lnovo, 'M  V30 ', 1 + allbond, ' ', 2, ' ', 18 + atommer, ' ', 41 + atommer]])
        allbond = allbond + 1
        bondblock = ''.join([str(x) for x in [bondblock, lnovo, 'M  V30 ', 1 + allbond, ' ', 1, ' ', 19 + atommer, ' ', 20 + atommer]])
        allbond = allbond + 1
        bondblock = ''.join([str(x) for x in [bondblock, lnovo, 'M  V30 ', 1 + allbond, ' ', 1, ' ', 19 + atommer, ' ', 30 + atommer]])
        allbond = allbond + 1
        bondblock = ''.join([str(x) for x in [bondblock, lnovo, 'M  V30 ', 1 + allbond, ' ', 1, ' ', 20 + atommer, ' ', 21 + atommer]])
        allbond = allbond + 1
        bondblock = ''.join([str(x) for x in [bondblock, lnovo, 'M  V30 ', 1 + allbond, ' ', 1, ' ', 21 + atommer, ' ', 22 + atommer]])
        allbond = allbond + 1
        bondblock = ''.join([str(x) for x in [bondblock, lnovo, 'M  V30 ', 1 + allbond, ' ', 1, ' ', 21 + atommer, ' ', 25 + atommer]])
        allbond = allbond + 1
        bondblock = ''.join([str(x) for x in [bondblock, lnovo, 'M  V30 ', 1 + allbond, ' ', 1, ' ', 22 + atommer, ' ', 23 + atommer]])
        allbond = allbond + 1
        bondblock = ''.join([str(x) for x in [bondblock, lnovo, 'M  V30 ', 1 + allbond, ' ', 1, ' ', 22 + atommer, ' ', 24 + atommer]])
        allbond = allbond + 1
        bondblock = ''.join([str(x) for x in [bondblock, lnovo, 'M  V30 ', 1 + allbond, ' ', 1, ' ', 23 + atommer, ' ', 19 + atommer]])
        allbond = allbond + 1
        bondblock = ''.join([str(x) for x in [bondblock, lnovo, 'M  V30 ', 1 + allbond, ' ', 1, ' ', 25 + atommer, ' ', 26 + atommer]])
        allbond = allbond + 1
        bondblock = ''.join([str(x) for x in [bondblock, lnovo, 'M  V30 ', 1 + allbond, ' ', 1, ' ', 26 + atommer, ' ', 27 + atommer]])
        allbond = allbond + 1
        bondblock = ''.join([str(x) for x in [bondblock, lnovo, 'M  V30 ', 1 + allbond, ' ', 2, ' ', 26 + atommer, ' ', 28 + atommer]])
        allbond = allbond + 1
        bondblock = ''.join([str(x) for x in [bondblock, lnovo, 'M  V30 ', 1 + allbond, ' ', 1, ' ', 29 + atommer, ' ', 26 + atommer]])
        allbond = allbond + 1
        bondblock = ''.join([str(x) for x in [bondblock, lnovo, 'M  V30 ', 1 + allbond, ' ', 1, ' ', 30 + atommer, ' ', 31 + atommer]])
        allbond = allbond + 1
        bondblock = ''.join([str(x) for x in [bondblock, lnovo, 'M  V30 ', 1 + allbond, ' ', 2, ' ', 31 + atommer, ' ', 32 + atommer]])
        allbond = allbond + 1
        bondblock = ''.join([str(x) for x in [bondblock, lnovo, 'M  V30 ', 1 + allbond, ' ', 1, ' ', 32 + atommer, ' ', 33 + atommer]])
        allbond = allbond + 1
        bondblock = ''.join([str(x) for x in [bondblock, lnovo, 'M  V30 ', 1 + allbond, ' ', 2, ' ', 33 + atommer, ' ', 34 + atommer]])
        allbond = allbond + 1
        bondblock = ''.join([str(x) for x in [bondblock, lnovo, 'M  V30 ', 1 + allbond, ' ', 1, ' ', 33 + atommer, ' ', 35 + atommer]])
        allbond = allbond + 1
        bondblock = ''.join([str(x) for x in [bondblock, lnovo, 'M  V30 ', 1 + allbond, ' ', 1, ' ', 34 + atommer, ' ', 30 + atommer]])
        allbond = allbond + 1
        bondblock = ''.join([str(x) for x in [bondblock, lnovo, 'M  V30 ', 1 + allbond, ' ', 2, ' ', 35 + atommer, ' ', 36 + atommer]])
        allbond = allbond + 1
        bondblock = ''.join([str(x) for x in [bondblock, lnovo, 'M  V30 ', 1 + allbond, ' ', 1, ' ', 35 + atommer, ' ', 39 + atommer]])
        allbond = allbond + 1
        bondblock = ''.join([str(x) for x in [bondblock, lnovo, 'M  V30 ', 1 + allbond, ' ', 1, ' ', 36 + atommer, ' ', 37 + atommer]])
        allbond = allbond + 1
        bondblock = ''.join([str(x) for x in [bondblock, lnovo, 'M  V30 ', 1 + allbond, ' ', 2, ' ', 37 + atommer, ' ', 38 + atommer]])
        allbond = allbond + 1
        bondblock = ''.join([str(x) for x in [bondblock, lnovo, 'M  V30 ', 1 + allbond, ' ', 1, ' ', 38 + atommer, ' ', 34 + atommer]])
        allbond = allbond + 1
      elif monomer == 'T':
        atomblock = ''.join([str(x) for x in [atomblock, lnovo, 'M  V30 ', allatom + 1, ' ', 'C', ' ', 0.8418 + 6 * bpairmono, ' ', -4.9677, ' ', 0, ' ', 0]])
        allatom = allatom + 1
        atomblock = ''.join([str(x) for x in [atomblock, lnovo, 'M  V30 ', allatom + 1, ' ', 'C', ' ', 1.6508 + 6 * bpairmono, ' ', -5.5555, ' ', 0, ' ', 0]])
        allatom = allatom + 1
        atomblock = ''.join([str(x) for x in [atomblock, lnovo, 'M  V30 ', allatom + 1, ' ', 'C', ' ', 1.3418 + 6 * bpairmono, ' ', -6.5066, ' ', 0, ' ', 0]])
        allatom = allatom + 1
        atomblock = ''.join([str(x) for x in [atomblock, lnovo, 'M  V30 ', allatom + 1, ' ', 'C', ' ', 0.3418 + 6 * bpairmono, ' ', -6.5066, ' ', 0, ' ', 0]])
        allatom = allatom + 1
        atomblock = ''.join([str(x) for x in [atomblock, lnovo, 'M  V30 ', allatom + 1, ' ', 'O', ' ', 0.0327 + 6 * bpairmono, ' ', -5.5555, ' ', 0, ' ', 0]])
        allatom = allatom + 1
        atomblock = ''.join([str(x) for x in [atomblock, lnovo, 'M  V30 ', allatom + 1, ' ', 'O', ' ', 2.3418 + 6 * bpairmono, ' ', -6.5066, ' ', 0, ' ', 0]])
        allatom = allatom + 1
        atomblock = ''.join([str(x) for x in [atomblock, lnovo, 'M  V30 ', allatom + 1, ' ', 'C', ' ', -0.6582 + 6 * bpairmono, ' ', -6.5066, ' ', 0, ' ', 0]])
        allatom = allatom + 1
        atomblock = ''.join([str(x) for x in [atomblock, lnovo, 'M  V30 ', allatom + 1, ' ', 'O', ' ', -1.6582 + 6 * bpairmono, ' ', -6.5066, ' ', 0, ' ', 0]])
        allatom = allatom + 1
        atomblock = ''.join([str(x) for x in [atomblock, lnovo, 'M  V30 ', allatom + 1, ' ', 'P', ' ', -2.6582 + 6 * bpairmono, ' ', -6.5066, ' ', 0, ' ', 0]])
        allatom = allatom + 1
        if bpairmono > 0:
          bondblock = ''.join([str(x) for x in [bondblock, lnovo, 'M  V30 ', 1 + allbond, ' ', 1, ' ', merj1, ' ', allatom]])
          allbond = allbond + 1
        merj1 = allatom - 3
        atomblock = ''.join([str(x) for x in [atomblock, lnovo, 'M  V30 ', allatom + 1, ' ', 'C', ' ', 1.4262 + 6 * bpairmono, ' ', 4.9677, ' ', 0, ' ', 0]])
        allatom = allatom + 1
        atomblock = ''.join([str(x) for x in [atomblock, lnovo, 'M  V30 ', allatom + 1, ' ', 'C', ' ', 0.6172 + 6 * bpairmono, ' ', 5.5555, ' ', 0, ' ', 0]])
        allatom = allatom + 1
        atomblock = ''.join([str(x) for x in [atomblock, lnovo, 'M  V30 ', allatom + 1, ' ', 'C', ' ', 0.9262 + 6 * bpairmono, ' ', 6.5066, ' ', 0, ' ', 0]])
        allatom = allatom + 1
        atomblock = ''.join([str(x) for x in [atomblock, lnovo, 'M  V30 ', allatom + 1, ' ', 'C', ' ', 1.9262 + 6 * bpairmono, ' ', 6.5066, ' ', 0, ' ', 0]])
        allatom = allatom + 1
        atomblock = ''.join([str(x) for x in [atomblock, lnovo, 'M  V30 ', allatom + 1, ' ', 'O', ' ', 2.2352 + 6 * bpairmono, ' ', 5.5555, ' ', 0, ' ', 0]])
        allatom = allatom + 1
        atomblock = ''.join([str(x) for x in [atomblock, lnovo, 'M  V30 ', allatom + 1, ' ', 'C', ' ', 2.9262 + 6 * bpairmono, ' ', 6.5066, ' ', 0, ' ', 0]])
        allatom = allatom + 1
        atomblock = ''.join([str(x) for x in [atomblock, lnovo, 'M  V30 ', allatom + 1, ' ', 'O', ' ', -0.0738 + 6 * bpairmono, ' ', 6.5066, ' ', 0, ' ', 0]])
        allatom = allatom + 1
        atomblock = ''.join([str(x) for x in [atomblock, lnovo, 'M  V30 ', allatom + 1, ' ', 'P', ' ', -1.0738 + 6 * bpairmono, ' ', 6.5066, ' ', 0, ' ', 0]])
        allatom = allatom + 1
        atomblock = ''.join([str(x) for x in [atomblock, lnovo, 'M  V30 ', allatom + 1, ' ', 'O', ' ', -2.0738 + 6 * bpairmono, ' ', 6.5066, ' ', 0, ' ', 0]])
        allatom = allatom + 1
        if bpairmono > 0:
          bondblock = ''.join([str(x) for x in [bondblock, lnovo, 'M  V30 ', 1 + allbond, ' ', 1, ' ', merj2, ' ', allatom]])
          allbond = allbond + 1
        merj2 = allatom - 3
        atomblock = ''.join([str(x) for x in [atomblock, lnovo, 'M  V30 ', allatom + 1, ' ', 'O', ' ', -0.5738 + 6 * bpairmono, ' ', 7.3726, ' ', 0, ' ', 0]])
        allatom = allatom + 1
        atomblock = ''.join([str(x) for x in [atomblock, lnovo, 'M  V30 ', allatom + 1, ' ', 'O', ' ', -1.5738 + 6 * bpairmono, ' ', 7.3726, ' ', 0, ' ', 0, ' ', 'CHG=-1']])
        allatom = allatom + 1
        atomblock = ''.join([str(x) for x in [atomblock, lnovo, 'M  V30 ', allatom + 1, ' ', 'O', ' ', -2.1582 + 6 * bpairmono, ' ', -7.3726, ' ', 0, ' ', 0, ' ', 'CHG=-1']])
        allatom = allatom + 1
        atomblock = ''.join([str(x) for x in [atomblock, lnovo, 'M  V30 ', allatom + 1, ' ', 'O', ' ', -3.1582 + 6 * bpairmono, ' ', -7.3726, ' ', 0, ' ', 0]])
        allatom = allatom + 1
        atomblock = ''.join([str(x) for x in [atomblock, lnovo, 'M  V30 ', allatom + 1, ' ', 'N', ' ', 1.4262 + 6 * bpairmono, ' ', 3.9677, ' ', 0, ' ', 0]])
        allatom = allatom + 1
        atomblock = ''.join([str(x) for x in [atomblock, lnovo, 'M  V30 ', allatom + 1, ' ', 'C', ' ', 2.2922 + 6 * bpairmono, ' ', 3.4677, ' ', 0, ' ', 0]])
        allatom = allatom + 1
        atomblock = ''.join([str(x) for x in [atomblock, lnovo, 'M  V30 ', allatom + 1, ' ', 'C', ' ', 2.2922 + 6 * bpairmono, ' ', 2.4677, ' ', 0, ' ', 0]])
        allatom = allatom + 1
        atomblock = ''.join([str(x) for x in [atomblock, lnovo, 'M  V30 ', allatom + 1, ' ', 'C', ' ', 1.4262 + 6 * bpairmono, ' ', 1.9677, ' ', 0, ' ', 0]])
        allatom = allatom + 1
        atomblock = ''.join([str(x) for x in [atomblock, lnovo, 'M  V30 ', allatom + 1, ' ', 'N', ' ', 0.5602 + 6 * bpairmono, ' ', 2.4677, ' ', 0, ' ', 0]])
        allatom = allatom + 1
        atomblock = ''.join([str(x) for x in [atomblock, lnovo, 'M  V30 ', allatom + 1, ' ', 'C', ' ', 0.5602 + 6 * bpairmono, ' ', 3.4677, ' ', 0, ' ', 0]])
        allatom = allatom + 1
        atomblock = ''.join([str(x) for x in [atomblock, lnovo, 'M  V30 ', allatom + 1, ' ', 'O', ' ', 1.4262 + 6 * bpairmono, ' ', 0.9677, ' ', 0, ' ', 0]])
        allatom = allatom + 1
        atomblock = ''.join([str(x) for x in [atomblock, lnovo, 'M  V30 ', allatom + 1, ' ', 'O', ' ', -0.3059 + 6 * bpairmono, ' ', 3.9677, ' ', 0, ' ', 0]])
        allatom = allatom + 1
        atomblock = ''.join([str(x) for x in [atomblock, lnovo, 'M  V30 ', allatom + 1, ' ', 'C', ' ', 3.1582 + 6 * bpairmono, ' ', 1.9677, ' ', 0, ' ', 0]])
        allatom = allatom + 1
        atomblock = ''.join([str(x) for x in [atomblock, lnovo, 'M  V30 ', allatom + 1, ' ', 'N', ' ', 0.8418 + 6 * bpairmono, ' ', -3.9677, ' ', 0, ' ', 0]])
        allatom = allatom + 1
        atomblock = ''.join([str(x) for x in [atomblock, lnovo, 'M  V30 ', allatom + 1, ' ', 'C', ' ', 0.0327 + 6 * bpairmono, ' ', -3.3799, ' ', 0, ' ', 0]])
        allatom = allatom + 1
        atomblock = ''.join([str(x) for x in [atomblock, lnovo, 'M  V30 ', allatom + 1, ' ', 'C', ' ', 0.3418 + 6 * bpairmono, ' ', -2.4289, ' ', 0, ' ', 0]])
        allatom = allatom + 1
        atomblock = ''.join([str(x) for x in [atomblock, lnovo, 'M  V30 ', allatom + 1, ' ', 'N', ' ', 1.3418 + 6 * bpairmono, ' ', -2.4289, ' ', 0, ' ', 0]])
        allatom = allatom + 1
        atomblock = ''.join([str(x) for x in [atomblock, lnovo, 'M  V30 ', allatom + 1, ' ', 'C', ' ', 1.6508 + 6 * bpairmono, ' ', -3.3799, ' ', 0, ' ', 0]])
        allatom = allatom + 1
        atomblock = ''.join([str(x) for x in [atomblock, lnovo, 'M  V30 ', allatom + 1, ' ', 'N', ' ', -0.9454 + 6 * bpairmono, ' ', -3.5878, ' ', 0, ' ', 0]])
        allatom = allatom + 1
        atomblock = ''.join([str(x) for x in [atomblock, lnovo, 'M  V30 ', allatom + 1, ' ', 'C', ' ', -1.6145 + 6 * bpairmono, ' ', -2.8447, ' ', 0, ' ', 0]])
        allatom = allatom + 1
        atomblock = ''.join([str(x) for x in [atomblock, lnovo, 'M  V30 ', allatom + 1, ' ', 'N', ' ', -1.3055 + 6 * bpairmono, ' ', -1.8936, ' ', 0, ' ', 0]])
        allatom = allatom + 1
        atomblock = ''.join([str(x) for x in [atomblock, lnovo, 'M  V30 ', allatom + 1, ' ', 'C', ' ', -0.3274 + 6 * bpairmono, ' ', -1.6857, ' ', 0, ' ', 0]])
        allatom = allatom + 1
        atomblock = ''.join([str(x) for x in [atomblock, lnovo, 'M  V30 ', allatom + 1, ' ', 'N', ' ', 0.1726 + 6 * bpairmono, ' ', -0.8197, ' ', 0, ' ', 0]])
        allatom = allatom + 1
        bondblock = ''.join([str(x) for x in [bondblock, lnovo, 'M  V30 ', 1 + allbond, ' ', 1, ' ', 1 + atommer, ' ', 2 + atommer]])
        allbond = allbond + 1
        bondblock = ''.join([str(x) for x in [bondblock, lnovo, 'M  V30 ', 1 + allbond, ' ', 1, ' ', 1 + atommer, ' ', 32 + atommer]])
        allbond = allbond + 1
        bondblock = ''.join([str(x) for x in [bondblock, lnovo, 'M  V30 ', 1 + allbond, ' ', 1, ' ', 2 + atommer, ' ', 3 + atommer]])
        allbond = allbond + 1
        bondblock = ''.join([str(x) for x in [bondblock, lnovo, 'M  V30 ', 1 + allbond, ' ', 1, ' ', 3 + atommer, ' ', 4 + atommer]])
        allbond = allbond + 1
        bondblock = ''.join([str(x) for x in [bondblock, lnovo, 'M  V30 ', 1 + allbond, ' ', 1, ' ', 3 + atommer, ' ', 6 + atommer]])
        allbond = allbond + 1
        bondblock = ''.join([str(x) for x in [bondblock, lnovo, 'M  V30 ', 1 + allbond, ' ', 1, ' ', 4 + atommer, ' ', 5 + atommer]])
        allbond = allbond + 1
        bondblock = ''.join([str(x) for x in [bondblock, lnovo, 'M  V30 ', 1 + allbond, ' ', 1, ' ', 4 + atommer, ' ', 7 + atommer]])
        allbond = allbond + 1
        bondblock = ''.join([str(x) for x in [bondblock, lnovo, 'M  V30 ', 1 + allbond, ' ', 1, ' ', 5 + atommer, ' ', 1 + atommer]])
        allbond = allbond + 1
        bondblock = ''.join([str(x) for x in [bondblock, lnovo, 'M  V30 ', 1 + allbond, ' ', 1, ' ', 7 + atommer, ' ', 8 + atommer]])
        allbond = allbond + 1
        bondblock = ''.join([str(x) for x in [bondblock, lnovo, 'M  V30 ', 1 + allbond, ' ', 1, ' ', 8 + atommer, ' ', 9 + atommer]])
        allbond = allbond + 1
        bondblock = ''.join([str(x) for x in [bondblock, lnovo, 'M  V30 ', 1 + allbond, ' ', 1, ' ', 9 + atommer, ' ', 21 + atommer]])
        allbond = allbond + 1
        bondblock = ''.join([str(x) for x in [bondblock, lnovo, 'M  V30 ', 1 + allbond, ' ', 2, ' ', 9 + atommer, ' ', 22 + atommer]])
        allbond = allbond + 1
        bondblock = ''.join([str(x) for x in [bondblock, lnovo, 'M  V30 ', 1 + allbond, ' ', 1, ' ', 10 + atommer, ' ', 11 + atommer]])
        allbond = allbond + 1
        bondblock = ''.join([str(x) for x in [bondblock, lnovo, 'M  V30 ', 1 + allbond, ' ', 1, ' ', 10 + atommer, ' ', 23 + atommer]])
        allbond = allbond + 1
        bondblock = ''.join([str(x) for x in [bondblock, lnovo, 'M  V30 ', 1 + allbond, ' ', 1, ' ', 11 + atommer, ' ', 12 + atommer]])
        allbond = allbond + 1
        bondblock = ''.join([str(x) for x in [bondblock, lnovo, 'M  V30 ', 1 + allbond, ' ', 1, ' ', 12 + atommer, ' ', 13 + atommer]])
        allbond = allbond + 1
        bondblock = ''.join([str(x) for x in [bondblock, lnovo, 'M  V30 ', 1 + allbond, ' ', 1, ' ', 12 + atommer, ' ', 16 + atommer]])
        allbond = allbond + 1
        bondblock = ''.join([str(x) for x in [bondblock, lnovo, 'M  V30 ', 1 + allbond, ' ', 1, ' ', 13 + atommer, ' ', 14 + atommer]])
        allbond = allbond + 1
        bondblock = ''.join([str(x) for x in [bondblock, lnovo, 'M  V30 ', 1 + allbond, ' ', 1, ' ', 13 + atommer, ' ', 15 + atommer]])
        allbond = allbond + 1
        bondblock = ''.join([str(x) for x in [bondblock, lnovo, 'M  V30 ', 1 + allbond, ' ', 1, ' ', 14 + atommer, ' ', 10 + atommer]])
        allbond = allbond + 1
        bondblock = ''.join([str(x) for x in [bondblock, lnovo, 'M  V30 ', 1 + allbond, ' ', 1, ' ', 16 + atommer, ' ', 17 + atommer]])
        allbond = allbond + 1
        bondblock = ''.join([str(x) for x in [bondblock, lnovo, 'M  V30 ', 1 + allbond, ' ', 1, ' ', 17 + atommer, ' ', 18 + atommer]])
        allbond = allbond + 1
        bondblock = ''.join([str(x) for x in [bondblock, lnovo, 'M  V30 ', 1 + allbond, ' ', 2, ' ', 17 + atommer, ' ', 19 + atommer]])
        allbond = allbond + 1
        bondblock = ''.join([str(x) for x in [bondblock, lnovo, 'M  V30 ', 1 + allbond, ' ', 1, ' ', 20 + atommer, ' ', 17 + atommer]])
        allbond = allbond + 1
        bondblock = ''.join([str(x) for x in [bondblock, lnovo, 'M  V30 ', 1 + allbond, ' ', 1, ' ', 23 + atommer, ' ', 24 + atommer]])
        allbond = allbond + 1
        bondblock = ''.join([str(x) for x in [bondblock, lnovo, 'M  V30 ', 1 + allbond, ' ', 2, ' ', 24 + atommer, ' ', 25 + atommer]])
        allbond = allbond + 1
        bondblock = ''.join([str(x) for x in [bondblock, lnovo, 'M  V30 ', 1 + allbond, ' ', 1, ' ', 25 + atommer, ' ', 26 + atommer]])
        allbond = allbond + 1
        bondblock = ''.join([str(x) for x in [bondblock, lnovo, 'M  V30 ', 1 + allbond, ' ', 1, ' ', 25 + atommer, ' ', 31 + atommer]])
        allbond = allbond + 1
        bondblock = ''.join([str(x) for x in [bondblock, lnovo, 'M  V30 ', 1 + allbond, ' ', 1, ' ', 26 + atommer, ' ', 27 + atommer]])
        allbond = allbond + 1
        bondblock = ''.join([str(x) for x in [bondblock, lnovo, 'M  V30 ', 1 + allbond, ' ', 2, ' ', 26 + atommer, ' ', 29 + atommer]])
        allbond = allbond + 1
        bondblock = ''.join([str(x) for x in [bondblock, lnovo, 'M  V30 ', 1 + allbond, ' ', 1, ' ', 27 + atommer, ' ', 28 + atommer]])
        allbond = allbond + 1
        bondblock = ''.join([str(x) for x in [bondblock, lnovo, 'M  V30 ', 1 + allbond, ' ', 1, ' ', 28 + atommer, ' ', 23 + atommer]])
        allbond = allbond + 1
        bondblock = ''.join([str(x) for x in [bondblock, lnovo, 'M  V30 ', 1 + allbond, ' ', 2, ' ', 28 + atommer, ' ', 30 + atommer]])
        allbond = allbond + 1
        bondblock = ''.join([str(x) for x in [bondblock, lnovo, 'M  V30 ', 1 + allbond, ' ', 1, ' ', 32 + atommer, ' ', 33 + atommer]])
        allbond = allbond + 1
        bondblock = ''.join([str(x) for x in [bondblock, lnovo, 'M  V30 ', 1 + allbond, ' ', 2, ' ', 33 + atommer, ' ', 34 + atommer]])
        allbond = allbond + 1
        bondblock = ''.join([str(x) for x in [bondblock, lnovo, 'M  V30 ', 1 + allbond, ' ', 1, ' ', 33 + atommer, ' ', 37 + atommer]])
        allbond = allbond + 1
        bondblock = ''.join([str(x) for x in [bondblock, lnovo, 'M  V30 ', 1 + allbond, ' ', 1, ' ', 34 + atommer, ' ', 35 + atommer]])
        allbond = allbond + 1
        bondblock = ''.join([str(x) for x in [bondblock, lnovo, 'M  V30 ', 1 + allbond, ' ', 2, ' ', 35 + atommer, ' ', 36 + atommer]])
        allbond = allbond + 1
        bondblock = ''.join([str(x) for x in [bondblock, lnovo, 'M  V30 ', 1 + allbond, ' ', 1, ' ', 36 + atommer, ' ', 32 + atommer]])
        allbond = allbond + 1
        bondblock = ''.join([str(x) for x in [bondblock, lnovo, 'M  V30 ', 1 + allbond, ' ', 2, ' ', 37 + atommer, ' ', 38 + atommer]])
        allbond = allbond + 1
        bondblock = ''.join([str(x) for x in [bondblock, lnovo, 'M  V30 ', 1 + allbond, ' ', 1, ' ', 38 + atommer, ' ', 39 + atommer]])
        allbond = allbond + 1
        bondblock = ''.join([str(x) for x in [bondblock, lnovo, 'M  V30 ', 1 + allbond, ' ', 2, ' ', 39 + atommer, ' ', 40 + atommer]])
        allbond = allbond + 1
        bondblock = ''.join([str(x) for x in [bondblock, lnovo, 'M  V30 ', 1 + allbond, ' ', 1, ' ', 40 + atommer, ' ', 34 + atommer]])
        allbond = allbond + 1
        bondblock = ''.join([str(x) for x in [bondblock, lnovo, 'M  V30 ', 1 + allbond, ' ', 1, ' ', 40 + atommer, ' ', 41 + atommer]])
        allbond = allbond + 1
      elif monomer == 'C':
        atomblock = ''.join([str(x) for x in [atomblock, lnovo, 'M  V30 ', allatom + 1, ' ', 'C', ' ', 0.9578 + 6 * bpairmono, ' ', -4.9677, ' ', 0, ' ', 0]])
        allatom = allatom + 1
        atomblock = ''.join([str(x) for x in [atomblock, lnovo, 'M  V30 ', allatom + 1, ' ', 'C', ' ', 1.7668 + 6 * bpairmono, ' ', -5.5555, ' ', 0, ' ', 0]])
        allatom = allatom + 1
        atomblock = ''.join([str(x) for x in [atomblock, lnovo, 'M  V30 ', allatom + 1, ' ', 'C', ' ', 1.4578 + 6 * bpairmono, ' ', -6.5066, ' ', 0, ' ', 0]])
        allatom = allatom + 1
        atomblock = ''.join([str(x) for x in [atomblock, lnovo, 'M  V30 ', allatom + 1, ' ', 'C', ' ', 0.4578 + 6 * bpairmono, ' ', -6.5066, ' ', 0, ' ', 0]])
        allatom = allatom + 1
        atomblock = ''.join([str(x) for x in [atomblock, lnovo, 'M  V30 ', allatom + 1, ' ', 'O', ' ', 0.1488 + 6 * bpairmono, ' ', -5.5555, ' ', 0, ' ', 0]])
        allatom = allatom + 1
        atomblock = ''.join([str(x) for x in [atomblock, lnovo, 'M  V30 ', allatom + 1, ' ', 'O', ' ', 2.4578 + 6 * bpairmono, ' ', -6.5066, ' ', 0, ' ', 0]])
        allatom = allatom + 1
        atomblock = ''.join([str(x) for x in [atomblock, lnovo, 'M  V30 ', allatom + 1, ' ', 'C', ' ', -0.5422 + 6 * bpairmono, ' ', -6.5066, ' ', 0, ' ', 0]])
        allatom = allatom + 1
        atomblock = ''.join([str(x) for x in [atomblock, lnovo, 'M  V30 ', allatom + 1, ' ', 'O', ' ', -1.5422 + 6 * bpairmono, ' ', -6.5066, ' ', 0, ' ', 0]])
        allatom = allatom + 1
        atomblock = ''.join([str(x) for x in [atomblock, lnovo, 'M  V30 ', allatom + 1, ' ', 'P', ' ', -2.5422 + 6 * bpairmono, ' ', -6.5066, ' ', 0, ' ', 0]])
        allatom = allatom + 1
        if bpairmono > 0:
          bondblock = ''.join([str(x) for x in [bondblock, lnovo, 'M  V30 ', 1 + allbond, ' ', 1, ' ', merj1, ' ', allatom]])
          allbond = allbond + 1
        merj1 = allatom - 3
        atomblock = ''.join([str(x) for x in [atomblock, lnovo, 'M  V30 ', allatom + 1, ' ', 'C', ' ', 1.5422 + 6 * bpairmono, ' ', 4.9677, ' ', 0, ' ', 0]])
        allatom = allatom + 1
        atomblock = ''.join([str(x) for x in [atomblock, lnovo, 'M  V30 ', allatom + 1, ' ', 'C', ' ', 0.7332 + 6 * bpairmono, ' ', 5.5555, ' ', 0, ' ', 0]])
        allatom = allatom + 1
        atomblock = ''.join([str(x) for x in [atomblock, lnovo, 'M  V30 ', allatom + 1, ' ', 'C', ' ', 1.0422 + 6 * bpairmono, ' ', 6.5066, ' ', 0, ' ', 0]])
        allatom = allatom + 1
        atomblock = ''.join([str(x) for x in [atomblock, lnovo, 'M  V30 ', allatom + 1, ' ', 'C', ' ', 2.0422 + 6 * bpairmono, ' ', 6.5066, ' ', 0, ' ', 0]])
        allatom = allatom + 1
        atomblock = ''.join([str(x) for x in [atomblock, lnovo, 'M  V30 ', allatom + 1, ' ', 'O', ' ', 2.3512 + 6 * bpairmono, ' ', 5.5555, ' ', 0, ' ', 0]])
        allatom = allatom + 1
        atomblock = ''.join([str(x) for x in [atomblock, lnovo, 'M  V30 ', allatom + 1, ' ', 'C', ' ', 3.0422 + 6 * bpairmono, ' ', 6.5066, ' ', 0, ' ', 0]])
        allatom = allatom + 1
        atomblock = ''.join([str(x) for x in [atomblock, lnovo, 'M  V30 ', allatom + 1, ' ', 'O', ' ', 0.0422 + 6 * bpairmono, ' ', 6.5066, ' ', 0, ' ', 0]])
        allatom = allatom + 1
        atomblock = ''.join([str(x) for x in [atomblock, lnovo, 'M  V30 ', allatom + 1, ' ', 'P', ' ', -0.9578 + 6 * bpairmono, ' ', 6.5066, ' ', 0, ' ', 0]])
        allatom = allatom + 1
        atomblock = ''.join([str(x) for x in [atomblock, lnovo, 'M  V30 ', allatom + 1, ' ', 'O', ' ', -1.9578 + 6 * bpairmono, ' ', 6.5066, ' ', 0, ' ', 0]])
        allatom = allatom + 1
        if bpairmono > 0:
          bondblock = ''.join([str(x) for x in [bondblock, lnovo, 'M  V30 ', 1 + allbond, ' ', 1, ' ', merj2, ' ', allatom]])
          allbond = allbond + 1
        merj2 = allatom - 3
        atomblock = ''.join([str(x) for x in [atomblock, lnovo, 'M  V30 ', allatom + 1, ' ', 'O', ' ', -0.4578 + 6 * bpairmono, ' ', 7.3726, ' ', 0, ' ', 0]])
        allatom = allatom + 1
        atomblock = ''.join([str(x) for x in [atomblock, lnovo, 'M  V30 ', allatom + 1, ' ', 'O', ' ', -1.4578 + 6 * bpairmono, ' ', 7.3726, ' ', 0, ' ', 0, ' ', 'CHG=-1']])
        allatom = allatom + 1
        atomblock = ''.join([str(x) for x in [atomblock, lnovo, 'M  V30 ', allatom + 1, ' ', 'O', ' ', -2.0422 + 6 * bpairmono, ' ', -7.3726, ' ', 0, ' ', 0, ' ', 'CHG=-1']])
        allatom = allatom + 1
        atomblock = ''.join([str(x) for x in [atomblock, lnovo, 'M  V30 ', allatom + 1, ' ', 'O', ' ', -3.0422 + 6 * bpairmono, ' ', -7.3726, ' ', 0, ' ', 0]])
        allatom = allatom + 1
        atomblock = ''.join([str(x) for x in [atomblock, lnovo, 'M  V30 ', allatom + 1, ' ', 'N', ' ', 1.5422 + 6 * bpairmono, ' ', 3.9677, ' ', 0, ' ', 0]])
        allatom = allatom + 1
        atomblock = ''.join([str(x) for x in [atomblock, lnovo, 'M  V30 ', allatom + 1, ' ', 'C', ' ', 2.4082 + 6 * bpairmono, ' ', 3.4677, ' ', 0, ' ', 0]])
        allatom = allatom + 1
        atomblock = ''.join([str(x) for x in [atomblock, lnovo, 'M  V30 ', allatom + 1, ' ', 'C', ' ', 2.4082 + 6 * bpairmono, ' ', 2.4677, ' ', 0, ' ', 0]])
        allatom = allatom + 1
        atomblock = ''.join([str(x) for x in [atomblock, lnovo, 'M  V30 ', allatom + 1, ' ', 'C', ' ', 1.5422 + 6 * bpairmono, ' ', 1.9677, ' ', 0, ' ', 0]])
        allatom = allatom + 1
        atomblock = ''.join([str(x) for x in [atomblock, lnovo, 'M  V30 ', allatom + 1, ' ', 'N', ' ', 0.6762 + 6 * bpairmono, ' ', 2.4677, ' ', 0, ' ', 0]])
        allatom = allatom + 1
        atomblock = ''.join([str(x) for x in [atomblock, lnovo, 'M  V30 ', allatom + 1, ' ', 'C', ' ', 0.6762 + 6 * bpairmono, ' ', 3.4677, ' ', 0, ' ', 0]])
        allatom = allatom + 1
        atomblock = ''.join([str(x) for x in [atomblock, lnovo, 'M  V30 ', allatom + 1, ' ', 'N', ' ', 1.5422 + 6 * bpairmono, ' ', 0.9677, ' ', 0, ' ', 0]])
        allatom = allatom + 1
        atomblock = ''.join([str(x) for x in [atomblock, lnovo, 'M  V30 ', allatom + 1, ' ', 'N', ' ', 0.9578 + 6 * bpairmono, ' ', -3.9677, ' ', 0, ' ', 0]])
        allatom = allatom + 1
        atomblock = ''.join([str(x) for x in [atomblock, lnovo, 'M  V30 ', allatom + 1, ' ', 'C', ' ', 0.1488 + 6 * bpairmono, ' ', -3.3799, ' ', 0, ' ', 0]])
        allatom = allatom + 1
        atomblock = ''.join([str(x) for x in [atomblock, lnovo, 'M  V30 ', allatom + 1, ' ', 'C', ' ', 0.4578 + 6 * bpairmono, ' ', -2.4289, ' ', 0, ' ', 0]])
        allatom = allatom + 1
        atomblock = ''.join([str(x) for x in [atomblock, lnovo, 'M  V30 ', allatom + 1, ' ', 'N', ' ', 1.4578 + 6 * bpairmono, ' ', -2.4289, ' ', 0, ' ', 0]])
        allatom = allatom + 1
        atomblock = ''.join([str(x) for x in [atomblock, lnovo, 'M  V30 ', allatom + 1, ' ', 'C', ' ', 1.7668 + 6 * bpairmono, ' ', -3.3799, ' ', 0, ' ', 0]])
        allatom = allatom + 1
        atomblock = ''.join([str(x) for x in [atomblock, lnovo, 'M  V30 ', allatom + 1, ' ', 'N', ' ', -0.8294 + 6 * bpairmono, ' ', -3.5878, ' ', 0, ' ', 0]])
        allatom = allatom + 1
        atomblock = ''.join([str(x) for x in [atomblock, lnovo, 'M  V30 ', allatom + 1, ' ', 'C', ' ', -1.4985 + 6 * bpairmono, ' ', -2.8447, ' ', 0, ' ', 0]])
        allatom = allatom + 1
        atomblock = ''.join([str(x) for x in [atomblock, lnovo, 'M  V30 ', allatom + 1, ' ', 'N', ' ', -1.1895 + 6 * bpairmono, ' ', -1.8936, ' ', 0, ' ', 0]])
        allatom = allatom + 1
        atomblock = ''.join([str(x) for x in [atomblock, lnovo, 'M  V30 ', allatom + 1, ' ', 'C', ' ', -0.2113 + 6 * bpairmono, ' ', -1.6857, ' ', 0, ' ', 0]])
        allatom = allatom + 1
        atomblock = ''.join([str(x) for x in [atomblock, lnovo, 'M  V30 ', allatom + 1, ' ', 'O', ' ', 0.2887 + 6 * bpairmono, ' ', -0.8197, ' ', 0, ' ', 0]])
        allatom = allatom + 1
        atomblock = ''.join([str(x) for x in [atomblock, lnovo, 'M  V30 ', allatom + 1, ' ', 'O', ' ', -0.1898 + 6 * bpairmono, ' ', 3.9677, ' ', 0, ' ', 0]])
        allatom = allatom + 1
        atomblock = ''.join([str(x) for x in [atomblock, lnovo, 'M  V30 ', allatom + 1, ' ', 'N', ' ', -1.7741 + 6 * bpairmono, ' ', -3.8768, ' ', 0, ' ', 0]])
        allatom = allatom + 1
        bondblock = ''.join([str(x) for x in [bondblock, lnovo, 'M  V30 ', 1 + allbond, ' ', 1, ' ', 1 + atommer, ' ', 2 + atommer]])
        allbond = allbond + 1
        bondblock = ''.join([str(x) for x in [bondblock, lnovo, 'M  V30 ', 1 + allbond, ' ', 1, ' ', 1 + atommer, ' ', 30 + atommer]])
        allbond = allbond + 1
        bondblock = ''.join([str(x) for x in [bondblock, lnovo, 'M  V30 ', 1 + allbond, ' ', 1, ' ', 2 + atommer, ' ', 3 + atommer]])
        allbond = allbond + 1
        bondblock = ''.join([str(x) for x in [bondblock, lnovo, 'M  V30 ', 1 + allbond, ' ', 1, ' ', 3 + atommer, ' ', 4 + atommer]])
        allbond = allbond + 1
        bondblock = ''.join([str(x) for x in [bondblock, lnovo, 'M  V30 ', 1 + allbond, ' ', 1, ' ', 3 + atommer, ' ', 6 + atommer]])
        allbond = allbond + 1
        bondblock = ''.join([str(x) for x in [bondblock, lnovo, 'M  V30 ', 1 + allbond, ' ', 1, ' ', 4 + atommer, ' ', 5 + atommer]])
        allbond = allbond + 1
        bondblock = ''.join([str(x) for x in [bondblock, lnovo, 'M  V30 ', 1 + allbond, ' ', 1, ' ', 4 + atommer, ' ', 7 + atommer]])
        allbond = allbond + 1
        bondblock = ''.join([str(x) for x in [bondblock, lnovo, 'M  V30 ', 1 + allbond, ' ', 1, ' ', 5 + atommer, ' ', 1 + atommer]])
        allbond = allbond + 1
        bondblock = ''.join([str(x) for x in [bondblock, lnovo, 'M  V30 ', 1 + allbond, ' ', 1, ' ', 7 + atommer, ' ', 8 + atommer]])
        allbond = allbond + 1
        bondblock = ''.join([str(x) for x in [bondblock, lnovo, 'M  V30 ', 1 + allbond, ' ', 1, ' ', 8 + atommer, ' ', 9 + atommer]])
        allbond = allbond + 1
        bondblock = ''.join([str(x) for x in [bondblock, lnovo, 'M  V30 ', 1 + allbond, ' ', 1, ' ', 9 + atommer, ' ', 21 + atommer]])
        allbond = allbond + 1
        bondblock = ''.join([str(x) for x in [bondblock, lnovo, 'M  V30 ', 1 + allbond, ' ', 2, ' ', 9 + atommer, ' ', 22 + atommer]])
        allbond = allbond + 1
        bondblock = ''.join([str(x) for x in [bondblock, lnovo, 'M  V30 ', 1 + allbond, ' ', 1, ' ', 10 + atommer, ' ', 11 + atommer]])
        allbond = allbond + 1
        bondblock = ''.join([str(x) for x in [bondblock, lnovo, 'M  V30 ', 1 + allbond, ' ', 1, ' ', 10 + atommer, ' ', 23 + atommer]])
        allbond = allbond + 1
        bondblock = ''.join([str(x) for x in [bondblock, lnovo, 'M  V30 ', 1 + allbond, ' ', 1, ' ', 11 + atommer, ' ', 12 + atommer]])
        allbond = allbond + 1
        bondblock = ''.join([str(x) for x in [bondblock, lnovo, 'M  V30 ', 1 + allbond, ' ', 1, ' ', 12 + atommer, ' ', 13 + atommer]])
        allbond = allbond + 1
        bondblock = ''.join([str(x) for x in [bondblock, lnovo, 'M  V30 ', 1 + allbond, ' ', 1, ' ', 12 + atommer, ' ', 16 + atommer]])
        allbond = allbond + 1
        bondblock = ''.join([str(x) for x in [bondblock, lnovo, 'M  V30 ', 1 + allbond, ' ', 1, ' ', 13 + atommer, ' ', 14 + atommer]])
        allbond = allbond + 1
        bondblock = ''.join([str(x) for x in [bondblock, lnovo, 'M  V30 ', 1 + allbond, ' ', 1, ' ', 13 + atommer, ' ', 15 + atommer]])
        allbond = allbond + 1
        bondblock = ''.join([str(x) for x in [bondblock, lnovo, 'M  V30 ', 1 + allbond, ' ', 1, ' ', 14 + atommer, ' ', 10 + atommer]])
        allbond = allbond + 1
        bondblock = ''.join([str(x) for x in [bondblock, lnovo, 'M  V30 ', 1 + allbond, ' ', 1, ' ', 16 + atommer, ' ', 17 + atommer]])
        allbond = allbond + 1
        bondblock = ''.join([str(x) for x in [bondblock, lnovo, 'M  V30 ', 1 + allbond, ' ', 1, ' ', 17 + atommer, ' ', 18 + atommer]])
        allbond = allbond + 1
        bondblock = ''.join([str(x) for x in [bondblock, lnovo, 'M  V30 ', 1 + allbond, ' ', 2, ' ', 17 + atommer, ' ', 19 + atommer]])
        allbond = allbond + 1
        bondblock = ''.join([str(x) for x in [bondblock, lnovo, 'M  V30 ', 1 + allbond, ' ', 1, ' ', 20 + atommer, ' ', 17 + atommer]])
        allbond = allbond + 1
        bondblock = ''.join([str(x) for x in [bondblock, lnovo, 'M  V30 ', 1 + allbond, ' ', 1, ' ', 23 + atommer, ' ', 24 + atommer]])
        allbond = allbond + 1
        bondblock = ''.join([str(x) for x in [bondblock, lnovo, 'M  V30 ', 1 + allbond, ' ', 2, ' ', 24 + atommer, ' ', 25 + atommer]])
        allbond = allbond + 1
        bondblock = ''.join([str(x) for x in [bondblock, lnovo, 'M  V30 ', 1 + allbond, ' ', 1, ' ', 25 + atommer, ' ', 26 + atommer]])
        allbond = allbond + 1
        bondblock = ''.join([str(x) for x in [bondblock, lnovo, 'M  V30 ', 1 + allbond, ' ', 2, ' ', 26 + atommer, ' ', 27 + atommer]])
        allbond = allbond + 1
        bondblock = ''.join([str(x) for x in [bondblock, lnovo, 'M  V30 ', 1 + allbond, ' ', 1, ' ', 26 + atommer, ' ', 29 + atommer]])
        allbond = allbond + 1
        bondblock = ''.join([str(x) for x in [bondblock, lnovo, 'M  V30 ', 1 + allbond, ' ', 1, ' ', 27 + atommer, ' ', 28 + atommer]])
        allbond = allbond + 1
        bondblock = ''.join([str(x) for x in [bondblock, lnovo, 'M  V30 ', 1 + allbond, ' ', 1, ' ', 28 + atommer, ' ', 23 + atommer]])
        allbond = allbond + 1
        bondblock = ''.join([str(x) for x in [bondblock, lnovo, 'M  V30 ', 1 + allbond, ' ', 2, ' ', 28 + atommer, ' ', 40 + atommer]])
        allbond = allbond + 1
        bondblock = ''.join([str(x) for x in [bondblock, lnovo, 'M  V30 ', 1 + allbond, ' ', 1, ' ', 30 + atommer, ' ', 31 + atommer]])
        allbond = allbond + 1
        bondblock = ''.join([str(x) for x in [bondblock, lnovo, 'M  V30 ', 1 + allbond, ' ', 2, ' ', 31 + atommer, ' ', 32 + atommer]])
        allbond = allbond + 1
        bondblock = ''.join([str(x) for x in [bondblock, lnovo, 'M  V30 ', 1 + allbond, ' ', 1, ' ', 31 + atommer, ' ', 35 + atommer]])
        allbond = allbond + 1
        bondblock = ''.join([str(x) for x in [bondblock, lnovo, 'M  V30 ', 1 + allbond, ' ', 1, ' ', 32 + atommer, ' ', 33 + atommer]])
        allbond = allbond + 1
        bondblock = ''.join([str(x) for x in [bondblock, lnovo, 'M  V30 ', 1 + allbond, ' ', 2, ' ', 33 + atommer, ' ', 34 + atommer]])
        allbond = allbond + 1
        bondblock = ''.join([str(x) for x in [bondblock, lnovo, 'M  V30 ', 1 + allbond, ' ', 1, ' ', 34 + atommer, ' ', 30 + atommer]])
        allbond = allbond + 1
        bondblock = ''.join([str(x) for x in [bondblock, lnovo, 'M  V30 ', 1 + allbond, ' ', 2, ' ', 35 + atommer, ' ', 36 + atommer]])
        allbond = allbond + 1
        bondblock = ''.join([str(x) for x in [bondblock, lnovo, 'M  V30 ', 1 + allbond, ' ', 1, ' ', 36 + atommer, ' ', 37 + atommer]])
        allbond = allbond + 1
        bondblock = ''.join([str(x) for x in [bondblock, lnovo, 'M  V30 ', 1 + allbond, ' ', 1, ' ', 36 + atommer, ' ', 41 + atommer]])
        allbond = allbond + 1
        bondblock = ''.join([str(x) for x in [bondblock, lnovo, 'M  V30 ', 1 + allbond, ' ', 1, ' ', 37 + atommer, ' ', 38 + atommer]])
        allbond = allbond + 1
        bondblock = ''.join([str(x) for x in [bondblock, lnovo, 'M  V30 ', 1 + allbond, ' ', 1, ' ', 38 + atommer, ' ', 32 + atommer]])
        allbond = allbond + 1
        bondblock = ''.join([str(x) for x in [bondblock, lnovo, 'M  V30 ', 1 + allbond, ' ', 2, ' ', 38 + atommer, ' ', 39 + atommer]])
        allbond = allbond + 1
      else:
        atomblock = ''.join([str(x) for x in [atomblock, lnovo, 'M  V30 ', allatom + 1, ' ', 'C', ' ', 0.9578 + 6 * bpairmono, ' ', -1.9677, ' ', 0, ' ', 0]])
        allatom = allatom + 1
        atomblock = ''.join([str(x) for x in [atomblock, lnovo, 'M  V30 ', allatom + 1, ' ', 'C', ' ', 1.8238 + 6 * bpairmono, ' ', -2.4677, ' ', 0, ' ', 0]])
        allatom = allatom + 1
        atomblock = ''.join([str(x) for x in [atomblock, lnovo, 'M  V30 ', allatom + 1, ' ', 'C', ' ', 1.8238 + 6 * bpairmono, ' ', -3.4677, ' ', 0, ' ', 0]])
        allatom = allatom + 1
        atomblock = ''.join([str(x) for x in [atomblock, lnovo, 'M  V30 ', allatom + 1, ' ', 'N', ' ', 0.9578 + 6 * bpairmono, ' ', -3.9677, ' ', 0, ' ', 0]])
        allatom = allatom + 1
        atomblock = ''.join([str(x) for x in [atomblock, lnovo, 'M  V30 ', allatom + 1, ' ', 'C', ' ', 0.0918 + 6 * bpairmono, ' ', -3.4677, ' ', 0, ' ', 0]])
        allatom = allatom + 1
        atomblock = ''.join([str(x) for x in [atomblock, lnovo, 'M  V30 ', allatom + 1, ' ', 'N', ' ', 0.0918 + 6 * bpairmono, ' ', -2.4677, ' ', 0, ' ', 0]])
        allatom = allatom + 1
        atomblock = ''.join([str(x) for x in [atomblock, lnovo, 'M  V30 ', allatom + 1, ' ', 'N', ' ', 0.9578 + 6 * bpairmono, ' ', -0.9677, ' ', 0, ' ', 0]])
        allatom = allatom + 1
        atomblock = ''.join([str(x) for x in [atomblock, lnovo, 'M  V30 ', allatom + 1, ' ', 'O', ' ', -0.7743 + 6 * bpairmono, ' ', -3.9677, ' ', 0, ' ', 0]])
        allatom = allatom + 1
        atomblock = ''.join([str(x) for x in [atomblock, lnovo, 'M  V30 ', allatom + 1, ' ', 'C', ' ', 0.9578 + 6 * bpairmono, ' ', -4.9677, ' ', 0, ' ', 0]])
        allatom = allatom + 1
        atomblock = ''.join([str(x) for x in [atomblock, lnovo, 'M  V30 ', allatom + 1, ' ', 'C', ' ', 1.7668 + 6 * bpairmono, ' ', -5.5555, ' ', 0, ' ', 0]])
        allatom = allatom + 1
        atomblock = ''.join([str(x) for x in [atomblock, lnovo, 'M  V30 ', allatom + 1, ' ', 'C', ' ', 1.4578 + 6 * bpairmono, ' ', -6.5066, ' ', 0, ' ', 0]])
        allatom = allatom + 1
        atomblock = ''.join([str(x) for x in [atomblock, lnovo, 'M  V30 ', allatom + 1, ' ', 'C', ' ', 0.4578 + 6 * bpairmono, ' ', -6.5066, ' ', 0, ' ', 0]])
        allatom = allatom + 1
        atomblock = ''.join([str(x) for x in [atomblock, lnovo, 'M  V30 ', allatom + 1, ' ', 'O', ' ', 0.1488 + 6 * bpairmono, ' ', -5.5555, ' ', 0, ' ', 0]])
        allatom = allatom + 1
        atomblock = ''.join([str(x) for x in [atomblock, lnovo, 'M  V30 ', allatom + 1, ' ', 'O', ' ', 2.4578 + 6 * bpairmono, ' ', -6.5066, ' ', 0, ' ', 0]])
        allatom = allatom + 1
        atomblock = ''.join([str(x) for x in [atomblock, lnovo, 'M  V30 ', allatom + 1, ' ', 'C', ' ', -0.5422 + 6 * bpairmono, ' ', -6.5066, ' ', 0, ' ', 0]])
        allatom = allatom + 1
        atomblock = ''.join([str(x) for x in [atomblock, lnovo, 'M  V30 ', allatom + 1, ' ', 'O', ' ', -1.5422 + 6 * bpairmono, ' ', -6.5066, ' ', 0, ' ', 0]])
        allatom = allatom + 1
        atomblock = ''.join([str(x) for x in [atomblock, lnovo, 'M  V30 ', allatom + 1, ' ', 'P', ' ', -2.5422 + 6 * bpairmono, ' ', -6.5066, ' ', 0, ' ', 0]])
        allatom = allatom + 1
        if bpairmono > 0:
          bondblock = ''.join([str(x) for x in [bondblock, lnovo, 'M  V30 ', 1 + allbond, ' ', 1, ' ', merj1, ' ', allatom]])
          allbond = allbond + 1
        merj1 = allatom - 3
        atomblock = ''.join([str(x) for x in [atomblock, lnovo, 'M  V30 ', allatom + 1, ' ', 'C', ' ', 1.5422 + 6 * bpairmono, ' ', 4.9677, ' ', 0, ' ', 0]])
        allatom = allatom + 1
        atomblock = ''.join([str(x) for x in [atomblock, lnovo, 'M  V30 ', allatom + 1, ' ', 'C', ' ', 0.7332 + 6 * bpairmono, ' ', 5.5555, ' ', 0, ' ', 0]])
        allatom = allatom + 1
        atomblock = ''.join([str(x) for x in [atomblock, lnovo, 'M  V30 ', allatom + 1, ' ', 'C', ' ', 1.0422 + 6 * bpairmono, ' ', 6.5066, ' ', 0, ' ', 0]])
        allatom = allatom + 1
        atomblock = ''.join([str(x) for x in [atomblock, lnovo, 'M  V30 ', allatom + 1, ' ', 'C', ' ', 2.0422 + 6 * bpairmono, ' ', 6.5066, ' ', 0, ' ', 0]])
        allatom = allatom + 1
        atomblock = ''.join([str(x) for x in [atomblock, lnovo, 'M  V30 ', allatom + 1, ' ', 'O', ' ', 2.3512 + 6 * bpairmono, ' ', 5.5555, ' ', 0, ' ', 0]])
        allatom = allatom + 1
        atomblock = ''.join([str(x) for x in [atomblock, lnovo, 'M  V30 ', allatom + 1, ' ', 'C', ' ', 3.0422 + 6 * bpairmono, ' ', 6.5066, ' ', 0, ' ', 0]])
        allatom = allatom + 1
        atomblock = ''.join([str(x) for x in [atomblock, lnovo, 'M  V30 ', allatom + 1, ' ', 'O', ' ', 0.0422 + 6 * bpairmono, ' ', 6.5066, ' ', 0, ' ', 0]])
        allatom = allatom + 1
        atomblock = ''.join([str(x) for x in [atomblock, lnovo, 'M  V30 ', allatom + 1, ' ', 'P', ' ', -0.9578 + 6 * bpairmono, ' ', 6.5066, ' ', 0, ' ', 0]])
        allatom = allatom + 1
        atomblock = ''.join([str(x) for x in [atomblock, lnovo, 'M  V30 ', allatom + 1, ' ', 'O', ' ', -1.9578 + 6 * bpairmono, ' ', 6.5066, ' ', 0, ' ', 0]])
        allatom = allatom + 1
        if bpairmono > 0:
          bondblock = ''.join([str(x) for x in [bondblock, lnovo, 'M  V30 ', 1 + allbond, ' ', 1, ' ', merj2, ' ', allatom]])
          allbond = allbond + 1
        merj2 = allatom - 3
        atomblock = ''.join([str(x) for x in [atomblock, lnovo, 'M  V30 ', allatom + 1, ' ', 'O', ' ', -0.4578 + 6 * bpairmono, ' ', 7.3726, ' ', 0, ' ', 0]])
        allatom = allatom + 1
        atomblock = ''.join([str(x) for x in [atomblock, lnovo, 'M  V30 ', allatom + 1, ' ', 'O', ' ', -1.4578 + 6 * bpairmono, ' ', 7.3726, ' ', 0, ' ', 0, ' ', 'CHG=-1']])
        allatom = allatom + 1
        atomblock = ''.join([str(x) for x in [atomblock, lnovo, 'M  V30 ', allatom + 1, ' ', 'N', ' ', 1.5422 + 6 * bpairmono, ' ', 3.9677, ' ', 0, ' ', 0]])
        allatom = allatom + 1
        atomblock = ''.join([str(x) for x in [atomblock, lnovo, 'M  V30 ', allatom + 1, ' ', 'C', ' ', 2.3512 + 6 * bpairmono, ' ', 3.3799, ' ', 0, ' ', 0]])
        allatom = allatom + 1
        atomblock = ''.join([str(x) for x in [atomblock, lnovo, 'M  V30 ', allatom + 1, ' ', 'N', ' ', 2.0422 + 6 * bpairmono, ' ', 2.4289, ' ', 0, ' ', 0]])
        allatom = allatom + 1
        atomblock = ''.join([str(x) for x in [atomblock, lnovo, 'M  V30 ', allatom + 1, ' ', 'C', ' ', 1.0422 + 6 * bpairmono, ' ', 2.4289, ' ', 0, ' ', 0]])
        allatom = allatom + 1
        atomblock = ''.join([str(x) for x in [atomblock, lnovo, 'M  V30 ', allatom + 1, ' ', 'C', ' ', 0.7332 + 6 * bpairmono, ' ', 3.3799, ' ', 0, ' ', 0]])
        allatom = allatom + 1
        atomblock = ''.join([str(x) for x in [atomblock, lnovo, 'M  V30 ', allatom + 1, ' ', 'C', ' ', 0.3731 + 6 * bpairmono, ' ', 1.6857, ' ', 0, ' ', 0]])
        allatom = allatom + 1
        atomblock = ''.join([str(x) for x in [atomblock, lnovo, 'M  V30 ', allatom + 1, ' ', 'N', ' ', -0.6051 + 6 * bpairmono, ' ', 1.8936, ' ', 0, ' ', 0]])
        allatom = allatom + 1
        atomblock = ''.join([str(x) for x in [atomblock, lnovo, 'M  V30 ', allatom + 1, ' ', 'C', ' ', -0.9141 + 6 * bpairmono, ' ', 2.8447, ' ', 0, ' ', 0]])
        allatom = allatom + 1
        atomblock = ''.join([str(x) for x in [atomblock, lnovo, 'M  V30 ', allatom + 1, ' ', 'N', ' ', -0.2450 + 6 * bpairmono, ' ', 3.5878, ' ', 0, ' ', 0]])
        allatom = allatom + 1
        atomblock = ''.join([str(x) for x in [atomblock, lnovo, 'M  V30 ', allatom + 1, ' ', 'O', ' ', 0.8731 + 6 * bpairmono, ' ', 0.8197, ' ', 0, ' ', 0]])
        allatom = allatom + 1
        atomblock = ''.join([str(x) for x in [atomblock, lnovo, 'M  V30 ', allatom + 1, ' ', 'O', ' ', -2.0422 + 6 * bpairmono, ' ', -7.3726, ' ', 0, ' ', 0, ' ', 'CHG=-1']])
        allatom = allatom + 1
        atomblock = ''.join([str(x) for x in [atomblock, lnovo, 'M  V30 ', allatom + 1, ' ', 'O', ' ', -3.0422 + 6 * bpairmono, ' ', -7.3726, ' ', 0, ' ', 0]])
        allatom = allatom + 1
        atomblock = ''.join([str(x) for x in [atomblock, lnovo, 'M  V30 ', allatom + 1, ' ', 'N', ' ', -1.4141 + 6 * bpairmono, ' ', 3.7107, ' ', 0, ' ', 0]])
        allatom = allatom + 1
        bondblock = ''.join([str(x) for x in [bondblock, lnovo, 'M  V30 ', 1 + allbond, ' ', 1, ' ', 1 + atommer, ' ', 2 + atommer]])
        allbond = allbond + 1
        bondblock = ''.join([str(x) for x in [bondblock, lnovo, 'M  V30 ', 1 + allbond, ' ', 1, ' ', 1 + atommer, ' ', 7 + atommer]])
        allbond = allbond + 1
        bondblock = ''.join([str(x) for x in [bondblock, lnovo, 'M  V30 ', 1 + allbond, ' ', 2, ' ', 2 + atommer, ' ', 3 + atommer]])
        allbond = allbond + 1
        bondblock = ''.join([str(x) for x in [bondblock, lnovo, 'M  V30 ', 1 + allbond, ' ', 1, ' ', 3 + atommer, ' ', 4 + atommer]])
        allbond = allbond + 1
        bondblock = ''.join([str(x) for x in [bondblock, lnovo, 'M  V30 ', 1 + allbond, ' ', 1, ' ', 4 + atommer, ' ', 5 + atommer]])
        allbond = allbond + 1
        bondblock = ''.join([str(x) for x in [bondblock, lnovo, 'M  V30 ', 1 + allbond, ' ', 1, ' ', 5 + atommer, ' ', 6 + atommer]])
        allbond = allbond + 1
        bondblock = ''.join([str(x) for x in [bondblock, lnovo, 'M  V30 ', 1 + allbond, ' ', 2, ' ', 5 + atommer, ' ', 8 + atommer]])
        allbond = allbond + 1
        bondblock = ''.join([str(x) for x in [bondblock, lnovo, 'M  V30 ', 1 + allbond, ' ', 2, ' ', 6 + atommer, ' ', 1 + atommer]])
        allbond = allbond + 1
        bondblock = ''.join([str(x) for x in [bondblock, lnovo, 'M  V30 ', 1 + allbond, ' ', 1, ' ', 9 + atommer, ' ', 4 + atommer]])
        allbond = allbond + 1
        bondblock = ''.join([str(x) for x in [bondblock, lnovo, 'M  V30 ', 1 + allbond, ' ', 1, ' ', 9 + atommer, ' ', 10 + atommer]])
        allbond = allbond + 1
        bondblock = ''.join([str(x) for x in [bondblock, lnovo, 'M  V30 ', 1 + allbond, ' ', 1, ' ', 10 + atommer, ' ', 11 + atommer]])
        allbond = allbond + 1
        bondblock = ''.join([str(x) for x in [bondblock, lnovo, 'M  V30 ', 1 + allbond, ' ', 1, ' ', 11 + atommer, ' ', 12 + atommer]])
        allbond = allbond + 1
        bondblock = ''.join([str(x) for x in [bondblock, lnovo, 'M  V30 ', 1 + allbond, ' ', 1, ' ', 11 + atommer, ' ', 14 + atommer]])
        allbond = allbond + 1
        bondblock = ''.join([str(x) for x in [bondblock, lnovo, 'M  V30 ', 1 + allbond, ' ', 1, ' ', 12 + atommer, ' ', 13 + atommer]])
        allbond = allbond + 1
        bondblock = ''.join([str(x) for x in [bondblock, lnovo, 'M  V30 ', 1 + allbond, ' ', 1, ' ', 12 + atommer, ' ', 15 + atommer]])
        allbond = allbond + 1
        bondblock = ''.join([str(x) for x in [bondblock, lnovo, 'M  V30 ', 1 + allbond, ' ', 1, ' ', 13 + atommer, ' ', 9 + atommer]])
        allbond = allbond + 1
        bondblock = ''.join([str(x) for x in [bondblock, lnovo, 'M  V30 ', 1 + allbond, ' ', 1, ' ', 15 + atommer, ' ', 16 + atommer]])
        allbond = allbond + 1
        bondblock = ''.join([str(x) for x in [bondblock, lnovo, 'M  V30 ', 1 + allbond, ' ', 1, ' ', 16 + atommer, ' ', 17 + atommer]])
        allbond = allbond + 1
        bondblock = ''.join([str(x) for x in [bondblock, lnovo, 'M  V30 ', 1 + allbond, ' ', 1, ' ', 17 + atommer, ' ', 39 + atommer]])
        allbond = allbond + 1
        bondblock = ''.join([str(x) for x in [bondblock, lnovo, 'M  V30 ', 1 + allbond, ' ', 2, ' ', 17 + atommer, ' ', 40 + atommer]])
        allbond = allbond + 1
        bondblock = ''.join([str(x) for x in [bondblock, lnovo, 'M  V30 ', 1 + allbond, ' ', 1, ' ', 18 + atommer, ' ', 19 + atommer]])
        allbond = allbond + 1
        bondblock = ''.join([str(x) for x in [bondblock, lnovo, 'M  V30 ', 1 + allbond, ' ', 1, ' ', 18 + atommer, ' ', 29 + atommer]])
        allbond = allbond + 1
        bondblock = ''.join([str(x) for x in [bondblock, lnovo, 'M  V30 ', 1 + allbond, ' ', 1, ' ', 19 + atommer, ' ', 20 + atommer]])
        allbond = allbond + 1
        bondblock = ''.join([str(x) for x in [bondblock, lnovo, 'M  V30 ', 1 + allbond, ' ', 1, ' ', 20 + atommer, ' ', 21 + atommer]])
        allbond = allbond + 1
        bondblock = ''.join([str(x) for x in [bondblock, lnovo, 'M  V30 ', 1 + allbond, ' ', 1, ' ', 20 + atommer, ' ', 24 + atommer]])
        allbond = allbond + 1
        bondblock = ''.join([str(x) for x in [bondblock, lnovo, 'M  V30 ', 1 + allbond, ' ', 1, ' ', 21 + atommer, ' ', 22 + atommer]])
        allbond = allbond + 1
        bondblock = ''.join([str(x) for x in [bondblock, lnovo, 'M  V30 ', 1 + allbond, ' ', 1, ' ', 21 + atommer, ' ', 23 + atommer]])
        allbond = allbond + 1
        bondblock = ''.join([str(x) for x in [bondblock, lnovo, 'M  V30 ', 1 + allbond, ' ', 1, ' ', 22 + atommer, ' ', 18 + atommer]])
        allbond = allbond + 1
        bondblock = ''.join([str(x) for x in [bondblock, lnovo, 'M  V30 ', 1 + allbond, ' ', 1, ' ', 24 + atommer, ' ', 25 + atommer]])
        allbond = allbond + 1
        bondblock = ''.join([str(x) for x in [bondblock, lnovo, 'M  V30 ', 1 + allbond, ' ', 1, ' ', 25 + atommer, ' ', 26 + atommer]])
        allbond = allbond + 1
        bondblock = ''.join([str(x) for x in [bondblock, lnovo, 'M  V30 ', 1 + allbond, ' ', 2, ' ', 25 + atommer, ' ', 27 + atommer]])
        allbond = allbond + 1
        bondblock = ''.join([str(x) for x in [bondblock, lnovo, 'M  V30 ', 1 + allbond, ' ', 1, ' ', 28 + atommer, ' ', 25 + atommer]])
        allbond = allbond + 1
        bondblock = ''.join([str(x) for x in [bondblock, lnovo, 'M  V30 ', 1 + allbond, ' ', 1, ' ', 29 + atommer, ' ', 30 + atommer]])
        allbond = allbond + 1
        bondblock = ''.join([str(x) for x in [bondblock, lnovo, 'M  V30 ', 1 + allbond, ' ', 2, ' ', 30 + atommer, ' ', 31 + atommer]])
        allbond = allbond + 1
        bondblock = ''.join([str(x) for x in [bondblock, lnovo, 'M  V30 ', 1 + allbond, ' ', 1, ' ', 31 + atommer, ' ', 32 + atommer]])
        allbond = allbond + 1
        bondblock = ''.join([str(x) for x in [bondblock, lnovo, 'M  V30 ', 1 + allbond, ' ', 2, ' ', 32 + atommer, ' ', 33 + atommer]])
        allbond = allbond + 1
        bondblock = ''.join([str(x) for x in [bondblock, lnovo, 'M  V30 ', 1 + allbond, ' ', 1, ' ', 32 + atommer, ' ', 34 + atommer]])
        allbond = allbond + 1
        bondblock = ''.join([str(x) for x in [bondblock, lnovo, 'M  V30 ', 1 + allbond, ' ', 1, ' ', 33 + atommer, ' ', 29 + atommer]])
        allbond = allbond + 1
        bondblock = ''.join([str(x) for x in [bondblock, lnovo, 'M  V30 ', 1 + allbond, ' ', 1, ' ', 34 + atommer, ' ', 35 + atommer]])
        allbond = allbond + 1
        bondblock = ''.join([str(x) for x in [bondblock, lnovo, 'M  V30 ', 1 + allbond, ' ', 2, ' ', 34 + atommer, ' ', 38 + atommer]])
        allbond = allbond + 1
        bondblock = ''.join([str(x) for x in [bondblock, lnovo, 'M  V30 ', 1 + allbond, ' ', 1, ' ', 35 + atommer, ' ', 36 + atommer]])
        allbond = allbond + 1
        bondblock = ''.join([str(x) for x in [bondblock, lnovo, 'M  V30 ', 1 + allbond, ' ', 2, ' ', 36 + atommer, ' ', 37 + atommer]])
        allbond = allbond + 1
        bondblock = ''.join([str(x) for x in [bondblock, lnovo, 'M  V30 ', 1 + allbond, ' ', 1, ' ', 36 + atommer, ' ', 41 + atommer]])
        allbond = allbond + 1
        bondblock = ''.join([str(x) for x in [bondblock, lnovo, 'M  V30 ', 1 + allbond, ' ', 1, ' ', 37 + atommer, ' ', 33 + atommer]])
        allbond = allbond + 1
      atommer = allatom
      bpairmono = bpairmono + 1
    atomblock = atomblock + lnovo + """M  V30 END ATOM
"""
    closer = """M  V30 END BOND
M  V30 END CTAB
M  END"""
    vitinfo = """
M  V30 BEGIN CTAB
M  V30 COUNTS """ + str(allatom) + " " + str(allbond) + " 0 0 1"
    dnamdlfile = ''.join([str(x2) for x2 in [dnamdlfile, ''.join([str(x) for x in [lnovo, ' ', 0, ' ', 0, '  0  0  0  0  0  0  0  0999 V3000']]), vitinfo, atomblock, bondblock, lnovo, closer]])
    naname = "dsdna" + str(dname) + ".mol"
    vmdl = open(naname, "w")
    vmdl.write(dnamdlfile)
    vmdl.close()
    txtid = naname
    convertformat()
    spename = "dsdna" + str(dname) + ".txt"





# Describe this function...
def vitamin():
  global randomizer, molid, command, worldfile, avatarx, avatary, avatarz, newline
  randomizer = random.randint(1, 120000)
  if randomizer <= 30:
    molid = 31
  elif randomizer <= 32:
    molid = 32
  elif randomizer <= 47:
    molid = 33
  elif randomizer <= 447:
    molid = 34
  elif randomizer <= 8447:
    molid = 35
  elif randomizer <= 8492:
    molid = 36
  elif randomizer <= 24000:
    molid = 37
  elif randomizer <= 24110:
    molid = 38
  elif randomizer <= 29110:
    molid = 39
  elif randomizer <= 30110:
    molid = 40
  elif randomizer <= 31000:
    molid = 41
  elif randomizer <= 32300:
    molid = 42
  elif randomizer <= 33400:
    molid = 43
  elif randomizer <= 47400:
    molid = 44
  else:
    molid = 30
  command = ''.join([str(x) for x in ['create_atoms 0 single ', avatarx, ' ', avatary, ' ', avatarz, ' mol ', molid, ' 3 rotate 0 1 0 0']])
  worldfile = ''.join([str(x2) for x2 in [worldfile, newline, command]])
  part_str()

# Describe this function...
def mineral():
  global randomizer, command, molid, worldfile, avatarx, avatary, avatarz, newline
  randomizer = random.randint(1, 2000000)
  if randomizer <= 35:
    molid = 51
  elif randomizer <= 45:
    molid = 52
  elif randomizer <= 195:
    molid = 53
  elif randomizer <= 1000195:
    molid = 54
  elif randomizer <= 1000197:
    molid = 55
  elif randomizer <= 1000242:
    molid = 56
  elif randomizer <= 1000253:
    molid = 59
  else:
    molid = 50
  command = ''.join([str(x) for x in ['create_atoms 0 single ', avatarx, ' ', avatary, ' ', avatarz, ' mol ', molid, ' 3 rotate 0 1 0 0']])
  worldfile = ''.join([str(x2) for x2 in [worldfile, newline, command]])
  part_str()


# Describe this function...
def nucgendna():
    global slime, genpure, finalform, command, worldfile, regnumr, originx, originy, originz, newline, blocksize, repeatx, repeaty, repeatz, avatarx, avatary, avatarz
    randomizer = random.randint(1, 60)
    if randomizer <= 20:
      finalform = random.choice(genpure)
    elif randomizer <= 40:
      finalform = slime
    else:
      dna1 = random.choice(genpure)
      dna2 = random.choice(genpure)
      seqs = {"seq1": dna1, "seq2": dna2}
      seqs = make_unaligned_seqs(data=seqs, moltype="dna")
      aln, tree = TreeAlign("HKY85", seqs, show_progress=False)
      newick = None
      newick = ''.join([str(x) for x in ['(seq1:', 1, ',seq2:', random.randint(1, 2000000000) / 1000000000, ');']])
      tree = make_tree(newick)
      sm = get_model("GN")
      lf = sm.make_likelihood_function(tree, digits=3, space=2)
      lf.set_alignment(aln)
      lf.optimise(show_progress=False)
      ancestors = lf.likely_ancestral_seqs()
      elder = str(ancestors)
      fnatext = elder
      lineloc = 1
      puregentext = ''
      textlength = len(fnatext)
      for count in range(int(textlength)):
        if fnatext[int(lineloc - 1)] == 'a' or fnatext[int(lineloc - 1)] == 'A':
          puregentext = str(puregentext) + 'A'
        elif fnatext[int(lineloc - 1)] == 'g' or fnatext[int(lineloc - 1)] == 'G':
            puregentext = str(puregentext) + 'G'
        elif fnatext[int(lineloc - 1)] == 'c' or fnatext[int(lineloc - 1)] == 'C':
            puregentext = str(puregentext) + 'C'
        elif fnatext[int(lineloc - 1)] == 't' or fnatext[int(lineloc - 1)] == 'T' or fnatext[int(lineloc - 1)] == 'u' or fnatext[int(lineloc - 1)] == 'U':
           puregentext = str(puregentext) + 'T'
        lineloc = lineloc + 1
      finalform = puregentext
    makestrand()

def protein():
  global molid, command, worldfile, avatarx, avatary, avatarz, newline
  molid = random.randint(200, 269)
  command = ''.join([str(x) for x in ['create_atoms 0 single ', avatarx, ' ', avatary, ' ', avatarz, ' mol ', molid, ' 3 rotate 0 1 0 0']])
  worldfile = ''.join([str(x2) for x2 in [worldfile, newline, command]])
  part_str()

# Describe this function...
def segmentation():
  global blockx, blocky, blockz, originx, originy, originz, repeatx, repeaty, repeatz, blocksize, avatary, avatarx, avatarz, molid, command, worldfile, newline, randomizer
  blockx = originx + 1000
  blocky = (originy + 970000000000) + 1000
  blockz = originz + 1000
  for count4 in range(int(math.floor(blocksize / 100000000000) - 10)):
    repeatx = math.floor(50000000000 / 100)
    repeaty = math.floor(1900000000 / 100)
    repeatz = math.floor(blocksize / 100) - 1000
    avatarx = blockx
    avatary = blocky
    avatarz = blockz
    for count3 in range(int(repeatx)):
      for count2 in range(int(repeaty)):
        for count in range(int(repeatz)):
          vitamin()
          avatarz = avatarz + 100
        avatarz = blockz
        avatary = avatary + 100
      avatary = blocky
      avatarx = avatarx + 100
    blockx = blockx + 100000000000
  blockx = originx + 1000
  blocky = (originy + 972000000000) + 1000
  blockz = originz + 1000
  for count8 in range(int(math.floor(blocksize / 70000000000) - 10)):
    repeatx = math.floor(50000000000 / 100)
    repeaty = math.floor(1900000000 / 100)
    repeatz = math.floor(blocksize / 100) - 1000
    avatarx = blockx
    avatary = blocky
    avatarz = blockz
    for count7 in range(int(repeatx)):
      for count6 in range(int(repeaty)):
        for count5 in range(int(repeatz)):
          molid = random.randint(60, 70)
          command = ''.join([str(x) for x in ['create_atoms 0 single ', avatarx, ' ', avatary, ' ', avatarz, ' mol ', molid, ' 3 rotate 0 1 0 0']])
          worldfile = ''.join([str(x2) for x2 in [worldfile, newline, command]])
          part_str()
          avatarz = avatarz + 100
        avatarz = blockz
        avatary = avatary + 100
      avatary = blocky
      avatarx = avatarx + 100
    blockx = blockx + 70000000000
  blockx = originx + 1000
  blocky = (originy + 974000000000) + 1000
  blockz = originz + 1000
  for count12 in range(int(math.floor(blocksize / 130000000000) - 10)):
    repeatx = math.floor(50000000000 / 100)
    repeaty = math.floor(1900000000 / 100)
    repeatz = math.floor(blocksize / 100) - 1000
    avatarx = blockx
    avatary = blocky
    avatarz = blockz
    for count11 in range(int(repeatx)):
      for count10 in range(int(repeaty)):
        for count9 in range(int(repeatz)):
          randomizer = random.randint(1, 10000000)
          if randomizer < 10000:
            molid = 100
          elif randomizer == 17001:
            molid = 101
          else:
            molid = random.randint(102, 121)
          command = ''.join([str(x3) for x3 in ['create_atoms 0 single ', avatarx, ' ', avatary, ' ', avatarz, ' mol ', molid, ' 3 rotate 0 1 0 0']])
          worldfile = ''.join([str(x4) for x4 in [worldfile, newline, command]])
          part_str()
          avatarz = avatarz + 100
        avatarz = blockz
        avatary = avatary + 100
      avatary = blocky
      avatarx = avatarx + 100
    blockx = blockx + 130000000000
  blockx = originx + 1000
  blocky = (originy + 976000000000) + 1000
  blockz = originz + 1000
  for count16 in range(int(math.floor(blocksize / 170000000000) - 10)):
    repeatx = math.floor(50000000000 / 100)
    repeaty = math.floor(1900000000 / 100)
    repeatz = math.floor(blocksize / 100) - 1000
    avatarx = blockx
    avatary = blocky
    avatarz = blockz
    for count15 in range(int(repeatx)):
      for count14 in range(int(repeaty)):
        for count13 in range(int(repeatz)):
          mineral()
          avatarz = avatarz + 100
        avatarz = blockz
        avatary = avatary + 100
      avatary = blocky
      avatarx = avatarx + 100
    blockx = blockx + 170000000000

def vitalparts():
  global command, worldfile, regnumr, originx, originy, originz, newline, blocksize, repeatx, repeaty, repeatz, avatarx, avatary, avatarz
  command = ''.join([str(x2) for x2 in ['region ', ''.join([str(x) for x in ['reg', regnumr, 'd']]), ' block ', originx + 1000, ' ', originx + (blocksize - 1000), ' ', originy + (960000000000 + 1000), ' ', originy + (970000000000 - 1000), ' ', originz + 1000, ' ', originz + (blocksize - 1000)]])
  worldfile = ''.join([str(x3) for x3 in [worldfile, newline, command]])
  worldfile = ''.join([str(x4) for x4 in [worldfile, newline, 'lattice sc 3.1061591792']])
  command = ''.join([str(x6) for x6 in ['create_atoms 0 region ', ''.join([str(x5) for x5 in ['reg', regnumr, 'd']]), ' mol 11 3 rotate 0 1 0 0']])
  worldfile = ''.join([str(x7) for x7 in [worldfile, newline, command]])
  command = ''.join([str(x9) for x9 in ['region ', ''.join([str(x8) for x8 in ['reg', regnumr, 'e']]), ' block ', originx + 1000, ' ', originx + (blocksize - 1000), ' ', originy + (980000000000 + 1000), ' ', originy + (990000000000 - 1000), ' ', originz + 1000, ' ', originz + (blocksize - 1000)]])
  worldfile = ''.join([str(x10) for x10 in [worldfile, newline, command]])
  worldfile = ''.join([str(x11) for x11 in [worldfile, newline, 'lattice sc 3.1061591792']])
  command = ''.join([str(x13) for x13 in ['create_atoms 0 region ', ''.join([str(x12) for x12 in ['reg', regnumr, 'e']]), ' mol 11 3 rotate 0 1 0 0']])
  worldfile = ''.join([str(x14) for x14 in [worldfile, newline, command]])
  part_str()
  repeatx = math.floor(((blocksize - 1000) - 2000) / 1000) - 5
  repeaty = math.floor((10000000000 - 2000) / 1000) - 5
  repeatz = math.floor(((blocksize - 1000) - 2000) / 1000) - 5
  avatarx = originx + 1000
  avatary = originy + (950000000000 + 1000)
  avatarz = originz + 1000
  for count3 in range(int(repeatx)):
    for count2 in range(int(repeaty)):
      for count in range(int(repeatz)):
        protein()
        avatarz = avatarz + 1000
      avatarz = originz + 1000
      avatary = avatary + 1000
    avatary = originy + (950000000000 + 1000)
    avatarx = avatarx + 1000
  repeatx = math.floor(((blocksize - 1000) - 2000) / 990000000000) - 5
  repeaty = math.floor((10000000000 - 2000) / 100) - 5
  repeatz = math.floor(((blocksize - 1000) - 2000) / 100) - 5
  avatarx = originx + 1000000000000
  avatary = originy + (990000000000 + 1000)
  avatarz = originz + 1000
  for count6 in range(int(repeatx)):
    for count5 in range(int(repeaty)):
      for count4 in range(int(repeatz)):
        worldfile = ''.join([str(x2) for x2 in [worldfile, newline, "dna ", avatarx, " ", avatary, " ", avatarz]])
        part_str()
        avatarz = avatarz + 100
      avatarz = originz + 1000
      avatary = avatary + 100
    avatary = originy + (990000000000 + 1000)
    avatarx = avatarx + 990000000000
  segmentation()


finalform = "ATCG"
blockx = 0
blocky = 0
blockz = 0
newline = """
"""
worldsize = 100000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
avatarx = 50000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000 + 100
avatary = 10000000000000 + 100
avatarz = 100
blocksize = 50000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000 - 200
repeatx = math.floor(blocksize / 20)
repeaty = math.floor((110000000000000 - 300) / 20)
repeatz = math.floor(blocksize / 20)
for count3 in range(int(repeatx)):
  for count2 in range(int(repeaty)):
    for count in range(int(repeatz)):
      avatarz = avatarz + 20
      randomizer = random.randint(1, 1000)
      molid = '1'
      if randomizer == 977:
        randomizer = random.randint(1, 1000)
        if randomizer == 977:
          randomizer = random.randint(1, 1000)
          if randomizer == 977:
            randomizer = random.randint(1, 1000)
            if randomizer == 977:
              randomizer = random.randint(1, 1000000)
              if randomizer == 977:
                randomizer = random.randint(1, 1000000)
                if randomizer == 977:
                  molid = '9'
                else:
                  molid = '8'
              else:
                molid = '7'
            else:
              molid = '6'
          else:
            molid = '5'
        else:
          molid = '4'
      elif randomizer < 850:
        molid = '1'
      elif randomizer < 900:
        molid = '3'
      else:
        molid = '2'
      command = ''.join([str(x) for x in ['create_atoms 0 single ', avatarx, ' ', avatary, ' ', avatarz, ' mol ', molid, ' 3 rotate 0 1 0 0']])
      worldfile = ''.join([str(x2) for x2 in [worldfile, newline, command]])
      part_str()
    avatary = avatary + 20
    avatarz = 100
  avatary = 10000000000000 + 100
  avatarx = avatarx + 20



originx = 100
originy = 10000000000000 + 100
originz = 100
regnumr = 1
for count in range(110):
  command = ''.join([str(x2) for x2 in ['region ', ''.join([str(x) for x in ['reg', regnumr, 'a']]), ' block ', originx + 1000, ' ', originx + (blocksize - 1000), ' ', originy + 1000, ' ', originy + (8300000000 - 1000), ' ', originz + 1000, ' ', originz + (blocksize - 1000)]])
  worldfile = ''.join([str(x3) for x3 in [worldfile, newline, command]])
  worldfile = ''.join([str(x4) for x4 in [worldfile, newline, 'lattice fcc 6.38']])
  command = ''.join([str(x6) for x6 in ['create_atoms 0 region ', ''.join([str(x5) for x5 in ['reg', regnumr, 'a']]), ' mol 58 3 rotate 0 1 0 0']])
  worldfile = ''.join([str(x3) for x3 in [worldfile, newline, command]])
  command = ''.join([str(x8) for x8 in ['region ', ''.join([str(x7) for x7 in ['reg', regnumr, 'b']]), ' block ', originx + 1000, ' ', originx + (blocksize - 1000), ' ', originy + (8300000000 + 1000), ' ', originy + (9000000000 - 1000), ' ', originz + 1000, ' ', originz + (blocksize - 1000)]])
  worldfile = ''.join([str(x9) for x9 in [worldfile, newline, command]])
  worldfile = ''.join([str(x10) for x10 in [worldfile, newline, 'lattice fcc 5.63']])
  command = ''.join([str(x12) for x12 in ['create_atoms 0 region ', ''.join([str(x11) for x11 in ['reg', regnumr, 'b']]), ' mol 57 3 rotate 0 1 0 0']])
  worldfile = ''.join([str(x3) for x3 in [worldfile, newline, command]])
  command = ''.join([str(x14) for x14 in ['region ', ''.join([str(x13) for x13 in ['reg', regnumr, 'c']]), ' block ', originx + 1000, ' ', originx + (blocksize - 1000), ' ', originy + (9000000000 + 1000), ' ', originy + (950000000000 - 1000), ' ', originz + 1000, ' ', originz + (blocksize - 1000)]])
  worldfile = ''.join([str(x15) for x15 in [worldfile, newline, command]])
  worldfile = ''.join([str(x16) for x16 in [worldfile, newline, 'lattice sc 3.1061591792']])
  command = ''.join([str(x18) for x18 in ['create_atoms 0 region ', ''.join([str(x17) for x17 in ['reg', regnumr, 'c']]), ' mol 11 3 rotate 0 1 0 0']])
  worldfile = ''.join([str(x3) for x3 in [worldfile, newline, command]])
  part_str()
  vitalparts()
  regnumr = regnumr + 1
  originy = originy + 1000000000000








command = """mass 1 1.008
mass 2 4.003
mass 3 6.941
mass 4 9.012
mass 5 10.811
mass 6 12.011
mass 7 14.007
mass 8 15.999
mass 9 18.998
mass 10 20.180
mass 11 22.990
mass 12 24.305
mass 13 26.982
mass 14 28.086
mass 15 30.974
mass 16 32.065
mass 17 35.453
mass 18 39.948
mass 19 39.098
mass 20 40.078
mass 21 44.956
mass 22 47.867
mass 23 50.942
mass 24 51.996
mass 25 54.938
mass 26 55.845
mass 27 58.933
mass 28 58.693
mass 29 63.546
mass 30 65.390
mass 31 69.723
mass 32 72.640
mass 33 74.922
mass 34 78.960
mass 35 79.904
mass 36 83.800
mass 37 85.468
mass 38 87.620
mass 39 88.906
mass 40 91.224
mass 41 92.906
mass 42 95.940
mass 43 98.000
mass 44 101.070
mass 45 102.906
mass 46 106.420
mass 47 107.868
mass 48 112.411
mass 49 114.818
mass 50 118.710
mass 51 121.760
mass 52 127.600
mass 53 126.905
mass 54 131.293
mass 55 132.906
mass 56 137.327
mass 57 138.906
mass 58 140.116
mass 59 140.908
mass 60 144.240
mass 61 145.000
mass 62 150.360
mass 63 151.964
mass 64 157.250
mass 65 158.925
mass 66 162.500
mass 67 164.930
mass 68 167.259
mass 69 168.934
mass 70 173.040
mass 71 174.967
mass 72 178.490
mass 73 180.948
mass 74 183.840
mass 75 186.207
mass 76 190.230
mass 77 192.217
mass 78 195.078
mass 79 196.967
mass 80 200.590
mass 81 204.383
mass 82 207.200
mass 83 208.980
mass 84 209.000
mass 85 210.000
mass 86 222.000
mass 87 223.000
mass 88 226.000
mass 89 227.000
mass 90 232.038
mass 91 231.036
mass 92 238.029
mass 93 237.000
mass 94 244.000
mass 95 243.000
mass 96 247.000
mass 97 247.000
mass 98 251.000
mass 99 252.000
mass 100 257.000
fix ensemble all nvt temp 308.15 308.15 100 tchain 1
fix xwalls all wall/reflect xlo 10 xhi 1000000000000000900 ylo 10 yhi 109990000000000000 zlo 10 zhi 1000000000000000900 units box
fix 1 all gravity 0.000000000000000234230942 vector 0 -1 0
fix reactions all acks2/reaxff 1 0.0 10.0 1.0e-6 reaxff
delete_atoms overlap 0.3 all all bond yes

timestep 1.0

# Output
thermo_style one
thermo 0

# Run the simulation
run 1000000000
"""

worldfile = ''.join([str(x) for x in [worldfile, newline, command]])

command = ''.join([str(x) for x in ['write_data era', era, '.txt nocoeff nofix']])

worldfile = ''.join([str(x) for x in [worldfile, newline, command]])


worldfile_list = [worldfile, worldfile_list]

from decimal import *



def check_h(jb, hb):
    if (((Decimal(jb[2]) - Decimal(hb[2])) ** Decimal(2)) + ((Decimal(jb[3]) - Decimal(hb[3])) ** Decimal(2)) + ((Decimal(jb[4]) - Decimal(hb[4])) ** Decimal(2))) ** Decimal(0.5) < Decimal(6):
        if jb[0] != hb[0]:
            return True
        else: return False
    else: return False





def h_find(ja):
    global atom_list, list_neighbors
    browselist2 = True
    list_sel2 = atom_list
    while browselist2 == True:
        if len(list_sel2[-1]) == 1:
            browselist2 = False
        cur_l = list_sel2[0]
        for h in cur_l:
            if check_h(ja, h) == True and len(list_neighbors) < 10000:
                list_neighbors.append(h)
        if len(list_sel2[-1]) != 1:
            list_sel2 = list_sel2[-1]





def run_lammps():
    global key_args, list_molecules, list_regions, lattice_size, lattice_type, posx, posy, posz, molpath, all_particles, world_structure, list_j, list_neighbors, atom_list
    for j in list_j:
        gravityvalue = random.randint(1, 200)
        if gravityvalue == 177:
            j[8] = str(Decimal(j[8]) + Decimal("-0.000000000000000019600445"))
        if Decimal(j[7]) > Decimal(499):
            j[7] = "499"
        if Decimal(j[7]) < Decimal(-499):
            j[7] = "-499"
        if Decimal(j[8]) > Decimal(499):
            j[8] = "499"
        if Decimal(j[8]) < Decimal(-499):
            j[8] = "-499"
        if Decimal(j[9]) > Decimal(499):
            j[9] = "499"
        if Decimal(j[9]) < Decimal(-499):
            j[9] = "-499"
        if Decimal(j[2]) < Decimal(0):
            j[2] = str(Decimal(j[2]) * Decimal(-1))
            j[7] = str(Decimal(j[7]) * Decimal(-1))
        if Decimal(j[2]) > Decimal(worldsize):
            j[2] = str(Decimal(worldsize) - (Decimal(j[2]) - Decimal(worldsize)))
            j[7] = str(Decimal(j[7]) * Decimal(-1))
        if Decimal(j[3]) < Decimal(0):
            j[3] = str(Decimal(j[3]) * Decimal(-1))
            j[8] = str(Decimal(j[8]) * Decimal(-1))
        if Decimal(j[3]) > Decimal("110000000000000000"):
            j[3] = str(Decimal("110000000000000000") - (Decimal(j[3]) - Decimal("110000000000000000")))
            j[8] = str(Decimal(j[8]) * Decimal(-1))
        if Decimal(j[4]) < Decimal(0):
            j[4] = str(Decimal(j[4]) * Decimal(-1))
            j[9] = str(Decimal(j[9]) * Decimal(-1))
        if Decimal(j[4]) > Decimal(worldsize):
            j[4] = str(Decimal(worldsize) - (Decimal(j[4]) - Decimal(worldsize)))
            j[9] = str(Decimal(j[9]) * Decimal(-1))
        list_neighbors = []
        h_find(j)
        list_neighbors = list_neighbors[0:10000]
        world_text = """LAMMPS data file via write_data, version 22 Dec 2022, timestep = 1000000

"""
        world_text = world_text + str(len(list_neighbors) + 1)
        text_input = """ atoms
100 atom types

-500 500 xlo xhi
-500 500 ylo yhi
-500 500 zlo zhi

Masses

1 1.008
2 4.003
3 6.941
4 9.012
5 10.811
6 12.011
7 14.007
8 15.999
9 18.998
10 20.18
11 22.99
12 24.305
13 26.982
14 28.086
15 30.974
16 32.065
17 35.453
18 39.948
19 39.098
20 40.078
21 44.956
22 47.867
23 50.942
24 51.996
25 54.938
26 55.845
27 58.933
28 58.693
29 63.546
30 65.39
31 69.723
32 72.64
33 74.922
34 78.95999999999999
35 79.904
36 83.8
37 85.468
38 87.62
39 88.90600000000001
40 91.224
41 92.90600000000001
42 95.94
43 98
44 101.07
45 102.906
46 106.42
47 107.868
48 112.411
49 114.818
50 118.71
51 121.76
52 127.6
53 126.905
54 131.293
55 132.906
56 137.327
57 138.906
58 140.116
59 140.908
60 144.24
61 145
62 150.36
63 151.964
64 157.25
65 158.925
66 162.5
67 164.93
68 167.259
69 168.934
70 173.04
71 174.967
72 178.49
73 180.948
74 183.84
75 186.207
76 190.23
77 192.217
78 195.078
79 196.967
80 200.59
81 204.383
82 207.2
83 208.98
84 209
85 210
86 222
87 223
88 226
89 227
90 232.038
91 231.036
92 238.029
93 237
94 244
95 243
96 247
97 247
98 251
99 252
100 257

Atoms # full

"""
        world_text = world_text + text_input
        lneo = """
"""
        atom_count = 1
        list_v = []
        list_a = []
        txtl_a = []
        txtl_v = []
        list_v.append([str(atom_count), j[7], j[8], j[9]])
        list_a.append([str(atom_count), "0", j[1], j[5], "0", "0", "0"])
        for hj in list_neighbors:
            atom_count = atom_count + 1
            rel_x = str(Decimal(hj[2]) - Decimal(j[2]))
            rel_y = str(Decimal(hj[3]) - Decimal(j[3]))
            rel_z = str(Decimal(hj[4]) - Decimal(j[4]))
            list_v.append([str(atom_count), hj[7], hj[8], hj[9]])
            list_a.append([str(atom_count), "0", hj[1], hj[5], rel_x, rel_y, rel_z])
        for xy in list_a:
            txtl_a.append(' '.join([str(elem) for elem in xy]))
        for zy in list_v:
            txtl_v.append(' '.join([str(elem) for elem in zy]))
        txt_a = lneo.join([str(elem) for elem in txtl_a])
        txt_v = lneo.join([str(elem) for elem in txtl_v])
        world_text = world_text + txt_a
        text_input = """

Velocities

"""
        world_text = world_text + text_input
        world_text = world_text + txt_v
        j_file = open("iworld.txt", "w")
        j_file.write(world_text)
        j_file.close()
        for ch in world_text:
            unused_var = 1
        os.system("chmod u+x /workspace/sudoku/lmp")
        os.system("#!/bin/sh")
        os.system("/workspace/sudoku/lmp -in u9.lmp")
        l_file = open("i1.txt", "r")
        l_txt = l_file.read()
        l_lmp = l_txt.splitlines()
        if "Velocities" in l_lmp and "Atoms # full" in l_lmp:
            v_idx = l_lmp.index("Velocities")
            a_idx = l_lmp.index("Atoms # full")
            l_vel = l_lmp[v_idx:]
            l_ats = l_lmp[a_idx:v_idx]
            for ai in l_ats:
                a2 = ai.split(" ")
                if a2[0] == "1":
                    j[2] = str(Decimal(j[2]) + Decimal(a2[4]))
                    j[3] = str(Decimal(j[3]) + Decimal(a2[5]))
                    j[4] = str(Decimal(j[4]) + Decimal(a2[6]))
                    j[5] = a2[3]
            for vi in l_vel:
                v2 = vi.split(" ")
                if v2[0] == "1":
                    j[7] = v2[1]
                    j[8] = v2[2]
                    j[9] = v2[3]

                


def onestep():
    global key_args, list_molecules, list_regions, lattice_size, lattice_type, posx, posy, posz, molpath, all_particles, world_structure, list_j, list_neighbors, atom_list
    llevl = 1
    lidx = 0
    browselist = True
    list_sel = atom_list
    while browselist == True:
        if len(list_sel[-1]) == 1:
            browselist = False
        list_j = list_sel[0]
        run_lammps()
        if len(list_sel[-1]) != 1:
            list_sel = list_sel[-1]
            llevl = llevl + 1





def monomer():
    global key_args, list_molecules, list_regions, lattice_size, lattice_type, posx, posy, posz, molpath, all_particles, world_structure, list_j
    selfile = open(molpath, "r")
    return_data = selfile.read()
    list_atoms = return_data.splitlines()
    size_molecular = int(list_atoms[2].split(" ")[0])
    list_masses = list_atoms[(list_atoms.index("Masses")):(list_atoms.index("Masses") + size_molecular + 2)]
    list_charges = list_atoms[(list_atoms.index("Charges")):(list_atoms.index("Charges") + size_molecular + 2)]
    list_coords = list_atoms[(list_atoms.index("Coords")):(list_atoms.index("Coords") + size_molecular + 2)]
    list_types = list_atoms[(list_atoms.index("Types")):(list_atoms.index("Types") + size_molecular + 2)]
    for mycount in range(int(size_molecular)):
        current_atom = []
        current_atom.append(all_particles)
        all_particles = all_particles + 1
        current_atom.append(list_types[mycount + 2].split(" ")[3])
        current_atom.append(list_coords[mycount + 2].split(" ")[3])
        current_atom.append(list_coords[mycount + 2].split(" ")[4])
        current_atom.append(list_coords[mycount + 2].split(" ")[5])
        current_atom.append(list_charges[mycount + 2].split(" ")[3])
        current_atom.append(list_masses[mycount + 2].split(" ")[3])
        current_atom[2] = str(Decimal(current_atom[2]) + posx)
        current_atom[3] = str(Decimal(current_atom[3]) + posy)
        current_atom[4] = str(Decimal(current_atom[4]) + posz)
        current_atom.append("0")
        current_atom.append("0")
        current_atom.append("0")
        list_j.append(current_atom)
        part_list()







def mol_single():
    global key_args, list_molecules, list_regions, lattice_size, lattice_type, posx, posy, posz, molpath, all_particles, world_structure, list_j
    for x in list_molecules:
        if x[0] == key_args[7]:
            molpath = x[1]
    posx = Decimal(key_args[3])
    posy = Decimal(key_args[4])
    posz = Decimal(key_args[5])
    monomer()



def mol_lattice():
    global key_args, list_molecules, list_regions, lattice_size, lattice_type, posx, posy, posz, molpath, all_particles, world_structure, list_j
    for x in list_molecules:
        if x[0] == key_args[5]:
            molpath = x[1]
    for r in range(int(len(list_regions))):
        if (list_regions[r])[0] == key_args[3]:
            selreg = list_regions[r]
    repeat_x = math.floor((selreg[2] - selreg[1]) / lattice_size) - 1
    repeat_y = math.floor((selreg[4] - selreg[3]) / lattice_size) - 1
    repeat_z = math.floor((selreg[6] - selreg[5]) / lattice_size) - 1
    if lattice_type == 1:
        posx = selreg[1]
        posy = selreg[3]
        posz = selreg[5]
        for Count1 in range(int(repeat_x)):
            for Count2 in range(int(repeat_y)):
                for Count3 in range(int(repeat_z)):
                    posz = posz + lattice_size
                    monomer()
                posz = selreg[5]
                posy = posy + lattice_size
            posy = selreg[3]
            posx = posx + lattice_size
    if lattice_type == 2:
        posx = selreg[1]
        posy = selreg[3]
        posz = selreg[5]
        for Count1 in range(int(repeat_x)):
            for Count2 in range(int(repeat_y)):
                for Count3 in range(int(repeat_z)):
                    posz = posz + lattice_size
                    monomer()
                posz = selreg[5]
                posy = posy + lattice_size
            posy = selreg[3]
            posx = posx + lattice_size
        posx = selreg[1] + lattice_size / 2
        posy = selreg[3] + lattice_size / 2
        posz = selreg[5] + lattice_size / 2
        for Count1 in range(int(repeat_x)):
            for Count2 in range(int(repeat_y)):
                for Count3 in range(int(repeat_z)):
                    posz = posz + lattice_size
                    monomer()
                posz = selreg[5] + lattice_size / 2
                posy = posy + lattice_size
            posy = selreg[3] + lattice_size / 2
            posx = posx + lattice_size
    if lattice_type == 3:
        posx = selreg[1]
        posy = selreg[3]
        posz = selreg[5]
        for Count1 in range(int(repeat_x)):
            for Count2 in range(int(repeat_y)):
                for Count3 in range(int(repeat_z)):
                    posz = posz + lattice_size
                    monomer()
                posz = selreg[5]
                posy = posy + lattice_size
            posy = selreg[3]
            posx = posx + lattice_size
        posx = selreg[1] + lattice_size / 2
        posy = selreg[3]
        posz = selreg[5] + lattice_size / 2
        for Count1 in range(int(repeat_x)):
            for Count2 in range(int(repeat_y)):
                for Count3 in range(int(repeat_z)):
                    posz = posz + lattice_size
                    monomer()
                posz = selreg[5] + lattice_size / 2
                posy = posy + lattice_size
            posy = selreg[3]
            posx = posx + lattice_size
        posx = selreg[1] + lattice_size / 2
        posy = selreg[3] + lattice_size / 2
        posz = selreg[5]
        for Count1 in range(int(repeat_x)):
            for Count2 in range(int(repeat_y)):
                for Count3 in range(int(repeat_z)):
                    posz = posz + lattice_size
                    monomer()
                posz = selreg[5]
                posy = posy + lattice_size
            posy = selreg[3] + lattice_size / 2
            posx = posx + lattice_size
        posx = selreg[1]
        posy = selreg[3] + lattice_size / 2
        posz = selreg[5] + lattice_size / 2
        for Count1 in range(int(repeat_x)):
            for Count2 in range(int(repeat_y)):
                for Count3 in range(int(repeat_z)):
                    posz = posz + lattice_size
                    monomer()
                posz = selreg[5] + lattice_size / 2
                posy = posy + lattice_size
            posy = selreg[3] + lattice_size / 2
            posx = posx + lattice_size


molpath = ""
posx = 0
posy = 0
posz = 0

list_j = []
all_particles = 1
world_structure = []
list_molecules = []
list_regions = []
lattice_size = 1
lattice_type = 1

def read_world():
    global worldfile, key_args, lammps_input, lattice_size, lattice_type, list_regions, list_molecules, posx, posy, posz, molpath
    lammps_input = worldfile.splitlines()
    for i in lammps_input:
        key_args = i.split(" ")
        if key_args[0] == "create_atoms":
            if key_args[2] == "single":
                mol_single()
            if key_args[2] == "region":
                mol_lattice()
        if key_args[0] == "lattice":
            lattice_size = Decimal(key_args[2])
            if key_args[1] == "sc":
                lattice_type = 1
            if key_args[1] == "bcc":
                lattice_type = 2
            if key_args[1] == "fcc":
                lattice_type = 3
        if key_args[0] == "region":
            list_regions.append([key_args[1], Decimal(key_args[3]), Decimal(key_args[4]), Decimal(key_args[5]), Decimal(key_args[6]), Decimal(key_args[7]), Decimal(key_args[8])])
        if key_args[0] == "molecule":
            list_molecules.append([key_args[1], key_args[2]])
        if key_args[0] == "dna":
            nucgendna()
            posx = Decimal(key_args[1])
            posy = Decimal(key_args[2])
            posz = Decimal(key_args[3])
            molpath = "dsdna" + str(dname) + ".txt"
            monomer()



browselist3 = True
llevl3 = 0
list_sel3 = worldfile_list
while browselist3 == True:
    if type(list_sel3[-1]) != type([]):
        browselist3 = False
    if type(list_sel3[-1]) == type([]):
        list_sel3 = list_sel3[-1]
        llevl3 = llevl3 + 1

llevl3 = llevl3 - 1

while llevl3 >= 0:
    list_sel3 = worldfile_list
    for ll_c3 in range(int(llevl3)):
        if type(list_sel3[-1]) == type([]):
            list_sel3 = list_sel3[-1]
    worldfile = list_sel3[0]
    read_world()
    llevl3 = llevl3 - 1


atom_list = [list_j, atom_list]

split_line = """
"""
while True:
    onestep()






