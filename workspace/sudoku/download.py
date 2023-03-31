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



import os
os.system("git clone -b release https://github.com/lammps/lammps.git mylammps")

xfile = open("a1.txt", "w")
xfile.write("placeholder")
xfile.close
xfile = open("b1.txt", "w")
xfile.write("placeholder")
xfile.close
xfile = open("c1.txt", "w")
xfile.write("placeholder")
xfile.close
xfile = open("d1.txt", "w")
xfile.write("placeholder")
xfile.close
xfile = open("e1.txt", "w")
xfile.write("placeholder")
xfile.close
xfile = open("f1.txt", "w")
xfile.write("placeholder")
xfile.close
xfile = open("g1.txt", "w")
xfile.write("placeholder")
xfile.close
xfile = open("h1.txt", "w")
xfile.write("placeholder")
xfile.close
xfile = open("i1.txt", "w")
xfile.write("placeholder")
xfile.close
xfile = open("j1.txt", "w")
xfile.write("placeholder")
xfile.close

xfile = open("aworld.txt", "w")
xfile.write("placeholder")
xfile.close
xfile = open("bworld.txt", "w")
xfile.write("placeholder")
xfile.close
xfile = open("cworld.txt", "w")
xfile.write("placeholder")
xfile.close
xfile = open("dworld.txt", "w")
xfile.write("placeholder")
xfile.close
xfile = open("eworld.txt", "w")
xfile.write("placeholder")
xfile.close
xfile = open("fworld.txt", "w")
xfile.write("placeholder")
xfile.close
xfile = open("gworld.txt", "w")
xfile.write("placeholder")
xfile.close
xfile = open("hworld.txt", "w")
xfile.write("placeholder")
xfile.close
xfile = open("iworld.txt", "w")
xfile.write("placeholder")
xfile.close
xfile = open("jworld.txt", "w")
xfile.write("placeholder")
xfile.close



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
os.system("pip install rdkit")
import rdkit
from rdkit import Chem
from os import listdir