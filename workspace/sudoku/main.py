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
import dimod
from dimod import ConstrainedQuadraticModel, Binary
from dwave.system import DWaveSampler, EmbeddingComposite ### imports
from dwave.system import LeapHybridDQMSampler
import urllib
import urllib3
from urllib import request
import socket
import socketserver
import requests
import os

############ create image of the system

os.system("sudo add-apt-repository ppa:brandonsnider/cdrtools")
os.system("sudo apt-get update")
os.system("sudo apt-get install cdda2wav cdrecord mkisofs")
os.system("sudo mkisofs -D -o ubuntu.iso /")

######### convert to current linux format filesystem

os.system("dd if=/workspace/sudoku/ubuntu.iso of=/workspace/sudoku/ubuntu.ext4")

######## create linkable binary file for the ubuntu filesystem
os.system("ld --relocatable --format=binary --output=ubuntu.o ubuntu.ext4")
os.system("nm ubuntu.o")

#################
os.system("gcc start.c  -o begin.o") # create a .o file for the program innitiator

# link files into executable code

os.system("ld --output=final begin.o ubuntu.o --entry=0000000000401000 --no-gc-sections")

############ upload to quantum computer


solver1 = DWaveSampler(solver={'topology__type': 'pegasus'})
x = solver1.client
x.permissive_ssl = True
ses = x.create_session()
fileq = open("/workspace/sudoku/final", "rb")
y = x.upload_problem_encoded(fileq.read())

while True:
    var_unused = 1


# Due to the code containing references to D-Wave,
#The following copyright notice is included just in case:



# Copyright 2019 D-Wave Systems Inc.
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