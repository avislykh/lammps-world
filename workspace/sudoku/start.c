// Copyright 2023 Alexander Vislykh
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.



#include <stdlib.h>
#include <stdio.h>
#include <string.h>

extern const unsigned char _binary_ubuntu_ext4_start[];
extern const unsigned char _binary_ubuntu_ext4_end[];
extern const size_t *_binary_ubuntu_ext4_size;

int main () {
   char command[50];

   strcpy( command, "python3 /workspace/sudoku/startup.py" );
   system(command);

   return(0);
} 