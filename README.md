# DumpyOS: A Data-Adaptive Multi-ary Index for Scalable Data Series Similarity Search (VLDB Journal in 2024)

DumpyOS is an enhanced version of Dumpy, focusing on building efficiency, search accuracy, and pruning-based search efficiency.

# License
This archive is free for use for academic and non-profit purposes, but if you use it, please reference it properly.

# Reference
Wang, Z., Wang, Q., Wang, P. et al. DumpyOS: A data-adaptive multi-ary index for scalable data series similarity search. The VLDB Journal (2024). https://doi.org/10.1007/s00778-024-00874-9

# Disclaimer
The code is provided without warranty of any kind. While we thoroughly tested all code bases on Ubuntu 20.04 LTS (Windows Subsystem of Linux), we do not guarantee that they are exempt from bugs, nor that they will work on other platforms. If you encounter any issues with the code, please feel free to propose them on the ISSUE page of this repo. We will do our best to address your concerns but do not promise to resolve all issues.

# Installation

## Libraries

2 libraries are required: GSL and boost:serialization.
GSL can be omitted if you don't need to generate random walk dataset on your own. Just remove related line in CMakeLists.txt, with the file name "src/Expr/RandDataGenerator.cpp" "include/Expr/RandDataGenerator.h" in COMMAND add_executable (10th line in CMakeLists.txt)
boost:serialization is used to serialize the memory index structure, so it is necessary for NOW, though an ad-hoc serialization version might be developed in the future.

## Build

Here is a Cmake project. All the information about compiling is in CMakeList.txt.
We test this project in Cmake 3.20 with C++ language standard 23, however, we SUBJECTIVELY think it can run on relatively lower versions.

*If you can use IDE like Clion*, just open this project use that and build it automatically and then run it.

*Else*, please refer to the following instructions.

1. create a "build" directory under the project

2. cd build

3. cmake ..

4. make

5. cd .. && ./bin/FADAS

## Run

All the configuration is written on the config.ini, including the information about dataset.
To finish your task, please READ and UPDATE config.ini IN DETAIL.

