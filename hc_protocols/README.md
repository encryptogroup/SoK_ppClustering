# Hierarchical Clustering

This repository contains an implementation of the PCA and OPT algorithms presented by Meng et al. in [MPO19](https://arxiv.org/abs/1904.04475).

## External Dependencies
The following libraries need to be installed separately and should be available to the build system and compiler.

- [ABY](https://github.com/encryptogroup/ABY)
- [GMP](https://gmplib.org/)
- [Boost](https://www.boost.org/) (1.72.0 or later)
- [LibPaillier](http://hms.isi.jhu.edu/acsc/libpaillier/)
- [RELIC](https://github.com/relic-toolkit/relic)
- [CNPY](https://github.com/rogersce/cnpy)
- [Nlohmann JSON](https://github.com/nlohmann/json)

## Compilation
The project uses [CMake](https://cmake.org/) for building the source code.
To compile, run the following from the `hc_protocols` directory:

```sh
mkdir build && cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make
```

## Usage
The `hc_benchmark` program can be used to benchmark the performance of the PCA and OPT protocols.

Execute the following commands from the `build` directory created during compilation to run the programs on the sample dataset:
```sh
# Run unit tests. Can be skipped.
ctest

# Benchmark performance of PCA on input dataset.
# The command below runs both parties with the first party running in the background. '--party2' is used to indicate the second party.
# Parties can also be run on separate machines.
# Outputs files `party1.json` and `party2.json` containing the benchmark data corresponding to each party.
./hc_benchmark --output ./party1.json --dataset ../../data/dataset.npy --clusters 2 --protocol pca & \
  ./hc_benchmark --output ./party2.json --dataset ../../data/dataset.npy --clusters 2 --protocol pca --party2

# Changing the argument to '--protocol' to opt will run the OPT protocol.
```

## Directory Structure
A summary of the source code structure is given below:

- `src`: This directory contains the implementation of the PCA and OPT protocols.
  It also provides a simple Python wrapper for easy interop between the Python and the C++ implementation (check description of `python` directory below).
- `include/mpo19`: This directory contains the header files for the implementation in `src`.
- `tests`: This directory contains extensive unit tests for the code in `src` and `include/mpo19`.
- `python`: This directory contains an example on utilizing the python interop provided by the core library. 
  This python source is utilized in `src/mpo19.cpp`.
- `benchmark`: This directory contains the programs used for benchmarking the protocols.
  It essentially wraps the implementation provided by `src` into a program with an easy to use interface.
- `cmake`: This directory consists of CMake configuration files for using external libraries.
