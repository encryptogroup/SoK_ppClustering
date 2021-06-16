# HE-Meanshift

This subdirectory contains an implementation of the HE-Meanshift clustering algorithm presented by Cheon et al. in [CKP19](https://eprint.iacr.org/2019/465).
We use the [HEAAN](https://github.com/snucrypto/HEAAN) library for homomorphic encryption.
The source code of the HEAAN library has been copied into the `extern/HEAAN` directory for ease of compilation and usage.

## External Dependencies
The following libraries need to be installed separately and should be available to the build system and compiler.

- [GMP](https://gmplib.org/)
- [NTL](https://www.shoup.net/ntl/) (11.0.0 or later)
- [Boost](https://www.boost.org/) (1.72.0 or later)
- [CNPY](https://github.com/rogersce/cnpy)
- [Nlohmann JSON](https://github.com/nlohmann/json)

## Compilation
The project uses [CMake](https://cmake.org/) for building the source code.
To compile, run the following from the `he_meanshift` directory:

```sh
mkdir build && cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make
```

## Usage
A short description of the compiled programs is given below.
All of them provide detailed usage description on using the `--help` option.

- `benchmarks/cluster`: Run and benchmark the HE-Meanshift clustering algorithm.
- `benchmarks/plaintext_cluster`: Run the plaintext variant of the protocol. Refer to [Directory Structure](#directory-structure) for more details.
- `tests/ckp19_test`: Run unit tests.

Execute the following commands from the `build` directory created during compilation to run the programs on the sample dataset:

```sh
# Run unit tests. Can be skipped.
./tests/ckp19_test --show-progress

# Run the HE-Meanshift protocol to cluster the sample dataset.
# Configuration options are given via sample_config.ini but can also be given through command line arguments
# Upon completion, stores the output in the sec_output directory. 
# sec_output/info.json contains benchmark data and sec_output/rep_0.npy contains the clustering output.
./benchmarks/cluster --output sec_output --dusts 2 --dataset ../../data/dataset.npy -c ../sample_config.ini

# Evaluate clustering output on external and internal indices.
python3 ../../utils/he_meanshift_evaluate.py ../../data dataset sec_output
```

**NOTE**: Generating the bootstrap keys might require significant time and memory with the default (recommended) HEAAN parameters.
Temporarily setting the `logN` parameter to a lower value, say `logN=10`, in `he_meanshift/extern/HEAAN/src/Params.h` might serve better for testing purposes.

## Directory Structure
A summary of the source code structure is given below:

- `src/ckp19`: This directory contains the implementation of the HE-Meanshift protocol which can be accessed through an instance of the `MeanShift` class.
  It also contains a plaintext implementation of the HE-Meanshift protocol.
  The main difference between the plaintext and secure variant is that the former computes directly on plaintext without using any secure computation.
  The output of the plaintext variant will be similar to the secure variant within an error bound (since HEAAN is an approximate HE scheme).
- `test`: This directory contains unit tests for the code in `src/ckp19`.
- `benchmark`: This directory contains the source for benchmarking the protocol implementation provided by `src/ckp19`.
  It essentially wraps the implementation provided by `src/ckp19` into individual programs with an easy to use interface.
- `extern`: This directory contains a clone of the HEAAN library for ease of compilation.
- `cmake`: This directory contains CMake configurations for using external libraries.
