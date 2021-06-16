# SoK: Efficient Privacy-preserving Clustering

This repository contains the source code for the PETS'21 paper [HMS+21] by Aditya Hegde, Helen Möllering, Thomas Schneider, and Hossein Yalame.

## Components

A brief description of the subdirectories in the codebase is given below. The README in each subdirectory provides more information on compilation and usage.

- `he_meanshift`: An implementation of the HE-Meanshift protocol presented by Cheon et al. in [CKP19](https://eprint.iacr.org/2019/465).
- `hc_protocols`: An implementation of the hierarchical clustering protocols of Meng et al. in [MPO19](https://arxiv.org/abs/1904.04475).
- `utils`: Scripts to automate simple tasks and aid in analysis.
- `data`: A sample dataset to use as input. See the [Datasets](#datasets) section for more details.

## Building the Project

### Docker

All required dependencies to compile and run the project are available through the docker image.
To use docker run the following from the root of the repository:

```sh
docker build -t sokppcluster .
docker run -t -i -v $PWD:/code sokppcluster

# The command below should be run inside the container.
cd /code
```

## Datasets

The datasets we use for evaluating clustering quality are available at the public GitHub repository [gagolews/clustering\_benchmarks\_v1](https://github.com/gagolews/clustering_benchmarks_v1).
While the above repository provides datasets in text format saved as `.gz` files, the C++ benchmark programs require the input dataset to be in Numpy's `.npy` format.
The `utils/transform_data.py` program can be used to convert the `.gz` file into `.npy` format.
Please refer to the README in the `utils` directory for usage information.

A sample dataset in the above formats along with the corresponding ground truth as created using [Sci-kit learn's](https://scikit-learn.org/stable/modules/generated/sklearn.datasets.make_blobs.html) `make_blobs` function is available in the `data` directory.
It consists of 128 data records each having 1 attribute and consists of 2 clusters.
