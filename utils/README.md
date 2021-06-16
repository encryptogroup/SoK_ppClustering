# Utils

This directory contains a few helpful Python scripts for performing a number of tasks.

## External Dependencies
In addition to a Python3 installation, the following packages are required for running the scripts.
All of them can be installed using `pip`.

- [Scikit-learn](https://scikit-learn.org/stable/)
- [Numpy](https://numpy.org/)
- [kneed](https://pypi.org/project/kneed/)

## Usage
The following scripts are available, each of them provides a detailed usage description on using the `--help` flag.

- `transform_dataset.py`: Transform datasets saved in text format as `.gz` files to `.npy` format.
- `ds_info.py`: Summarize the dataset and clusters using any available ground truth labels.
- `sklearn_cluster.py`: Run KMeans, KMeans++, Meanshift, DBSCAN, and Hierarchical clustering with single and complete linkage on a given dataset and output a summary of the scores of each algorithm.
  In order to automate running several clustering algorithms, the script chooses certain parameters based on heuristics e.g., the epsilon parameter for DBSCAN (refer `cluster/algo.py`).
  However, manually setting the parameter value might lead to significant improvement in clustering performance.

Execute the following commands from the `utils` directory to run the above scripts on the sample dataset:

```sh
# Change dataset format.
python3 transform_dataset.py ../data/dataset.data.gz ./transformed.npy

# Summarize the dataset.
python3 ds_info.py ../data dataset

# Run common clustering algorithms on input dataset.
# Outputs a file named `dataset.json` with a summary of the performance of different clustering algorithms.
python3 sklearn_cluster.py ../data dataset
```
