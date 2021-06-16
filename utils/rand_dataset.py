# MIT License
#
# Copyright (c) 2021 Aditya Shridhar Hegde
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

import sklearn.datasets as skdatasets
import numpy as np
import argparse
import os


def generate_and_store_dataset(n_clusters, n_points, n_dim, seed, data_dir):
    X, y = skdatasets.make_blobs(
        n_samples=n_points, n_features=n_dim, centers=n_clusters, random_state=seed
    )
    np.save(os.path.join(data_dir, f"n{n_points}_d{n_dim}_k{n_clusters}.npy"), X)


def cli_args():
    parser = argparse.ArgumentParser(description="Generate random datasets")
    parser.add_argument("clusters", help="Number of clusters.", type=int)
    parser.add_argument("points", help="Number of data points.", type=int)
    parser.add_argument("dim", help="Dimension of each point.", type=int)
    parser.add_argument("--seed", default=2603, help="Seed for RNG.")
    parser.add_argument(
        "--data-dir", default="", help="Path to directory where dataset will be saved."
    )

    args = parser.parse_args()
    return args


if __name__ == "__main__":
    args = cli_args()
    generate_and_store_dataset(
        args.clusters, args.points, args.dim, args.seed, args.data_dir
    )
