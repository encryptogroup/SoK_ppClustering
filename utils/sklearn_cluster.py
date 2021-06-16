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

import json
import cluster.algo as cluster
from cluster.dataset import Dataset
import argparse
import os


def run_clustering_algorithms(ds, seeds=None):
    if seeds is None:
        seeds = list(range(200, 210))

    perf_data = {}

    perf_data["kmeans"] = []
    perf_data["kmeans++"] = []
    km = cluster.KMeans(ds)
    kmpp = cluster.KMeansPP(ds)
    for n_cluster in ds.n_clusters:
        perf_data["kmeans"].append(km(n_cluster, seeds=seeds))
        perf_data["kmeans++"].append(kmpp(n_cluster, seeds=seeds))
    print("Finished KMeans.")

    ms = cluster.Meanshift(ds)
    perf_data["meanshift"] = ms()
    print("Finished MeanShift.")

    db = cluster.DBSCAN(ds)
    perf_data["dbscan"] = db()
    print("Finished DBSCAN.")

    perf_data["hs_single_linkage"] = []
    perf_data["hs_complete_linkage"] = []
    hs = cluster.Hierarchical(ds)
    for n_cluster in ds.n_clusters:
        perf_data["hs_single_linkage"].append(hs("single", n_cluster))
        perf_data["hs_complete_linkage"].append(hs("complete", n_cluster))
    print("Finished Hierarchical.")

    return perf_data


def cli_args():
    parser = argparse.ArgumentParser(
        description="Run sklearn clustering algorithms on a dataset"
    )
    parser.add_argument(
        "data_dir", help="Path to directory where dataset and labels are stored."
    )
    parser.add_argument(
        "dataset",
        help="Name of dataset. Expects a file by name 'dataset.gz' and 'labels*.gz' containing numpy matrices in text format.",
    )
    parser.add_argument(
        "--output-dir", default="", help="Path to output directory where data is saved."
    )
    parser.add_argument("--seed", default=2603, type=int, help="Seed for RNG.")

    args = parser.parse_args()
    return args


if __name__ == "__main__":
    args = cli_args()
    ds = Dataset.load_gz(args.data_dir, args.dataset)
    ds.print_stats()
    print()

    perf_data = run_clustering_algorithms(ds, list(range(args.seed, args.seed + 10)))
    perf_data["groundtruth"] = ds.groundtruth_external_metrics()

    with open(os.path.join(args.output_dir, f"{args.dataset}.json"), "w") as f:
        json.dump(perf_data, f)
