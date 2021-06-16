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
from pathlib import Path
import numpy as np
from .score import Score


class Dataset:
    @staticmethod
    def load_gz(data_dir, dname):
        ds = Dataset()

        data_dir = Path(data_dir)
        ds.data = np.loadtxt((data_dir / f"{dname}.data.gz").resolve(), ndmin=2)
        ds.labels = []

        for f in data_dir.glob(f"{dname}.labels*.gz"):
            # Subtract 1 to ensure label for noise is -1
            label = np.loadtxt(f.resolve(), dtype=np.intc) - 1
            ds.labels.append(label)

        # Number of clusters is the derived from the first label set
        ds.n_clusters = set([len(np.unique(label[label != -1])) for label in ds.labels])

        return ds

    @staticmethod
    def random(n_samples, n_features, centers, random_state):
        ds = Dataset()
        X, y = skdatasets.make_blobs(
            n_samples=n_samples,
            n_features=n_features,
            centers=centers,
            random_state=random_state,
        )

        ds.data = X
        ds.labels = [y]
        ds.n_clusters = centers

        return ds

    def print_stats(self):
        npoints, dim = self.data.shape
        print("Number of data points:", npoints)
        print("Dimension of data point:", dim)
        print("Number of labels (groundtruth):", len(self.labels))

        print("Clusters per label:", end=" ")
        for c in self.n_clusters:
            print(c, end=" ")
        print()

        print("Noise points per label:", end=" ")
        for lbl in self.labels:
            print((lbl == -1).sum(), end=" ")
        print()

    def groundtruth_external_metrics(self):
        return [
            Score.evaluate(self.data, label, label).external() for label in self.labels
        ]
