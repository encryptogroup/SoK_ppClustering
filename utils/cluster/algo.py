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

import sklearn.cluster as skcluster
from sklearn.neighbors import NearestNeighbors
from sklearn.preprocessing import MinMaxScaler
import numpy as np
import timeit
import kneed
from .score import Score


class KMeans:
    def __init__(self, dataset):
        self.ds = dataset

    def __call__(self, n_clusters, *, seeds=200, max_iter=20, data_range=(0, 100)):
        if isinstance(seeds, int):
            seeds = [seeds]

        tol = 1e-4
        if max_iter is not None:
            tol = 0

        data = MinMaxScaler(data_range).fit_transform(self.ds.data)

        score_list = []
        time = 0
        iterations = 0

        for seed in seeds:
            rng = np.random.default_rng(seed)

            start = timeit.default_timer()
            cluster_init = rng.uniform(
                data_range[0], data_range[1], (n_clusters, data.shape[1])
            )
            km = skcluster.KMeans(
                n_clusters,
                init=cluster_init,
                n_init=1,
                max_iter=max_iter,
                tol=tol,
                random_state=seed,
                algorithm="full",
            )
            pred = km.fit_predict(data)
            end = timeit.default_timer()

            time += end - start
            iterations += km.n_iter_
            score_list.append(Score.evaluate_on_dataset(self.ds, pred))

        mean_score_list = [Score.mean(i).all() for i in zip(*score_list)]

        return {
            "time": time / len(seeds),
            "iterations": iterations / len(seeds),
            "score": mean_score_list,
            "seeds": seeds,
            "n_clusters": n_clusters,
        }


# KMeans++
class KMeansPP:
    def __init__(self, dataset):
        self.ds = dataset

    def __call__(self, n_clusters, *, seeds=200, max_iter=20, data_range=(0, 100)):
        if isinstance(seeds, int):
            seeds = [seeds]

        tol = 1e-4
        if max_iter is not None:
            tol = 0

        data = MinMaxScaler(data_range).fit_transform(self.ds.data)

        score_list = []
        time = 0
        iterations = 0

        for seed in seeds:
            start = timeit.default_timer()
            km = skcluster.KMeans(
                n_clusters,
                init="k-means++",
                n_init=1,
                max_iter=max_iter,
                tol=tol,
                random_state=seed,
                algorithm="full",
            )
            pred = km.fit_predict(data)
            end = timeit.default_timer()

            time += end - start
            iterations += km.n_iter_
            score_list.append(Score.evaluate_on_dataset(self.ds, pred))

        mean_score_list = [Score.mean(i).all() for i in zip(*score_list)]

        return {
            "time": time / len(seeds),
            "iterations": iterations / len(seeds),
            "score": mean_score_list,
            "seeds": seeds,
            "n_clusters": n_clusters,
        }


class Meanshift:
    def __init__(self, dataset):
        self.ds = dataset

    def __call__(self, *, max_iter=20, bandwidth=None):
        start = timeit.default_timer()
        ms = skcluster.MeanShift(bandwidth=bandwidth, max_iter=max_iter)
        pred = ms.fit_predict(self.ds.data)
        end = timeit.default_timer()

        time = end - start
        iterations = ms.n_iter_
        score = Score.evaluate_on_dataset(self.ds, pred)

        return {
            "time": time,
            "iterations": iterations,
            "score": [s.all() for s in score],
            "bandwidth": bandwidth,
        }


class DBSCAN:
    def __init__(self, dataset):
        self.ds = dataset

    def estimate_eps(self):
        nn = NearestNeighbors(n_neighbors=2, metric="sqeuclidean")
        neighbours = nn.fit(self.ds.data)
        distances, indices = neighbours.kneighbors(self.ds.data)
        distances = np.sort(distances, axis=0)
        distances = distances[:, 1]
        kneedle = kneed.KneeLocator(
            range(len(self.ds.data)),
            distances,
            S=1.0,
            curve="convex",
            direction="increasing",
            interp_method="polynomial",
        )
        return kneedle.elbow_y

    def __call__(self, eps=None, *, min_samples=5):
        if eps is None:
            eps = self.estimate_eps()

        start = timeit.default_timer()
        db = skcluster.DBSCAN(
            eps=eps, min_samples=min_samples, metric="sqeuclidean", algorithm="brute"
        )
        pred = db.fit_predict(self.ds.data)
        end = timeit.default_timer()

        time = end - start
        score = Score.evaluate_on_dataset(self.ds, pred)

        return {
            "time": time,
            "eps": eps,
            "min_samples": min_samples,
            "score": [s.all() for s in score],
        }


class Hierarchical:
    def __init__(self, dataset):
        self.ds = dataset

    def __call__(self, linkage, n_clusters):
        start = timeit.default_timer()
        hc = skcluster.AgglomerativeClustering(
            n_clusters=n_clusters,
            linkage=linkage,
            affinity="sqeuclidean",
        )
        pred = hc.fit_predict(self.ds.data)
        end = timeit.default_timer()

        time = end - start
        score = Score.evaluate_on_dataset(self.ds, pred)

        return {
            "time": time,
            "score": [s.all() for s in score],
            "linkage": linkage,
            "n_clusters": n_clusters,
        }
