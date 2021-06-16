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

import sklearn.metrics as skmetrics
import numpy as np


class Score:
    DEFAULT = {
        "adjusted_rand_index": -1,
        "adjusted_mutual_info": 0,
        "silhouette": -1,
        "calinski_harabasz": 0,
        "generalized_dunn_index_gd33": 0,
    }

    EXTERNAL = ["silhouette", "calinski_harabasz", "generalized_dunn_index_gd33"]
    INTERNAL = ["adjusted_rand_index", "adjusted_mutual_info"]

    def __init__(self):
        self.scores_ = {**Score.DEFAULT}
        self.scale_factor_ = 1.0

    def external(self):
        if self.scores_ is None:
            return {k: None for k in Score.EXTERNAL}

        res = {k: self.scores_[k] for k in Score.EXTERNAL}
        res["scale_factor"] = self.scale_factor_
        return res

    def internal(self):
        if self.scores_ is None:
            return {k: None for k in Score.INTERNAL}

        res = {k: self.scores_[k] for k in Score.INTERNAL}
        res["scale_factor"] = self.scale_factor_
        return res

    def all(self):
        if self.scores_ is None:
            return {k: None for k in Score.DEFAULT}

        res = {**self.scores_}
        res["scale_factor"] = self.scale_factor_
        return res

    @staticmethod
    def mean(iterable):
        score_list = list(filter(lambda x: x.scores_ is not None, iterable))
        mean = score_list[0]

        for score in score_list[1:]:
            for k in mean.scores_:
                mean.scores_[k] += score.scores_[k]

        for k in mean.scores_:
            mean.scores_[k] /= len(score_list)

        return mean

    @staticmethod
    def evaluate(data, groundtruth, pred):
        assert groundtruth.shape == pred.shape
        assert data.shape[0] == groundtruth.shape[0]

        num_total = len(pred[groundtruth != -1])
        num_noise_pred = ((groundtruth != -1) & (pred == -1)).sum()
        scale_factor = (num_total - num_noise_pred) / (num_total)

        selector = (groundtruth != -1) & (pred != -1)
        data = data[selector]
        pred = pred[selector]
        groundtruth = groundtruth[selector]

        scores = {}
        if len(np.unique(pred)) == 1:
            # if all points are mapped to a single label then assign score to None since
            # the metric don't make much sense in this case.
            scores = None
        else:
            scores["adjusted_rand_index"] = Score.adjusted_rand_index(groundtruth, pred)
            scores["adjusted_mutual_info"] = Score.adjusted_mutual_info(groundtruth, pred)
            scores["silhouette"] = Score.silhouette_score(data, pred)
            scores["calinski_harabasz"] = Score.calinski_harabasz_score(data, pred)
            scores["generalized_dunn_index_gd33"] = Score.generalized_dunn_index(data, pred)

            for k in scores:
                scores[k] *= scale_factor

        res = Score()
        res.scores_ = scores
        res.scale_factor_ = scale_factor
        return res

    @staticmethod
    def evaluate_on_dataset(ds, pred):
        return [Score.evaluate(ds.data, label, pred) for label in ds.labels]

    @staticmethod
    def adjusted_rand_index(groundtruth, pred):
        return skmetrics.adjusted_rand_score(groundtruth, pred)

    @staticmethod
    def adjusted_mutual_info(groundtruth, pred):
        return skmetrics.adjusted_mutual_info_score(groundtruth, pred)

    @staticmethod
    def silhouette_score(data, pred, metric="sqeuclidean"):
        if len(set(pred)) == 1:
            return None

        return skmetrics.silhouette_score(data, pred, metric=metric)

    @staticmethod
    def calinski_harabasz_score(data, pred):
        if len(set(pred)) == 1:
            return None

        return skmetrics.calinski_harabasz_score(data, pred)

    @staticmethod
    def generalized_dunn_index(
        data, pred, seperation_metric=None, cohesion_metric=None
    ):
        if len(set(pred)) == 1:
            return None

        def gd3_seperation_metric(c1, c2):
            numerator = skmetrics.pairwise_distances(c1, c2, metric="sqeuclidean").sum()
            denominator = c1.shape[0] * c2.shape[0]
            return numerator / denominator

        def gd3_cohesion_metric(cpoints):
            mean = cpoints.mean(0)
            numerator = (
                2
                * skmetrics.pairwise_distances(
                    cpoints, [mean], metric="sqeuclidean"
                ).sum()
            )
            return numerator / 2

        if seperation_metric is None:
            seperation_metric = gd3_seperation_metric
        if cohesion_metric is None:
            cohesion_metric = gd3_cohesion_metric

        labels = list(set(pred))

        numerator = float("inf")
        denominator = float("-inf")
        for i in range(len(labels)):
            denominator = max(denominator, cohesion_metric(data[pred == labels[i]]))
            for j in range(i + 1, len(labels)):
                if i == j:
                    continue
                numerator = min(
                    numerator,
                    seperation_metric(data[pred == labels[i]], data[pred == labels[j]]),
                )

        return numerator / denominator
