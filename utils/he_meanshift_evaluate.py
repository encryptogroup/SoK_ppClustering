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
from cluster.score import Score
from cluster.dataset import Dataset
import argparse
import os
import numpy as np
from pathlib import Path


def encode_label(v):
    res = 0
    for i in range(len(v)):
        if v[i] == 1:
            res += (2**i)
    return res


def evaluate(ds, pred_dir, threshold):
    pred_dir = Path(pred_dir)
    score_list = []

    for rep in pred_dir.glob("rep_*.npy"):
        pred = np.load(rep.resolve())
        pred[pred >= threshold] = 1
        pred[pred < threshold] = 0
        labels = np.array([encode_label(i) for i in pred])
        score_list.append(Score.evaluate_on_dataset(ds, labels))

    mean_score_list = [Score.mean(i).all() for i in zip(*score_list)]
    return mean_score_list


def cli_args():
    parser = argparse.ArgumentParser(
        description="Evaluate output of CKP19 meanshift."
    )
    parser.add_argument(
        "data_dir",
        help="Path to directory where dataset and labels are stored.",
    )
    parser.add_argument(
        "dataset",
        help="Name of dataset. Expects a file by name 'dataset.gz' and 'labels*.gz' containing numpy matrices in text format.",
    )
    parser.add_argument(
        "pred_dir",
        help="Path to directory containing output of clustering.",
    )
    parser.add_argument(
        "--output_dir", default="", help="Path to output directory."
    )
    parser.add_argument(
        "--threshold", default=0.1, type=float, help="Threshold for converting labels to 0 and 1."
    )

    args = parser.parse_args()
    return args


if __name__ == "__main__":
    args = cli_args()
    ds = Dataset.load_gz(args.data_dir, args.dataset)
    ds.print_stats()
    print()

    output = {}
    output["scores"] = evaluate(ds, args.pred_dir, args.threshold)
    pred_dir = Path(args.pred_dir)

    for i, score in enumerate(output["scores"]):
        print("--- Mean Scores for groundtruth label", i, "---")
        for k in score:
            print(k, ":", score[k])
        print()

    with open(os.path.join(args.output_dir, f"{pred_dir.name}.json"), "w") as f:
        json.dump(output, f)
