// MIT License
//
// Copyright (c) 2021 Aditya Shridhar Hegde
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

#include "plaintext_clustering.hpp"

#include <random>

using Point = PlainMeanShift::Point;

PlainMeanShift::PlainMeanShift() { seed = std::rand(); }

PlainMeanShift::PlainMeanShift(long seed) : seed(seed) {}

double PlainMeanShift::inverse(double inp, double m, long steps) const {
  inp *= 2.0 / m;
  double a = 2 - inp;
  double b = 1 - inp;

  for (int i = 0; i < steps; ++i) {
    b *= b;
    a *= (1 + b);
  }

  return (a * 2.0) / m;
}

double PlainMeanShift::kernel(const Point &x, const Point &y,
                              long degree) const {
  double res = 0;
  for (size_t i = 0; i < x.size(); ++i) {
    res += (x[i] - y[i]) * (x[i] - y[i]);
  }

  res = 1 - res;

  for (long i = 0; i < degree; ++i) res *= res;

  return res;
}

void PlainMeanShift::minIdx(std::vector<double> &res,
                            const std::vector<double> &x, long t,
                            long inv_steps) const {
  res.resize(x.size());

  double sum = 0;

  for (size_t i = 0; i < x.size(); ++i) {
    res[i] = 1 - x[i];
    for (long j = 0; j < t; ++j) {
      res[i] *= res[i];
    }

    sum += res[i];
  }

  for (auto &i : res) {
    i = i * inverse(sum, x.size(), inv_steps);
  }
}

void PlainMeanShift::modeSeeking(std::vector<Point> &res,
                                 const std::vector<Point> &points, long d,
                                 long steps, long kdegree,
                                 long inv_steps) const {
  res.resize(d);

  std::mt19937 gen(seed);
  std::uniform_int_distribution<> distrib(0, points.size() - 1);
  for (long i = 0; i < d; ++i) {
    res[i] = points[distrib(gen)];
  }

  for (long iter = 0; iter < steps; ++iter) {
    for (auto &dust : res) {
      double sum = 0;
      auto A = std::vector<double>(dust.size(), 0);

      for (auto &point : points) {
        double a = kernel(dust, point, kdegree);
        sum += a;
        for (size_t i = 0; i < A.size(); ++i) {
          A[i] += a * point[i];
        }
      }

      for (size_t i = 0; i < A.size(); ++i) {
        dust[i] = A[i] * inverse(sum, points.size(), inv_steps);
      }
    }
  }
}

void PlainMeanShift::pointLabeling(std::vector<std::vector<double>> &res,
                                   const std::vector<Point> &points,
                                   const std::vector<Point> &dusts,
                                   long kdegree, long inv_steps,
                                   long minidx_t) const {
  std::vector<double> nhbd(dusts.size(), 0);

  for (size_t i = 0; i < dusts.size(); ++i) {
    for (auto &dust : dusts) {
      nhbd[i] += kernel(dusts[i], dust, kdegree);
    }
  }

  res.resize(points.size());
  for (size_t i = 0; i < points.size(); ++i) {
    std::vector<double> temp(dusts.size());

    for (size_t j = 0; j < dusts.size(); ++j) {
      for (size_t k = 0; k < points[i].size(); ++k) {
        temp[j] += (points[i][k] - dusts[j][k]) * (points[i][k] - dusts[j][k]);
      }
    }

    minIdx(res[i], temp, minidx_t, inv_steps);

    for (size_t j = 0; j < dusts.size(); ++j) {
      res[i][j] *= nhbd[j];
    }
  }
}

void PlainMeanShift::cluster(std::vector<std::vector<double>> &res,
                             const std::vector<Point> &points, long num_dusts,
                             long iterations, long mode_kdegree,
                             long label_kdegree, long mode_inv_steps,
                             long label_inv_steps, long minidx_t) {
  std::vector<Point> dusts;
  modeSeeking(dusts, points, num_dusts, iterations, mode_kdegree,
              mode_inv_steps);
  pointLabeling(res, points, dusts, label_kdegree, label_inv_steps, minidx_t);
}
