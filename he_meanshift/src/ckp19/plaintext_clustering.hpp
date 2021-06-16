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

#pragma once

#include <random>
#include <vector>

class PlainMeanShift {
 public:
  using Point = std::vector<double>;
  long seed;

  PlainMeanShift();
  PlainMeanShift(long seed);

  double inverse(double inp, double m, long steps) const;
  double kernel(const Point &x, const Point &y, long degree) const;
  void minIdx(std::vector<double> &res, const std::vector<double> &x, long t,
              long inv_steps) const;

  void modeSeeking(std::vector<Point> &res, const std::vector<Point> &points,
                   long d, long steps, long kdegree, long inv_steps) const;
  void pointLabeling(std::vector<std::vector<double>> &res,
                     const std::vector<Point> &points,
                     const std::vector<Point> &dusts, long kdegree,
                     long inv_steps, long minidx_t) const;
  void cluster(std::vector<std::vector<double>> &res,
               const std::vector<Point> &points, long num_dusts,
               long iterations, long mode_kdegree, long label_kdegree,
               long mode_inv_steps, long label_inv_steps, long minidx_t);
};
