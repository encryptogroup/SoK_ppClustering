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

#include "dataset.hpp"

#include <cnpy.h>

#include <algorithm>
#include <numeric>
#include <random>
#include <stdexcept>
#include <string>

using Point = Dataset::Point;

Dataset::Dataset(long dim, long npoints, long seed)
    : data(npoints, Point(dim)) {
  std::mt19937 gen(seed);
  std::uniform_real_distribution<> dist(-1.0, 1.0);

  auto next_val = [&]() -> double { return dist(gen); };

  for (auto &i : data) {
    std::generate(i.begin(), i.end(), next_val);
  }
}

Dataset Dataset::loadNpy(const std::string &fpath) {
  auto npy_arr = cnpy::npy_load(fpath);

  if (npy_arr.shape.size() != 2) {
    throw std::logic_error("Expected 2D matrix for dataset.");
  }

  size_t rows = npy_arr.shape[0];
  size_t cols = npy_arr.shape[1];

  Dataset ds;
  ds.data.resize(rows, Point(cols));
  auto npy_data = npy_arr.data<double>();

  for (size_t i = 0; i < rows; ++i) {
    for (size_t j = 0; j < cols; ++j) {
      ds.data[i][j] = npy_data[i * cols + j];
    }
  }

  return ds;
}

const Point &Dataset::operator[](std::size_t idx) const { return data[idx]; }

const std::vector<Point> &Dataset::getPoints() const { return data; }

size_t Dataset::size() const { return data.size(); }

double Dataset::l2Norm(const Point &p) const {
  double res = 0;
  for (auto const &i : p) {
    res += i * i;
  }

  return std::sqrt(res);
}

void Dataset::rescaleDatasetInPlace(double r) {
  double max_norm = -1 * std::numeric_limits<double>::infinity();
  for (auto const &p : data) {
    max_norm = std::max(l2Norm(p), max_norm);
  }

  max_norm /= r;

  for (auto &p : data) {
    for (auto &i : p) {
      i /= max_norm;
    }
  }
}

Dataset Dataset::rescaleDataset(double r) const {
  Dataset scaled(*this);
  scaled.rescaleDatasetInPlace(r);

  return scaled;
}

void Dataset::encryptSIMD(EncryptedDataset &res, Scheme &scheme, long num_dusts,
                          long slots, long logp, long logq) const {
  size_t num_points = data.size();
  size_t dim = data[0].size();

  // Number of points in a ciphertext, excluding repition of data for each dust
  // i.e., number of unique points.
  long npoints = slots / (dim * num_dusts);

  // All points can be packed into one ciphertext
  if (num_points < npoints) npoints = num_points;

  // Current implementation assumes that there are exactly the same number of
  // points per ciphertext. We can either find closet divisor of `num_points`
  // lesser than `npoints` or modify the mean shift implementation to allow
  // different number of plaintex points per ciphertext. Both might lead to
  // performance losses. Currently, this case is not supported.
  if (num_points % npoints != 0) {
    std::string err =
        "Incompatible dataset (num_points=" + std::to_string(num_points) +
        ", npoints=" + std::to_string(npoints) + ")";
    throw std::runtime_error(err);
  }

  long num_cts = num_points / npoints;
  std::vector<Ciphertext> cpoints(num_cts);

  // Ideally, we need only `dim * npoints * num_dusts` slots per
  // ciphertext. However, HEAAN allows encryption only if number of slots
  // are a power of 2.
  long log_slots = std::ceil(std::log2((double)dim * npoints * num_dusts));
  long req_slots = 1 << log_slots;

  std::vector<double> plaintext(req_slots, 0);

  for (long i = 0; i < num_cts; ++i) {
    for (long j = 0; j < npoints; ++j) {
      for (long k = 0; k < dim; ++k) {
        plaintext[j * dim + k] = data[i * npoints + j][k];
      }
    }

    // repeat points for each dust
    for (long j = 1; j < num_dusts; ++j) {
      std::copy_n(plaintext.begin(), npoints * dim,
                  plaintext.begin() + j * npoints * dim);
    }

    scheme.encrypt(cpoints[i], plaintext.data(), req_slots, logp, logq);
  }

  res.dim = dim;
  res.npoints = npoints;
  res.log_slots = log_slots;
  res.points = std::move(cpoints);
}
