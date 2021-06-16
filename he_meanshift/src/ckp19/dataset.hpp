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

#include <HEAAN.h>

#include <string>
#include <vector>

// Packed ciphertext with associated data
struct EncryptedDataset {
  // Dimension of each point
  long dim;

  // Number of unique points per ciphertext, i.e., excludes repetition for
  // dusts
  long npoints;

  // Log of number of slots used per ciphertext. Power of 2 to greater than
  // or equal to `dim * npoints * num_dusts`
  long log_slots;
  std::vector<Ciphertext> points;
};

class Dataset {
 public:
  using Point = std::vector<double>;

  Dataset() = default;

  // Generate a dataset with `npoints` random points eaching having `dim`
  // attributes.
  Dataset(long dim, long npoints, long seed = 0);

  // Read matrix from a `npy` file.
  static Dataset loadNpy(const std::string &fpath);

  const Point &operator[](std::size_t idx) const;
  const std::vector<Point> &getPoints() const;

  size_t size() const;

  double l2Norm(const Point &p) const;
  Dataset rescaleDataset(double r) const;
  void rescaleDatasetInPlace(double r);
  void encryptSIMD(EncryptedDataset &res, Scheme &scheme, long num_dusts,
                   long slots, long logp, long logq) const;

 private:
  std::vector<Point> data;
  long dim;
};
