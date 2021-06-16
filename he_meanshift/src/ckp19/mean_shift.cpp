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

#include "mean_shift.hpp"

#include <NTL/ZZ.h>

#include <algorithm>
#include <numeric>
#include <random>
#include <stdexcept>

MeanShift::MeanShift(Scheme &scheme, long logn, long logq, long logq_boot,
                     long logp, long logt)
    : scheme(scheme),
      logn(logn),
      logq(logq),
      logq_boot(logq_boot),
      logp(logp),
      logt(logt) {
  slots = 1 << logn;
  refreshSeed();
}

MeanShift::MeanShift(Scheme &scheme, long logn, long logq, long logq_boot,
                     long logp, long logt, long seed)
    : scheme(scheme),
      logn(logn),
      logq(logq),
      logq_boot(logq_boot),
      logp(logp),
      logt(logt),
      seed(seed) {
  slots = 1 << logn;
}

void MeanShift::refreshSeed(long s) { seed = s; }

void MeanShift::refreshSeed() { seed = std::rand(); }

void MeanShift::leftRotateAndEqual(Ciphertext &x, long r) {
  while (r & (-r)) {
    scheme.leftRotateFastAndEqual(x, r & (-r));
    r -= (r & (-r));
  }
}

void MeanShift::rightRotateAndEqual(Ciphertext &x, long r) {
  while (r & (-r)) {
    scheme.rightRotateFastAndEqual(x, r & (-r));
    r -= (r & (-r));
  }
}

// x.logq := x.logq
void MeanShift::repeatSlotRight(Ciphertext &x, long rep, long step) {
  Ciphertext orig;
  Ciphertext temp;
  orig.copy(x);

  long rot = 1;
  for (long i = NumBits(rep) - 2; i >= 0; --i) {
    temp.copy(x);
    rightRotateAndEqual(temp, step * rot);
    scheme.addAndEqual(x, temp);
    rot *= 2;

    if (bit(rep, i)) {
      temp.copy(orig);
      rightRotateAndEqual(temp, rot * step);
      scheme.addAndEqual(x, temp);
      rot += 1;
    }
  }
}

// x.logq := x.logq
void MeanShift::sumSlotsLeft(Ciphertext &x, long num, long step) {
  Ciphertext orig;
  Ciphertext temp;
  orig.copy(x);

  long rot = 1;
  for (long i = NumBits(num) - 2; i >= 0; --i) {
    temp.copy(x);
    leftRotateAndEqual(temp, rot * step);
    scheme.addAndEqual(x, temp);
    rot *= 2;

    if (bit(num, i)) {
      temp.copy(orig);
      leftRotateAndEqual(temp, rot * step);
      scheme.addAndEqual(x, temp);
      rot += 1;
    }
  }
}

// x.logq := x.logq - logp
void MeanShift::clearSlotsRight(Ciphertext &x, long ret_size, long skip,
                                long max_slot) {
  std::vector<complex<double>> selector(slots, complex<double>{0, 0});

  for (long i = 0; i + ret_size <= max_slot; i += skip) {
    std::fill(selector.begin() + i, selector.begin() + i + ret_size,
              complex<double>{1, 0});
  }

  scheme.multByConstVecAndEqual(x, selector.data(), logp);
  scheme.reScaleByAndEqual(x, logp);
}

// res.logq = x.logq - (steps + 2) * logp
void MeanShift::inverse(Ciphertext &res, Ciphertext &x, double m, long steps) {
  Ciphertext scaled;
  scheme.multByConst(scaled, x, 2.0 / m, logp);
  scheme.reScaleByAndEqual(scaled, logp);

  SchemeAlgo schemeAlgo(scheme);
  schemeAlgo.inverse(res, scaled, scaled.logp, steps);

  scheme.multByConstAndEqual(res, 2.0 / m, res.logp);
  scheme.reScaleByAndEqual(res, logp);
}

// res[i].logq = x[i].logq - (inv_steps + t + 3) * logp
void MeanShift::minIdx(std::vector<Ciphertext> &res, std::vector<Ciphertext> &x,
                       long t, long inv_steps) {
  Ciphertext sum;
  res.resize(x.size());

  for (size_t i = 0; i < x.size(); ++i) {
    // res[i] := 1 - x[i]
    scheme.negate(res[i], x[i]);
    scheme.addConstAndEqual(res[i], 1.0, logp);

    // res[i] := res[i]^(2^t)
    for (long j = 0; j < t; ++j) {
      scheme.squareAndEqual(res[i]);
      scheme.reScaleByAndEqual(res[i], logp);
    }

    // sum := sum + res[i]
    if (i == 0)
      sum.copy(res[0]);
    else {
      scheme.addAndEqual(sum, res[i]);
    }
  }

  // inv := 1/sum
  Ciphertext inv;
  inverse(inv, sum, x.size(), inv_steps);

  // res[i] := res[i] / sum
  for (size_t i = 0; i < res.size(); ++i) {
    scheme.modDownToAndEqual(res[i], inv.logq);
    scheme.multAndEqual(res[i], inv);
    scheme.reScaleByAndEqual(res[i], logp);
  }
}

// res.logq := x.logq - (inv_steps + t + 4) logp
void MeanShift::minIdxSIMD(Ciphertext &res, Ciphertext &x, long dim,
                           long npoints, long t, long inv_steps,
                           long batchsize) {
  scheme.negate(res, x);
  scheme.addConstAndEqual(res, 1.0, logp);

  for (long j = 0; j < t; ++j) {
    scheme.squareAndEqual(res);
    scheme.reScaleByAndEqual(res, logp);
  }

  Ciphertext sum;
  sum.copy(res);
  sumSlotsLeft(sum, batchsize, npoints * dim);
  clearSlotsRight(sum, npoints * dim, slots, slots);
  repeatSlotRight(sum, batchsize, npoints * dim);

  Ciphertext inv;
  inverse(inv, sum, batchsize, inv_steps);

  scheme.modDownToAndEqual(res, inv.logq);
  scheme.multAndEqual(res, inv);
  scheme.reScaleByAndEqual(res, logp);
}

// res.logq = x.logq - logp
void MeanShift::l2NormSquared(Ciphertext &res, Ciphertext &x) {
  Ciphertext rot;
  scheme.square(res, x);
  scheme.reScaleByAndEqual(res, logp);
  rot.copy(res);
  for (long i = 1; i < slots; ++i) {
    scheme.leftRotateFastAndEqual(rot, 1);
    scheme.addAndEqual(res, rot);
  }
}

// res.logq := x.logq - 2 logp
void MeanShift::l2NormSquaredSIMD(Ciphertext &res, Ciphertext &x, long dim,
                                  long num_points) {
  scheme.square(res, x);
  scheme.reScaleByAndEqual(res, logp);
  sumSlotsLeft(res, dim, 1);
  clearSlotsRight(res, 1, dim, num_points * dim);
  repeatSlotRight(res, dim, 1);
}

// res.logq := x.logq - (degree + 1) * logp
void MeanShift::kernel(Ciphertext &res, Ciphertext &x, Ciphertext &y,
                       long degree) {
  // diff = (x - y)^2
  Ciphertext diff;
  scheme.sub(diff, x, y);

  l2NormSquared(res, diff);

  // res = (1 - res)^(2^degree)
  scheme.negateAndEqual(res);
  scheme.addConstAndEqual(res, 1.0, logp);
  for (long i = 0; i < degree; ++i) {
    scheme.squareAndEqual(res);
    scheme.reScaleByAndEqual(res, logp);
  }
}

// res.logp := x.logp - (degree + 2) logp
void MeanShift::kernelSIMD(Ciphertext &res, Ciphertext &x, Ciphertext &y,
                           long dim, long num_points, long degree) {
  // diff = (x - y)^2
  Ciphertext diff;
  scheme.sub(diff, x, y);

  l2NormSquaredSIMD(res, diff, dim, num_points);

  // res = (1 - res)^(2^degree)
  scheme.negateAndEqual(res);
  scheme.addConstAndEqual(res, 1.0, logp);
  for (long i = 0; i < degree; ++i) {
    scheme.squareAndEqual(res);
    scheme.reScaleByAndEqual(res, logp);
  }
}

void MeanShift::modeSeeking(std::vector<Ciphertext> &res,
                            std::vector<Ciphertext> &points, long d, long steps,
                            long kdegree, long inv_steps) {
  res.resize(d);

  // Choose dusts
  std::mt19937 gen(seed);
  std::uniform_int_distribution<> distrib(0, points.size() - 1);
  for (long i = 0; i < d; ++i) {
    res[i].copy(points[distrib(gen)]);
  }

  // Gradient ascent
  for (long iter = 0; iter < steps; ++iter) {
    // dust.logq := dust.logq - (kdegree + inv_steps + 4) * logp
    for (auto &dust : res) {
      // Bootstrap if modulus not large enough for one more iteration
      if (dust.logq < (kdegree + inv_steps + 4) * logp + 1) {
        scheme.bootstrapAndEqual(dust, logq_boot, logQ, logt);
      }

      Ciphertext sum, A, a, tpoint, temp;

      for (size_t k = 0; k < points.size(); ++k) {
        if (points[k].logq != dust.logq)
          scheme.modDownTo(tpoint, points[k], dust.logq);
        else
          tpoint.copy(points[k]);

        kernel(a, tpoint, dust, kdegree);

        scheme.modDownTo(temp, tpoint, a.logq);
        scheme.multAndEqual(temp, a);
        scheme.reScaleByAndEqual(temp, logp);

        if (k == 0) {
          sum.copy(a);
          A.copy(temp);
        } else {
          scheme.addAndEqual(sum, a);
          scheme.addAndEqual(A, temp);
        }
      }

      Ciphertext inv;
      inverse(inv, sum, points.size(), inv_steps);

      scheme.modDownToAndEqual(A, inv.logq);
      scheme.mult(dust, inv, A);
      scheme.reScaleByAndEqual(dust, logp);
    }
  }
}

void MeanShift::refreshDusts(Ciphertext &dusts, long dim, long npoints,
                             long d) {
  Ciphertext temp1, temp2;

  // The cost of bootstrapping is proportional to the number of underlying
  // plaintext slots. Bootstrapping constitutes one of the most expensive
  // operations in the algorithm.
  //
  // Currently dusts are of the form (d1, d1, ..., d1, d2, d2, ..., dd).
  // We want convert it to (d1, d2, d3, ..., dd) to reduce the number of
  // plaintext slots.
  std::vector<complex<double>> selector(slots, complex<double>{0, 0});
  for (long i = 0; i < d; ++i) {
    std::fill(selector.begin() + i * npoints * dim,
              selector.begin() + i * npoints * dim + dim,
              complex<double>{1, 0});
    scheme.multByConstVec(temp1, dusts, selector.data(), logp);
    scheme.reScaleByAndEqual(temp1, logp);
    std::fill(selector.begin() + i * npoints * dim,
              selector.begin() + i * npoints * dim + dim,
              complex<double>{0, 0});

    if (i == 0) {
      temp2.copy(temp1);
    } else {
      leftRotateAndEqual(temp1, i * dim * (npoints - 1));
      scheme.addAndEqual(temp2, temp1);
    }
  }

  long step = 1 << ((long)std::ceil(std::log2(dim * d)));
  repeatSlotRight(temp2, slots / step, step);
  auto temp = temp2.n;
  temp2.n = step;

  scheme.bootstrapAndEqual(temp2, logq_boot, logQ, logt);
  temp2.n = temp;

  for (long i = 0; i < d; ++i) {
    std::fill(selector.begin() + i * dim, selector.begin() + (i + 1) * dim,
              complex<double>{1, 0});
    scheme.multByConstVec(temp1, temp2, selector.data(), logp);
    scheme.reScaleByAndEqual(temp1, logp);
    std::fill(selector.begin() + i * dim, selector.begin() + (i + 1) * dim,
              complex<double>{0, 0});

    if (i == 0) {
      dusts.copy(temp1);
    } else {
      rightRotateAndEqual(temp1, i * dim * (npoints - 1));
      scheme.addAndEqual(dusts, temp1);
    }
  }
}

void MeanShift::modeSeekingSIMD(Ciphertext &dusts,
                                std::vector<Ciphertext> &points, long dim,
                                long npoints, long d, long steps, long kdegree,
                                long inv_steps) {
  Ciphertext temp;

  // Choose dusts
  std::mt19937 gen(seed);
  // Can also choose a ciphertext at random and then a point in the ciphertext
  // at random. However, this approach allows easy verification when computed
  // in plaintext (check plaintext_clustering).
  std::uniform_int_distribution<> distrib(0, npoints * points.size() - 1);
  std::vector<complex<double>> selector(slots, complex<double>{0, 0});

  // Assume d < slots / (dim * npoints)
  for (long i = 0; i < d; ++i) {
    long r = distrib(gen);
    long pi = r % npoints;
    long ci = r / npoints;

    std::fill(selector.begin() + pi * dim, selector.begin() + (pi + 1) * dim,
              complex<double>{1, 0});
    scheme.multByConstVec(temp, points[ci], selector.data(), logp);
    scheme.reScaleByAndEqual(temp, logp);
    std::fill(selector.begin() + pi * dim, selector.begin() + (pi + 1) * dim,
              complex<double>{0, 0});

    long beg = i * npoints;
    if (beg < pi) {
      leftRotateAndEqual(temp, (pi - beg) * dim);
    } else if (beg > pi) {
      rightRotateAndEqual(temp, (beg - pi) * dim);
    }

    if (i == 0)
      dusts.copy(temp);
    else
      scheme.addAndEqual(dusts, temp);
  }

  // dusts will be repeated in loop below. The current dust format is expected
  // at the start of each iteration.

  // Gradient ascent
  // Now, dusts.logq = points.logq - logp

  Ciphertext grad, sum, A, temp2, tpoint;

  for (long iter = 0; iter < steps; ++iter) {
    // Bootstrap if modulus not large enough for one more iteration
    if (dusts.logq - (kdegree + inv_steps + 6) * logp < logq_boot) {
      refreshDusts(dusts, dim, npoints, d);
    }

    if (dusts.logq - (kdegree + inv_steps + 6) * logp < logq_boot) {
      throw std::runtime_error(
          "Modulus not large enough for mode seeking after bootstrap.");
    }

    // In the beginning of each iteration, only the first `dim` slots of dusts
    // have the right value and everything else has 0.
    repeatSlotRight(dusts, npoints, dim);

    for (size_t i = 0; i < points.size(); ++i) {
      scheme.modDownTo(tpoint, points[i], dusts.logq);
      kernelSIMD(grad, dusts, tpoint, dim, npoints * d, kdegree);

      scheme.modDownTo(temp2, tpoint, grad.logq);
      scheme.mult(temp, grad, temp2);
      scheme.reScaleByAndEqual(temp, logp);
      sumSlotsLeft(temp, npoints, dim);
      if (i == 0) {
        A.copy(temp);
      } else {
        scheme.addAndEqual(A, temp);
      }

      sumSlotsLeft(grad, npoints, dim);
      if (i == 0) {
        sum.copy(grad);
      } else {
        scheme.addAndEqual(sum, grad);
      }
    }

    Ciphertext inv;
    inverse(inv, sum, npoints * points.size(), inv_steps);

    scheme.modDownToAndEqual(A, inv.logq);
    scheme.mult(dusts, inv, A);
    scheme.reScaleByAndEqual(dusts, logp);
    clearSlotsRight(dusts, dim, dim * npoints, d * npoints * dim);
  }

  repeatSlotRight(dusts, npoints, dim);
}

// z := max(kdegree + 2, inv_steps + minidx_t + 5)
// res.logq = dusts.logq - z * logp
void MeanShift::pointLabeling(std::vector<std::vector<Ciphertext>> &res,
                              std::vector<Ciphertext> &dusts,
                              std::vector<Ciphertext> &points, long kdegree,
                              long inv_steps, long minidx_t) {
  std::vector<Ciphertext> nbhd(dusts.size());
  Ciphertext temp;

  for (size_t i = 0; i < dusts.size(); ++i) {
    for (size_t j = 0; j < dusts.size(); ++j) {
      kernel(temp, dusts[i], dusts[j], kdegree);

      if (j == 0) {
        nbhd[i].copy(temp);
      } else {
        scheme.addAndEqual(nbhd[i], temp);
      }
    }
  }

  std::vector<Ciphertext> Cbar;
  std::vector<Ciphertext> norms(dusts.size());

  res.resize(points.size());

  for (size_t i = 0; i < points.size(); ++i) {
    for (size_t j = 0; j < dusts.size(); ++j) {
      scheme.sub(temp, points[i], dusts[j]);
      l2NormSquared(norms[j], temp);
    }

    minIdx(Cbar, norms, minidx_t, inv_steps);
    res[i].resize(dusts.size());

    for (size_t j = 0; j < dusts.size(); ++j) {
      // make compatible for multiplication
      if (Cbar[j].logq < nbhd[j].logq) {
        // this will be run once for all points
        scheme.modDownToAndEqual(nbhd[j], Cbar[j].logq);
      } else if (Cbar[j].logq > nbhd[j].logq) {
        scheme.modDownToAndEqual(Cbar[j], nbhd[j].logq);
      }

      scheme.mult(res[i][j], Cbar[j], nbhd[j]);
      scheme.reScaleByAndEqual(res[i][j], logp);
    }
  }
}

// res.logq := dusts.logq - max(kdegree + 5, inv_steps + minidx_t + 7) logp
void MeanShift::pointLabelingSIMD(std::vector<Ciphertext> &res,
                                  Ciphertext &dusts,
                                  std::vector<Ciphertext> &points, long dim,
                                  long npoints, long d, long kdegree,
                                  long inv_steps, long minidx_t) {
  Ciphertext temp1, temp2, nbhd;

  // Currently dusts are of the form (d1, d1, ..., d1, d2, d2, ..., dd).
  // We want to now convert it to
  // (d1, d2, d3, ..., dd, 0, .., 0, d1, d2, d3, ... (repeated d times), dd, 0,
  // .., 0) so that we can easily compute pairwise kernel on the dusts.
  std::vector<complex<double>> selector(slots, complex<double>{0, 0});
  for (long i = 0; i < d; ++i) {
    std::fill(selector.begin() + i * npoints * dim,
              selector.begin() + i * npoints * dim + dim,
              complex<double>{1, 0});
    scheme.multByConstVec(temp1, dusts, selector.data(), logp);
    scheme.reScaleByAndEqual(temp1, logp);
    std::fill(selector.begin() + i * npoints * dim,
              selector.begin() + i * npoints * dim + dim,
              complex<double>{0, 0});

    if (i == 0) {
      temp2.copy(temp1);
    } else {
      leftRotateAndEqual(temp1, i * dim * (npoints - 1));
      scheme.addAndEqual(temp2, temp1);
    }
  }

  // `temp2` will now contain the dusts in the required format.
  repeatSlotRight(temp2, d, npoints * dim);

  // Compute `nbhd` which is the sum of the kernel value of a dust to every
  // other dust. (Refer paper for details)
  scheme.modDownBy(temp1, dusts, logp);
  kernelSIMD(nbhd, temp2, temp1, dim, npoints * d, kdegree);
  sumSlotsLeft(nbhd, d, dim);
  clearSlotsRight(nbhd, dim, npoints * dim, d * npoints * dim);
  repeatSlotRight(nbhd, npoints, dim);

  // Compute labels for each point ciphertext and store result in `res`.
  res.resize(points.size());
  for (size_t i = 0; i < points.size(); ++i) {
    scheme.sub(temp1, dusts, points[i]);
    l2NormSquaredSIMD(temp2, temp1, dim, npoints * d);
    minIdxSIMD(temp1, temp2, dim, npoints, minidx_t, inv_steps, d);

    if (nbhd.logq < temp1.logq)
      scheme.modDownToAndEqual(temp1, nbhd.logq);
    else if (nbhd.logq > temp1.logq)
      scheme.modDownToAndEqual(nbhd, temp1.logq);

    scheme.mult(res[i], temp1, nbhd);
    scheme.reScaleByAndEqual(res[i], logp);
  }
}

void MeanShift::clusterSIMD(std::vector<Ciphertext> &res,
                            std::vector<Ciphertext> &points, long dim,
                            long npoints, long num_dusts, long iterations,
                            long mode_kdegree, long label_kdegree,
                            long mode_inv_steps, long label_inv_steps,
                            long minidx_t) {
  Ciphertext dusts;
  modeSeekingSIMD(dusts, points, dim, npoints, num_dusts, iterations,
                  mode_kdegree, mode_inv_steps);

  // `dusts` has a lower modulus
  // Check if modulus suffices for point labeling
  if (dusts.logq < logq_boot)
    throw std::runtime_error(
        "Modulus lower than logl after mode seeking. Cannot bootstrap.");

  if (dusts.logq <
      std::max(label_kdegree + 5, label_inv_steps + minidx_t + 7) * logp) {
    refreshDusts(dusts, dim, npoints, num_dusts);
    repeatSlotRight(dusts, npoints, dim);
  }

  if (dusts.logq <
      std::max(label_kdegree + 5, label_inv_steps + minidx_t + 7) * logp)
    throw std::runtime_error(
        "Modulus not large enough for point labelling after bootstrap.");

  if (dusts.logq <
      std::max(label_kdegree + 5, label_inv_steps + minidx_t + 7) * logp)
    throw std::runtime_error(
        "Modulus not large enough for point labelling after bootstrap.");

  // `pointLabelingSIMD` expects `dusts` and `points` to have same modulus.
  std::vector<Ciphertext> mod_points(points.size());
  for (size_t i = 0; i < points.size(); ++i) {
    scheme.modDownTo(mod_points[i], points[i], dusts.logq);
  }

  pointLabelingSIMD(res, dusts, mod_points, dim, npoints, num_dusts,
                    label_kdegree, label_inv_steps, minidx_t);
}
