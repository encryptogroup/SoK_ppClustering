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

#include <random>
#include <vector>

// Perform secure mean shift clustering.
class MeanShift {
 public:
  long logp;  // Log of noise parameter
  long logn;  // Log of total usable slots
  long logq;  // Log of ciphertext modulus
  long logt;
  long logq_boot;  // Log of modulus to perform bootstrap
  long slots;      // Total usable slots. Equals 1 << logn.
  long seed;       // Random seed used for choosing dusts.
  Scheme &scheme;  // Reference to scheme for performing HE operations

  MeanShift(Scheme &scheme, long logn, long logq, long logq_boot, long logp,
            long logt);
  MeanShift(Scheme &scheme, long logn, long logq, long logq_boot, long logp,
            long logt, long seed);

  void refreshSeed(long seed);
  void refreshSeed();

  // Using HEAAN's `addRightRotKeys`/`addLeftRotKeys` only allows
  // rotation in powers of 2. This allows rotation by arbitrary amounts by
  // decomposing it into powers of 2.
  void leftRotateAndEqual(Ciphertext &x, long r);
  void rightRotateAndEqual(Ciphertext &x, long r);

  // Utility functions to implement SIMD operations

  // Sum `num` slots, each of size `step` into the first, left-most set of
  // slots. This is done in a SIMD manner.
  void sumSlotsLeft(Ciphertext &x, long num, long step);

  // Clear slots towards the right. Retain `ret_size` slots for every `skip`
  // slots until `max_slot`, after which everything is cleared.
  void clearSlotsRight(Ciphertext &x, long ret_size, long skip, long max_slot);

  // Assumes that the slots to which the value has to be copied have
  // already been cleared (i.e. they have a 0).
  void repeatSlotRight(Ciphertext &x, long rep, long step);

  // Compute inverse of each slot. `m` is maximum magnitude of
  // any slot.
  void inverse(Ciphertext &res, Ciphertext &x, double m, long steps);

  // The non-SIMD variants treat each point as a seperate Ciphertext with
  // `slots` number of dimensions.

  void l2NormSquared(Ciphertext &res, Ciphertext &x);

  // `num_points` is the total number of points in the ciphertext for which the
  // l2Norm has to be computed. Note that this is not the same as `npoints` used
  // with respect to dataset size in mean shift clustering.
  void l2NormSquaredSIMD(Ciphertext &res, Ciphertext &x, long dim,
                         long num_points);

  void minIdx(std::vector<Ciphertext> &res, std::vector<Ciphertext> &x, long t,
              long inv_steps);
  void minIdxSIMD(Ciphertext &res, Ciphertext &x, long dim, long npoints,
                  long t, long inv_steps, long batchsize);

  // Note that this computes the *derivative* of the original kernel within a
  // multiplicative constant (Algorithm 4). Assumes x.logq == y.logq.
  void kernel(Ciphertext &res, Ciphertext &x, Ciphertext &y, long degree);

  // `num_points` is the total number of points in the ciphertext for which the
  // kernel has to be computed. Note that this is not the same as `npoints` used
  // with respect to dataset size in mean shift clustering.
  void kernelSIMD(Ciphertext &res, Ciphertext &x, Ciphertext &y, long dim,
                  long num_points, long degree);

  // Utility function to efficiently bootstrap dusts.
  void refreshDusts(Ciphertext &dusts, long dim, long npoints, long d);

  void modeSeeking(std::vector<Ciphertext> &res,
                   std::vector<Ciphertext> &points, long d, long steps,
                   long kdegree, long inv_steps);

  // Each ciphertext in point consists of the set of npoint points repeated
  // d times i.e., (p_1, p_2, ..., p_npoints, p_1, ...)
  void modeSeekingSIMD(Ciphertext &res, std::vector<Ciphertext> &points,
                       long dim, long npoints, long d, long steps, long kdegree,
                       long inv_steps);

  // The result has the same value for each dimension for the corresponding
  // point in both variants.
  void pointLabeling(std::vector<std::vector<Ciphertext>> &res,
                     std::vector<Ciphertext> &dusts,
                     std::vector<Ciphertext> &points, long kdegree,
                     long inv_steps, long minidx_t);

  // Assumes dusts.logq = points.logq
  void pointLabelingSIMD(std::vector<Ciphertext> &res, Ciphertext &dusts,
                         std::vector<Ciphertext> &points, long dim,
                         long npoints, long d, long kdegree, long inv_steps,
                         long minidx_t);

  // Run full clustering algorithm.
  void clusterSIMD(std::vector<Ciphertext> &res,
                   std::vector<Ciphertext> &points, long dim, long npoints,
                   long num_dusts, long iterations, long mode_kdegree,
                   long label_kdegree, long mode_inv_steps,
                   long label_inv_steps, long minidx_t);

  // Create necessary keys for rotation and bootstrapping.
  void setupScheme(SecretKey &secretKey);
};
