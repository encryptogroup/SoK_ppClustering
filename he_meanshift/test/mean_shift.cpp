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

#include <HEAAN.h>
#include <limits.h>
#include <stdlib.h>

#include <boost/test/data/test_case.hpp>
#include <boost/test/unit_test.hpp>
#include <ckp19/dataset.hpp>
#include <ckp19/mean_shift.hpp>
#include <ckp19/plaintext_clustering.hpp>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <vector>

#include "common.hpp"

namespace utf = boost::unit_test;

// Test fixture
struct MeanShiftFixture {
  Ring ring;
  SecretKey secretKey;
  Scheme scheme;
  MeanShift meanshift;

  MeanShiftFixture()
      : secretKey(ring),
        scheme(secretKey, ring),
        meanshift(scheme, SchemeConstants.logn, SchemeConstants.logq,
                  SchemeConstants.logq_boot, SchemeConstants.logp,
                  SchemeConstants.logt, SchemeConstants.seed) {
    BOOST_TEST_REQUIRE(SchemeConstants.logn < logN - 1);

    scheme.addLeftRotKeys(secretKey);
    scheme.addRightRotKeys(secretKey);

    // Tests will involve bootstrap operations on ciphertext with 2 or 4 slots.
    scheme.addBootKey(secretKey, 3, SchemeConstants.logq_boot + 4);
    scheme.addBootKey(secretKey, 4, SchemeConstants.logq_boot + 4);
  }
};

// Utility functions

void scale_vector(double *v, long n, double r) {
  double sum;
  for (int i = 0; i < n; ++i) sum += v[i] * v[i];

  sum = std::sqrt(sum);
  sum /= r;

  for (int i = 0; i < n; ++i) v[i] /= sum;
};

void getDustCiphertextFromPlaintext(
    Ciphertext &dusts, MeanShift &meanshift,
    const std::vector<std::vector<double>> &points, long dim, long npoints,
    long num_dusts, long slots) {
  std::vector<complex<double>> pdusts(slots, 0);
  for (long i = 0; i < num_dusts; ++i) {
    for (long j = 0; j < npoints; ++j) {
      for (long k = 0; k < dim; ++k) {
        pdusts[i * npoints * dim + j * dim + k] = points[i % points.size()][k];
      }
    }
  }

  meanshift.scheme.encrypt(dusts, pdusts.data(), slots, meanshift.logp,
                           meanshift.logq);
}

// Since most of these algorithms take signfificant time to execute if `logN`
// is large, the tests are quite minimal and don't verify across a lot of
// different parameter values.
BOOST_FIXTURE_TEST_SUITE(mean_shift, MeanShiftFixture)

BOOST_AUTO_TEST_CASE(inverse, *utf::tolerance(1e-4)) {
  srand(SchemeConstants.seed);
  meanshift.refreshSeed(SchemeConstants.seed);

  double m = 20;
  double steps = 11;

  std::unique_ptr<double[]> pval(
      EvaluatorUtils::randomRealArray(meanshift.slots, m));

  Ciphertext cval;
  scheme.encrypt(cval, pval.get(), meanshift.slots, meanshift.logp,
                 meanshift.logq);
  Ciphertext cinv;
  meanshift.inverse(cinv, cval, m, steps);

  BOOST_TEST(cinv.logq == cval.logq - (steps + 2) * meanshift.logp,
             "cinv.logq: " << cinv.logq << ", cval.logq: " << cval.logq);
  BOOST_TEST(cinv.logp == cval.logp,
             "cinv.logp: " << cinv.logp << ", cval.logp: " << cval.logp);

  std::unique_ptr<complex<double>[]> pinv(scheme.decrypt(secretKey, cinv));

  for (long i = 0; i < meanshift.slots; ++i) {
    double exp = 1.0 / pval[i];
    double output = pinv[i].real();
    BOOST_TEST(exp == output,
               "output: " << pinv[i].real() << ", expected: " << exp);
  }
}

BOOST_TEST_DECORATOR(*utf::tolerance(1e-5))
BOOST_DATA_TEST_CASE(kernelSIMD, utf::data::xrange<long>(2, 4), dim) {
  srand(SchemeConstants.seed);
  meanshift.refreshSeed(SchemeConstants.seed);

  long degree = 5;
  // Deliberate choice to not make `npoints` a power of 2. Having buffer values
  // in plaintext slots can help catch errors in implementation.
  long npoints = 6;

  // `npoints` points should be encrypted within the same ciphertext.
  BOOST_TEST_REQUIRE(npoints * dim < meanshift.slots);

  long req_slots = npoints * dim;
  std::unique_ptr<double[]> x(EvaluatorUtils::randomRealArray(meanshift.slots));
  std::unique_ptr<double[]> y(EvaluatorUtils::randomRealArray(meanshift.slots));
  scale_vector(x.get(), req_slots, 0.5);
  scale_vector(y.get(), req_slots, 0.5);

  Ciphertext cx, cy;
  scheme.encrypt(cx, x.get(), meanshift.slots, meanshift.logp, meanshift.logq);
  scheme.encrypt(cy, y.get(), meanshift.slots, meanshift.logp, meanshift.logq);

  Ciphertext ckernel;
  meanshift.kernelSIMD(ckernel, cx, cy, dim, npoints, degree);

  BOOST_TEST(ckernel.logq == cx.logq - (degree + 2) * meanshift.logp,
             "ckernel.logq: " << ckernel.logq << ", cx.logq: " << cx.logq);
  BOOST_TEST(ckernel.logp == cx.logp,
             "ckernel.logp: " << ckernel.logp << ", cx.logp: " << cx.logp);

  std::unique_ptr<complex<double>[]> output(
      meanshift.scheme.decrypt(secretKey, ckernel));

  PlainMeanShift pmeanshift(meanshift.seed);
  for (long i = 0; i < npoints; ++i) {
    double kernel = pmeanshift.kernel(
        std::vector<double>(x.get() + i * dim, x.get() + (i + 1) * dim),
        std::vector<double>(y.get() + i * dim, y.get() + (i + 1) * dim),
        degree);

    for (long j = 0; j < dim; ++j) {
      BOOST_TEST(
          output[i * dim + j].real() == kernel,
          "output: " << output[i * dim + j].real() << ", expected: " << kernel);
    }
  }
}

BOOST_TEST_DECORATOR(*utf::tolerance(5e-4))
BOOST_DATA_TEST_CASE(refreshDusts, utf::data::xrange(2, 4), dim) {
  srand(SchemeConstants.seed);
  meanshift.refreshSeed(SchemeConstants.seed);

  long num_dusts = 4;
  // Compute `npoints` (data points per ciphertext) based on given params
  long npoints = meanshift.slots / (dim * num_dusts);

  BOOST_TEST_REQUIRE(npoints >= 1);

  Dataset ds(dim, npoints, SchemeConstants.seed);
  Ciphertext cdusts;
  getDustCiphertextFromPlaintext(cdusts, meanshift, ds.getPoints(), dim,
                                 npoints, num_dusts, meanshift.slots);
  std::unique_ptr<complex<double>[]> orig(scheme.decrypt(secretKey, cdusts));

  meanshift.refreshDusts(cdusts, dim, npoints, num_dusts);
  std::unique_ptr<complex<double>[]> bootstrap(
      scheme.decrypt(secretKey, cdusts));

  BOOST_TEST(
      cdusts.logq > cdusts.logp,
      "cdusts.logq: " << cdusts.logq << ", cdusts.logp: " << cdusts.logp);
  BOOST_TEST(cdusts.logq > SchemeConstants.logq_boot,
             "cdusts.logq: " << cdusts.logq
                             << ", logq_boot: " << SchemeConstants.logq_boot);

  for (long i = 0; i < num_dusts; ++i) {
    for (long j = 0; j < dim; ++j) {
      // `refreshDusts` sets only the first instance of the dust
      long pos = i * dim * npoints + j;
      BOOST_TEST(orig[pos].real() == bootstrap[pos].real(),
                 "output: " << bootstrap[pos].real()
                            << ", expected: " << orig[pos].real());
    }
  }
}

BOOST_TEST_DECORATOR(*utf::tolerance(5e-4))
BOOST_DATA_TEST_CASE(modeSeekingSIMD,
                     utf::data::xrange(2, 4) *
                         utf::data::make(std::vector<long>{1, 5}),
                     dim, steps) {
  srand(SchemeConstants.seed);
  meanshift.refreshSeed(SchemeConstants.seed);

  long num_point_cts = 2;
  long num_dusts = 4;
  long kdegree = 5;
  long inv_steps = 8;

  // Compute `npoints` (data points per ciphertext) based on given params
  long npoints = meanshift.slots / (dim * num_dusts);

  BOOST_TEST_REQUIRE(npoints >= 1);

  Dataset ds(dim, num_point_cts * npoints, SchemeConstants.seed);
  ds.rescaleDatasetInPlace(0.5);

  EncryptedDataset eds;
  ds.encryptSIMD(eds, scheme, num_dusts, meanshift.slots, meanshift.logp,
                 meanshift.logq);

  // Meanshift object should ideally be constructed with `eds.log_slots`. For
  // testing, we work our way backwards and set `eds.log_slots` to be equal to
  // `meanshift.logn`.
  BOOST_TEST_REQUIRE(eds.log_slots == meanshift.logn);

  std::vector<std::vector<double>> dusts;
  PlainMeanShift pmeanshift(meanshift.seed);
  pmeanshift.modeSeeking(dusts, ds.getPoints(), num_dusts, steps, kdegree,
                         inv_steps);

  Ciphertext cdusts;
  meanshift.modeSeekingSIMD(cdusts, eds.points, eds.dim, eds.npoints, num_dusts,
                            steps, kdegree, inv_steps);

  if (steps == 1) {
    BOOST_TEST(cdusts.logq == (eds.points[0].logq -
                               (kdegree + inv_steps + 7) * meanshift.logp),
               "cdusts.logq: " << cdusts.logq
                               << ", dataset.logq: " << eds.points[0].logq);
    BOOST_TEST(cdusts.logp == eds.points[0].logp,
               "cdusts.logp: " << cdusts.logp
                               << ", dataset.logp: " << eds.points[0].logp);
  } else {
    BOOST_TEST(
        cdusts.logq > cdusts.logp,
        "cdusts.logq: " << cdusts.logq << ", cdusts.logp: " << cdusts.logp);
    BOOST_TEST(cdusts.logq > SchemeConstants.logq_boot,
               "cdusts.logq: " << cdusts.logq
                               << ", logq_boot: " << SchemeConstants.logq_boot);
  }

  std::unique_ptr<complex<double>[]> pdusts(
      meanshift.scheme.decrypt(secretKey, cdusts));

  for (long i = 0; i < num_dusts; ++i) {
    for (long j = 0; j < npoints; ++j) {
      for (long k = 0; k < dim; ++k) {
        long pos = i * npoints * dim + j * dim + k;
        BOOST_TEST(
            pdusts[pos].real() == dusts[i][k],
            "output: " << pdusts[pos].real() << ", expected: " << dusts[i][k]);
      }
    }
  }
}

BOOST_TEST_DECORATOR(*utf::tolerance(5e-1))
BOOST_DATA_TEST_CASE(pointLabelingSIMD, utf::data::xrange(2, 4), dim) {
  srand(SchemeConstants.seed);
  meanshift.refreshSeed(SchemeConstants.seed);

  long num_point_cts = 2;
  long num_dusts = 4;
  long kdegree = 4;
  long inv_steps = 13;
  long minidx_t = 4;
  long npoints = meanshift.slots / (dim * num_dusts);

  BOOST_TEST_REQUIRE(npoints >= 1);

  Dataset ds(dim, num_point_cts * npoints, SchemeConstants.seed);
  ds.rescaleDatasetInPlace(0.4);

  EncryptedDataset eds;
  ds.encryptSIMD(eds, scheme, num_dusts, meanshift.slots, meanshift.logp,
                 meanshift.logq);

  // Meanshift object should ideally be constructed with `eds.log_slots`. For
  // testing, we work our way backwards and set `eds.log_slots` to be equal to
  // `meanshift.logn`.
  BOOST_TEST_REQUIRE(eds.log_slots == meanshift.logn);

  std::vector<std::vector<double>> dusts(ds.getPoints().begin(),
                                         ds.getPoints().begin() + num_dusts);

  Ciphertext cdusts;
  getDustCiphertextFromPlaintext(cdusts, meanshift, dusts, dim, npoints,
                                 num_dusts, meanshift.slots);

  BOOST_TEST_REQUIRE(cdusts.logq == eds.points[0].logq);
  BOOST_TEST_REQUIRE(cdusts.logp == eds.points[0].logp);

  std::vector<std::vector<double>> labels;
  PlainMeanShift pmeanshift(meanshift.seed);
  pmeanshift.pointLabeling(labels, ds.getPoints(), dusts, kdegree, inv_steps,
                           minidx_t);

  std::vector<Ciphertext> clabels;
  meanshift.pointLabelingSIMD(clabels, cdusts, eds.points, dim, npoints,
                              num_dusts, kdegree, inv_steps, minidx_t);

  for (long p = 0; p < num_point_cts; ++p) {
    std::unique_ptr<complex<double>[]> plabel(
        meanshift.scheme.decrypt(secretKey, clabels[p]));

    for (long i = 0; i < num_dusts; ++i) {
      for (long j = 0; j < npoints; ++j) {
        for (long k = 0; k < dim; ++k) {
          double output = plabel[i * npoints * dim + j * dim + k].real();
          double exp = labels[p * npoints + j][i];
          BOOST_TEST(output == exp,
                     "output: " << output << ", expected: " << exp);
        }
      }
    }
  }
}

BOOST_AUTO_TEST_SUITE_END()
