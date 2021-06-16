// MIT License
//
// Copyright (c) 2021 Oliver Schick
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

#include <iostream>
#include <iomanip>

#include <mpo19/clustering.hpp>
#include <mpo19/typedefs.hpp>
#include <mpo19/point.hpp>
#include <mpo19/linkage.hpp>
#include <mpo19/python.hpp>
#include <boost/hana.hpp>
#include <numeric>
#include <utility>
#include <iostream>
#include <iomanip>
#include <random>
#include <unordered_set>
#include <sstream>

#define BOOST_TEST_MODULE clustering_test
#include <boost/test/included/unit_test.hpp>
#include <boost/test/data/test_case.hpp>
#include <boost/test/data/monomorphic.hpp>

using namespace mpo19;
namespace bdata = boost::unit_test::data;

constexpr int m = 1'000'000;

BOOST_DATA_TEST_CASE(clustering_test
                     , bdata::random(0, m) ^ bdata::xrange(10)
                     , seed, idx)
{
    (void) idx;
    constexpr size_t n = 101;
    size_t t = 1;
    std::mt19937 mt(seed);
    std::uniform_int_distribution<int> distr{1, m};

    static_assert(n*(n-1) < m, "");

    std::unordered_set<size_t> randoms;

    while(randoms.size() < n*(n-1)/2) {
        randoms.insert(distr(mt));
    }

    auto rand_it = randoms.begin();
    table_t<size_t> precomputed{boost::extents[n][n]};
    for(index_t i = 0; i != index_t(n); ++i) {
        for(index_t j = 0; j != index_t(n); ++j) {
            if(i > j) continue;
            if(i == j) {
                precomputed[i][j] = 0;
            } else {
                precomputed[i][j] = *rand_it;
                precomputed[j][i] = precomputed[i][j];
                ++rand_it;
            }
        }
    }
    
    std::stringstream ss;

    auto python_dendrogram = python{}.hierarchical_clustering(precomputed, t, "single");
    for(auto const& c : python_dendrogram) {
        ss << c;
    }
    std::string python_dendrogram_str = ss.str();
    ss.str(std::string{});

    hierarchical_clustering<opt_min_linkage<table_t<size_t>>> hc{precomputed};
    auto hc_dendrogram = hc(t);
    for(auto const& c : hc_dendrogram) {
        ss << c;
    }
    std::string hc_dendrogram_str = ss.str();
    ss.str(std::string{});
    
    BOOST_TEST(python_dendrogram_str == hc_dendrogram_str);

}
