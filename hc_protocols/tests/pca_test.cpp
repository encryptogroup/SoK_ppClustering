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

#include <pca_test.hpp>
#include <mpo19/pca.hpp>
#include <mpo19/typedefs.hpp>
#include <mpo19/point.hpp>
#include <mpo19/random.hpp>
#include <mpo19/cluster.hpp>
#include <mpo19/clustering.hpp>
#include <mpo19/python.hpp>
#include <mpo19/linkage.hpp>
#include <mpo19/config.hpp>

#include <vector>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <chrono>
#include <thread>
#include <numeric>
#include <unordered_set>
#include <exception>

#include <boost/range/algorithm_ext.hpp>
#include <boost/range/adaptors.hpp>

#define BOOST_TEST_MODULE pca_test
#include <boost/test/included/unit_test.hpp>
#include <boost/test/data/test_case.hpp>
#include <boost/test/data/monomorphic.hpp>

table_t<size_t> pca_internal_state::Sigma;

using namespace mpo19;
namespace bdata = boost::unit_test::data;

constexpr int m = 10'000;
constexpr int min_seed = 0;
constexpr int max_seed = 1'000'000;

template<typename RandomDistribution>
point<int64_t> random_point(size_t point_dimension, std::mt19937& mt, RandomDistribution&& rand)
{
    point<int64_t> p(point_dimension);
    for(size_t i = 0; i != point_dimension; ++i) {
        p[i] = rand(mt);
    }
    return p;
}

template<typename Dist>
index_t check_unique_dist_matrix(std::vector<point<int64_t>>& input, Dist&& dist)
{
    std::unordered_set<size_t> distances;
    size_t n = input.size();
    for(index_t i = 0; i != index_t(n); ++i) {
        for(index_t j = i + 1; j != index_t(n); ++j) {
            size_t d = dist(input[i], input[j]);
            if(!distances.insert(d).second || d == 0) {
                return j;
            }
        }
    }
    return n;
}

template<typename RandomDistribution, typename Dist>
std::vector<point<int64_t>> generate_random_no_tie_input(
                             RandomDistribution&& rand, std::seed_seq seed
                             , size_t n, size_t point_dimension, Dist&& dist)
{
    constexpr size_t max_trials = 1'000;
    std::unordered_set<size_t> randoms;
    std::vector<point<int64_t>> input;
    input.reserve(n);
    std::mt19937 mt(seed);

    for(index_t i = 0; i != index_t(n); ++i) {
        input.emplace_back(random_point(point_dimension, mt, rand));
    }

    int counter = 0;
    while(true) {
        index_t check = check_unique_dist_matrix(input, dist);
        if(check == index_t(n)) break;
        input[check] = random_point(point_dimension, mt, rand);
        if(counter == max_trials) {
            throw std::runtime_error{"Could not generate random matrix."};
        }
        ++counter;
    }

    return input;
}

bool operator==(std::vector<cluster> const& lhs, std::vector<cluster> const& rhs)
{
    std::stringstream ss;
    for(auto const& c : lhs) {
        ss << c;
    }
    std::string lhs_str = ss.str();
    ss.str(std::string{});

    for(auto const& c : rhs) {
        ss << c;
    }

    std::string rhs_str = ss.str();
    ss.str(std::string{});

    return lhs_str == rhs_str;
}

size_t dist(point<int64_t> const& p1, point<int64_t> const& p2)
{
    return apply_on_elements([](int64_t element) {
        return size_t(element * element);
    }, p1 - p2).fold();
}

template<typename T>
bool is_no_tie_dist_matrix(table_t<T>& dist_matrix)
{
    std::unordered_set<size_t> distances;
    size_t n = dist_matrix.size();
    dist_matrix.reindex(0);
    for(index_t i = 0; i != index_t(n); ++i) {
        for(index_t j = 0; j != index_t(n); ++j) {
            size_t d = dist_matrix[i][j];
            if(i == j) {
                assert(d == 0);
                continue;
            } else {
                assert(d != 0);
            }
            distances.insert(d);
        }
    }
    return distances.size() == n*(n-1)/2;
}

template<typename T>
std::ostream& operator<<(std::ostream& os, table_t<T> const& tbl)
{

    for(auto&& row : tbl) {
        os << "(";
        for(auto&& e : row) {
            os << std::setw(6) << e;
        }
        os << ")";
        os << "\n";
    }
    return os;
}

BOOST_DATA_TEST_CASE(pcure0_test
                     , bdata::random(min_seed, max_seed) ^ bdata::xrange(1)
                     , seed, idx)
{
    try {
        using namespace boost::python;
        using boost::adaptors::sliced;
        using namespace std::chrono;
        using namespace std::literals::chrono_literals;

        size_t n = 80;

        std::vector<point<int64_t>> input =
                                     generate_random_no_tie_input(
                                         std::uniform_int_distribution<int64_t> {-m, m}
                                         , {seed}, n, 3, dist
                                     );

        std::stringstream ss;

        size_t s = n/4;
        size_t p = 2;
        size_t q = 3;
        size_t t1 = 3;
        size_t t2 = 5;
        size_t t = 2;

        std::cout << "PCURE0: " << std::endl;

        do {
            try {
                std::cout << "seed=" << seed << ", idx=" << idx << std::endl;
                std::thread other_party{[&]{
                        P1 p1{7766};
                        try
                        {
                            p1.pcure0(input | sliced(input.size()/2, input.size()), s, p, q, t1, t2, t);
                        } catch(mpo19::aby_error const& e)
                        {
                            return;
                        }
                    }};
                other_party.detach();


                P2 p2{"localhost", 7766};
                auto start = steady_clock::now();
                auto out = p2.pcure0(input | sliced(0, input.size()/2), s, p, q, t1, t2, t);
                auto end = steady_clock::now();
                std::cout << "\ntook: " << duration<double> {end - start} .count() << " s\n" << std::endl;

                std::cout << "\n";
                for(auto const& p : input) {
                    std::cout << p << " ";
                }
                std::cout << std::endl;

                std::cout << "\n";
                for(auto const& c : out) {
                    std::cout << c << "\n";
                }
                std::cout << std::endl;

                break;
            } catch(mpo19::aby_error const& e) {
                std::cout << "Caught an aby error:\n";
                std::cout << e.what() << "\n";
                std::cout << "Retrying..." << std::endl;
            }
        } while(true);


    } catch(boost::python::error_already_set const&) {
        PyErr_Print();
    }

}

BOOST_DATA_TEST_CASE(pcure0_opt_test
                     , bdata::random(min_seed, max_seed) ^ bdata::xrange(1)
                     , seed, idx)
{
    try {
        using namespace boost::python;
        using boost::adaptors::sliced;
        using namespace std::chrono;
        using namespace std::literals::chrono_literals;

        size_t n = 80;

        std::vector<point<int64_t>> input =
                                     generate_random_no_tie_input(
                                         std::uniform_int_distribution<int64_t> {-m, m}
                                         , {seed}, n, 3, dist
                                     );

        std::stringstream ss;

        size_t s = n/4;
        size_t p = 2;
        size_t q = 3;
        size_t t1 = 3;
        size_t t2 = 5;
        size_t t = 2;

        std::cout << "PCURE0 OPT: " << std::endl;

        do {
            try {
                std::cout << "seed=" << seed << ", idx=" << idx << std::endl;
                std::thread other_party{[&]{
                        P1 p1{7766};
                        try
                        {
                            p1.pcure0(input | sliced(input.size()/2, input.size()), s, p, q, t1, t2, t);
                        } catch(mpo19::aby_error const& e)
                        {
                            return;
                        }
                    }};
                other_party.detach();


                P2 p2{"localhost", 7766};
                auto start = steady_clock::now();
                auto out = p2.pcure0(input | sliced(0, input.size()/2), s, p, q, t1, t2, t);
                auto end = steady_clock::now();
                std::cout << "\ntook: " << duration<double> {end - start} .count() << " s\n" << std::endl;

                std::cout << "\n";
                for(auto const& p : input) {
                    std::cout << p << " ";
                }
                std::cout << std::endl;

                std::cout << "\n";
                for(auto const& c : out) {
                    std::cout << c << "\n";
                }
                std::cout << std::endl;

                break;
            } catch(mpo19::aby_error const& e) {
                std::cout << "Caught an aby error:\n";
                std::cout << e.what() << "\n";
                std::cout << "Retrying..." << std::endl;
            }
        } while(true);


    } catch(boost::python::error_already_set const&) {
        PyErr_Print();
    }
}

BOOST_DATA_TEST_CASE(pcure1_test
                     , bdata::random(min_seed, max_seed) ^ bdata::xrange(1)
                     , seed, idx)
{
    using namespace boost::python;
    using boost::adaptors::sliced;
    using namespace std::chrono;
    using namespace std::literals::chrono_literals;

    (void) idx;

    size_t n = 80;

    std::vector<point<int64_t>> input =
                                 generate_random_no_tie_input(
                                     std::uniform_int_distribution<int64_t> {-m, m}
                                     , {seed}, n, 3, dist
                                 );

    size_t s = n/4;
    size_t p = 2;
    size_t q = 3;
    size_t t1 = 3;
    size_t t2 = 5;
    size_t t = 2;

    std::cout << "PCURE1: " << std::endl;

    std::thread other_party{[&]{
            P1 p1{7766};
            try
            {
                p1.pcure1(input | sliced(input.size()/2, input.size()), s, p, q, t1, t2, t);
            } catch(mpo19::aby_error const& e)
            {
                return;
            }
        }};
    other_party.detach();


    P2 p2{"localhost", 7766};
    try {
        auto out = p2.pcure1(input | sliced(0, input.size()/2), s, p, q, t1, t2, t);
        for(cluster const& c : out) {
            std::cout << c << std::endl;
        }
    } catch(mpo19::aby_error const& e) {
        return;
    }
}

BOOST_DATA_TEST_CASE(pcure1_opt_test
                     , bdata::random(min_seed, max_seed) ^ bdata::xrange(1)
                     , seed, idx)
{
    using namespace boost::python;
    using boost::adaptors::sliced;
    using namespace std::chrono;
    using namespace std::literals::chrono_literals;

    (void) idx;

    size_t n = 80;

    std::vector<point<int64_t>> input =
                                 generate_random_no_tie_input(
                                     std::uniform_int_distribution<int64_t> {-m, m}
                                     , {seed}, n, 3, dist
                                 );

    size_t s = n/4;
    size_t p = 2;
    size_t q = 3;
    size_t t1 = 3;
    size_t t2 = 5;
    size_t t = 2;

    std::cout << "PCURE1 OPT: " << std::endl;

    std::thread other_party{[&]{
            P1 p1{7766};
            try
            {
                p1.pcure1_opt(input | sliced(input.size()/2, input.size()), s, p, q, t1, t2, t);
            } catch(mpo19::aby_error const& e)
            {
                return;
            }
        }};
    other_party.detach();


    P2 p2{"localhost", 7766};
    try {
        auto out = p2.pcure1_opt(input | sliced(0, input.size()/2), s, p, q, t1, t2, t);
        for(cluster const& c : out) {
            std::cout << c << std::endl;
        }
        std::cout << std::endl;
    } catch(mpo19::aby_error const& e) {
        return;
    }
}

BOOST_DATA_TEST_CASE(pcure2_test
                     , bdata::random(min_seed, max_seed) ^ bdata::xrange(1)
                     , seed, idx)
{
    using namespace boost::python;
    using boost::adaptors::sliced;
    using namespace std::chrono;
    using namespace std::literals::chrono_literals;

    (void) idx;

    size_t n = 80;

    std::vector<point<int64_t>> input =
                                 generate_random_no_tie_input(
                                     std::uniform_int_distribution<int64_t> {-m, m}
                                     , {seed}, n, 3, dist
                                 );

    size_t s = n/4;
    size_t q = 3;
    size_t t1 = 3;
    size_t t2 = 5;
    size_t t = 2;

    std::cout << "PCURE2: " << std::endl;

    std::thread other_party{[&]{
            P1 p1{7766};
            try
            {
                p1.pcure2(input | sliced(input.size()/2, input.size()), s, q, t1, t2, t);
            } catch(mpo19::aby_error const& e)
            {
                return;
            }
        }};
    other_party.detach();


    P2 p2{"localhost", 7766};
    try {
        auto out = p2.pcure2(input | sliced(0, input.size()/2), s, q, t1, t2, t);
        for(cluster const& c : out) {
            std::cout << c << std::endl;
        }
    } catch(mpo19::aby_error const& e) {
        return;
    }
}

BOOST_DATA_TEST_CASE(pcure2_opt_test
                     , bdata::random(min_seed, max_seed) ^ bdata::xrange(1)
                     , seed, idx)
{
    using namespace boost::python;
    using boost::adaptors::sliced;
    using namespace std::chrono;
    using namespace std::literals::chrono_literals;

    (void) idx;

    size_t n = 80;

    std::vector<point<int64_t>> input =
                                 generate_random_no_tie_input(
                                     std::uniform_int_distribution<int64_t> {-m, m}
                                     , {seed}, n, 3, dist
                                 );

    size_t s = n/4;
    size_t q = 3;
    size_t t1 = 3;
    size_t t2 = 5;
    size_t t = 2;

    std::cout << "PCURE2 OPT: " << std::endl;

    std::thread other_party{[&]{
            P1 p1{7766};
            try
            {
                p1.pcure2_opt(input | sliced(input.size()/2, input.size()), s, q, t1, t2, t);
            } catch(mpo19::aby_error const& e)
            {
                return;
            }
        }};
    other_party.detach();


    P2 p2{"localhost", 7766};
    try {
        auto out = p2.pcure2_opt(input | sliced(0, input.size()/2), s, q, t1, t2, t);
        for(cluster const& c : out) {
            std::cout << c << std::endl;
        }
        std::cout << std::endl;
    } catch(mpo19::aby_error const& e) {
        return;
    }
}

constexpr size_t n = 11;
constexpr size_t t = 1;

BOOST_DATA_TEST_CASE(opt_test
                     , bdata::random(min_seed, max_seed) ^ bdata::xrange(10)
                     , seed, idx)
{
    try {
        using namespace boost::python;
        using boost::adaptors::sliced;
        using namespace std::chrono;
        using namespace std::literals::chrono_literals;


        std::vector<point<int64_t>> input =
                                     generate_random_no_tie_input(
                                         std::uniform_int_distribution<int64_t> {-m, m}
                                         , {seed}, n, 3, dist
                                     );

        std::stringstream ss;

        do {
            try {
                std::cout << "seed=" << seed << ", idx=" << idx << std::endl;
                std::thread other_party{[&]{
                        P1 p1{7766};
                        try
                        {
                            p1.execute_opt(input | sliced(input.size()/2, input.size()), t);
                            throw aby_error{"test"};
                        } catch(mpo19::aby_error const& e)
                        {
                            return;
                        }
                    }};
                other_party.detach();


                P2 p2{"localhost", 7766};
                auto start = steady_clock::now();
                auto out = p2.execute_opt(input | sliced(0, input.size()/2), t);
                auto end = steady_clock::now();
                std::cout << "\ntook: " << duration<double> {end - start} .count() << " s\n" << std::endl;

                BOOST_TEST(is_no_tie_dist_matrix(pca_internal_state::Sigma));
                BOOST_TEST(std::get<0>(out).size() == t);
                BOOST_TEST(std::get<1>(out).size() == t);
                BOOST_TEST(std::get<2>(out).size() == t);
                BOOST_TEST(python{} .hierarchical_clustering(pca_internal_state::Sigma, t, "single") == std::get<0>(out));
                break;
            } catch(mpo19::aby_error const& e) {
                std::cout << "Caught an aby error:\n";
                std::cout << e.what() << "\n";
                std::cout << "Retrying..." << std::endl;
            }
        } while(true);


    } catch(boost::python::error_already_set const&) {
        PyErr_Print();
    }
}

BOOST_DATA_TEST_CASE(random_test
                     , bdata::random(min_seed, max_seed) ^ bdata::xrange(10)
                     , seed, idx)
{
    try {
        using namespace boost::python;
        using boost::adaptors::sliced;
        using namespace std::chrono;
        using namespace std::literals::chrono_literals;


        std::vector<point<int64_t>> input =
                                     generate_random_no_tie_input(
                                         std::uniform_int_distribution<int64_t> {-m, m}
                                         , {seed}, n, 3, dist
                                     );

        std::stringstream ss;

        do {
            try {
                std::cout << "seed=" << seed << ", idx=" << idx << std::endl;
                std::thread other_party{[&]{
                        P1 p1{7766};
                        try
                        {
                            p1.execute_pca(input | sliced(input.size()/2, input.size()), t);
                            throw aby_error{"test"};
                        } catch(mpo19::aby_error const& e)
                        {
                            return;
                        }
                    }};
                other_party.detach();


                P2 p2{"localhost", 7766};
                auto start = steady_clock::now();
                auto out = p2.execute_pca(input | sliced(0, input.size()/2), t);
                auto end = steady_clock::now();
                std::cout << "\ntook: " << duration<double> {end - start} .count() << " s\n" << std::endl;

                BOOST_TEST(is_no_tie_dist_matrix(pca_internal_state::Sigma));
                BOOST_TEST(std::get<0>(out).size() == t);
                BOOST_TEST(std::get<1>(out).size() == t);
                BOOST_TEST(std::get<2>(out).size() == t);
                BOOST_TEST(python{} .hierarchical_clustering(pca_internal_state::Sigma, t, "complete") == std::get<0>(out));
                break;
            } catch(mpo19::aby_error const& e) {
                std::cout << "Caught an aby error:\n";
                std::cout << e.what() << "\n";
                std::cout << "Retrying..." << std::endl;
            }
        } while(true);


    } catch(boost::python::error_already_set const&) {
        PyErr_Print();
    }
}

/*
BOOST_AUTO_TEST_CASE(iris_test)
{
    using namespace boost::python;
    using boost::adaptors::sliced;
    using namespace std::chrono;
    using namespace std::literals::chrono_literals;

    constexpr size_t t = 3;

    try {
        python py;
        std::vector<point<int>> iris = py.load_iris();
        size_t n = iris.size();

        for(auto const& p : iris)
            std::cout << p << ", ";
        std::cout << std::endl;

        random_permutation p(iris.size());

        std::thread other_party{[&]{
                P1 p1{7766};
                p1.execute_pca(iris | sliced(iris.size()/2, iris.size()), t);
            }};
        other_party.detach();

        P2 p2{"localhost", 7766};

        auto start = steady_clock::now();
        auto out = p2.execute_pca(iris | sliced(0, iris.size()/2), t);
        auto end = steady_clock::now();

        std::cout << "\ntook: " << duration<double> {end - start} .count() << " s\n" << std::endl;

        std::vector<cluster>& pca_dendrogram = std::get<0>(out);

        BOOST_TEST(pca_dendrogram.size() == 3);

        table_t<size_t> precomputed{boost::extents[n][n]};
        for(size_t i = 0; i != n; ++i) {
            for(size_t j = 0; j != n; ++j) {
                precomputed[i][j] = static_cast<size_t>(pca_internal_state::Sigma[i][j]);
            }
        }

        hierarchical_clustering<
        linkage_precomputed<
        table_t<size_t>, max_linkage_precomputed
        > > hc{precomputed};

        BOOST_TEST(hc(t) == pca_dendrogram);

        for(point<double> const& p : std::get<1>(out)) {
            std::cout << p << " ";
        }
        std::cout << std::endl;

        for(int i : std::get<2>(out)) {
            std::cout << i << ", ";
        }
        std::cout << std::endl;

    } catch(boost::python::error_already_set const&) {
        PyErr_Print();
    }
}
*/
