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

#include <mpo19/utility.hpp>
#include <mpo19/random.hpp>
#include <algorithm>
#include <boost/range/algorithm.hpp>

#define BOOST_TEST_MODULE point_test
#include <boost/test/included/unit_test.hpp>
#include <boost/test/data/test_case.hpp>
#include <boost/test/data/monomorphic.hpp>
#include <boost/multi_array.hpp>

namespace bdata = boost::unit_test::data;

BOOST_AUTO_TEST_CASE(buffered_rng_test)
{
    using rng_t = buffered_secure_rng;
    size_t test_sz = 2*rng_t::default_buffer_size + 1;
    rng_t rng_empty
    , rng_default_size(rng_t::default_buffer_size)
    , rng_filled(test_sz);

    BOOST_TEST(rng_default_size.remaining() == rng_t::default_buffer_size);
    BOOST_TEST(rng_filled.remaining() == test_sz);
    BOOST_TEST(rng_empty.remaining() == 0);

    rng_t::result_type r1 = rng_empty()
                            , r2 = rng_default_size() + rng_default_size()
                            , r3 = rng_filled();

    BOOST_TEST(rng_empty.remaining() == rng_t::default_buffer_size - 1);
    BOOST_TEST(rng_default_size.remaining() == rng_t::default_buffer_size - 2);
    BOOST_TEST(rng_filled.remaining() == test_sz - 1);

    rng_t::result_pointer rs = rng_filled.get_rand_numbers(test_sz-3);
    BOOST_TEST(rng_filled.remaining() == 2);

    rs = rng_filled.get_rand_numbers(3*test_sz);
    //Buffer should be filled according to size and then completely consumed
    BOOST_TEST(rng_filled.remaining() == 0);
    for(int i = 0; i < 1000; ++i) {
        rs = rng_filled.get_rand_numbers(3*test_sz);
        std::vector<rng_t::result_type> tmp;
        tmp.reserve(3*test_sz);
        std::copy(rs, rs + 3*test_sz, std::back_inserter(tmp));

        //tmp holds a copy of rs_dash, so it must be equal
        BOOST_TEST(std::equal(rs, rs + 3*test_sz, tmp.cbegin()));

        buffered_secure_rng::buffer_t rs_copy = rng_filled.move_buffer();
        BOOST_TEST(rng_filled.remaining() == 0);

        rng_t::result_pointer rs_dash = rng_filled.get_rand_numbers(3*test_sz);
        BOOST_TEST(rng_filled.remaining() == 0);

        //Probabiliy of rs == rs_dash is so small, that test should never fail
        //if rs and rs_dash are really random.
        BOOST_TEST(!std::equal(rs_dash, rs_dash + 3*test_sz, tmp.cbegin()));
        //rs_copy holds the buffer before it got reallocated, so it must be equal to tmp
        BOOST_TEST(std::equal(rs_copy.begin(), rs_copy.end(), tmp.cbegin()));
        std::pair<rng_t::result_pointer, size_t> rbits;
        rbits = rng_filled.get_rand_bits(41);
        BOOST_TEST(rbits.first[rbits.second - 1] < 2);

        rbits = rng_filled.get_rand_bits(42);
        BOOST_TEST(rbits.first[rbits.second - 1] < 4);

        rbits = rng_filled.get_rand_bits(43);
        BOOST_TEST(rbits.first[rbits.second - 1] < 8);

        rbits = rng_filled.get_rand_bits(44);
        BOOST_TEST(rbits.first[rbits.second - 1] < 16);

        rbits = rng_filled.get_rand_bits(45);
        BOOST_TEST(rbits.first[rbits.second - 1] < 32);

        rbits = rng_filled.get_rand_bits(46);
        BOOST_TEST(rbits.first[rbits.second - 1] < 64);

        rbits = rng_filled.get_rand_bits(47);
        BOOST_TEST(rbits.first[rbits.second - 1] < 128);
    }
    //Remove unused variable warnings
    (void) r1;
    (void) r2;
    (void) r3;
}

BOOST_AUTO_TEST_CASE(random_permutation_test)
{
    using rng_t = buffered_secure_rng;
    size_t n = 100;
    random_permutation const pi_1(n);

    std::vector<size_t> m(n);
    std::iota(m.begin(), m.end(), 0);

    auto mapping = pi_1.permute(m);

    //=================================================================
    // Test 1-dimensional case
    //=================================================================

    std::vector<int> v, v_copy;
    v.reserve(n);
    std::generate_n(std::back_inserter(v), n, rng_t{n});
    boost::range::copy(v, std::back_inserter(v_copy));
    auto permuted_v = pi_1.permute(v);

    BOOST_TEST(std::is_permutation(v.begin(), v.end(), permuted_v.begin(), permuted_v.end()));
    //Probability of v and permuted_v being equal is negligible, if permutation is random
    BOOST_TEST(!std::equal(v.begin(), v.end(), permuted_v.begin()));

    //Check Consistency with mapping
    for(size_t i = 0; i != n; ++i) {
        BOOST_TEST( permuted_v(i) == v[mapping(i)]);
    }

    //Same permutation applied on identical ranges
    BOOST_TEST(std::equal(permuted_v.begin(), permuted_v.end(), pi_1.permute(v_copy).begin()));
    
    //permute_inplace has the same effect as permute
    pi_1.permute_inplace(v_copy);
    BOOST_TEST(std::equal(permuted_v.begin(), permuted_v.end(), v_copy.begin()));
    

    //=================================================================
    // Test 2-dimensional case
    //=================================================================

    rng_t rng;

    boost::multi_array<int, 2> a{boost::extents[n][n]};
    for(size_t i = 0; i != n; ++i) {
        for(size_t j = 0; j != n; ++j) {
            a[i][j] = rng();
        }
    }

    auto permuted_a = pi_1.permute(a);
    boost::multi_array<int, 2> a_copy = a;
    
    pi_1.permute_inplace(a_copy);

    for(size_t i = 0; i != n; ++i) {
        for(size_t j = 0; j != n; ++j) {
            BOOST_TEST(permuted_a(i)(j) == a[mapping(i)][mapping(j)]);
            BOOST_TEST(permuted_a(i)(j) == a_copy[i][j]);
        }
    }

    //=================================================================
    // Test 3-dimensional case
    //=================================================================

    boost::multi_array<int, 3> b{boost::extents[n][n][n]};
    for(size_t i = 0; i != n; ++i) {
        for(size_t j = 0; j != n; ++j) {
            for(size_t k = 0; k != n; ++k) {
                b[i][j][k] = rng();
            }
        }
    }

    auto permuted_b = pi_1.permute(b);
    boost::multi_array<int, 3> b_copy = b;

    pi_1.permute_inplace(b_copy);

    for(size_t i = 0; i != n; ++i) {
        for(size_t j = 0; j != n; ++j) {
            for(size_t k = 0; k != n; ++k) {
                BOOST_TEST(permuted_b(i)(j)(k) == b[mapping(i)][mapping(j)][mapping(k)]);
                BOOST_TEST(permuted_b(i)(j)(k) == b_copy[i][j][k]);
            }
        }
    }

}
