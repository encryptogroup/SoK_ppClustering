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

#include <vector>
#include <list>
#include <map>
#include <random>
#include <functional>
#include <string>
#include <numeric>
#include <algorithm>

#define BOOST_TEST_MODULE circuits_test
#include <boost/test/included/unit_test.hpp>
#include <boost/test/data/test_case.hpp>
#include <boost/test/data/monomorphic.hpp>

namespace bdata = boost::unit_test::data;
constexpr int m = 10000;

BOOST_DATA_TEST_CASE(utility_log_test
                     , bdata::random(1, m) ^ bdata::xrange(10000)
                     , s1, idx)
{
    (void) idx;
    static_assert(ceil_log2(0) == 0, "");
    static_assert(ceil_log2(1) == 0, "");
    static_assert(floor_log2(0) == 0, "");
    static_assert(floor_log2(1) == 0, "");
    uint64_t clog = ceil_log2(s1), flog = floor_log2(s1);
    BOOST_TEST( (1u << flog) <= uint64_t(s1) );
    BOOST_TEST( (1u << clog) >= uint64_t(s1) );
    if(clog == flog) BOOST_TEST( (1u << flog) == uint64_t(s1) );
    else BOOST_TEST( (1u << clog) > uint64_t(s1) );

}

BOOST_DATA_TEST_CASE(nothrow_swap_test
                     , bdata::random(0, m) ^ bdata::random(0, m)
                     ^ bdata::xrange(10000)
                     , s1, s2, idx)
{
    (void) idx;
    auto s1_copy = s1;
    auto s2_copy = s2;

    static_assert(noexcept(nothrow_swap), "nothrow_swap throws an exception");
    nothrow_swap(s1_copy, s2_copy);
    BOOST_TEST(s1_copy == s2);
    BOOST_TEST(s2_copy == s1);

}

BOOST_DATA_TEST_CASE(utility_words_to_hold_bitlen_test
                     , bdata::random(1, m) ^ bdata::xrange(10000)
                     , s1, idx)
{
    (void) idx;

    BOOST_TEST(bytes_to_hold_bitlen(s1) == words_to_hold_bitlen(s1, CHAR_BIT));
    BOOST_TEST(bytes_to_hold_bitlen(s1 * CHAR_BIT) == s1);
    BOOST_TEST(words_to_hold_bitlen(s1 * 16, 16) == s1);
    BOOST_TEST(words_to_hold_bitlen(s1 * 32, 16) == 2*s1);
    BOOST_TEST(words_to_hold_bitlen(s1 * 16, 32) == s1/2 + s1%2);

}

BOOST_DATA_TEST_CASE(utility_code_len_test
                     , bdata::random(1, std::numeric_limits<int>::max()) ^ bdata::xrange(10000)
                     , len, idx)
{
    (void) idx;
    char buf[encode_len_buffer_size];
    size_t num_chars = encode_len(buf, len);
    size_t dec_len = decode_len([&, i = 0u](unsigned char* c, size_t len) mutable{
        assert(i + len <= encode_len_buffer_size);
        std::copy(buf + i, buf + i + len, c);
        i += len;
    });

    BOOST_TEST(dec_len == len);
    BOOST_TEST(num_chars == words_to_hold_bitlen(ceil_log2(len), CHAR_BIT - 1));
}

BOOST_DATA_TEST_CASE(utility_test
                     , bdata::random(0, m) ^ bdata::xrange(10)
                     , s1, idx)
{
    (void) idx;

    static_assert(is_range<std::vector<int>>(), "");
    static_assert(is_range<std::list<int>>(), "");
    static_assert(is_range<std::map<int, double>>(), "");
    static_assert(!is_range<int>(), "");
    static_assert(!is_range<std::map<int, double>::iterator>(), "");

    static_assert(is_serializable<int>(), "");
    static_assert(!is_serializable<std::vector<int>>(), "");

    static_assert(is_deserializable<int>(), "");
    static_assert(!is_deserializable<std::vector<int>>(), "");

    auto f = [](int, double) {};
    auto any_arg = [](auto) {};
    auto any = [](auto...) {};

    static_assert(is_invocable<decltype(f)(int, double)>::value, "");
    static_assert(!is_invocable<decltype(f)(int, double, float)>::value, "");
    static_assert(is_invocable<decltype(any_arg)(int)>::value, "");
    static_assert(is_invocable<decltype(any_arg)(double)>::value, "");
    static_assert(!is_invocable<decltype(any_arg)()>::value, "");
    static_assert(is_invocable<decltype(any)(int, double, char, float)>::value, "");
    static_assert(is_invocable<decltype(any)(double, double, std::vector<int>)>::value, "");
    static_assert(is_invocable<decltype(any)()>::value, "");

    std::vector<int> v(100 + idx), v_copy;
    std::vector<std::vector<int>> vv(10 + idx), vv_copy;

    std::mt19937 gen(s1 + idx);
    std::uniform_int_distribution<int> d{0, m};

    std::generate(v.begin(), v.end(), [&] {return d(gen);});

    for(auto& vec : vv) {
        vec.resize(10 + idx);
        std::generate(vec.begin(), vec.end(), [&] {return d(gen);});
    }

    v_copy = v;
    vv_copy = vv;

    multi_for_each(v, [&](auto& el) {
        el += s1;
    });

    std::transform(v_copy.begin(), v_copy.end(), v_copy.begin(), [&](auto const& el) {
        return el + s1;
    });

    BOOST_TEST(std::equal(v.begin(), v.end(), v_copy.begin()));

    multi_for_each(vv, vv_copy, [](auto& a, auto& b) {
        a += 2*b;
    });

    for(size_t i = 0; i != vv.size(); ++i) {
        for(size_t j = 0; j != vv[i].size(); ++j) {
            BOOST_TEST(vv[i][j]/3 == vv_copy[i][j]);
        }
    }

    BOOST_TEST(tree_accumulate(v, std::plus<> {}) == std::accumulate(v.begin(), v.end(), 0, std::plus<> {}));
    auto res1 = tree_accumulate(v, std::plus<> {}, [](auto&& a, size_t idx) {
        return 2*a + idx;
    });

    int sum = 0;
    for(size_t i = 0; i != v.size(); ++i) {
        sum += 2*v[i] + i;
    }

    BOOST_TEST(res1 == sum);
}
