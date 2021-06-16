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

#include <mpo19/paillier_wrapper.hpp>
#include <cstring>

#define BOOST_TEST_MODULE paillier_wrapper_test
#include <boost/test/included/unit_test.hpp>
#include <boost/test/data/test_case.hpp>
#include <boost/test/data/monomorphic.hpp>

namespace bdata = boost::unit_test::data;
using namespace mpo19;

constexpr int m = 100;

constexpr key_id global = "Global";

bool equal(plaintext const& lhs, plaintext const& rhs)
{
    std::vector<char> dec_lhs = lhs;
    std::vector<char> dec_rhs = rhs;
    assert(!dec_lhs.empty() && !dec_rhs.empty());
    return strcmp(dec_lhs.data(), dec_rhs.data()) == 0;
}

bool equal(paillier<global> const& lhs, paillier<global> const& rhs)
{
    return equal(decrypt(lhs), decrypt(rhs));
}

BOOST_DATA_TEST_CASE(paillier_operation_test
                     , bdata::random(0, m) ^ bdata::random(0, m)
                     ^ bdata::xrange(10)
                     , s1, s2, idx)
{
    (void) idx;
    BOOST_TEST(paillier<global>::get_keys().has_public_key());
    BOOST_TEST(paillier<global>::get_keys().has_private_key());
    paillier<global>  p1 = paillier<global>::gen_random(40);
    paillier<global>  p2 = paillier<global>::gen_random(40);
    paillier<global>  p3 = paillier<global>::gen_random(40);


    //Associativity
    BOOST_TEST( equal( (p1 + p2) + p3, p1 + (p2 + p3) ) );

    BOOST_TEST( equal( (s1 * p1) * s2, s1 * (p1 * s2) ) );
    BOOST_TEST( equal(p1 * (s1 * s2), (p1 * s1) * s2) );

    //Commutativity
    BOOST_TEST( equal(p1 + p2, p2 + p1) );
    BOOST_TEST( equal(s1 * p1, p1 * s1) );

    //Neutral Element
    BOOST_TEST( equal(p1 + paillier<global>::encrypt_zero(), p1) );
    BOOST_TEST( equal(p1 * 1, p1) );

    //Distributivity
    BOOST_TEST( equal(s1 * (p1 + p2), s1 * p1 + s1 * p2) );

    //Plus is consistent
    BOOST_TEST( equal(p1 + p1, 2*p1) );
    BOOST_TEST( equal(p1 + p1 + p1, 3*p1) );

    //Consistency with +=, *=
    paillier<global>  enc_s1 = paillier<global>::encrypt(s1);
    paillier<global>  enc_s2 = paillier<global>::encrypt(s2);

    enc_s1 += enc_s2;
    BOOST_TEST( equal(enc_s1, paillier<global>::encrypt(s1+s2) ) );
    enc_s2 += enc_s1;
    BOOST_TEST( equal(enc_s2, paillier<global>::encrypt(s2+s1+s2) ) );
    enc_s1 += enc_s2;
    BOOST_TEST( equal(enc_s1, paillier<global>::encrypt(s1+s2+s2+s1+s2) ) );
    enc_s1 *= s2;
    BOOST_TEST( equal(enc_s1, paillier<global>::encrypt( (s1+s2+s2+s1+s2)*s2) ) );
    enc_s1 *= s1;
    BOOST_TEST( equal(enc_s1, paillier<global>::encrypt( (s1+s2+s2+s1+s2)*s2*s1) ) );
}

#include <iostream>
#include <fstream>

BOOST_DATA_TEST_CASE(paillier_serialization_test
                     , bdata::random(0, m) ^ bdata::xrange(10), s1, idx)
{
    (void) idx;
    paillier<global> p1{paillier<global>::gen_random(40)};
    paillier<global> p2{paillier<global>::gen_random(40)};

    std::ofstream ofs;
    ofs.exceptions(std::ios::failbit | std::ios::badbit);
    ofs.open("test.bin", std::ios::binary | std::ios::trunc | std::ios::out);

    ofs << plaintext{} << paillier<global>::get_keys().get_public() << " "
        << plaintext{} << p1 << p2 << s1 << " "
        << p1 + p2 << p2 * s1
        << decrypt(p1) << plaintext{};
    ofs.close();

    std::ifstream ifs;
    ifs.exceptions(std::ios::failbit | std::ios::badbit);
    ifs.open("test.bin", std::ios::binary | std::ios::in);
    paillier_keys::public_key pub;
    paillier<global>  sp1, sp2, sp1p2, sp2s1;
    plaintext dec_p1, inv_p1, inv_p2, inv_p3;
    int ss1;
    ifs >> inv_p1 >> pub;
    char skip = ifs.get();
    BOOST_TEST(' ' == skip);
    ifs >> inv_p2 >> sp1 >> sp2 >> ss1;
    skip = ifs.get();
    BOOST_TEST(' ' == skip);
    ifs >> sp1p2 >> sp2s1 >> dec_p1 >> inv_p3;
    ifs.close();

    BOOST_TEST(equal(p1, sp1));
    BOOST_TEST(equal(p2, sp2));
    BOOST_TEST(s1 == ss1);
    BOOST_TEST(equal(p1 + p2, sp1p2));
    BOOST_TEST(equal(p2 * s1, sp2s1));
    BOOST_TEST(equal(decrypt(p1), dec_p1));
    BOOST_TEST(!inv_p1.valid());
    BOOST_TEST(!inv_p2.valid());
    BOOST_TEST(!inv_p3.valid());
}
