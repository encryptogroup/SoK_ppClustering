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
#include <mpo19/point.hpp>

#define BOOST_TEST_MODULE point_test
#include <boost/test/included/unit_test.hpp>
#include <boost/test/data/test_case.hpp>
#include <boost/test/data/monomorphic.hpp>

namespace bdata = boost::unit_test::data;

BOOST_AUTO_TEST_CASE(apply_on_elements_test)
{
    using namespace boost;
    point<uint32_t> uut1{4, 5, 6}, uut2{2, 3, 4}, uut3{9, 8, 7};
    point<double> res1(uut1.dimension());
    apply_on_elements(res1
        , [](auto const& a, auto const& b)
          {
            return (double) a * b;
          }
        , uut1, uut2);
    BOOST_TEST(range::equal(res1.elements()
               , std::vector<double>{4*2, 5*3, 6*4}));

    point<std::string> res2 = apply_on_elements(
        [](auto const& a, auto const& b, auto const& c)
        {
            std::string res{};
            res += std::to_string(a);
            res += std::string{" * "};
            res += std::to_string(b);
            res += std::string{" + "};
            res += std::to_string(c);
            res += std::string{" = "};
            res += std::to_string(a*b+c);
            return res;
        }, uut1, uut2, uut3);

    BOOST_TEST(range::equal(res2.elements(), std::vector<std::string>
    {
        std::string{"4 * 2 + 9 = "} + std::to_string(4*2+9),
        std::string{"5 * 3 + 8 = "} + std::to_string(5*3+8),
        std::string{"6 * 4 + 7 = "} + std::to_string(6*4+7)
    }));
}

BOOST_AUTO_TEST_CASE(point_test_examples)
{
    using namespace boost;
    point<uint32_t> uut1(3);
    BOOST_TEST(uut1.dimension() == 3);

    std::vector<uint32_t> v(3);
    point<uint32_t> uut2(v.begin(), v.end());
    BOOST_TEST(uut2.dimension() == 3);

    BOOST_TEST( (point<double>{4.0, 5.0, 6.0}) == (point<uint32_t>{4, 5, 6}) );

    BOOST_TEST( (point<double>{1.2, 1.3, 1.4} + point<int>{3, 5, 4}) 
                    == 
                (point<double>{3+1.2, 5+1.3, 4+1.4}) );
    BOOST_TEST( (point<double>{1.2, 1.3, 1.4} + point<int>{3, 5, 4}) 
                    != 
                 (point<int>{4, 6, 5}) );

    BOOST_TEST( (point<double>{1.2, -1.3, -1.4} - point<int>{-2, 1, -14}) 
                    == 
                (point<double>{1.2-(-2), -1.3-1, -1.4-(-14)}) );

    BOOST_TEST( (point<double>{1.2, 1.3, 1.4} - point<int>{-2, 1, -14}) 
                    != 
                (point<int>{3, 0, 15}) );

    BOOST_TEST( (point<double>{1.2, -1.3, -1.4} * 3) 
                    == 
                (point<double>{1.2*3, -1.3*3, -1.4*3}) );

    BOOST_TEST( (point<double>{1.2, 1.3, 1.4} * 3) != (point<int>{3, 3, 3}) );

    BOOST_TEST( (-point<int>{1, -2, 3}) == (point<int>{-1, 2, -3}) );

    BOOST_TEST( (point<int>{3, 4, 5}.fold()) == 12 );

    BOOST_TEST( (point<int>{2, 5, 7}.fold(std::multiplies<>{})) == 70);

    struct implicit_int_string : public std::string{
        using base = std::string;
        implicit_int_string(int i) : base{std::to_string(i)}{}
        implicit_int_string(base const& b) : base{b}{}
        implicit_int_string(base&& b) : base{std::forward<base>(b)}{}
    };

    BOOST_TEST( (point<int>{2, 5, 7}.fold([](implicit_int_string&& res, int arg){
                     using namespace std::string_literals;
                     return implicit_int_string("f("s + std::move(res) + ", "s + std::to_string(arg) + ")"s);
                 })) 
                 == 
                 "f(f(2, 5), 7)");

    point<uint32_t> uut3{4, 5, 6};

    BOOST_TEST( (uut3.fold()) == 15 );

    point<uint32_t> const uut4{10, 9, 8, 7, 6, 5};

    BOOST_TEST( (uut4.fold(std::multiplies<>{})) == (10*9*8*7*6*5) );
}

constexpr int m = 1000;
auto random_3d_point = 
    bdata::random(0, m) ^ bdata::random(0, m) ^ bdata::random(0, m);

BOOST_DATA_TEST_CASE(operator_consistency_test
                     , random_3d_point ^ random_3d_point ^ random_3d_point 
                       ^ bdata::random(0, m) ^ bdata::random(0, m) 
                       ^ bdata::xrange(100)
                     , x1, y1, z1
                     , x2, y2, z2
                     , x3, y3, z3
                     , s1, s2, idx)
{
    (void) idx;
    point<int> p1{x1, y1, z1}, p2{x2, y2, z2}, p3{x3, y3, z3};

    //Associativity
    BOOST_TEST( ((p1 + p2) + p3) == (p1 + (p2 + p3)) );

    BOOST_TEST( ((s1 * p1) * s2) == (s1 * (p1 * s2)) );
    BOOST_TEST( (p1 * (s1 * s2)) == ((p1 * s1) * s2) );

    //Commutativity
    BOOST_TEST( (p1 + p2) == (p2 + p1) );
    BOOST_TEST( (s1 * p1) == (p1 * s1) );

    //Neutral Element
    BOOST_TEST( (p1 + point<int>{0, 0, 0}) == p1 );
    BOOST_TEST( (p1 * 1) == p1 );

    //Distributivity
    BOOST_TEST( (s1 * (p1 + p2)) == (s1 * p1 + s1 * p2) );
    
    //Plus is consistent
    BOOST_TEST( (p1 + p1) == (2*p1) );
    BOOST_TEST( (p1 + p1 + p1) == (3*p1) );

    //Subtraction and Negation are consistent
    BOOST_TEST( (-p1) == (p1*(-1)) );
    BOOST_TEST( (p1 + (-p2)) == (p1 - p2) );
    BOOST_TEST( (p1 - p2) == (-(p2 - p1)) ); 
    BOOST_TEST( (-(p1 - p2)) == ((-p1) + p2) );
    BOOST_TEST( (p1 - p1) == (point<int>{0, 0, 0}) );
    BOOST_TEST( (point<int>{0, 0, 0} - p1) == (-p1) );
    BOOST_TEST( ((-p1) - p1) == ((-2) * p1) );
    BOOST_TEST( ((-p1) - p1 - p1) == ((-3) * p1) );

    //Consistency with =, +=, -=, *=
    point<int> a_p1 = p1, a_p2 = p2, a_p3{p3};
     BOOST_TEST( (a_p1) == (p1) );
    BOOST_TEST( (a_p2) == (p2) );
    BOOST_TEST( (a_p3) == (p3) );
    a_p1 += p2;
    BOOST_TEST( a_p1 == (p1 + p2) );
    BOOST_TEST( (a_p1 += a_p1 += p3) == (2*(p1 + p2 + p3)) );
    a_p2 *= s1;
    BOOST_TEST( a_p2 == (p2 * s1) );
    BOOST_TEST( (a_p2 *= s2) == (p2 * s1 * s2) );
    a_p3 -= p1;
    BOOST_TEST( a_p3 == (p3 - p1) );
    BOOST_TEST( (a_p3 -= p2) == (p3 - p1 - p2) );
    BOOST_TEST( (a_p3 -= a_p3 -= p3) == (point<int>{0, 0, 0}) );
    
}

