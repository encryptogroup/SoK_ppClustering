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

#include <mpo19/cluster.hpp>
#include <iostream>

#define BOOST_TEST_MODULE linkage_test
#include <boost/test/included/unit_test.hpp>
#include <boost/test/data/test_case.hpp>
#include <boost/test/data/monomorphic.hpp>

namespace bdata = boost::unit_test::data;

BOOST_DATA_TEST_CASE(linkage_test
                     , bdata::xrange(1)
                     , idx)
{
    (void) idx;
    constexpr size_t len = 10;
    std::vector<cluster> dendrogram;
    dendrogram.reserve(len + idx);
    for(size_t i = 0; i < len; ++i){
        dendrogram.emplace_back(i);
    }
    
    merge(dendrogram[0], dendrogram[1], 1);
    merge(dendrogram[2], dendrogram[3], 2);
    merge(dendrogram[0], dendrogram[2], 3);
    merge(dendrogram[4], dendrogram[8], 4);
    merge(dendrogram[0], dendrogram[8], 5);
    merge(dendrogram[5], dendrogram[6], 6);
    merge(dendrogram[5], dendrogram[7], 7);
    
    
    for(auto&& d : dendrogram){
        std::cout << d << std::endl;
    }

}
