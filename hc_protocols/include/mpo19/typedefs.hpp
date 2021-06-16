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

#ifndef MPO19_TYPEDEFS_30062020
#define MPO19_TYPEDEFS_30062020

#include <mpo19/point.hpp>

#include <memory>
#include <scoped_allocator>
#include <type_traits>

#include <boost/multi_array.hpp>

template<typename T, typename Alloc = std::allocator<T>>
struct allocator_adaptor {
    using type = Alloc;
};

template<typename T, typename Alloc>
struct allocator_adaptor<point<T>, Alloc> {
    using type = std::scoped_allocator_adaptor<Alloc>;
};

template<typename T>
using table_t = boost::multi_array<T, 2, typename allocator_adaptor<T>::type>;
template<typename T>
using array_t = boost::multi_array<T, 1, typename allocator_adaptor<T>::type>;
template<typename T>
using row_t = typename boost::multi_array<T, 2>::template array_view<1>::type;
using index_range_t = boost::multi_array_types::index_range;
using index_t = boost::multi_array_types::index;
using extent_range_t = boost::multi_array_types::extent_range;

//Alias for ABY roles
//Whether P1 assumes the role of client or server doesn't matter.
//It may even be inconsistent with the setup phase.
//However, it is important that P2 is the party holding the blinded
//values (matrix B), while P1 is the one holding the blinds (matrix R).
#define ROLE_P1 SERVER
#define ROLE_P2 CLIENT

#endif //MPO19_TYPEDEFS_30062020
