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

#ifndef MPO19_UTILITY_20062020
#define MPO19_UTILITY_20062020

#include <gmp.h>

extern "C" {
#include <paillier.h>
}

#include <boost/hana.hpp>

#include <boost/preprocessor/repetition.hpp>
#include <boost/preprocessor/arithmetic/sub.hpp>
#include <boost/preprocessor/punctuation/comma_if.hpp>

#include <cinttypes>
#include <vector>
#include <limits>
#include <cassert>
#include <memory>
#include  <experimental/propagate_const>

#define FOR_EACH_VARIADIC(expr)  using expander = int[]; expander{ ( (void) (expr), 0 )... }

template<typename T, typename U>
inline void nothrow_swap(T& t, U& u) noexcept
{
    using std::swap;
    static_assert(noexcept(swap(t, u))
                  , "swap may throw an exception");
    swap(t, u);
}

constexpr size_t words_to_hold_bitlen(size_t bits, size_t bits_in_word)
{
    assert(0 != bits_in_word);
    return bits/bits_in_word + (bits % bits_in_word == 0 ? 0 : 1);
}

constexpr size_t bytes_to_hold_bitlen(size_t bits)
{
    return words_to_hold_bitlen(bits, CHAR_BIT);
}

namespace detail
{

constexpr uint64_t integer_log2(uint64_t x, bool ceil)
{
    if(x == 0) return 0;

    constexpr size_t num_bits = sizeof(decltype(x)) * CHAR_BIT;
    uint64_t a = x;
    uint64_t res = 0;
    for(size_t i = num_bits/2; i != 0; i /=2) {
        //If right-shifting a results in a number not equal zero, then a
        //has between i and prev(i) bits, where prev(i) is defined as:
        //prev(0) -> num_bits
        //prev(i) -> value of i in previous iteration
        uint64_t tmp = a >> i;
        if(tmp != 0) {
            //We set a to the right shifted value
            a = tmp;
            //We add i to res, as a has at least res + i bits
            res = res + i;
        }
        //else a has between res and i bits, so we do nothing here
    }
    assert(a == 1);
    //res now contains ⌊log2(x)⌋
    assert( (1u << res) <= x && x < (1u << (res + 1)) );
    if(ceil && (1u << res) < x)
        return res + 1;
    else
        return res;
}

} //detail

//Calculates ceil_log2(a) = ⌈log2(a)⌉, for a ≠ 0
//and ceil_log2(a) = 1, for a = 0
constexpr uint64_t ceil_log2(uint64_t a)
{
    return detail::integer_log2(a, true);
}

//Calculates ceil_log2(a) = ⌊log2(a)⌋, for a ≠ 0
//and floor_log2(a) = 1, for a = 0
constexpr uint64_t floor_log2(uint64_t a)
{
    return detail::integer_log2(a, false);
}

constexpr unsigned char clear_msb = 0x7F;
constexpr unsigned char select_msb = 0x80;

template<typename CharT>
size_t encode_len(CharT* buf, size_t len)
{
    size_t i = 0;
    for(; len != 0; ++i) {
        unsigned char c =  static_cast<unsigned char>(len & clear_msb);
        len >>= CHAR_BIT - 1;
        if(len != 0) c |= select_msb;
        buf[i] = static_cast<CharT>(c);
    }
    assert(len == 0);
    return i;
}

template<typename Reader>
size_t decode_len(Reader&& reader)
{
    size_t res = 0;
    //End condition evaluated in body
    for(int i = 0;; i += (CHAR_BIT - 1) ) {
        assert(res >> i == 0);
        unsigned char c = 0;
        reader(&c, sizeof(c));
        res |= (c & clear_msb) << i;
        if((c & select_msb) == 0) break;
    }
    return res;
}

namespace detail
{
constexpr size_t calculate_encode_len_buffer_size()
{
    size_t res = (sizeof(size_t) * CHAR_BIT) / (CHAR_BIT - 1);
    if(res * (CHAR_BIT - 1) < sizeof(size_t) * CHAR_BIT) res += 1;
    return res;
}
}

constexpr size_t encode_len_buffer_size =
    detail::calculate_encode_len_buffer_size();


//Check if T is a range. A range is defined as any instance t of type T,
//for which both, boost::begin(t) and boost::end(t) form valid expressions
template<typename T>
constexpr bool is_range()
{
    using boost::hana::is_valid;
    //Function to check if std::begin(r) forms a valid expression
    auto test_begin = [](auto const& r) -> decltype((void) std::begin(r)) {};
    //Function to check if std::end(r) forms a valid expression
    auto test_end = [](auto const& r) -> decltype((void) std::end(r)) {};
    //We use decltype so that std::declval<T>() is in unevaluated context
    //The return type is either hana::true_ or hana::false_, which is implicitly
    //convertible to bool, holding the expected boolean value
    using bool_t = decltype(is_valid(test_begin, std::declval<T>())
                            && is_valid(test_end, std::declval<T>()));
    return bool_t{};
}

template<typename T1, typename T2, typename... Ts>
constexpr bool is_range()
{
    return boost::hana::and_(is_range<T1>(), is_range<T2>(), is_range<Ts>()...);
}

//Check if T serializable. An instance t of type T is serializable if, the expression
//os << t is valid, with os being an instance derived from std::ostream
template<typename T>
constexpr bool is_serializable()
{
    using boost::hana::is_valid;
    auto test_serialize = [](auto&& r) -> decltype((void) (std::declval<std::ostream>() << r)) {};
    using bool_t = decltype(is_valid(test_serialize, std::declval<T>()));
    return  bool_t{};
}

//Check if T deserializable. An instance t of type T is deserializable if, the expression
//is >> t is valid, with is being an instance derived from std::istream
template<typename T>
constexpr bool is_deserializable()
{
    using boost::hana::is_valid;
    auto test_deserialize = [](auto&& r) -> decltype((void) (std::declval<std::istream>() >> r)) {
        static_assert(std::is_same<decltype(r), T&>::value, "");
    };
    using bool_t = decltype(is_valid(test_deserialize, std::declval<T>()));
    return  bool_t{};
}

//Preprocessor will generate the following set of functions:
//template <typename Range0, ..., typename RangeN-1, typename Op>
//void multi_for_each(Range0&& range0, ..., RangeN-1&& rangeN-1, Op const& op);
//
//Applies N-Ary operation to every element in the input ranges, which are iterated in
//a parallel fashion, i.e. in iteration i, the following holds:
//op(range0[i], ..., rangeN-1[i]).
//Furthermore if the input ranges have multiple dimension, the ranges are iterated,
//as if they were flattened before.
//All ranges must be of equal size and equal dimension. In case of different sizes the
//behavior is unspecified and in case of different dimensions compilation fails.


//Using variadic templates to generate above functions, would yield more complex code

#define MPO19_INIT_ITERATORS(z, n, unused) \
    auto it ## n = std::begin(range ## n);

#define MPO19_DEREF_ITERATORS(z, n, unused) \
    auto&&  r ## n = *it ## n;


#define MPO19_MULTI_FOR_EACH(z, n, unused) \
    template <BOOST_PP_ENUM_PARAMS(n, typename Range), typename Op, \
              std::enable_if_t<is_range<BOOST_PP_ENUM_PARAMS(n, Range)>(), int> = 0> \
    void multi_for_each(BOOST_PP_ENUM_BINARY_PARAMS(n, Range, && range), Op&& op) \
    { \
        /*Initialize it0, ..., itN-1 with begin of range0, ..., rangeN-1*/ \
        BOOST_PP_REPEAT(n, MPO19_INIT_ITERATORS, unused) \
        \
        auto end0 = std::end(range0); \
        for(; it0 != end0; BOOST_PP_ENUM_PARAMS(n, ++it)/* == ++it0, ..., ++itN-1*/) \
        { \
            /*Dereference iterators and hold an rvalue reference, to prevent pass by value later on*/ \
            BOOST_PP_REPEAT(n, MPO19_DEREF_ITERATORS, unused) \
            \
            /*Recurse as the dereferenced values might be ranges themselves */ \
            multi_for_each(BOOST_PP_ENUM_PARAMS(n, r), op); \
            \
        } \
    } \
    /*Passed a non-range, so we simply apply op on every argument*/ \
    template <BOOST_PP_ENUM_PARAMS(n, typename T), typename Op, \
              std::enable_if_t<!is_range<BOOST_PP_ENUM_PARAMS(n, T)>(), int> = 0> \
    void multi_for_each(BOOST_PP_ENUM_BINARY_PARAMS(n, T, && t), Op&& op) \
    { \
        op(BOOST_PP_ENUM_PARAMS(n, t)); \
    }

BOOST_PP_REPEAT_FROM_TO(1, 4, MPO19_MULTI_FOR_EACH, ~)

#undef MPO19_INIT_ITERATORS
#undef MPO19_DEREF_ITERATORS
#undef MPO19_MULTI_FOR_EACH

namespace detail
{

template<typename T>
struct always_void {
    using type = void;
};

template<typename T>
using always_void_t = typename always_void<T>::type;

} //detail

template<typename T, typename = void>
struct is_invocable : std::false_type {};

template<typename F, typename... Args>
struct is_invocable<F(Args...)
, detail::always_void_t<decltype(
std::declval<F>()(std::declval<Args>()...))>> : std::true_type {};


constexpr struct identity_t {
    template<typename T>
    inline T operator()(T&& t) const
    {
        return t;
    }
} identity;

namespace detail
{

template<typename F, typename Arg,
         std::enable_if_t<
             is_invocable<F(Arg&&)>::value &&
             !is_invocable<F(Arg&&, size_t)>::value, int> = 0>
auto invoke_maybe_with_idx(F&& f, Arg&& arg, size_t)
{
    return f(std::forward<Arg>(arg));
}

template<typename F, typename Arg,
         std::enable_if_t<is_invocable<F(Arg&&, size_t)>::value, int> = 0>
auto invoke_maybe_with_idx(F&& f, Arg&& arg, size_t s)
{
    return f(std::forward<Arg>(arg), s);
}

template<typename RandomIt, typename BinaryOperation, typename LeafTransformation>
auto tree_accumulate_impl(RandomIt first
                          , RandomIt last
                          , BinaryOperation&& op
                          , LeafTransformation&& lt
                          , RandomIt const& base)
{
    if(std::distance(first, last) == 1)
        return invoke_maybe_with_idx(lt, *first, size_t(std::distance(base, first)));
    else
        return op(
                   tree_accumulate_impl(first, first + std::distance(first, last)/2, op, lt, base),
                   tree_accumulate_impl(first + std::distance(first, last)/2, last, op, lt, base)
               );
}

} //detail

template<typename RandomIt, typename BinaryOperation, typename LeafTransformation = decltype(identity)&>
auto tree_accumulate(RandomIt first, RandomIt last, BinaryOperation&& op, LeafTransformation&& lt = identity)
{
    return detail::tree_accumulate_impl(first, last, op, lt, first);
}

template<typename RandomAccessRange, typename BinaryOperation, typename LeafTransformation = decltype(identity)&>
auto tree_accumulate(RandomAccessRange&& range, BinaryOperation&& op, LeafTransformation&& lt = identity)
{
    return tree_accumulate(range.begin(), range.end(), op, lt);
}

#endif //MPO19_UTILITY_20062020
