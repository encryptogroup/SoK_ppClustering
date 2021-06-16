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

#ifndef MPO19_POINT_12062020
#define MPO19_POINT_12062020

#include <cassert>
#include <vector>
#include <type_traits>
#include <iostream>
#include <memory>
#include <scoped_allocator>

#include <boost/range/algorithm.hpp>

template<typename, typename> struct point;

template<bool Cond, typename T, std::enable_if_t<Cond, int> = 0>
T&& move_if(T& t)
{
    return static_cast<T&&>(t);
}

template<bool Cond, typename T, std::enable_if_t<!Cond, int> = 0>
T& move_if(T& t)
{
    return t;
}

constexpr struct apply_on_elements_t {
    template<typename R, typename RAlloc
             , typename NAryOp
             , typename T, typename Alloc
             ,  typename... Ts, typename... Allocs>
    void operator()(point<R, RAlloc>& res
                    , NAryOp const& op
                    , point<T, Alloc> const& p
                    , point<Ts, Allocs> const&... ps) const
    {
        size_t d = p.dimension();
        assert(res.dimension() == d);
        for(size_t i = 0; i != d; ++i) {
            res[i] = op(p[i], ps[i]...);
        }
    }

    template<typename NAryOp, typename T, typename Alloc, typename... Ts, typename... Allocs>
    auto operator()(NAryOp const& op
                    , point<T, Alloc> const& p
                    , point<Ts, Allocs> const&... ps) const
    {
        using result_element_t = decltype(op(p[0], ps[0]...));
        using result_allocator_t = typename std::allocator_traits<Alloc>::template rebind_alloc<result_element_t>;
        using result_t = point<result_element_t, result_allocator_t>;
        size_t d = p.dimension();
        result_t res;
        res.reserve(d);
        for(size_t i = 0; i != d; ++i) {
            res.push_back(op(p[i], ps[i]...));
        }
        assert(res.dimension() == d);
        return res;
    }

} apply_on_elements;

template<typename T, typename Alloc = std::allocator<T>>
struct point : private std::vector<T, Alloc> {
private:
    using Base_ = std::vector<T, Alloc>;
public:
    using Base_::Base_;
    using Base_::data;
    using element_type = typename Base_::value_type;

    using Base_::operator[];

    point(point const&) = default;
    point(point&&) = default;
    ~point() = default;
    point& operator=(point const&) = default;
    point& operator=(point&&) = default;

    point& operator+=(point const& rhs)
    {
        apply_on_elements(*this, std::plus<> {}, *this, rhs);
        return *this;
    }

    point& operator-=(point const& rhs)
    {
        apply_on_elements(*this, std::minus<> {}, *this, rhs);
        return *this;
    }

    point& operator*=(element_type rhs)
    {
        auto scalar_multiply = [&rhs](auto&& a) {
            return std::forward<decltype(a)>(a) * rhs;
        };
        apply_on_elements(*this, scalar_multiply, *this);
        return *this;
    }

    size_t dimension() const
    {
        return Base_::size();
    }

    Base_ const& elements() const
    {
        return static_cast<Base_ const&>(*this);
    }

    template<typename BinaryOp = std::plus<>>
    auto fold(BinaryOp&& op = std::plus<> {}) const&
    {
        return fold_impl(*this, std::forward<BinaryOp>(op));
    }

    template<typename BinaryOp = std::plus<>>
    auto fold(BinaryOp&& op = std::plus<> {}) && {
        return fold_impl(std::move(*this), std::forward<BinaryOp>(op));
    }

private:

    template<typename Point, typename BinaryOp>
    static auto fold_impl(Point&& p, BinaryOp&& op)
    {
        using res_t = decltype(op(std::declval<T>(), std::declval<T>()));
        constexpr bool p_is_rvalue = std::is_rvalue_reference<Point&&>::value;
        size_t d = p.dimension();
        assert(1 <= d);
        res_t res{move_if<p_is_rvalue>(p[0])};
        for(size_t i = 1; i != d; ++i) {
            res = op(std::move(res), move_if<p_is_rvalue>(p[i]));
        }
        return res;
    }

    friend void swap(point& lhs, point& rhs)
    {
        using std::swap;
        using base_t = point::Base_;
        swap(static_cast<base_t&>(lhs), static_cast<base_t&>(rhs));
    }

    friend std::istream& operator>>(std::istream& is, point& p)
    {
        is >> std::ws;
        char c;
        is >> c;
        if('/' == c) return is;
        assert('(' == c);
        size_t i = 0;
        size_t dim = p.dimension();
        do {
            //point does have a value at that position, so insert
            if(i < dim) {
                is >> p[i];
                //point does not have a value at that position, so push_back
            } else {
                T t;
                is >> t;
                p.push_back(std::move(t));
            }
            ++i;
            is >> c;
        } while(c != ')');

        return is;
    }

    friend struct apply_on_elements_t;
};

template<typename T>
std::ostream& operator<<(std::ostream& os, point<T> const& p)
{
    size_t dim = p.dimension();
    if(0 == dim) return os << '/';
    os << '(';
    for(size_t i = 0; i != dim; ++i) {
        os << p[i];
        if(i != dim-1) os << ',';
    }
    os << ')';
    return os;
}

template<typename A, typename AllocA, typename B, typename AllocB>
inline auto operator+(point<A, AllocA> const& lhs, point<B, AllocB> const& rhs)
{
    return apply_on_elements(std::plus<> {}, lhs, rhs);
}

template<typename A, typename AllocA, typename B, typename AllocB>
inline auto operator-(point<A, AllocA> const& lhs, point<B, AllocB> const& rhs)
{
    return apply_on_elements(std::minus<> {}, lhs, rhs);
}

template<typename A, typename AllocA>
inline auto operator-(point<A, AllocA> const& p)
{
    return apply_on_elements(std::negate<> {}, p);
}

template<typename A, typename AllocA, typename ArithT
         , std::enable_if_t<std::is_arithmetic<ArithT>::value, int> = 0>
inline auto operator*(point<A, AllocA> const& lhs, ArithT const& rhs)
{
    auto scalar_multiply = [&rhs](auto&& a) {
        return std::forward<decltype(a)>(a) * rhs;
    };
    return apply_on_elements(scalar_multiply, lhs);
}

template<typename A, typename AllocA, typename ArithT
         , std::enable_if_t<std::is_arithmetic<ArithT>::value, int> = 0>
inline auto operator*(ArithT const& lhs, point<A, AllocA> const& rhs)
{
    return rhs * lhs;
}

template<typename A, typename AllocA, typename B, typename AllocB>
inline auto operator==(point<A, AllocA> const& lhs, point<B, AllocB> const& rhs)
{
    using boost::range::equal;
    return equal(lhs.elements(), rhs.elements());
}

template<typename A, typename AllocA, typename B, typename AllocB>
inline auto operator!=(point<A, AllocA> const& lhs, point<B, AllocB> const& rhs)
{
    return !(lhs == rhs);
}

#endif //MPO19_POINT_12062020
