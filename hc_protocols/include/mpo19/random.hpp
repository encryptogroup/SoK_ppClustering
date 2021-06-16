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

#ifndef MPO19_RANDOM_22062020
#define MPO19_RANDOM_22062020

#include <mpo19/utility.hpp>

#include <boost/iterator/permutation_iterator.hpp>
#include <boost/range/adaptors.hpp>
#include <boost/range/iterator_range_core.hpp>

#include <cassert>
#include <algorithm>
#include <random>
#include <memory>
#include <utility>
#include <iterator>
#include <functional>
#include <algorithm>
#include <type_traits>

struct buffered_secure_rng {
    using result_type = uint8_t;
    using result_pointer = std::vector<result_type>::pointer;
    using result_reference = uint8_t&;
    using buffer_t = std::vector<result_type>;

    static constexpr size_t default_buffer_size = 1024;

    buffered_secure_rng();

    explicit buffered_secure_rng(size_t buffer_size);

    static constexpr result_type min()
    {
        return std::numeric_limits<result_type>::min();
    }

    static constexpr result_type max()
    {
        return std::numeric_limits<result_type>::max();
    }

    //Get one random result_type
    //After calling this function every pointer to random data
    //might be invalidated afterwards.
    result_type operator()();

    //Getting a pointer to random data containing at least amount data bytes
    //The resulting pointer is not owned and may be invalidated after
    //succeeding calls to non-const functions. It is the caller's responsibility
    //to secure the returned data before it gets invalidated.
    result_pointer get_rand_numbers(size_t amount);

    //Getting a pointer to random data holding exactly amount data bits
    //If bits is not a multiple of CHAR_BIT, then the remaining bits of the last
    //byte will be set to a random value, while the other bits will be set to 0
    std::pair<result_pointer, size_t> get_rand_bits(size_t amount);

    size_t remaining() const
    {
        return std::distance(pos_, buffer_.crend());
    }

    //Buffer is filled to contain at least amount random numbers
    //May invalidate every pointer to random data
    void fill_buffer(size_t amount);

    //We move the buffer to another location andre initialize our buffer.
    //Allows a user to safely securre random numbers without needing
    //to reallocate memory
    buffer_t move_buffer()
    {
        buffer_t res = std::move(buffer_);
        buffer_ = buffer_t{};
        pos_ = buffer_.crbegin();
        return res;
    }

private:
    buffer_t buffer_;
    buffer_t::const_reverse_iterator pos_;
};

struct random_permutation {
private:
    using permutation_container_t_ = std::vector<size_t>;

public:

    random_permutation(size_t n);

    //Permute a range lazily
    //Note that accessing the elements of the range through operator[] yields
    //undefined behaviour. Access the elements through operator(). For example
    //instead of writing:
    //auto e = permuted[0];
    //write this:
    //auto e = permute(0);
    template<typename Range, std::enable_if_t<is_range<Range>(), int> = 0>
    auto permute(Range const& range) const
    {
        return permute_(range);
    }

    //Permute a range greedily and in-place
    //The input range is modified in the process of doing so.
    template<typename... Args>
    void permute_inplace(Args&... args) const
    {
        FOR_EACH_VARIADIC(permute_inplace_(args));
    }

private:

    template<typename T>
    static constexpr bool is_range_of_range_()
    {
        return is_range<T>()
               && is_range<typename std::iterator_traits<
               decltype(std::begin(std::declval<T>()))>::reference>();
    }

    template<typename Range, std::enable_if_t<is_range_of_range_<Range>(), int> = 0>
    auto permute_(Range const& range) const
    {
        using boost::adaptors::transformed;
        //lazily transform this range, so that every subrange is permuted
        auto rec = range | transformed([this](auto const& range) {
            return this->permute_(range);
        });
        //lazily permute this range
        return permute_impl_(rec);
    }

    //Base case: Range has no sub-ranges, so only lazily permute this range
    template<typename Range, std::enable_if_t<!is_range_of_range_<Range>(), int> = 0>
    auto permute_(Range const& range) const
    {
        return permute_impl_(range);
    }

    //Perform the lazy permutation
    template<typename Range>
    auto permute_impl_(Range const& range) const
    {
        auto first = boost::make_permutation_iterator(std::begin(range), permutation_.begin());
        auto last = boost::make_permutation_iterator(std::begin(range), permutation_.end());
        return boost::make_iterator_range(first, last);
    }

    //Permute range greedily and inplace
    template<typename Range, std::enable_if_t<is_range<Range>(), int> = 0>
    void permute_inplace_(Range& range) const
    {
        using range_iterator_t = typename Range::iterator;
        using ref_t = typename range_iterator_t::reference;

        range_iterator_t first = std::begin(range);
        range_iterator_t last = std::end(range);

        assert((size_t) std::distance(first, last) == permutation_inplace_.size() + 1);

        size_t pos = 0;
        //last is decremented by 2, as we iterate until n-2
        --last;
        for(range_iterator_t it = first; it != last; ++it, ++pos) {
            //pos should never be greater than the size of permutation
            assert(pos < permutation_inplace_.size());
            size_t target = permutation_inplace_[pos];
            //The target index to permute with should be between the current position and n-1
            assert(pos <= target );
            assert(target <= (size_t) std::distance(first, last));
            //Permute-swap current position with target position
            ref_t lhs_ref = *it;
            if(pos != target) {
                ref_t rhs_ref =  *(first + target);
                multi_for_each(lhs_ref, rhs_ref, [](auto& lhs,auto& rhs) {
                    using std::swap;
                    swap(lhs, rhs);
                });
            }
            permute_inplace_(lhs_ref);
        }
        ref_t ref = *last;
        permute_inplace_(ref);
    }

    //Permuting a value is a noop
    template<typename T, std::enable_if_t<!is_range<T>(), int> = 0>
    void permute_inplace_(T&) const {}

    //permutation_ contains the indices used for a lazy permutation
    //permutation_inplace_ contains the indices to perform a permutation inplace
    permutation_container_t_ permutation_, permutation_inplace_;
};

template<typename Range>
auto random_sample(Range& input, size_t sample_size)
{
    using std::swap;
    using rng_t  = buffered_secure_rng;
    std::vector<std::decay_t<decltype(input[0])>> res;
    res.reserve(sample_size);
    
    rng_t r;
    for(size_t i = 0; i != sample_size; ++i){
        size_t selection = std::uniform_int_distribution<std::size_t>{0, input.size() - (i + 1)}(r);
        assert(0 <= selection && selection < input.size() - i);
        res.emplace_back(input[selection]);
        swap(input[selection], input[input.size() - (i + 1)]);
    }
    
    return res;
}

#endif //MPO19_RANDOM_22062020
