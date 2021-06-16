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

#include <mpo19/random.hpp>
#include <iostream>
#include <fstream>

using result_type = buffered_secure_rng::result_type;
using result_pointer = buffered_secure_rng::result_pointer;
using result_reference = buffered_secure_rng::result_reference;
using buffer_t = buffered_secure_rng::buffer_t;

constexpr size_t buffered_secure_rng::default_buffer_size;

buffered_secure_rng::buffered_secure_rng()
    : buffer_{}, pos_{buffer_.crbegin()} {}

buffered_secure_rng::buffered_secure_rng(size_t buffer_size)
    : buffer_{}, pos_{buffer_.crbegin()}
{
    fill_buffer(buffer_size);
}

result_type buffered_secure_rng::operator()()
{
    if(pos_ == buffer_.crend()) {
        fill_buffer(default_buffer_size);
    }
    result_type res = *pos_;
    ++pos_;
    return res;
}

result_pointer buffered_secure_rng::get_rand_numbers(size_t amount)
{
    assert(0 < amount);
    fill_buffer(amount);
    size_t remaining = this->remaining();
    assert(remaining >= amount);
    result_pointer res = buffer_.data() + (remaining - amount);
    pos_ += amount;
    assert(this->remaining() == remaining - amount);
    return res;
}

std::pair<result_pointer, size_t> buffered_secure_rng::get_rand_bits(size_t amount)
{
    using uresult_type = std::make_unsigned_t<result_type>;
    constexpr size_t word_bitlen = sizeof(result_type) * CHAR_BIT;
    //Get minimal amount of result_types to hold amount bits
    size_t result_amount = words_to_hold_bitlen(amount, word_bitlen);
    //Get result_amount number of rand_numbers
    result_pointer nums = get_rand_numbers(result_amount);
    //Reference to last element in returned rand_numbers
    result_reference last = nums[result_amount-1];
    //Number of spare bits in last element (these will be set to 0)
    size_t spare_bits = result_amount * word_bitlen  - amount;

    //mask must be an unsigned type, so that the right shift operator is not implementation defined
    uresult_type mask = 0;
    mask = ~mask;

    //bitwise operators should be applied only on unsigned types.
    last = static_cast<result_type>(static_cast<uresult_type>(last) & (mask >> spare_bits));
    return std::make_pair(nums, result_amount);
}

void buffered_secure_rng::fill_buffer(size_t amount)
{
    size_t remaining = this->remaining();
    if(remaining < amount) {
        amount = std::max(amount, default_buffer_size);
        //Do not free previously allocated elements, but allocate if needed
        if(buffer_.size() < amount) buffer_.resize(amount);

        assert(buffer_.size() >= amount);

        std::ifstream rand;
        rand.exceptions(std::ios::failbit | std::ios::badbit | std::ios::eofbit);
        rand.open("/dev/urandom", std::ios::binary);
        rand.read( reinterpret_cast<char*>(buffer_.data() + remaining)
                                      ,  (amount - remaining) * sizeof(result_type));

        //Advance pos_ to point to the first element of newly created random data
        pos_ = buffer_.crbegin() + (buffer_.size() - amount);
        assert(this->remaining() == amount);
        
        assert(std::addressof(*pos_) == buffer_.data() + amount - 1);
    }
}

random_permutation::random_permutation(size_t n)
{
    using rng_t  = buffered_secure_rng;
    rng_t rng;

    for(size_t i = 0; i != n; ++i) {
        permutation_.push_back(i);
    }
    permutation_inplace_.reserve(n-1);

    //We store generated random values in permutation_inplace_
    //When permuting inplace we then perform the Fisher-Yates shuffle,
    //using the stored random values, to yield the same permutation
    //for every invokation.
    for(size_t i = 0; i != n-1; ++i) {
        permutation_inplace_.push_back(std::uniform_int_distribution<size_t> {i, n-1}(rng));
    }
    //We shuffle the vector inplace to match the semantic of permute_inplace
    permute_inplace(permutation_);
}
