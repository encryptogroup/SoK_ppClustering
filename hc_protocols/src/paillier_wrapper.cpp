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
#include <iostream>
#include <algorithm>

namespace mpo19{

void init_plaintext(paillier_plaintext_t*& pt)
{
    pt = static_cast<paillier_plaintext_t*>(malloc(sizeof(paillier_plaintext_t)));
    assert(pt != nullptr);
    mpz_init(pt->m);
}

void import_into_plaintext(paillier_plaintext_t* pt, void const* m, size_t len)
{
    assert(pt != nullptr);
    assert(m != nullptr);
    //Least significant word first (ABY returns little endian values)
    constexpr int order = -1;
    //ABY returns a pointer to uint8_t values
    constexpr size_t word_size = sizeof(uint8_t);
    //Use local endian for the words
    constexpr int endian = 0;
    constexpr int nails = 0;
    mpz_import(pt->m, len, order, word_size, endian, nails, m);
}

plaintext::plaintext(void const* m, size_t len)
{
    init_plaintext(plain_);
    assert(plain_ != nullptr);
    import_into_plaintext(plain_, m, len);
}

plaintext operator-(plaintext const& plain)
{
    assert(plain.valid());
    paillier_plaintext_t* res;
    //res must be allocated with malloc, as it is deleted using free in destructor
    res = static_cast<paillier_plaintext_t*>(malloc(sizeof(paillier_plaintext_t)));
    assert(nullptr != res);
    mpz_init(res->m);
    mpz_neg(res->m, plain.plain_->m);
    return plaintext{res};
}

std::ostream& operator<<(std::ostream& os, plaintext const& plain)
{
    assert(os.good() && "Stream should not be in error state");
    //Using a size of 0, means invalid plaintext
    if(!plain.valid()) {
        assert(!plain.valid());
        os.put(0);
        return os;
    }
    assert(plain.valid());
    //Encode length of plaintext
    char buf[encode_len_buffer_size];
    size_t len = plain.size();
    assert(1 <= len && "Length of plaintext must at least be 1");
    size_t s = encode_len(buf, len);

    assert(1 <= s && "Length encoding must at least take one byte");
    assert(s == words_to_hold_bitlen(ceil_log2(len), CHAR_BIT - 1)
           && "Amount of bytes length encoding must be taking");
    assert(std::none_of(buf, buf + s, [](auto const& e) {
        return e == 0;
    }) && "There should be no 0 value in range [buf, buf + s)");

    os.write(buf, s);
    assert(os.good() && "Stream should not go in error state here");

    //Write plaintext
    std::vector<char> exp = plain;
    assert(!exp.empty() && "plain is valid, so exp should not be empty");
    os.write(exp.data(), len);
    return os;
}

std::istream& operator>>(std::istream& is, plaintext& plain)
{
    assert(is.good() && "Stream should not be in error state");
    size_t len = decode_len([&](unsigned char* buf, size_t l) {
        is.read(reinterpret_cast<char*>(buf), l);
        assert(is.good() && "Stream should not go in error state here");
        assert(size_t(is.gcount()) == l && "exactly l characters must have been read");
    });
    //len is 0, so received plaintext is invalid
    if(len == 0) {
        plain = plaintext{};
        assert(len == 0);
        assert(!plain.valid() && "len was 0, so plaintext should not be valid");
        return is;
    }
    assert(1 <= len && "len should be 1 here");

    auto mem = std::make_unique<char[]>(len);
    assert(nullptr != mem && "mem must not be nullptr");
    is.read(mem.get(), len);
    assert(is.good() && "Stream should not go in error state here");
    assert(size_t(is.gcount()) == len && "exactly len characters must have been read");

    plain = plaintext{static_cast<void const*>(mem.get()), len};
    assert(plain.valid() && "plain must be valid now");
    return is;
}

} //namespace mpo19
