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

#ifndef MPO19_PAILLIER_WRAPPER_12062020
#define MPO19_PAILLIER_WRAPPER_12062020


//utility.hpp also includes paillier.h, which is not
//guarded, therefore we must not include it directly
//here or else we might get multiple inclusion of paillier.h
#include <mpo19/config.hpp>
#include <mpo19/utility.hpp>
#include <mpo19/random.hpp>

#include <cstring>
#include <cassert>
#include <memory>
#include <type_traits>
#include <iostream>

namespace mpo19
{

typedef const char(key_id)[];

//==========================================================================
//struct plaintext
//
//Struct describing a plaintext value that may be encrypted to a paillier
//instance or returned by decrypting a paillier instance.
//==========================================================================
struct plaintext {

    plaintext() = default;


    plaintext(plaintext const& rhs)
        : plain_{static_cast<paillier_plaintext_t*>(malloc(sizeof(paillier_plaintext_t)))}
    {
        assert(plain_ != nullptr);
        mpz_init_set(plain_->m, rhs.plain_->m);
    }

    plaintext(plaintext&& rhs) noexcept
    {
        swap(*this, rhs);
    }

    plaintext& operator=(plaintext rhs) noexcept
    {
        swap(*this, rhs);
        return *this;
    }

    plaintext& operator+=(plaintext const& rhs)
    {
        assert(valid() && rhs.valid());
        mpz_add(plain_->m, plain_->m, rhs.plain_->m);
        return *this;
    }

    inline friend plaintext operator+(plaintext&& lhs, plaintext const& rhs)
    {
        assert(lhs.valid() && rhs.valid());
        lhs += rhs;
        return std::move(lhs);
    }

    inline friend plaintext operator+(plaintext const& lhs, plaintext&& rhs)
    {
        return std::move(rhs) + lhs;
    }

    inline friend plaintext operator+(plaintext&& lhs, plaintext&& rhs)
    {
        return std::move(lhs) + rhs;
    }

    friend plaintext operator-(plaintext const& plain);

    ~plaintext()
    {
        if(nullptr != plain_)
            paillier_freeplaintext(plain_);
    }

    plaintext(paillier_plaintext_t* plain) : plain_{plain} {}

    plaintext(unsigned long ui) : plain_{paillier_plaintext_from_ui(ui)} {}

    plaintext(void const* m, size_t len);

    plaintext(char const* str)
        : plaintext{static_cast<void const*>(str), strlen(str)} {}

    bool valid() const
    {
        return nullptr != plain_;
    }

    size_t size() const
    {
        return bytes_to_hold_bitlen(mpz_sizeinbase(plain_->m, 2));
    }

    operator paillier_plaintext_t*() const
    {
        return plain_;
    }

    //Implicit conversion to std::vector<T> with T being an arithmetic type
    //In case of std::vector<cv char> an additional char is allocated and
    //set to '\0', so that it can be used as cstring. Note that only
    //std::vector<cv char> will be null-terminated, while
    //std::vector<cv unsigned char> will not.
    //If plaintext is not in a valid state, a default constructed std::vector<T> will be returned.
    template<typename T
             , std::enable_if_t<
                 std::is_arithmetic<T>::value
                 && !std::is_floating_point<T>::value, int> = 0>
    operator std::vector<T>() const
    {
        constexpr bool null_termination = std::is_same<T, std::remove_cv_t<char>>::value;
        if(!valid()) return std::vector<T> {};

        size_t req_size = words_to_hold_bitlen(mpz_sizeinbase(plain_->m, 2), sizeof(T) * CHAR_BIT);
        assert(1 <= req_size && "At least one word should be reserved for the result");
        //Allocate an additional item for null termination, if T = char
        if(null_termination) ++req_size;

        std::vector<T> res;
        res.resize(req_size);
        assert(res.size() == req_size);

        size_t num_read;
        //Least significant word first (ABY expects little endian ordering)
        constexpr int order = -1;
        //Use local endian for the words
        constexpr int endian = 0;
        constexpr int nails = 0;
        mpz_export(res.data(), &num_read, order, sizeof(T), endian, nails, plain_->m);

        assert(num_read <= req_size && "bytes read should not exceed requested size");
        assert( (!null_termination || num_read <= res.size() - 1)
                && "at least one element must not have been set for null termination");
        if(null_termination) res[num_read] = '\0';
        return res;
    }

    template<typename T, std::enable_if_t<std::is_arithmetic<T>::value, int> = 0>
    operator T() const
    {
        if(std::is_floating_point<T>::value) {
            return static_cast<T>(mpz_get_d(plain_->m));
        } else {
            return static_cast<T>(mpz_get_ui(plain_->m));
        }
    }

    friend void swap(plaintext& lhs, plaintext& rhs) noexcept
    {
        nothrow_swap(lhs.plain_, rhs.plain_);
    }

private:

    paillier_plaintext_t* plain_ = nullptr;

};

//======================================================================
// Random generator
//======================================================================

constexpr struct {
    plaintext operator()(size_t kappa) const
    {
        static buffered_secure_rng rng;
        auto r = rng.get_rand_bits(kappa);
        return {static_cast<void const*>(r.first), r.second};
    }
} gen_random_plain;

//==========================================================================
//struct paillier_keys
//
//Struct holding a public key for encryption and optionally a private key
//for decryption.
//Default constructed paillier_keys are in an unspecified state.Only the following
//operations are safe on unspecified paillier_keys:
// - move-assign
// - move-from (the move-to object becomes unspecified)
// - destructor call
// - destroy_keys()
// - has_public_key()
// - has_private_key()
//==========================================================================
struct paillier_keys {
    //serializable wrapper around paillier_pubkey_t
    struct public_key {
        operator paillier_pubkey_t* () const
        {
            return key;
        }
        paillier_pubkey_t* key = nullptr;

        friend std::ostream& operator<<(std::ostream& os, public_key const& pubkey)
        {
            os << paillier_pubkey_to_hex(pubkey);
            assert(os.good() && "Stream should not go into error state here");
            return os;
        }

        friend std::istream& operator>>(std::istream& is, public_key& pubkey)
        {
            std::string s;
            is >> s;
            assert(is.good() && "Stream should not go into error state here");
            //Const cast is save, as paillier_pubkey_from_hex does not modify its input
            pubkey.key = paillier_pubkey_from_hex(const_cast<char*>(s.c_str()));
            assert(pubkey.key != nullptr);
            return is;
        }
    };

    //serializable wrapper around paillier_privkey_t
    struct private_key {
        operator paillier_prvkey_t* () const
        {
            return key;
        }
        paillier_prvkey_t* key = nullptr;
    };

    paillier_keys() = default;

    paillier_keys(paillier_keys const&) = delete;

    paillier_keys(paillier_keys&& rhs) noexcept
    {
        swap(*this, rhs);
    }

    paillier_keys& operator=(paillier_keys const&) = delete;

    paillier_keys& operator=(paillier_keys&& rhs) noexcept
    {
        swap(*this, rhs);
        return *this;
    }

    ~paillier_keys()
    {
        destroy_keys();
    }

    paillier_keys(public_key pubkey) : pubkey_{pubkey} {}

    explicit paillier_keys(size_t bits)
    {
        generate_keys(bits);
    }

    public_key const& get_public() const
    {
        return pubkey_;
    }

    private_key const & get_private() const
    {
        return privkey_;
    }

    void generate_keys(size_t bits = paillier_bits)
    {
        assert(!has_public_key() && !has_private_key());
        paillier_keygen(bits, &pubkey_.key, &privkey_.key
                        , &paillier_get_rand_devurandom);

    }

    //destroy paillier_keys object and leave it in an unspecified state
    void destroy_keys()
    {
        if(nullptr != pubkey_.key) {
            paillier_freepubkey(pubkey_.key);
            pubkey_.key = nullptr;
        }
        if(nullptr != privkey_.key) {
            paillier_freeprvkey(privkey_.key);
            privkey_.key = nullptr;
        }
    }

    //check if this instance has a valid public key
    bool has_public_key() const
    {
        return nullptr != pubkey_.key;
    }

    //check if this instance has a valid private key
    bool has_private_key() const
    {
        return nullptr != privkey_.key;
    }

private:
    public_key pubkey_{nullptr};
    private_key privkey_{nullptr};

    friend void swap(paillier_keys& lhs, paillier_keys& rhs) noexcept
    {
        nothrow_swap(lhs.pubkey_.key, rhs.pubkey_.key);
        nothrow_swap(lhs.privkey_.key, rhs.privkey_.key);
    }

    template<key_id> friend struct paillier;

};

static_assert(std::is_nothrow_move_constructible<paillier_keys>::value
              , "paillier::paillier_keys's move constructor should not throw");

static_assert(std::is_nothrow_move_assignable<paillier_keys>::value
              , "paillier::paillier_keys's move assignment should not throw");

//Basic class for functionality that is not dependent on key id
struct basic_paillier {

    //Check if a paillier value is operable
    bool operable() const
    {
        return  nullptr != cipher_;
    }

    //Get the length in bytes needed to store a paillier cipher
    size_t get_cipher_len() const
    {
        assert(operable());
        return bytes_to_hold_bitlen(mpz_sizeinbase(cipher_->c, 2));
    }

protected:

    basic_paillier() = default;

    basic_paillier(basic_paillier const&) = default;

    basic_paillier(basic_paillier&&) = default;

    basic_paillier& operator=(basic_paillier const&) = default;

    basic_paillier& operator=(basic_paillier&&) = default;

    basic_paillier(paillier_ciphertext_t* cipher)
        : cipher_{cipher} {}

    ~basic_paillier()
    {
        if(nullptr != cipher_ ) paillier_freeciphertext(cipher_);
    }

    paillier_ciphertext_t* cipher_ = nullptr;
};

//==============================================================================
//struct paillier
//
//Wrapper around a Paillier Cryptosystem.
//The template parameter ID of type key_id is used to differentiate paillier
//values encrypted with possibly different paillier_keys at compile-time, such that
//addition of those values results in a compile-time error.
//The template parameter must be a constexpr cstring variable that also has
//linkage, like a global variable or a constexpr static class member:
//  constexpr key_id foo = "foo";   //global variable
//  void func()
//  {
//      paillier<foo> ...
//   }
//Be also aware that constexpr global variables are not guaranteed to retain
//the same address accross different translation units, so it is possible that
//paillier<foo> in file bar.cpp and paillier<foo> in file baz.cpp might be two
//different types. If a key_id foo needs to be used accross several translation
//units one should define foo as a static member variable of a class or as inline
//variable (since C++17).
//The key of a paillier<ID> is lazily initialized, i.e. the first time an operation
//needs a valid key. Examples include the operation paillier<ID>::encrypt and
//paillier<ID>::get_keys using a randomly generated keypair using
//paillier_bits as the security parameter.
//It is also possible to explicitly initialize the keypair by calling
//paillier<ID>::init_keys using an already generated paillier<ID>::paillier_keys as
//input parameter, a custom security parameter or a public key (in the latter
//case no decryption is possible).
//The key can be deinitialized by calling paillier<ID>::deinit(). In that case
//all living paillier<ID> values transition to an unspecified state. Note that an
//unspecified state is different to the inoperable state in that not even
//operable() might safely be called.
//Freshly encrypted paillier instances are always operable, reflected by the
//method paillier::operable() returning true. All provided methods are
//always safe on operable instances.
//Default constructed instances - i.e. those that were not returned by
//either overload of method paillier<ID>::encrypt - are inoperable,
//reflected by the method operable() returning false
//When moving an instance, the instance moved from will become
//inoperable (and not unspecified).
//Only methods safe on paillier values in unspecified stat,e plus the method
//operable() are safe on an inoperable instances.
//==============================================================================
template<key_id ID>
struct paillier : basic_paillier {

    paillier() = default;

    paillier(paillier const& rhs) = delete;

    paillier(paillier&& rhs) noexcept
    {
        swap(*this, rhs);
    }

    paillier& operator=(paillier const& rhs) = delete;

    paillier& operator=(paillier&& rhs) noexcept
    {
        swap(*this, rhs);
        return *this;
    }

    ~paillier() = default;

    paillier(plaintext plain)
        : paillier{paillier::encrypt(plain)} {}

    paillier& operator+=(paillier const& rhs)
    {
        assert(operable() && rhs.operable());
        paillier_mul(k_->get_public(), cipher_, cipher_, rhs.cipher_);
        return *this;
    }

    paillier& operator*=(plaintext const& rhs)
    {
        assert(operable());
        paillier_exp(k_->get_public(), cipher_, cipher_, rhs);
        return *this;
    }

    //Reserve space for a paillier value.
    static paillier encrypt_zero()
    {
        return {paillier_create_enc_zero()};
    }

    //decrypt a paillier value to a plaintext
    friend plaintext decrypt(paillier const& p)
    {
        assert(nullptr != k_
               && "k_ must be initialized");
        assert(k_->has_private_key()
               && "k_ must have a valid private key when decrypting");
        assert(p.operable()
               && "input paillier value must be operable to decrypt");
        return paillier_dec(nullptr, k_->get_public()
                            , k_->get_private(), p.cipher_);
    }

    //Encrypt a plaintext to a paillier value
    static paillier encrypt(plaintext const& plain)
    {
        init_keys_();
        assert(nullptr != k_ && "k_ should be initialized");
        assert(k_->has_public_key());
        return paillier{paillier_enc(
                            nullptr
                            , paillier::k_->get_public()
                            , plain
                            , paillier_get_rand_devurandom
                        )};
    }

    //Generate a random paillier value
    static paillier gen_random(size_t bitlen)
    {
        size_t len = bytes_to_hold_bitlen(bitlen);
        auto buf = std::make_unique<unsigned char[]>(len);
        paillier_get_rand_devurandom(static_cast<void*>(buf.get()), len);
        return paillier::encrypt(plaintext{static_cast<void const*>(buf.get()), len});
    }

    static paillier_keys const& init_keys(size_t bits = paillier_bits)
    {
        return init_keys(paillier_keys{bits});
    }

    static paillier_keys const& init_keys(paillier_keys&& k)
    {
        assert(k_ == nullptr && "Should not initialize paillier_keys twice");
        return init_keys_(std::addressof(k));
    }

    static void deinit_keys()
    {
        assert(k_ != nullptr && "should not deinit an uninitialized key");
        k_->destroy_keys();
        assert(!k_->has_public_key()
               && "k_ should be successfully destructed at this point");
        k_ = nullptr;
        serialization_length_  = 0;
        assert(k_ == nullptr && "k_ should be set to nullptr at end of function");
    }

    static paillier_keys const& get_keys()
    {
        init_keys_();
        assert(nullptr != k_);
        assert(k_->has_public_key());
        return *k_;
    }

private:

    using basic_paillier::basic_paillier;

    static paillier_keys* k_;
    static size_t serialization_length_;

    static inline paillier_keys get_init_keys_(paillier_keys* movable)
    {
        return nullptr != movable ? std::move(*movable) : paillier_keys{paillier_bits};
    }

    static paillier_keys const& init_keys_(paillier_keys* movable = nullptr)
    {
        static paillier_keys k = get_init_keys_(movable);
        //If k is not valid (e.g. because it was destroyed)
        if(!k.has_public_key()) {
            k = get_init_keys_(movable);
        }
        assert(k.has_public_key()
               && "k should have been correctly initialized here"
               " (i.e. k should have at least a valid public key)");
        k_ = std::addressof(k);
        serialization_length_ = bytes_to_hold_bitlen(2*k.get_public().key->bits);
        return k;
    }

    //We keep swap in derived class, as we dont want to swap
    //paillier values with different keys
    friend void swap(paillier& lhs, paillier& rhs) noexcept
    {
        nothrow_swap(lhs.cipher_, rhs.cipher_);
    }

    inline friend paillier operator+(paillier const& lhs, paillier const& rhs)
    {
        assert(lhs.operable() && rhs.operable());
        paillier res = encrypt_zero();
        paillier_mul(k_->get_public(), res.cipher_,  lhs.cipher_, rhs.cipher_);
        return res;
    }

    inline friend paillier operator*(paillier const& lhs, plaintext const& rhs)
    {
        assert(lhs.operable());
        paillier res = encrypt_zero();
        paillier_exp(k_->get_public(), res.cipher_, lhs.cipher_, rhs);
        return res;
    }

    friend std::ostream& operator<<(std::ostream& os, paillier const& p)
    {
        //Do not serialize inoperable values.
        //Serializing inoperable paillier values is usually not safe in terms
        //of security as it often leaks information about the local control flow
        if(!p.operable()) return os;
        size_t len = serialization_length_;
        assert(p.operable() && "p must be operable here");
        assert(0 < len && "The size of the serialization must not be 0");
        void* cipher_mem = paillier_ciphertext_to_bytes(len, p.cipher_);
        assert(cipher_mem != nullptr);
        os.write(static_cast<char const*>(cipher_mem), len);
        assert(os.good() && "os should not go into error state here");
        free(cipher_mem);
        return os;
    }

    friend std::istream& operator>>(std::istream& is, paillier& p)
    {
        assert(nullptr != k_
               && "k_ must be initialized before deserializing a paillier value");
        size_t len = serialization_length_;
        assert(0 < len && "The size of the serialization must not be 0");

        auto mem = std::make_unique<char[]>(len);
        assert(nullptr != mem && "mem must not be nullptr");
        is.read(mem.get(), len);
        assert(is.good() && "Stream should not go in error state");
        assert(size_t(is.gcount()) == len && "Exactly len bytes must have been read");

        paillier_ciphertext_t* c =
            paillier_ciphertext_from_bytes(static_cast<void*>(mem.get()), len);
        assert(c != nullptr);
        p = paillier{c};
        assert(p.operable());
        return is;
    }

};

//======================================================================
// Implementation
//======================================================================

template<key_id ID>
paillier_keys* paillier<ID>::k_ = nullptr;

template<key_id ID>
size_t paillier<ID>::serialization_length_ = 0;

static_assert(std::is_nothrow_move_constructible<plaintext>::value
              , "plaintext's move constructor should not throw");

static_assert(std::is_nothrow_move_assignable<plaintext>::value
              , "plaintext's move assignment should not throw");

template<key_id ID>
paillier<ID> operator+(paillier<ID> const&, paillier<ID> const&);

template<key_id ID>
paillier<ID> operator+(paillier<ID>&&, paillier<ID> const&);

template<key_id ID>
paillier<ID> operator+(paillier<ID> const&, paillier<ID>&&);

template<key_id ID>
paillier<ID> operator+(paillier<ID>&&, paillier<ID>&&);


template<key_id ID>
paillier<ID> operator*(paillier<ID> const&, plaintext const&);

template<key_id ID>
paillier<ID> operator*(plaintext const& lhs, paillier<ID> const& rhs);

template<key_id ID>
paillier<ID> operator*(paillier<ID>&& lhs, plaintext const& rhs);

template<key_id ID>
paillier<ID> operator*(plaintext const& lhs, paillier<ID>&& rhs);

std::ostream& operator<<(std::ostream&, plaintext const&);
std::istream& operator>>(std::istream&, plaintext&);

template<key_id ID>
std::ostream& operator<<(std::ostream&, paillier<ID> const&);
template<key_id ID>
std::istream& operator>>(std::istream&, paillier<ID>&);


template<key_id ID>
inline paillier<ID> operator+(paillier<ID>&& lhs, paillier<ID> const& rhs)
{
    lhs += rhs;
    return std::move(lhs);
}

template<key_id ID>
inline paillier<ID> operator+(paillier<ID> const& lhs, paillier<ID>&& rhs)
{
    return std::move(rhs) + lhs;
}

template<key_id ID>
inline paillier<ID> operator+(paillier<ID>&& lhs, paillier<ID>&& rhs)
{
    return std::move(lhs) + rhs;
}

template<key_id ID>
inline paillier<ID> operator*(plaintext const& lhs, paillier<ID> const& rhs)
{
    return rhs * lhs;
}

template<key_id ID>
inline paillier<ID> operator*(paillier<ID>&& lhs, plaintext const& rhs)
{
    lhs *= rhs;
    return std::move(lhs);
}

template<key_id ID>
inline paillier<ID> operator*(plaintext const& lhs, paillier<ID>&& rhs)
{
    return std::move(rhs) * lhs;
}

} //namespace mpo19
#endif  //MPO19_PAILLIER_WRAPPER_12062020
