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

#include <mpo19/channel.hpp>
#include <mpo19/paillier_wrapper.hpp>
#include <mpo19/point.hpp>
#include <chrono>
#include <thread>
#include <exception>

#define BOOST_TEST_MODULE channel_test
#include <boost/test/included/unit_test.hpp>
#include <boost/test/data/test_case.hpp>
#include <boost/test/data/monomorphic.hpp>

using namespace mpo19;
namespace bdata = boost::unit_test::data;

constexpr int m = 1000;
constexpr key_id client_key = "CLIENT";
constexpr key_id server_peer_key = "SERVER";

bool equal(plaintext const& lhs, plaintext const& rhs)
{
    std::vector<char> dec_lhs = lhs;
    std::vector<char>dec_rhs = rhs;
    assert(!dec_lhs.empty() && !dec_rhs.empty());
    return strcmp(dec_lhs.data(), dec_rhs.data()) == 0;
}

template<key_id ID>
bool equal(paillier<ID> const& lhs, paillier<ID> const& rhs)
{
    return equal(decrypt(lhs), decrypt(rhs));
}

BOOST_DATA_TEST_CASE(channel_paillier_test
                     , bdata::random(0, m) ^ bdata::random(0, m) ^ bdata::xrange(10)
                     , s1, s2, idx)
{
    using namespace std::chrono_literals;
    (void) idx;

    boost::asio::io_context io_context;
    channel chan{io_context};

    std::thread server{[]()
    {
        boost::asio::io_context io_context;
        channel chan{io_context};
        auto peer = chan.open(7766);
        if(peer.address().to_string() != "127.0.0.1")
            throw std::runtime_error("Not connected to localhost as expected!");
        paillier_keys::public_key pub;
        chan >> pub >> channel::seperator;
        paillier<server_peer_key>::init_keys(pub);

        paillier<server_peer_key> p1, p2;
        chan >> p1 >> p2;
        chan << p1 + p2 << channel::send;

        int s1;
        chan.get(s1);
        int s2 = chan.get<int>();
        chan.put(s1 + s2);
        chan.transmit();
        chan.put(chan.get<paillier<server_peer_key>>()
                 + chan.get<paillier<server_peer_key>>())
        .transmit();
        std::vector<char> hello = chan.get(strlen("Hello")+1);
        hello.resize(hello.size() + strlen(" World!") + 1);
        chan.put(strcat(hello.data(), " World!")).transmit();

        paillier<server_peer_key>::deinit_keys();
    }};

    //Try to connect to server
    for(int i = 0; i != 1000; ++i) {
        boost::system::error_code ec;
        auto peer = chan.open("localhost", "7766", ec);
        if(ec) {
            std::this_thread::sleep_for(10ms);
        }
        //break and skip else statement (python-like for-else)
        else {
            BOOST_TEST(peer.address().to_string() == "127.0.0.1");
            goto CONNECTED;
        }
    }
    /*else*/{
        throw std::runtime_error("Could not connect to server");
    }
CONNECTED:

    paillier<client_key> p1 = paillier<client_key>::gen_random(40);
    paillier<client_key> p2 = paillier<client_key>::gen_random(40);
    paillier<client_key> result;

    chan << paillier<client_key>::get_keys().get_public()<< channel::seperator
         << p1 << p2 << channel::send;
    chan >> result;

    BOOST_TEST(equal(p1 + p2, result));

    chan.put(s1);
    chan.put(s2);
    chan.transmit();
    BOOST_TEST(s1 + s2 == chan.get<int>());

    paillier<client_key> p_s1, p_s2;
    p_s1 = p_s1.encrypt(s1);
    p_s2 = p_s2.encrypt(s2);
    chan.put(p_s1).put(p_s2).transmit();
    paillier<client_key> p_res;
    chan.get(p_res);
    BOOST_TEST(equal(p_res, p_s1 + p_s2));

    chan.put("Hello").transmit();
    std::vector<char> hello_world(strlen("Hello World!") + 1);
    chan.get(hello_world.data(), hello_world.size());

    BOOST_TEST(strcmp(hello_world.data(), "Hello World!") == 0);

    server.join();
    paillier<client_key>::deinit_keys();
}

BOOST_DATA_TEST_CASE(channel_point_test
                     , bdata::xrange(10)
                     , idx)
{
    using namespace std::chrono_literals;
    (void) idx;

    boost::asio::io_context io_context;
    channel chan{io_context};

    std::thread server{[]()
    {
        boost::asio::io_context io_context;
        channel chan{io_context};
        auto peer = chan.open(7766);
        if(peer.address().to_string() != "127.0.0.1")
            throw std::runtime_error("Not connected to localhost as expected!");
        paillier_keys::public_key pub;
        chan >> pub >> channel::seperator;
        paillier<server_peer_key>::init_keys(pub);
        point<paillier<server_peer_key>> p, q;
        chan >> p >> q;
        chan << p + q << channel::send;
        paillier<server_peer_key>::deinit_keys();
    }};

    //Try to connect to server
    for(int i = 0; i != 1000; ++i) {
        boost::system::error_code ec;
        auto peer = chan.open("localhost", "7766", ec);
        if(ec) {
            std::this_thread::sleep_for(10ms);
        }
        //break and skip else statement (python-like for-else)
        else {
            BOOST_TEST(peer.address().to_string() == "127.0.0.1");
            goto CONNECTED;
        }
    }
    /*else*/ {
        throw std::runtime_error("Could not connect to server");
    }
CONNECTED:

    size_t dim = 3 + idx;
    point<paillier<client_key>> p(dim);
    point<paillier<client_key>> q(dim);

    for(size_t i = 0; i != dim; ++i) {
        p[i] = paillier<client_key>::gen_random(40);
        q[i] = paillier<client_key>::gen_random(40);
    }
    point<paillier<client_key>> result, sum = p + q;
    chan << paillier<client_key>::get_keys().get_public()
         << channel::seperator << p << q << channel::send;
    chan >> result;
    server.join();

    for(size_t i = 0; i != dim; ++i) {
        equal(sum[i], result[i]);
    }

}
