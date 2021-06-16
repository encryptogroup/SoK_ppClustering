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

#ifndef MPO19_CHANNEL_16062020
#define  MPO19_CHANNEL_16062020

#include <mpo19/utility.hpp>

#include <boost/asio/io_context.hpp>
#include <boost/asio/buffered_stream.hpp>
#include <boost/asio/streambuf.hpp>
#include <boost/asio/ip/tcp.hpp>
#include <boost/asio/connect.hpp>
#include <boost/asio/read.hpp>
#include <boost/asio/write.hpp>
#include <boost/endian/conversion.hpp>

#include <iosfwd>
#include <type_traits>

namespace mpo19
{

// Having these as static class functions in `channel` leads to issues due to
// seperate translation units across source files.
uint64_t bytes_sent();
uint64_t bytes_received();
void reset_counters();

struct channel {

    channel(boost::asio::io_context& io_context);

    boost::asio::ip::tcp::endpoint open(std::string const& destination
                                        , std::string const& service, boost::system::error_code& ec);

    boost::asio::ip::tcp::endpoint open(unsigned short port);

    void close();

private:
    struct send_t {};
    struct seperator_t {
        static constexpr char val = ' ';
        constexpr operator char() const
        {
            return val;
        }
    };
public:

    channel& put(char const* s, std::streamsize n)
    {
        os_.write(s, n);
        return *this;
    }

    channel& get(char* s, std::streamsize n)
    {
        prepare_input_buffer_();
        is_.read(s, n);
        return *this;
    }

    std::vector<char> get(std::streamsize n)
    {
        std::vector<char> res(n);
        get(res.data(), n);
        return res;
    }

    template<typename T, std::enable_if_t<std::is_arithmetic<T>::value, int> = 0>
    channel& put(T const& t)
    {
        T network_order_t = boost::endian::native_to_big(t);
        return put(reinterpret_cast<char*>(&network_order_t), sizeof(T));
    }

    template<typename T, std::enable_if_t<std::is_arithmetic<T>::value, int> = 0>
    channel& get(T& t)
    {
        get(reinterpret_cast<char*>(&t), sizeof(T));
        boost::endian::big_to_native_inplace(t);
        return *this;
    }

    template<typename T
             , std::enable_if_t<
                 !std::is_arithmetic<T>::value
                 && is_serializable<T>(), int> = 0>
    channel& put(T const& t)
    {
        return (*this) << t;
    }

    template<typename T
             , std::enable_if_t<
                 !std::is_arithmetic<T>::value
                 && is_deserializable<T>(), int> = 0>
    channel& get(T& t)
    {
        return (*this) >> t;
    }

    template<typename T>
    T get()
    {
        T t;
        get(t);
        return t;
    }

    void transmit();

    static constexpr send_t send{};
    static constexpr seperator_t seperator{};

    template<typename T>
    friend channel& operator<<(channel& chan, T&& t)
    {
        chan.os_ << t;
        return chan;
    }

    friend channel& operator<<(channel& chan, seperator_t)
    {
        chan.os_ << seperator_t::val;
        return chan;
    }

    friend channel& operator<<(channel& chan, send_t)
    {
        chan.transmit();
        return chan;
    }

    template<typename T>
    friend channel& operator>>(channel& chan, T& t)
    {
        using namespace boost::asio;
        chan.prepare_input_buffer_();
        chan.is_ >> t;
        return chan;
    }

    friend channel& operator>>(channel& chan, seperator_t)
    {
        //Consume seperator
        chan.buf_.consume(1);
        return chan;
    }

private:
    boost::asio::io_context& io_context_;
    boost::asio::ip::tcp::socket socket_;
    boost::asio::streambuf buf_;
    std::istream is_;
    std::ostream os_;

    void send_length_encoding_(size_t len);
    void prepare_input_buffer_();
    void read_some_length_encoding();
};
}

#endif // MPO19_CHANNEL_16062020
