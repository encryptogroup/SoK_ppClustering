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
#include <iostream>

namespace mpo19
{

uint64_t BYTES_SENT = 0;
uint64_t BYTES_RECV = 0;

uint64_t bytes_sent() {
    return BYTES_SENT;
}

uint64_t bytes_received() {
    return BYTES_RECV;
}

void reset_counters() {
    BYTES_RECV = 0;
    BYTES_SENT = 0;
}

channel::channel(boost::asio::io_context& io_context)
    : io_context_{io_context}
    , socket_{io_context}
    , is_{&buf_}, os_{&buf_} {}


boost::asio::ip::tcp::endpoint channel::open(std::string const& destination
        , std::string const& service, boost::system::error_code& ec)
{
    using namespace boost::asio;
    using ip::tcp;
    tcp::resolver resolver(io_context_);
    tcp::resolver::results_type endpoints =
        resolver.resolve(destination, service, ec);
    if(!ec)
        return connect(socket_, endpoints, ec);
    else
        return tcp::endpoint{};
}

boost::asio::ip::tcp::endpoint channel::open(unsigned short port)
{
    using boost::asio::ip::tcp;
    tcp::endpoint peer{tcp::v4(), port};
    tcp::acceptor acceptor{io_context_, peer};
    acceptor.set_option(tcp::acceptor::reuse_address(true));
    acceptor.accept(socket_, peer);
    return peer;
}

void channel::close()
{
    if(0 < buf_.size()) *this << channel::send;
    socket_.close();
}

void channel::send_length_encoding_(size_t len)
{
    using namespace boost::asio;
    unsigned char len_buf[encode_len_buffer_size];
    //last is last element written by encode_len_
    size_t last = encode_len(len_buf, len);
    //Send encode of length of input to be sent
    write(socket_, buffer(len_buf, last));
}

void channel::read_some_length_encoding()
{
    size_t bytes_read = socket_.read_some(buf_.prepare(encode_len_buffer_size));
    buf_.commit(bytes_read);
    assert(bytes_read == buf_.size());
}

void channel::prepare_input_buffer_()
{
    using namespace boost::asio;

    if(0 == buf_.size()) {
        read_some_length_encoding();
        size_t len = decode_len([this](unsigned char* buf, size_t len) {
            if(0 == buf_.size()) read_some_length_encoding();
            assert(0 < buf_.size() && "buf_ should never be empty here");
            buf_.sgetn(reinterpret_cast<char*>(buf), len);
            assert(buf_.size() < encode_len_buffer_size);
        });
        assert(buf_.size() <=
               encode_len_buffer_size - words_to_hold_bitlen(ceil_log2(len), CHAR_BIT-1)
               && "buf_ size should be smaller or equal to encode_len_buffer_size minus"
               " the bytes needed needed to encode len");
        size_t bytes_read = read(socket_, buf_.prepare(len - buf_.size()), transfer_all());
        //Move received data into input buffer
        buf_.commit(bytes_read);
        BYTES_RECV += len;
        assert(len == buf_.size() && "input sequence must contain exactly len bytes");
    }
}

void channel::transmit()
{
    using boost::asio::transfer_all;

    // NOTE: We are not including the length of encoding the
    // message-length since this is an implementatio detail.
    BYTES_SENT += buf_.size();

    //First send length encoding
    send_length_encoding_(buf_.size());
    //Send input
    size_t len = boost::asio::write(socket_, buf_.data(), transfer_all());
    //Remove written data from input buffer
    buf_.consume(len);
    assert(buf_.size() == 0
           && "There should be nothing in the input sequence now");
}
} //namespace mpo19
