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

#ifdef NDEBUG
#define BOOST_DISABLE_ASSERTS
#endif

#include <mpo19/typedefs.hpp>
#include <mpo19/utility.hpp>
#include <mpo19/point.hpp>
#include <mpo19/paillier_wrapper.hpp>
#include <mpo19/linkage.hpp>
#include <mpo19/channel.hpp>
#include <mpo19/random.hpp>
#include <mpo19/cluster.hpp>
#include <mpo19/pca.hpp>
#include <mpo19/clustering.hpp>
#include <mpo19/pcure.hpp>

#include <cassert>
#include <algorithm>
#include <numeric>
#include <type_traits>
#include <cstdlib>
#include <thread>
#include <chrono>
#include <ratio>
#include <iomanip>

#include <boost/range/adaptors.hpp>
#include <boost/range/numeric.hpp>
#include <boost/range/algorithm_ext.hpp>
#include <boost/multi_array.hpp>
#include <boost/phoenix.hpp>
#include <boost/hana.hpp>


#ifdef MPO19_TESTING_PCA_
#include <pca_test.hpp>
#endif //MPO19_TESTING_PCA_

using namespace mpo19;

//======================================================================
//Helpers for  encryption/decryption
//======================================================================
template<key_id ID>
struct enc_t {

    paillier<ID> operator()(plaintext const& plain) const
    {
        return paillier<ID>::encrypt(plain);
    }

    point<paillier<ID>> operator()(point<int64_t> const& in) const
    {
        return apply_on_elements([&](auto const& plain) {
            return paillier<ID>::encrypt(plain);
        }, in);
    }

    point<paillier<ID>> operator()(point<plaintext> const& in) const
    {
        return apply_on_elements([&](auto const& plain) {
            return paillier<ID>::encrypt(plain);
        }, in);
    }

    template<typename Range, std::enable_if_t<is_range<Range>(), int> = 0>
    auto operator()(Range const& range) const
    {
        using boost::adaptors::transformed;
        return range | transformed([&](auto const& e) {
            return (*this)(e);
        });
    }

};

template<key_id ID>
constexpr enc_t<ID> enc{};

constexpr struct dec_t {

    template<key_id ID>
    plaintext operator()(paillier<ID> const& plain) const
    {
        return decrypt(plain);
    }

    template<key_id ID>
    point<plaintext> operator()(point<paillier<ID>> const& in) const
    {
        return apply_on_elements([&](auto const& plain) {
            return decrypt(plain);
        }, in);
    }

    template<typename Range, std::enable_if_t<is_range<Range>(), int> = 0>
    auto operator()(Range const& range) const
    {
        using boost::adaptors::transformed;
        return range | transformed([this](auto const& e) {
            return (*this)(e);
        });
    }

    template<typename InputIt>
    auto operator()(InputIt const& first, InputIt const& last) const
    {
        return (*this)(boost::make_iterator_range(first, last));
    }

} dec;

//======================================================================
// Dist function
//======================================================================

struct dist_t {
    size_t operator()(point<int64_t> const& lhs
                      , point<int64_t> const& rhs) const
    {
        using namespace boost::phoenix::placeholders;
        assert(lhs.dimension() == rhs.dimension());
        return apply_on_elements(arg1*arg1, lhs - rhs).fold();
    }

    template<key_id ID>
    paillier<ID> operator()(point<int64_t> const& p
                            , std::array<point<paillier<ID>>*, 3> const& qs) const
    {
        namespace phoenix = boost::phoenix;
        using namespace boost::phoenix::placeholders;
        using boost::hana::compose;
        point<paillier<ID>>& minus_2q = *(qs[1]);
        point<paillier<ID>>& q_sqr = *(qs[2]);

        auto square_encrypt_k = compose(enc<ID>, arg1 * arg1);

        point<paillier<ID>> p_sqr = apply_on_elements(square_encrypt_k, p);
        //Compute -2q*p
        point<paillier<ID>> minus_2qp = apply_on_elements(arg1 * arg2, p, minus_2q);
        //Compute p² + (-2q)p + q²
        auto res = (std::move(p_sqr) + std::move(minus_2qp) + q_sqr).fold();
        return res;
    }
};

//======================================================================
// Communicator between P1 and P2
//======================================================================

struct communicator {

    communicator(unsigned short port)
        : endpoint_{chan_.open(port)} {}

    communicator(std::string const& destination
                 , unsigned short port
                 , int trials = 10
                 , std::chrono::milliseconds interval =
                 std::chrono::milliseconds{100})
    {
        for(int i = 0; i != trials; ++i) {
            boost::system::error_code ec;
            auto endpoint = chan_.open(destination, std::to_string(port), ec);
            if(ec) {
                std::this_thread::sleep_for(interval);
            } else {
                endpoint_ = std::move(endpoint);
                return;
            }
        }
        throw std::runtime_error("Could not establish a connection");
    }

    size_t exchange(size_t nx)
    {
        using mpo19::channel;
        chan_ << nx << channel::seperator << channel::send;
        size_t res_nx;
        chan_ >> res_nx >> channel::seperator;
        return res_nx;
    }

    paillier_keys exchange(paillier_keys const& pkx)
    {
        using mpo19::channel;
        chan_ << pkx.get_public() << channel::seperator << channel::send;
        paillier_keys::public_key res_pkx;
        chan_ >> res_pkx >> channel::seperator;
        return res_pkx;
    }

    void synchronize()
    {
        chan_.put('\0').transmit();
        chan_.get<char>();
    }

    void close()
    {
        chan_.close();
    }

public:

    template<typename... Ts>
    void send(Ts const&... ts)
    {
        using mpo19::channel;
        FOR_EACH_VARIADIC( (send_impl_(ts)) );
        chan_ << channel::send;
    }

    template<typename Range, std::enable_if_t<is_range<Range>(), int> = 0>
    void send_impl_(Range const& range)
    {
        for(auto const& e : range) {
            send_impl_(e);
        };
    }

    template<typename T>
    void send_impl_(std::vector<T> const& vec)
    {
        chan_.put(vec.size());
        for(T const& t : vec) {
            send_impl_(t);
        }
    }

    template<typename T, std::enable_if_t<!is_range<T>(), int> = 0>
    void send_impl_(T const& t)
    {
        chan_.put(t);
    }

    template<typename... Ts>
    void receive(Ts&... ts)
    {
        FOR_EACH_VARIADIC( (receive_impl_(ts)) );
    }

    template<typename Range, std::enable_if_t<is_range<Range>(), int> = 0>
    void receive_impl_(Range& range)
    {
        for(auto it = range.begin(); it != range.end(); ++it) {
            auto&& e = *it;
            receive_impl_(e);
        };
    }

    template<typename T>
    void receive_impl_(std::vector<T>& vec)
    {
        size_t sz;
        chan_.get(sz);
        vec.clear();
        vec.reserve(sz);
        for(size_t i = 0; i != sz; ++i) {
            T t;
            receive_impl_(t);
            vec.emplace_back(std::move(t));
        }
    }

    template<typename T, std::enable_if_t<!is_range<T>(), int> = 0>
    void receive_impl_(T& t)
    {
        chan_.get(t);
    }

    std::string get_peer_ip() const
    {
        return endpoint_.address().to_string();
    }

    unsigned short get_peer_port() const
    {
        return endpoint_.port();
    }

private:
    boost::asio::io_context io_context_;
    mpo19::channel chan_{io_context_};
    boost::asio::ip::tcp::endpoint endpoint_;

};

template<typename T, typename InterLinkagePolicy, typename ClearClusterDist>
struct intercluster_dist_t : public InterLinkagePolicy {

    template<typename... Args>
    intercluster_dist_t(
        std::vector<point<T>> const& values
        , communicator& comm
        , size_t kappa
        , ClearClusterDist& clear_cluster_dist
        , Args&&... args)
        : InterLinkagePolicy{std::forward<Args>(args)...}
        , values_{values}
        , dist_{}
        , comm_{comm}
        , kappa_{kappa}
        , clear_cluster_dist_{clear_cluster_dist} {}

    auto operator()(cluster const& c1, cluster const& c2)
    {
        return clear_cluster_dist_(c1, c2);
    }

    template<key_id ID>
    auto operator()(cluster const& c1, std::array<std::vector<point<paillier<ID>>>*, 3> c2)
    {
        assert(std::all_of(std::begin(c2), std::end(c2), [&](std::vector<point<paillier<ID>>> const* v) {
            return v != nullptr;
        }));
        auto c1_indices = get_indices(c1);

        size_t c1_indices_size = c1_indices.size();
        size_t c2_indices_size = c2[0]->size();

        assert(0 < c1_indices_size);
        assert(0 < c2_indices_size);
        assert(std::all_of(std::begin(c2), std::end(c2), [&](std::vector<point<paillier<ID>>> const* v) {
            return c2_indices_size== v->size();
        }));

        plaintext res;
        InterLinkagePolicy::invalidate(res);

        comm_.send(c1_indices_size * c2_indices_size);

        for(size_t i = 0; i != c1_indices_size; ++i) {
            for(size_t j = 0; j != c2_indices_size; ++j) {
                plaintext r = gen_random_plain(kappa_);
                std::array<point<paillier<ID>>*, 3>
                                            arr{std::addressof((*c2[0])[j]), std::addressof((*c2[1])[j]), std::addressof((*c2[2])[j])};
                paillier<ID> p_enc = dist_(values_[c1_indices[i]], arr);
                p_enc += enc<ID>(r);
                comm_.send(p_enc);
                res = InterLinkagePolicy::cluster_distance(res, r);
            }
        }
        paillier<ID> res_enc;
        comm_.receive(res_enc);
        res_enc += -res;
        return res_enc;
    }

    template<key_id ID>
    void calc_intercluster_distance_with_peer()
    {
        plaintext res;
        InterLinkagePolicy::invalidate(res);

        size_t n;
        comm_.receive(n);

        for(size_t i = 0; i != n; ++i) {
            paillier<ID> p_enc;
            comm_.receive(p_enc);
            plaintext p = dec(p_enc);
            res = InterLinkagePolicy::cluster_distance(res, p);
        }

        comm_.send(enc<ID>(res));
    }

private:
    std::vector<point<T>> const& values_;
    dist_t dist_;
    communicator& comm_;
    size_t kappa_;
    ClearClusterDist& clear_cluster_dist_;
};


//Port used for setup, clustering and output phase
constexpr unsigned short default_port = 7766;
//Port used for ArgMin and MaxDist calculation
constexpr unsigned short default_aby_port = 6677;
//Ports must be different, else socket cannot be bound
static_assert(default_port != default_aby_port);

std::string default_destination{"localhost"};

struct P1::impl_ {

    //pk for P1, pk_dash for P2
    static constexpr key_id pk = "P1 - pk ";
    static constexpr key_id pk_dash = "P1 - pk' ";
    static constexpr e_role role = ROLE_P1;

    impl_(unsigned short port)
        : comm{port}, port{port}
    {
        paillier<pk_dash>::init_keys(comm.exchange(paillier<pk>::init_keys(paillier_bits)));
    }

    ~impl_()
    {
        paillier<pk_dash>::deinit_keys();
        paillier<pk>::deinit_keys();
    }

    // NOTE: Temporary fix to get API working
    std::string get_aby_ip() {
        // Party 1 is server (dependent of definition in typedefs.h)
        return "0.0.0.0";
    }

    void init_L(array_t<point<int64_t>>& input)
    {
        namespace phoenix = boost::phoenix;
        namespace range = boost::range;
        using namespace phoenix::placeholders;

        using boost::indices;
        using phoenix::bind;
        using boost::extents;

        index_t n1 = in->n1, n = in->n;
        size_t d = get_dimension(input);

        assert(index_t(input.size()) == n1);

        //We encrypt our input and store it into L
        //Compute array L : L_i= [[p_i]],    i ∈ [1 : n1]
        auto& L = in->L;
        L.reindex(1);

        assert(index_t(L.size()) == n);

        range::copy(enc<pk_dash>(input), L.begin());
        auto other_input = L[indices[n1 < index_range_t{} <= n]];
        comm.receive(other_input);

        //Compute array S: S_i = s_i,    s_i <--$-- {0,1}^κ,    i ∈ [1:n]
        array_t<point<plaintext> > S{extents[n]};
        std::generate(S.begin(), S.end(), [&] {
            point<plaintext> res(d);
            apply_on_elements(res, bind(gen_random_plain, kappa), res);
            return res;
        });
        S.reindex(1);
        //Blind L as: L_i := L_i · [[S_i]]
        multi_for_each(L, enc<pk_dash>(S), arg1 += arg2);
        comm.send(in->pi_1.permute(enc<pk>(S)), in->pi_1.permute(L));
        array_t<point<paillier<pk>>> S_enc{extents[n]};
        comm.receive(S_enc, L);
        //Decrypt S as: S_i := < S_i >,    i ∈ [1 : n]
        range::copy(dec(S_enc), S.begin());
        //Unblind L as: L_i := L_i · [[S_i]]^(-1)
        multi_for_each(L, S, arg1 += bind(enc<pk_dash>, -arg2));
    }

    void init_L(pcure_clusters<int64_t> const& input)
    {
        index_t n1 = in->n1;
        array_t<point<int64_t>> arr{boost::extents[n1]};
        size_t d = get_dimension(input.values);

        assert(index_t(input.partition_dendrograms.size()) == n1);

        for(index_t i = 0; i != n1; ++i) {
            point<int64_t> p(d);
            auto indices = get_indices(input.partition_dendrograms[i]);
            for(index_t idx : indices) {
                p += input.values[idx];
            }
            arr[i] = std::move(p);
        }
        init_L(arr);
    }

    template<typename Dist>
    void init_B(array_t<point<int64_t>>& input
                , table_t<paillier<pk_dash>>& B
                , Dist&& dist)
    {
        namespace phoenix = boost::phoenix;
        namespace range = boost::range;
        namespace adaptors = boost::adaptors;
        using namespace phoenix::placeholders;

        using boost::extents;
        using boost::indices;
        using phoenix::bind;
        using phoenix::cref;

        index_t n1 = in->n1, n2 = in->n2, n = in->n;
        input.reindex(1);

        init_L(input);

        //Receive Matrices H and D from P2
        table_t<point<paillier<pk_dash>>> H{extents[precomputed][n2]};
        //We reindex H to 2, as H[1] won't be inserted in H, but in L[n1+1..n]
        H.reindex(1);

        //We insert H[1][1..n2] directly into L[n1+1..n]
        auto H_1 = H[indices[1][index_range_t{}]];
        auto H_2 = H[indices[2][index_range_t{}]];
        auto H_3 = H[indices[3][index_range_t{}]];

        //H_1[1..n2] maps to L[n1+1..n], so we reindex H_1 to H_3
        //to reflect that fact.
        H_1.reindex(n1+1);
        H_2.reindex(n1+1);
        H_3.reindex(n1+1);

        //We read D directly into B
        auto D = B[indices[n1 < index_range_t{} <= n][n1 < index_range_t{} <= n]];

        //Receive H, D from other party, with the following effects:
        //L_i = H_1,i-n1 ,    i ∈ [n1+1 : n]
        //B_i,j = D_i-n1,j-n1,    i,j ∈ [n1+1 : n]
        comm.receive(H_1, H_2, H_3, D);

        //Compute matrix B:
        assert(index_t(B.size()) == n);
        for(index_t i = 1; i <= index_t(B.size()); ++i) {
            assert(index_t(B[i].size()) == n);
            for(index_t j = 1; j <= index_t(B[i].size()); ++j) {
                //NOTE Optimize by only calculating the lower half of the symmetric matrix
                if(i <= n1 && j <= n1) {
                    //B_i,j = [[dist(p_i, p_j)]],     if i,j ≤ n1
                    B[i][j] = enc<pk_dash>(dist(input[i], input[j]));
                } else if(n1 < i && n1 < j) {
                    //D[1..n2][1..n2] is an alias of B[n1+1..n][n1+1..n], so nothing to do here
                    //if you delete this else if-clause, think about modifying the else clause in that case
                }
                //Compute dist(p_i, q_j), i in {1, n1}, j ∈ {1, n2}
                //It always holds that i ≤ n1 ⊻ j ≤ n1
                else {
                    //P1's input is represented by indices between 1..n1, so
                    //we take the minimum between i and j, where ii ≤ n1 holds
                    index_t ii = std::min(i, j);
                    assert(ii <= n1);
                    //H is represented by indices n1+1..n, so we always take the maximum
                    //where n1 < jj ≤ n holds.
                    //Note that we make use of the reindexing of H to base n1+1 above.
                    index_t jj = std::max(i, j);
                    assert(n1 < jj && j <= n);
                    assert(ii < jj);
                    //We don't need to take care of argument order, as dist(a, b) == dist(b, a)
                    B[i][j] = dist(input[ii]
                    , std::array<point<paillier<pk_dash>>*, precomputed> {
                        &H_1[jj], &H_2[jj], &H_3[jj]
                    });
                }
            }
        }
    }

    template<typename Dist>
    void init_B(pcure_clusters<int64_t>& input
                , table_t<paillier<pk_dash>>& B
                , Dist&& dist)
    {
        using boost::extents;
        using boost::indices;

        index_t n1 = in->n1, n2 = in->n2, n = in->n;
        init_L(input);

        table_t<std::vector<point<paillier<pk_dash>>>> H{extents[precomputed][n2]};
        H.reindex(std::array<size_t, 2> {size_t(1), size_t(n1+1)});

        //We read D directly into B
        auto D = B[indices[n1 < index_range_t{} <= n][n1 < index_range_t{} <= n]];

        comm.receive(H, D);

        //Compute matrix B:
        assert(index_t(B.size()) == n);
        for(index_t i = 1; i <= n; ++i) {
            assert(index_t(B[i].size()) == n);
            for(index_t j = 1; j <= n; ++j) {
                //NOTE Optimize by only calculating the lower half of the symmetric matrix
                if(i <= n1 && j <= n1) {
                    //B_i,j = [[dist(p_i, p_j)]],     if i,j ≤ n1
                    B[i][j] = enc<pk_dash>(
                                  dist(input.partition_dendrograms[i-1]
                                       , input.partition_dendrograms[j-1]));
                } else if(n1 < i && n1 < j) {
                    //D[1..n2][1..n2] is an alias of B[n1+1..n][n1+1..n], so nothing to do here
                    //if you delete this else if-clause, think about modifying the else clause in that case
                }
                //Compute dist(p_i, q_j), i in {1, n1}, j ∈ {1, n2}
                //It always holds that i ≤ n1 ⊻ j ≤ n1
                else {
                    assert(i != j);
                    //P1's input is represented by indices between 1..n1, so
                    //we take the minimum between i and j, where ii ≤ n1 holds
                    index_t ii = std::min(i, j);
                    assert(ii <= n1);
                    //H is represented by indices n1+1..n, so we always take the maximum
                    //where n1 < jj ≤ n holds.
                    //Note that we make use of the reindexing of H to base n1+1 above.
                    index_t jj = std::max(i, j);
                    assert(n1 < jj && j <= n);
                    assert(ii < jj);
                    //We don't need to take care of argument order, as dist(a, b) == dist(b, a)
                    std::array<std::vector<point<paillier<pk_dash>>>*, 3>
                    arr{std::addressof(H[1][jj]), std::addressof(H[2][jj]), std::addressof(H[3][jj])};
                    B[i][j] = dist(input.partition_dendrograms[ii-1], arr);
                }
            }
        }
    }

    template<typename Range>
    size_t get_dimension(Range& input) const
    {
        return (*input.begin()).dimension();
    }

    template<typename Input, typename Dist = dist_t>
    table_t<plaintext>& setup(Input& input, Dist&& dist = dist_t{})
    {
        namespace phoenix = boost::phoenix;
        namespace range = boost::range;
        namespace adaptors = boost::adaptors;
        using namespace phoenix::placeholders;

        using boost::extents;
        using boost::indices;
        using phoenix::bind;
        using phoenix::cref;

        static_assert(precomputed == 3
                      , "Any other amount than 3 precomputed values is not supported");

        index_t n = in->n;

        //================================================================
        // First round
        //
        // Get H, D from P1, then generate S, R, L, B and send these to P2.
        //================================================================

        //R is an alias for Sigma for P1
        table_t<plaintext>& R = in->Sigma;
        //Compute matrix R: R_i,j = r_i,j,    r_i,j  <--$-- {0,1}^κ,    i,j ∈ [1 : n]
        R.reindex(1);
        std::generate(R.data(), R.data() + R.num_elements(), bind(gen_random_plain, kappa));

        table_t<paillier<pk_dash>> B{extents[n][n]};
        B.reindex(1);

        init_B(input, B, dist);

        multi_for_each(B, R, [](paillier<pk_dash>& b, plaintext const& r) {
            b += enc<pk_dash>(r);
        });

        //Encrypt S as: S_i := [S_i]
        //Encrypt R as:R_i,j := [R_i,j]
        //Permute S, L, R and B via a random permutation π_1(n)
        //Send {S, L, R, B}
        comm.send(in->pi_1.permute(enc<pk>(R)), in->pi_1.permute(B));

        table_t<paillier<pk>> R_enc{extents[n][n]};

        comm.receive(R_enc);

        //Decrypt R as: R_i,j := < R_i,j >,    i,j ∈ [1 : n]
        multi_for_each(R, R_enc, arg1 = bind(dec, arg2));

        return in->Sigma;
    }

    output_t output()
    {
        using dendrogram_t = std::vector<cluster>;
        using E_dec_t = std::vector<point<double>>;
        using J_t = std::vector<size_t>;

        constexpr size_t dendrogram_idx = 0;
        constexpr size_t E_dec_idx = 1;
        constexpr size_t J_idx = 2;

        auto& Sigma = in->Sigma;
        auto& L = in->L;
        auto& dendrogram_ = in->dendrogram_;

        auto res = std::make_tuple(dendrogram_t{}, E_dec_t{}, J_t{});

        assert(Sigma.index_bases()[0] == Sigma.index_bases()[1]
               && "Index base for both dimensions of Sigma must be equal");

        index_t n = in->n;
        L.reindex(0);
        assert(L.size() == size_t(n) && "L must be of size n");

        dendrogram_t& dendro = std::get<dendrogram_idx>(res);
        E_dec_t& E_dec = std::get<E_dec_idx>(res);

        //Initialize arrays E,J: E_i= J_i =⊥    i ∈ [1 : n]
        std::vector<point<paillier<pk_dash>>> E;
        E.reserve(n);
        J_t& J = std::get<J_idx>(res);
        J.reserve(n);
        //foreach i=1, ..., n do:
        for(cluster& R_ii : dendrogram_) {
            //if R_i,i encodes a target cluster C_i then
            if(R_ii.is_cluster()) {
                //R_ii is not a leaf, so that means the cluster is a target cluster
                //Only clusters that were merged at least once are target clusters

                //find the index set I_i of points in cluster C_i
                std::vector<index_t> I_i = get_indices(R_ii);
                assert(std::all_of(I_i.begin(), I_i.end(), [&](index_t idx) {
                    return 0 <= idx && idx < n;
                })
                && "All indices idx must be between 0 <= idx < n");

                point<paillier<pk_dash>> r;
                //Set E_i = ∏_{j ∈ I_i} L_j
                if(1 < I_i.size()) {
                    for(size_t j = 1; j != I_i.size(); ++j) {
                        if(j == 1) {
                            r = L[I_i[j - 1]] + L[I_i[j]];
                        } else {
                            r += L[I_i[j]];
                        }
                    }
                } else if(1 == I_i.size()) {
                    //Move is save here, as L[I_i[0]] will be used only once
                    r = std::move(L[I_i[0]]);
                } else {
                    assert(false && "cluster::is_cluster() should not return true on an empty cluster");
                }
                assert(0 < r.dimension());
                E.emplace_back(std::move(r));

                //J_i = |I_i|
                J.emplace_back(I_i.size());
                dendro.emplace_back(std::move(R_ii));
            }
        }

        assert(E.size() == dendro.size());
        assert(J.size() == dendro.size());

        comm.send(E, J);
        E_dec.resize(dendro.size());

        assert(E_dec.size() == dendro.size());

        comm.receive(E_dec);

        return res;
    }

    //Perform initial phase for P1, i.e. exchanging input sizes
    //and public paillier_keys with P2
    struct initial {

        initial(communicator& comm, size_t n1)
            : n1{n1}, n2{comm.exchange(n1)}, n{n1 + n2} {}
        //n1 for P1, n2 for P2, n = n1 + n2
        size_t n1, n2, n;

        table_t<plaintext> Sigma{boost::extents[n][n]};
        array_t<point<paillier<pk_dash>>> L{boost::extents[n]};
        random_permutation pi_1{n};
        std::vector<cluster> dendrogram_;
    };

    mutable communicator comm;
    std::unique_ptr<initial> in;
    unsigned short port;
};

key_id P1::impl_::pk;
key_id P1::impl_::pk_dash;

struct P2::impl_ {

    static constexpr key_id pk = "P2 - pk ";
    static constexpr key_id pk_dash = "P2 - pk' ";
    static constexpr e_role role = ROLE_P2;

    impl_(std::string const& destination, unsigned short port)
        : comm{destination, port}, port{port}
    {
        paillier<pk>::init_keys(comm.exchange(paillier<pk_dash>::init_keys(paillier_bits)));
    }

    ~impl_()
    {

        paillier<pk>::deinit_keys();
        paillier<pk_dash>::deinit_keys();
    }

    // NOTE: Temporary fix to get API working
    std::string get_aby_ip() const {
        // Party 2 is client (dependent of definition in typedefs.h)
        return comm.get_peer_ip();
    }


    void init_L(array_t<point<int64_t>>& input)
    {
        namespace phoenix = boost::phoenix;
        using namespace phoenix::placeholders;

        using boost::indices;
        using phoenix::bind;
        using boost::extents;

        index_t n = in->n;
        size_t d = get_dimension(input);

        assert(input.size() == in->n2);

        comm.send(enc<pk_dash>(input));

        array_t<point<plaintext>> S_dash{extents[n]};

        //Compute array S': S'_i = s'_i,    s'_i <--$-- {0,1}^κ,    i ∈ [1:n]
        std::generate(S_dash.data(), S_dash.data() + S_dash.num_elements(), [&] {
            point<plaintext> res(d);
            apply_on_elements(res, bind(gen_random_plain, kappa), res);
            return res;
        });

        array_t<point<paillier<pk>>> S{extents[n]};
        array_t<point<paillier<pk_dash>>> L{extents[n]};
        comm.receive(S, L);
        //Blind L_i = L_i · [[S'_i]]
        multi_for_each(L, enc<pk_dash>(S_dash), arg1 += arg2);
        //Blind S_i = S_i · [S'_i]
        multi_for_each(S, enc<pk>(S_dash), arg1 += arg2);
        comm.send(in->pi_2.permute(S), in->pi_2.permute(L));
    }

    void init_L(pcure_clusters<int64_t> const& input)
    {
        index_t n2 = in->n2;
        array_t<point<int64_t>> arr{boost::extents[n2]};
        size_t d = get_dimension(input.values);

        assert(index_t(input.partition_dendrograms.size()) == n2);

        for(index_t i = 0; i != n2; ++i) {
            point<int64_t> p(d);
            auto indices = get_indices(input.partition_dendrograms[i]);
            for(index_t idx : indices) {
                p += input.values[idx];
            }
            arr[i] = std::move(p);
        }
        init_L(arr);
    }

    template<typename Dist>
    void init_Sigma(array_t<point<int64_t>>& input, Dist&& dist)
    {
        namespace phoenix = boost::phoenix;
        namespace range = boost::range;
        namespace adaptors = boost::adaptors;
        using namespace phoenix::placeholders;

        using boost::indices;
        using phoenix::bind;
        using phoenix::cref;
        using boost::extents;

        index_t n2 = in->n2;
        input.reindex(1);

        init_L(input);

        //Function to square all elements of a point
        auto square_elements = [](point<int64_t> const& p) {
            using namespace boost::phoenix::placeholders;
            return apply_on_elements(arg1 * arg1, p);
        };

        table_t<point<paillier<pk_dash>>> H{extents[precomputed][n2]};
        H.reindex(1);

        //Initialize H
        //H_1,i = [[q_i]],    i ∈ [1 : n2]
        range::copy(enc<pk_dash>(input)
                    , H[indices[1][index_range_t{}]].begin());
        //H_2,i = [[-2*q_i]],    i ∈ [1 : n2]
        range::copy(enc<pk_dash>(input | adaptors::transformed(-2*arg1))
                    , H[indices[2][index_range_t{}]].begin());
        //H_3,i = [[q_i²]],    i ∈ [1 : n2]
        range::copy(enc<pk_dash>(input | adaptors::transformed(square_elements))
                    , H[indices[3][index_range_t{}]].begin());

        table_t<paillier<pk_dash>> D{extents[n2][n2]};
        D.reindex(1);

        //D_i,j = [[dist(q_i,q_j)]],    i, j ∈ [1 : n2]
        for(index_t i = 1; i <= index_t(D.size()); ++i) {
            for(index_t j = 1; j <= index_t(D[i].size()); ++j) {
                D[i][j] = enc<pk_dash>(dist(input[i], input[j]));
            }
        }

        comm.send(H, D);
    }


    template<typename Dist>
    void init_Sigma(pcure_clusters<int64_t> const& input, Dist&& dist)
    {
        namespace phoenix = boost::phoenix;
        namespace range = boost::range;
        namespace adaptors = boost::adaptors;
        using namespace phoenix::placeholders;

        using boost::extents;
        using boost::indices;
        using phoenix::bind;
        using phoenix::cref;

        index_t n1 = in->n1, n2 = in->n2, n = in->n;
        init_L(input);

        //Function to square all elements of a point
        auto square_elements = [](point<int64_t> const& p) {
            using namespace boost::phoenix::placeholders;
            return apply_on_elements(arg1 * arg1, p);
        };

        auto make_H_cluster = [&](auto&& f, cluster const& c) {
            auto indices = get_indices(c);
            std::vector<point<paillier<pk_dash>>> res;
            res.reserve(indices.size());
            for(index_t i : indices) {
                res.emplace_back(f(input.values[i]));
            }
            return res;
        };

        table_t<std::vector<point<paillier<pk_dash>>>> H{extents[precomputed][n2]};
        H.reindex(1);
        //Initialize H
        //H_1,i = [[q_i]],    i ∈ [1 : n2]
        range::copy(input.partition_dendrograms | adaptors::transformed(
        [&](cluster const& c) {
            return make_H_cluster(enc<pk_dash>, c);
        })
        , H[indices[1][index_range_t{}]].begin());

        //H_2,i = [[-2*q_i]],    i ∈ [1 : n2]
        range::copy(input.partition_dendrograms | adaptors::transformed(
        [&](cluster const& c) {
            return make_H_cluster(bind(enc<pk_dash>, -2*arg1), c);
        })
        , H[indices[2][index_range_t{}]].begin());

        //H_3,i = [[q_i²]],    i ∈ [1 : n2]
        range::copy(input.partition_dendrograms | adaptors::transformed(
        [&](cluster const& c) {
            return make_H_cluster(bind(enc<pk_dash>, bind(square_elements, arg1)), c);
        })
        , H[indices[3][index_range_t{}]].begin());

        table_t<paillier<pk_dash>> D{extents[n2][n2]};
        D.reindex(1);

        //D_i,j = [[dist(q_i,q_j)]],    i, j ∈ [1 : n2]
        for(index_t i = 1; i <= index_t(D.size()); ++i) {
            for(index_t j = 1; j <= index_t(D[i].size()); ++j) {
                D[i][j] = enc<pk_dash>(dist(
                                           input.partition_dendrograms[i-1]
                                           , input.partition_dendrograms[j-1]));
            }
        }

        comm.send(H, D);

        for(index_t i = 1; i <= n; ++i) {
            for(index_t j = 1; j <= n; ++j) {
                //NOTE Optimize by only calculating the lower half of the symmetric matrix
                if(i <= n1 && j <= n1) {}
                else if(n1 < i && n1 < j) {}
                else {
                    dist.template calc_intercluster_distance_with_peer<pk_dash>();
                }
            }
        }
    }

    template<typename Range>
    size_t get_dimension(Range& input) const
    {
        return (*input.begin()).dimension();
    }

    template<typename Input, typename Dist = dist_t>
    table_t<plaintext>& setup(Input& input, Dist&& dist = dist_t{})
    {
        namespace phoenix = boost::phoenix;
        namespace range = boost::range;
        namespace adaptors = boost::adaptors;
        using namespace phoenix::placeholders;

        using boost::extents;
        using boost::indices;
        using phoenix::bind;
        using phoenix::cref;

        static_assert(precomputed == 3
                      , "Any other amount than 3 precomputed values is not supported");

        index_t n = in->n;

        init_Sigma(input, std::forward<Dist>(dist));

        table_t<plaintext> R_dash{extents[n][n]};
        R_dash.reindex(1);

        //Compute matrix R': R'_i,j = r'_i,j,    r'_i,j  <--$-- {0,1}^κ,    i,j ∈ [1 : n]
        std::generate(R_dash.data()
                      , R_dash.data() + R_dash.num_elements()
                      , phoenix::bind(gen_random_plain, kappa));

        table_t<paillier<pk>> R{extents[n][n]};
        table_t<paillier<pk_dash>> B_enc{extents[n][n]};

        comm.receive(R, B_enc);

        //B is an alias for Sigma for P2
        table_t<plaintext>& B = in->Sigma;
        B.reindex(1);

        //Decrypt and re-blind matrix B: B_i,j = <B_i,j> + R_i,j
        multi_for_each(B, B_enc, R_dash, arg1 = bind(dec, arg2) + arg3);

        //Blind R_i,j = R_i,j · [R'_i,j]
        multi_for_each(R, R_dash, arg1 += bind(enc<pk>, arg2));

        //Permute S, L and R via a random permutation π_2(n)
        //Send {R} to P1
        comm.send(in->pi_2.permute(R));

        in->pi_2.permute_inplace(B);

        return in->Sigma;
    }

    output_t output()
    {
        using boost::adaptors::filtered;

        using dendrogram_t = std::vector<cluster>;
        using E_t = std::vector<point<double>>;
        using J_t = std::vector<size_t>;

        constexpr size_t dendrogram_idx = 0;
        constexpr size_t E_idx = 1;
        constexpr size_t J_idx = 2;

        auto& dendrogram_ = in->dendrogram_;

        auto res = std::make_tuple(dendrogram_t{}, E_t{}, J_t{});

        for(cluster& c : dendrogram_) {
            if(c.is_cluster()) {
                std::get<dendrogram_idx>(res).emplace_back(std::move(c));
            }
        }

        size_t num_clusters = std::get<dendrogram_idx>(res).size();

        std::vector<point<paillier<pk_dash>>> E_enc(num_clusters);
        J_t& J = std::get<J_idx>(res);
        J.resize(num_clusters);
        comm.receive(E_enc, J);

        E_t& E = std::get<E_idx>(res);
        E.reserve(num_clusters);

        multi_for_each(E_enc, J, [&](auto const& e, auto const& j) {
            assert(0 < e.dimension());
            E.emplace_back(apply_on_elements([&](auto const& p) {
                return static_cast<int64_t>(dec(p))/static_cast<double>(j);
            }, e));
        });

        assert(E.size() == num_clusters);

        comm.send(E);

        return res;
    }

    //Perform initial phase for P2, i.e. connecting with P1 and
    //exchanging input sizes and public paillier_keys
    struct initial {

        initial(communicator& comm, size_t n2)
            : n1{comm.exchange(n2)}, n2{n2}, n{n1 + n2} {}

        //n1 for P1, n2 for P2, n = n1 + n2
        size_t n1, n2, n;

        table_t<plaintext> Sigma{boost::extents[n][n]};
        random_permutation pi_2{n};
        std::vector<cluster> dendrogram_;

    };

    mutable communicator comm;
    //pk for P1, pk_dash for P2
    std::unique_ptr<initial> in;
    unsigned short port;
};

key_id P2::impl_::pk;
key_id P2::impl_::pk_dash;

//======================================================================
// Common functionality for both parties
//======================================================================

template<typename Party>
void basic_party<Party>::clustering(size_t t)
{
    using clustering_t =
        hierarchical_clustering<
        aby_max_linkage<table_t<plaintext> > >;
    constexpr uint32_t nthreads = 1;
    constexpr e_mt_gen_alg mt_alg = MT_OT;
    constexpr e_role role = derived_::impl_::role;

    auto& pimpl =  static_cast<derived_&>(*this).pimpl_;

    uint32_t bitlen = kappa + 1;
    unsigned short aby_port = pimpl->port + 1;

    //Synchronize with other party to prevent timeout in ABY
    pimpl->comm.synchronize();
    pimpl->in->dendrogram_ = clustering_t{
        pimpl->in->Sigma
        , role
        , pimpl->get_aby_ip()
        , aby_port
        , get_sec_lvl(symmetric_bits)
        , bitlen
        , nthreads
        , mt_alg
    }(t);

    assert(pimpl->in->dendrogram_.size() == size_t(t)
           && "dendrogram_ should be of size t");
}

template<typename Party>
void basic_party<Party>::clustering_opt(size_t t)
{
    using clustering_t =
        hierarchical_clustering<
        opt_aby_min_linkage<table_t<plaintext> > >;
    constexpr uint32_t nthreads = 1;
    constexpr e_mt_gen_alg mt_alg = MT_OT;
    constexpr e_role role = derived_::impl_::role;

    auto& pimpl =  static_cast<derived_&>(*this).pimpl_;

    uint32_t bitlen = kappa + 1;
    unsigned short aby_port = pimpl->port + 1;

    //Synchronize with other party to prevent timeout in ABY
    pimpl->comm.synchronize();
    pimpl->in->dendrogram_ = clustering_t{
        pimpl->in->Sigma
        , role
        , pimpl->get_aby_ip()
        , aby_port
        , get_sec_lvl(symmetric_bits)
        , bitlen
        , nthreads
        , mt_alg
    }(t);

    assert(pimpl->in->dendrogram_.size() == size_t(t)
           && "dendrogram_ should be of size t");
}

#define MPO19_STORE_SIGMA \
    auto& in = static_cast<derived_ &>(*this).pimpl_->in; \
    auto& Sigma = static_cast<derived_ &>(*this).pimpl_->in->Sigma; \
    auto& comm = static_cast<derived_ const&>(*this).pimpl_->comm; \
    constexpr e_role role = derived_::impl_::role; \
    \
    table_t<size_t>& sig = pca_internal_state::Sigma; \
    sig.resize(boost::extents[in->n][in->n]); \
    if(ROLE_P1 == role) { \
        comm.send(Sigma); \
    } else if(ROLE_P2 == role) { \
        table_t<plaintext> blinds{boost::extents[in->n][in->n]}; \
        comm.receive(blinds); \
        \
        assert(Sigma.num_elements() == blinds.num_elements() \
               && "Sigma and blinds must have the same number of elements"); \
        \
        assert(Sigma.shape()[0] ==blinds.shape()[0] \
               && "Sigma and blinds must have the same shape"); \
        assert(Sigma.shape()[1] == blinds.shape()[1] \
               && "Sigma and blinds must have the same shape"); \
        \
        multi_for_each(sig, Sigma, blinds \
        , [](size_t& r, plaintext const& s, plaintext const& b) { \
            assert( (s.valid() && b.valid()) \
                    && "s and b must be valid"); \
            r = static_cast<size_t>(s + -b); \
        }); \
    }

template<typename Party>
void basic_party<Party>::print_sigma() const
{
    auto const& in = static_cast<derived_ const&>(*this).pimpl_->in;
    auto const& Sigma = static_cast<derived_ const&>(*this).pimpl_->in->Sigma;
    auto& comm = static_cast<derived_ const&>(*this).pimpl_->comm;
    e_role role = derived_::impl_::role;

    if(ROLE_P1 == role) {
        comm.send(Sigma);
    } else if(ROLE_P2 == role) {
        table_t<plaintext> blinds{boost::extents[in->n][in->n]};
        comm.receive(blinds);

        assert(Sigma.num_elements() == blinds.num_elements()
               && "Sigma and blinds must have the same number of elements");


        assert(Sigma.shape()[0] ==blinds.shape()[0]
               && "Sigma and blinds must have the same shape");
        assert(Sigma.shape()[1] == blinds.shape()[1]
               && "Sigma and blinds must have the same shape");

        std::cout << "Sigma:\n";
        multi_for_each(Sigma, blinds
        , [&, counter=0u](plaintext const& s, plaintext const& b) mutable{
            assert( ((s.valid() && b.valid())
                     || (!s.valid() && !b.valid()))
                    && "s and b are valid or they are both not valid");

            std::cout << std::setw(6);

            if(s.valid())
            {
                uint64_t p = s + -b;
                std::cout << p;
            } else
            {
                std::cout << '/';
            }

            if(++counter == Sigma.shape()[0])
            {
                std::cout << "\n";
                counter = 0;
            } else
            {
                std::cout << ", ";
            }
        });
    }
}

template<typename Party>
typename basic_party<Party>::output_t
basic_party<Party>::execute_pca(array_t<point<int64_t>>& input, size_t t)
{
    typename derived_::impl_& p = *static_cast<derived_&>(*this).pimpl_;
    p.in = std::make_unique<typename derived_::impl_::initial>(p.comm, input.size());
    assert(nullptr != p.in);
    p.setup(input);
#ifdef MPO19_TESTING_PCA_
    MPO19_STORE_SIGMA
#endif //MPO19_TESTING_PCA_

    clustering(t);
    return p.output();
}

template<typename Party>
typename basic_party<Party>::output_t
basic_party<Party>::execute_opt(array_t<point<int64_t>>& input, size_t t)
{
    typename derived_::impl_& p = *static_cast<derived_&>(*this).pimpl_;
    p.in = std::make_unique<typename derived_::impl_::initial>(p.comm, input.size());
    assert(nullptr != p.in);

    p.setup(input);

#ifdef MPO19_TESTING_PCA_
    MPO19_STORE_SIGMA
#endif //MPO19_TESTING_PCA_

    clustering_opt(t);
    return p.output();
}

template<typename Party>
void basic_party<Party>::insert_other_parties_centroids(std::vector<point<double>>& centroids)
{
    auto& comm = static_cast<derived_&>(*this).pimpl_->comm;

    comm.send(centroids);

    std::vector<point<double>> other_centroids;
    comm.receive(other_centroids);

    centroids.reserve(centroids.size() + other_centroids.size());
    centroids.insert(centroids.end(), other_centroids.begin(), other_centroids.end());
}

template<typename Party>
std::vector<cluster>
basic_party<Party>::pcure0(std::vector<point<int64_t>> const& input
                           , size_t s, size_t p, size_t q, size_t t1, size_t t2, size_t t)
{
    constexpr e_role role = derived_::impl_::role;

    typename derived_::impl_& impl = *static_cast<derived_&>(*this).pimpl_;
    impl.in = std::make_unique<typename derived_::impl_::initial>(impl.comm, input.size());
    assert(nullptr != impl.in);

    pcure<clear_clustering_t<max_linkage>> pcure{s, p, q, t1, t2, t/2 + (role == ROLE_P1 ? t%2 : 0)};
    auto centroids = pcure.get_cluster_centroids(input);
    insert_other_parties_centroids(centroids);
    return pcure.assign_clusters(input, centroids);
}

template<typename Party>
std::vector<cluster>
basic_party<Party>::pcure0_opt(std::vector<point<int64_t>> const& input
                               , size_t s, size_t p, size_t q, size_t t1, size_t t2, size_t t)
{
    constexpr e_role role = derived_::impl_::role;

    typename derived_::impl_& impl = *static_cast<derived_&>(*this).pimpl_;
    impl.in = std::make_unique<typename derived_::impl_::initial>(impl.comm, input.size());
    assert(nullptr != impl.in);

    pcure<clear_clustering_t<max_linkage>> pcure{s, p, q, t1, t2, t/2 + (role == ROLE_P1 ? t%2 : 0)};
    auto centroids = pcure.get_cluster_centroids(input);
    insert_other_parties_centroids(centroids);
    return pcure.assign_clusters(input, centroids);
}

template<typename Party>
std::vector<point<double>> basic_party<Party>::remove_outliers(output_t& out, pcure_clusters<int64_t>& clusters, size_t t2)
{
    using boost::accumulate;

    typename derived_::impl_& impl = *static_cast<derived_&>(*this).pimpl_;
    constexpr e_role role = derived_::impl_::role;

    size_t offset = (role == ROLE_P1 ? 0 : impl.in->n1);
    std::vector<size_t> cluster_sizes;
    cluster_sizes.reserve(std::get<0>(out).size());
    for(cluster const& c : std::get<0>(out)) {
        size_t real_size = accumulate(get_indices(c), 0, [&](size_t sum, index_t i) {
            if( (role == ROLE_P1 && i >= index_t(clusters.partition_dendrograms.size())) ||
                (role == ROLE_P2 && i < index_t(offset))) {
                return sum;
            }
            i = i - offset;
            assert(0 <= i && i < index_t(clusters.partition_dendrograms.size()));
            return sum + get_indices(clusters.partition_dendrograms[i]).size();
        });
        cluster_sizes.emplace_back(real_size);
    }
    impl.comm.send(cluster_sizes);
    std::vector<size_t> cluster_sizes_peer;
    impl.comm.receive(cluster_sizes_peer);
    assert(cluster_sizes.size() == cluster_sizes_peer.size());
    for(size_t i = 0; i != cluster_sizes.size(); ++i) {
        cluster_sizes[i] += cluster_sizes_peer[i];
    }

    auto& centroids = std::get<1>(out);

    centroids.erase(std::remove_if(centroids.begin(), centroids.end(), [&, i = size_t(0)](point<double> const& p) mutable {
        assert(i < std::get<0>(out).size());
        return p.dimension() == 0 || cluster_sizes[i++] < t2;
    })
    , centroids.end());

    return centroids;
}

template<typename Party>
std::vector<point<double>> basic_party<Party>::pca_on_A_clusters(pcure_clusters<int64_t>& clusters, size_t t2, size_t t)
{
    assert(!clusters.partition_dendrograms.empty());

    typename derived_::impl_& impl = *static_cast<derived_&>(*this).pimpl_;
    impl.in = std::make_unique<typename derived_::impl_::initial>(impl.comm, clusters.partition_dendrograms.size());
    assert(nullptr != impl.in);

    constexpr uint32_t nthreads = 1;
    constexpr e_mt_gen_alg mt_alg = MT_OT;
    constexpr e_role role = derived_::impl_::role;
    uint32_t bitlen = kappa + 1;
    unsigned short aby_port = impl.port + 1;

    cluster_dist_t<int64_t, max_linkage_policy> clear_dist{clusters.values};

    {
        using dist_t = intercluster_dist_t<int64_t, basic_aby_linkage_policy<aby_max_linkage_policy>,  decltype(clear_dist)>;

        dist_t dist{
            clusters.values, impl.comm, kappa, clear_dist
            , role, impl.get_aby_ip(), aby_port
            , get_sec_lvl(symmetric_bits), bitlen, nthreads, mt_alg
        };

        impl.setup(clusters, dist);
    } //destroy dist, to release aby_port

    clustering(t);
    auto out = impl.output();

    return remove_outliers(out, clusters, t2);
}

template<typename Party>
std::vector<point<double>> basic_party<Party>::opt_on_A_clusters(pcure_clusters<int64_t>& clusters, size_t t2, size_t t)
{
    assert(!clusters.partition_dendrograms.empty());

    typename derived_::impl_& impl = *static_cast<derived_&>(*this).pimpl_;
    impl.in = std::make_unique<typename derived_::impl_::initial>(impl.comm, clusters.partition_dendrograms.size());
    assert(nullptr != impl.in);

    constexpr uint32_t nthreads = 1;
    constexpr e_mt_gen_alg mt_alg = MT_OT;
    constexpr e_role role = derived_::impl_::role;
    uint32_t bitlen = kappa + 1;
    unsigned short aby_port = impl.port + 1;

    cluster_dist_t<int64_t, min_linkage_policy> clear_dist{clusters.values};

    {
        using dist_t = intercluster_dist_t<int64_t, basic_aby_linkage_policy<aby_min_linkage_policy>,  decltype(clear_dist)>;

        dist_t dist{
            clusters.values, impl.comm, kappa, clear_dist
            , role, impl.get_aby_ip(), aby_port
            , get_sec_lvl(symmetric_bits), bitlen, nthreads, mt_alg
        };

        impl.setup(clusters, dist);
    } //destroy dist, to release aby_port

    clustering_opt(t);
    auto out = impl.output();

    return remove_outliers(out, clusters, t2);
}

template<typename Party>
std::vector<cluster>
basic_party<Party>::pcure1(std::vector<point<int64_t>> const& input
                           , size_t s, size_t p, size_t q, size_t t1, size_t t2, size_t t)
{

    pcure<clear_clustering_t<max_linkage>> pcure{s, p, q, t1, t2, t};
    std::vector<point<int64_t>> in{input};
    pcure_clusters<int64_t> clusters = pcure.A_clustering(in);

    return assign_clusters(input, pca_on_A_clusters(clusters, t2, t));
}

template<typename Party>
std::vector<cluster>
basic_party<Party>::pcure1_opt(std::vector<point<int64_t>> const& input
                               , size_t s, size_t p, size_t q, size_t t1, size_t t2, size_t t)
{

    pcure<clear_clustering_t<max_linkage>> pcure{s, p, q, t1, t2, t};
    std::vector<point<int64_t>> in{input};
    pcure_clusters<int64_t> clusters = pcure.A_clustering(in);

    return assign_clusters(input, opt_on_A_clusters(clusters, t2, t));
}

template<typename Party>
std::vector<cluster>
basic_party<Party>::pcure2(array_t<point<int64_t>>& input
                           , size_t s, size_t q, size_t t1, size_t t2, size_t t)
{
    assert(input.size() > s &&
           "Cannot take a sample of size s with s being bigger than the input size");
    assert(s/q > 0 &&
           "The target clusters shall be grater than zero, i.e. s/q > 0");

    size_t const target_clusters = s/q;

    std::vector<point<int64_t>> in{input.begin(), input.end()};

    auto sample = random_sample(in, s);
    array_t<point<int64_t>> sample_arr{boost::extents[sample.size()]};
    boost::copy(sample, sample_arr.begin());

    typename derived_::impl_& impl = *static_cast<derived_&>(*this).pimpl_;
    impl.in = std::make_unique<typename derived_::impl_::initial>(impl.comm, sample_arr.size());
    assert(nullptr != impl.in);

    impl.setup(sample_arr);

    using clustering_t =
        hierarchical_clustering<
        aby_max_linkage<table_t<plaintext> > >;

    constexpr uint32_t nthreads = 1;
    constexpr e_mt_gen_alg mt_alg = MT_OT;
    constexpr e_role role = derived_::impl_::role;

    uint32_t bitlen = kappa + 1;
    unsigned short aby_port = impl.port + 1;

    //Synchronize with other party to prevent timeout in ABY
    impl.comm.synchronize();
    clustering_t clustering{
        impl.in->Sigma
        , role
        , impl.get_aby_ip()
        , aby_port
        , get_sec_lvl(symmetric_bits)
        , bitlen
        , nthreads
        , mt_alg
    };

    impl.in->dendrogram_ = clustering.partial_evaluation(target_clusters);

    //Invalidate clusters of size < t1
    for(index_t i = 0; i != index_t(impl.in->dendrogram_.size()); ++i) {
        if(impl.in->dendrogram_[i].is_cluster() && get_indices(impl.in->dendrogram_[i]).size() < t1) {
            impl.in->dendrogram_[i] = cluster{};
            //Invalidate cluster distances associated with the invalidated cluster

            assert(impl.in->Sigma.index_bases()[0] == impl.in->Sigma.index_bases()[1]
                   && "Index base for both dimensions of Sigma must be equal");
            index_t base = impl.in->Sigma.index_bases()[0];
            for(index_t j = 0; j != index_t(impl.in->Sigma.size()); ++j) {
                if(i != j) {
                    clustering.invalidate(impl.in->Sigma[i + base][j + base]);
                    clustering.invalidate(impl.in->Sigma[j + base][i + base]);
                }
            }
        }
    }

    impl.in->dendrogram_ = clustering.partial_evaluation(t, std::move(impl.in->dendrogram_));

    auto out = impl.output();

    auto& dendrogram = std::get<0>(out);
    auto& centroids = std::get<1>(out);
    auto& cluster_sizes = std::get<2>(out);

    assert(dendrogram.size() == cluster_sizes.size());
    assert(dendrogram.size() == centroids.size());
    for(size_t i = 0; i != dendrogram.size(); ++i) {
        if(cluster_sizes[i] < t2) {
            centroids.erase(centroids.begin() + i);
        }
    }

    return assign_clusters(input, centroids);
}

template<typename Party>
std::vector<cluster>
basic_party<Party>::pcure2_opt(array_t<point<int64_t>>& input
                               , size_t s, size_t q, size_t t1, size_t t2, size_t t)
{
    assert(input.size() > s &&
           "Cannot take a sample of size s with s being bigger than the input size");
    assert(s/q > 0 &&
           "The target clusters shall be grater than zero, i.e. s/q > 0");

    size_t const target_clusters = s/q;

    std::vector<point<int64_t>> in{input.begin(), input.end()};

    auto sample = random_sample(in, s);
    array_t<point<int64_t>> sample_arr{boost::extents[sample.size()]};
    boost::copy(sample, sample_arr.begin());

    typename derived_::impl_& impl = *static_cast<derived_&>(*this).pimpl_;
    impl.in = std::make_unique<typename derived_::impl_::initial>(impl.comm, sample_arr.size());
    assert(nullptr != impl.in);

    impl.setup(sample_arr);

    using clustering_t =
        hierarchical_clustering<
        opt_aby_min_linkage<table_t<plaintext> > >;

    constexpr uint32_t nthreads = 1;
    constexpr e_mt_gen_alg mt_alg = MT_OT;
    constexpr e_role role = derived_::impl_::role;

    uint32_t bitlen = kappa + 1;
    unsigned short aby_port = impl.port + 1;

    //Synchronize with other party to prevent timeout in ABY
    impl.comm.synchronize();
    clustering_t clustering{
        impl.in->Sigma
        , role
        , impl.get_aby_ip()
        , aby_port
        , get_sec_lvl(symmetric_bits)
        , bitlen
        , nthreads
        , mt_alg
    };

    impl.in->dendrogram_ = clustering.partial_evaluation(target_clusters);

    //Invalidate clusters of size < t1
    for(index_t i = 0; i != index_t(impl.in->dendrogram_.size()); ++i) {
        if(impl.in->dendrogram_[i].is_cluster() && get_indices(impl.in->dendrogram_[i]).size() < t1) {
            impl.in->dendrogram_[i] = cluster{};
            //Invalidate cluster distances associated with the invalidated cluster
            for(index_t j = 0; j != index_t(impl.in->Sigma.size()); ++j) {

                if(i != j) {
                    clustering.invalidate(impl.in->Sigma[i][j]);
                    clustering.invalidate(impl.in->Sigma[j][i]);
                }
            }
            clustering.invalidate(clustering.precomputed_merge_values[i]);
        }
    }

    impl.in->dendrogram_ = clustering.partial_evaluation(t, std::move(impl.in->dendrogram_));

    auto out = impl.output();

    auto& dendrogram = std::get<0>(out);
    auto& centroids = std::get<1>(out);
    auto& cluster_sizes = std::get<2>(out);

    assert(dendrogram.size() == cluster_sizes.size());
    assert(dendrogram.size() == centroids.size());
    for(size_t i = 0; i != dendrogram.size(); ++i) {
        if(cluster_sizes[i] < t2) {
            centroids.erase(centroids.begin() + i);
        }
    }

    return assign_clusters(input, centroids);
}

P1::P1(unsigned short port)
    : pimpl_{std::make_unique<impl_>(port)} {}

P1::~P1() = default;

P2::P2(std::string const& destination, unsigned short port)
    : pimpl_{std::make_unique<impl_>(destination, port)} {}

P2::~P2() = default;

void P1::print_L()
{
    pimpl_->comm.send(pimpl_->in->L);
}

void P2::print_L()
{
    std::vector<point<paillier<impl_::pk_dash>>> L(pimpl_->in->n);
    pimpl_->comm.receive(L);
    for(auto const& pnt : L) {
        std::cout << apply_on_elements([](auto&& p) {
            return static_cast<int64_t>(dec(p));
        }, pnt) << " ";
    }
    std::cout << std::endl;
}

template class basic_party<P1>;
template class basic_party<P2>;
