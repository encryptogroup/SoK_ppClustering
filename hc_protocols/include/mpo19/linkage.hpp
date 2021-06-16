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

#ifndef MPO19_LINKAGE_30062020
#define MPO19_LINKAGE_30062020

#include <mpo19/typedefs.hpp>
#include <mpo19/utility.hpp>
#include <mpo19/paillier_wrapper.hpp>
#include <mpo19/cluster.hpp>

#include <abycore/aby/abyparty.h>
#include <abycore/circuit/booleancircuits.h>
#include <abycore/sharing/sharing.h>
#include <abycore/ABY_utils/ABYconstants.h>
#include <ENCRYPTO_utils/crypto/crypto.h>

#include <boost/hana.hpp>

#include <utility>
#include <cstddef>
#include <numeric>
#include <exception>
#include <sstream>

namespace mpo19
{

uint64_t aby_bytes_sent();
uint64_t aby_bytes_received();
void aby_reset_counters();
void aby_update_counters();

//template<T, U> struct linkage_precomputed
//
//Provided specialization for T = table_t<T'>, U = LinkagePolicy
//See specialization below for more details
template<typename, typename>
struct linkage_precomputed;


/* template struct linkage_precomputed<table_t<T>, LinkagePolicy>
 *
 * Uses a precomputed distance table to compute the indices of the clusters to merge
 * Forwards certain details of the computation of the merge indices to LinkagePolicy.
 *
 * Requirements on LinkagePolicy:
 *
 *     Contructors:
 *         Any
 *             additional constructor arguments than those needed by linkage_precomputed
 *             are forwarded to LinkagePolicy's constructor
 *
 *     Methods:
 *         let lifted_t be any type fulfilling the requirements describd ine method lift:
 *
 *         lifted_t lift(T const& t, index_t i, index_t j)
 *              lifts the passed entry t in the table at position i,j to a lifted_t.
 *              lifted_t must store enough information so that two values
 *              of lifted_t can be compared and the corresponding indices can
 *              be retrieved afterwards.
 *
 *         lifted_t select(lifted_t lhs, lifted_t rhs)
 *              selects a lifted_t from two input lifted_t by comparing the stored
 *              table entry value and selecting one according to some comparison.
 *              The two lifted arguments are moved into this function, therefore it
 *              is also valid to provide overloads for lifted_t&& or lifted const& instead
 *
 *         std::pair<index_t, index_t> retrieve_indices(lifted_t l)
 *              returns the two indices of the corresponding table entry stored in l
 *
 *         T cluster_distance(T&, T&)
 *              returns the cluster distance between two clusters.
 *              linkage_precomputed updates the distance table such that the cluster
 *              distance between two clusters can be computed by simply using the
 *              two input entries, if the cluster_distance operation is associative
 *
 *         void invalidate(T& entry)
 *              Invalidates the input entry e.g. by setting it to a special value.
 *              T counts as invalidad iff for every lifted_t li and every lifted_t lv
 *              with lv not being invalidated : select(li, lv) == lv
 *
 *         bool is_valid(T const& t)
 *              checks if t denotes a valid entry in the table. If t was invalidated before
 *              then is_valid(t) == true, else is_valid(t) == false
 *
 * This struct assumes that the cluster_distance method is associative to yield
 * correct results, which is the case for e.g. min and max linkage
 */
template<typename T, typename LinkagePolicy>
struct linkage_precomputed<table_t<T>, LinkagePolicy> : public LinkagePolicy {

    using LinkagePolicy::lift;
    using LinkagePolicy::select;
    using LinkagePolicy::retrieve_indices;
    using LinkagePolicy::cluster_distance;
    using LinkagePolicy::is_valid;
    using LinkagePolicy::invalidate;


    template<typename... Args>
    linkage_precomputed(table_t<T>& table, Args&&... args)
        : LinkagePolicy{std::forward<Args>(args)...}, table_{table}
    {
        assert(table_.shape()[0] == table_.shape()[1]
               && "table_ should be a square matrix");
        assert(table_.index_bases()[0] == table_.index_bases()[1]
               && "index base of each dimesion should be equal");
    }

    std::pair<index_t, index_t> get_merge_indices();

    size_t get_input_size() const
    {
        return table_.size();
    }

    table_t<T>& table_;
};


template<typename, typename>
struct linkage_opt;
/*
 * This struct assumes that the select operation is associative with the
 * cluster_distance operation and the cluster_distance operation is associative
 * to yield correct results.
 * More precisely all of the following statements must be true:
 * Notation: a <cd> b means cluster_distance(a, b)
 *                   a <sel> b means select(a, b)
 *
 *     1. (a <cd> b) <cd> c = a <cd> (b <cd> c)
 *     2. (a <sel> b) <cd> c = a <sel> (b <cd> c)
 *
 */
template<typename T, typename LinkagePolicy>
struct linkage_opt<table_t<T>, LinkagePolicy> : LinkagePolicy {

    using LinkagePolicy::lift;
    using LinkagePolicy::select;
    using LinkagePolicy::retrieve_indices;
    using LinkagePolicy::cluster_distance;
    using LinkagePolicy::is_valid;
    using LinkagePolicy::invalidate;

    template<typename... Args>
    linkage_opt(table_t<T>& table, Args&&... args);

    std::pair<index_t, index_t> get_merge_indices();

    size_t get_input_size() const
    {
        return table_.size();
    }

    table_t<T>& table_;
    std::vector<index_t> precomputed_merge_indices{};
    std::vector<T> precomputed_merge_values{};
};

template<typename LinkagePolicy>
struct basic_linkage_policy : LinkagePolicy {

    template<typename T>
    auto lift(T const& val, index_t i, index_t j)
    {
        return boost::hana::make_tuple(val, i, j);
    }

    template<typename Lifted>
    Lifted select(Lifted&& lhs, Lifted&& rhs)
    {
        using namespace boost::hana::literals;
        if(lhs[0_c] <= rhs[0_c]) {
            return std::move(lhs);
        } else {
            return std::move(rhs);
        }
    }

    template<typename Lifted>
    std::pair<index_t, index_t> retrieve_indices(Lifted const& l)
    {
        using namespace boost::hana::literals;
        return std::make_pair(l[1_c], l[2_c]);
    }

    using LinkagePolicy::cluster_distance;

    template<typename T>
    bool is_valid(T const& entry)
    {
        return std::numeric_limits<T>::max() != entry;
    }

    template<typename T>
    void invalidate(T& entry)
    {
        //A max value will never be selected by select,
        //which tries to minimize the distance
        entry = std::numeric_limits<T>::max();
    }

};

struct max_linkage_policy {
    template<typename T>
    T cluster_distance(T const& i_dist, T const& j_dist)
    {
        return std::max(i_dist, j_dist);
    }
};

struct min_linkage_policy {
    template<typename T>
    T cluster_distance(T const& i_dist, T const& j_dist)
    {
        return std::min(i_dist, j_dist);
    }
};

template<typename ABYLinkagePolicy>
struct basic_aby_linkage_policy :  ABYLinkagePolicy {

    using lifted_t = boost::hana::tuple<share*, share*, share*>;

    basic_aby_linkage_policy(e_role role
                             , const std::string& addr = "127.0.0.1"
                             , uint16_t port = 7766
                             , seclvl seclvl = LT
                             , uint32_t bitlen = 32
                             , uint32_t nthreads = 2
                             , e_mt_gen_alg mg_algo = MT_OT
                             , uint32_t reservegates = 65536
                             , const std::string& abycircdir = ABY_CIRCUIT_DIR);

    lifted_t lift(plaintext const& p, index_t i, index_t j);

    lifted_t select(lifted_t&& lhs, lifted_t && rhs);

    plaintext cluster_distance(plaintext& lhs, plaintext& rhs);

    std::pair<index_t, index_t> retrieve_indices(lifted_t const& l);

    bool is_valid(plaintext const& entry)
    {
        return entry.valid();
    }

    void invalidate(plaintext& entry)
    {
        entry = plaintext{};
    }

    share* put_idx_(index_t idx);

    //Put plaintext into ABY Circuit and unblind it using ABY.
    share* put_unblinded_input_(plaintext const& input);

    e_role role_;
    uint32_t bitlen_;
    ABYParty party_;
    Circuit* circ_;
};

struct aby_max_linkage_policy {
    share* cluster_distance(Circuit* circ, share* left_cluster_distance, share* right_cluster_distance)
    {
        share* gt = circ->PutGTGate(left_cluster_distance, right_cluster_distance);
        share* max = circ->PutMUXGate(left_cluster_distance, right_cluster_distance, gt);
        return max;
    }
};

struct aby_min_linkage_policy {
    share* cluster_distance(Circuit* circ, share* left_cluster_distance, share* right_cluster_distance)
    {
        share* gt = circ->PutGTGate(left_cluster_distance, right_cluster_distance);
        share* min = circ->PutMUXGate(right_cluster_distance, left_cluster_distance, gt);
        return min;
    }
};


template<typename T>
using max_linkage = linkage_precomputed<T, basic_linkage_policy<max_linkage_policy>>;
template<typename T>
using opt_min_linkage = linkage_opt<T, basic_linkage_policy<min_linkage_policy>>;
template<typename T>
using aby_max_linkage = linkage_precomputed<T, basic_aby_linkage_policy<aby_max_linkage_policy>>;
template<typename T>
using opt_aby_min_linkage = linkage_opt<T, basic_aby_linkage_policy<aby_min_linkage_policy>>;

//==============================================================================
//Implementation
//==============================================================================

struct aby_error : std::runtime_error {
    using std::runtime_error::runtime_error;
};

template<typename T, typename LinkagePolicy>
std::pair<index_t, index_t> linkage_precomputed<table_t<T>, LinkagePolicy>::get_merge_indices()
{
    index_t base = table_.index_bases()[0];
    index_t i = 0, j = 1;
    size_t n = get_input_size();

    assert(i != j &&
           "i should never be equal to j, as we don't want to"
           " assign an element on the diagonal into val");
    auto val = lift(table_[i + base][j + base], i, j);
    ++j;
    for(; i != index_t(n); ++i, j = 0) {
        for(; j != index_t(n); ++j) {
            //Skip if i == j, to skip elements on the diagonal
            if(i == j) continue;
            assert(!(i == 0 && j == 1) &&
                   "We already checked (0, 1) so we don't want to check it again");
            val = select(std::move(val), lift(table_[i + base][j + base], i, j));
        }
    }
    auto p = retrieve_indices(std::move(val));

    if(!(0 <= p.first && p.first < index_t(n) &&
         0 <= p.second && p.second < index_t(n))) {
        std::stringstream ss;
        ss << "ABY returned garbage indices: (" << p.first << ", " << p.second << ")";
        throw aby_error{ss.str()};
    }

    assert(p.first != p.second);

    if(p.first > p.second)
        std::swap(p.first, p.second);

    i = p.first;
    j = p.second;

    for(index_t k = 0; k != index_t(n); ++k) {

        //Update linkages for newly merged cluster
        if(k != i && k != j
           && is_valid(table_[i + base][k + base])
           && is_valid(table_[j + base][k + base])) {

            T c_dist = cluster_distance(
                           table_[i + base][k + base]
                           , table_[j + base][k + base]);

            table_[i + base][k + base] = c_dist;
            table_[k + base][i + base] = std::move(c_dist);
        }

        //Cluster j only refers to i after merge,
        //so we invalidate entries related to j
        if(k != j) {
            invalidate(table_[k + base][j + base]);
            invalidate(table_[j + base][k + base]);
        }
    }

    return p;
}

template<typename T, typename LinkagePolicy>
template<typename... Args>
linkage_opt<table_t<T>, LinkagePolicy>::linkage_opt(table_t<T>& table, Args&&... args)
    : LinkagePolicy{std::forward<Args>(args)...}, table_{table}
{
    precomputed_merge_indices.reserve(get_input_size());
    precomputed_merge_values.reserve(get_input_size());

    index_t base = table_.index_bases()[0];
    size_t n = get_input_size();

    for(index_t i = 0, j = 1; i != index_t(n); ++i) {
        assert(i != j &&
               "i should never be equal to j, as we don't want to"
               " assign an element on the diagonal to val");

        auto val = lift(table_[i + base][j + base], i, j);
        ++j;

        for(; j != index_t(n); ++j) {
            //Skip if i == j, to skip elements on the diagonal
            if(i != j) {
                assert(!( (i == 0 && j == 1) || (i != 0 && j == 0) ) &&
                       "val is set to the element in (0, 1) in the first"
                       "iteration and to (i, 0) in subsequent iterations, "
                       "so we shouldn't try to select val a second time "
                       "in this loop");
                val = select(std::move(val), lift(table_[i + base][j + base], i, j));
            }
        }
        auto p = retrieve_indices(std::move(val));

        if(!(0 <= p.first && p.first < index_t(n) &&
             0 <= p.second && p.second < index_t(n))) {
            std::stringstream ss;
            ss << "ABY returned garbage indices: (" << p.first << ", " << p.second << ")";
            throw aby_error{ss.str()};
        }

        assert(p.first == i &&
               "We only input the ith row, so the selected row index "
               "should always be i");
        assert(p.first != p.second &&
               "A diagonal element should never be chosen");

        precomputed_merge_values.emplace_back(table_[p.first + base][p.second + base]);
        //We insert the min of each row i into precomputed_merge_indices
        precomputed_merge_indices.emplace_back(p.second);
        j = 0;
    }
    assert(precomputed_merge_values.size() == n);
    assert(precomputed_merge_indices.size() == n);
}


template<typename T, typename LinkagePolicy>
std::pair<index_t, index_t>
linkage_opt<table_t<T>, LinkagePolicy>::get_merge_indices()
{
    index_t base = table_.index_bases()[0];
    size_t n = get_input_size();

    index_t i = 0, j = precomputed_merge_indices[i];
    auto val = lift(precomputed_merge_values[i], i, j);
    for(i = 1; i != index_t(n); ++i) {
        if(std::numeric_limits<index_t>::max() != j) {
            j = precomputed_merge_indices[i];
            val = select(std::move(val), lift(precomputed_merge_values[i], i, j));
        }
    }

    auto p = retrieve_indices(std::move(val));
    if(p.first > p.second)
        std::swap(p.first, p.second);

    i = p.first;
    j = p.second;

    //As i and j got merged, we invalidate the value in i, as
    //row j will effectively stop to exist after updating, so the
    //value stored in precomputed_merge_values[i] is
    //invalid and should be replaced by another value
    invalidate(precomputed_merge_values[i]);

    //Update
    for(index_t k = 0; k != index_t(n); ++k) {
        if(k != i && k != j
           && is_valid(table_[i + base][k + base])
           && is_valid(table_[j + base][k + base])) {

            T c_dist = cluster_distance(
                           table_[i + base][k + base]
                           , table_[j + base][k + base]);

            auto& R_i = precomputed_merge_values[i];
            auto sel = select(lift(c_dist, 1, 1), lift(R_i, 2, 2));

            auto indices = retrieve_indices(std::move(sel));
            if(indices.first == 1) {
                R_i = c_dist;
                precomputed_merge_indices[i] = k;
            }

            table_[i + base][k + base] = c_dist;
            table_[k + base][i + base] = c_dist;

            auto& R_k = precomputed_merge_values[k];
            sel = select(lift(c_dist, 1, 1), lift(R_k, 2, 2));

            indices = retrieve_indices(std::move(sel));
            if(indices.first == 1) {
                R_k = std::move(c_dist);
                precomputed_merge_indices[k] = i;
            }

        }

        if(k != j) {
            invalidate(table_[j + base][k + base]);
            invalidate(table_[k + base][j + base]);
        }
    }
    invalidate(precomputed_merge_values[j]);

    return p;
}

template<typename ABYLinkagePolicy>
basic_aby_linkage_policy<ABYLinkagePolicy>::basic_aby_linkage_policy(e_role role
        , const std::string& addr
        , uint16_t port
        , seclvl seclvl
        , uint32_t bitlen
        , uint32_t nthreads
        , e_mt_gen_alg mg_algo
        , uint32_t reservegates
        , const std::string& abycircdir)
    : role_{role}, bitlen_{bitlen}
    , party_{role, addr, port, seclvl, bitlen, nthreads, mg_algo, reservegates, abycircdir}
    , circ_{party_.GetSharings()[S_YAO]->GetCircuitBuildRoutine()} {}

template<typename ABYLinkagePolicy>
typename basic_aby_linkage_policy<ABYLinkagePolicy>::lifted_t
basic_aby_linkage_policy<ABYLinkagePolicy>::lift(plaintext const& p, index_t i, index_t j)
{
    namespace hana = boost::hana;

    //We skip invalid values
    if(!p.valid()) {
        return {nullptr, nullptr,  nullptr};
    }

    share* s_p = put_unblinded_input_(p);
    share* s_i = put_idx_(i);
    share* s_j = put_idx_(j);

    return {s_p, s_i, s_j};
}

template<typename ABYLinkagePolicy>
typename basic_aby_linkage_policy<ABYLinkagePolicy>::lifted_t
basic_aby_linkage_policy<ABYLinkagePolicy>::select(lifted_t&& lhs, lifted_t&& rhs)
{
    namespace hana = boost::hana;
    using namespace hana::literals;

    //either lhs or rhs are invalid, so we skip them
    //by returning the other tuple as winner
    if(lhs[0_c] == nullptr) {
        return std::move(rhs);
    }
    if(rhs[0_c] == nullptr) {
        return std::move(lhs);
    }

    share* gt = circ_->PutGTGate(lhs[0_c], rhs[0_c]);

    lifted_t res = hana::transform(hana::zip(rhs, lhs)
    , [&, this](hana::tuple<share*, share*> zipped) {
        assert(hana::contains(rhs, zipped[0_c])
               && "zipped[0_c] should contain a share taken from rhs");
        assert(hana::contains(lhs, zipped[1_c])
               && "zipped[1_c] should contain a share taken from lhs");
        //If gt is 1, then lhs[0_c] > rhs[0_c], so we select rhs.
        //If gt is 0, then lhs[0_c] <= rhs[0_c] ,  so we select lhs.
        return circ_->PutMUXGate(zipped[0_c], zipped[1_c], gt);
    });
    return res;
}

template<typename ABYLinkagePolicy>
plaintext basic_aby_linkage_policy<ABYLinkagePolicy>::cluster_distance(plaintext& lhs, plaintext& rhs)
{

    if(!lhs.valid()){
        return std::move(rhs);
    }
    BooleanCircuit* circ = static_cast<BooleanCircuit*>(
                               party_.GetSharings()[S_YAO]->GetCircuitBuildRoutine());

    //B_i,k − R_i,k
    share* unblinded_lhs = put_unblinded_input_(lhs);
    //B_j,k − R_j,k
    share* unblinded_rhs = put_unblinded_input_(rhs);

    share* linkage = ABYLinkagePolicy::cluster_distance(circ_, unblinded_lhs, unblinded_rhs);

    plaintext X;
    //Put X in Circuit
    share* blind = nullptr;
    if(ROLE_P1 == role_) {
        //bitlen_ is equals kappa + 1
        X = gen_random_plain(bitlen_ - 1);
        std::vector<uint_fast32_t> X_vec = X;
        assert(!X_vec.empty() && "X_vec must not be empty");
        blind = circ->PutINGate(X_vec.data(), bitlen_, ROLE_P1);
    } else {
        blind = circ->PutDummyINGate(bitlen_);
    }
    assert(blind != nullptr && "blind must have been initialized here");

    //max{B_i,k − R_i,k , B_j,k − R_j,k} + X
    share* reblind = circ->PutADDGate(linkage, blind);

    share* s_out = circ->PutOUTGate(reblind, ROLE_P2);

    party_.ExecCircuit();
    aby_update_counters();

    uint8_t* out = nullptr;
    if(ROLE_P2 == role_) {
        assert(!X.valid() && "X must not be valid for P2");
        //P2 gets output Y (which is stored in X)
        out = s_out ->get_clear_value_ptr();
        X = plaintext{static_cast<void const*>(out)
                      , bytes_to_hold_bitlen(bitlen_)};
    }

    party_.Reset();
    //The result (⊥; Y) is translated to (X; Y) where P1
    //gets X from MaxDist instead of generating it himself
    return X;
}

template<typename ABYLinkagePolicy>
std::pair<index_t, index_t> basic_aby_linkage_policy<ABYLinkagePolicy>::retrieve_indices(lifted_t const& l)
{
    namespace hana = boost::hana;
    using namespace boost::hana::literals;

    assert(hana::all_of(l, [](share* s) {
        return s != nullptr;
    }) &&
    "l must be a tuple pointing to a valid value "
    "and two valid indices");

    share* out1 = circ_->PutOUTGate(l[1_c], ALL);
    share* out2 = circ_->PutOUTGate(l[2_c], ALL);

    party_.ExecCircuit();
    aby_update_counters();

    index_t r1 = out1->get_clear_value<index_t>();
    index_t r2 = out2->get_clear_value<index_t>();

    party_.Reset();

    return std::make_pair(r1, r2);
}

template<typename ABYLinkagePolicy>
share* basic_aby_linkage_policy<ABYLinkagePolicy>::put_idx_(index_t idx)
{
    return circ_->PutCONSGate(idx, bitlen_);
};

//Put plaintext into ABY Circuit and unblind it using ABY.
template<typename ABYLinkagePolicy>
share* basic_aby_linkage_policy<ABYLinkagePolicy>::put_unblinded_input_(plaintext const& input)
{
    assert(input.valid() && "input must be valid here");
    using std::swap;
    //Own input
    std::vector<uint_fast32_t> in = input;
    assert(!in.empty() && "input is valid, so in should not be empty");
    share* lhs = circ_->PutINGate(in.data(), bitlen_, role_);

    //Other party's input
    share* rhs = circ_->PutDummyINGate(bitlen_);

    //ROLE_P2 holds the blinded values and
    //ROLE_P1 holds the blinds, so input of
    //ROLE_P1 should be on the rhs, as we calculate
    //lhs - rhs <=> (x + blind) - blind
    if(ROLE_P1 == role_) swap(lhs, rhs);

    return circ_->PutSUBGate(lhs, rhs);
}


} //namespace mpo19
#endif //MPO19_LINKAGE_30062020
