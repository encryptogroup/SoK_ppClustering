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

#ifndef MPO19_PCA_17072020
#define MPO19_PCA_17072020

#include <mpo19/cluster.hpp>
#include <mpo19/point.hpp>
#include <mpo19/utility.hpp>
#include <mpo19/typedefs.hpp>
#include <mpo19/config.hpp>
#include <mpo19/pcure.hpp>
#include <tuple>
#include <vector>
#include <string>
#include <experimental/propagate_const>
#include <boost/range/algorithm.hpp>

struct P1;
struct P2;

template<typename Party>
struct basic_party {
private:
    using self_ = basic_party;
    using derived_ = Party;

public:
    using output_t = std::tuple<
                     std::vector<cluster>
                     , std::vector<point<double>>
                     , std::vector<size_t>>;

    void print_sigma() const;

    output_t execute_pca(array_t<point<int64_t>>& input, size_t t);
    output_t execute_opt(array_t<point<int64_t>>& input, size_t t);

    template<typename Range>
    output_t execute_pca(Range const& range, size_t t)
    {
        array_t<point<int64_t>> in{boost::extents[range.size()]};
        boost::copy(range, in.begin());
        return execute_pca(in, t);
    }

    template<typename Range>
    output_t execute_opt(Range const& range, size_t t)
    {
        array_t<point<int64_t>> in{boost::extents[range.size()]};
        boost::copy(range, in.begin());
        return execute_opt(in, t);
    }

    std::vector<cluster> pcure0(std::vector<point<int64_t>> const& input
                                , size_t s, size_t p, size_t q, size_t t1, size_t t2, size_t t);
    std::vector<cluster> pcure0_opt(std::vector<point<int64_t>> const& input
                                    , size_t s, size_t p, size_t q, size_t t1, size_t t2, size_t t);

    template<typename Range>
    std::vector<cluster> pcure0(Range const& input
                                , size_t s, size_t p, size_t q, size_t t1, size_t t2, size_t t)
    {
        std::vector<point<int64_t>> in{input.begin(), input.end()};
        return pcure0(in, s, p, q, t1, t2, t);
    }

    template<typename Range>
    std::vector<cluster> pcure0_opt(Range const& input
                                    , size_t s, size_t p, size_t q, size_t t1, size_t t2, size_t t)
    {
        std::vector<point<int64_t>> in{input.begin(), input.end()};
        return pcure0_opt(in, s, p, q, t1, t2, t);
    }


    std::vector<cluster> pcure1(std::vector<point<int64_t>> const& input
                                , size_t s, size_t p, size_t q, size_t t1, size_t t2, size_t t);
    std::vector<cluster> pcure1_opt(std::vector<point<int64_t>> const& input
                                    , size_t s, size_t p, size_t q, size_t t1, size_t t2, size_t t);

    template<typename Range>
    std::vector<cluster> pcure1(Range const& input
                                , size_t s, size_t p, size_t q, size_t t1, size_t t2, size_t t)
    {
        std::vector<point<int64_t>> in{input.begin(), input.end()};
        return pcure1(in, s, p, q, t1, t2, t);
    }

    template<typename Range>
    std::vector<cluster> pcure1_opt(Range const& input
                                    , size_t s, size_t p, size_t q, size_t t1, size_t t2, size_t t)
    {
        std::vector<point<int64_t>> in{input.begin(), input.end()};
        return pcure1_opt(in, s, p, q, t1, t2, t);
    }

    std::vector<cluster> pcure2(array_t<point<int64_t>>& input
                                , size_t s, size_t q, size_t t1, size_t t2, size_t t);
    std::vector<cluster> pcure2_opt(array_t<point<int64_t>>& input
                                    , size_t s, size_t q, size_t t1, size_t t2, size_t t);

    template<typename Range>
    std::vector<cluster> pcure2(Range const& input
                                , size_t s, size_t q, size_t t1, size_t t2, size_t t)
    {
        array_t<point<int64_t>> in{boost::extents[input.size()]};
        boost::copy(input, in.begin());
        return pcure2(in, s, q, t1, t2, t);
    }
    
    template<typename Range>
    std::vector<cluster> pcure2_opt(Range const& input
                                    , size_t s, size_t q, size_t t1, size_t t2, size_t t)
    {
        array_t<point<int64_t>> in{boost::extents[input.size()]};
        boost::copy(input, in.begin());
        return pcure2(in, s, q, t1, t2, t);
    }

protected:

    ~basic_party() = default;

private:

    void clustering(size_t t);
    void clustering_opt(size_t t);
    void insert_other_parties_centroids(std::vector<point<double>>& centroids);
    std::vector<point<double>> pca_on_A_clusters(mpo19::pcure_clusters<int64_t>& clusters, size_t t2, size_t t);
    std::vector<point<double>> opt_on_A_clusters(mpo19::pcure_clusters<int64_t>& clusters, size_t t2, size_t t);
    std::vector<point<double>> remove_outliers(output_t& out, mpo19::pcure_clusters<int64_t>& clusters, size_t t2);
};

struct P1 : basic_party<P1> {
private:
    using self_ = P1;
    using base_ = basic_party<self_>;
    friend base_;

public:
    P1(unsigned short port);

    void print_L();

    ~P1();

private:

    struct impl_;
    std::experimental::propagate_const<
    std::unique_ptr<impl_>> pimpl_;
};

struct P2  : basic_party<P2> {
private:
    using self_ = P2;
    using base_ = basic_party<self_>;
    friend base_;

public:

    P2(std::string const& destination, unsigned short port);

    ~P2();

    void print_L();

private:

    struct impl_;
    std::experimental::propagate_const<
    std::unique_ptr<impl_>> pimpl_;
};

#endif //MPO19_PCA_17072020
