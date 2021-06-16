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

#ifndef MPO19_PCURE_02082020
#define MPO19_PCURE_02082020

#include <mpo19/typedefs.hpp>
#include <mpo19/random.hpp>
#include <mpo19/clustering.hpp>
#include <mpo19/linkage.hpp>
#include <mpo19/paillier_wrapper.hpp>
#include <boost/phoenix.hpp>

#include <iostream>
#include <type_traits>

namespace mpo19
{

constexpr struct calc_distance_t {
    size_t operator()(point<int64_t> const& lhs
                      , point<int64_t> const& rhs) const
    {
        using namespace boost::phoenix::placeholders;
        assert(lhs.dimension() == rhs.dimension());
        return apply_on_elements(arg1 * arg1, lhs - rhs).fold();
    }

    double operator()(point<double> const& lhs
                      , point<int64_t> const& rhs) const
    {
        using namespace boost::phoenix::placeholders;
        assert(lhs.dimension() == rhs.dimension());
        point<double> res(lhs.dimension());
        apply_on_elements(res, arg1 - arg2, lhs, rhs);
        apply_on_elements(res, arg1 * arg1, res);
        return res.fold();
    }

} calc_distance;

template<typename T, typename LinkagePolicy>
struct cluster_dist_t : public LinkagePolicy {

    template<typename... Args>
    cluster_dist_t(std::vector<point<T>> const& values, Args&&... args)
        : LinkagePolicy{std::forward<Args>(args)...}, values_{values} {}

    auto operator()(cluster const& c1, cluster const& c2)
    {
        auto c1_indices = get_indices(c1);
        auto c2_indices = get_indices(c2);

        size_t c1_indices_size = c1_indices.size();
        size_t c2_indices_size = c2_indices.size();

        assert(0 < c1_indices_size);
        assert(0 < c2_indices_size);

        auto res = calc_distance(values_[c1_indices[0]], values_[c2_indices[0]]);

        size_t i = 0, j = 1;
        for(; i != c1_indices_size; ++i, j = 0) {
            for(; j != c2_indices_size; ++j) {
                assert( !(i == 0 && j == 0) &&
                        "We already calculated the distance for (0, 0)");
                res =
                    LinkagePolicy::cluster_distance(
                        res, calc_distance(values_[c1_indices[i]], values_[c2_indices[j]]));
            }
        }

        return res;
    }

private:
    std::vector<point<T>> const& values_;
};

template<typename T>
struct pcure_clusters {
    //We use a vector<vector>, as the rows might have different sizes
    std::vector<cluster> partition_dendrograms;
    std::vector<point<T>> values;
};

template<typename T>
std::ostream& operator<<(std::ostream& os, pcure_clusters<T> const& clusters)
{
    os << '(';
    for(cluster const& c : clusters.partition_dendrograms) {
        os << c;
    }
    os << ')';
    os << '(';
    for(point<T> const& p : clusters.values) {
        os << p;
    }
    os << ')';
    return os;
}

template<typename T>
std::istream& operator>>(std::istream& is, pcure_clusters<T>& clusters)
{
    clusters.partition_dendrograms.clear();
    clusters.values.clear();

    char c = is.get();
    assert(c == '(');

    while(is.peek() != ')') {
        cluster clust;
        is >> clust;
        clusters.partition_dendrograms.emplace_back(std::move(clust));
    }

    c = is.get();
    assert(c == ')');
    c = is.get();
    assert(c == '(');

    while(is.peek() != ')') {
        point<T> p;
        is >> p;
        clusters.values.emplace_back(std::move(p));
    }

    return is;
}

template<template<class> class LinkagePolicy>
struct clear_clustering_t {

    template<typename RandomAccessRange>
    auto operator()(RandomAccessRange const& partition, size_t t)
    {
        using dist_table_t = table_t<decltype(calc_distance(partition[0], partition[0]))>;
        size_t n = partition.size();
        dist_table_t dist_matrix{boost::extents[n][n]};

        for(index_t  i = 0; i != index_t(n); ++i) {
            for(index_t j = 0; j != index_t(n); ++j) {
                dist_matrix[i][j] = calc_distance(partition[i], partition[j]);
            }
        }

        return hierarchical_clustering<LinkagePolicy<dist_table_t>> {dist_matrix}(t);
    }

    template<typename T>
    auto operator()(pcure_clusters<T> const& clusters, size_t t)
    {
        size_t n = clusters.partition_dendrograms.size();
        table_t<T> dist_matrix{boost::extents[n][n]};

        hierarchical_clustering<cluster_dist_t<T, LinkagePolicy<table_t<T>>>>
        hc{clusters.values, dist_matrix};

        cluster_dist_t<T, LinkagePolicy<table_t<T>>>& cluster_dist = hc;

        for(index_t i = 0; i != index_t(n); ++i) {
            for(index_t j = 0; j != index_t(n); ++j) {
                if(i == j) {
                    dist_matrix[i][j] = 0;
                } else {
                    dist_matrix[i][j] = cluster_dist(
                                            clusters.partition_dendrograms[i]
                                            , clusters.partition_dendrograms[j]);
                }
            }
        }

        return hc(t);
    }
};

template<typename Range>
std::vector<cluster> assign_clusters(
    Range const& input
    , std::vector<point<double>> const& centroids)
{
    assert(!centroids.empty());
    std::vector<cluster> res(centroids.size());
    size_t rank = 0;

    for(size_t i = 0; i != input.size(); ++i) {
        size_t min_index = 0;
        double min_distance = calc_distance(centroids[0], input[i]);
        for(size_t j = 1; j != centroids.size(); ++j) {
            double distance = calc_distance(centroids[j], input[i]);
            if(distance < min_distance) {
                min_index = j;
                min_distance = distance;
            }
        }
        if(res[min_index].is_cluster()) ++rank;
        cluster singleton{cluster::leaf_t(i)};
        merge(res[min_index], singleton, rank);
        assert(!singleton.is_cluster());
    }
    
    res.erase(std::remove_if(res.begin(), res.end(), [](cluster const& c){
        return !c.is_cluster();
    }), res.end());

    return res;
}

template<typename T>
void transform_pcure_dendrogram(
    pcure_clusters<T>& clusters
    , std::vector<cluster>& dendrogram
    , size_t base_rank)
{
    for(cluster& c : dendrogram) {
        assert(c.is_cluster());
        std::vector<index_t> indices = get_indices(c);
        assert(clusters.partition_dendrograms[indices[0]].is_cluster());
        c = std::move(clusters.partition_dendrograms[indices[0]]);

        for(size_t i = 1, rank = base_rank; i != indices.size(); ++i, ++rank) {
            index_t idx = indices[i];
            assert(idx < clusters.partition_dendrograms.size());
            assert(clusters.partition_dendrograms[idx].is_cluster());
            merge(c, clusters.partition_dendrograms[idx], rank);
        }

        assert(c.is_cluster());
    }
}

template<typename ClusteringPolicy>
struct pcure : public ClusteringPolicy {

    template<typename... Args>
    pcure(size_t s, size_t p, size_t q, size_t t1, size_t t2, size_t t, Args&&... args)
        : ClusteringPolicy{std::forward<Args>(args)...}, s{s}, p{p}, q{q}, t1{t1}, t2{t2}, t{t} {}

    template<typename T>
    pcure_clusters<T> A_clustering(std::vector<point<T>>& input)
    {
        using boost::adaptors::sliced;
        using boost::adaptors::filtered;
        using boost::range::copy;

        struct add_partition_idx_t : boost::static_visitor<void> {

            add_partition_idx_t(size_t partition_idx)
                : partition_idx_{partition_idx} {}

            result_type operator()(cluster::union_t& u) const
            {
                u.get_left().visit(*this);
                u.get_right().visit(*this);
            }

            result_type operator()(cluster::leaf_t& l) const
            {
                l += partition_idx_;
            }

        private:
            size_t partition_idx_;
        };

        assert(input.size() > s &&
               "Cannot take a sample of size s with s being bigger than the input size");
        assert(s > p &&
               "Cannot partition s into more than s partitions");
        assert(s/(p*q) > 0 &&
               "The target clusters shall be grater than zero, i.e. s/(p*q) > 0");

        size_t const target_clusters = s/(p*q);

        pcure_clusters<T> res;

        res.partition_dendrograms.reserve(s/q);

        //if s is not a multiple of p, then the last s % p items won't be included in any partition
        auto sample = random_sample(input, s);

        for(size_t i = 0; i != p; ++i) {
            size_t const slice_first = i * s/p;
            size_t const slice_last = (i+1) * s/p;
            assert(slice_first < slice_last);

            auto partition = sample | sliced(slice_first, slice_last);

            auto dendrogram = static_cast<ClusteringPolicy&>(*this)(partition, target_clusters);
            auto filtered_dendrogram = dendrogram | filtered([&](cluster const& c) -> bool {
                return get_indices(c).size() >= t1;
            });

            for(cluster& c : filtered_dendrogram) {
                if(slice_first != 0)
                    c.visit(add_partition_idx_t{slice_first});
                res.partition_dendrograms.emplace_back(std::move(c));
            }
        }
        res.values = std::move(sample);

        return res;
    }

    template<typename T>
    pcure_clusters<T> B_clustering(pcure_clusters<T>& clusters)
    {
        using boost::range::remove_if;

        assert(clusters.partition_dendrograms.size() <= s/q);

        //The result in this context is the merge history of input clusters,
        //therefore we need to convert it to a merge history of indices
        auto dendrogram = static_cast<ClusteringPolicy&>(*this)(clusters, t);
        size_t base_rank = s/p - s/(p*q) + 1;

        for(cluster& c : dendrogram) {
            assert(c.is_cluster());
            std::vector<index_t> indices = get_indices(c);
            assert(clusters.partition_dendrograms[indices[0]].is_cluster());
            c = std::move(clusters.partition_dendrograms[indices[0]]);

            for(size_t i = 1, rank = base_rank; i != indices.size(); ++i, ++rank) {
                assert(clusters.partition_dendrograms[indices[i]].is_cluster());
                merge(c, clusters.partition_dendrograms[indices[i]], rank);
            }

            assert(c.is_cluster());
        }

        dendrogram.erase(remove_if(dendrogram, [&](cluster const& c) {
            return get_indices(c).size() < t2;
        })
        , dendrogram.end());

        clusters.partition_dendrograms = std::move(dendrogram);

        return std::move(clusters);
    }

    template<typename T>
    std::vector<point<double>> get_cluster_centroids(std::vector<point<T>> const& input)
    {
        using namespace boost::phoenix::placeholders;

        std::vector<point<T>> in{input};
        auto C_A = A_clustering(in);
        auto C_B = B_clustering(C_A);

        std::vector<point<double>> centroids;
        centroids.reserve(C_B.partition_dendrograms.size());
        for(cluster const& c : C_B.partition_dendrograms) {
            std::vector<index_t> indices = get_indices(c);
            point<double> sum(C_B.values.front().dimension(), 0.0);
            for(index_t i : indices) {
                apply_on_elements(sum, arg1 + arg2, sum, C_B.values[i]);
            }
            apply_on_elements(sum, arg1/indices.size(), sum);
            centroids.emplace_back(std::move(sum));
        }

        return centroids;
    }

    template<typename T>
    std::vector<cluster> assign_clusters(
        std::vector<point<T>> const& input
        , std::vector<point<double>> const& centroids)
    {
        return ::mpo19::assign_clusters(input, centroids);
    }

private:
    size_t s, p, q, t1, t2, t;
};

}

#endif //MPO19_PCURE_02082020
