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

#ifndef MPO19_CLUSTERING_22072020
#define MPO19_CLUSTERING_22072020

#include <mpo19/typedefs.hpp>
#include <mpo19/cluster.hpp>
#include <boost/range/algorithm.hpp>

namespace mpo19{
template<typename Linkage>
struct hierarchical_clustering : public Linkage {

    using dendrogram_t  = std::vector<cluster>;

    using Linkage::Linkage;

    dendrogram_t operator()(size_t t)
    {
        using std::swap;

        size_t n = Linkage::get_input_size();

        assert(t < n && "t must be strictly smaller than n");

        dendrogram_t dendrogram;

        //Initialize dendrogram
        for(index_t i = 0; i != index_t(n); ++i) {
            dendrogram.emplace_back(cluster::leaf_t(i));
        }

        for(size_t l = 0; l < n - t; ++l) {
            //get indices of clusters to be merged
            auto p = get_merge_indices();

            size_t i = std::get<0>(p);
            size_t j = std::get<1>(p);

            assert(i < j && "By convention i should be smaller than j");
            assert(0 <= i && i < n && "i should be between 0 and n");
            assert(0 <= j && j < n && "j should be between 0 and n");

            //merge the two clusters at position i and j
            merge(dendrogram[i], dendrogram[j], l + 1);
        }
        
        assert(dendrogram.size() == n && "dendrogram should be of size n here");
        //Erase all clusters that only refer to other clusters through a pointer
        dendrogram.erase(
        boost::remove_if(dendrogram, [](cluster const& c) {
            return !c.is_cluster();
        })
        , dendrogram.end());
        
        assert(dendrogram.size() == t && "dendrogram should be of size t here");

        return dendrogram;
    }
    
    dendrogram_t partial_evaluation(size_t t, dendrogram_t dendrogram = dendrogram_t{})
    {
        using std::swap;
        using boost::range::count_if;

        size_t n = Linkage::get_input_size();

        assert(t < n && "t must be strictly smaller than n");

        size_t start;
        if(dendrogram.empty()) {
            //Initialize dendrogram
            for(index_t i = 0; i != index_t(n); ++i) {
                dendrogram.emplace_back(cluster::leaf_t(i));
            }
            start = 0;
        }
        else {
            assert(dendrogram.size() == n);
            start = n - count_if(dendrogram, [](cluster const& c){
                return c.is_cluster();
            });
        }

        for(size_t l = start; l < n - t; ++l) {
            //get indices of clusters to be merged
            auto p = get_merge_indices();

            size_t i = std::get<0>(p);
            size_t j = std::get<1>(p);

            assert(i < j && "By convention i should be smaller than j");
            assert(0 <= i && i < n && "i should be between 0 and n");
            assert(0 <= j && j < n && "j should be between 0 and n");

            //merge the two clusters at position i and j
            merge(dendrogram[i], dendrogram[j], l + 1);
        }

        assert(dendrogram.size() == n && "dendrogram should be of size n here");

        return dendrogram;
    }

private:
    using Linkage::get_merge_indices;
};

}
#endif //MPO19_CLUSTERING_22072020
