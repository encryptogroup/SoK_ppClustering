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

#include <mpo19/python.hpp>

namespace mpo19
{

typedef char const python_code[];

python::python()
{
    using namespace boost::python;
    static object main_module = (Py_Initialize(), import("__main__"));
    static object main_namespace = main_module.attr("__dict__");
    main_namespace_ = std::addressof(main_namespace);
}

std::vector<cluster> python::hierarchical_clustering(
    table_t<size_t>& distance_matrix, size_t t, char const* alg)
{
    using namespace boost::python;

    assert(distance_matrix.shape()[0] == distance_matrix.shape()[1]
           && "distance_matrix should be a squared matrix");

    size_t n = distance_matrix.shape()[0];

    std::stringstream py;
    py <<
       "import numpy as np\n"
       "from sklearn.cluster import AgglomerativeClustering\n"
       "distance_matrix = np.array([";

    for(size_t i = 0; i < n; ++i) {
        if(i != 0) py << ", ";
        py << "[";
        for(size_t j = 0; j < n; ++j) {
            if(j != 0) py << ", ";
            py << distance_matrix[i][j];
        }
        py << "]";
    }
    py << "])\n";
    py <<
       "cluster = AgglomerativeClustering("
       "                     n_clusters = " << t << ","
       "                     affinity = 'precomputed',"
       "                     linkage = '" << alg << "').fit(distance_matrix)\n"
       "dendrogram = cluster.children_\n";

    object ignored = exec(py.str().c_str(), (*main_namespace_));

    object py_dendrogram = (*main_namespace_)["dendrogram"];
    auto py_dendrogram_shape = py_dendrogram.attr("shape");

    assert(n - 1 == extract<size_t>(py_dendrogram_shape[0])
           && "cluster.children_ : array-like of shape (n_samples-1, 2)");

    assert(2 == extract<size_t>(py_dendrogram_shape[1])
           && "cluster.children_ : array-like of shape (n_samples-1, 2)");

    std::vector<cluster> dendrogram;
    dendrogram.reserve(n);
    std::cout << std::endl;
    for(size_t i = 0; i != n; ++i) {
        dendrogram.emplace_back(cluster::leaf_t(i));
    }
    assert(dendrogram.size() == n);

    std::vector<size_t> merged_clusters_idices(n - 1);

    for(size_t i = 0; i != n - t; ++i) {
        size_t left = extract<size_t>(long_(py_dendrogram[i][0]));
        size_t right = extract<size_t>(long_(py_dendrogram[i][1]));

        assert(left != right);

        //Get clusters to merge either from dendrogram (left|right < n)
        //or from already merged clusters (left|right >= n)
        if(left >= n) {
            left = left - n;
            assert(left < i);
            left = merged_clusters_idices[left];
        }

        if(right >= n) {
            right = right - n;
            assert(right < i);
            right = merged_clusters_idices[right];
        }
        assert(left != right);

        if(left > right) {
            std::swap(left,  right);
        }
        assert(left < right);

        merge(dendrogram[left], dendrogram[right], i + 1);
        merged_clusters_idices[i] = left;
    }
    auto new_end = boost::remove_if(dendrogram, [](cluster const& c) {
        return !c.is_cluster();
    });

    dendrogram.resize(std::distance(dendrogram.begin(), new_end));

    return dendrogram;
}

std::vector<point<int>> python::load_iris()
{
    using namespace boost::python;
    constexpr python_code load_iris =
        "from sklearn.datasets import load_iris\n"
        "iris = load_iris()\n"
        "X = iris.data";

    object ignored = exec(load_iris, *main_namespace_);

    object py_iris = (*main_namespace_)["X"];
    auto py_iris_shape = py_iris.attr("shape");

    std::vector<point<int>> res;
    size_t iris_size = extract<size_t>(py_iris_shape[0]);
    res.reserve(iris_size);

    for(size_t i = 0; i != iris_size; ++i) {
        size_t d = extract<size_t>(py_iris_shape[1]);
        point<int> p(d);
        for(size_t j = 0; j != d; ++j) {
            p[j] = cm_to_mm(extract<float>(py_iris[i][j]));
        }
        res.emplace_back(std::move(p));
    }
    return res;
}

int python::cm_to_mm(float cm)
{
    float mm = cm * 10;
    assert(mm - floor(mm) < 0.01f
           && "Centimeter should not have more than one decimal place");
    return static_cast<int>(mm);
}

}
