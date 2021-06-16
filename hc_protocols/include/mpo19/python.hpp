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

#ifndef MPO19_PYTHON_22072020
#define MPO19_PYTHON_22072020

#define BOOST_BIND_GLOBAL_PLACEHOLDERS
#include <boost/python.hpp>
#include <mpo19/cluster.hpp>

namespace mpo19
{

namespace detail
{


template<typename T>
T convert(float f)
{
    return static_cast<T>(round(f));
}

}

struct python {

    python();

    std::vector<cluster> hierarchical_clustering(
        table_t<size_t>& distance_matrix, size_t t, char const* alg = "complete");

    std::vector<point<int>> load_iris();

    python& exec_file(const char* filename)
    {
        boost::python::exec_file(filename, *main_namespace_, *main_namespace_);
        return *this;
    }

    struct python_object_wrapper {

        template<typename T>
        operator std::vector<point<T>>() const
        {
            using namespace boost::python;

            std::vector<point<T>> res;

            auto shape = obj_.attr("shape");
            size_t vec_size = extract<size_t>(shape[0]);
            res.reserve(vec_size);

            size_t point_dim = extract<size_t>(shape[1]);

            for(size_t i = 0; i != vec_size; ++i) {
                point<T> p(point_dim);
                for(size_t j = 0; j != point_dim; ++j) {
                    float f = extract<float>(obj_[i][j]);
                    p[j] = detail::convert<T>(f);
                }
                res.emplace_back(std::move(p));
            }

            return res;
        }
        
        template<typename T>
        operator T() const
        {
            return boost::python::extract<T>(obj_);
        }

        boost::python::object obj_;
    };

    python_object_wrapper operator[](const char* name)
    {
        return {(*main_namespace_)[name]};
    }


private:
    static int cm_to_mm(float cm);
    boost::python::object* main_namespace_;

};
}
#endif //MPO19_PYTHON_22072020
