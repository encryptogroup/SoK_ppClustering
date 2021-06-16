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

#include <mpo19/paillier_wrapper.hpp>
#include <mpo19/pca.hpp>
#include <mpo19/python.hpp>

#include <chrono>
#include <thread>
#include <iostream>
#include <sstream>

#include <boost/range/adaptors.hpp>

using namespace mpo19;

void benchmark_pca(std::vector<point<int64_t>> const& input, size_t target_clusters, std::string const& ip, unsigned short port)
{
    using boost::adaptors::sliced;
    using namespace std::chrono;
    using namespace std::literals::chrono_literals;
    std::cout << "PCA on input size: " << input.size() << std::endl;
    random_permutation p(input.size());
    std::thread other_party{[&]{
            P1 p1{port};
            p1.execute_pca(input | sliced(input.size()/2, input.size()), target_clusters);
        }};
    other_party.detach();

    P2 p2{ip, port};

    auto start = steady_clock::now();
    p2.execute_pca(input | sliced(0, input.size()/2), target_clusters);
    auto end = steady_clock::now();

    std::cout << "\nPCA for input of size "<< input.size() << " took: " << duration<double> {end - start} .count() << " s\n" << std::endl;
}

void benchmark_opt(std::vector<point<int64_t>> const& input, size_t target_clusters, std::string const& ip, unsigned short port)
{
    using boost::adaptors::sliced;
    using namespace std::chrono;
    using namespace std::literals::chrono_literals;
    std::cout << "OPT on input size: " << input.size() << std::endl;
    random_permutation p(input.size());
    std::thread other_party{[&]{
            P1 p1{port};
            p1.execute_opt(input | sliced(input.size()/2, input.size()), target_clusters);
        }};
    other_party.detach();

    P2 p2{ip, port};

    auto start = steady_clock::now();
    p2.execute_opt(input | sliced(0, input.size()/2), target_clusters);
    auto end = steady_clock::now();

    std::cout << "\nOPT for input of size "<< input.size() << " took: " << duration<double> {end - start} .count() << " s\n" << std::endl;
}

#include <mpo19/pcure.hpp>

#include <iostream>
#include <sstream>
#include <iomanip>

int main()
{
    using namespace boost::python;

    try {
        python py;
        //This path assumes that the executable is in mpo19/<some_dir>
        //you might want to update this path
        py.exec_file("../python/read_data.py");

        //read_data.py sets the variables iris, wine, ip, port

        //Conversion to std::vector<point<int64_t>> asserts that
        //the python variable contains a python object obj, such that
        //obj[i][j] is a number.
        std::vector<point<int64_t>> iris = py["iris"];
        std::vector<point<int64_t>> wine = py["wine"];
        //ip is set to a string
        std::string ip = py["ip"];
        //port is set to an integer number <= 65535
        unsigned short port = py["port"];
        (void) port;

        constexpr size_t s = 75;
        constexpr size_t p = 5;
        constexpr size_t q = 3;
        constexpr size_t t1 = 3;
        constexpr size_t t2 = 5;
        constexpr size_t t = 6;

        for(auto const& p : iris){
            std::cout << p << " ";
        }
        std::cout << "\n" << std::endl;


        pcure<clear_clustering_t<max_linkage>> pcure{s, p, q, t1, t2, t};
        auto centroids = pcure.get_cluster_centroids(iris);
        auto clusters = pcure.assign_clusters(iris, centroids);

        for(auto const& p : centroids) {
            std::cout << p << " ";
        }
        std::cout << std::endl;

        for(auto const& p : clusters) {
            std::cout << p << "\n";
        }
        std::cout << std::endl;

        std::cout << "\n\n";

        pcure_clusters<int64_t> p_c;
        p_c.partition_dendrograms = std::move(clusters);
        p_c.values = std::move(iris);

        std::cout << p_c << std::endl;

        std::stringstream ss;
        ss << p_c;
        ss >> p_c;

        std::cout << "\n\n" << p_c << std::endl;

        return 0;

    } catch(boost::python::error_already_set const&) {
        PyErr_Print();
    }
    try {
        python py;
        //This path assumes that the executable is in mpo19/<some_dir>
        //you might want to update this path
        py.exec_file("../python/read_data.py");

        //read_data.py sets the variables iris, wine, ip, port

        //Conversion to std::vector<point<int64_t>> asserts that
        //the python variable contains a python object obj, such that
        //obj[i][j] is a number.
        std::vector<point<int64_t>> iris = py["iris"];
        std::vector<point<int64_t>> wine = py["wine"];
        //ip is set to a string
        std::string ip = py["ip"];
        //port is set to an integer number <= 65535
        unsigned short port = py["port"];

        benchmark_opt(iris, 5, ip, port);
        benchmark_opt(wine, 5, ip, port);
        benchmark_pca(iris, 5, ip, port);
        benchmark_pca(wine, 5, ip, port);


    } catch(boost::python::error_already_set const&) {
        PyErr_Print();
    }
}
