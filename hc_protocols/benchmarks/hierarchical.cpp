// MIT License
//
// Copyright (c) 2021 Aditya Shridhar Hegde
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

#include <ENCRYPTO_utils/timer.h>

#include <algorithm>
#include <boost/program_options.hpp>
#include <boost/range/adaptors.hpp>
#include <fstream>
#include <limits>
#include <mpo19/channel.hpp>
#include <mpo19/paillier_wrapper.hpp>
#include <mpo19/pca.hpp>
#include <mpo19/python.hpp>
#include <nlohmann/json.hpp>
#include <sstream>
#include <stdexcept>

#include "utils.hpp"

using namespace mpo19;
using json = nlohmann::json;
namespace bpo = boost::program_options;
using boost::adaptors::sliced;

// clang-format off
bpo::options_description programOptions() {
  bpo::options_description desc("Following options are supported by config file too.");
  desc.add_options()
    ("ip", bpo::value<std::string>()->default_value("localhost"), "IP address of party1.")
    ("port", bpo::value<unsigned short>()->default_value(8000), "Port at which party1 is listening.")
    ("output,o", bpo::value<std::string>()->required(), "File to save benchmarks.")
    ("repeat,r", bpo::value<unsigned int>()->default_value(1), "Number of times to run benchmarks.")
    ("dataset", bpo::value<std::string>()->required(), "Path to dataset in npy format.")
    ("clusters", bpo::value<size_t>()->required(), "Number of clusters.")
    ("party2", "Default is party1. Using this switch denotes party 2.")
    ("protocol", bpo::value<std::string>()->required(), "Protocol name (pca | opt).");

  return desc;
}
// clang-format on

int main(int argc, char* argv[]) {
  auto prog_opts(programOptions());

  bpo::options_description cmdline(
      "Benchmark Meng et al. clustering protocols.");
  cmdline.add(prog_opts);
  cmdline.add_options()(
      "config,c", bpo::value<std::string>(),
      "configuration file for easy specification of cmd line arguments")(
      "help,h", "produce help message");

  bpo::variables_map opts;
  bpo::store(bpo::command_line_parser(argc, argv).options(cmdline).run(), opts);

  if (opts.count("help") != 0) {
    std::cout << cmdline << std::endl;
    return 0;
  }

  if (opts.count("config") > 0) {
    std::string cpath(opts["config"].as<std::string>());
    std::ifstream fin(cpath.c_str());

    if (fin.fail()) {
      std::cerr << "Could not open configuration file at " << cpath << "\n";
      return 1;
    }

    bpo::store(bpo::parse_config_file(fin, prog_opts), opts);
  }

  // Validate program options
  try {
    bpo::notify(opts);

    // Check if output file already exists
    std::ifstream ftemp(opts["output"].as<std::string>());
    if (ftemp.good()) {
      ftemp.close();
      throw std::runtime_error("Output file aready exists.");
    }
    ftemp.close();

    auto protocol = opts["protocol"].as<std::string>();
    if (protocol != "pca" && protocol != "opt") {
      throw std::runtime_error("Invalid protocol option.");
    }
  } catch (const std::exception& ex) {
    std::cerr << ex.what() << std::endl;
    return 1;
  }

  try {
    auto dataset_path = opts["dataset"].as<std::string>();
    auto num_clusters = opts["clusters"].as<size_t>();
    auto repeat = opts["repeat"].as<unsigned int>();
    auto ip = opts["ip"].as<std::string>();
    auto port = opts["port"].as<unsigned short>();
    auto party2 = opts.count("party2") != 0;
    auto output_file = opts["output"].as<std::string>();
    auto protocol = opts["protocol"].as<std::string>();

    auto dataset = load_dataset(dataset_path);

    json bench;
    bench["details"] = {{"dataset", dataset_path},
                        {"repeat", repeat},
                        {"clusters", num_clusters}};

    std::cout << "--- Details ---" << std::endl;
    print_json(bench["details"]);

    bench["benchmarks"] = json::array();
    auto &bench_data = bench["benchmarks"];

    for (unsigned int i = 0; i < repeat; ++i) {
      random_permutation perm(dataset.size());
      perm.permute_inplace(dataset);

      if (party2) {
        P2 p2{ip, port};
        auto input = dataset | sliced(0, dataset.size() / 2);

        if (protocol == "pca") {
          BENCHMARK(bench_data, p2.execute_pca, input, num_clusters);
        } else {
          BENCHMARK(bench_data, p2.execute_opt, input, num_clusters);
        }

      } else {
        P1 p1{port};
        auto input = dataset | sliced(dataset.size() / 2, dataset.size());

        if (protocol == "pca") {
          BENCHMARK(bench_data, p1.execute_pca, input, num_clusters);
        } else {
          BENCHMARK(bench_data, p1.execute_opt, input, num_clusters);
        }
      }

      save_json(bench, output_file);
      std::cout << std::endl;
    }

    std::cout << "\n--- Stats ---\n";
    bench["stats"]["vmpeak"] = getProcStatus("VmPeak:");
    bench["stats"]["vmhwm"] = getProcStatus("VmHWM:");
    save_json(bench, output_file);
    print_json(bench["stats"]);

  } catch (const std::exception& ex) {
    std::cerr << ex.what() << "\nFatal error" << std::endl;
    return 1;
  }

  return 0;
}
