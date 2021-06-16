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

#include <HEAAN.h>
#include <cnpy.h>
#include <errno.h>
#include <sys/stat.h>

#include <boost/program_options.hpp>
#include <chrono>
#include <ckp19/dataset.hpp>
#include <ckp19/mean_shift.hpp>
#include <ckp19/plaintext_clustering.hpp>
#include <filesystem>
#include <fstream>
#include <nlohmann/json.hpp>
#include <numeric>
#include <string>

#include "utils.hpp"

namespace fs = std::filesystem;
namespace bpo = boost::program_options;
using json = nlohmann::json;
using TimeUnit = std::chrono::system_clock::time_point;

void benchmarkClustering(const bpo::variables_map &opts) {
  long seed = opts["seed"].as<long>();
  long repeat = opts["repeat"].as<long>();
  long threads = opts["threads"].as<long>();
  std::string dataset = opts["dataset"].as<std::string>();

  srand(seed);
  SetNumThreads(threads);

  bool save_data = false;
  fs::path output_dir;
  if (opts.count("output") > 0) {
    save_data = true;
    output_dir = fs::path(opts["output"].as<std::string>());
  }

  long logq = opts["logq"].as<long>();
  long logp = opts["logp"].as<long>();
  long logt = opts["logt"].as<long>();
  long logq_boot = opts["logq-boot"].as<long>();

  long num_dusts = opts["dusts"].as<long>();
  long iterations = opts["iterations"].as<long>();
  long mode_kdegree = opts["mode-kdegree"].as<long>();
  long label_kdegree = opts["label-kdegree"].as<long>();
  long mode_inv_steps = opts["mode-invsteps"].as<long>();
  long label_inv_steps = opts["label-invsteps"].as<long>();
  long minidx_t = opts["min-idx-t"].as<long>();

  json info;

  TimePoint start;
  Ring ring;
  SecretKey secretKey(ring);
  Scheme scheme(secretKey, ring);
  scheme.addLeftRotKeys(secretKey);
  scheme.addRightRotKeys(secretKey);
  TimePoint end;
  info["stats"]["init_scheme"] = end - start;

  Dataset ds = Dataset::loadNpy(dataset);
  ds.rescaleDatasetInPlace(0.5);

  long logl = std::ceil(std::log2(ds[0].size() * num_dusts));
  start = TimePoint();
  scheme.addBootKey(secretKey, logl, logq_boot + logt);
  end = TimePoint();
  info["stats"]["gen_bootstrap_key"] = end - start;

  info["details"] = {{"seed", seed},
                     {"repeat", repeat},
                     {"threads", threads},
                     {"dataset", dataset},
                     {"logN", logN},
                     {"logQ", logQ},
                     {"logq", logq},
                     {"logp", logp},
                     {"logt", logt},
                     {"logq-boot", logq_boot},
                     {"logl", logl},
                     {"dusts", num_dusts},
                     {"iterations", iterations},
                     {"mode-kdegree", mode_kdegree},
                     {"label-kdegree", label_kdegree},
                     {"mode-invsteps", mode_inv_steps},
                     {"label-invsteps", label_inv_steps},
                     {"min-idx-t", minidx_t}};

  std::cout << "--- Information ---\n";
  for (json::iterator it = info["details"].begin(); it != info["details"].end();
       ++it) {
    std::cout << it.key() << " : " << it.value() << "\n";
  }
  std::cout << std::endl;

  std::cout << "--- Setup ---\n";
  for (json::iterator it = info["stats"].begin(); it != info["stats"].end();
       ++it) {
    std::cout << it.key() << " : " << it.value() << "\n";
  }
  std::cout << std::endl;

  info["benchmarks"] = json::array();

  fs::create_directories(output_dir);
  auto info_file = output_dir / "info.json";

  for (long i = 0; i < repeat; ++i) {
    std::cout << "--- Repetition " << i + 1 << " ---\n";

    json iter;

    EncryptedDataset eds;
    start = TimePoint();
    ds.encryptSIMD(eds, scheme, num_dusts, Nh, logp, logq);
    end = TimePoint();
    iter["encrypt_dataset_time"] = end - start;

    if (i == 0) {
      info["details"]["plaintext_points_per_ciphertext"] = eds.npoints;
      info["details"]["num_ciphertexts"] = eds.points.size();
      info["details"]["slots_used"] = eds.log_slots;
    }

    MeanShift meanshift(scheme, eds.log_slots, logq, logq_boot, logp, logt,
                        seed + i);

    std::vector<Ciphertext> clabels;

    TimePoint start;
    meanshift.clusterSIMD(clabels, eds.points, eds.dim, eds.npoints, num_dusts,
                          iterations, mode_kdegree, label_kdegree,
                          mode_inv_steps, label_inv_steps, minidx_t);
    TimePoint end;
    iter["cluster_time"] = end - start;

    std::vector<double> olabels(ds.size() * num_dusts);
    double decrypt_time = 0;
    for (size_t p = 0; p < eds.points.size(); ++p) {
      TimePoint start;
      std::unique_ptr<complex<double>[]> output(
          scheme.decrypt(secretKey, clabels[p]));
      TimePoint end;
      decrypt_time += end - start;

      for (int i = 0; i < num_dusts; ++i) {
        for (int j = 0; j < eds.npoints; ++j) {
          // Store the final label as the value in the first dimension
          olabels[p * eds.npoints * num_dusts + j * num_dusts + i] =
              output[i * eds.npoints * eds.dim + j * eds.dim].real();
        }
      }
    }
    iter["decrypt_output_time"] = decrypt_time;

    if (save_data) {
      info["benchmarks"].push_back(iter);
      saveJson(info, info_file);

      auto data_file =
          (output_dir / ("rep_" + std::to_string(i) + ".npy")).string();
      cnpy::npy_save(data_file, olabels.data(),
                     {ds.size(), static_cast<unsigned long>(num_dusts)}, "w");
    }
  }

  if (save_data) {
    // Highest virutal memory usage and physical memory usage
    info["stats"]["vmpeak"] = getProcStatus("VmPeak:");
    info["stats"]["vmhwm"] = getProcStatus("VmHWM:");
    saveJson(info, info_file);
  }
}

// clang-format off
bpo::options_description generic_program_options() {
  bpo::options_description desc("Followig options are supported by config file too");
  desc.add_options()
    ("output,o", bpo::value<std::string>(), "Directory to save benchmarks and other data.")
    ("seed", bpo::value<long>()->default_value(2602), "Seed used for RNG.")
    ("repeat,r", bpo::value<long>()->default_value(1), "Number of times to run benchmarks.")
    ("threads", bpo::value<long>()->default_value(1), "Number of threads to use.")
    ("dataset", bpo::value<std::string>()->required(), "Path to dataset in npy format.");

  desc.add_options()
    ("logq", bpo::value<long>()->default_value(logQ), "Log of ciphertext modulus.")
    ("logp", bpo::value<long>()->default_value(30), "Log of noise parameter.")
    ("logt", bpo::value<long>()->default_value(4), "Log of bootstrap noise parameter.")
    ("logq-boot", bpo::value<long>()->default_value(40), "Log of bootstrap ciphertext modulus.");

  desc.add_options()
    ("dusts", bpo::value<long>()->required(), "Number of dusts.")
    ("iterations", bpo::value<long>()->required(), "Number of iterations for mode seeking.")
    ("mode-kdegree", bpo::value<long>()->required(), "Log of kernel degree in mode seeking (\\Gamma_1).")
    ("label-kdegree", bpo::value<long>()->required(), "Log of kernel degree in point labeling (\\Gamma_2).")
    ("mode-invsteps", bpo::value<long>()->required(), "Number of iterations for inverse in mode seeking (\\zeta_1).")
    ("label-invsteps", bpo::value<long>()->required(), "Number of iterations for inverse in point labeling (\\zeta_2).")
    ("min-idx-t", bpo::value<long>()->required(), "Log of degree for minimum index algorithm (t).");

  return desc;
}
// clang-format on

int main(int argc, char *argv[]) {
  auto generic(generic_program_options());

  bpo::options_description cmdline("Run CKP19 mean shift clustering.");
  cmdline.add(generic);
  cmdline.add_options()(
      "config,c", bpo::value<std::string>(),
      "Configuration file for easy specification of cmd line arguments.")(
      "help,h", "Produce help message.");

  bpo::variables_map opts;
  bpo::store(bpo::command_line_parser(argc, argv).options(cmdline).run(), opts);

  if (opts.count("help") != 0) {
    std::cout << cmdline << std::endl;
    return 0;
  }

  if (opts.count("config") != 0) {
    std::string cpath(opts["config"].as<std::string>());
    std::ifstream fin(cpath.c_str());

    if (fin.fail()) {
      std::cerr << "Could not open configuration file at " << cpath << "\n";
      return 1;
    }

    bpo::store(bpo::parse_config_file(fin, generic), opts);
  }

  try {
    bpo::notify(opts);

    // Check if output file already exists
    if (fs::exists(opts["output"].as<std::string>())) {
      throw std::runtime_error("Output directory aready exists.");
    }
  } catch (const std::exception &ex) {
    std::cerr << ex.what() << std::endl;
    return 1;
  }

  try {
    benchmarkClustering(opts);
  } catch (const std::exception &ex) {
    std::cerr << ex.what() << "\nFatal error" << std::endl;
    return 1;
  }

  return 0;
}
