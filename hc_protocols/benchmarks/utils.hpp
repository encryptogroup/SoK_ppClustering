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

#include <cnpy.h>
#include <sys/types.h>
#include <unistd.h>

#include <chrono>
#include <iostream>
#include <nlohmann/json.hpp>

#define BENCHMARK(bench, f, ...)                      \
  ({                                                  \
    reset_counters();                                 \
    aby_reset_counters();                             \
    TimePoint start;                                  \
    f(__VA_ARGS__);                                   \
    TimePoint end;                                    \
    json res;                                         \
    res["time"] = end - start;                        \
    res["setup_bytes_sent"] = bytes_sent();           \
    res["setup_bytes_received"] = bytes_received();   \
    res["aby_bytes_sent"] = aby_bytes_sent();         \
    res["aby_bytes_received"] = aby_bytes_received(); \
    bench.push_back(res);                             \
    std::cout << "--- " << #f << " ---\n";            \
    print_json(res);                                  \
  })

using json = nlohmann::json;

class TimePoint {
 public:
  using ChronoTimePoint = std::chrono::high_resolution_clock::time_point;
  using TimeUnit = std::chrono::duration<double, std::milli>;

  TimePoint() : time_(ChronoTimePoint::clock::now()) {}

  double operator-(const TimePoint& rhs) {
    return std::chrono::duration_cast<TimeUnit>(time_ - rhs.time_).count();
  }

 private:
  ChronoTimePoint time_;
};

bool save_json(const json& data, const std::string& fpath) {
  std::ofstream fout;
  fout.open(fpath, std::fstream::out);
  if (!fout.is_open()) {
    std::cerr << "Could not open save file at " << fpath << "\n";
    return false;
  }

  fout << data;
  fout.close();

  std::cout << "Saved data in " << fpath << std::endl;

  return true;
}

void print_json(const json& data) {
  for (json::const_iterator it = data.begin(); it != data.end(); ++it) {
    std::cout << it.key() << " : " << it.value() << "\n";
  }
  std::cout << std::endl;
}

std::vector<point<int64_t>> load_dataset(std::string dataset) {
  auto npy_arr = cnpy::npy_load(dataset);
  auto npy_data = npy_arr.data<double>();

  if (npy_arr.shape.size() != 2) {
    throw std::logic_error("Expected 2D matrix for dataset.");
  }

  size_t rows = npy_arr.shape[0];
  size_t cols = npy_arr.shape[1];

  std::vector<double> vmin(cols, std::numeric_limits<double>::infinity());

  for (size_t i = 0; i < rows; ++i) {
    for (size_t j = 0; j < cols; ++j) {
      vmin[j] = std::min(vmin[j], npy_data[i * cols + j]);
    }
  }

  std::vector<point<int64_t>> res(rows, point<int64_t>(cols));

  for (size_t i = 0; i < rows; ++i) {
    for (size_t j = 0; j < cols; ++j) {
      double v = npy_data[i * cols + j] - vmin[j];
      v *= 1000;
      res[i][j] = static_cast<int64_t>(v);
    }
  }

  return res;
}

// Reference: https://gist.github.com/k3vur/4169316
int64_t getProcStatus(const std::string& key) {
  int64_t value = 0;

  std::string temp = "/proc/" + std::to_string(getpid()) + "/status";
  const char* filename = temp.c_str();

  std::ifstream procfile(filename);
  std::string word;
  while (procfile.good()) {
    procfile >> word;
    if (word == key) {
      procfile >> value;
      break;
    }

    // Skip to end of line.
    procfile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
  }

  if (procfile.fail()) {
    return -1;
  }

  return value;
}
