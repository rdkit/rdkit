//
// Copyright (C) David Cosgrove 2025.
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
// This implements a very basic thread-safe progress bar, based on code
// here
// https://www.cppstories.com/2020/02/inidicators.html/
// which is written by the author of https://github.com/p-ranav/indicators
// which is released under an MIT license.

#ifndef PROGRESSBAR_H
#define PROGRESSBAR_H

#include <atomic>
#include <chrono>
#include <ctime>
#include <mutex>
#include <iostream>

namespace RDKit {
class ProgressBar {
 public:
  ProgressBar() : d_startTime(std::chrono::system_clock::now()) {}
  ProgressBar(unsigned int width, std::uint64_t numToDo)
      : d_bar_width(width),
        d_numToDo(numToDo),
        d_startTime(std::chrono::system_clock::now()) {}
  ProgressBar(const ProgressBar &) = delete;
  ProgressBar &operator=(const ProgressBar &) = delete;
  ProgressBar(ProgressBar &&) = delete;
  ProgressBar &operator=(ProgressBar &&) = delete;
  ~ProgressBar() = default;

  void update(float value, std::ostream &os = std::cout) {
    set_progress(value);
    write_progress(os);
  }
  void increment(std::uint64_t incr = 1, std::ostream &os = std::cout) {
    d_done += incr;
    float newProgress =
        100.0 * static_cast<float>(d_done) / static_cast<float>(d_numToDo);
    update(newProgress, os);
  }

  void set_progress(float value) {
    std::unique_lock lock{d_mutex};
    d_progress = value;
  }

  void write_progress(std::ostream &os = std::cout) {
    std::unique_lock lock{d_mutex};
    if (d_progress > 100.0f) {
      return;
    }
    os << "\r" << std::flush;
    os << "[";
    const auto completed = static_cast<size_t>(
        d_progress * static_cast<float>(d_bar_width) / 100.0);
    for (size_t i = 0; i < d_bar_width; ++i) {
      if (i <= completed)
        os << "*";
      else
        os << " ";
    }

    // End bar
    os << "]";

    auto currTime = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed = currTime - d_startTime;
    std::chrono::duration<double> totToDo = elapsed / (d_progress / 100.0);
    auto projFinish =
        d_startTime +
        std::chrono::duration_cast<std::chrono::system_clock::duration>(
            totToDo);

    // Convert to time_t for localtime
    std::time_t finish_time_t =
        std::chrono::system_clock::to_time_t(projFinish);

    // Convert to local time structure
    std::tm local_tm = *std::localtime(&finish_time_t);
    std::string fmt = "%H:%M %a";
    static constexpr std::int64_t weekSecs = 7 * 24 * 3600;
    if (std::chrono::duration_cast<std::chrono::seconds>(totToDo).count() >
        weekSecs) {
      fmt = "%d %m %y";
    }
    // Write progress percentage and actual values
    os << " " << std::min(static_cast<size_t>(d_progress), size_t(100)) << "%"
       << " " << d_done << "/" << d_numToDo << " "
       << std::put_time(&local_tm, fmt.c_str()) << std::flush;
  }

 private:
  std::mutex d_mutex;
  float d_progress{0.0};
  unsigned int d_bar_width{60};
  std::uint64_t d_numToDo{0};
  std::atomic<std::uint64_t> d_done{0};
  std::chrono::time_point<std::chrono::system_clock> d_startTime;
};
}  // namespace RDKit

#endif  // PROGRESSBAR_H
