#pragma once

#include <common_types.hpp>
#include <json/json.hpp>
#include <util/debug.hpp>

namespace arb {

/// Description class for a regular spike source: a cell that generates
/// spikes with a fixed period over a given time interval.

struct ipss_cell_description {
    time_type start_time;
    time_type stop_time;

    // Every sample_delta we sample if we should emit a spike (in ms)
    double sample_delta;

    // vector of spike_rates starting at times
    // The vector needs at least a single entry
    std::vector<std::pair<time_type, double>> rates_per_time;
    bool interpolate;

    // rates_per_time: A vector of spike_rates each starting at the supplied time
    // The first time-rate pair should have a time before the start_time
    ipss_cell_description(time_type start_time, time_type stop_time, time_type sample_delta,
        std::vector<std::pair<time_type, double>> rates_per_time, bool interpolate = true) :
        start_time(start_time),
        stop_time(stop_time),
        sample_delta(sample_delta),
        rates_per_time(std::move(rates_per_time)),
        interpolate(interpolate)
    {}

    // Collect all the cell parameters from its json description
    ipss_cell_description(nlohmann::json const& cell_options) {
        start_time = cell_options["start_time"];
        stop_time = cell_options["stop_time"];
        sample_delta = cell_options["sample_delta"];
        interpolate = cell_options["interpolate"];

        auto times = cell_options["times"];
        auto rates = cell_options["rates"];

        EXPECTS(times.size() == rates.size());
        for (unsigned idx = 0; idx < times.size(); ++idx) {
            rates_per_time.push_back({ times[idx], rates[idx] });
        }
    }

};

} // namespace arb
