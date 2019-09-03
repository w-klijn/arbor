#include "../gtest.h"

#include <arbor/common_types.hpp>

#include "epoch.hpp"
#include "fvm_lowered_cell.hpp"
#include "mc_cell_group.hpp"
#include "util/rangeutil.hpp"

#include "common.hpp"
#include "../common_cells.hpp"
#include "../simple_recipes.hpp"

using namespace arb;

namespace {
    execution_context context;

    fvm_lowered_cell_ptr lowered_cell() {
        return make_fvm_lowered_cell(backend_kind::multicore, context);
    }

    cable_cell make_cell() {
        auto c = make_cell_ball_and_stick();

        c.add_detector({0, 0}, 0);
        c.segment(1)->set_compartments(101);

        return c;
    }
}

ACCESS_BIND(
    std::vector<cell_member_type> mc_cell_group::*,
    private_spike_sources_ptr,
    &mc_cell_group::spike_sources_)

TEST(mc_cell_group, get_kind) {
    mc_cell_group group{{0}, cable1d_recipe(make_cell()), lowered_cell()};

    EXPECT_EQ(cell_kind::cable, group.get_cell_kind());
}

TEST(mc_cell_group, test) {
    auto rec = cable1d_recipe(make_cell());
    rec.nernst_ion("na");
    rec.nernst_ion("ca");
    rec.nernst_ion("k");

    mc_cell_group group{{0}, rec, lowered_cell()};
    group.advance(epoch(0, 50), 0.01, {});

    // Model is expected to generate 4 spikes as a result of the
    // fixed stimulus over 50 ms.
    EXPECT_EQ(4u, group.spikes().size());
}

TEST(mc_cell_group, sources) {
    // Make twenty cells, with an extra detector on gids 0, 3 and 17
    // to make things more interesting.
    std::vector<cable_cell> cells;

    for (int i=0; i<20; ++i) {
        cells.push_back(make_cell());
        if (i==0 || i==3 || i==17) {
            cells.back().add_detector({1, 0.3}, 2.3);
        }

        EXPECT_EQ(1u + (i==0 || i==3 || i==17), cells.back().detectors().size());
    }

    std::vector<cell_gid_type> gids = {3u, 4u, 10u, 16u, 17u, 18u};
    auto rec = cable1d_recipe(cells);
    rec.nernst_ion("na");
    rec.nernst_ion("ca");
    rec.nernst_ion("k");

    mc_cell_group group{gids, rec, lowered_cell()};

    // Expect group sources to be lexicographically sorted by source id
    // with gids in cell group's range and indices starting from zero.

    const auto& sources = group.*private_spike_sources_ptr;
    for (unsigned j = 0; j<sources.size(); ++j) {
        auto id = sources[j];
        if (j==0) {
            EXPECT_EQ(id.gid, gids[0]);
            EXPECT_EQ(id.index, 0u);
        }
        else {
            auto prev = sources[j-1];
            EXPECT_GT(id, prev);
            EXPECT_EQ(id.index, id.gid==prev.gid? prev.index+1: 0u);
        }
    }
}

