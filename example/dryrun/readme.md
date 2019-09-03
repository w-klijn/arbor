# Dryrun Example

A miniapp that demonstrates how to use dry-run mode on a simple network
duplicated across `num_ranks`.

It uses the `arb::tile` to build a network of `num_cells_per_rank` cells of type
`arb::cable_cell`. The network is translated over `num_ranks` domains using
`arb::symmetric_recipe`.

Example:
```
num_cells_per_rank = 4;
num_ranks = 2;

gids = {0, ..., 7}

tile connections:
    dest <- src:

    0 <- 1;
    1 <- 3;
    2 <- 5;
    3 <- 7;

symmetric_recipe inferred connections:
    translate tile connections for second doimain:
    dest <- src:

    0 <- 1;
    1 <- 3;
    2 <- 5;
    3 <- 7;

    4 <- 5;
    5 <- 7;
    6 <- 1;
    7 <- 3;

```

The model of the *tile* can be configured using a json configuration file:

```
./bench.exe params.json
```

An example parameter file for a dry-run is:
```
{
    "name": "dry run test",
    "dry-run": true,
    "num-cells-per-rank": 1000,
    "num-ranks": 2,
    "duration": 100,
    "min-delay": 1,
    "depth": 5,
    "branch-probs": [1.0, 0.5],
    "compartments": [100, 2]
}

```

The parameter file for the equivalent MPI run is:
```
{
    "name": "MPI test",
    "dry-run": false,
    "num-cells-per-rank": 1000,
    "duration": 100,
    "min-delay": 1,
    "depth": 5,
    "branch-probs": [1.0, 0.5],
    "compartments": [100, 2]
}

```
These 2 files should provide exactly the same spike.gdf files.


The parameters in the file:
  * `name`: a string with a name for the benchmark.
  * `dry-run`: a bool indicating whether or not to use dry-run mode.
    if false, use MPI if available, and local distribution if not.
  * `num-ranks`: the number of domains to simulate. Only for dry-run
    mode.
  * `num-cells-per-rank`: the number of cells on a single tile.
    The total number of cells in the model = num-cells-per-rank *
    num-ranks.
  * `duration`: the length of the simulated time interval, in ms.
  * `fan-in`: the number of incoming connections on each cell.
  * `min-delay`: the minimum delay of the network.
  * `spike-frequency`: frequency of the independent Poisson processes that
    generate spikes for each cell.
  * `realtime-ratio`: the ratio between time taken to advance a single cell in
    the simulation and the simulated time. For example, a value of 1 indicates
    that the cell is simulated in real time, while a value of 0.1 indicates
    that 10s can be simulated in a single second.

The network is randomly connected with no self-connections and `fan-in`
incoming connections on each cell, with every connection having delay of
`min-delay`.