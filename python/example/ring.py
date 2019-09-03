import sys
import arbor

class ring_recipe (arbor.recipe):

    def __init__(self, n=4):
        # The base C++ class constructor must be called first, to ensure that
        # all memory in the C++ class is initialized correctly.
        arbor.recipe.__init__(self)
        self.ncells = n
        self.params = arbor.cell_parameters()

    # The num_cells method that returns the total number of cells in the model
    # must be implemented.
    def num_cells(self):
        return self.ncells

    # The cell_description method returns a cell
    def cell_description(self, gid):
        return arbor.make_cable_cell(gid, self.params)

    def num_targets(self, gid):
        return 1

    def num_sources(self, gid):
        return 1

    # The kind method returns the type of cell with gid.
    # Note: this must agree with the type returned by cell_description.
    def cell_kind(self, gid):
        return arbor.cell_kind.cable

    # Make a ring network
    def connections_on(self, gid):
        src = (gid-1)%self.ncells
        w = 0.01
        d = 10
        return [arbor.connection(arbor.cell_member(src,0), arbor.cell_member(gid,0), w, d)]

    # Attach a generator to the first cell in the ring.
    def event_generators(self, gid):
        if gid==0:
            sched = arbor.explicit_schedule([1])
            return [arbor.event_generator(arbor.cell_member(0,0), 0.1, sched)]
        return []


context = arbor.context(threads=4, gpu_id=None)
print(context)

meters = arbor.meter_manager()
meters.start(context)

recipe = ring_recipe(100)
print(f'{recipe}')

meters.checkpoint('recipe-create', context)

decomp = arbor.partition_load_balance(recipe, context)
print(f'{decomp}')

hint = arbor.partition_hint()
hint.prefer_gpu = True
hint.gpu_group_size = 1000
print(f'{hint}')

hints = dict([(arbor.cell_kind.cable, hint)])
decomp = arbor.partition_load_balance(recipe, context, hints)
print(f'{decomp}')

meters.checkpoint('load-balance', context)

sim = arbor.simulation(recipe, decomp, context)
print(f'{sim} finished')

meters.checkpoint('simulation-init', context)

recorder = arbor.attach_spike_recorder(sim)

sim.run(1000)

meters.checkpoint('simulation-run', context)

print(f'{arbor.meter_report(meters, context)}')

for s in recorder.spikes:
    print(s)
