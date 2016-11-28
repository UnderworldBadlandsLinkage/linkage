import numpy as np
from scipy.interpolate import griddata
import underworld as uw


class LinkageModel(object):
    """
    A LinkageModel joins an Underworld and a Badlands model. Underworld models
    the domain in 3D and Badlands models the surface processes.

    LinkageModel hides the details of how the two models communicate so you can
    concentrate on building them.

    To define a linked model, instantiate a LinkageModel object and set, at
    least, the following members:

    mesh: the FeMesh that the model is defined over
    swarm: TODO
    material_index: an array of material indices. These will be changed when Badlands runs a tick.
    velocity_field: an Underworld velocity field. It must be defined over the same domain as the Badlands surface.

        self.velocityField is used to advect the Badlands surface, so it should
        be representative of the movement of the surface you are modelling.

        On each iteration, the material types of self.mesh will be altered
        according to Badland's processes. Particles which have transitioned
        from air to sediment will be given the self.sedimentIndex material
        type. Particles which have transitioned from sediment to air will be
        given the self.airIndex material type. You can override these if you
        wish.

        NOTE: the Badlands surface and Underworld mesh must be defined over the
        same coordinate system.


    badlands_model: the Badlands model which you would like to use. It must be initalised (XML and DEM loaded). It must be at time=0.
    update_function: a function which is called on each iteration. It must perform all of the per-iteration processing required for your Underworld model (usually, solving the Stokes equation and advecting the system.) It has the interface func(linkage, update_seconds)
    You can get the current time from linkage.t
    You should not run more than max_seconds or your checkpoints will no longer be synchronised.
    You must return the number of seconds that you *did* run (usually advector.get_max_dt())

        You might also use this to output interesting data to disk or to modify
        the behaviour of the model at a specific time.
        The update_function member can be changed at any time (say, if you're running in a Jupyter Notebook and want to adjust the model mid-run).

    checkpoint_function: a function which is called at the end of each checkpoint interval. You should write any relevant Underworld state to disk. Badlands state is written to disk automatically.
        The function is called with arguments 'linkage' and 'checkpoint_number'
        linkage: the linkage model object
        checkpoint_number: integer increasing by 1 per checkpoint. The initial state of the system is checkpointed as number 0, and subsequent states are 1..n.
        time_years: current time of the model (in years) - note that this happens AFTER advection. Usually you will write this to disk so output can be loaded synchronised with Badlands.

    air_material_indices
    sediment_material_indices

    You can then run the model by calling the run_for_years() function.

    Output will be saved as files on disk.
    """

    def __init__(self):
        # These MUST be set after you instantiate the object
        self.velocity_field = None
        self.material_index = None
        self.update_function = None
        self.badlands_model = None

        # You probably want to override some of these settings

        # Badlands and Underworld will write synchronised output this many
        # years apart
        self.checkpoint_interval = 10000

        # When Badlands deposits sediment at a location, its material type is
        # changed to deposited_material_index. If Badlands erodes material from
        # a location, the material type is changed to eroded_material_indexits
        # material type is changed to deposited_material_index. If Badlands
        # erodes material from a location, the material type is changed to
        # eroded_material_index. You can change these to anything you like.
        self.deposited_material_index = 103
        self.eroded_material_index = 104

        # --- You don't need to modify any settings below this line ---

        self.time_years = 0.  # Simulation time in years. We start at year 0.
        self._model_started = False  # Have we performed one-time initialisation yet?
        self._disp_inserted = False

        self.SECONDS_PER_YEAR = float(365 * 24 * 60 * 60)
        self._checkpoint_number = 0  # used to give each checkpoint a unique index
        self._next_checkpoint_years = None

        # Set this variable to stop Badlands changes from modifying the
        # Underworld material types. This is usually used for open-loop testing
        # to ensure that BL and UW do not desync after many timesteps.
        self.disable_material_changes = False

    def run_for_years(self, years):
        """
        Run the model for a number of years.

        If you pass 0 for the number of years, all of the initial states will
        be set up but the model will not be run.
        """
        if not self._model_started:
            self._startup()

        end_years = self.time_years + years
        while self.time_years < end_years:
            # What's the longest we can run before we have to write a checkpoint?
            max_seconds = (self._next_checkpoint_years - self.time_years) * self.SECONDS_PER_YEAR

            # Ask the Underworld model to update
            dt_seconds = self.update_function(self, max_seconds)
            assert dt_seconds <= max_seconds, "Maximum dt (seconds) for the update function was %s, but it ran for more than that (%s seconds)" % (max_seconds, dt_seconds)

            # Do we need to write a checkpoint later?
            # TODO: make sure floating point imperfections don't desync the seconds/years counters on both sides, especially around writing checkpoints
            write_checkpoint = False
            if dt_seconds == max_seconds:
                write_checkpoint = True

            dt_years = dt_seconds / self.SECONDS_PER_YEAR

            # Advect the Badlands interface surface
            print 'advect %s seconds' % dt_seconds
            self._surface_advector.integrate(dt_seconds)

            # Synchronise the velocity field across nodes
            print "TODO ENSURE SINGLE CORE"
            # assert uw.size() == 1, "MPI is not supported yet"
            # TODO implement sync across MPI nodes

            # the entire velocity vector on each particle in METERS PER SECOND
            tracer_velocity_mps = self.velocity_field.evaluate(self._surface_tracers)

            ### INTERFACE PART 1: UW->BL
            # Use the tracer vertical velocities to deform the Badlands TIN
            # convert from meters per second to meters displacement over the whole iteration
            tracer_disp = tracer_velocity_mps * self.SECONDS_PER_YEAR * dt_years
            self._inject_badlands_displacement(self.time_years, dt_years, tracer_disp)

            # Run the Badlands model to the same time point
            self.badlands_model.run_to_time(self.time_years + dt_years)

            ### INTERFACE PART 2: BL->UW
            self._update_material_types()

            if write_checkpoint:
                # checkpoint fields and swarm
                self.checkpoint_function(self, self._checkpoint_number, self.time_years)
                self._checkpoint_number += 1
                self._next_checkpoint_years += self.checkpoint_interval

            # Advance time
            self.time_years += dt_years

    @staticmethod
    def generate_flat_dem(minCoord, maxCoord, resolution, elevation):
        """
        Generate a flat DEM. This can be used as the initial Badlands state.

        minCoord: tuple of (X, Y, Z) coordinates defining the minimum bounds of
                  the DEM. Only the X and Y coordinates are used.
        maxCoord: tuple of (X, Y, Z) coordinates defining the maximum bounds of
                  the DEM. Only the X and Y coordinates are used.
        resolution: tuple of (X resolution, Y resolution) defining the number
                    of points dividing the DEM along the X and Y axes.
        elevation: the Z parameter that each point is created at

        For your convenience, minCoord and maxCoord are designed to have the
        same formatting as the Underworld FeMesh_Cartesian minCoord and
        maxCoord parameters.

        Note that the initial elevation of the Badlands surface should coincide
        with the material transition in Underworld.
        """
        items = []
        # FIXME: there should be a fast numpy way to do this
        for y in np.linspace(minCoord[1], maxCoord[1], resolution[1]):
            for x in np.linspace(minCoord[0], maxCoord[0], resolution[0]):
                items.append([x, y, elevation])

        # NOTE: Badlands uses the difference in X coord of the first two points to determine the resolution.
        # This is something we should fix.
        # This is why we loop in y/x order instead of x/y order.
        return np.array(items)

    def _startup(self):
        """
        Perform one-time initialisation of the models.

        We load everything, transfer the initial surface elevation from
        Badlands to Underworld, then write the initial state to disk.
        """

        assert not self._model_started

        # Make sure UW and BL are operating over the same XY domain
        rg = self.badlands_model.recGrid
        bl_xy = (rg.rectX.min(), rg.rectX.max(), rg.rectY.min(), rg.rectY.max())
        uw_xy = (self.mesh.minCoord[0], self.mesh.maxCoord[0], self.mesh.minCoord[1], self.mesh.maxCoord[1])
        assert bl_xy == uw_xy, "Badlands and Underworld must operate over the same domain (Badlands has %s, but Underworld has %s)" % (bl_xy, uw_xy)

        self.badlands_model.input.disp3d = True  # enable 3D displacements
        self.badlands_model.input.region = 0  # TODO: check what this does

        # Override the checkpoint/display interval in the Badlands model to
        # ensure BL and UW are synced
        self.badlands_model.input.tDisplay = self.checkpoint_interval

        # create tracers used to track movement of the Badlands surface
        bl_tracers = uw.swarm.Swarm(self.mesh)
        rg = self.badlands_model.recGrid
        dem_coords = np.column_stack((rg.rectX, rg.rectY, rg.rectZ))
        bl_tracers.add_particles_with_coordinates(dem_coords)

        self._surface_advector = uw.systems.SwarmAdvector(swarm=bl_tracers, velocityField=self.velocity_field, order=2)
        self._surface_tracers = bl_tracers

        # Transfer the initial DEM state to Underworld
        self._update_material_types()

        self._next_checkpoint_years = self.checkpoint_interval

        # Perform an initial checkpoint
        self.checkpoint_function(self, self._checkpoint_number, self.time_years)
        self._checkpoint_number += 1

        self._model_started = True

    def _determine_particle_state(self):
        # Given Badlands' mesh, determine if each particle in 'volume' is above
        # (False) or below (True) it.

        # To do this, for each X/Y pair in 'volume', we interpolate its Z value
        # relative to the mesh in blModel. Then, if the interpolated Z is
        # greater than the supplied Z (i.e. Badlands mesh is above particle
        # elevation) it's sediment (True). Else, it's air (False).

        known_xy = self.badlands_model.recGrid.tinMesh['vertices']  # points that we have known elevation for
        known_z = self.badlands_model.elevation  # elevation for those points

        volume = self.swarm.particleCoordinates.data

        interpolate_xy = volume[:, [0, 1]]
        # linear interpolation should be plenty as we're running Badlands at
        # higher resolution than Underworld
        interpolate_z = griddata(points=known_xy, values=known_z, xi=interpolate_xy, method='linear')

        # True for sediment, False for air
        flags = volume[:, 2] < interpolate_z

        return flags

    def load_badlands_dem_file(self, filename):
        self.badlands_model.build_mesh(filename, verbose=False)

    def load_badlands_dem_array(self, array):
        # for now, write it out to a temp file and load that into badlands
        np.savetxt('/tmp/dem.csv', array)
        self.load_badlands_dem_file('/tmp/dem.csv')

    def _inject_badlands_displacement(self, time, dt, disp):#tracers, disp, blModel, time, dt, n, min_coord, max_coord, display_interval=None):
        """
        Takes a plane of tracer points and their DISPLACEMENTS in 3D over time
        period dt. Injects it into Badlands as 3D tectonic movement.
        """

        # TODO: what does this do? Can it be removed?
        self.badlands_model.force.merge3d = self.badlands_model.input.Afactor * self.badlands_model.recGrid.resEdges * 0.5

        # The Badlands 3D interpolation map is the displacement of each DEM
        # node at the end of the time period relative to its starting position.
        # If you start a new displacement file, it is treated as starting at
        # the DEM starting points (and interpolated onto the TIN as it was at
        # that tNow).

        # kludge; don't keep adding new entries
        if self._disp_inserted:
            self.badlands_model.force.T_disp[0, 0] = time
            self.badlands_model.force.T_disp[0, 1] = (time + dt)
        else:
            self.badlands_model.force.T_disp = np.vstack(([time, time + dt], self.badlands_model.force.T_disp))
            self._disp_inserted = True

        self.badlands_model.force.injected_disps = disp

    def _update_material_types(self):
        if self.disable_material_changes:
            return

        # What do the materials (in air/sediment terms) look like now?
        material_flags = self._determine_particle_state()

        # If any materials changed state, update the Underworld material types
        # TODO(perf): vectorise
        for index, material in enumerate(self.material_index.data):
            # convert air to sediment
            if int(material) in self.air_material_indices and material_flags[index]:
                self.material_index.data[index] = self.deposited_material_index

            # convert sediment to air
            if int(material) in self.sediment_material_indices and not material_flags[index]:
                self.material_index.data[index] = self.eroded_material_index
