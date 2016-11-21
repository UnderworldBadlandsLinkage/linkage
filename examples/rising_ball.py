# Rising Ball demo
# We start with a flat surface and a low-density sphere underground. As the
# sphere rises due to the density difference, the surface is deformed.

from linkage import LinkageModel


class RisingBall(LinkageModel):
    def badlands_init(self):
        # describe how the surface looks at the start
        dem = self.generate_flat_dem(x_range=(0, 100e3), y_range=(0, 100e3), x_resolution=180, y_resolution=180)
        self.badlands.load_dem(array=dem)

        # self.badlands.load_dem('nodes.csv')

    def underworld_init(self):
        # describe how the model initially looks in 3D
        # Set size and position of cube
        radius = 10e3
        centre = (50e3, 50e3, -20e3)  # make sure the cube starts underground. We can't model non-flat starting surfaces in Badlands yet.


        mesh = uw.mesh.FeMesh_Cartesian( elementType = ("Q1/dQ0"), 
                                         elementRes  = [UNDERWORLD_RESOLUTION] * 3, 
                                         minCoord    = MIN_COORD, 
                                         maxCoord    = MAX_COORD)

        velocityField    = uw.mesh.MeshVariable( mesh=mesh,         nodeDofCount=mesh.dim )
        pressureField    = uw.mesh.MeshVariable( mesh=mesh.subMesh, nodeDofCount=1 )

        velocityField.data[:] = [0.,0.,0.]
        pressureField.data[:] = 0.

        # Create the swarm, material index variable and swarm advector

        todo: how much of this is plumbing and should be hidden?
        swarm = uw.swarm.Swarm(mesh=mesh)
        materialIndex = swarm.add_variable(dataType="int", count=1)

        swarmLayout = uw.swarm.layouts.GlobalSpaceFillerLayout(swarm=swarm, particlesPerCell=20)
        swarm.populate_using_layout(layout=swarmLayout)

        advector = uw.systems.SwarmAdvector(swarm=swarm, velocityField=velocityField, order=2)

        # Let's initialise the 'materialVariable' data to represent two different materials. 
airIndex = 0
heavyIndex = 1
lightIndex = 2
sedimentIndex = 3
erodedIndex = 4

for index, coord in enumerate(swarm.particleCoordinates.data):
    offset = coord - centre
    if (offset[0]**2 + offset[1]**2 + offset[2]**2 < radius**2):
        materialIndex.data[index] = lightIndex
    elif coord[2] > AIR_ELEVATION:
        materialIndex.data[index] = airIndex
    else:
        materialIndex.data[index] = heavyIndex

# Create tracers points for interpolating velocity, used to comm. with BL
bl_tracers = uw.swarm.Swarm(mesh)

# grid between x,y:[5,15] regularly spaced 10 points in each direction
# TODO: pull points from badlands and use them to determine velocity
plane = np.vstack([np.mgrid[MIN_COORD[0]:MAX_COORD[0]:10j, MIN_COORD[0]:MAX_COORD[1]:10j]]).reshape(2, -1).T
# TODO: we need to know exactly how the 3d displacement mesh in badlands works to know how best to do this

# initial height just below 'air'
z = np.full(shape=(len(plane), 1), fill_value=AIR_ELEVATION)
some = np.c_[plane, z]  # way to join to column vectors together in numpy - efficiency unknown

# create swarm
bl_tracers.add_particles_with_coordinates(some)
advector2 = uw.systems.SwarmAdvector(swarm=bl_tracers, velocityField=velocityField, order=2)


# Set viscosities and densities of the model.
viscosityMapFn = 1e19

# Here we set a density of '0.' for the lightMaterial, and '1.' for the heavymaterial.
mappingDictDensity = {airIndex: 3300.0,
                      lightIndex: 3240.0,
                      heavyIndex: 3300.0,
                      sedimentIndex: 3240.0,
                      erodedIndex: 3300.0}
densityFn = fn.branching.map(fn_key=materialIndex, mapping=mappingDictDensity)

# And the final buoyancy force function.
buoyancyFn = densityFn * 9.8 * [0.0, 0.0, -1.0]

# wall velocity boundary conditions - free slip on all walls

iWalls = mesh.specialSets["MinI_VertexSet"] + mesh.specialSets["MaxI_VertexSet"]
jWalls = mesh.specialSets["MinJ_VertexSet"] + mesh.specialSets["MaxJ_VertexSet"]
kWalls = mesh.specialSets["MinK_VertexSet"] + mesh.specialSets["MaxK_VertexSet"]

velocityBC = uw.conditions.DirichletCondition(variable=velocityField,
                                              indexSetsPerDof=(iWalls, jWalls, kWalls))

# combine all the above into Stokes system and get solver
stokesPIC = uw.systems.Stokes(velocityField = velocityField,
                              pressureField = pressureField,
                              swarm         = swarm,
                              conditions    = [velocityBC,],
                              fn_viscosity  = viscosityMapFn,
                              fn_bodyforce  = buoyancyFn)
solver = uw.systems.Solver(stokesPIC)

# some analytics

vdotv = fn.math.dot(velocityField, velocityField)
v2sum_integral = uw.utils.Integral(mesh=mesh, fn=vdotv)
volume_integral = uw.utils.Integral(mesh=mesh, fn=1.)


    def on_tick(self, years):
        pass  # no time-dependent behaviour for this model
        # you might want to write custom output or metrics from here


model = RisingBall()
model.checkpoint_interval = 10000
model.run_for_years(100000)
