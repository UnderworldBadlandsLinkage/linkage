# "No Erosion Clock"
#
# This model tests two aspects of the linkage:
# * the coordinate systems between UW and BL are synchronised
# * the scaling between the two is the same
#
# Normally, you run this model for a long period (at least 20 timesteps) and make sure both sides are synchronised.


import os

import underworld as uw
from underworld import function as fn

from pyBadlands.model import Model as BadlandsModel

from linkagemodel.linkage import LinkageModel

# parameters for our model composition
MIN_COORD = (0., 0., -80e3)
MAX_COORD = (100e3, 70e3, 20e3)  # this crashes at DEM load? do we have X/Y swapped?
# MAX_COORD = (100e3, 100e3, 20e3)
INITIAL_AIR_ELEVATION = 0.0  # the height at which we transition from sediment to air

# Put a rectangular prism underground
cubeSize = (10e3, 20e3, 20e3)
centre = (50e3, 25e3, -30e3)

linkage = LinkageModel()


### SET UP THE BADLANDS MODEL
badlands_model = BadlandsModel()
badlands_model.load_xml('no_erosion_clock.xml')
linkage.badlands_model = badlands_model

# IMPORTANT: disable BL->UW transfers
# Because erosion is disabled in the Badlands config, this completely uncouples
# the two models. They independently advect the surface. In this way, we can
# ensure that they advect identically.
linkage.disable_material_changes = True

# Load a flat DEM
bl_dem = linkage.generate_flat_dem(minCoord=MIN_COORD, maxCoord=MAX_COORD, resolution=(180, 180), elevation=INITIAL_AIR_ELEVATION)
linkage.load_badlands_dem_array(bl_dem)


### SET UP THE UNDERWORLD MODEL

# All output will go to the 'uwout' directory, which we will create
uw_output_path = 'uwout'
try:
    os.mkdir(uw_output_path)
except OSError:
    # probably already exists
    pass

# Underworld models normally run at a much lower resolution than
# Badlands models in order to keep the computation time reasonable.
UNDERWORLD_RESOLUTION = 20

# This is the mesh whose material types will be changed by Badlands
mesh = uw.mesh.FeMesh_Cartesian(elementType=("Q1/dQ0"),
                                elementRes =[UNDERWORLD_RESOLUTION] * 3,
                                minCoord   =MIN_COORD,
                                maxCoord   =MAX_COORD)

# We want to track velocity and pressure.
velocityField = uw.mesh.MeshVariable(mesh=mesh, nodeDofCount=mesh.dim)
pressureField = uw.mesh.MeshVariable(mesh=mesh.subMesh, nodeDofCount=1)

# Set initial states
velocityField.data[:] = [0., 0., 0.]
pressureField.data[:] = 0.

# Create the swarm, material index variable and swarm advector
swarm = uw.swarm.Swarm(mesh=mesh)
materialIndex = swarm.add_variable(dataType="int", count=1)

swarmLayout = uw.swarm.layouts.GlobalSpaceFillerLayout(swarm=swarm, particlesPerCell=20)
swarm.populate_using_layout(layout=swarmLayout)

advector = uw.systems.SwarmAdvector(swarm=swarm, velocityField=velocityField, order=2)

airIndex = 0
heavyIndex = 1
lightIndex = 2
sedimentIndex = 3
erodedIndex = 4

# configure material mappings
linkage.deposited_material_index = sedimentIndex
linkage.eroded_material_index = erodedIndex
linkage.air_material_indices = [airIndex, erodedIndex]
linkage.sediment_material_indices = [heavyIndex, lightIndex, sedimentIndex]

# Define a rectangular prism
for index, coord in enumerate(swarm.particleCoordinates.data):
    offset = coord - centre
    if offset[0]**2 < cubeSize[0]**2 and offset[1]**2 < cubeSize[1]**2 and offset[2]**2 < cubeSize[2]**2:
        materialIndex.data[index] = lightIndex
    elif coord[2] > INITIAL_AIR_ELEVATION:
        materialIndex.data[index] = airIndex
    else:
        materialIndex.data[index] = heavyIndex

# Set viscosities and densities of the model.
viscosityMapFn = 1e19

# TODO update below comment
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
stokesPIC = uw.systems.Stokes(velocityField=velocityField,
                              pressureField=pressureField,
                              swarm        =swarm,
                              conditions   =[velocityBC, ],
                              fn_viscosity =viscosityMapFn,
                              fn_bodyforce =buoyancyFn)
solver = uw.systems.Solver(stokesPIC)

# Linkage needs to know about the mesh and the swarm the model is defined over
# TODO: can you eliminate one of these variables?
linkage.mesh = mesh
linkage.swarm = swarm
# This velocity field will be used to deform the Badlands surface
linkage.velocity_field = velocityField
# This array will store the updated material types from Badlands
linkage.material_index = materialIndex


def update_function(linkage, max_seconds):
    # Get solution for initial configuration.
    solver.solve()

    # Determine the maximum possible timestep for the advection system.
    dtmax_seconds = advector.get_max_dt()
    dt_seconds = min(max_seconds, dtmax_seconds)

    advector.integrate(dt_seconds)

    return dt_seconds
linkage.update_function = update_function


def checkpoint_function(linkage, checkpoint_number, time_years):
    mH = mesh.save(os.path.join(uw_output_path, "mesh.h5"))

    file_prefix = os.path.join(uw_output_path, 'velocity-%s' % checkpoint_number)
    handle = velocityField.save('%s.h5' % file_prefix)
    velocityField.xdmf('%s.xdmf' % file_prefix, handle, 'velocity', mH, 'mesh', modeltime=time_years)

    file_prefix = os.path.join(uw_output_path, 'pressure-%s' % checkpoint_number)
    handle = pressureField.save('%s.h5' % file_prefix)
    pressureField.xdmf('%s.xdmf' % file_prefix, handle, 'pressure', mH, 'mesh', modeltime=time_years)

    sH = swarm.save(os.path.join(uw_output_path, 'swarm-%s.h5' % checkpoint_number))

    file_prefix = os.path.join(uw_output_path, 'material-%s' % checkpoint_number)
    handle = materialIndex.save('%s.h5' % file_prefix)
    materialIndex.xdmf('%s.xdmf' % file_prefix, handle, 'material', sH, 'swarm', modeltime=time_years)
linkage.checkpoint_function = checkpoint_function


### RUN THE MODEL
linkage.checkpoint_interval = 10000
linkage.run_for_years(200000)
