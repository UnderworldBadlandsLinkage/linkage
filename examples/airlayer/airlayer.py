# Demonstrate higher Underworld resolution near Badlands surface and use of
# low-density air

import os

from pyBadlands.model import Model as BadlandsModel
import underworld as uw
from underworld import function as fn
from linkagemodel.linkage import LinkageModel

# Create mesh and associated mesh variables
MIN_COORD = (0., 0., -80e3)
MAX_COORD = (100e3, 100e3, 20e3)
UNDERWORLD_RESOLUTION = [20, 20, 32]
BADLANDS_RESOLUTION = (180, 180)
# BADLANDS_RESOLUTION = (120, 120)
AIR_ELEVATION = 0.0  # the height at which we transition from sediment to air

# Set size and position of sphere
radius = 10e3
centre = (50e3, 50e3, -20e3)

# Let's initialise the 'materialVariable' data to represent two different materials. 
airIndex = 0
heavyIndex = 1
lightIndex = 2
sedimentIndex = 3
erodedIndex = 4

### SET UP THE LINKAGE
linkage = LinkageModel()


### SET UP THE UNDERWORLD MODEL

# All output will go to the 'uwout' directory, which we will create
uw_output_path = 'uwout'
try:
    os.mkdir(uw_output_path)
except OSError:
    # probably already exists
    pass

mesh = uw.mesh.FeMesh_Cartesian( elementType = ("Q1/dQ0"), 
                                 elementRes  = UNDERWORLD_RESOLUTION, 
                                 minCoord    = MIN_COORD, 
                                 maxCoord    = MAX_COORD)

velocityField    = uw.mesh.MeshVariable( mesh=mesh,         nodeDofCount=mesh.dim )
pressureField    = uw.mesh.MeshVariable( mesh=mesh.subMesh, nodeDofCount=1 )

velocityField.data[:] = [0.,0.,0.]
pressureField.data[:] = 0.
linkage.velocity_field = velocityField


# MEXICO_2
####  HOW TO DEFORM THE MESH FOR FOCUSSED RESOLUTION ####

# IAN FIXME this breaks add_particles; no idea why. Only 20% of the particles actually get added.

# we use the poisson equation, varying diffusivity to solve for
# new distribution of points, in the vertical axis

'''
yField    = uw.mesh.MeshVariable( mesh=mesh, nodeDofCount=1 )
coord = fn.coord()
yBCs = uw.conditions.DirichletCondition( variable    = yField, 
                                         indexSetsPerDof = ( mesh.specialSets["MinJ_VertexSet"]+mesh.specialSets["MaxJ_VertexSet"] ) )

yField.data[mesh.specialSets["MinJ_VertexSet"].data] = MIN_COORD[2]
yField.data[mesh.specialSets["MaxJ_VertexSet"].data] = MAX_COORD[2]

intensityFn = 1.0+6.0*uw.function.math.exp(-10.0*uw.function.math.pow((coord[1]-0.4),2.0))

nodeDiffuse = uw.systems.SteadyStateHeat(yField, 
                                     fn_diffusivity=intensityFn,
                                     fn_heating=0.0,
                                     conditions=[yBCs])

nodeSolve = uw.systems.Solver(nodeDiffuse)
nodeSolve.solve()
with mesh.deform_mesh():
     mesh.data[:,1] = yField.data[:,0]

# If the above is used then it must be repeated for the
# 'mesh0' created further on in this script
########################################################
'''


### SET UP THE BADLANDS MODEL

badlands_model = BadlandsModel()
badlands_model.load_xml('airlayer.xml')
linkage.badlands_model = badlands_model

linkage.material_map = [
    [erodedIndex, airIndex],
    [sedimentIndex, heavyIndex, lightIndex]
]

dem = linkage.generate_flat_dem(minCoord=MIN_COORD, maxCoord=MAX_COORD, resolution=BADLANDS_RESOLUTION, elevation=AIR_ELEVATION)
print 'dem shape %s' % (dem.shape,)
linkage.load_badlands_dem_array(dem)


# Create the swarm, material index variable and swarm advector
swarm = uw.swarm.Swarm(mesh=mesh)
materialIndex = swarm.add_variable(dataType="int", count=1)

swarmLayout = uw.swarm.layouts.GlobalSpaceFillerLayout(swarm=swarm, particlesPerCell=20)
swarm.populate_using_layout(layout=swarmLayout)

advector = uw.systems.SwarmAdvector(swarm=swarm, velocityField=velocityField, order=2)


for index, coord in enumerate(swarm.particleCoordinates.data):
    offset = coord - centre
    if (offset[0]**2 + offset[1]**2 + offset[2]**2 < radius**2):
        materialIndex.data[index] = lightIndex
    elif coord[2] > AIR_ELEVATION:
        materialIndex.data[index] = airIndex
    else:
        materialIndex.data[index] = heavyIndex

# Set viscosities and densities of the model.
viscosityMapFn = 1e19

# Here we set a density of '0.' for the lightMaterial, and '1.' for the heavymaterial.
mappingDictDensity = {#airIndex: 3300.0, 
                      # MEXICO_3
                      airIndex: 3.,
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
                              voronoi_swarm = swarm,
                              conditions    = [velocityBC,],
                              fn_viscosity  = viscosityMapFn,
                              fn_bodyforce  = buoyancyFn)
solver = uw.systems.Solver(stokesPIC)

# some analytics

vdotv = fn.math.dot(velocityField, velocityField)
v2sum_integral = uw.utils.Integral(mesh=mesh, fn=vdotv)
volume_integral = uw.utils.Integral(mesh=mesh, fn=1.)

# configure linkage
linkage.mesh = mesh
linkage.swarm = swarm
linkage.velocity_field = velocityField
linkage.material_index = materialIndex


global first_update
first_update = True


def update_function(linkage, max_seconds):
    global first_update
    if first_update:
        # MEXICO_2
        '''
        with linkage.np_mesh.deform_mesh():
            linkage.np_mesh.data[:, 1] = yField.data[:, 0]
        '''
        first_update = False

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


# run the model
linkage.checkpoint_interval = 10000
linkage.run_for_years(100000)
