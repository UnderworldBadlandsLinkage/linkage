
# coding: utf-8

# 3D raising blob
# run in parallel with: 'mpirun -np 4 python parallelTest.py'
# --------
# -----
# 
# * Simple 3D isoviscous mechanical model
# * A negatively buoyant blob ascends from 'mantle' into 'air'. 
# * This model is for testing the initial coupling of Underworld and Badlands (non representative units used)
# 
# ** DEV notes **
# * See the main loop for coords_bl, v_bl and vz numpy arrays:
#   coords_bl - the coordinates of the tracer location in underworld
#   v_bl      - the entire velocity vector at those locations
#   vz        - only the vertical component at those location (commented out)
# * 'plotly' is needed to visualise the final cell, download it with pip


## TODO: you probably want to scale the badlands mesh size; units are meters there, so we want to see 10kunits on each axis

# In[ ]:

import bridge

import numpy as np
from pyBadlands.model import Model as badlandsModel
import underworld as uw
import math
from underworld import function as fn
import glucifer
from mpl_toolkits.mplot3d import Axes3D
from mpi4py import MPI
from scipy.interpolate import griddata

comm = MPI.COMM_WORLD
rank = uw.rank()


# In[ ]:

# Create output directory, parallel safe

output_dir="./raisingSphere_parallel/" 

if uw.rank() == 0: # only rank 0 creates directory
    try:
        import os
        if os.path.exists('./'+output_dir+'/'):
            # ALWAYS DELETE directory if found!
            import shutil
            shutil.rmtree('./'+output_dir+'/')
        
        os.makedirs("./"+output_dir+"/")
    except OSError:
        raise

comm.Barrier()  # wait for all


# Create mesh and associated mesh variables
MIN_COORD = (0., 0., -80e3)
MAX_COORD = (100e3, 100e3, 20e3)
# MEXICO_1
UNDERWORLD_RESOLUTION = [20,20,32]
# BADLANDS_RESOLUTION = 60
BADLANDS_RESOLUTION = 180
AIR_ELEVATION = 0.0  # the height at which we transition from sediment to air

# Set size and position of cube
radius = 10e3
centre = (50e3, 50e3, -20e3)  # make sure the cube starts underground. We can't model non-flat starting surfaces in Badlands yet.
# centre = (50e3, 50e3, -30e3)  # make sure the cube starts underground. We can't model non-flat starting surfaces in Badlands yet.

# Uncomment this to verify that coordinate systems match up
# cubeSize = (1., 4., 4.)
# centre = (10., 14., 8.)  # make sure the cube starts underground. We can't model non-flat starting surfaces in Badlands yet.


mesh = uw.mesh.FeMesh_Cartesian( elementType = ("Q1/dQ0"), 
                                 elementRes  = UNDERWORLD_RESOLUTION, 
                                 minCoord    = MIN_COORD, 
                                 maxCoord    = MAX_COORD)

velocityField    = uw.mesh.MeshVariable( mesh=mesh,         nodeDofCount=mesh.dim )
pressureField    = uw.mesh.MeshVariable( mesh=mesh.subMesh, nodeDofCount=1 )

velocityField.data[:] = [0.,0.,0.]
pressureField.data[:] = 0.


# MEXICO_2
###  HOW TO DEFORM THE MESH FOR FOCUSSED RESOLUTION ####

## we use the poisson equation, varying diffusivity to solve for
## new distribution of points, in the vertical axis

# yField    = uw.mesh.MeshVariable( mesh=mesh, nodeDofCount=1 )
# coord = fn.coord()
# yBCs = uw.conditions.DirichletCondition( variable    = yField, 
#                                          indexSetsPerDof = ( mesh.specialSets["MinJ_VertexSet"]+mesh.specialSets["MaxJ_VertexSet"] ) )

# yField.data[mesh.specialSets["MinJ_VertexSet"].data] = MIN_COORD[2]
# yField.data[mesh.specialSets["MaxJ_VertexSet"].data] = MAX_COORD[2]

# intensityFn = 1.0+6.0*uw.function.math.exp(-10.0*uw.function.math.pow((coord[1]-0.4),2.0))

# nodeDiffuse = uw.systems.SteadyStateHeat(yField, 
#                                      fn_diffusivity=intensityFn,
#                                      fn_heating=0.0,
#                                      conditions=[yBCs])

# nodeSolve = uw.systems.Solver(nodeDiffuse)
# nodeSolve.solve()
# with mesh.deform_mesh():
#      mesh.data[:,1] = yField.data[:,0]

## If the above is used then it must be repeated for the
## 'mesh0' created further on in this script
#########################################################


# In[ ]:

# init badlands

# TODO: we need to set up a mesh in badlands with configurable resolution and using the same coordinate system as above

# NOTE: Badlands currently assumes that the mesh is square and has equal resolution on each axis

blModel = badlandsModel()

# Build the initial mesh for Badlands
# FIXME: I haven't translated the coord system from UW->BL yet

items = []
# FIXME: there should be a fast numpy way to do this
for y in np.linspace(MIN_COORD[1], MAX_COORD[1], BADLANDS_RESOLUTION):
    for x in np.linspace(MIN_COORD[0], MAX_COORD[0], BADLANDS_RESOLUTION):
        items.append([x, y, AIR_ELEVATION])

# NOTE: Badlands uses the difference in X coord of the first two points to determine the resolution.
# This is something we should fix.
# This is why we loop in y/x order instead of x/y order.
blMesh = np.array(items)
np.savetxt('flat/nodes.csv', blMesh)
blModel.load_xml('flat/flat.xml')


# Create the swarm, material index variable and swarm advector

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
			# MEXICO_3
			# airIndex: 3
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

### FOR OUTPUT ###

# create checkpoint function
def checkpoint(mesh, fieldDict, swarm, swarmDict, index, modeltime=0., prefix=None, enable_xdmf=True):
    import os
    # Check the prefix is valid
    if prefix is not None:
        if not os.path.exists(prefix):
            raise ValueError("prefix given '{}' doesn't exist".format(prefix))

    if not isinstance(index, int):
        raise TypeError("'index' is not of type int")
    ii = str(index)

    if mesh is not None:
        # Error check the mesh and fields
        if not isinstance(mesh, uw.mesh.FeMesh):
            raise TypeError("'mesh' is not of type uw.mesh.FeMesh")
        if not isinstance(fieldDict, dict):
            raise TypeError("'fieldDict' is not of type dict")
        for key, value in fieldDict.iteritems():
            if not isinstance(value, uw.mesh.MeshVariable):
                raise TypeError("'fieldDict' must contain uw.mesh.MeshVariable elements")

        # see if we have already saved the mesh. It only needs to be saved once
        if not hasattr(checkpoint, 'mH'):
            checkpoint.mH = mesh.save(prefix + "mesh.h5")
        mh = checkpoint.mH

        for key, value in fieldDict.iteritems():
            filename = prefix + key + '-' + ii
            handle = value.save(filename + '.h5')
            if enable_xdmf:
                value.xdmf(filename, handle, key, mh, 'mesh', modeltime=modeltime)

    # is there a swarm
    if swarm is not None:
        # Error check the swarms
        if not isinstance(swarm, uw.swarm.Swarm):
            raise TypeError("'swarm' is not of type uw.swarm.Swarm")
        if not isinstance(swarmDict, dict):
            raise TypeError("'swarmDict' is not of type dict")
        for key, value in swarmDict.iteritems():
            if not isinstance( value, uw.swarm.SwarmVariable ):
                raise TypeError("'fieldDict' must contain uw.swarm.SwarmVariable elements")
    
        sH = swarm.save(prefix+"swarm-"+ii+".h5")
        for key,value in swarmDict.iteritems():
            filename = prefix+key+'-'+ii
            handle = value.save(filename+'.h5')
            if enable_xdmf: value.xdmf(filename, handle, key, sH, 'swarm', modeltime=modeltime)
                

# setup of checkpoint dictionaries and initial checkpoint

fields = {'velocity':velocityField, 'pressure':pressureField}
swarmV = {'material':materialIndex}

checkpoint(None, None, swarm, swarmV, 0, prefix=output_dir)


def determine_particle_state(volume, blModel):
    # given badlands' mesh, determine if each particle in 'volume' is above (False) or below (True) it.
    
    # To do this, for each X/Y pair in 'volume', we interpolate its Z value relative to the mesh in blModel. Then,
    # if the interpolated Z is greater than the supplied Z (i.e. Badlands mesh is above particle elevation) it's 
    # sediment (True). Else, it's air (False).
    
    # This is probably a really slow way to go about this problem, but it should be quick to implement.

    # Technically, tinMesh could represent a mesh that overlaps itself in the Z
    # axis. I'm going to assume that it doesn't. This will very slightly reduce
    # the accuracy of the interpolation but it shouldn't matter.


    known_xy = blModel.recGrid.tinMesh['vertices']  # points that we have known elevation for
    known_z = blModel.elevation  # elevation for those points

    interpolate_xy = volume[:, [0, 1]]
    # linear interpolation should be plenty as we're running badlands at higher resolution than underworld
    interpolate_z = griddata(points=known_xy, values=known_z, xi=interpolate_xy, method='linear')

    # True for sediment, False for air
    flags = volume[:, 2] < interpolate_z

    return flags


# Stepping. Initialise time and timestep.
SECONDS_PER_YEAR = float(365 * 24 * 60 * 60)

# checkpoint number
step = 0

# Finish the simulation after this many years
END_TIME_YEARS = 100000.
# TODO: patch the badlands simulation end time to match this so you don't run over the end of the xml-specified runtime

# Current time
time_years = 0.

# We will interchange data between UW and BL on every iteration
# We only generate output (checkpoints) at specified intervals
CHECKPOINT_INTERVAL_YEARS = 10000.

# At which year will we produce another checkpoint?
next_checkpoint_years = CHECKPOINT_INTERVAL_YEARS

# build a non-partitioned mesh with same box size
mesh0 = uw.mesh.FeMesh_Cartesian( elementType = ("Q1/dQ0"), 
			      elementRes  = UNDERWORLD_RESOLUTION,
			      minCoord    = MIN_COORD, 
			      maxCoord    = MAX_COORD,
			      partitioned = False)

# MEXICO_2
# with mesh0.deform_mesh():
#      mesh0.data[:,1] = yField.data[:,0]


# Main loop
while time_years < END_TIME_YEARS:
    vrms = math.sqrt(v2sum_integral.evaluate()[0] / volume_integral.evaluate()[0])
    # Get solution for initial configuration.
    solver.solve()

    # Retrieve the maximum possible timestep for the advection system.
    dtmax_seconds = advector.get_max_dt()
    assert dtmax_seconds != np.inf

    dtmax_years = dtmax_seconds / SECONDS_PER_YEAR

    if next_checkpoint_years == time_years:
        # checkpoint fields and swarm
        # UW writes a checkpoint at the START of the timestep
        checkpoint(velocityField.mesh, fields, swarm, swarmV, step, prefix=output_dir, modeltime=time_years)
        step += 1

        next_checkpoint_years += CHECKPOINT_INTERVAL_YEARS

    if dtmax_years + time_years > next_checkpoint_years:
        # we're going to pass the checkpoint mark
        dt_years = next_checkpoint_years - time_years
        endtime_years = next_checkpoint_years  # avoid floating point issues by doing this here
        # on the next iteration, we will write the checkpoint and advance next_checkpoint_years
    else:
        dt_years = dtmax_years
        endtime_years += dt_years

    dt_seconds = dt_years * SECONDS_PER_YEAR


    # Advect using this timestep size.
    advector.integrate(dt_seconds)
    advector2.integrate(dt_seconds)

    # Each CPU saves its view of the velocity field so it can be reconstructed everywhere
    mesh.save('temp/mesh-%s.h5' % step)
    velocityField.save('temp/velfield-%s.h5' % step)
    bl_tracers.save('temp/tracers-%s.h5' % step)

    # load previous mesh coordinate data onto new non-partitioned mesh
    mesh0.load('temp/mesh-%s.h5' % step)

    velocityField0 = uw.mesh.MeshVariable(mesh=mesh0, nodeDofCount=mesh0.dim)
    velocityField0.load('temp/velfield-%s.h5' % step)

    tracers0 = uw.swarm.Swarm(mesh0)
    tracers0.load('temp/tracers-%s.h5' % step)

    # tracers0 contains the tracers across all nodes
    # FIXME badlands might want to define the tracer locations to coincide with the surface mesh
    tracer_velocity_mps = velocityField0.evaluate(tracers0)  # the entire velocity vector on each particle in METERS PER SECOND

    ### INTERFACE PART 1: UW->BL
    # Use the tracer vertical velocities to deform the Badlands TIN
    # from meters per second to meters displacement over the whole iteration
    tracer_disp = tracer_velocity_mps * SECONDS_PER_YEAR * dt_years

    bridge.inject_mesh_to_badlands_3d(tracers0.particleCoordinates.data, tracer_disp, blModel, time_years, dt_years, n=BADLANDS_RESOLUTION, min_coord=MIN_COORD, max_coord=MAX_COORD, display_interval=CHECKPOINT_INTERVAL_YEARS)
    # FIXME: you may want to put in the tin nodes as tracers so you can get new velocity for them right out of underworld

    blModel.run_to_time(endtime_years)
    #####################

    ### INTERFACE PART 2: BL->UW

    # UW switch materials
    flags = determine_particle_state(swarm.particleCoordinates.data, blModel)
    for index, material in enumerate(materialIndex.data):
        # convert air to sediment
        if material in [airIndex, erodedIndex] and flags[index]:
            materialIndex.data[index] = sedimentIndex
        # convert everything-but-air to air
        if material not in [airIndex, erodedIndex] and not flags[index]:
            # materialIndex.data[index] = airIndex
            materialIndex.data[index] = erodedIndex

    ### ADVANCE TIME
    time_years = endtime_years

# write a final checkpoint from uw
checkpoint(velocityField.mesh, fields, swarm, swarmV, step, prefix=output_dir, modeltime=time_years)

# TODO: write an initial state from BL (before any displacement from UW)
