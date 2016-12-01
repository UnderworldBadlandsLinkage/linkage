# DEM demo
# Load a DEM into Badlands and verify that it is transferred to Underworld

# After you run this, you can open material-0.xdmf and tin.series.xdmf in
# ParaView. Setting the Z scale to 10 will make things clearer. You should also
# change the colour values on the material

import numpy
import os
import pandas

from pyBadlands.model import Model as BadlandsModel
import underworld as uw
from underworld import function as fn
from linkagemodel.linkage import LinkageModel

# We're going to load from a DEM; let's just load it to determine the Underworld mesh size
dem = pandas.read_csv('nodes.csv', sep=' ', header=None, na_filter=False, dtype=numpy.float, low_memory=False)

# set up min/max coord for UW
min_coord = [dem[0].min(), dem[1].min(), -20e3]
max_coord = [dem[0].max(), dem[1].max(), +20e3]

print 'The model has bounds %s to %s' % (min_coord, max_coord)

linkage = LinkageModel()


### SET UP THE BADLANDS MODEL
badlands_model = BadlandsModel()
badlands_model.load_xml('dem.xml')  # The XML file specifies the DEM to load
linkage.badlands_model = badlands_model


### SET UP THE UNDERWORLD MODEL

# All output will go to the 'uwout' directory, which we will create
uw_output_path = 'uwout'
try:
    os.mkdir(uw_output_path)
except OSError:
    # probably already exists
    pass

UNDERWORLD_RESOLUTION = 20

# This is the mesh whose material types will be changed by Badlands
mesh = uw.mesh.FeMesh_Cartesian(elementType=("Q1/dQ0"),
                                elementRes =[UNDERWORLD_RESOLUTION] * 3,
                                minCoord   =min_coord,
                                maxCoord   =max_coord)

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

# Set the whole thing to be 'sediment'. At initialisation, Badlands will dividie it into sediment/air.
materialIndex.data[:] = linkage.material_map[1][0]

advector = uw.systems.SwarmAdvector(swarm=swarm, velocityField=velocityField, order=2)

# Set viscosities and densities of the model.
viscosityMapFn = 1e19
mappingDictDensity = {#0: 3300.0,  # FIXME: there shouldn't be a 0
                      linkage.deposited_material_index: 3300.0,
                      linkage.eroded_material_index: 3300.0}
densityFn = fn.branching.map(fn_key=materialIndex, mapping=mappingDictDensity)
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
# Just initialise. Don't run.
# linkage.run_for_years(0)

# Run a tiny bit so Badlands produces a checkpoint.
linkage.run_for_years(1)
# material-0.xdmf should show the landscape, and it should look identical to the Badlands tin.xdmf.
