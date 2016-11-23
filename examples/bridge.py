import numpy as np
from scipy.interpolate import griddata

'''
# uncomment if you want plotting
import matplotlib
matplotlib.use('Agg')  # noninteractive
import matplotlib.pyplot as plt
'''

# Underworld uses coord system:
#     x is horizontal (right positive)
#     y is vertical (up positive)
#     z is into the screen (positive is closer to the user)
# Badlands uses coord system (?? TODO CHECK THIS)
#     x is horizontal (right positive)
#     y is into the screen (positive away from user)
#     z is vertical (up positive)
# so we need to exchange y and z
# the direction of underworld's z might be backwards

# TODO: the coord systems appear to be the same now. Can we remove this?
def uw_to_bl_coords(arr):
    '''Convert a three-column numpy array from Underworld to Badlands coordinate systems'''
    return np.column_stack([arr[:, 0], arr[:, 2], arr[:, 1]])


def write_initial_badlands_mesh(mesh, output_filename, transition_elev):
    # LATER: would be nice if we could do this all in-memory
    # FIXME: this function is a kludgey mess
    
    epsilon = 0.01

    bl_mesh = uw_to_bl_coords(mesh.data)

    # FIXME YUCKY - we're trying to take a 2D slice of the mesh
    # this is probably the wrong approach; rethink this
    # this also breaks if the resolution isn't quite right
    mask = np.logical_and(bl_mesh[:, 2] > (transition_elev - epsilon), bl_mesh[:, 2] < (transition_elev + epsilon))

    masked_bl_mesh = bl_mesh[mask]

    # this should be shape (res + 1, res + 1) if the mask operation was done correctly
    print masked_bl_mesh.shape

    np.savetxt(output_filename, masked_bl_mesh)


global disp_inserted
disp_inserted = False


def inject_mesh_to_badlands_3d(tracers, disp, blModel, time, dt, n, min_coord, max_coord, display_interval=None):
    """
    Takes a plane of tracer points and their DISPLACEMENTS in 3D over time
    period dt. Injects it into Badlands as 3D tectonic movement.
    """

    # enable 3D displacements
    blModel.input.disp3d = True

    # TODO: check what this does
    blModel.input.region = 0

    if display_interval is not None:
        blModel.input.tDisplay = display_interval

    blModel.force.merge3d = blModel.input.Afactor * blModel.recGrid.resEdges * 0.5

    # NOTE: we (Badlands) are going to cheat a bit here. Instead of
    # interpolating on top of the TIN, we're going to interpolate on top
    # of a DEM. This will cost a little accuracy. I suspect it will not
    # be noticeable.
    #
    # We can interpolate on the TIN later; it just takes more time to 
    # implement.

    # Interpolate back to grid
    grid_y, grid_x = np.mgrid[min_coord[0]:max_coord[0]:n * 1j, min_coord[1]:max_coord[1]:n * 1j] 

    # I'm not entirely certain, but I think the interpretation of the 3d
    # displacement map is the displacement of each DEM node at the end of the
    # time period relative to its starting position. If you start a new
    # displacement file, it is treated as starting at the DEM starting points
    # (and interpolated onto the TIN as it as at that tNow).

    # see forceSim.py:472 if you want to interpolate directly onto the DEM
    tracer_xy = np.column_stack([tracers[:, 0], tracers[:, 1]])

    grid_dx = griddata(tracer_xy, disp[:, 0], (grid_x, grid_y), method='cubic')
    grid_dy = griddata(tracer_xy, disp[:, 1], (grid_x, grid_y), method='cubic')
    grid_dz = griddata(tracer_xy, disp[:, 2], (grid_x, grid_y), method='cubic')

    # FIXME kludge; don't keep added new entries
    global disp_inserted
    if disp_inserted:
        print 'INITIAL INSERT'

        blModel.force.T_disp[0, 0] = time
        blModel.force.T_disp[0, 1] = (time + dt)
    else:
        print 'SUBSEQUENT INSERT'
        blModel.force.T_disp = np.vstack(([time, time + dt], blModel.force.T_disp))
        disp_inserted = True

    # grid_z is the new displacement
    # make a 3D displacement array
    disp = np.column_stack((grid_dx.flatten(), grid_dy.flatten(), grid_dz.flatten()))

    blModel.force.injected_disps = disp

    # It looks like badlands doesn't apply interpolation on every tick; rather,
    # just on the interval specified here. This should be smaller than the
    # timesteps it's taking but large so your run doesn't take too long.
    # TODO: talk to Tristan about this, seems weird; can't we just update on
    # every timestep?

    # COMMENT THIS OUT TO DISABLE INTERPOLATION
    # which is not ideal but that's what we're testing
    # blModel.force.time3d = 100   # FIXME: does time3d refer to dt or bl_endtime?
