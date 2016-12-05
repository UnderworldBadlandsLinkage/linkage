from linkagemodel.smm import SimpleMaterialModel

INITIAL_AIR_ELEVATION = 0.

# material indices
airIndex = 0
heavyIndex = 1
lightIndex = 2
sedimentIndex = 3
erodedIndex = 4

# specify material densities
# NOTE: the first material assigned to a Badlands layer will be the default erosion/deposition type
materials = [
    {'uw_index': erodedIndex, 'density': 3300.0, 'bl_layer': 0},  # bl_layer 0 is air
    {'uw_index': airIndex, 'density': 3300.0, 'bl_layer': 0},
    {'uw_index': sedimentIndex, 'density': 3300.0, 'bl_layer': 1},  # bl_layer 1 is top sediment layer
    {'uw_index': heavyIndex, 'density': 3300.0, 'bl_layer': 1},
    {'uw_index': lightIndex, 'density': 3000.0, 'bl_layer': 1},
]

smm = SimpleMaterialModel(badlands_config='rising_smm.xml',
                          underworld_resolution=[20, 20, 32],
                          materials=materials,
                          # The following are only needed if you do not load a DEM
                          min_coords=(0., 0., -80e3),
                          max_coords=(100e3, 100e3, 20e3),
                          badlands_resolution=(180, 180),
                          surface_elevation=INITIAL_AIR_ELEVATION)

'''
smm = SimpleMaterialModel(badlands_config='rising_smm.xml',
                          underworld_resolution=[20, 20, 32],
                          materials=materials,
                          # If you load a DEM, you need to specify the elevation (Z) range
                          elev_range=(-1e3, 1e3))
'''

SPHERE_RADIUS = 12e3
SPHERE_CENTRE = (50e3, 50e3, -18e3)

# Set up the initial material state
# Define a sphere
for index, coord in enumerate(smm.particle_coordinates.data):
    offset = coord - SPHERE_CENTRE
    if (offset[0]**2 + offset[1]**2 + offset[2]**2 < SPHERE_RADIUS**2):
        smm.material_index.data[index] = lightIndex
    elif coord[2] > INITIAL_AIR_ELEVATION:
        smm.material_index.data[index] = airIndex
    else:
        smm.material_index.data[index] = heavyIndex

smm.checkpoint_interval = 10000
smm.run_for_years(100000)
