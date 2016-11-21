class LinkageModel(object):
    """
    A LinkageModel joins an Underworld and a Badlands model. Underworld models
    the domain in 3D and Badlands models the surface processes.

    LinkageModel hides the details of how the two models communicate so you can
    concentrate on building them.

    To define a linked model, create a subclass of LinkageModel and override
    the badlands_init() and underworld_init() functions. You can then run the
    model by calling the run_for_years() function.

    Output will be saved as files on disk.

    The raw Badlands model object can be accessed as the 'badlands' member. The
    raw Underworld model object can be accessed as the 'underworld' member.
    """

    def __init__(self):
        self.t = 0.  # Simulation time in years. We start at year 0.
        self.badlands = 
        self.underworld = 

    def badlands_init(self):
        """
        Define the initial state of the Badlands model. Many models will load a
        DEM (using self.badlands.load_dem(filename)). You can create a flat
        DEM with the self.generate_flat_dem() function and then load it using
        load_dem().

        You can also load an XML configuration file with something like
        self.badlands.load_xml(filename).

        You can also override settings made in the XML input file at this time.
        """
        raise NotImplementedError("The badlands_init() function must be overridden.")

    def underworld_init(self):
        """
        Define the initial state of the Underworld model.
        """
        raise NotImplementedError("The underworld_init() function must be overridden.")

    def on_tick(self, years):
        """
        If you want to modify the model as it runs, override this function. It
        will be called at each iteration. The simulation time (in years) will
        be passed into the 'years' parameter.

        If you don't need this, you don't need to override it.
        """
        pass

    def run_for_years(self, years):
        """
        Run the model for a number of years.
        """
        todo

    def generate_flat_dem(self< paramstodo):
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
        review this function
