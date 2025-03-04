.. _usersguide_depletion:

===========================
Depletion and Transmutation
===========================

OpenMC supports coupled depletion, or burnup, calculations through the
:mod:`openmc.deplete` Python module. OpenMC solves the transport equation to
obtain transmutation reaction rates, and then the reaction rates are used to
solve a set of transmutation equations that determine the evolution of nuclide
densities within a material. The nuclide densities predicted as some future time
are then used to determine updated reaction rates, and the process is repeated
for as many timesteps as are requested.

The depletion module is designed such that the flux/reaction rate solution (the
transport "operator") is completely isolated from the solution of the
transmutation equations and the method used for advancing time. At present, the
:mod:`openmc.deplete` module offers a single transport operator,
:class:`openmc.deplete.Operator` (which uses the OpenMC transport solver), but
in principle additional operator classes based on other transport codes could be
implemented and no changes to the depletion solver itself would be needed. The
operator class requires a :class:`openmc.Geometry` instance and a
:class:`openmc.Settings` instance::

    geom = openmc.Geometry()
    settings = openmc.Settings()
    ...

    op = openmc.deplete.Operator(geom, settings)

Any material that contains a fissionable nuclide is depleted by default, but
this can behavior can be changed with the :attr:`Material.depletable` attribute.

.. important:: The volume must be specified for each material that is depleted by
               setting the :attr:`Material.volume` attribute. This is necessary
               in order to calculate the proper normalization of tally results
               based on the source rate.

:mod:`openmc.deplete` supports multiple time-integration methods for determining
material compositions over time. Each method appears as a different class.
For example, :class:`openmc.deplete.CECMIntegrator` runs a depletion calculation
using the CE/CM algorithm (deplete over a timestep using the middle-of-step
reaction rates). An instance of :class:`openmc.deplete.Operator` is passed to
one of these functions along with the timesteps and power level::

    power = 1200.0e6  # watts
    timesteps = [10.0, 10.0, 10.0]  # days
    openmc.deplete.CECMIntegrator(op, timesteps, power, timestep_units='d').integrate()

The coupled transport-depletion problem is executed, and once it is done a
``depletion_results.h5`` file is written. The results can be analyzed using the
:class:`openmc.deplete.Results` class. This class has methods that allow for
easy retrieval of k-effective, nuclide concentrations, and reaction rates over
time::

    results = openmc.deplete.Results("depletion_results.h5")
    time, keff = results.get_keff()

Note that the coupling between the transport solver and the transmutation solver
happens in-memory rather than by reading/writing files on disk.

Fixed-Source Transmutation
==========================

When the ``power`` or ``power_density`` argument is used for one of the
Integrator classes, it is assumed that OpenMC is running in k-eigenvalue mode,
and normalization of tally results is performed based on energy deposition. It
is also possible to run a fixed-source simulation and perform normalization
based on a known source rate. First, as with all fixed-source calculations, we
need to set the run mode::

    settings.run_mode = 'fixed source'

Additionally, all materials that you wish to deplete need to be marked as such
using the :attr:`Material.depletable` attribute::

    mat = openmc.Material()
    mat.depletable = True

When constructing the :class:`~openmc.deplete.Operator`, you should indicate
that normalization of tally results will be done based on the source rate rather
than a power or power density::

    op = openmc.deplete.Operator(geometry, settings, normalization_mode='source-rate')

Finally, when creating a depletion integrator, use the ``source_rates`` argument::

    integrator = openmc.deplete.PredictorIntegrator(op, timesteps, sources_rates=...)

As with the ``power`` argument, you can provide a different source rate for each
timestep in the calculation. A zero source rate for a given timestep will result
in a decay-only step, where all reaction rates are zero.

Caveats
=======

Energy Deposition
-----------------

The default energy deposition mode, ``"fission-q"``, instructs the
:class:`openmc.deplete.Operator` to normalize reaction rates using the product
of fission reaction rates and fission Q values taken from the depletion chain.
This approach does not consider indirect contributions to energy deposition,
such as neutron heating and energy from secondary photons. In doing this, the
energy deposited during a transport calculation will be lower than expected.
This causes the reaction rates to be over-adjusted to hit the user-specific
power, or power density, leading to an over-depletion of burnable materials.

There are some remedies. First, the fission Q values can be directly set in a
variety of ways. This requires knowing what the total fission energy release
should be, including indirect components. Some examples are provided below::

    # use a dictionary of fission_q values
    fission_q = {"U235": 202e+6}  # energy in eV

    # create a modified chain and write it to a new file
    chain = openmc.deplete.Chain.from_xml("chain.xml", fission_q)
    chain.export_to_xml("chain_mod_q.xml")
    op = openmc.deplete.Operator(geometry, setting, "chain_mod_q.xml")

    # alternatively, pass the modified fission Q directly to the operator
    op = openmc.deplete.Operator(geometry, setting, "chain.xml",
        fission_q=fission_q)


A more complete way to model the energy deposition is to use the modified
heating reactions described in :ref:`methods_heating`.  These values can be used
to normalize reaction rates instead of using the fission reaction rates with::

    op = openmc.deplete.Operator(geometry, settings, "chain.xml",
        normalization_mode="energy-deposition")

These modified heating libraries can be generated by running the latest version
of :meth:`openmc.data.IncidentNeutron.from_njoy`, and will eventually be bundled
into the distributed libraries.

Local Spectra and Repeated Materials
------------------------------------

It is not uncommon to explicitly create a single burnable material across many
locations. From a pure transport perspective, there is nothing wrong with
creating a single 3.5 wt.% enriched fuel ``fuel_3``, and placing that fuel in
every fuel pin in an assembly or even full core problem. This certainly
expedites the model making process, but can pose issues with depletion. Under
this setup, :mod:`openmc.deplete` will deplete a single ``fuel_3`` material
using a single set of reaction rates, and produce a single new composition for
the next time step. This can be problematic if the same ``fuel_3`` is used in
very different regions of the problem.

As an example, consider a full-scale power reactor core with vacuum boundary
conditions, and with fuel pins solely composed of the same ``fuel_3`` material.
The fuel pins towards the center of the problem will surely experience a more
intense neutron flux and greater reaction rates than those towards the edge of
the domain. This indicates that the fuel in the center should be at a more
depleted state than periphery pins, at least for the fist depletion step.
However, without any other instructions, OpenMC will deplete ``fuel_3`` as a
single material, and all of the fuel pins will have an identical composition at
the next transport step.

This can be countered by instructing the operator to treat repeated instances
of the same material as a unique material definition with::

    op = openmc.deplete.Operator(geometry, settings, chain_file,
        diff_burnable_mats=True)

For our example problem, this would deplete fuel on the outer region of the
problem with different reaction rates than those in the center. Materials will
be depleted corresponding to their local neutron spectra, and have unique
compositions at each transport step.  The volume of the original ``fuel_3``
material must represent the volume of **all** the ``fuel_3`` in the problem.
When creating the unique materials, this volume will be equally distributed
across all material instances.


.. note::

    This will increase the total memory usage and run time due to an increased
    number of tallies and material definitions.

