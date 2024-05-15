#!/usr/bin/env python
"""
MASC creates a model star cluster, which can then be used in N-body simulations
or for other purposes.

It requires AMUSE, which can be downloaded from http://amusecode.org or
https://github.com/amusecode/amuse.

Currently not feature-complete yet, and function/argument names are
subject to change.
"""

import logging
import numpy

from amuse.units import (
    units,
    nbody_system,
    generic_unit_converter,
)
from amuse.units.trigo import sin, cos
from amuse.datamodel.particles import Particles
from amuse.ic.plummer import new_plummer_sphere
from amuse.ic.kingmodel import new_king_model
from .stars import new_masses
from .binaries import new_binary_distribution
try:
    from amuse.ic.fractalcluster import new_fractal_cluster_model
except ImportError:
    new_fractal_cluster_model = None


def new_star_cluster(
    stellar_mass=None,
    initial_mass_function="salpeter",
    upper_mass_limit=125. | units.MSun,
    lower_mass_limit=0.1 | units.MSun,
    number_of_stars=1024,
    number_of_binaries=0,
    effective_radius=3.0 | units.parsec,
    star_distribution="plummer",
    star_distribution_w0=7.0,
    star_distribution_fd=2.0,
    star_metallicity=0.01,
    initial_binary_mass_fraction=0,
    binary_semi_major_axis_distribution=None,
    eccentricity_distribution=None,
    return_binaries=False,
    **kwargs
):
    """
    Takes properties for a star cluster and returns the stars.
    When using an IMF, either the stellar mass is fixed (within
    stochastic error) or the number of stars is fixed. When using
    equal-mass stars, both are fixed.
    stellar_mass takes precedence over number_of_stars!
    """

    # sanity checks
    if stellar_mass is None and number_of_stars is None:
        print("Error: must provide either stellar_mass or number_of_stars")
        return -1
    if stellar_mass is None:
        if number_of_binaries < 0:
            print("number_of_binaries must be >= 0")
            return -1
        if number_of_binaries > number_of_stars:
            print("number_of_binaries must be <= number_of_stars")
            return -1

    # First, generate the masses of the systems (single stars and binary stars)
    def generate_single_stars():
        if stellar_mass is not None:
            initial_single_star_mass_fraction = (
                1
                - initial_binary_mass_fraction
            )
            single_stars_mass = (
                initial_single_star_mass_fraction * stellar_mass
            )
        elif number_of_stars is not None:
            number_of_single_stars = number_of_stars - number_of_binaries
            single_stars_mass = None

        single_star_masses = new_masses(
            stellar_mass=single_stars_mass,
            initial_mass_function=initial_mass_function,
            upper_mass_limit=upper_mass_limit,
            lower_mass_limit=lower_mass_limit,
            number_of_stars=number_of_single_stars,
        )
        stars = Particles(mass=single_star_masses)
        stars.position = [0, 0, 0] | units.pc
        stars.velocity = [0, 0, 0] | units.kms
        return stars

    def generate_binary_stars():
        if stellar_mass is not None:
            total_mass = initial_binary_mass_fraction * stellar_mass
            primary_masses = None
        else:
            total_mass = None
            primary_masses = new_masses(
                number_of_stars=number_of_binaries,
                initial_mass_function=initial_mass_function,
                upper_mass_limit=upper_mass_limit,
                lower_mass_limit=lower_mass_limit,
            )

        stars, binaries = new_binary_distribution(
            primary_masses=primary_masses,
            total_mass=total_mass,
            initial_mass_function=initial_mass_function,
            upper_mass_limit=upper_mass_limit,
            lower_mass_limit=lower_mass_limit,
            semi_major_axis_distribution=binary_semi_major_axis_distribution,
            eccentricity_distribution=eccentricity_distribution,
        )
        return stars, binaries

    single_stars = generate_single_stars()
    binary_component_stars, binary_pairs = generate_binary_stars()
    all_stars = single_stars | binary_component_stars
    all_systems = single_stars | binary_pairs

    number_of_single_stars = len(single_stars)
    number_of_stars = len(all_stars)
    number_of_systems = len(all_systems)
    number_of_binaries = len(binary_pairs)

    if number_of_binaries > 0:
        total_mass = single_stars.total_mass() + binary_pairs.total_mass()
    else:
        total_mass = single_stars.total_mass()

    converter = generic_unit_converter.ConvertBetweenGenericAndSiUnits(
        total_mass,
        1. | units.kms,
        effective_radius,
    )
    # Give stars a position and velocity, based on the distribution model.
    if star_distribution == "plummer":
        spatial_distribution = new_plummer_sphere(
            number_of_systems,
            convert_nbody=converter,
        )
    elif star_distribution == "king":
        spatial_distribution = new_king_model(
            number_of_systems,
            star_distribution_w0,
            convert_nbody=converter,
        )
    elif star_distribution == "fractal":
        spatial_distribution = new_fractal_cluster_model(
            number_of_systems,
            fractal_dimension=star_distribution_fd,
            convert_nbody=converter,
        )
    else:
        return -1, "No stellar distribution"

    # randomise spatial_distribution (needed for fractal distribution in
    # particular)
    random_indices = numpy.arange(len(spatial_distribution))
    numpy.random.shuffle(random_indices)

    # move systems (not stars!) to their positions
    all_systems.position += spatial_distribution[random_indices].position
    all_systems.velocity += spatial_distribution[random_indices].velocity

    # set other stellar parameters.
    all_stars.metallicity = star_metallicity

    # Virialize the star cluster if > 1 star
    if number_of_systems > 1:
        all_systems.move_to_center()
        all_systems.scale_to_standard(
            convert_nbody=converter,
            # virial_ratio=virial_ratio,
            # smoothing_length_squared= ...,
        )
    if number_of_binaries > 0:
        # make sure the component stars of the binary systems are also placed
        # correctly
        binary_component_stars[0:number_of_binaries].position += \
            binary_pairs.position
        binary_component_stars[number_of_binaries:].position += \
            binary_pairs.position
        binary_component_stars[0:number_of_binaries].velocity += \
            binary_pairs.velocity
        binary_component_stars[number_of_binaries:].velocity += \
            binary_pairs.velocity

    # Record the cluster's initial parameters to the particle distribution
    all_stars.collection_attributes.initial_mass_function = \
        initial_mass_function.lower()
    all_stars.collection_attributes.upper_mass_limit = upper_mass_limit
    all_stars.collection_attributes.lower_mass_limit = lower_mass_limit
    all_stars.collection_attributes.number_of_stars = number_of_stars
    all_stars.collection_attributes.effective_radius = effective_radius
    all_stars.collection_attributes.star_distribution = star_distribution
    all_stars.collection_attributes.star_distribution_w0 = star_distribution_w0
    all_stars.collection_attributes.star_distribution_fd = star_distribution_fd
    all_stars.collection_attributes.star_metallicity = star_metallicity

    # Derived/legacy values
    all_stars.collection_attributes.converter_mass = \
        converter.to_si(1 | nbody_system.mass)
    all_stars.collection_attributes.converter_length =\
        converter.to_si(1 | nbody_system.length)
    all_stars.collection_attributes.converter_speed =\
        converter.to_si(1 | nbody_system.speed)

    if return_binaries:
        return single_stars, binary_component_stars, binary_pairs
    return all_stars


def new_stars_from_sink(
        origin,
        upper_mass_limit=125 | units.MSun,
        lower_mass_limit=0.1 | units.MSun,
        default_radius=0.25 | units.pc,
        velocity_dispersion=1 | units.kms,
        logger=None,
        initial_mass_function="kroupa",
        distribution="random",
        randomseed=None,
        **keyword_arguments
):
    """
    Form stars from an origin particle that keeps track of the properties of
    this region.
    """
    logger = logger or logging.getLogger(__name__)
    if randomseed is not None:
        logger.info("setting random seed to %i", randomseed)
        numpy.random.seed(randomseed)

    try:
        initialised = origin.initialised
    except AttributeError:
        initialised = False
    if not initialised:
        logger.debug(
            "Initialising origin particle %i for star formation",
            origin.key
        )
        next_mass = new_star_cluster(
            initial_mass_function=initial_mass_function,
            upper_mass_limit=upper_mass_limit,
            lower_mass_limit=lower_mass_limit,
            number_of_stars=1,
            **keyword_arguments
        )
        origin.next_primary_mass = next_mass[0].mass
        origin.initialised = True

    if origin.mass < origin.next_primary_mass:
        logger.debug(
            "Not enough in star forming region %i to form the next star",
            origin.key
        )
        return Particles()

    mass_reservoir = origin.mass - origin.next_primary_mass
    stellar_masses = new_star_cluster(
        stellar_mass=mass_reservoir,
        upper_mass_limit=upper_mass_limit,
        lower_mass_limit=lower_mass_limit,
        imf=initial_mass_function,
    ).mass
    number_of_stars = len(stellar_masses)

    new_stars = Particles(number_of_stars)
    new_stars.age = 0 | units.yr
    new_stars[0].mass = origin.next_primary_mass
    new_stars[1:].mass = stellar_masses[:-1]
    origin.next_primary_mass = stellar_masses[-1]
    new_stars.position = origin.position
    new_stars.velocity = origin.velocity

    try:
        radius = origin.radius
    except AttributeError:
        radius = default_radius
    rho = numpy.random.random(number_of_stars) * radius
    theta = (
        numpy.random.random(number_of_stars)
        * (2 * numpy.pi | units.rad)
    )
    phi = (
        numpy.random.random(number_of_stars) * numpy.pi | units.rad
    )
    x = rho * sin(phi) * cos(theta)
    y = rho * sin(phi) * sin(theta)
    z = rho * cos(phi)
    new_stars.x += x
    new_stars.y += y
    new_stars.z += z

    velocity_magnitude = numpy.random.normal(
        scale=velocity_dispersion.value_in(units.kms),
        size=number_of_stars,
    ) | units.kms
    velocity_theta = (
        numpy.random.random(number_of_stars)
        * (2 * numpy.pi | units.rad)
    )
    velocity_phi = (
        numpy.random.random(number_of_stars)
        * (numpy.pi | units.rad)
    )
    vx = velocity_magnitude * sin(velocity_phi) * cos(velocity_theta)
    vy = velocity_magnitude * sin(velocity_phi) * sin(velocity_theta)
    vz = velocity_magnitude * cos(velocity_phi)
    new_stars.vx += vx
    new_stars.vy += vy
    new_stars.vz += vz

    new_stars.origin = origin.key
    origin.mass -= new_stars.total_mass()

    return new_stars
