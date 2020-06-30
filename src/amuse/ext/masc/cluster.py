#!/usr/bin/env python
"""
make_a_star_cluster.py creates a model star cluster,
which can then be used in N-body simulations or for other purposes.

It requires AMUSE, which can be downloaded from http://amusecode.org or
https://github.com/amusecode/amuse.

Currently not feature-complete yet, and function/argument names are
subject to change.

-- Steven Rieder steven at rieder punt nl
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
try:
    from amuse.ic.fractalcluster import new_fractal_cluster_model
except ImportError:
    new_fractal_cluster_model = None


def new_star_cluster(
        stellar_mass=False,
        initial_mass_function="salpeter",
        upper_mass_limit=125. | units.MSun,
        lower_mass_limit=0.1 | units.MSun,
        number_of_stars=1024,
        effective_radius=3.0 | units.parsec,
        star_distribution="plummer",
        star_distribution_w0=7.0,
        star_distribution_fd=2.0,
        star_metallicity=0.01,
        # initial_binary_fraction=0,
        **kwargs
):
    """
    Create stars.
    When using an IMF, either the stellar mass is fixed (within
    stochastic error) or the number of stars is fixed. When using
    equal-mass stars, both are fixed.
    """

    imf_name = initial_mass_function.lower()
    if imf_name == "salpeter":
        from amuse.ic.salpeter import new_salpeter_mass_distribution
        initial_mass_function = new_salpeter_mass_distribution
    elif imf_name == "kroupa":
        from amuse.ic.brokenimf import new_kroupa_mass_distribution
        initial_mass_function = new_kroupa_mass_distribution
    elif imf_name == "flat":
        from amuse.ic.flatimf import new_flat_mass_distribution
        initial_mass_function = new_flat_mass_distribution
    elif imf_name == "fixed":
        from amuse.ic.flatimf import new_flat_mass_distribution

        def new_fixed_mass_distribution(
                number_of_particles, *list_arguments, **keyword_arguments
        ):
            return new_flat_mass_distribution(
                number_of_particles,
                mass_min=stellar_mass/number_of_stars,
                mass_max=stellar_mass/number_of_stars,
            )
        initial_mass_function = new_fixed_mass_distribution

    if stellar_mass:
        # Add stars to cluster, until mass limit reached (inclusive!)
        mass = initial_mass_function(
            0,
            mass_min=lower_mass_limit,
            mass_max=upper_mass_limit,
        )
        while mass.sum() < stellar_mass:
            mass.append(
                initial_mass_function(
                    1,
                    mass_min=lower_mass_limit,
                    mass_max=upper_mass_limit,
                )[0]
            )
        total_mass = mass.sum()
        number_of_stars = len(mass)
    else:
        # Give stars their mass
        mass = initial_mass_function(
            number_of_stars,
            mass_min=lower_mass_limit,
            mass_max=upper_mass_limit,
        )
        total_mass = mass.sum()

    converter = generic_unit_converter.ConvertBetweenGenericAndSiUnits(
        total_mass,
        1. | units.kms,
        effective_radius,
    )
    # Give stars a position and velocity, based on the distribution model.
    if star_distribution == "plummer":
        stars = new_plummer_sphere(
            number_of_stars,
            convert_nbody=converter,
        )
    elif star_distribution == "king":
        stars = new_king_model(
            number_of_stars,
            star_distribution_w0,
            convert_nbody=converter,
        )
    elif star_distribution == "fractal":
        stars = new_fractal_cluster_model(
            number_of_stars,
            fractal_dimension=star_distribution_fd,
            convert_nbody=converter,
        )
    else:
        return -1, "No stellar distribution"

    # set the stellar mass.
    stars.mass = mass

    # set other stellar parameters.
    stars.metallicity = star_metallicity

    # Virialize the star cluster if > 1 star
    if len(stars) > 1:
        stars.move_to_center()
        stars.scale_to_standard(
            convert_nbody=converter,
            # virial_ratio=virial_ratio,
            # smoothing_length_squared= ...,
        )

    # Record the cluster's initial parameters to the particle distribution
    stars.collection_attributes.initial_mass_function = imf_name
    stars.collection_attributes.upper_mass_limit = upper_mass_limit
    stars.collection_attributes.lower_mass_limit = lower_mass_limit
    stars.collection_attributes.number_of_stars = number_of_stars

    stars.collection_attributes.effective_radius = effective_radius

    stars.collection_attributes.star_distribution = star_distribution
    stars.collection_attributes.star_distribution_w0 = star_distribution_w0
    stars.collection_attributes.star_distribution_fd = star_distribution_fd

    stars.collection_attributes.star_metallicity = star_metallicity

    # Derived/legacy values
    stars.collection_attributes.converter_mass = \
        converter.to_si(1 | nbody_system.mass)
    stars.collection_attributes.converter_length =\
        converter.to_si(1 | nbody_system.length)
    stars.collection_attributes.converter_speed =\
        converter.to_si(1 | nbody_system.speed)

    return stars


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
