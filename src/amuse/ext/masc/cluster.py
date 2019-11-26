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

from amuse.units import (
    units,
    nbody_system,
    generic_unit_converter,
)


def new_star_cluster(
        stellar_mass=False,
        initial_mass_function="salpeter",
        upper_mass_limit=125. | units.MSun,
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

    if stellar_mass:
        # Add stars to cluster, until mass limit reached (inclusive!)
        if initial_mass_function == "kroupa":
            from amuse.ic.brokenimf import new_kroupa_mass_distribution
            mass = new_kroupa_mass_distribution(0)
            while mass.sum() < stellar_mass:
                mass.append(
                    new_kroupa_mass_distribution(
                        1,
                        mass_max=upper_mass_limit,
                    )[0]
                )
            total_mass = mass.sum()
            number_of_stars = len(mass)
        elif initial_mass_function == "salpeter":
            from amuse.ic.salpeter import new_salpeter_mass_distribution
            mass = new_salpeter_mass_distribution(0)
            while mass.sum() < stellar_mass:
                mass.append(
                    new_salpeter_mass_distribution(
                        1,
                        mass_max=upper_mass_limit,
                    )[0]
                )
            total_mass = mass.sum()
            number_of_stars = len(mass)
        elif initial_mass_function == "fixed":
            mass_of_individual_star = stellar_mass / number_of_stars
            mass = mass_of_individual_star
            total_mass = stellar_mass
        else:
            return -1, "No mass function"
    else:
        # Give stars their mass
        if initial_mass_function == "kroupa":
            from amuse.ic.brokenimf import new_kroupa_mass_distribution
            mass = new_kroupa_mass_distribution(
                number_of_stars,
                mass_max=upper_mass_limit,
            )
            total_mass = mass.sum()
        elif initial_mass_function == "salpeter":
            from amuse.ic.salpeter import new_salpeter_mass_distribution
            mass = new_salpeter_mass_distribution(
                number_of_stars,
                mass_max=upper_mass_limit,
            )
            total_mass = mass.sum()
        elif initial_mass_function == "fixed":
            mass = mass_of_individual_star
            total_mass = number_of_stars * mass
        else:
            return -1, "No mass function"

    converter = generic_unit_converter.ConvertBetweenGenericAndSiUnits(
        total_mass,
        1. | units.kms,
        effective_radius,
    )
    # Give stars a position and velocity, based on the distribution model.
    if star_distribution == "plummer":
        from amuse.ic.plummer import new_plummer_sphere
        stars = new_plummer_sphere(
            number_of_stars,
            convert_nbody=converter,
        )
    elif star_distribution == "king":
        from amuse.ic.kingmodel import new_king_model
        stars = new_king_model(
            number_of_stars,
            star_distribution_w0,
            convert_nbody=converter,
        )
    elif star_distribution == "fractal":
        from amuse.ic.fractalcluster import new_fractal_cluster_model
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
    stars.collection_attributes.initial_mass_function = initial_mass_function
    stars.collection_attributes.upper_mass_limit = upper_mass_limit
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
