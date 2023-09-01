#!/usr/bin/env python
"""
Create new stars / stellar properties
"""

import logging
import numpy

from amuse.units import units
from amuse.ic.salpeter import new_salpeter_mass_distribution
from amuse.ic.brokenimf import new_kroupa_mass_distribution
from amuse.ic.flatimf import new_flat_mass_distribution


def new_masses(
    stellar_mass=False,
    initial_mass_function="kroupa",
    upper_mass_limit=125. | units.MSun,
    lower_mass_limit=0.1 | units.MSun,
    number_of_stars=None,
    exceed_mass=True,
    sort_by_mass=False,
):
    """
    Creates new stellar masses.
    """
    imf_name = initial_mass_function.lower()
    if imf_name == "salpeter":
        initial_mass_function = new_salpeter_mass_distribution
    elif imf_name == "kroupa":
        initial_mass_function = new_kroupa_mass_distribution
    elif imf_name == "flat":
        initial_mass_function = new_flat_mass_distribution
    elif imf_name == "fixed":

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
        # best underestimate mean_mass a bit for faster results
        mean_mass = max(0.25 | units.MSun, lower_mass_limit)
        mass = initial_mass_function(
            max(1, int(stellar_mass / mean_mass)),
            mass_min=lower_mass_limit,
            mass_max=upper_mass_limit,
        )
        previous_number_of_stars = len(mass)
        if exceed_mass:
            # Allow one final star to exceed stellar_mass
            final_star = 1+numpy.argmax(mass.cumsum() > stellar_mass)
            if (final_star > 1 and final_star < len(mass)):
                mass = mass[:final_star]
        else:
            # Limit to stars not exceeding stellar_mass
            mass = mass[mass.cumsum() < stellar_mass]

        additional_mass = [] | units.MSun
        while True:
            if stellar_mass < mass.sum():
                break
            if previous_number_of_stars + len(additional_mass) > len(mass):
                break
            # We don't have enough stars yet, or at least not tested this
            additional_mass = initial_mass_function(
                max(1, int(stellar_mass / mean_mass)),
                mass_min=lower_mass_limit,
                mass_max=upper_mass_limit,
            )
            if exceed_mass:
                # Allow one final star to exceed stellar_mass
                final_star = 1+numpy.argmax(
                    mass.sum() + additional_mass.cumsum() > stellar_mass
                )
                if (final_star > 1 and final_star < len(mass)):
                    additional_mass = additional_mass[:final_star]
                mass.append(additional_mass)
            else:
                # Limit to stars not exceeding stellar_mass
                additional_mass_used = additional_mass[
                    mass.sum() + additional_mass.cumsum() < stellar_mass
                ]
                mass.append(additional_mass_used)
                if len(additional_mass_used) < len(additional_mass):
                    break
        number_of_stars = len(mass)
    else:
        # Give stars their mass
        mass = initial_mass_function(
            number_of_stars,
            mass_min=lower_mass_limit,
            mass_max=upper_mass_limit,
        )

    if sort_by_mass:
        mass = mass.sorted()[::-1]
        if exceed_mass:
            final_star = 1+numpy.argmax(mass.cumsum() > stellar_mass)
            if (final_star > 1 and final_star < len(mass)):
                mass = mass[:final_star]

    return mass
