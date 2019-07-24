from __future__ import print_function

# from math import cos, sin, pi
# from amuse.ext.orbital_elements import (
#         new_binary_from_orbital_elements,
#         # orbital_elements_from_binary,
#         )

from amuse.lab import (
        units, constants,
        Particles,
        # new_salpeter_mass_distribution,
        write_set_to_file,
        )
from amuse.ext.orbital_elements import (
        generate_binaries,
        )

import numpy as np
from numpy import pi
from numpy.random import (
        random,
        uniform,
        )


def circular_velocity(
        primary,
        secondary,
        ):
        v_circ = (
            constants.G * (primary.mass + secondary.mass)
            / (secondary.position-primary.position).length()
            )**0.5
        return v_circ


def orbital_period_to_semi_major_axis(
        period, mass1, mass2=0 | units.MSun, G=constants.G):
    """
    returns semi-major axis for given period and masses
    """
    mu = G * (mass1 + mass2)
    a3 = (period/(2*np.pi))**2 * mu
    return a3**(1./3)


def new_binary_distribution(
        primary_mass,
        secondary_mass=None,
        binaries=None,
        min_mass=0.08 | units.MSun,
        ):
    """
    Takes primary masses, and returns a set of stars and a set of binaries
    formed by these stars.
    Secondary masses are given by a uniform random mass ratio with the primary
    masses. If optional secondary masses are given, these are used instead.
    binaries is an optional particleset used for the binary pairs, with given
    positions and velocities. Other parameters are ignored.
    """
    N = len(primary_mass)
    if binaries is None:
        binaries = Particles(N)
        # Should give some position/velocity as well?
    elif len(binaries) != N:
        print("binaries must be None or have the same lenght as primary_mass")
        return -1
    if secondary_mass is None:
        # Now, we need to specify the mass ratio in the binaries.
        # A flat distribution seems to be OK.
        mass_ratio = uniform(N)
        # This gives us the secondaries' masses
        secondary_mass = mass_ratio * primary_mass
        # secondaries are min_mass at least
        secondary_mass = np.maximum(secondary_mass, min_mass)
    elif len(secondary_mass) != N:
        print("Number of secondaries is unequal to number of primaries!")
        return -1
    else:
        # Make sure primary_mass is the larger of the two, and secondary_mass
        # the smaller.
        pm = primary_mass.maximum(secondary_mass)
        sm = primary_mass.minimum(secondary_mass)
        primary_mass = pm
        secondary_mass = sm
        del(pm, sm)
    # Now, we need to calculate the semi-major axes for the binaries. Since the
    # observed quantity is orbital periods, we start from there.
    mean_log_orbital_period = 5  # 10log of the period in days, (Duchene&Kraus)
    sigma_log_orbital_period = 2.3
    orbital_period = np.random.lognormal(
            size=N,
            mean=np.log(10) * mean_log_orbital_period,
            sigma=np.log(10) * sigma_log_orbital_period,
            ) | units.day
    # We need the masses to calculate the corresponding semi-major axes.
    semi_major_axis = orbital_period_to_semi_major_axis(
            orbital_period,
            primary_mass,
            secondary_mass,
            )
    # Eccentricity: square root of random value
    eccentricity = np.sqrt(random(N))
    # Other orbital elements at random
    inclination = pi * random(N) | units.rad
    true_anomaly = 2 * pi * random(N) | units.rad
    longitude_of_the_ascending_node = 2 * pi * random(N) | units.rad
    argument_of_periapsis = 2 * pi * random(N) | units.rad
    primaries, secondaries = generate_binaries(
            primary_mass,
            secondary_mass,
            semi_major_axis,
            eccentricity=eccentricity,
            inclination=inclination,
            true_anomaly=true_anomaly,
            longitude_of_the_ascending_node=longitude_of_the_ascending_node,
            argument_of_periapsis=argument_of_periapsis,
            G=constants.G,
            )
    stars = Particles()

    primaries.position += binaries.position
    secondaries.position += binaries.position
    primaries.velocity += binaries.velocity
    secondaries.velocity += binaries.velocity
    primaries = stars.add_particles(primaries)
    secondaries = stars.add_particles(secondaries)

    binaries.eccentricity = eccentricity
    binaries.semi_major_axis = semi_major_axis
    for i in range(len(primaries)):
        binaries[i].child1 = primaries[i]
        binaries[i].child2 = secondaries[i]
        # Probably needed
        binaries[i].mass = primaries[i].mass + secondaries[i].mass
    return stars, binaries


if __name__ == "__main__":
    N = 4
    stars, binaries = new_binary_distribution(N, model="default")
    write_set_to_file(stars, "stars.hdf5", "amuse")
    write_set_to_file(binaries, "stars-binarypairs.hdf5", "amuse")
