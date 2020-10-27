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

import os
import argparse

import numpy as np

from amuse.io import write_set_to_file, read_set_from_file
from amuse.datamodel import Particles
from amuse.units import units
from amuse.support.console import set_printing_strategy
from amuse.ext.masc import new_star_cluster, new_stars_from_sink


def new_argument_parser():
    "Parse command line arguments"
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--sinks',
        dest='sinks',
        default=None,
        help="Generate stars from sink particles in this file (EXPERIMENTAL)",
    )
    parser.add_argument(
        '-o',
        dest='clustername',
        default="auto",
        help="Outputfile name [auto]",
    )
    parser.add_argument(
        '-type',
        dest='filetype',
        default="amuse",
        help="Output file type ([amuse]/txt/starlab/nemo)",
    )
    parser.add_argument(
        '-N',
        dest='number_of_stars',
        default=1024,
        type=int,
        help='Number of stars [1024]',
    )
    parser.add_argument(
        '-M',
        dest='cluster_mass',
        type=float,
        default=0,
        help='Cluster mass (in MSun) (takes precedence over number of stars)',
    )
    parser.add_argument(
        '-dist',
        dest='star_distribution',
        default="plummer",
        help="Star distribution ([plummer]/king/fractal)",
    )
    parser.add_argument(
        '-gasdist',
        dest='gas_distribution',
        default="none",
        help="Gas distribution ([none]/plummer/king/fractal) \
                NOT IMPLEMENTED YET",
    )
    parser.add_argument(
        '-imf',
        dest='initial_mass_function',
        default="kroupa",
        help="IMF ([kroupa]/salpeter/fixed)",
    )
    parser.add_argument(
        '-modelnr',
        dest='cluster_model_number',
        type=int,
        default=30483,
        help='Model number',
    )
    parser.add_argument(
        '-kingw',
        dest='king_parameter_w0',
        type=float,
        default=7.0,
        help='King W0 parameter [7.0]',
    )
    parser.add_argument(
        '-fractald',
        dest='fractal_parameter_fd',
        type=float,
        default=2.0,
        help='Fractal distribution [2.0]',
    )
    parser.add_argument(
        '-uml',
        dest='upper_mass_limit',
        type=float,
        default=125.,
        help='Upper mass limit for stars (in MSun) [125.]',
    )
    parser.add_argument(
        '-lml',
        dest='lower_mass_limit',
        type=float,
        default=0.1,
        help='Lower mass limit for stars (in MSun) [0.1]',
    )
    parser.add_argument(
        '-Z',
        dest='metallicity',
        type=float,
        default=0.01,
        help='Metallicity [0.01]',
    )
    parser.add_argument(
        '-ibf',
        dest='initial_binary_fraction',
        type=float,
        default=0.,
        help='Initial binary fraction (0-1) [0.] NOT IMPLEMENTED YET',
    )
    parser.add_argument(
        '-R',
        dest='effective_radius',
        type=float,
        default=3.0,
        help=(
            "Cluster radius / default radius when forming stars from sinks"
            " [3.0] (in parsec)"
        ),
    )
    parser.add_argument(
        '--vel',
        dest='velocity_dispersion',
        type=float,
        default=1.0,
        help=(
            "Default velocity dispersion when forming stars from sinks"
            " [1.0] (in km/s)"
        ),
    )
    # parser.add_argument(
    #     '-vr',
    #     dest='virial_ratio',
    #     type=float,
    #     default=0.5,
    #     help="Virial ration [0.5], 0.5=stable, 0.75=just expelled gas, \
    #             0.1=collapsing",
    # )
    args = parser.parse_args()
    return args


def main():
    "Make a star cluster"
    set_printing_strategy(
        "custom",
        preferred_units=[units.MSun, units.parsec, units.yr, units.kms],
        precision=5,
    )
    clustertemplate = "TESTCluster_%08i"

    args = new_argument_parser()
    sinks = args.sinks
    cluster_model_number = args.cluster_model_number
    star_distribution = args.star_distribution
    # gas_distribution = args.gas_distribution
    king_parameter_w0 = args.king_parameter_w0
    fractal_parameter_fd = args.fractal_parameter_fd
    initial_mass_function = args.initial_mass_function.lower()

    number_of_stars = args.number_of_stars
    if args.cluster_mass != 0:
        cluster_mass = args.cluster_mass | units.MSun
    else:
        cluster_mass = False
    upper_mass_limit = args.upper_mass_limit | units.MSun
    lower_mass_limit = args.lower_mass_limit | units.MSun
    effective_radius = args.effective_radius | units.parsec
    metallicity = args.metallicity
    # virial_ratio = args.virial_ratio
    filetype = args.filetype

    # not implemented yet
    # initial_binary_fraction = args.initial_binary_fraction

    np.random.seed(cluster_model_number)

    if not (number_of_stars or cluster_mass or sinks):
        print(
            "no number of stars, cluster mass or origin sinks given, exiting"
        )
        exit()

    if sinks is not None:
        sinks = read_set_from_file(sinks, "amuse")
        stars = Particles()
        for sink in sinks:
            try:
                velocity_dispersion = sink.u.sqrt()
            except AttributeError:
                velocity_dispersion = args.velocity_dispersion | units.kms
            new_stars = new_stars_from_sink(
                sink,
                upper_mass_limit=upper_mass_limit,
                lower_mass_limit=lower_mass_limit,
                default_radius=effective_radius,
                velocity_dispersion=velocity_dispersion,
                initial_mass_function=initial_mass_function,
                # logger=logger,
            )
            stars.add_particles(
                new_stars
            )
    else:
        stars = new_star_cluster(
            stellar_mass=cluster_mass,
            initial_mass_function=initial_mass_function,
            upper_mass_limit=upper_mass_limit,
            lower_mass_limit=lower_mass_limit,
            number_of_stars=number_of_stars,
            effective_radius=effective_radius,
            star_distribution=star_distribution,
            star_distribution_w0=king_parameter_w0,
            star_distribution_fd=fractal_parameter_fd,
            star_metallicity=metallicity,
        )

    print(
        "%i stars generated (%s)"
        % (len(stars), stars.total_mass().in_(units.MSun))
    )

    if args.clustername != "auto":
        clustertemplate = args.clustername + "%s"

    stars_file_exists = True
    sinks_file_exists = True
    N = -1
    while (stars_file_exists or sinks_file_exists):
        N += 1
        starsfilename = (
            clustertemplate % N
            + "-stars." + filetype
        )
        stars_file_exists = os.path.isfile(starsfilename)
        sinksfilename = (
            clustertemplate % N
            + "-sinks." + filetype
        )
        sinks_file_exists = os.path.isfile(sinksfilename)

    write_set_to_file(stars, starsfilename, filetype)
    if sinks is not None:
        write_set_to_file(sinks, sinksfilename, filetype)


if __name__ in ["__main__"]:
    main()
