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

from amuse.io import write_set_to_file
from amuse.units import units
from amuse.support.console import set_printing_strategy
from amuse.ext.masc import new_star_cluster


def new_argument_parser():
    "Parse command line arguments"
    parser = argparse.ArgumentParser()
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
        help='Cluster mass (takes precedence over number of stars)',
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
        help='Upper mass limit (in MSun) [125.]',
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
        help='Cluster radius [3.0] (parsec)',
    )
    parser.add_argument(
        '-vr',
        dest='virial_ratio',
        type=float,
        default=0.5,
        help="Virial ration [0.5], 0.5=stable, 0.75=just expelled gas, \
                0.1=collapsing",
    )
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
    cluster_model_number = args.cluster_model_number
    star_distribution = args.star_distribution
    # gas_distribution = args.gas_distribution
    king_parameter_w0 = args.king_parameter_w0
    fractal_parameter_fd = args.fractal_parameter_fd
    initial_mass_function = args.initial_mass_function
    number_of_stars = args.number_of_stars
    if args.cluster_mass != 0:
        cluster_mass = args.cluster_mass | units.MSun
    else:
        cluster_mass = False
    upper_mass_limit = args.upper_mass_limit | units.MSun
    effective_radius = args.effective_radius | units.parsec
    metallicity = args.metallicity
    # virial_ratio = args.virial_ratio
    filetype = args.filetype

    # not implemented yet
    # initial_binary_fraction = args.initial_binary_fraction

    np.random.seed(cluster_model_number)

    if not (number_of_stars or cluster_mass):
        print("no number of stars or cluster mass given, exiting")
        exit()

    stars = new_star_cluster(
        stellar_mass=cluster_mass,
        initial_mass_function=initial_mass_function,
        upper_mass_limit=upper_mass_limit,
        number_of_stars=number_of_stars,
        effective_radius=effective_radius,
        star_distribution=star_distribution,
        star_distribution_w0=king_parameter_w0,
        star_distribution_fd=fractal_parameter_fd,
        star_metallicity=metallicity,
    )

    print(stars.mass.sum())
    print(len(stars))

    if args.clustername != "auto":
        clustername = args.clustername
    else:
        cluster_file_exists = True
        N = -1
        while cluster_file_exists:
            N += 1
            clustername = (
                clustertemplate % N
                + "." + filetype
            )
            cluster_file_exists = os.path.isfile(clustername)

    write_set_to_file(stars, clustername, filetype)


if __name__ in ["__main__"]:
    main()
