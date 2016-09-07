#!/usr/bin/env python

"""
make_a_star_cluster.py creates a model star cluster,
which can then be used in N-body simulations or for other purposes.

It requires AMUSE, which can be downloaded from http://amusecode.org or
https://github.com/amusecode/amuse.
"""

import os,sys
import numpy as np

from amuse.lab import *
from amuse.io import write_set_to_file
from amuse.units import units
from amuse.units import generic_unit_converter
from amuse.support.console import set_printing_strategy

import argparse

def new_argument_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument(
            '-o',
            dest='clustername',
            default="auto",
            help="Outputfile name [auto]",
            )
    parser.add_argument(
            '-type',
            dest='file type',
            default="amuse",
            help="Output file type ([amuse]/ascii/starlab/nemo)",
            )
    parser.add_argument(
            '-N', 
            dest='number_of_stars', 
            type=int, 
            help='Number of stars',
            )
    parser.add_argument(
            '-M', 
            dest='cluster_mass', 
            type=float, 
            default = 0,
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
            help="Gas distribution ([none]/plummer/king/fractal) NOT IMPLEMENTED YET",
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
            help='Virial ration [0.5], 0.5=stable, 0.75=just expelled gas, 0.1=collapsing',
            )
    args = parser.parse_args()
    return args

def new_binary(
        mass1, 
        mass2, 
        semi_major_axis,       
        eccentricity = 0, 
        ):

    mass_fraction_particle_1 = mass1 / (mass1 + mass2)

    binary = datamodel.Particles(2)
        
    binary[0].mass = mass1
    binary[1].mass = mass2

    mu = constants.G * total_mass

    velocity_perihelion = np.sqrt( 
            mu / semi_major_axis  * 
            ( 
                (1.0 + eccentricity) / 
                (1.0 - eccentricity)
                )
            )
    radius_perihelion = (
            semi_major_axis * 
            (1.0 - eccentricity)
            )
    
    binary[0].position = (
            (1.0 - mass_fraction_particle_1) * 
            radius_perihelion * 
            [1.0,0.0,0.0]
            )
    binary[1].position = -(
            mass_fraction_particle_1 * 
            radius_perihelion * 
            [1.0,0.0,0.0]
            )

    binary[0].velocity = (
            (1.0 - mass_fraction_particle_1) * 
            velocity_perihelion * 
            [0.0,1.0,0.0]
            )
    binary[1].velocity = -(
            mass_fraction_particle_1 * 
            velocity_perihelion * 
            [0.0,1.0,0.0]
            )

    return binary
    

def create_binaries(
        center_of_mass_particles,
        ):
    singles_in_binaries = datamodel.Particles()
    for i in range(len(center_of_mass_particles)):
        binary = center_of_mass_particles[i]
        particles_in_binary = self.new_binary(
            binary.binary_mass1,
            binary.binary_mass2,
            binary.binary_semi_major_axis,
            binary.binary_eccentricity,
        )
            
        particles_in_binary.radius = binary.binary_semi_major_axis
        
        binary.child1 = particles_in_binary[0]
        binary.child2 = particles_in_binary[1]
        binary.mass = binary.binary_mass1 + binary.binary_mass2
       
        particles_in_binary.position += binary.position
        particles_in_binary.velocity += binary.velocity
        singles_in_binaries.add_particles(particles_in_binary)
        ##FIXME Rotate binary around phi, theta axes

    return center_of_mass_particles, singles_in_binaries


def setup_binaries(self):
    center_of_mass_particles = new_plummer_sphere(self.N_binaries, self.converter)
    ## Generate semi-major axes based on models
    center_of_mass_particles.binary_semi_major_axis = 30. | units.AU 
    ## Generate eccentricities based on models
    center_of_mass_particles.binary_eccentricity = 0.0
    ## Generate angles at random?
    center_of_mass_particles.binary_theta = np.random.random()*np.pi
    center_of_mass_particles.binary_phi   = np.random.random()*np.pi*2


    center_of_mass_particles.binary_mass1 = new_salpeter_mass_distribution(self.N_binaries)
    center_of_mass_particles.binary_mass2 = new_salpeter_mass_distribution(self.N_binaries)
    center_of_mass_particles.mass = center_of_mass_particles.binary_mass1 + center_of_mass_particles.binary_mass2
    center_of_mass_particles.radius = self.interaction_radius
    self.binaries, self.singles_in_binaries = self.create_binaries(
            center_of_mass_particles,
            )

    
if __name__ in ["__main__"]:
    set_printing_strategy(
            "custom", 
            preferred_units = [units.MSun, units.parsec, units.yr, units.kms], 
            precision       = 5
            )
    clustertemplate = "TESTCluster_%08i.hdf5"

    args = new_argument_parser()
    cluster_model_number    = args.cluster_model_number
    star_distribution       = args.star_distribution
    gas_distribution        = args.gas_distribution
    king_parameter_w0       = args.king_parameter_w0
    fractal_parameter_fd    = args.fractal_parameter_fd
    initial_mass_function   = args.initial_mass_function
    number_of_stars         = args.number_of_stars
    if args.cluster_mass != 0:
        cluster_mass        = args.cluster_mass | units.MSun
    else:
        cluster_mass        = False
    upper_mass_limit        = args.upper_mass_limit | units.MSun
    effective_radius        = args.effective_radius | units.parsec
    metallicity             = args.metallicity
    virial_ratio            = args.virial_ratio
    filetype                = args.filetype

    ## not implemented yet
    initial_binary_fraction = args.initial_binary_fraction

    np.random.seed(cluster_model_number)

    if not (number_of_stars or cluster_mass):
        print "no number of stars or cluster mass given, exiting"
        exit()

    if cluster_mass:
        ## Add stars to cluster, until mass limit reached (inclusive!)
        if initial_mass_function == "kroupa":
            from amuse.ic.brokenimf import new_kroupa_mass_distribution
            mass = new_kroupa_mass_distribution(0)
            while mass.sum() < cluster_mass:
                mass.append(
                        new_kroupa_mass_distribution(
                            1,
                            mass_max = upper_mass_limit,
                            )
                        )
            total_mass = mass.sum()
            number_of_stars = len(mass)
        elif initial_mass_function == "salpeter":
            from amuse.ic.salpeter import new_salpeter_mass_distribution
            mass = new_salpeter_mass_distribution(0)
            while mass.sum() < cluster_mass:
                mass.append(
                        new_salpeter_mass_distribution(
                            1,
                            mass_max = upper_mass_limit,
                            )
                        )
            total_mass = mass.sum()
            number_of_stars = len(mass)
        elif initial_mass_function == "fixed":
            mass_of_each_star = cluster_mass / number_of_stars
            mass = mass_of_each_star
            total_mass = cluster_mass
        else:
            print "No mass function"
            exit()
    else:
        ## Give stars their mass       
        if initial_mass_function == "kroupa":
            from amuse.ic.brokenimf import new_kroupa_mass_distribution
            mass = new_kroupa_mass_distribution(
                    number_of_stars,
                    mass_max = upper_mass_limit,
                    )
            total_mass = mass.sum()
        elif initial_mass_function == "salpeter":
            from amuse.ic.salpeter import new_salpeter_mass_distribution
            mass = new_salpeter_mass_distribution(
                    number_of_stars,
                    mass_max = upper_mass_limit,
                    )
            total_mass = mass.sum()
        elif initial_mass_function == "fixed":
            mass = mass_of_each_star
            total_mass = number_of_stars * mass
        else:
            print "No mass function"
            exit()

    converter = generic_unit_converter.ConvertBetweenGenericAndSiUnits(
            total_mass,
            1|units.kms,
            effective_radius,
            )

    ## Distribute stars (position + velocity)
    if star_distribution == "plummer":
        from amuse.ic.plummer import new_plummer_sphere
        stars = new_plummer_sphere(
                number_of_stars, 
                convert_nbody = converter,
                )
    elif star_distribution == "king":
        from amuse.ic.kingmodel import new_king_model
        stars = new_king_model(
                number_of_stars,
                king_parameter_w0,
                convert_nbody = converter,
                )
    elif star_distribution == "fractal":
        from amuse.ic.fractalcluster import new_fractal_cluster_model
        stars = new_fractal_cluster_model(
                number_of_stars,
                fractal_dimension = fractal_parameter_fd,
                convert_nbody = converter,
                )
    else:
        print "No stellar distribution"
        exit()

    ## set mass of the stars
    stars.mass = mass

    ## Other parameters
    stars.metallicity = metallicity

    ## Virialize cluster!
    stars.move_to_center()
    stars.scale_to_standard(
            convert_nbody = converter,
            #virial_ratio = virial_ratio,
            #smoothing_length_squared = ...,
            )

    ## Write initial parameters to distribution
    stars.collection_attributes.cluster_model_number    = cluster_model_number
    stars.collection_attributes.star_distribution       = star_distribution
    stars.collection_attributes.gas_distribution        = gas_distribution
    stars.collection_attributes.king_parameter_w0       = king_parameter_w0
    stars.collection_attributes.fractal_parameter_fd    = fractal_parameter_fd
    stars.collection_attributes.initial_mass_function   = initial_mass_function
    stars.collection_attributes.number_of_stars         = number_of_stars
    stars.collection_attributes.upper_mass_limit        = upper_mass_limit
    stars.collection_attributes.effective_radius        = effective_radius
    stars.collection_attributes.metallicity             = metallicity
    stars.collection_attributes.initial_binary_fraction = initial_binary_fraction
    stars.collection_attributes.virial_ratio            = virial_ratio

    ## Derived/legacy values
    stars.collection_attributes.Mcluster = stars.mass.sum()
    stars.collection_attributes.Rcluster = effective_radius
    stars.collection_attributes.N        = len(stars)
    stars.collection_attributes.seed     = cluster_model_number

    print stars.mass.sum()
    print len(stars)

    if args.clustername != "auto":
        clustername = args.clustername
    else:
        cluster_file_exists = True
        N=-1
        while cluster_file_exists:
            N+=1
            clustername = clustertemplate%N
            cluster_file_exists = os.path.isfile(clustername)

    write_set_to_file(stars,clustername,filetype)
