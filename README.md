# MASC
MASC (Make a star cluster) is an AMUSE-based script to generate a star cluster, to be used in N-body simulations or for other purposes. It aims to include all popular distribution models, initial mass functions, and many more features.

## Supported features
- Stellar distributions: Plummer, King, fractal-based
- Initial mass function: Salpeter, Kroupa, equal-mass
- Metallicity
- Virial ratio

## Requirements
MASC requires AMUSE (https://github.com/amusecode/amuse) to work.
The GUI version requires QT4 or QT5 to be installed.

## Usage:
Package version:

    from amuse.units import units
    from amuse_masc import make_a_star_cluster
    cluster = make_a_star_cluster.new_cluster(
        number_of_stars=1024,
        # stellar_mass=False,
        initial_mass_function='salpeter',
        upper_mass_limit=125.0 | units.MSun,
        effective_radius=3.0 | units.parsec,
        star_distribution='plummer',
        # star_distribution_w0=7.0,
        # star_distribution_fd=2.0,
        star_metallicity=0.01,
    )

Standalone version:

    $ python make_a_star_cluster.py -h
    usage: make_a_star_cluster.py [-h] [-o CLUSTERNAME] [-type FILE TYPE]
                                  [-N NUMBER_OF_STARS] [-M CLUSTER_MASS]
                                  [-dist STAR_DISTRIBUTION]
                                  [-gasdist GAS_DISTRIBUTION]
                                  [-imf INITIAL_MASS_FUNCTION]
                                  [-modelnr CLUSTER_MODEL_NUMBER]
                                  [-kingw KING_PARAMETER_W0]
                                  [-fractald FRACTAL_PARAMETER_FD]
                                  [-uml UPPER_MASS_LIMIT] [-Z METALLICITY]
                                  [-ibf INITIAL_BINARY_FRACTION]
                                  [-R EFFECTIVE_RADIUS] [-vr VIRIAL_RATIO]
    
    optional arguments:
      -h, --help            show this help message and exit
      -o CLUSTERNAME        Outputfile name [auto]
      -type FILE TYPE       Output file type ([amuse]/ascii/starlab/nemo)
      -N NUMBER_OF_STARS    Number of stars
      -M CLUSTER_MASS       Cluster mass (takes precedence over number of stars)
      -dist STAR_DISTRIBUTION
                            Star distribution ([plummer]/king/fractal)
      -gasdist GAS_DISTRIBUTION
                            Gas distribution ([none]/plummer/king/fractal) NOT
                            IMPLEMENTED YET
      -imf INITIAL_MASS_FUNCTION
                            IMF ([kroupa]/salpeter/fixed)
      -modelnr CLUSTER_MODEL_NUMBER
                            Model number
      -kingw KING_PARAMETER_W0
                            King W0 parameter [7.0]
      -fractald FRACTAL_PARAMETER_FD
                            Fractal distribution [2.0]
      -uml UPPER_MASS_LIMIT
                            Upper mass limit (in MSun) [125.]
      -Z METALLICITY        Metallicity [0.01]
      -ibf INITIAL_BINARY_FRACTION
                            Initial binary fraction (0-1) [0.] NOT IMPLEMENTED YET
      -R EFFECTIVE_RADIUS   Cluster radius [3.0] (parsec)
      -vr VIRIAL_RATIO      Virial ration [0.5], 0.5=stable, 0.75=just expelled
                            gas, 0.1=collapsing

## COPYRIGHT/LICENSE
See LICENSE
