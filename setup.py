#!/usr/bin/env python3
from setuptools import setup


name = 'amuse-masc'
author = 'Steven Rieder'
author_email = 'steven+masc@rieder.nl'
license_ = "Apache License 2.0"
url = 'https://github.com/rieder/masc'


install_requires = [
    'wheel>=0.32',
    'amuse-framework>=2022.6.0',
]
setup_requires = []
description = 'AMUSE tool to create star cluster initial conditions'
with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()
long_description_content_type = "text/markdown"

extensions = []

all_data_files = []

packages = [
    'amuse.ext.masc',
]

package_data = {
}

classifiers = [
    # Maturity of amuse_masc
    "Development Status :: 4 - Beta",
    # Intended audience
    "Intended Audience :: Science/Research",
    # License
    "License :: OSI Approved :: Apache Software License",
    # OS support
    "Operating System :: OS Independent",
    # Python versions supported by amuse-masc
    "Programming Language :: Python :: 3",
    "Programming Language :: Python :: 3.7",
    "Programming Language :: Python :: 3.8",
    "Programming Language :: Python :: 3.9",
    "Programming Language :: Python :: 3.10",
    'Topic :: Scientific/Engineering :: Astronomy',    
]

try:
    from src.amuse.ext.masc.version import version
    use_scm_version = False
except ImportError:
    version = False
    setup_requires += ['setuptools_scm',]
    use_scm_version = {
        "root": ".",
        "relative_to": __file__,
        "write_to": "src/amuse/ext/masc/version.py",
    }

setup(
    name=name,
    use_scm_version=use_scm_version,
    setup_requires=setup_requires,
    version=version,
    classifiers=classifiers,
    url=url,
    project_urls={
        "Bug Tracker": "https://github.com/rieder/masc/issues",
    },
    author_email=author_email,
    author=author,
    license=license_,
    description=description,
    long_description=long_description,
    long_description_content_type=long_description_content_type,
    install_requires=install_requires,
    ext_modules=extensions,
    package_dir={
        'amuse.ext.masc': 'src/amuse/ext/masc',
    },
    packages=packages,
    package_data=package_data,
    data_files=all_data_files,
    scripts=["bin/masc.py", ],
    python_requires=">=3.7, <4",
)
