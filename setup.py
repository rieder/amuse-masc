from setuptools import setup

version = "0.4.0"
name = 'amuse-masc'
author = 'Steven Rieder'
author_email = 'steven@rieder.nl'
license_ = "MIT"
url = 'http://amusecode.org'

classifiers = [
    # Python versions supported by amuse-masc
    "Programming Language :: Python :: 3",
    "Programming Language :: Python :: 3.5",
    "Programming Language :: Python :: 3.6",
    "Programming Language :: Python :: 3.7",
    "Programming Language :: Python :: 3.8",

    # License
    "License :: OSI Approved :: MIT License",

    # OS support
    "Operating System :: OS Independent",

    # Maturity of amuse_masc
    "Development Status :: 4 - Beta",

    # Intended audience
    "Intended Audience :: Science/Research",
]

install_requires = [
    'wheel>=0.32',
    'amuse-framework>=13.0.0',
]
description = 'Make a star cluster'
with open("README.md", "r") as fh:
    long_description = fh.read()
long_description_content_type = "text/markdown"

extensions = []

all_data_files = []

packages = ['amuse.ext.masc']

package_dir = {
    'amuse.ext.masc': 'src/amuse/ext/masc'
}

package_data = {
}

setup(
    name=name,
    version=version,
    classifiers=classifiers,
    url=url,
    author_email=author_email,
    author=author,
    license=license_,
    description=description,
    long_description=long_description,
    long_description_content_type=long_description_content_type,
    install_requires=install_requires,
    ext_modules=extensions,
    package_dir=package_dir,
    packages=packages,
    package_data=package_data,
    data_files=all_data_files,
    python_requires='>=3.5, <4',
    scripts=["masc.py"],
)
