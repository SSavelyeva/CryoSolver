# -*- coding: utf-8 -*-
"""
Setup file for CryoSolver
"""

import setuptools
from codecs import open
from os import path

here = path.abspath(path.dirname(__file__))

with open(path.join(here, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()
	
setuptools.setup(
    name="CryoSolver", 
    version="1.0.0",
    author="Sofiya Savelyeva, Steffen Klöppel",
    author_email="sofiya.savelyeva@tu-dresden.de, steffen.kloeppel@tu-dresden.de",
    description="Package for cryogenic cycle simulation",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/SSavelyeva/CryoSolver ",
    download_url="https://github.com/SSavelyeva/CryoSolver/archive/v1.0.0.tar.gz",
    packages=setuptools.find_packages(exclude=['examples']),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: BSD-3-Clause",
        "Operating System :: Windows",
        "Topic :: Scientific/Engineering"
    ],
    install_requires=[
    'ctREFPROP'
    ],
    python_requires='>=3.6',
    zip_safe=True,
	license='BSD-3-Clause',
    keywords=['cryogenics', 'refrigeration', 'cycle', 'simulation'],
	setup_requires=['setuptools>=38.6.0'],
)
