#!/usr/bin/env python

"""
setup.py for linkagemodel
"""
from numpy.distutils.core import setup

setup(
    name="linkagemodel",
    version="0.1",
    author="Ian Howson",
    author_email="",
    description=("Model that links Underworld and Badlands"),
    long_description=open('README.md').read(),
    classifiers=[
        "Development Status :: 3 - Alpha",
    ],
    packages=['linkagemodel'],
    scripts=[],
)
