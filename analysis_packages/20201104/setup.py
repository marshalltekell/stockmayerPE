#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun  2 17:44:00 2020

@author: marshalltekell
"""

from setuptools import setup
from Cython.Build import cythonize

setup(
    ext_modules = cythonize("cython_functions_C1.pyx", annotate=True)
)