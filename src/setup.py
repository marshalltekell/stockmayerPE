#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 13 15:23:14 2020

@author: marshalltekell


from http://www.swig.org/Doc1.3/Python.html#Python_nn2
"""

from distutils.core import setup, Extension
import numpy as np


# Obtain the numpy include directory.  This logic works across numpy versions.
try:
    numpy_include = np.get_include()
except AttributeError:
    numpy_include = np.get_numpy_include()


cfunctions_module = Extension('_cfunctions',
                           sources=['cfunctions_wrap.c', 'cfunctions.c'],
                           include_dirs = [numpy_include],
                           )

setup (name = 'cfunctions',
       version = '0.1',
       author      = "Marshall Tekell",
       description = """Analysis functions for molecular dynamics simulation of Stockmayer polymer electrolyte""",
       ext_modules = [cfunctions_module],
       py_modules = ["cfunctions"],
       )