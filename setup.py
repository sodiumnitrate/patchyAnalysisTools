from setuptools import setup, find_packages
import patchyAnalysisTools

import sys
from pybind11 import get_cmake_dir
from pybind11.setup_helpers import Pybind11Extension, build_ext


ext_modules = [
	Pybind11Extension("RDF", 
    ["src/RDF.cpp"],),

    Pybind11Extension("sq",
    ["src/sq.cpp"],),
]

setup(
    name='patchyAnalysisTools',
    version=patchyAnalysisTools.__version__,
    description='simple tools to analyze output of patchy particle MC sims',
    url='https://github.com/sodiumnitrate/patchyAnalysisTools.git',
    author='Irem Altan',
    author_email='irem.altan@yale.edu',
    license='',
    packages=find_packages(),
    install_requires=['numpy','pybind11','matplotlib'],
    python_requires='>=3.6'
)
