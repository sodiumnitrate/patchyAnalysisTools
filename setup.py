from setuptools import setup, find_packages
import patchyAnalysisTools

setup(
    name='patchyAnalysisTools',
    version=patchyAnalysisTools.__version__,
    description='simple tools to analyze output of patchy particle MC sims',
    url='https://github.com/sodiumnitrate/patchyAnalysisTools.git',
    author='Irem Altan',
    author_email='irem.altan@yale.edu',
    license='',
    packages=find_packages(),
    install_requires=['numpy'],
    python_requires='>=3.6'
)
