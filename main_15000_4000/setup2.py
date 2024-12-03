from setuptools import setup
from Cython.Build import cythonize
setup(name="FullSim", ext_modules=cythonize("more_cythonized.pyx"))
