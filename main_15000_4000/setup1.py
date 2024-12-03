from setuptools import setup
from Cython.Build import cythonize
setup(name="FasterSim", ext_modules=cythonize("fastest_cython.pyx"))
