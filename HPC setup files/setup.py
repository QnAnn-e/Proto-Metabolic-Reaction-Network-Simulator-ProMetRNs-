from setuptools import setup
from Cython.Build import cythonize
setup(name="FasterSim", ext_modules=cythonize("cythonized_faster_sim.pyx"))
