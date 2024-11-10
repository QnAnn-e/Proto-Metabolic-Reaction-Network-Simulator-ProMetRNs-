from setuptools import setup
from Cython.Build import cythonize
setup(name="FullSim", ext_modules=cythonize("cythonized_faster_sim.pyx"))
