from setuptools import setup
from Cython.Build import cythonize

setup(
    ext_modules=cythonize("freg.pyx"),
)

#/Users/tom/opt/anaconda3/envs/tpm/bin/python