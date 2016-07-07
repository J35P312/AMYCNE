from distutils.core import setup
from Cython.Build import cythonize

setup(
  name = 'AMYCNE',
  ext_modules = cythonize("common.py"),
)
