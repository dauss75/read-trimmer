from distutils.core import setup
from Cython.Build import cythonize

setup(name="Fastq Utils",
      ext_modules=cythonize("trimmer/_utils.pyx"))
