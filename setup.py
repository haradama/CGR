from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
import numpy

setup(
    name = "pycgr",
    version = "0.0.1",
    author = "Masafumi Harada",
    cmdclass = {"build_ext": build_ext},
    ext_modules = [Extension("core", ["pycgr.pyx"])],
    include_dirs = [numpy.get_include()]
)
