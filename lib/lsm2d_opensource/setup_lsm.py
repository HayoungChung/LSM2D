from distutils.core import setup, Extension
from Cython.Build import cythonize
import numpy as np

sources = [
    'LSM2d/LS_toolbox.cpp',
    'LSM2d/LSM2d.pyx',
    'M2DO_LSM/src/hole.cpp',
    'M2DO_LSM/src/mesh.cpp',
    'M2DO_LSM/src/boundary.cpp',
    'M2DO_LSM/src/fast_marching_method.cpp',
    'M2DO_LSM/src/heap.cpp',
    'M2DO_LSM/src/level_set.cpp',
    'M2DO_LSM/src/optimise.cpp',
]

setup(
    ext_modules = cythonize(Extension(
        "LSM2d", sources=sources,
        language="c++", extra_compile_args=['-std=c++11'],
        include_dirs = [np.get_include(),"./M2DO_LSM/include"],
	    library_dirs = ['usr/lib','usr/',"./M2DO_LSM/include"], 
    )),
)
