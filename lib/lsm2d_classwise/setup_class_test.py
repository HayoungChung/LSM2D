from distutils.core import setup, Extension
from Cython.Build import cythonize
import numpy as np

sources = [
    'lsm_classes.pyx',
    'M2DO_LSM/src/hole.cpp',
    'M2DO_LSM/src/mesh.cpp',
    'M2DO_LSM/src/level_set.cpp',
    'M2DO_LSM/src/boundary.cpp',
    'M2DO_LSM/src/fast_marching_method.cpp',
    'M2DO_LSM/src/heap.cpp',
    'M2DO_LSM/src/optimise_noNLOPT.cpp',
]

setup(
    ext_modules = cythonize(Extension(
        "lsm_classes", sources=sources,
        language="c++", extra_compile_args=['-std=c++11','-g'],
        include_dirs = [np.get_include(),"./M2DO_LSM/include"],
	    library_dirs = ['usr/lib','usr/',"./M2DO_LSM/include"], 
	extra_link_args=["-g"],
    ), gdb_debug=True),
)
