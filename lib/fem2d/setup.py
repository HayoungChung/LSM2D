from distutils.core import setup, Extension
from Cython.Build import cythonize
import numpy as np

sources = [
    'fem2d/fem_solver.cpp',
    'fem2d/fem2d.pyx'
]

setup(
    ext_modules = cythonize(Extension(
        "fem2d.fem2d", sources=sources,
        language="c++", extra_compile_args=['-std=c++11'],
        include_dirs=[np.get_include()]
    )),
    packages=[
        'fem2d',
        'fem2d.openmdao',
        'fem2d.utils',
    ]
)
