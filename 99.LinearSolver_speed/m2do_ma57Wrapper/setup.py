from distutils.core import setup, Extension
from Cython.Build import cythonize
import numpy as np

sources = [
    'ma57_solver.cpp',
    'ma57_pyWrap.pyx'
]

setup(
    ext_modules = cythonize(Extension(
        "ma57_pyWrap", sources=sources,
        language="c++", extra_compile_args=['-std=c++11'],
        include_dirs=[np.get_include()],
        libraries = ['lapack','ma57','gfortran'],   
        library_dirs=['/usr/local/lib'],#,'/usr/lib/x86_64-linux-gnu/'],
        # runtime_library_dirs=['/usr/local/lib','/usr/lib/x86_64-linux-gnu/'],
        extra_objects=['/usr/lib/libmetis.a']
    )),
)
