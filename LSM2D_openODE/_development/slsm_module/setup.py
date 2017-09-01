from distutils.core import setup, Extension

module_slsm_Module = Extension('slsm_Module', \
include_dirs = ['usr/include','home/m2do/Dropbox/0.Working/1.FEA_LSM_toPYTHON_m2do/3.mylinE/slsm_module/src'],\
libraries = ['boost_python','nlopt','metis','ma57','blas','cblas','lapack','lapacke','tmglib','fakemetis','m'], \
library_dirs = ['usr/lib','usr/'], \
sources = ['slsm_Module.cpp','src/Boundary.cpp','src/FastMarchingMethod.cpp','src/Heap.cpp',\
'src/Hole.cpp','src/InputOutput.cpp','src/LevelSet.cpp',\
'src/Mesh.cpp','src/Optimise.cpp','src/Sensitivity.cpp'], \
extra_compile_args=['-std=c++11','-g','-ggdb'])

setup(name = "slsm_Module",version="0.1",description="this is a test",ext_modules =[module_slsm_Module])
