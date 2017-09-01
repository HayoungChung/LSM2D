from distutils.core import setup, Extension

module_Opt_Module = Extension('Opt_Module', \
include_dirs = ['usr/include','home/m2do/Dropbox/0.Working/1.FEA_LSM_toPYTHON_m2do/3.mylinE/Opt_Module/src'],\
libraries = ['boost_python','nlopt','metis','ma57','blas','cblas','lapack','lapacke','tmglib','fakemetis','m'], \
library_dirs = ['usr/lib','usr/'], \
sources = ['Optimise_boost.cpp','Boundary.cpp','FastMarchingMethod.cpp','Heap.cpp',\
'Hole.cpp','InputOutput.cpp','LevelSet.cpp',\
'Mesh.cpp','Sensitivity.cpp'], \
extra_compile_args=['-std=c++11','-g','-ggdb'])

setup(name = "Opt_Module",version="0.1",description="this is a test",ext_modules =[module_Opt_Module])
