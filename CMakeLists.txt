
# to build with cmake
# create a build directory and move into it
# $ mkdir build
# $ cd build
# generate the makefile (to do only ones, if we don't add files or change makefiles)
# don't forget the two points at the end of the command '..'.
# It runs cmake in the 'build' directory
# but with the data from the '..' directory.
# $ cmake ..

cmake_minimum_required (VERSION 3.10 FATAL_ERROR)

# check command line args, user can specify SM to override manual setting
set(SM "0" CACHE STRING "GPU SM value")


project(helios-r C CXX)


#some C++/C flags
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_CXX_FLAGS "-Wall -pedantic -MMD")
set(CMAKE_C_FLAGS "-Wall -pedantic -MMD")

#Release codegen: keep the default -O3 -DNDEBUG and add loop unrolling for the
#tight radiative-transfer kernels. Do NOT add -march=native here: FastChem's
#Eigen breaks under it (see project memory). -ffast-math is also intentionally
#omitted to keep results bit-stable.
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG -funroll-loops")
set(CMAKE_C_FLAGS_RELEASE "-O3 -DNDEBUG -funroll-loops")

set(CMAKE_BUILD_TYPE Release)

#output directories
set(CMAKE_BINARY_DIR ${CMAKE_SOURCE_DIR})
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_DEBUG ${CMAKE_BINARY_DIR})
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_RELEASE ${CMAKE_BINARY_DIR})



#define the source files 
set(SRC
  src/model_main/model_main.cpp
  src/agb_model/agb_model.cpp
  src/config/config.cpp
  src/additional/aux_functions.cpp
  src/spectral_grid/spectral_grid.cpp
  src/spectral_grid/convert.cpp
  src/spectral_grid/spectral_grid_interpolate.cpp
  src/atmosphere/atmosphere.cpp
  src/atmosphere/atmosphere_read_write.cpp
  src/chemistry/fastchem_chemistry.cpp
  src/transport_coefficients/cross_section_file.cpp
  src/transport_coefficients/opacity_species.cpp
  src/transport_coefficients/sampled_data.cpp
  src/transport_coefficients/species_rayleigh_cross_sections.cpp
  src/transport_coefficients/transport_coeff.cpp
  src/transport_coefficients/h_m.cpp
  src/transport_coefficients/h2_m.cpp
  src/radiative_transfer/radiative_transfer.cpp
  src/radiative_transfer/output.cpp
  src/radiative_transfer/impact_parameter.cpp
  src/radiative_transfer/radiation_field.cpp
  src/radiative_transfer/moment_system.cpp
  src/temperature_correction/temperature_correction.cpp
  src/dust/dust_species.cpp
  src/dust/analytic_dust.cpp
  src/dust/gail_sedlmayr_dust.cpp
  src/dust/gail_sedlmayr_thermo_functions.cpp
  src/hydrodynamics/hydrodynamics.cpp
  src/hydrodynamics/structure_solver.cpp
  _deps/lx_mie-src/mie/mie.cpp)
  


#check for OpenMP
find_package(OpenMP REQUIRED)


include(FetchContent)

#download FastChem 4.0.2 if necessary
FetchContent_Declare(
  fastchem
  GIT_REPOSITORY https://github.com/exoclime/fastchem
  GIT_TAG        v4.0.2
)


FetchContent_GetProperties(fastchem)
set(FASTCHEM_ONLY_LIBRARY ON)

if (NOT fastchem_POPULATED)
  FetchContent_Populate(fastchem)
  add_subdirectory(${fastchem_SOURCE_DIR} ${fastchem_BINARY_DIR})
endif()



#download LX-MIE if necessary
FetchContent_Declare(
  lx_mie
  GIT_REPOSITORY https://github.com/daniel-kitzmann/LX-MIE
  GIT_TAG        bf62d51017ad46dd82218274c2d1b571d2ccae5f
)

FetchContent_GetProperties(lx_mie)
if (NOT lx_mie_POPULATED)
  FetchContent_Populate(lx_mie)
endif()



#download Eigen (header-only) for the linear algebra in the Henyey structure solver
FetchContent_Declare(
  eigen
  GIT_REPOSITORY https://gitlab.com/libeigen/eigen
  GIT_TAG        3.4.0
)

FetchContent_GetProperties(eigen)
if (NOT eigen_POPULATED)
  FetchContent_Populate(eigen)
endif()



#build CppAD (provides autodiff for the Henyey structure-solver Jacobian).
#CppAD's CMake generates configure.hpp and forces CMAKE_INSTALL_PREFIX, so it is
#built/installed in isolation via ExternalProject rather than add_subdirectory.
include(ExternalProject)
set(CPPAD_INSTALL_DIR ${CMAKE_BINARY_DIR}/cppad-install)

ExternalProject_Add(cppad_external
  GIT_REPOSITORY https://github.com/coin-or/CppAD
  GIT_TAG        20250000.3
  GIT_SHALLOW    TRUE
  UPDATE_DISCONNECTED TRUE
  CMAKE_ARGS
    -Dcppad_prefix=${CPPAD_INSTALL_DIR}
    -Dcppad_testvector=std
    -DCMAKE_BUILD_TYPE=Release
  BUILD_BYPRODUCTS ${CPPAD_INSTALL_DIR}/lib/libcppad_lib.so
)

#the include dir must exist at configure time for target_include_directories
file(MAKE_DIRECTORY ${CPPAD_INSTALL_DIR}/include)



# build application
add_executable(agb_wind ${SRC})
add_dependencies(agb_wind cppad_external)
target_include_directories(agb_wind PUBLIC ${CPPAD_INSTALL_DIR}/include ${eigen_SOURCE_DIR})

#cppad_lib is a shared library in a non-standard location; add it to the rpath
#so the runtime loader finds it (the executable is placed in the source dir)
set_target_properties(agb_wind PROPERTIES BUILD_RPATH ${CPPAD_INSTALL_DIR}/lib)

#link files for the final helios-r code
target_link_libraries(agb_wind PUBLIC fastchem_lib OpenMP::OpenMP_CXX ${CPPAD_INSTALL_DIR}/lib/libcppad_lib.so)
