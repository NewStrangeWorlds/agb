
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
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_CXX_FLAGS "-Wall -pedantic -MMD")
set(CMAKE_C_FLAGS "-Wall -pedantic -MMD")

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
  src/transport_coefficients/transport_coeff.cpp)
  


#check for OpenMP
find_package(OpenMP REQUIRED)


include(FetchContent)

#download the FastChem 3.1 if necessary
FetchContent_Declare(
  fastchem
  GIT_REPOSITORY https://github.com/exoclime/fastchem
  GIT_TAG        e8dbeb7ec7719dbd0021d672a0b5d330e924a764
)

FetchContent_GetProperties(fastchem)
set(FASTCHEM_ONLY_LIBRARY ON)

if (NOT fastchem_POPULATED)
  FetchContent_Populate(fastchem)
  add_subdirectory(${fastchem_SOURCE_DIR} ${fastchem_BINARY_DIR})
endif()


# build application
add_executable(agb_wind ${SRC})
target_include_directories(agb_wind PUBLIC)

#link files for the final helios-r code
target_link_libraries(agb_wind PUBLIC fastchem_lib OpenMP::OpenMP_CXX)
