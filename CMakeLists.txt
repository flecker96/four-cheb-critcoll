cmake_minimum_required(VERSION 3.16)

project(CritCollCheb LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# -------------------------
# Dependencies
# -------------------------

# FFTW (robust method via pkg-config)
find_package(PkgConfig REQUIRED)
pkg_check_modules(FFTW REQUIRED fftw3)

# HDF5
find_package(HDF5 REQUIRED COMPONENTS CXX)

# JSON
find_package(nlohmann_json REQUIRED)

# LAPACK
find_package(LAPACK REQUIRED)
find_library(LAPACKE_LIB lapacke)  # C interface

# -------------------------
# Source files
# -------------------------

set(SOURCES
    src/main.cpp
    src/common.cpp
    src/NewtonSolver.cpp
    src/SpectralTransformer.cpp
    src/StatePacker.cpp
    src/EOMevaluator.cpp
    src/Sampler.cpp
)

add_executable(a ${SOURCES})


# -------------------------
# Includes
# -------------------------

target_include_directories(a PRIVATE
    include
    include/Eigen 
    include/unsupported/Eigen
    ${FFTW_INCLUDE_DIRS}
)

# -------------------------
# Link libraries
# -------------------------

target_link_libraries(a PRIVATE
    ${FFTW_LIBRARIES}
    HDF5::HDF5
    nlohmann_json::nlohmann_json
    LAPACK::LAPACK
    lapacke
    ${BLAS_LIBRARIES}
)

# -------------------------
# Optimization flags
# -------------------------

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

target_compile_options(a PRIVATE
    -O3
    -march=native
    -Wall
)