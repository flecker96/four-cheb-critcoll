#pragma once
/**
 * @file common.hpp
 * @brief Common type aliases, utility functions, and third-party includes
 *        for the critical collapse solver.
 *
 * @details
 * This header centralizes:
 *  - Standard library and third-party includes.
 *  - Type aliases for reals, complex numbers, and vectors/matrices.
 *  - Parallelism headers (OpenMP/MPI) enabled via compile-time flags.
 *  - Shared numerical/IO utility functions (approximate equality checks,
 *    printing, writing vectors/matrices, small fitting routines).
 *
 * It is intended to be included across the project for consistent types
 * and helper functions.
 */

// ========== Standard Library ==========
#include <iostream>
#include <iomanip>
#include <cmath>
#include <complex>
#include <vector>
#include <array>
#include <algorithm>
#include <numeric>
#include <stdexcept>
#include <functional>
#include <memory>
#include <cstring>
#include <string>
#include <sstream>
#include <fstream>
#include <filesystem>
#include <cassert>
#include <thread>

// ========== Third-Party Libraries ==========
#include <nlohmann/json.hpp> ///< JSON for Modern C++

// ========== LAPACK ==========
#include <lapacke.h>        ///< LAPACK C interface

// ========== FFTW ==========
#include <fftw3.h>          ///< FFTW3 for spectral transforms

// =========== HDF5 =======================
#include "H5Cpp.h"

// ========== Parallelism ==========
#ifdef USE_OPENMP
#include <omp.h>            ///< OpenMP parallelism
#endif

#ifdef USE_MPI
#include <mpi.h>            ///< MPI parallelism
#endif

#ifdef USE_HYBRID
#include <mpi.h>
#include <omp.h>
#endif


// ========== Aliases ===============
using real_t     = double;                     ///< Floating point type used globally.
using complex_t  = std::complex<real_t>;       ///< Complex number type.
using vec_real   = std::vector<real_t>;        ///< Vector of real values.
using vec_complex= std::vector<complex_t>;     ///< Vector of complex values.
using mat_real   = std::vector<std::vector<real_t>>;   ///< Matrix of real values.
using mat_complex= std::vector<std::vector<complex_t>>;///< Matrix of complex values.
using json       = nlohmann::json;             ///< JSON type alias.

// ============ Common Functions =======

/**
 * @brief Check approximate equality of two complex numbers.
 * @param a First number.
 * @param b Second number.
 * @param tol Relative tolerance (default 1e-15).
 * @return true if |a-b| ≤ tol.
 */
bool almost_equal(complex_t a, complex_t b, double tol = 1e-15);

/**
 * @brief Check approximate equality of two real numbers.
 * @param a First number.
 * @param b Second number.
 * @param tol Relative tolerance (default 1e-15).
 * @return true if |a-b| ≤ tol.
 */
bool almost_equal(double a, double b, double tol = 1e-15);


bool isInteger(real_t x);

/**
 * @brief Print a real-valued vector to stdout.
 */
void print_vec(const vec_real& vec);

/**
 * @brief Print a vector of strings to stdout.
 */
void print_vec(const std::vector<std::string>& vec);

/**
 * @brief Print a complex-valued vector to stdout.
 */
void print_vec(const vec_complex& vec);

/**
 * @brief Write a matrix of reals to a text file.
 * @param filename Output file path.
 * @param mat      Matrix of real values.
 */
void write_mat(std::string filename, mat_real& mat);

/**
 * @brief Write a real vector to a text file.
 * @param filename Output file path.
 * @param vec      Vector of real values.
 */
void write_vec(std::string filename, vec_real& vec);

/**
 * @brief Write a complex vector to a text file.
 * @param filename Output file path.
 * @param vec      Vector of complex values.
 */
void write_vec(std::string filename, vec_complex& vec);

/**
 * @brief Construct a quadratic design matrix given three x-values.
 *
 * @details
 * Builds the Vandermonde-like design vector [1, x, x²] for each input
 * x ∈ {x1, x2, x3}, concatenated into a flat vector.
 *
 * @param x1 First x-value.
 * @param x2 Second x-value.
 * @param x3 Third x-value.
 * @return Flattened design matrix coefficients.
 */
vec_real build_design_matrix(real_t x1, real_t x2, real_t x3);

/**
 * @brief Fit a quadratic polynomial y ≈ a + bx + cx² in least squares sense.
 *
 * @param x_vals Vector of x-samples.
 * @param y_vals Vector of y-samples (same length as x_vals).
 * @return Coefficients [a, b, c] minimizing least squares error.
 *
 * @throws std::invalid_argument if input sizes mismatch or insufficient points.
 */
vec_real fit_quadratic_least_squares(const vec_real& x_vals, const vec_real& y_vals);

/**
 * @brief Extracts elements with stride n, ensuring the last element is always included.
 *
 * @param a Input vector to sample from.
 * @param n The stride interval (every n-th element).
 * @return A new vector containing the sampled elements, plus the final element of 'a'.
 */
vec_real every_nth_element(const vec_real& a, size_t n);
