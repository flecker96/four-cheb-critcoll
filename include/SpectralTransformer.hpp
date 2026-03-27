#pragma once
/**
 * @file SpectralTransformer.hpp
 * @brief Thin wrapper around FFTW to perform Fourier transforms and spectral operations.
 *
 * @details
 * Provides Fourier transforms (real/complex, forward/backward) and
 * additional utilities needed in the critical collapse solver:
 * - Differentiate in spectral space.
 * - Integrate with λ-shift.
 * - Interpolate Fourier series at arbitrary x.
 * - Mode truncation/doubling (halve/double resolution).
 * - Solve simple inhomogeneous equations in Fourier space.
 *
 * This class owns FFTW plans and work arrays. RAII ensures cleanup.
 */

#include "common.hpp"

/**
 * @class SpectralTransformer
 * @brief Encapsulates Fourier-based spectral operations.
 *
 * @section usage Usage
 * Construct with number of modes and period, then call forwardFFT/backwardFFT
 * to convert between time/τ domain and frequency domain. Use additional
 * helpers for differentiation, interpolation, etc.
 */
class SpectralTransformer
{
  private:
    size_t N, M;               ///< Number of real-space grid points in x and t directions.
    /*real_t period;*/          ///< Period (echoing period Δ).
    real_t k0;              ///< Fundamental frequency (2π/Δ).
    fftw_plan fft_forward_plan;
    fftw_plan fft_backward_plan;
    fftw_plan dct_plan;
    fftw_plan dct_plan_halved;
    fftw_complex *fft_data {nullptr}; ///< Work arrays.
    real_t *dct_in, *dct_out; ///< Work arrays.
    real_t *dct_in_halved, *dct_out_halved; ///< Work arrays.

  public:
    /**
     * @brief Construct transformer for N points with given period.
     * @param N      Number of real-space samples.
     * @param period Period of the signal (Δ).
     */
    explicit SpectralTransformer(size_t N, size_t M);

    /// Destructor: destroys FFTW plans and frees memory.
    ~SpectralTransformer();

    void forwardFFT(vec_complex& data);

    void backwardFFT(vec_complex& data);

    void forwardCheb(vec_complex& data);
    
    void backwardCheb(vec_complex& data);

    void forwardChebHalf(vec_complex& data);

    void backwardChebHalf(vec_complex& data);

    void differentiate_t(const vec_complex& in, vec_complex& out, real_t period_=0.0);

    void differentiate_x(const vec_complex& in, vec_complex& out);

    /**
     * @brief Double Fourier modes by zero-padding (double resolution).
     * @param in  Input coefficients.
     * @param out Output upsampled coefficients.
     */
    void doubleModes(const vec_complex& in, vec_complex& out);

    void halveModes(const vec_complex& in, vec_complex& out);
};
