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
    fftw_plan forward_plan {nullptr};   ///< FFTW plan: forward.
    fftw_plan backward_plan {nullptr};  ///< FFTW plan: backward.
    fftw_plan plan_dct {nullptr}; 
    fftw_plan forward_plan_many;
    fftw_plan forward_plan_many2;
    fftw_plan backward_plan_many;
    fftw_plan backward_plan_many2;
    fftw_plan plan_dct_many;
    fftw_plan plan_dct_many_halved;
    fftw_complex *forward_data {nullptr}, *backward_data {nullptr}, *fft_data {nullptr}; ///< Work arrays.
    real_t *in_data_dct {nullptr}, *out_data_dct {nullptr}; ///< Work arrays.
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

    void forwardFFT_time(const vec_complex& in, vec_complex& out);

    void forwardFFT_time2(vec_complex& data);

    void backwardFFT_time(const vec_complex& in, vec_complex& out);

    void backwardFFT_time2(vec_complex& data);

    void forwardDCT_space(const vec_complex& in, vec_complex& out);

    void forwardDCT_space2(vec_complex& data);

    void backwardDCT_space(const vec_complex& in, vec_complex& out);
    
    void backwardDCT_space2(vec_complex& data);

    void forwardChebHalf(vec_complex& data);

    void backwardChebHalf(vec_complex& data);

    /// Forward FFT (real → complex).
    void forwardFFT(const vec_real& in, vec_complex& out);

    /// Backward FFT (complex → real).
    void backwardFFT(const vec_complex& in, vec_real& out);

    /// Forward FFT (complex → complex).
    void forwardFFT(const vec_complex& in, vec_complex& out);

    /// Backward FFT (complex → complex).
    void backwardFFT(const vec_complex& in, vec_complex& out);

    void forwardDCT(const vec_real& in, vec_real& out);

    void backwardDCT(const vec_real& in, vec_real& out);

    /**
     * @brief Differentiate a Fourier series.
     * @param in      Input Fourier coefficients.
     * @param out     Output differentiated Fourier coefficients.
     * @param period_ Optional override of period (defaults to ctor value).
     */
    void differentiate(const vec_complex& in, vec_complex& out, real_t period_=0.0);

    void differentiate_t(const vec_complex& in, vec_complex& out, real_t period_=0.0);

    void differentiate_x(const vec_complex& in, vec_complex& out);

    /**
     * @brief Truncate Fourier modes by factor 2 (half resolution).
     * @param in  Input coefficients.
     * @param out Output truncated coefficients.
     */
    void halveModesOld(const vec_complex& in, vec_complex& out);

    void doubleModesOld(const vec_complex& in, vec_complex& out);
    /**
     * @brief Double Fourier modes by zero-padding (double resolution).
     * @param in  Input coefficients.
     * @param out Output upsampled coefficients.
     */
    void doubleModes(const vec_complex& in, vec_complex& out);

    void halveModes(const vec_complex& in, vec_complex& out);
};
