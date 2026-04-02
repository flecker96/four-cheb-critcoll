#pragma once
/**
 * @file SpectralTransformer.hpp
 * @brief Thin wrapper around FFTW to perform Fourier transforms and spectral operations.
 *
 * @details
 * This class owns FFTW plans and work arrays. RAII ensures cleanup.
 */

#include "common.hpp"

/**
 * @class SpectralTransformer
 * @brief Encapsulates Fourier-based spectral operations.
 * @section usage Usage
 * Construct with number of modes, then call forwardFFT/backwardFFT and Cheb
 * to convert between real domain and spectral domain. Use additional
 * helpers for differentiation, interpolation, etc.
 */
class SpectralTransformer
{
  private:
    size_t N, M;               ///< Number of real-space grid points in x and t directions.
    real_t k0;              ///< Fundamental frequency (2π/Δ).
    fftw_plan fft_forward_plan;
    fftw_plan fft_backward_plan;
    fftw_plan dct_plan;
    fftw_plan dct_plan_halved;
    fftw_complex *fft_data {nullptr}; ///< Work arrays.
    real_t *dct_in {nullptr}, *dct_out {nullptr}; ///< Work arrays.
    real_t *dct_in_halved {nullptr}, *dct_out_halved {nullptr}; ///< Work arrays.

  public:
    /**
     * @brief Construct transformer for NxM points.
     * @param N      Number of real-space samples in t direction.
     * @param M      Number of real-space samples in x direction.
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

    /**
     * @brief Compute spectral t-derivative.
     * @param in  Input coefficients.
     * @param out Output coefficients.
     * @param period Delta.
     */
    void differentiate_t(const vec_complex& in, vec_complex& out, real_t period_=0.0);

    /**
     * @brief Compute spectral x-derivative.
     * @param in  Input coefficients.
     * @param out Output coefficients.
     */
    void differentiate_x(const vec_complex& in, vec_complex& out);

    /**
     * @brief Double Fourier and Chebyshev modes by zero-padding (double resolution).
     * @param in  Input coefficients.
     * @param out Output upsampled coefficients.
     */
    void doubleModes(const vec_complex& in, vec_complex& out);

    /**
     * @brief Halve Fourier and Chebyshev modes by deleting UV part.
     * @param in  Input coefficients.
     * @param out Output upsampled coefficients.
     */
    void halveModes(const vec_complex& in, vec_complex& out);

    /**
     * @brief Double only Fourier modes by zero-padding (double resolution).
     * @param in  Input coefficients.
     * @param out Output upsampled coefficients.
     */

    //Increase by certain factor
    void increaseModes_spec(const vec_complex& in, vec_complex& out, real_t fact, real_t facx);

    //Decrease by certain factor
    void decreaseModes_spec(const vec_complex& in, vec_complex& out, real_t fact, real_t facx);
};
