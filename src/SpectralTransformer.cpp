//==============================================================================
// SpectralTransformer.cpp
// Thin wrapper around FFTW for periodic spectral transforms and utilities.
// Features:
//   • Forward/backward FFTs for complex data (FFTW complex-to-complex).
//   • Spectral differentiation in t with 2π/period wave numbers (odd Nyquist-safe).
//   • Spectral differentiation in x with Clenshaw recurrence pattern.
//   • Mode decimation/expansion (halveModes/doubleModes) for Newton packing.
//   • Inhomogeneous first-order ODE solve via integrating factor in τ.
// Notes:
//   • FFTW uses e^{-i k x} for FORWARD; we invert directions to match DFT
//     conventions used elsewhere, hence the BACKWARD/ FORWARD swap in plans.
//   • Keep the user's one-line loop/if braces `{}` intact.
//==============================================================================

#include "SpectralTransformer.hpp"

//------------------------------------------------------------------------------
// Ctor: allocate FFTW work buffers and build forward/backward plans.
// We deliberately swap BACKWARD/FORWARD flags to match our DFT sign convention.
//------------------------------------------------------------------------------
SpectralTransformer::SpectralTransformer(size_t N_, size_t M_)
    : N(N_), M(M_)
{
    
    fft_data = fftw_alloc_complex(N * M);

    dct_in  = fftw_alloc_real(N * M);
    dct_out = fftw_alloc_real(N * M);

    dct_in_halved  = fftw_alloc_real(N * M / 4);
    dct_out_halved = fftw_alloc_real(N * M / 4);

    // Plan for strided FFT on whole array
    const int rank = 1;
    int n[1] = { static_cast<int>(N) };

    const int howmany = static_cast<int>(M);

    const int istride = static_cast<int>(M);
    const int ostride = static_cast<int>(M);

    const int idist = 1;
    const int odist = 1;

    fft_forward_plan = fftw_plan_many_dft(
        rank, n, howmany,
        fft_data, nullptr,
        istride, idist,
        fft_data, nullptr,
        ostride, odist,
        FFTW_BACKWARD,
        FFTW_ESTIMATE
    );

    fft_backward_plan = fftw_plan_many_dft(
        rank, n, howmany,
        fft_data, nullptr,
        istride, idist,
        fft_data, nullptr,
        ostride, odist,
        FFTW_FORWARD,
        FFTW_ESTIMATE
    );


    // Plan for strided DCT on whole array
    int n_dct[1] = { static_cast<int>(M) };   // DCT length Nx
    const int howmany_dct = static_cast<int>(N);  // number of time slices

    const int stride_dct = 1;     // contiguous in x
    const int dist_dct   = static_cast<int>(M);

    fftw_r2r_kind kind[1] = { FFTW_REDFT00 };

    dct_plan = fftw_plan_many_r2r(
        rank, n_dct, howmany_dct,
        dct_in, nullptr,
        stride_dct, dist_dct,
        dct_out, nullptr,
        stride_dct, dist_dct,
        kind,
        FFTW_ESTIMATE
    );

    // plan for strided DCT of the dealiased array
    int n_dct_halved[1] = { static_cast<int>(M/2) };   // DCT length Nx
    const int howmany_dct_halved = static_cast<int>(N/2);  // number of time slices

    const int stride_dct_halved = 1;     // contiguous in x
    const int dist_dct_halved   = static_cast<int>(M/2);

    fftw_r2r_kind kind_halved[1] = { FFTW_REDFT00 };

   dct_plan_halved = fftw_plan_many_r2r(
        rank, n_dct_halved, howmany_dct_halved,
        dct_in_halved, nullptr,
        stride_dct_halved, dist_dct_halved,
        dct_out_halved, nullptr,
        stride_dct_halved, dist_dct_halved,
        kind_halved,
        FFTW_ESTIMATE
    );

}

//------------------------------------------------------------------------------
// Dtor: free plans and work arrays.
//------------------------------------------------------------------------------
SpectralTransformer::~SpectralTransformer()
{
    fftw_destroy_plan(dct_plan);
    fftw_destroy_plan(dct_plan_halved);
    fftw_destroy_plan(fft_forward_plan);
    fftw_destroy_plan(fft_backward_plan);
    fftw_free(fft_data);
    fftw_free(dct_in);
    fftw_free(dct_out);
    fftw_free(dct_in_halved);
    fftw_free(dct_out_halved);
}


void SpectralTransformer::forwardFFT(vec_complex& data)
{
    size_t size = N * M;

    // (optional safety check)
    if (data.size() != size)
        throw std::runtime_error("forwardFFT_time: wrong input size");

    // reinterpret std::complex as fftw_complex
    auto* ptr = reinterpret_cast<fftw_complex*>(data.data());

    // execute in-place FFT
    fftw_execute_dft(fft_forward_plan, ptr, ptr);

    // normalization
    for (size_t i = 0; i < size; ++i)
        data[i] /= static_cast<real_t>(N);
}

void SpectralTransformer::backwardFFT(vec_complex& data)
{
    size_t size = N * M;

    // (optional safety check)
    if (data.size() != size)
        throw std::runtime_error("forwardFFT_time: wrong input size");

    // reinterpret std::complex as fftw_complex
    auto* ptr = reinterpret_cast<fftw_complex*>(data.data());

    // execute in-place FFT
    fftw_execute_dft(fft_backward_plan, ptr, ptr);

}


void SpectralTransformer::forwardCheb(vec_complex& data)
{
    if (data.size() != N*M)
        throw std::runtime_error("forwardDCT_space: wrong size");

    size_t size = N * M;

    vec_real real_part(size);
    vec_real imag_part(size);

    for(size_t i=0;i<size;i++){
        real_part[i] = data[i].real();
        imag_part[i] = data[i].imag();
    }

    /* real part */

    std::copy(real_part.begin(), real_part.end(), dct_in);
    fftw_execute(dct_plan);
    std::copy(dct_out, dct_out + size, real_part.begin());

    /* imaginary part */

    std::copy(imag_part.begin(), imag_part.end(), dct_in);
    fftw_execute(dct_plan);
    std::copy(dct_out, dct_out + size, imag_part.begin());

    for(size_t i=0;i<size;i++){
        data[i] = complex_t(
            real_part[i] / (static_cast<real_t>(M) - 1.0),
            imag_part[i] / (static_cast<real_t>(M) - 1.0)
        );
    }
}

void SpectralTransformer::backwardCheb(vec_complex& data)
{
    if (data.size() != N*M)
        throw std::runtime_error("forwardDCT_space: wrong size");

    size_t size = N * M;

    vec_real real_part(size);
    vec_real imag_part(size);

    for(size_t i=0;i<size;i++){
        real_part[i] = data[i].real();
        imag_part[i] = data[i].imag();
    }

    /* real part */

    std::copy(real_part.begin(), real_part.end(), dct_in);
    fftw_execute(dct_plan);
    std::copy(dct_out, dct_out + size, real_part.begin());

    /* imaginary part */

    std::copy(imag_part.begin(), imag_part.end(), dct_in);
    fftw_execute(dct_plan);
    std::copy(dct_out, dct_out + size, imag_part.begin());

    for(size_t i=0;i<size;i++){
        data[i] = complex_t(
            real_part[i] / 2.0 ,
            imag_part[i] / 2.0
        );
    }
}

void SpectralTransformer::forwardChebHalf(vec_complex& data)
{
    if (data.size() != N*M/4)
        throw std::runtime_error("forwardDCT_space: wrong size");

    size_t size = N * M / 4;

    vec_real real_part(size);
    vec_real imag_part(size);

    for(size_t i=0;i<size;i++){
        real_part[i] = data[i].real();
        imag_part[i] = data[i].imag();
    }

    /* real part */

    std::copy(real_part.begin(), real_part.end(), dct_in_halved);
    fftw_execute(dct_plan_halved);
    std::copy(dct_out_halved, dct_out_halved + size, real_part.begin());

    /* imaginary part */

    std::copy(imag_part.begin(), imag_part.end(), dct_in_halved);
    fftw_execute(dct_plan_halved);
    std::copy(dct_out_halved, dct_out_halved + size, imag_part.begin());

    for(size_t i=0;i<size;i++){
        data[i] = complex_t(
            real_part[i] / (static_cast<real_t>(M) / 2.0 - 1.0),
            imag_part[i] / (static_cast<real_t>(M) / 2.0 - 1.0)
        );
    }
}

void SpectralTransformer::backwardChebHalf(vec_complex& data)
{
    if (data.size() != N*M/4)
        throw std::runtime_error("forwardDCT_space: wrong size");

    size_t size = N * M / 4;

    vec_real real_part(size);
    vec_real imag_part(size);

    for(size_t i=0;i<size;i++){
        real_part[i] = data[i].real();
        imag_part[i] = data[i].imag();
    }

    /* real part */

    std::copy(real_part.begin(), real_part.end(), dct_in_halved);
    fftw_execute(dct_plan_halved);
    std::copy(dct_out_halved, dct_out_halved + size, real_part.begin());

    /* imaginary part */

    std::copy(imag_part.begin(), imag_part.end(), dct_in_halved);
    fftw_execute(dct_plan_halved);
    std::copy(dct_out_halved, dct_out_halved + size, imag_part.begin());

    for(size_t i=0;i<size;i++){
        data[i] = complex_t(
            real_part[i] / 2.0,
            imag_part[i] / 2.0
        );
    }
}


void SpectralTransformer::differentiate_t(
    const vec_complex& in,
    vec_complex& out,
    real_t period)
{
    real_t k0c = (period == 0.0) ? k0 : 2*M_PI / period;

    if (out.size() != N*M)
        throw std::runtime_error("differentiate_t: wrong output size");

    for (size_t i = 0; i < M; ++i)  // loop over spatial index
    {
        for (size_t k = 0; k < N; ++k)  // Fourier mode index
        {
            size_t idx = k * M + i;

            if (k != N/2)
            {
                int m = (k < N/2)
                    ? static_cast<int>(k)
                    : static_cast<int>(k) - static_cast<int>(N);

                out[idx] = -complex_t(0.0, m * k0c) * in[idx];
            }
            else
            {
                out[idx] = complex_t(0.0);
            }
        }
    }
}

void SpectralTransformer::differentiate_x(
    const vec_complex& in,
    vec_complex& out)
{
    if (out.size() != N*M)
        throw std::runtime_error("differentiate_x: wrong output size");

    for (size_t i = 0; i < N; ++i)  // loop over temporal index
    {
        out[i*M + (M-1)] = complex_t(0.0, 0.0);
        out[i*M + (M-2)] = 2.0 * static_cast<int>(M - 1) * in[i*M + (M-1)];

        for (int k = static_cast<int>(M - 3); k >= 0; --k)  // Chebyshev mode index
        {
            size_t idx = i * M + k;

            out[idx] = out[idx + 2] + 2.0 * (k + 1) * in[idx + 1];
        }
    }
}

void SpectralTransformer::doubleModes(const vec_complex& in, vec_complex& out)
{ 
    std::fill(out.begin(), out.end(), complex_t(0.0, 0.0));

    //Now fill the vector with the IR modes
    for (size_t i=0; i<M/2; ++i)
    {
        //Positive frequencies
        for (size_t k=0; k<N/4; ++k)
        {
            out[M*k + i] = in[(M/2)*k + i]; 
        }
        
        //Nyquist
        out[M*(N/4) + i] = 0.5*in[(M/2)*N/4 + i];
        out[M*(3*N/4) + i] = 0.5*in[(M/2)*N/4 + i]; 

        //Negative frequencies
        for (size_t k=N/4+1; k<N/2; ++k)
        {
            out[M*(N/2 + k) + i] = in[(M/2)*k + i]; 
        }

        //Multiply old highest Chebyshev coefficient by 2 (normalization of series)
        if (i==(M/2-1)){
            for (size_t k=0; k<N; ++k)
                {
                    out[M*k + i] = 2.0*out[M*k + i]; 
                }
        }
    }
}

void SpectralTransformer::doubleModes_t(const vec_complex& in, vec_complex& out)
{ 
    size_t Ntnew = 2*N;
    size_t Nxnew = M; 

    if (out.size() != 2*N*M)
        throw std::runtime_error("doubleModes_t: wrong output size");

    std::fill(out.begin(), out.end(), complex_t(0.0, 0.0));

    //Now fill the vector with the IR modes
    for (size_t i=0; i<Nxnew; ++i)
    {
        //Positive frequencies
        for (size_t k=0; k<Ntnew/4; ++k)
        {
            out[Nxnew*k + i] = in[Nxnew*k + i]; 
        }
        
        //Nyquist
        out[Nxnew*(Ntnew/4) + i] = 0.5*in[Nxnew*Ntnew/4 + i];
        out[Nxnew*(3*Ntnew/4) + i] = 0.5*in[Nxnew*Ntnew/4 + i]; 

        //Negative frequencies
        for (size_t k=Ntnew/4+1; k<Ntnew/2; ++k)
        {
            out[Nxnew*(Ntnew/2 + k) + i] = in[Nxnew*k + i]; 
        }
    }
}

void SpectralTransformer::halveModes(const vec_complex& in, vec_complex& out)
{
    vec_complex tmp(N*M/4);

    //Read out IR modes
    for (size_t i=0; i<M/2; ++i)
    {
        //Positive frequencies
        for (size_t k=0; k<N/4; ++k)
        {
            tmp[(M/2)*k + i] = in[M*k + i]; 
        }
        
        //Nyquist
        tmp[(M/2)*N/4 + i] = in[M*(N/4) + i] + in[M*(3*N/4) + i];

        //Negative frequencies
        for (size_t k=N/4+1; k<N/2; ++k)
        {
            tmp[(M/2)*k + i] = in[M*(N/2 + k) + i];
        }

        //Multiply highest Chebyshev coefficient by 2 (normalization of series)
        if (i==(M/2-1)){
            for (size_t k=0; k<N/2; ++k)
                {
                    tmp[(M/2)*k + i] = tmp[(M/2)*k + i] / 2.0; 
                }
        }
    }

    out.swap(tmp);
}
