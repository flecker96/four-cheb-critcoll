//==============================================================================
// SpectralTransformer.cpp
// Thin wrapper around FFTW for periodic spectral transforms and utilities.
// Features:
//   • Forward/backward FFTs for real/complex data (FFTW complex-to-complex).
//   • Spectral differentiation with 2π/period wave numbers (odd Nyquist-safe).
//   • λ-integration in Fourier space: solve (λ - i m k0) Ĝ = ˆf mode-wise.
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

    forward_data = fftw_alloc_complex(N * M);
    backward_data = fftw_alloc_complex(N * M);
    in_data_dct = fftw_alloc_real(M);
    out_data_dct = fftw_alloc_real(M);

    // Changing FFT forward and backward flag due to DFT definition of FFTW3 (-exp(...) vs. exp(...))
    forward_plan  = fftw_plan_dft_1d(static_cast<int>(N), forward_data, backward_data, FFTW_BACKWARD, FFTW_ESTIMATE);
    backward_plan = fftw_plan_dft_1d(static_cast<int>(N), backward_data, forward_data, FFTW_FORWARD,  FFTW_ESTIMATE);

    // Plan for DCT is the same for forward and backward trafo
    plan_dct  = fftw_plan_r2r_1d(static_cast<int>(M), in_data_dct, out_data_dct, FFTW_REDFT00, FFTW_ESTIMATE);

    // Plan for strided FFT on whole array
    const int rank = 1;
    int n[1] = { static_cast<int>(N) };

    const int howmany = static_cast<int>(M);

    const int istride = static_cast<int>(M);
    const int ostride = static_cast<int>(M);

    const int idist = 1;
    const int odist = 1;

    forward_plan_many2 = fftw_plan_many_dft(
        rank, n, howmany,
        fft_data, nullptr,
        istride, idist,
        fft_data, nullptr,
        ostride, odist,
        FFTW_BACKWARD,
        FFTW_ESTIMATE
    );

    backward_plan_many2 = fftw_plan_many_dft(
        rank, n, howmany,
        fft_data, nullptr,
        istride, idist,
        fft_data, nullptr,
        ostride, odist,
        FFTW_FORWARD,
        FFTW_ESTIMATE
    );

    //Old versions
    forward_plan_many = fftw_plan_many_dft(
        rank, n, howmany,
        forward_data, nullptr,
        istride, idist,
        backward_data, nullptr,
        ostride, odist,
        FFTW_BACKWARD,
        FFTW_ESTIMATE
    );

    backward_plan_many = fftw_plan_many_dft(
        rank, n, howmany,
        backward_data, nullptr,
        istride, idist,
        forward_data, nullptr,
        ostride, odist,
        FFTW_FORWARD,
        FFTW_ESTIMATE
    );

    //Plan for DCT of big array (padded by zeros for AA)
    int n_dct[1] = { static_cast<int>(M) };   // DCT length Nx
    const int howmany_dct = static_cast<int>(N);  // number of time slices

    const int stride_dct = 1;     // contiguous in x
    const int dist_dct   = static_cast<int>(M);

    fftw_r2r_kind kind[1] = { FFTW_REDFT00 };

    plan_dct_many = fftw_plan_many_r2r(
        rank, n_dct, howmany_dct,
        dct_in, nullptr,
        stride_dct, dist_dct,
        dct_out, nullptr,
        stride_dct, dist_dct,
        kind,
        FFTW_ESTIMATE
    );

    // plan for DCT of the dealiased array
    int n_dct_halved[1] = { static_cast<int>(M/2) };   // DCT length Nx
    const int howmany_dct_halved = static_cast<int>(N/2);  // number of time slices

    const int stride_dct_halved = 1;     // contiguous in x
    const int dist_dct_halved   = static_cast<int>(M/2);

    fftw_r2r_kind kind_halved[1] = { FFTW_REDFT00 };

    plan_dct_many_halved = fftw_plan_many_r2r(
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
    fftw_destroy_plan(forward_plan_many);
    fftw_destroy_plan(backward_plan_many);
    fftw_destroy_plan(plan_dct_many);
    fftw_destroy_plan(plan_dct_many_halved);
    fftw_destroy_plan(forward_plan);
    fftw_destroy_plan(backward_plan);
    fftw_destroy_plan(plan_dct);
    fftw_free(forward_data);
    fftw_free(backward_data);
    fftw_free(dct_in);
    fftw_free(dct_out);
    fftw_free(dct_in_halved);
    fftw_free(dct_out_halved);
    fftw_free(in_data_dct);
    fftw_free(out_data_dct);
}

void SpectralTransformer::forwardFFT_time(const vec_complex& in,
                                          vec_complex& out)
{
    size_t size = N * M;

    for (size_t i = 0; i < size; ++i)
    {
        forward_data[i][0] = in[i].real();
        forward_data[i][1] = in[i].imag();
    }

    fftw_execute(forward_plan_many);

    if (out.size() != size) out.resize(size);

    for (size_t i = 0; i < size; ++i)
    {
        out[i] = complex_t(backward_data[i][0],
                           backward_data[i][1])
                 / static_cast<real_t>(N);
    }
}

void SpectralTransformer::forwardFFT_time2(vec_complex& data)
{
    size_t size = N * M;

    // (optional safety check)
    if (data.size() != size)
        throw std::runtime_error("forwardFFT_time: wrong input size");

    // reinterpret std::complex as fftw_complex
    auto* ptr = reinterpret_cast<fftw_complex*>(data.data());

    // execute in-place FFT
    fftw_execute_dft(forward_plan_many2, ptr, ptr);

    // normalization
    for (size_t i = 0; i < size; ++i)
        data[i] /= static_cast<real_t>(N);
}

void SpectralTransformer::backwardFFT_time(const vec_complex& in,
                                          vec_complex& out)
{
    size_t size = N * M;

    for (size_t i = 0; i < size; ++i)
    {
        backward_data[i][0] = in[i].real();
        backward_data[i][1] = in[i].imag();
    }

    fftw_execute(backward_plan_many);

    out.resize(size);

    for (size_t i = 0; i < size; ++i)
    {
        out[i] = complex_t(forward_data[i][0],
                           forward_data[i][1]);
    }
}

void SpectralTransformer::backwardFFT_time2(vec_complex& data)
{
    size_t size = N * M;

    // (optional safety check)
    if (data.size() != size)
        throw std::runtime_error("forwardFFT_time: wrong input size");

    // reinterpret std::complex as fftw_complex
    auto* ptr = reinterpret_cast<fftw_complex*>(data.data());

    // execute in-place FFT
    fftw_execute_dft(backward_plan_many2, ptr, ptr);

}


void SpectralTransformer::forwardDCT_space(
        const vec_complex& in,
        vec_complex& out)
{
    if (in.size() != N*M)
        throw std::runtime_error("forwardDCT_space: wrong size");

    size_t size = N * M;

    vec_real real_part(size);
    vec_real imag_part(size);

    for(size_t i=0;i<size;i++){
        real_part[i] = in[i].real();
        imag_part[i] = in[i].imag();
    }

    /* real part */

    std::copy(real_part.begin(), real_part.end(), dct_in);
    fftw_execute(plan_dct_many);
    std::copy(dct_out, dct_out + size, real_part.begin());

    /* imaginary part */

    std::copy(imag_part.begin(), imag_part.end(), dct_in);
    fftw_execute(plan_dct_many);
    std::copy(dct_out, dct_out + size, imag_part.begin());

    out.resize(size);

    for(size_t i=0;i<size;i++){
        out[i] = complex_t(
            real_part[i] / (static_cast<real_t>(M) - 1.0),
            imag_part[i] / (static_cast<real_t>(M) - 1.0)
        );
    }
}

void SpectralTransformer::forwardDCT_space2(vec_complex& data)
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
    fftw_execute(plan_dct_many);
    std::copy(dct_out, dct_out + size, real_part.begin());

    /* imaginary part */

    std::copy(imag_part.begin(), imag_part.end(), dct_in);
    fftw_execute(plan_dct_many);
    std::copy(dct_out, dct_out + size, imag_part.begin());

    for(size_t i=0;i<size;i++){
        data[i] = complex_t(
            real_part[i] / (static_cast<real_t>(M) - 1.0),
            imag_part[i] / (static_cast<real_t>(M) - 1.0)
        );
    }
}

void SpectralTransformer::backwardDCT_space(
        const vec_complex& in,
        vec_complex& out)
{
    if (in.size() != N*M)
        throw std::runtime_error("forwardDCT_space: wrong size");

    size_t size = N * M;

    vec_real real_part(size);
    vec_real imag_part(size);

    for(size_t i=0;i<size;i++){
        real_part[i] = in[i].real();
        imag_part[i] = in[i].imag();
    }

    /* real part */

    std::copy(real_part.begin(), real_part.end(), dct_in);
    fftw_execute(plan_dct_many);
    std::copy(dct_out, dct_out + size, real_part.begin());

    /* imaginary part */

    std::copy(imag_part.begin(), imag_part.end(), dct_in);
    fftw_execute(plan_dct_many);
    std::copy(dct_out, dct_out + size, imag_part.begin());

    out.resize(size);

    for(size_t i=0;i<size;i++){
        out[i] = complex_t(
            real_part[i] / 2.0 ,
            imag_part[i] / 2.0
        );
    }
}

void SpectralTransformer::backwardDCT_space2(vec_complex& data)
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
    fftw_execute(plan_dct_many);
    std::copy(dct_out, dct_out + size, real_part.begin());

    /* imaginary part */

    std::copy(imag_part.begin(), imag_part.end(), dct_in);
    fftw_execute(plan_dct_many);
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
    fftw_execute(plan_dct_many_halved);
    std::copy(dct_out_halved, dct_out_halved + size, real_part.begin());

    /* imaginary part */

    std::copy(imag_part.begin(), imag_part.end(), dct_in_halved);
    fftw_execute(plan_dct_many_halved);
    std::copy(dct_out_halved, dct_out_halved + size, imag_part.begin());

    for(size_t i=0;i<size;i++){
        data[i] = complex_t(
            real_part[i] / (static_cast<real_t>(M/2) - 1.0),
            imag_part[i] / (static_cast<real_t>(M/2) - 1.0)
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
    fftw_execute(plan_dct_many_halved);
    std::copy(dct_out_halved, dct_out_halved + size, real_part.begin());

    /* imaginary part */

    std::copy(imag_part.begin(), imag_part.end(), dct_in_halved);
    fftw_execute(plan_dct_many_halved);
    std::copy(dct_out_halved, dct_out_halved + size, imag_part.begin());

    for(size_t i=0;i<size;i++){
        data[i] = complex_t(
            real_part[i] / 2.0,
            imag_part[i] / 2.0
        );
    }
}

//------------------------------------------------------------------------------
// forwardFFT (real → complex): pack real input into complex buffer, execute,
// and scale by 1/N to obtain unitary-like convention.
//------------------------------------------------------------------------------
void SpectralTransformer::forwardFFT(const vec_real& in, vec_complex& out)
{
    for (size_t i=0; i<N; ++i)
    {
        forward_data[i][0] = in[i];
        forward_data[i][1] = 0.0;
    }

    fftw_execute(forward_plan);

    out.resize(N);
    for (size_t i=0; i<N; ++i)
    {
        out[i] = complex_t(backward_data[i][0], backward_data[i][1]) / static_cast<real_t>(N);
    }
}

//------------------------------------------------------------------------------
// backwardFFT (complex → real): inverse transform and extract real part.
//------------------------------------------------------------------------------
void SpectralTransformer::backwardFFT(const vec_complex& in, vec_real& out)
{
    for (size_t i=0; i<N; ++i)
    {
        backward_data[i][0] = in[i].real();
        backward_data[i][1] = in[i].imag();
    }

    fftw_execute(backward_plan);

    out.resize(N);
    for (size_t i=0; i<N; ++i)
    {
        out[i] = forward_data[i][0];
    }
}

//------------------------------------------------------------------------------
// forwardFFT (complex → complex): complex-to-complex forward with 1/N scaling.
//------------------------------------------------------------------------------
void SpectralTransformer::forwardFFT(const vec_complex& in, vec_complex& out)
{
    for (size_t i=0; i<N; ++i)
    {
        forward_data[i][0] = in[i].real();
        forward_data[i][1] = in[i].imag();
    }

    fftw_execute(forward_plan);

    out.resize(N);
    for (size_t i=0; i<N; ++i)
    {
        out[i] = complex_t(backward_data[i][0], backward_data[i][1]) / static_cast<real_t>(N);
    }
}

//------------------------------------------------------------------------------
// backwardFFT (complex → complex): inverse complex transform.
//------------------------------------------------------------------------------
void SpectralTransformer::backwardFFT(const vec_complex& in, vec_complex& out)
{
    for (size_t i=0; i<N; ++i)
    {
        backward_data[i][0] = in[i].real();
        backward_data[i][1] = in[i].imag();
    }

    fftw_execute(backward_plan);

    out.resize(N);
    for (size_t i=0; i<N; ++i)
    {
        out[i] = complex_t(forward_data[i][0], forward_data[i][1]);
    }
}

void SpectralTransformer::forwardDCT(const vec_real& in, vec_real& out)
{
    for (size_t i=0; i<N; ++i)
    {
        in_data_dct[i] = in[i];
    }

    fftw_execute(plan_dct);

    out.resize(N);
    for (size_t i=0; i<N; ++i)
    {
        out[i] = out_data_dct[i] / (static_cast<real_t>(N) - 1.0);
    }
}

void SpectralTransformer::backwardDCT(const vec_real& in, vec_real& out)
{
    for (size_t i=0; i<N; ++i)
    {
        in_data_dct[i] = in[i];
    }

    fftw_execute(plan_dct);

    out.resize(N);
    for (size_t i=0; i<N; ++i)
    {
        out[i] = out_data_dct[i] / 2.0;
    }
}

//------------------------------------------------------------------------------
// differentiate: spectral derivative in τ.
// Uses m = k for k < N/2 and m = k - N for k ≥ N/2 (wrap-around), and
// multiplies by -i m k0c. The Nyquist (k=N/2) is set to 0 for stability.
// Optionally override the period (changes k0c).
//------------------------------------------------------------------------------
void SpectralTransformer::differentiate(const vec_complex& in, vec_complex& out, real_t period)
{
    real_t k0c = (period == 0.0) ? k0 : 2*M_PI / period;
    size_t size = in.size();
        
    out.resize(size);
    for (size_t k=0; k<size; ++k)
    {
        if (k != size/2)
        {
            int m = (k<size/2) ? static_cast<int>(k) : static_cast<int>(k) - static_cast<int>(size);
            out[k] = -complex_t(0.0, m*k0c) * in[k];
        }
        else
        {
            out[size/2] = complex_t(0.0);
        }   
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
        out[i*M + (M-2)] = 2.0 * static_cast<int>(M - 1) * in[i*M + (M-2)];

        for (int k = static_cast<int>(M - 3); k >= 0; --k)  // Chebyshev mode index
        {
            size_t idx = i * M + k;

            out[idx] = out[idx + 2] + 2.0 * (k + 1) * in[idx + 1];
        }
    }
}

//------------------------------------------------------------------------------
// halveModes: downsample complex spectrum by 2, folding high/low halves and
// summing the two cosine Nyquist partners into the new center bin.
// Layout: [0 .. N/2-1 | N/2 | N/2+1 .. N-1]  → size→N/2 array.
//------------------------------------------------------------------------------
void SpectralTransformer::halveModesOld(const vec_complex& in, vec_complex& out)
{
    size_t N_in  = in.size();
    size_t N_out = N_in/2;

    vec_complex tmp(N_out);

    for (size_t k=0; k<N_out/2; ++k)
    {
        tmp[k] = in[k];
    }

    // Fold the two Nyquist partners
    tmp[N_out/2] = in[N_out/2] + in[3*N_out/2];

    for (size_t k=N_out/2+1; k<N_out; ++k)
    {
        tmp[k] = in[N_out+k];
    }

    out.swap(tmp);
}

/*void SpectralTransformer::halve_block(
    const vec_complex& in,
    vec_complex& out)
{
    size_t N_in  = in.size();
    size_t N_out = N_in/4;

    vec_complex tmp(N_out);
    
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
}*/

//------------------------------------------------------------------------------
// doubleModes: upsample complex spectrum by 2, splitting the cosine Nyquist
// into two half-amplitude copies and zeroing interstitial high-freq slots.
// Result size is 2*N_in with appropriate zero padding.
//------------------------------------------------------------------------------
void SpectralTransformer::doubleModesOld(const vec_complex& in, vec_complex& out)
{
    size_t N_in  = in.size();
    size_t N_out = 2*N_in;

    out.resize(N_out);

    // Low frequencies (strictly below Nyquist)
    for (size_t k=0; k<N_in/2; ++k)
    {
        out[k] = in[k];
    }

    // Split the Nyquist cosine into two symmetric bins
    out[N_in/2]   = 0.5*in[N_in/2];
    out[3*N_in/2] = 0.5*in[N_in/2];

    // Upper half copied to the high-frequency end
    for (size_t k=N_in/2; k<N_in; ++k)
    {
        out[N_in+k] = in[k];
    }

    // Zero the guard band between the split Nyquist bins
    for (size_t k=N_in/2+1; k<3*N_in/2; ++k)
    {
        out[k] = complex_t(0.0);
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

        //Divide old highest Chebyshev coefficient by 2 (normalization of series)
        if (i==(M/2-1)){
            for (size_t k=0; k<N; ++k)
                {
                    out[M*k + i] = 0.5*out[M*k + i]; 
                }
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
                    tmp[(M/2)*k + i] = 2.0*tmp[(M/2)*k + i]; 
                }
        }
    }

    out.swap(tmp);
}
