//==============================================================================
// StatePacker.cpp
// Builds spectral state vectors from boundary Taylor data for the DSS problem.
// - Left expansion (x≈0): 5th-order regular-center Taylor series.
// - Right expansion (x≈1): 3rd-order series near the self-similar horizon.
// - Spectral packing/unpacking of odd/even fields (halve/double modes).
// - Field <-> state-vector conversions and auxiliary constraint solve (IA2).
// All Fourier ops (FFT, differentiation, inhom solves) are delegated to `fft`.
//==============================================================================

#include "StatePacker.hpp"

//------------------------------------------------------------------------------
// Construct with:
//   Nt  : # grid points in τ (periodic).
//   Dim   : rational spacetime dimension d ∈ (3,4] in our use.
//   Delta : echoing period Δ (stored both here and inside FFT helper).
// Initializes the internal FFT helper with (Nt, Delta).
//------------------------------------------------------------------------------
StatePacker::StatePacker(size_t Nt_, size_t Nx_, real_t Dim_)
    : Nt(Nt_), Nx(Nx_), Nnewton(Nt_* Nx_ / 2), Dim(Dim_), fft(Nt_, Nx_)
    {
        //Work buffers
        FF.resize(Nt*Nx/4);
        OmF.resize(Nt*Nx/4);
        PiF.resize(Nt*Nx/4);
        PsiF.resize(Nt*Nx/4);
        Ftmp.resize(Nt*Nx);
        Omtmp.resize(Nt*Nx);
        Pitmp.resize(Nt*Nx);
        Psitmp.resize(Nt*Nx);

        tmp.resize(Nt*Nx/4);
        Y.resize(Nt*Nx);
        dtY.resize(Nt*Nx);
        dxY.resize(Nt*Nx);

    }

//------------------------------------------------------------------------------
// packSpectralFields
// Pack three *time-domain* series (Odd1, Odd2, Even) into a *real* vector Z
// in half-spectrum storage (after halving modes). Layout:
//   - Odd1: store positive odd modes' Re/Im interleaved.
//   - Odd2: idem, offset by Nnewton/3.
//   - Even: even modes' Re/Im interleaved, offset by 2*Nnewton/3,
//           plus a special high-frequency cosine scalar at Z[2*Nnewton/3+1].
//------------------------------------------------------------------------------
void StatePacker::pack(
    const vec_real& F, const vec_real& Om, const vec_real& Pi, const vec_real& Psi, vec_real& vec)
{
    for(size_t i=0; i<Nt*Nx; i++){
        Ftmp[i] = F[i];
        Pitmp[i] = Pi[i];
        Omtmp[i] = Om[i];
        Psitmp[i] = Psi[i];
    }

    fft.forwardFFT(Ftmp);
    fft.forwardCheb(Ftmp); 
    fft.halveModes(Ftmp, FF);
    fft.backwardChebHalf(FF);

    fft.forwardFFT(Pitmp);
    fft.forwardCheb(Pitmp);
    fft.halveModes(Pitmp, PiF);
    fft.backwardChebHalf(PiF);
    
    fft.forwardFFT(Omtmp);
    fft.forwardCheb(Omtmp);
    fft.halveModes(Omtmp, OmF);
    fft.backwardChebHalf(OmF);

    fft.forwardFFT(Psitmp);
    fft.forwardCheb(Psitmp);
    fft.halveModes(Psitmp, PsiF);
    fft.backwardChebHalf(PsiF);


    //now fill vec
    if (vec.size() != Nnewton)
        throw std::runtime_error("pack: wrong size of output vector");


    for (size_t i=0; i<Nx/2; ++i)
    {
        //Fill in
        for (size_t j=0; j<Nt/8; ++j)
        {
            //Pi
            vec[2*j + (Nt/4)*i] = PiF[(Nx/2)*(2*j + 1) + i].real();
            vec[2*j + 1 + (Nt/4)*i] = PiF[(Nx/2)*(2*j + 1) + i].imag();
            //Psi
            vec[2*j + (Nt/4)*i + Nt*Nx/8] = PsiF[(Nx/2)*(2*j + 1) + i].real();
            vec[2*j + 1 + (Nt/4)*i + Nt*Nx/8] = PsiF[(Nx/2)*(2*j + 1) + i].imag();
            //Om
            vec[2*j + (Nt/4)*i + Nt*Nx/4] = OmF[(Nx/2)*(2*j) + i].real();
            vec[2*j + 1 + (Nt/4)*i + Nt*Nx/4] = OmF[(Nx/2)*(2*j) + i].imag();
            //F
            vec[2*j + (Nt/4)*i + 3*Nt*Nx/8] = FF[(Nx/2)*(2*j) + i].real();
            vec[2*j + 1 + (Nt/4)*i + 3*Nt*Nx/8] = FF[(Nx/2)*(2*j) + i].imag();
        }

        //Nyquist
        vec[1 + (Nt/4)*i + Nt*Nx/4] = OmF[(Nx/2)*(Nt/4) + i].real();  //Store HF cosine in imaginary part of zero mode
        vec[1 + (Nt/4)*i + 3*Nt*Nx/8] = FF[(Nx/2)*(Nt/4) + i].real();
    }
    
}

// This takes a Newton inputvector vec, storing 
// (RePic1, ImPic1, RePic3, ImPic3, ..., RePi(Nx/2-1)(Nt/4-1), ImPi(Nx/2-1)(Nt/4-1), [Psi], [Om], Refc0, ReFcN/2, 0, Imfc2, Refc4, Imfc4, ... )
// in halved variables in x-omega space and builds the fields. build fields with bcs. pad with zeros for antialiasing, output Y in spectral space for building derivatives next. 
void StatePacker::unpack(const vec_real& vec, vec_complex& Y)
{
    //Fill halved vectors in x-omega space with zeros
    std::fill(PiF.begin(), PiF.end(), complex_t(0.0));
    std::fill(PsiF.begin(), PsiF.end(), complex_t(0.0));
    std::fill(OmF.begin(), OmF.end(), complex_t(0.0));
    std::fill(FF.begin(), FF.end(), complex_t(0.0));

    if (Y.size() != Nt*Nx)
        throw std::runtime_error("unpack: wrong size of output vector");
    
    for (size_t i=0; i<Nx/2; ++i){
            for (size_t j=0; j<Nt/8; ++j)
            {
                // Pi: take values of Z (in x-om space) and build contiguous vector Pi (Nt/2 joined blocks of length Nx/2 each)
                // i=0 (x=0) is special: use stored value of pic for this
                PiF[(2*j + 1)*Nx/2 + i]        = complex_t(vec[2*j + (Nt/4)*i], vec[2*j + 1 + (Nt/4)*i]);
                PiF[(Nt/2 - 2*j - 1)*Nx/2 + i]    = std::conj(PiF[(2*j + 1)*Nx/2 + i]);

                // i=0 (x=0) is special: use stored value of psic for this
                PsiF[(2*j + 1)*Nx/2 + i]        = complex_t(vec[2*j + (Nt/4)*i + Nt*Nx/8], vec[2*j+1 + (Nt/4)*i + Nt*Nx/8]);
                PsiF[(Nt/2 - 2*j - 1)*Nx/2 + i]    = std::conj(PsiF[(2*j + 1)*Nx/2 + i]);

                // Om can be done in one go, no seperate regularity condition needed
                OmF[(2*j)*Nx/2 + i]  = complex_t(vec[2*j + (Nt/4)*i + Nt*Nx/4], vec[2*j+1 + (Nt/4)*i + Nt*Nx/4]);
                if (j!=0) OmF[(Nt/2 - 2*j)*Nx/2 + i] = std::conj(OmF[(2*j)*Nx/2 + i]);

                // i=0 (x=0) is special: use stored value of fc for this
                FF[(2*j)*Nx/2 + i] = complex_t(vec[2*j + (Nt/4)*i + 3*Nt*Nx/8], vec[2*j+1 + (Nt/4)*i + 3*Nt*Nx/8]);
                if (j!=0) FF[(Nt/2 - 2*j)*Nx/2 + i] = std::conj(FF[(2*j)*Nx/2 + i]);
            }
        
        //Store Nyquist back where it belongs and make DC component purely real again
        OmF[(Nt/4)*Nx/2 + i] = complex_t(OmF[i].imag());
        OmF[i] = complex_t(OmF[i].real());

        FF[(Nt/4)*Nx/2 + i] = complex_t(FF[i].imag());
        FF[i] = complex_t(FF[i].real());        
       
    }

    //build state vector and double
    for (size_t i=0; i<Nt*Nx/4; ++i) {
        tmp[i] = FF[i] + PiF[i] + complex_t(0.0, 1.0)*(OmF[i] + PsiF[i]);
    }

    fft.forwardChebHalf(tmp);
    fft.doubleModes(tmp, Y);
}




void StatePacker::buildFields(const vec_complex& Yin, real_t Delta, vec_real& F, vec_real& Om, vec_real& Pi, vec_real& Psi, 
                            vec_real& dtF, vec_real& dtOm, vec_real& dtPi, vec_real& dtPsi, 
                            vec_real& dxF, vec_real& dxOm, vec_real& dxPi, vec_real& dxPsi)
{

    for (size_t i=0; i<Nt*Nx; ++i){
        Y[i] = Yin[i];
    }

    fft.differentiate_x(Y, dxY);
    fft.differentiate_t(Y, dtY, Delta);

    fft.backwardCheb(Y);
    fft.backwardCheb(dtY);
    fft.backwardCheb(dxY);

    fft.backwardFFT(Y);
    fft.backwardFFT(dtY);
    fft.backwardFFT(dxY);

    StateVectorToFields(Y, F, Om, Pi, Psi);
    StateVectorToFields(dtY, dtF, dtOm, dtPi, dtPsi);
    StateVectorToFields(dxY, dxF, dxOm, dxPi, dxPsi);
    
}

void StatePacker::NewtonToFields(const vec_real& vec, vec_real& F, vec_real& Om, vec_real& Pi, vec_real& Psi)
{
    unpack(vec, Y);
    fft.backwardCheb(Y);
    fft.backwardFFT(Y);

    StateVectorToFields(Y, F, Om, Pi, Psi);
}

//------------------------------------------------------------------------------
// StateVectorToFields 
// Takes Y = f + Pi + i(Om + Psi) and decomposes in to separate fields
//------------------------------------------------------------------------------
void StatePacker::StateVectorToFields(const vec_complex& Y, vec_real& Even1, vec_real& Even2,
     vec_real& Odd1, vec_real& Odd2)
{
     //take Y in (doubled) x-t space and compute fields
    for (size_t i=0; i<Nx; ++i)
    {
        for (size_t j=0; j<Nt/2; ++j)
        {
            Even1[Nx*j + i] = 0.5 * (Y[Nx*j + i].real() + Y[Nx*(j + Nt/2) + i].real());
            Odd1[Nx*j + i] = 0.5 * (Y[Nx*j + i].real() - Y[Nx*(j + Nt/2) + i].real());
            Even2[Nx*j + i] = 0.5 * (Y[Nx*j + i].imag() + Y[Nx*(j + Nt/2) + i].imag());
            Odd2[Nx*j + i] = 0.5 * (Y[Nx*j + i].imag() - Y[Nx*(j + Nt/2) + i].imag());
        }
        for (size_t j=0; j<Nt/2; ++j)
        {
            Even1[Nx*(j + Nt/2) + i] = Even1[Nx*j + i];
            Odd1[Nx*(j + Nt/2) + i] = - Odd1[Nx*j + i];
            Even2[Nx*(j + Nt/2) + i] = Even2[Nx*j + i];
            Odd2[Nx*(j + Nt/2) + i] = - Odd2[Nx*j + i];
        }
    }
}


