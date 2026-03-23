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
    : Nt(Nt_), Nx(Nx_), Nnewton(Nt_* Nx_ / 2 - Nt_ / 4), Dim(Dim_), fft(Nt_, Nx_)
    {
        //Work buffers
        fF.resize(Nt*Nx/4);
        OmF.resize(Nt*Nx/4);
        PiF.resize(Nt*Nx/4);
        PsiF.resize(Nt*Nx/4);
        ftmp.resize(Nt*Nx);
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
    const vec_real& f, const vec_real& Om, const vec_real& Pi, const vec_real& Psi,
    const vec_real& x2, vec_real& vec)
{
    for(size_t i=0; i<Nt*Nx; i++){
        ftmp[i] = f[i];
        Omtmp[i] = Om[i];
        Pitmp[i] = Pi[i];
        Psitmp[i] = Psi[i];
    }

    fft.forwardFFT_time2(ftmp);
    fft.forwardDCT_space2(ftmp);
    fft.halveModes(ftmp, fF);
    fft.backwardChebHalf(fF);

    fft.forwardFFT_time2(Omtmp);
    fft.forwardDCT_space2(Omtmp);
    fft.halveModes(Omtmp, OmF);
    fft.backwardChebHalf(OmF);

    fft.forwardFFT_time2(Pitmp);
    fft.forwardDCT_space2(Pitmp);
    fft.halveModes(Pitmp, PiF);
    fft.backwardChebHalf(PiF);

    fft.forwardFFT_time2(Psitmp);
    fft.forwardDCT_space2(Psitmp);
    fft.halveModes(Psitmp, PsiF);
    fft.backwardChebHalf(PsiF);

    

    //now fill vec
    if (vec.size() != Nnewton)
        throw std::runtime_error("pack: wrong size of output vector");

    for (size_t i=0; i<Nx/2; ++i)
    {
        if (i==0) 
        {
            //Fill in fc, Omc, Pic and Psic
            for (size_t j=0; j<Nt/8; ++j)
            {
                //Pic
                vec[2*j + (Nt/4)*i] = PiF[(Nx/2)*(2*j + 1) + i].real();
                vec[2*j + 1 + (Nt/4)*i] = PiF[(Nx/2)*(2*j + 1) + i].imag();
                //Psic
                vec[2*j + Nt*i/4 + Nt*Nx/8] = PsiF[(Nx/2)*(2*j + 1) + i].real();
                vec[2*j + 1 + Nt*i/4 + Nt*Nx/8] = PsiF[(Nx/2)*(2*j + 1) + i].imag();
                //Om
                vec[2*j + Nt*i/4 + Nt*Nx/4] = OmF[(Nx/2)*(2*j) + i].real();
                vec[2*j + 1 + Nt*i/4 + Nt*Nx/4] = OmF[(Nx/2)*(2*j) + i].imag();
                //f
                vec[2*j + Nt*i/4 + 3*Nt*Nx/8] = fF[(Nx/2)*(2*j) + i].real();
                vec[2*j + 1 + Nt*i/4 + 3*Nt*Nx/8] = fF[(Nx/2)*(2*j) + i].imag();
            }

            
            //Nyquist
            vec[1 + Nt*i/4 + Nt*Nx/4] = OmF[(Nx/2)*(2*Nt/8) + i].real();  //Store HF cosine in imaginary part of zero mode
            vec[1 + Nt*i/4 + 3*Nt*Nx/8] = fF[(Nx/2)*(2*Nt/8) + i].real();

        }
        else
        {
            for (size_t j=0; j<Nt/8; ++j)
                {
                    //PI = (Pi - Pic)/x^2 for all values except at x=0
                    vec[2*j + (Nt/4)*i] = (PiF[(Nx/2)*(2*j + 1) + i].real() - PiF[(Nx/2)*(2*j + 1)].real() ) / (x2[i]*x2[i]);
                    vec[2*j + 1 + (Nt/4)*i] = (PiF[(Nx/2)*(2*j + 1) + i].imag() - PiF[(Nx/2)*(2*j + 1)].imag() ) / (x2[i]*x2[i]);
                    //PSI (same)
                    vec[2*j + Nt*i/4 + Nt*Nx/8] = (PsiF[(Nx/2)*(2*j + 1) + i].real() - PsiF[(Nx/2)*(2*j + 1)].real() ) / (x2[i]*x2[i]);
                    vec[2*j + 1 + Nt*i/4 + Nt*Nx/8] = (PsiF[(Nx/2)*(2*j + 1) + i].imag() - PsiF[(Nx/2)*(2*j + 1)].imag() ) / (x2[i]*x2[i]);
                    //Om (no bc)
                    vec[2*j + Nt*i/4 + Nt*Nx/4] = OmF[(Nx/2)*(2*j) + i].real();
                    vec[2*j + 1 + Nt*i/4 + Nt*Nx/4] = OmF[(Nx/2)*(2*j) + i].imag();
                    //f ; the point at the SSH (x=1) is not stored
                    if (i != Nx/2-1)
                    {
                    vec[2*j + Nt*i/4 + 3*Nt*Nx/8] = ((fF[(Nx/2)*(2*j) + i].real() - 1.0)/((1.0 - x2[i])*exp(x2[i])) - fF[(Nx/2)*(2*j)].real() + 1.0) / (x2[i]*x2[i]);
                    vec[2*j + 1 + Nt*i/4 + 3*Nt*Nx/8] = ((fF[(Nx/2)*(2*j) + i].imag() - 1.0)/((1.0 - x2[i])*exp(x2[i])) - fF[(Nx/2)*(2*j)].imag() + 1.0) / (x2[i]*x2[i]);
                    }
                }
            //Nyquist
            vec[1 + Nt*i/4 + Nt*Nx/4] = OmF[(Nx/2)*(2*Nt/8) + i].real();
            if (i != Nx/2-1){ 
                vec[1 + Nt*i/4 + 3*Nt*Nx/8] = ((fF[(Nx/2)*(2*Nt/8) + i].real() - 1.0)/((1.0 - x2[i])*exp(x2[i])) - fF[(Nx/2)*(2*Nt/8)].real() + 1.0) / (x2[i]*x2[i]);
            }
        }
    }
    
}

// This takes a Newton inputvector vec, storing 
// (RePic1, ImPic1, RePic3, ImPic3, ..., RePi(Nx/2-1)(Nt/4-1), ImPi(Nx/2-1)(Nt/4-1), [Psi], [Om], Refc0, ReFcN/2, 0, Imfc2, Refc4, Imfc4, ... )
// in halved variables in x-omega space and builds the fields. build fields with bcs. pad with zeros for antialiasing, output Y in spectral space for building derivatives next. 
void StatePacker::unpack(const vec_real& vec, const vec_real& x2, vec_complex& Y)
{
    //Fill halved vectors in x-omega space with zeros
    std::fill(PiF.begin(), PiF.end(), complex_t(0.0));
    std::fill(PsiF.begin(), PsiF.end(), complex_t(0.0));
    std::fill(OmF.begin(), OmF.end(), complex_t(0.0));
    std::fill(fF.begin(), fF.end(), complex_t(0.0));

    if (Y.size() != Nt*Nx)
        throw std::runtime_error("unpack: wrong size of output vector");

    for (size_t i=0; i<Nx/2; ++i){
            if (i==0)
            {
               for (size_t j=0; j<Nt/8; ++j)
                {
                // Pi: take values of Z (in x-om space) and build contiguous vector Pi (Nt/2 joined blocks of length Nx/2 each)
                // i=0 (x=0) is special: use stored value of pic for this
                PiF[(2*j + 1)*Nx/2 + i]        = complex_t(vec[2*j + i*Nt/4], vec[2*j + 1 + i*Nt/4]);
                PiF[(Nt/2 - 2*j - 1)*Nx/2 + i]    = std::conj(PiF[(2*j + 1)*Nx/2 + i]);

                // i=0 (x=0) is special: use stored value of psic for this
                PsiF[(2*j + 1)*Nx/2 + i]        = complex_t(vec[2*j + i*Nt/4 + Nt*Nx/8], vec[2*j+1 + i*Nt/4 + Nt*Nx/8]);
                PsiF[(Nt/2 - 2*j - 1)*Nx/2 + i]    = std::conj(PsiF[(2*j + 1)*Nx/2 + i]);

                // Om can be done in one go, no seperate regularity condition needed
                OmF[(2*j)*Nx/2 + i]  = complex_t(vec[2*j + i*Nt/4 + Nt*Nx/4], vec[2*j+1 + i*Nt/4 + Nt*Nx/4]);
                if (j!=0) OmF[(Nt/2 - 2*j)*Nx/2 + i] = std::conj(OmF[(2*j)*Nx/2 + i]);

                // i=0 (x=0) is special: use stored value of fc for this
                fF[(2*j)*Nx/2 + i] = complex_t(vec[2*j + i*Nt/4 + 3*Nt*Nx/8], vec[2*j+1 + i*Nt/4 + 3*Nt*Nx/8]);
                if (j!=0) fF[(Nt/2 - 2*j)*Nx/2 + i] = std::conj(fF[(2*j)*Nx/2 + i]);
                }
            }
            else if (i==(Nx/2-1)) 
            {
                for (size_t j=0; j<Nt/8; ++j)
                {
                    // Pi: take values of vec (in x-om space) and build contiguous vector Pi (Nt/2 joined blocks of length Nx/2 each)
                    // for the other i>0 implement regularity condition (\partial _x pi =0)
                    PiF[(2*j + 1)*Nx/2 + i]        = PiF[(2*j + 1)*Nx/2] + x2[i] * x2[i] * complex_t(vec[2*j + i*Nt/4], vec[2*j+1 + i*Nt/4]);
                    PiF[(Nt/2 - 2*j - 1)*Nx/2 + i]    = PiF[(Nt/2 - 2*j - 1)*Nx/2] + x2[i] * x2[i] * std::conj(PiF[(2*j + 1)*Nx/2 + i]);
                    // for the other i>0 implement regularity condition (\partial _x psi =0)
                    PsiF[(2*j + 1)*Nx/2 + i]        = PsiF[(2*j + 1)*Nx/2] + x2[i] * x2[i] * complex_t(vec[2*j + i*Nt/4 + Nt*Nx/8], vec[2*j+1 + i*Nt/4 + Nt*Nx/8]);
                    PsiF[(Nt/2 - 2*j - 1)*Nx/2 + i]    = PsiF[(Nt/2 - 2*j - 1)*Nx/2] + x2[i] * x2[i] * std::conj(PsiF[(2*j + 1)*Nx/2 + i]);
                    // Om can be done in one go, no seperate regularity condition needed
                    OmF[(2*j)*Nx/2 + i]  = complex_t(vec[2*j + i*Nt/4 + Nt*Nx/4], vec[2*j+1 + i*Nt/4 + Nt*Nx/4]);
                    if (j!=0) OmF[(Nt/2 - 2*j)*Nx/2 + i]    = std::conj(OmF[(2*j)*Nx/2 + i]);
                    //Set f to 1 at SSH, i.e. only switch on zero mode
                    if (j==0) fF[(2*j)*Nx/2 + i] = complex_t(1.0, 0.0);
                }
            }
            else {
                for (size_t j=0; j<Nt/8; ++j)
                {
                    // Pi: take values of vec (in x-om space) and build contiguous vector Pi (Nt/2 joined blocks of length Nx/2 each)
                    // for the other i>0 implement regularity condition (\partial _x pi =0)
                    PiF[(2*j + 1)*Nx/2 + i]        = PiF[(2*j + 1)*Nx/2] + x2[i] * x2[i] * complex_t(vec[2*j + i*Nt/4], vec[2*j+1 + i*Nt/4]);
                    PiF[(Nt/2 - 2*j - 1)*Nx/2 + i]    = PiF[(Nt/2 - 2*j - 1)*Nx/2] + x2[i] * x2[i] * std::conj(PiF[(2*j + 1)*Nx/2 + i]);
                    // for the other i>0 implement regularity condition (\partial _x psi =0)
                    PsiF[(2*j + 1)*Nx/2 + i]        = PsiF[(2*j + 1)*Nx/2] + x2[i] * x2[i] * complex_t(vec[2*j + i*Nt/4 + Nt*Nx/8], vec[2*j+1 + i*Nt/4 + Nt*Nx/8]);
                    PsiF[(Nt/2 - 2*j - 1)*Nx/2 + i]    = PsiF[(Nt/2 - 2*j - 1)*Nx/2] + x2[i] * x2[i] * std::conj(PsiF[(2*j + 1)*Nx/2 + i]);
                    // Om can be done in one go, no seperate regularity condition needed
                    OmF[(2*j)*Nx/2 + i]  = complex_t(vec[2*j + i*Nt/4 + Nt*Nx/4], vec[2*j+1 + i*Nt/4 + Nt*Nx/4]);
                    if (j!=0) OmF[(Nt/2 - 2*j)*Nx/2 + i]    = std::conj(OmF[(2*j)*Nx/2 + i]);
                    //build f by f=1+(1-x) * e^x * (fc - 1 + x^2 F(x))
                    fF[(2*j)*Nx/2 + i] = 1.0 + (1.0 - x2[i])*exp(x2[i])*(fF[(2*j)*Nx/2] - 1.0 + x2[i]*x2[i]*complex_t(vec[2*j + i*Nt/4 + 3*Nt*Nx/8], vec[2*j+1 + i*Nt/4 + 3*Nt*Nx/8]));
                    if (j!=0) 
                        fF[(Nt/2 - 2*j)*Nx/2 + i] = std::conj(fF[(2*j)*Nx/2 + i]);
                }
            }
        

        //Store Nyquist back where it belongs and make DC component purely real again
        OmF[(Nt/4)*Nx/2 + i] = complex_t(OmF[i].imag());
        OmF[i] = complex_t(OmF[i].real());

        fF[(Nt/4)*Nx/2 + i] = complex_t(fF[i].imag());
        fF[i] = complex_t(fF[i].real());        
    }
    
    //build state vector
    for (size_t i=0; i<Nt*Nx/4; ++i) {
        tmp[i] = fF[i] + PiF[i] + complex_t(0.0, 1.0)*(OmF[i] + PsiF[i]);
    }

    fft.forwardChebHalf(tmp);
    fft.doubleModes(tmp, Y);
}


void StatePacker::condenseResidual(const vec_real& fRes, const vec_real& OmRes, const vec_real& PiRes, const vec_real& PsiRes, vec_real& vec)
{
    //convert into complex vectors
    for(size_t i=0; i<Nt*Nx; i++){
        ftmp[i] = fRes[i];
        Omtmp[i] = OmRes[i];
        Pitmp[i] = PiRes[i];
        Psitmp[i] = PsiRes[i];
    }

    fft.forwardFFT_time2(ftmp);
    fft.forwardDCT_space2(ftmp);
    fft.halveModes(ftmp, fF);
    fft.backwardChebHalf(fF);

    fft.forwardFFT_time2(Omtmp);
    fft.forwardDCT_space2(Omtmp);
    fft.halveModes(Omtmp, OmF);
    fft.backwardChebHalf(OmF);

    fft.forwardFFT_time2(Pitmp);
    fft.forwardDCT_space2(Pitmp);
    fft.halveModes(Pitmp, PiF);
    fft.backwardChebHalf(PiF);

    fft.forwardFFT_time2(Psitmp);
    fft.forwardDCT_space2(Psitmp);
    fft.halveModes(Psitmp, PsiF);
    fft.backwardChebHalf(PsiF);

    //now fill vec
    if (vec.size() != Nnewton)
        throw std::runtime_error("condenseResidual: wrong size of output vector");

    for (size_t i=0; i<Nx/2; ++i)
    {
        for (size_t j=0; j<Nt/8; ++j)
            {
                //PI
                vec[2*j + (Nt/4)*i] = PiF[(Nx/2)*(2*j + 1) + i].real();
                vec[2*j + 1 + (Nt/4)*i] = PiF[(Nx/2)*(2*j + 1) + i].imag();
                //PSI
                vec[2*j + Nt*i/4 + Nt*Nx/8] = PsiF[(Nx/2)*(2*j + 1) + i].real();
                vec[2*j + 1 + Nt*i/4 + Nt*Nx/8] = PsiF[(Nx/2)*(2*j + 1) + i].imag();
                //Om 
                vec[2*j + Nt*i/4 + Nt*Nx/4] = OmF[(Nx/2)*(2*j) + i].real();
                vec[2*j + 1 + Nt*i/4 + Nt*Nx/4] = OmF[(Nx/2)*(2*j) + i].imag();
                //f ; the eom at the center (x=0) is not used
                if (i != 0)
                {
                vec[2*j + (Nt/4)*(i-1) + 3*Nt*Nx/8] = fF[(Nx/2)*(2*j) + i].real();
                vec[2*j + 1 + (Nt/4)*(i-1) + 3*Nt*Nx/8] = fF[(Nx/2)*(2*j) + i].imag();
                }
            }
        //Nyquist
        vec[1 + Nt*i/4 + Nt*Nx/4] = OmF[(Nx/2)*(2*Nt/8) + i].real();
        if (i != 0){ 
            vec[1 + (Nt/4)*(i-1) + 3*Nt*Nx/8] = fF[(Nx/2)*(2*Nt/8) + i].real();
        }
    }

}

void StatePacker::buildFields(const vec_complex& Yin, real_t Delta, vec_real& f, vec_real& Om, vec_real& Pi, vec_real& Psi, 
                            vec_real& dtf, vec_real& dtOm, vec_real& dtPi, vec_real& dtPsi, 
                            vec_real& dxf, vec_real& dxOm, vec_real& dxPi, vec_real& dxPsi)
{

    for (size_t i=0; i<Nt*Nx; ++i){
        Y[i] = Yin[i];
    }

    fft.differentiate_x(Y, dxY);
    fft.differentiate_t(Y, dtY, Delta);


    fft.backwardDCT_space2(Y);
    fft.backwardDCT_space2(dtY);
    fft.backwardDCT_space2(dxY);

    fft.backwardFFT_time2(Y);
    fft.backwardFFT_time2(dtY);
    fft.backwardFFT_time2(dxY);

    StateVectorToFields(Y, f, Om, Pi, Psi);
    StateVectorToFields(dtY, dtf, dtOm, dtPi, dtPsi);
    StateVectorToFields(dxY, dxf, dxOm, dxPi, dxPsi);

}

void StatePacker::NewtonToFields(const vec_real& vec, const vec_real& x2, vec_real& f, vec_real& Om, vec_real& Pi, vec_real& Psi)
{
    unpack(vec, x2, Y);
    fft.backwardDCT_space2(Y);
    fft.backwardFFT_time2(Y);

    StateVectorToFields(Y, f, Om, Pi, Psi);
}

//------------------------------------------------------------------------------
// StateVectorToFields 
// Takes Y = f + Pi + i(Om + Psi) and decomposes in to deparate fields
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


//------------------------------------------------------------------------------
// unpackSpectralFields
// Inverse of packSpectralFields: rebuild full complex spectra from Z by
// populating the half-spectrum and mirroring (Hermitian) before inverse FFT.
// Restores the special high-frequency cosine to both endpoints.
//------------------------------------------------------------------------------
void StatePacker::unpackSpectralFields(const vec_real& Z,
    vec_real& Odd1, vec_real& Odd2, vec_real& Even)
{
    vec_complex Odd1F(Nt/2, complex_t(0.0));
    vec_complex Odd2F(Nt/2, complex_t(0.0));
    vec_complex EvenF(Nt/2, complex_t(0.0));

    for (size_t j=0; j<Nnewton/6; ++j)
    {
        // Odd1: positive odd indices, then conjugate mirror
        Odd1F[2*j+1]               = complex_t(Z[2*j], Z[2*j+1]);
        Odd1F[Nt/2 - 2*j - 1]    = std::conj(Odd1F[2*j+1]);

        // Odd2
        Odd2F[2*j+1]               = complex_t(Z[2*j + Nnewton/3], Z[2*j+1 + Nnewton/3]);
        Odd2F[Nt/2 - 2*j - 1]    = std::conj(Odd2F[2*j+1]);

        // Even: even indices; mirror except j=0 (DC)
        EvenF[2*j]                 = complex_t(Z[2*j + 2*Nnewton/3], Z[2*j+1 + 2*Nnewton/3]);
        if (j != 0)
            EvenF[Nt/2 - 2*j]    = std::conj(EvenF[2*j]);
    }

    // Enforce pure-real DC (imag=0)
    EvenF[0] = complex_t(EvenF[0].real());

    // Expand to full spectrum
    fft.doubleModes(Odd1F, Odd1F);
    fft.doubleModes(Odd2F, Odd2F);
    fft.doubleModes(EvenF, EvenF);

    // Restore high-frequency cosine shared on symmetric endpoints
    EvenF[Nnewton/3] = complex_t(Z[2*Nnewton/3+1]) / 2.0;
    EvenF[Nnewton]   = complex_t(Z[2*Nnewton/3+1]) / 2.0;

    // Back to time domain
    fft.backwardFFT(Odd1F, Odd1);
    fft.backwardFFT(Odd2F, Odd2);
    fft.backwardFFT(EvenF, Even);
}

//------------------------------------------------------------------------------
// FieldsToStateVector
// Compose the complex *state* Y from (U,V,F) in time domain as
//   Y_j = U_j + i (V_j + F_j)
// then forward FFT and halve modes for compact spectral storage.
//------------------------------------------------------------------------------
void StatePacker::FieldsToStateVector(const vec_real& U, const vec_real& V,
    const vec_real& F, vec_complex& Y)
{
    Y.resize(Nt);

    for (size_t j=0; j<Nt; ++j)
    {
        Y[j] = complex_t(U[j], V[j] + F[j]);
    }

    fft.forwardFFT(Y, Y);
    fft.halveModes(Y, Y);
}

/*//------------------------------------------------------------------------------
// StateVectorToFields (with derivatives & IA2)
// Expand spectral Y to time-domain fields (U,V,F) and their τ-derivatives,
// and solve the auxiliary constraint for IA2(τ) at given X.
// Steps:
//   1) Double modes (rebuild full spectrum), fix special mode to preserve
//      parity/phasing, then inverse FFT → composite vector.
//   2) Split into U,V,F by symmetry: U,V odd; F even.
//   3) Differentiate in spectral space → dUdt,dVdt,dFdt.
//   4) Solve inhom constraint for IA2 with coefficients from U,V,F.
//------------------------------------------------------------------------------
void StatePacker::StateVectorToFields(vec_complex& Y, vec_real& U, vec_real& V,
     vec_real& F, vec_real& IA2, vec_real& dUdt, vec_real& dVdt, vec_real& dFdt, real_t X)
{
    vec_complex compVec1, compVec2;
    vec_real Coeff1(Nt), Coeff2(Nt);

    // Rebuild full spectrum and enforce phase fix on special mode
    Delta = Y[2].real();
    fft.doubleModes(Y, compVec1);
    compVec1[2] = complex_t(-compVec1[compVec1.size()-2].real(), compVec1[2].imag());

    // Differentiate in spectral space, then go to time domain
    fft.differentiate(compVec1, compVec2, Delta);
    fft.backwardFFT(compVec1, compVec1);
    fft.backwardFFT(compVec2, compVec2);

   
    // Serial branch (same logic)
    Delta = Y[2].real();
    fft.doubleModes(Y, compVec1);
    compVec1[2] = complex_t(-compVec1[compVec1.size()-2].real(), compVec1[2].imag());

    fft.differentiate(compVec1, compVec2, Delta);
    fft.backwardFFT(compVec1, compVec1);

    for (size_t j=0; j<Nt/2; ++j)
    {
        U[j] = 0.5 * (compVec1[j].real() - compVec1[j+Nt/2].real());
        V[j] = 0.5 * (compVec1[j].imag() - compVec1[j+Nt/2].imag());
        F[j] = 0.5 * (compVec1[j].imag() + compVec1[j+Nt/2].imag());
    }
    for (size_t j=0; j<Nt/2; ++j)
    {
        U[j+Nt/2] = -U[j];
        V[j+Nt/2] = -V[j];
        F[j+Nt/2] =  F[j];
    }

    fft.backwardFFT(compVec2, compVec2);
    for (size_t j=0; j<Nt/2; ++j)
    {
        dUdt[j] = 0.5 * (compVec2[j].real() - compVec2[j+Nt/2].real());
        dVdt[j] = 0.5 * (compVec2[j].imag() - compVec2[j+Nt/2].imag());
        dFdt[j] = 0.5 * (compVec2[j].imag() + compVec2[j+Nt/2].imag());
    }
    for (size_t j=0; j<Nt/2; ++j)
    {
        dUdt[j+Nt/2] = -dUdt[j];
        dVdt[j+Nt/2] = -dVdt[j];
        dFdt[j+Nt/2] =  dFdt[j];
    }

    for (size_t j=0; j<Nt; ++j)
    {
        Coeff1[j] = - ( ((X+F[j])*U[j]*U[j] + (X-F[j])*V[j]*V[j]) * std::pow((Dim - 2.0),3) / (8.0*X)
                         + (Dim - 3.0) );
        Coeff2[j] =  (Dim - 3.0);
    }
    fft.solveInhom(Coeff1, Coeff2, IA2, Delta);
}*/

