#include "Sampler.hpp"

Sampler::Sampler(size_t Nt_, size_t Nx_)
    : Nt(Nt_), Nx(Nx_)
    {
        //Work buffers
        F.resize(Nt*Nx);
        Om.resize(Nt*Nx);
        Pi.resize(Nt*Nx);
        Psi.resize(Nt*Nx);
    }

void Sampler::double_Nt(SimulationConfig& configIn, SimulationConfig& configOut)
{
    //Work arrays
    Fnew.resize(2*Nt*Nx);
    Omnew.resize(2*Nt*Nx);
    Pinew.resize(2*Nt*Nx);
    Psinew.resize(2*Nt*Nx);

    SpectralTransformer fft(Nt, Nx);
    SpectralTransformer fftnew(2*Nt, Nx);

    for (size_t i=0; i<Nt*Nx; ++i)
    {
        F[i] = configIn.F[i];
        Om[i] = configIn.Om[i];
        Pi[i] = configIn.Pi[i];
        Psi[i] = configIn.Psi[i];
    }

    fft.forwardFFT(F);
    fft.forwardFFT(Om);
    fft.forwardFFT(Pi);
    fft.forwardFFT(Psi);

    fft.doubleModes_t(F, Fnew);
    fft.doubleModes_t(Om, Omnew);
    fft.doubleModes_t(Pi, Pinew);
    fft.doubleModes_t(Psi, Psinew);

    fftnew.backwardFFT(Fnew);
    fftnew.backwardFFT(Omnew);
    fftnew.backwardFFT(Pinew);
    fftnew.backwardFFT(Psinew);

    //Write to new Configuration file
    for (size_t i=0; i<2*Nt*Nx; ++i)
    {
        configOut.F[i] = Fnew[i].real();
        configOut.Om[i] = Omnew[i].real();
        configOut.Pi[i] = Pinew[i].real();
        configOut.Psi[i] = Psinew[i].real();
    }

    configOut.Delta = configIn.Delta;
    configOut.Nt = 2 * configIn.Nt;
    configOut.Nx = configIn.Nx;
    configOut.Dim = configIn.Dim;
    configOut.MaxIterNewton = configIn.MaxIterNewton;
    configOut.EpsNewton = configIn.EpsNewton;
    configOut.PrecisionNewton = configIn.PrecisionNewton;
    configOut.Converged = false;

}

/*
void Sampler::double_Nx(SimulationConfig& configIn)
{
    
}
*/