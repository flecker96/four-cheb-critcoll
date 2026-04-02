#include "Sampler.hpp"

//Currently this can only increase the modes (pad with zeros)
SimulationConfig changeModes(SimulationConfig& configIn, real_t fact, real_t facx)
{
    size_t Nt = configIn.Nt;
    size_t Nx = configIn.Nx;

    //Checks
    real_t Nt_test = configIn.Nt * fact;
    real_t Nx_test = configIn.Nx * facx; 

    if (!isInteger(Nt_test))
        throw std::runtime_error("changeModes: fact * Nt is not an integer!");
    if (!isInteger(Nx_test))
        throw std::runtime_error("changeModes: facx * Nx is not an integer!");
    
    size_t Ntnew = static_cast<size_t>(std::round(Nt_test));
    size_t Nxnew = static_cast<size_t>(std::round(Nx_test));

    if (Ntnew%8!=0)
        throw std::runtime_error("changeModes: New Nt must be divisible by 8.");
    if (Nxnew%2!=0)
        throw std::runtime_error("changeModes: New Nx must be divisible by 2.");

    std::cout << "Ntnew = " << Ntnew << ", Nxnew = " << Nxnew << std::endl;

    //Work arrays
    vec_complex F(Nt * Nx);
    vec_complex Om(Nt * Nx);
    vec_complex Pi(Nt * Nx);
    vec_complex Psi(Nt * Nx);
    vec_complex Fnew(Ntnew * Nxnew);
    vec_complex Omnew(Ntnew * Nxnew);
    vec_complex Pinew(Ntnew * Nxnew);
    vec_complex Psinew(Ntnew * Nxnew);

    SpectralTransformer fft(Nt, Nx);
    SpectralTransformer fftnew(Ntnew, Nxnew);

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

    fft.forwardCheb(F);
    fft.forwardCheb(Om);
    fft.forwardCheb(Pi);
    fft.forwardCheb(Psi);

    if (fact >= 1.0 && facx >= 1.0) 
    {
        fft.increaseModes_spec(F, Fnew, fact, facx);
        fft.increaseModes_spec(Om, Omnew, fact, facx);
        fft.increaseModes_spec(Pi, Pinew, fact, facx);
        fft.increaseModes_spec(Psi, Psinew, fact, facx);
    }
    else if (fact <= 1.0 && facx <= 1.0)
    {
        fft.decreaseModes_spec(F, Fnew, fact, facx);
        fft.decreaseModes_spec(Om, Omnew, fact, facx);
        fft.decreaseModes_spec(Pi, Pinew, fact, facx);
        fft.decreaseModes_spec(Psi, Psinew, fact, facx);
    }
    else 
    {
        throw std::runtime_error("changeModes: fact/facx should be both larger or smaller than 1.");
    }

    fftnew.backwardFFT(Fnew);
    fftnew.backwardFFT(Omnew);
    fftnew.backwardFFT(Pinew);
    fftnew.backwardFFT(Psinew);

    fftnew.backwardCheb(Fnew);
    fftnew.backwardCheb(Omnew);
    fftnew.backwardCheb(Pinew);
    fftnew.backwardCheb(Psinew);

    //Write to new Configuration file
    SimulationConfig configOut(Ntnew, Nxnew);

    for (size_t i=0; i<Ntnew*Nxnew; ++i)
    {
        configOut.F[i] = Fnew[i].real();
        configOut.Om[i] = Omnew[i].real();
        configOut.Pi[i] = Pinew[i].real();
        configOut.Psi[i] = Psinew[i].real();
    }

    configOut.Delta = configIn.Delta;
    configOut.Nt = Ntnew;
    configOut.Nx = Nxnew;
    configOut.Dim = configIn.Dim;
    configOut.MaxIterNewton = configIn.MaxIterNewton;
    configOut.EpsNewton = configIn.EpsNewton;
    configOut.PrecisionNewton = configIn.PrecisionNewton;
    configOut.Converged = false;

    return configOut;
}

