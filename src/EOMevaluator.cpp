#include "EOMevaluator.hpp"

//------------------------------------------------------------------------------
// Ctor: initialize IRK stepper (Gauss–Legendre s=2) with shared initGen.
//------------------------------------------------------------------------------
EOMevaluator::EOMevaluator(int Ntau_, int Nx_, real_t Dim_, StatePacker& packer_)
    : Nt(Ntau_), Nx(Nx_), Dim(Dim_), packer(packer_)
{
    // Work buffers
    f.resize(Nt*Nx);
    Om.resize(Nt*Nx);
    Pi.resize(Nt*Nx); 
    Psi.resize(Nt*Nx);

    dtf.resize(Nt*Nx);
    dtOm.resize(Nt*Nx);
    dtPi.resize(Nt*Nx);
    dtPsi.resize(Nt*Nx);

    dxf.resize(Nt*Nx);
    dxOm.resize(Nt*Nx);
    dxPi.resize(Nt*Nx);
    dxPsi.resize(Nt*Nx);

    fRes.resize(Nt*Nx);
    OmRes.resize(Nt*Nx);
    PiRes.resize(Nt*Nx);
    PsiRes.resize(Nt*Nx);

}

//Output Y in x-t space
void EOMevaluator::ComputeResidual(const vec_complex& Yin, const real_t& Delta, const vec_real& x, 
                        const vec_real& zprime, vec_real& outputVec)
{
    packer.buildFields(Yin, Delta, f, Om, Pi, Psi, dtf, dtOm, dtPi, dtPsi, dxf, dxOm, dxPi, dxPsi);           
    for (size_t j=0; j<Nx; ++j) std::cout << std::setprecision(15) << Pi[j] << "," << std::endl;
    exit(0);
    for (size_t i=0; i<Nx; ++i)
    {
        for (size_t j=0; j<Nt; ++j)
        {
            fRes[Nx*j + i] = (zprime[i] * dxf[Nx*j + i] / f[Nx*j + i]) - x[i] * Om[Nx*j + i];
            OmRes[Nx*j + i] = x[i] * zprime[i] * dxOm[Nx*j + i] 
                                - Om[Nx*j + i] * (1.0 - Dim + x[i]*x[i]*Pi[Nx*j + i]*Pi[Nx*j + i] + x[i]*x[i]*x[i]*x[i]*Psi[Nx*j + i]*Psi[Nx*j + i]) 
                                + x[i] * x[i] * Om[Nx*j + i] * Om[Nx*j + i] 
                                - (Dim - 3.0) * (Pi[Nx*j + i]*Pi[Nx*j + i] + x[i]*x[i]*Psi[Nx*j + i]*Psi[Nx*j + i]);
            PiRes[Nx*j + i] = x[i] * zprime[i] * dxPi[Nx*j + i] + dtPi[Nx*j + i] + Pi[Nx*j + i] - f[Nx*j + i] * (x[i]*zprime[i]*dxPsi[Nx*j + i] + (Dim - 1.0 + x[i]*x[i]*Om[Nx*j + i])*Psi[Nx*j + i]);
            PsiRes[Nx*j + i] = x[i]*x[i]*zprime[i]*dxPsi[Nx*j + i] + x[i]*dtPsi[Nx*j + i] + 2.0*x[i]*Psi[Nx*j + i] - f[Nx*j + i] * (zprime[i] * dxPi[Nx*j + i] + x[i]*Om[Nx*j + i]*Pi[Nx*j + i]);
        }
    }

    packer.condenseResidual(fRes, OmRes, PiRes, PsiRes, outputVec);

}

