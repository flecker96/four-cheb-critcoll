#include "EOMevaluator.hpp"

EOMevaluator::EOMevaluator(int Ntau_, int Nx_, real_t Dim_, StatePacker& packer_)
    : Nt(Ntau_), Nx(Nx_), Dim(Dim_), packer(packer_)
{
    // Work buffers
    F.resize(Nt*Nx);
    Om.resize(Nt*Nx);
    Pi.resize(Nt*Nx); 
    Psi.resize(Nt*Nx);

    dtF.resize(Nt*Nx);
    dtOm.resize(Nt*Nx);
    dtPi.resize(Nt*Nx);
    dtPsi.resize(Nt*Nx);

    dxF.resize(Nt*Nx);
    dxOm.resize(Nt*Nx);
    dxPi.resize(Nt*Nx);
    dxPsi.resize(Nt*Nx);

    FRes.resize(Nt*Nx);
    OmRes.resize(Nt*Nx);
    PiRes.resize(Nt*Nx);
    PsiRes.resize(Nt*Nx);
    constr.resize(Nt*Nx);
    monitor_vec.resize(Nt*Nx/2);

}


void EOMevaluator::ComputeResidual(const vec_complex& Yin, const real_t& Delta, const vec_real& x, 
                        const vec_real& zprime, bool monitor, vec_real& outputVec)
{
    packer.buildFields(Yin, Delta, F, Om, Pi, Psi, dtF, dtOm, dtPi, dtPsi, dxF, dxOm, dxPi, dxPsi);

    
    for (size_t i=0; i<Nx; ++i)
    {
        for (size_t j=0; j<Nt; ++j)
        {
            FRes[Nx*j + i] = (1.0 - x[i])*zprime[i]*dxF[Nx*j + i] - F[Nx*j + i] - (1.0 + (1.0 - x[i])*F[Nx*j + i])  *Om[Nx*j + i] / 2.0;
            PiRes[Nx*j + i] = 2.0*x[i]*zprime[i]*dxPi[Nx*j + i] + Pi[Nx*j + i] + dtPi[Nx*j + i] 
                                - (1.0 + (1.0 - x[i])*F[Nx*j + i]) * ((Dim - 1.0 + x[i]*Om[Nx*j + i])*Psi[Nx*j + i] + 2.0*x[i]*zprime[i]*dxPsi[Nx*j + i]);
            PsiRes[Nx*j + i] = 2.0*x[i]*zprime[i]*dxPsi[Nx*j + i] + 2.0*Psi[Nx*j + i] + dtPsi[Nx*j + i]
                                - (1.0 + (1.0 - x[i])*F[Nx*j + i]) * (Om[Nx*j + i]*Pi[Nx*j + i] + 2.0*zprime[i]*dxPi[Nx*j + i]);
            constr[Nx*j + i] = 2.0*x[i]*zprime[i]*dxOm[Nx*j + i] 
                            + Om[Nx*j + i] * (Dim - 1.0 - x[i]*(Pi[Nx*j + i]*Pi[Nx*j + i] + x[i]*Psi[Nx*j + i]*Psi[Nx*j + i]))
                            + x[i]*Om[Nx*j + i]*Om[Nx*j + i] - (Dim - 3.0)*(Pi[Nx*j + i]*Pi[Nx*j + i] + x[i]*Psi[Nx*j + i]*Psi[Nx*j + i]);
        }
    }

    if (monitor)
    {
        for (size_t i=0; i<Nx; ++i)
        {
            for (size_t j=0; j<Nt; ++j)
            {
            // last equation of motion, for monitoring error
            OmRes[Nx*j + i] = 2.0*x[i]*zprime[i]*dxOm[Nx*j + i] + dtOm[Nx*j + i] 
                            - 2.0*(Dim - 3.0)*(1.0 + (1.0 - x[i])*F[Nx*j + i])*Pi[Nx*j + i]*Psi[Nx*j + i] 
                            + 2.0*Om[Nx*j + i]*(1.0 - x[i]*(1.0 + (1.0 - x[i])*F[Nx*j + i]) * Pi[Nx*j + i] * Psi[Nx*j + i]);
            }
        }

        //Dealias and extract components associated to OmRes
        packer.pack(FRes, OmRes, PiRes, PsiRes, monitor_vec);
        vec_real OmRes_reduced(monitor_vec.begin() + Nt*Nx/4, monitor_vec.begin() + 3*Nt*Nx/8);

        //Compute l2 norm
        real_t sum = std::transform_reduce(OmRes_reduced.cbegin(), OmRes_reduced.cend(), 0.0, std::plus{},
                                       [](auto x){return x*x;});
        real_t err_constr = std::sqrt(sum / static_cast<real_t>(OmRes_reduced.size()));

        std::cout << "Error in OmRes = " << err_constr << std::endl;

    }

    packer.pack(FRes, constr, PiRes, PsiRes, outputVec);

    /*std::ofstream out("test.txt");
    for (real_t x : FRes) {
        out << std::setprecision(15) << x << "\n";
    }
    out.close();
    exit(0);*/

}

