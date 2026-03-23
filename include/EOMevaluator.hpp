#include "common.hpp"
#include "StatePacker.hpp"

class EOMevaluator
{
  private:
    // ===== Simulation state =====
    size_t Nt, Nx;                  ///< Number of τ samples per period.
    real_t Dim;

    // ===== Components =====
    StatePacker& packer;               ///< Generates near-boundary Taylor expansions.
    
    // ======== Auxiliary ==============
    vec_real f, Om, Pi, Psi;
    vec_real dtf, dtOm, dtPi, dtPsi;
    vec_real dxf, dxOm, dxPi, dxPsi;
    vec_real fRes, OmRes, PiRes, PsiRes;

  public:

    //Cstor
    EOMevaluator(int Nt_, int Nx_, real_t Dim_, StatePacker& packer_);

    void ComputeResidual(const vec_complex& Yin, const real_t& Delta, const vec_real& x, 
                        vec_real& outputVec);
   
};