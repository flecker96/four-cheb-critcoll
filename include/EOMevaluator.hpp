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
    vec_real F, Om, Pi, Psi;
    vec_real dtF, dtOm, dtPi, dtPsi;
    vec_real dxF, dxOm, dxPi, dxPsi;
    vec_real FRes, OmRes, PiRes, PsiRes, constr, monitor_vec;

  public:

    //Cstor
    EOMevaluator(int Nt_, int Nx_, real_t Dim_, StatePacker& packer_);

    void ComputeResidual(const vec_complex& Yin, const real_t& Delta, const vec_real& x, const vec_real& xprime,  
                        bool monitor, vec_real& outputVec);
   
};