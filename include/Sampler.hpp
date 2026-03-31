#include "common.hpp"
#include "SimulationConfig.hpp"
#include "SpectralTransformer.hpp"

class Sampler
{
  private:
    // ===== Simulation state =====
    size_t Nt, Nx;                        
    
    // ======== Auxiliary ==============
    vec_complex F, Om, Pi, Psi;
    vec_complex Fnew, Omnew, Pinew, Psinew;

  public:

    //Cstor
    Sampler(size_t Nt_, size_t Nx_);

    void double_Nt(SimulationConfig& configIn, SimulationConfig& configOut);

    //not yet implemented..
    //void double_Nx(SimulationConfig& configIn);
   
};