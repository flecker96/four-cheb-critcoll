#pragma once
/**
 * @file StatePacker.hpp
 * @brief Construction of near-boundary initial data and conversions between
 *        field representations for the critical collapse solver.
 *
 * @details
 * The StatePacker creates high-order Taylor expansions of the
 * DSS fields near the left boundary (x≈0) and right boundary (x≈1). It also
 * provides utilities to pack/unpack parity-separated spectral fields and to
 * convert between
 *   - real-valued field triples (U,V,F)(τ), and
 *   - the complex spectral state vector Y( k ) used by the ODE/Shooting solvers.
 *
 * Internally, derivatives and inhomogeneous solves in τ are performed via a
 * SpectralTransformer (FFTW-based). The echoing period Δ is carried alongside
 * the transforms; by convention the implementation stores Δ in a component of
 * Y (see @ref FieldsToStateVector / @ref StateVectorToFields notes).
 *
 * Conventions:
 *  - Ntau real samples represent periodic functions over one echoing period τ∈[0,Δ).
 *  - Vectors `fc`, `psic`, `Up` have length Ntau.
 *  - Parity: two odd fields (Odd1, Odd2) and one even field (Even) appear in
 *    several packing utilities; in the usual gauge these correspond to
 *    Odd1≡U, Odd2≡V, Even≡F.
 */

#include "common.hpp"
#include "SpectralTransformer.hpp"

/**
 * @class StatePacker
 * @brief Builds fields from spectral representations and vice versa.
 *
 * @section responsibilities Responsibilities
 * - Convert (Pi,Psi,F,Om) ↔ spectral state vector Y.
 * - Pack/unpack parity-split spectral arrays to flat storage.
 *
 * @section requirements Requirements
 * - Nt ≥ 8 and even (typical spectral constraints).
 */
class StatePacker
{
  private:
    /// Number of τ grid points per period.
    size_t Nt;

    /// Number of x grid points.
    size_t Nx;

    /// Internal Newton vector length used by the shooting/Newton scheme (typically 3*Ntau/4).
    size_t Nnewton;

    /// Physical (rational) spacetime dimension D.
    real_t Dim;

    /// FFTW-based spectral engine (derivatives, integrals, solves).
    SpectralTransformer fft;

    // buffers
    vec_complex FF, OmF, PiF, PsiF; 
    vec_complex Ftmp, Omtmp, Pitmp, Psitmp;
    vec_complex tmp, Y, dtY, dxY;


  public:
    /**
     * @brief Construct a StatePacker.
     * @param Nt_   Number of τ samples per period.
     * @param Nx_     Number of x samples.
     * @param Dim_    Physical dimension D.
     *
     * @note The internal SpectralTransformer is initialized with (Nt_, Nx_).
     */
    StatePacker(size_t Nt_, size_t Nx_, real_t Dim_);

    void pack(const vec_real& F, const vec_real& Om, const vec_real& Pi, const vec_real& Psi, vec_real& vec);

    void unpack(const vec_real& inputVec, vec_complex& Y);

    void buildFields(const vec_complex& Y, real_t Delta, vec_real& f, vec_real& Om, vec_real& Pi , vec_real& Psi, 
                    vec_real& dtf, vec_real& dtOm, vec_real& dtPi , vec_real& dtPsi,
                    vec_real& dxf, vec_real& dxOm, vec_real& dxPi , vec_real& dxPsi);

    void NewtonToFields(const vec_real& vec, vec_real& F, vec_real& Om, vec_real& Pi, vec_real& Psi);

    /**
     * @brief Unpack complex spectral state Y to fields (no derivatives).
     * @param[in,out] Y  Spectral state vector.
     * @param[out] U     U(τ).
     * @param[out] V     V(τ).
     * @param[out] F     F(τ).
     */
    void StateVectorToFields(const vec_complex& Y, vec_real& Even1, vec_real& Even2,
                            vec_real& Odd1, vec_real& Odd2);
};
