#pragma once
/**
 * @file NewtonSolver.hpp
 * @brief Nonlinear solver (Newton–Raphson) for the boundary value problem
 *        arising in discretely self-similar scalar field collapse.
 *
 * @details
 * The NewtonSolver iteratively adjusts the initial data vector to enforce
 * matching conditions between left- and right-boundary Taylor expansions.
 * It drives the mismatch residual to zero using a Jacobian-based Newton
 * method. Parallel variants exist with OpenMP, MPI, or hybrid execution.
 *
 * Responsibilities:
 *  - Initialize near-boundary data via InitialConditionGenerator.
 *  - Perform shooting integration across the domain via ShootingSolver.
 *  - Assemble mismatch residuals at the matching point XMid.
 *  - Construct and solve the Newton system J·dx = -F.
 *  - Track convergence and write results/output files.
 */

#include "common.hpp"
#include "StatePacker.hpp"
#include "SimulationConfig.hpp"
#include "EOMevaluator.hpp"

/**
 * @class NewtonSolver
 * @brief Newton–Raphson driver enforcing boundary matching in DSS collapse.
 *
 * @section workflow Workflow
 * - Construct from a SimulationConfig and output folder.
 * - Generate τ-grid and boundary data.
 * - Run Newton iterations until mismatch < TolNewton or maxIts reached.
 * - Return results as JSON; optionally benchmark execution.
 *
 * @section parallelism Parallelism
 * - OpenMP: replicate local data structures for thread-level parallel shooting.
 * - MPI: distribute Jacobian assembly across ranks, collective solves.
 * - Hybrid: combine MPI+OpenMP.
 */
class NewtonSolver
{
  private:
    // ===== Simulation state =====
    SimulationConfig config;      ///< Simulation parameters (dimension, Ntau, Δ, etc.)
    SimulationConfig& result;
    size_t Nt, Nx;                  ///< Number of τ samples per period.
    size_t Nnewton;               ///< Number of Newton unknowns (typically 3*Ntau/4).
    size_t maxIts;                ///< Maximum Newton iterations.
    real_t Dim;                   ///< Physical spacetime dimension D.
    real_t Delta;                 ///< Echoing period Δ.
    real_t EpsNewton;             ///< Step length damping or regularization.
    real_t TolNewton;             ///< Convergence tolerance on mismatch norm.
    bool Debug;                   ///< Debug flag for verbose diagnostics.
    bool Verbose;                 ///< Verbose printing flag.
    bool Converged;               ///< True if Newton converged.

    // ===== Components =====
    vec_complex Yin;         ///< Packed state vectors at left/right boundaries.
    bool benchmark;                    ///< If true, run in benchmark mode.
    vec_real F, Om, Pi, Psi;             ///< Boundary input functions (τ-dependent).
    vec_real fRes, OmRes, PiRes, PsiRes;
    vec_real xGrid, xGridHalf, z_prime;                    ///< Radial grid points (left → right).
    vec_real in0, out0;                ///< Working vectors for shooting I/O.
    std::filesystem::path baseFolder;  ///< Path where outputs are written.
    StatePacker packer;               ///< Generates near-boundary Taylor expansions.
    EOMevaluator evaluator;

    
    void generateGrid();
    /**
     * @brief Evaluate equations of motion for a given configuration.
     * @param inputVec  Newton unknown vector.
     * @param outputVec Output residual vector.
     * @param fieldVals Optional JSON container for storing fields along integration.
     */
    void EOM(vec_real& inputVec, vec_real& outputVec, json* fieldVals=nullptr);

    /**
     * @brief Assemble finite-difference Jacobian of residuals.
     * @param baseInput  Current Newton input vector.
     * @param baseOutput Corresponding mismatch vector.
     * @param[out] jacobian Jacobian matrix J = ∂F/∂x.
     */
    void assembleJacobian(const vec_real& baseInput, const vec_real& baseOutput,
                          mat_real& jacobian);

    /**
     * @brief Solve linear system A·dx = rhs.
     * @param[in,out] A   Matrix (modified in-place if LAPACK factorization used).
     * @param[in,out] rhs Right-hand side vector (overwritten with solution).
     * @param[out] dx     Solution increment.
     */
    void solveLinearSystem(const mat_real& A, vec_real& rhs, vec_real& dx);
    
    /**
     * @brief Compute L² norm of a vector.
     * @param vc Input vector.
     * @return sqrt(Σ vᵢ²).
     */
    real_t computeL2Norm(const vec_real& vc);

  public:
    /**
     * @brief Construct NewtonSolver with simulation config and output path.
     * @param configIn   SimulationConfig object (Ntau, Nx, Dim, Δ, etc.).
     * @param configOut  SimulationConfig object for output.
     * @param benchmarkIn If true, enable benchmark mode.
     */
    NewtonSolver(SimulationConfig configIn, SimulationConfig& configOut, bool benchmarkIn=false);

    /**
     * @brief Run Newton solver until convergence or max iterations.
     * @param benchmark_result Optional JSON to collect benchmark statistics. (deprecated)
     * @return JSON result dictionary with fields (Converged, NewtonIts, Delta, etc.).
     */
    void run(json* benchmark_result=nullptr);

    /**
     * @brief Write final converged output and metadata to files.
     * @param newtonIts     Number of Newton iterations performed.
     * @param mismatchNorm  Final mismatch L² norm.
     */
    void writeFinalOutput(size_t newtonIts, real_t mismatchNorm);
};
