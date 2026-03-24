//==============================================================================
// NewtonSolver.cpp
// Newton–Kantorovich solver for the critical-collapse boundary value problem.
// Responsibilities:
//   • Hold simulation config and working buffers (spectral/state/grids).
//   • Build left/right Taylor initial data and call the ShootingSolver.
//   • Assemble finite-difference Jacobian via repeated “shoot” evaluations.
//   • Solve J·dx = −res with LAPACK (row-major), apply line damping, iterate.
//   • Support Serial / OpenMP / MPI / Hybrid backends and benchmarking.
//==============================================================================

#include "NewtonSolver.hpp"

//------------------------------------------------------------------------------
// Ctor: cache config scalars & sizes, allocate work arrays, construct helpers.
// If MPI/Hybrid, also create a contiguous MPI type to send/recv vectors.
//------------------------------------------------------------------------------
NewtonSolver::NewtonSolver(SimulationConfig configIn, SimulationConfig& configOut, bool benchmarkIn)
    : config(configIn), result(configOut), Nt(configIn.Nt), Nx(configIn.Nx), Nnewton(configIn.Nt * configIn.Nx/2 - configIn.Nt/4), maxIts(configIn.MaxIterNewton), 
    Dim(configIn.Dim), Delta(configIn.Delta), slowErr(configIn.SlowError), EpsNewton(configIn.EpsNewton), TolNewton(configIn.PrecisionNewton),
    Debug(configIn.Debug), Verbose(configIn.Verbose), Converged(configIn.Converged), benchmark(benchmarkIn), 
    f(configIn.f), Om(configIn.Om), Pi(configIn.Pi), Psi(configIn.Psi),
    packer(configIn.Nt, configIn.Nx, configIn.Dim),
    evaluator(configIn.Nt, configIn.Nx, configIn.Dim, packer)
{
    // Work buffers
    in0.resize(Nnewton);
    out0.resize(Nnewton);
    Yin.resize(Nt * Nx);

    // xGrid
    xGrid.resize(Nx);
    z_prime.resize(Nx);
    xGridHalf.resize(Nx / 2);
}

//------------------------------------------------------------------------------
// run (Serial)
//------------------------------------------------------------------------------
void NewtonSolver::run(json* benchmark_result)
{
    if (!Converged)
    {
        real_t fac = 1.0, err = 1.0, errOld = 1.0;
        int mismatch_increase_count {};

        generateGrid();
        packer.pack(f, Om, Pi, Psi, xGridHalf, in0);            //Build the vector in0 
        in0[3*Nt*Nx/8 + 2] = Delta;                             //Store Delta in slot for Re(fc_2) (gauged to zero)


        vec_real in0old = in0;

        for (size_t its=0; its<maxIts; ++its)
        {
            
            std::cout << "Newton iteration: " << its+1 << std::endl;

            errOld = err;
            
            shoot(in0, out0);
           
            err = computeL2Norm(out0);
            std::cout << "Mismatch norm: " << err << std::endl;

            if (err<TolNewton)
            {
                Converged = true;
                std::cout << "The solution has converged!" << std::endl << std::endl;
                writeFinalOutput(its, err);
                break;
            }

            // Safeguard against divergence
            if (std::log10(err) > std::log10(errOld)) mismatch_increase_count++; 

            if (mismatch_increase_count >= 5)
            {
                Converged = true;
                Delta = in0old[3*Nt*Nx/8 + 2];
                in0old[3*Nt*Nx/8 + 2] = 0.0;
                packer.NewtonToFields(in0old, xGridHalf, f, Om, Pi, Psi);
                writeFinalOutput(its, errOld);
                std::cerr << "Mismatch increased – terminating Newton.\n";
                break;                    
            }

            // Assemble J via finite differences
            mat_real J(Nnewton, vec_real(Nnewton));
            assembleJacobian(in0, out0, J);
            

            // Solve J·dx = −res and update input with damping
            vec_real dx(Nnewton);
            vec_real rhs = out0;
            std::for_each(rhs.begin(), rhs.end(), [](auto& e){ e *= -1.0; });

            solveLinearSystem(J, rhs, dx);

            fac = std::min(1.0, slowErr / err);

            if (std::log10(err) < std::log10(errOld))
            {
                in0old = in0;                 // store trial point
                mismatch_increase_count = 0;
            }
            
            for (size_t i=0; i<Nnewton; ++i) in0[i] += fac * dx[i];

        }

        if (!Converged)
        {
            std::cerr << "Newton method did not converge in " << maxIts << " iterations. For dimension: " << Dim << std::endl;
            std::exit(EXIT_FAILURE);
        }
    }
    else
    {
        // Already converged: recompute residual for logging and store output
        generateGrid();
        packer.pack(f, Om, Pi, Psi, xGridHalf, in0);
        in0[3*Nt*Nx/8 + 2] = Delta;
        shoot(in0, out0);
        real_t err = computeL2Norm(out0);
        std::cout << "The solution is already marked as converged!" << std::endl;
        writeFinalOutput(0, err);
    }

}



//------------------------------------------------------------------------------
// shoot
// 1) Extract (Up, ψc, fc, Δ) from packed input vector (spectral).
// 2) Build left/right boundary data (Taylor) and call IRK shooter.
// 3) Transform mismatch state → fields → pack back to spectral output.
// If `fieldVals` is provided, the shooter exports intermediate field slices.
//------------------------------------------------------------------------------
void NewtonSolver::shoot(vec_real& inputVec, vec_real& outputVec, json* fieldVals)
{
    // Extract Δ stored in the slot of in0 where Re(fc_2) is
    Delta = inputVec[3*Nt*Nx/8 + 2];
    inputVec[3*Nt*Nx/8 + 2] = 0.0;

    packer.unpack(inputVec, xGridHalf, Yin); //Yin is in spectral space, doubled
    inputVec[3*Nt*Nx/8 + 2] = Delta;

    //Take Yin (in spectral space), evaluate EOM and condense
    evaluator.ComputeResidual(Yin, Delta, xGrid, z_prime, outputVec);
}



//------------------------------------------------------------------------------
// generateGrid
// Build a monotone grid in x by uniform spacing in the logit variable
//   z = log(x) - log(1-x).
// Left block: [XLeft, XMid] with NLeft segments; Right: [XMid, XRight].
//------------------------------------------------------------------------------
void NewtonSolver::generateGrid()
{
    // currently they are mapped to x in a linear fashion; may be changed according to needs
    for (size_t k = 0; k < Nx; ++k)
    {
        real_t z = std::cos(M_PI * k / (static_cast<real_t>(Nx) - 1.0));
        xGrid[k] = (1.0 - z) / 2.0;
        z_prime[k] = - 2.0;
    }

    for (size_t k = 0; k < Nx/2; ++k)
    {
        real_t z = std::cos(M_PI * k / (static_cast<real_t>(Nx/2) - 1.0));
        xGridHalf[k] = (1 - z) / 2;
    }

}


//------------------------------------------------------------------------------
// assembleJacobian (Serial/OpenMP)
// Serial: loop over columns; OpenMP: thread-local shooters/generators and
//         private buffers to build columns in parallel safely.
// Note: we fill J by columns (jacobian[j][i]) for LAPACK row-major layout later.
//------------------------------------------------------------------------------
void NewtonSolver::assembleJacobian(const vec_real& baseInput, const vec_real& baseOutput, mat_real& jacobian)
{
    Verbose = false;

    std::cout << "Starting to assemble Jacobian: " << std::endl << std::endl;

    auto toc_outer = std::chrono::high_resolution_clock::now();

    for (size_t i=0; i<Nnewton; ++i)
    {
        auto toc_inner = std::chrono::high_resolution_clock::now();

        vec_real perturbedInput = baseInput;
        perturbedInput[i] += EpsNewton;

        vec_real perturbedOutput(Nnewton);
        shoot(perturbedInput, perturbedOutput);

        for (size_t j=0; j<Nnewton; ++j)
        {
            jacobian[j][i] = (perturbedOutput[j] - baseOutput[j]) / EpsNewton;
        }

        //if (config.Verbose)
        //{
            auto tic_inner = std::chrono::high_resolution_clock::now();
            std::cout << "Varying parameter: " << i+1 << "/" << Nnewton
                      << " in " << static_cast<real_t>((tic_inner-toc_inner).count()) / 1e9
                      << " s." << std::endl;
        //}
    }

    auto tic_outer = std::chrono::high_resolution_clock::now();
    std::cout << "Time for Newton Iteration: " << static_cast<real_t>((tic_outer-toc_outer).count()) / 1e9
              << " s." << std::endl << std::endl;

    Verbose = config.Verbose;
}


//------------------------------------------------------------------------------
// solveLinearSystem (Serial/OpenMP)
// Same LAPACK call; row-major dense copy made by memcpy per row.
//------------------------------------------------------------------------------
void NewtonSolver::solveLinearSystem(const mat_real& A_in, vec_real& rhs, vec_real& dx)
{
    vec_real A_flat(Nnewton * Nnewton);
    for (size_t i=0; i<Nnewton; ++i)
    {
        std::memcpy(&A_flat[i * Nnewton], A_in[i].data(), Nnewton * sizeof(real_t));
    }

    std::vector<lapack_int> ipiv(Nnewton);

    lapack_int info = LAPACKE_dgesv(LAPACK_ROW_MAJOR, static_cast<int32_t>(Nnewton), 1, A_flat.data(),
                                    static_cast<int32_t>(Nnewton), ipiv.data(), rhs.data(), 1);
    if (info != 0)
    {
        std::cerr << "ERROR: LAPACKE_dgesv failed with error code " << info << std::endl;
        std::exit(EXIT_FAILURE);
    }
    dx = rhs;
}


//------------------------------------------------------------------------------
// computeL2Norm: sqrt( mean( v.^2 ) ). Used as mismatch metric.
//------------------------------------------------------------------------------
real_t NewtonSolver::computeL2Norm(const vec_real& vc)
{
    real_t sum = std::transform_reduce(vc.cbegin(), vc.cend(), 0.0, std::plus{},
                                       [](auto x){return x*x;});
    return std::sqrt(sum / static_cast<real_t>(vc.size()));
}

//------------------------------------------------------------------------------
// writeFinalOutput
// Fill `resultDict` with solver metadata, final mismatch, and initial data.
// On parallel builds only rank 0 prints the confirmation message.
//------------------------------------------------------------------------------
void NewtonSolver::writeFinalOutput(size_t newtonIts, real_t mismatchNorm)
{
    result.f = f;
    result.Om = Om;
    result.Pi = Pi;
    result.Psi = Psi;

    result.Nt = Nt;
    result.Nx = Nx;
    result.Dim = Dim;
    result.Converged = Converged;

    result.EpsNewton = EpsNewton;
    result.PrecisionNewton = TolNewton;
    result.SlowError = slowErr;
    result.MaxIterNewton = maxIts;
    result.IterNewton = newtonIts;
    result.ErrorNorm = mismatchNorm;

    std::cout << "Final result stored in simulation dictionary for dimension D = " << Dim << "." << std::endl << std::endl;
   
}
