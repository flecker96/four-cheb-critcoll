//==============================================================================
// NewtonSolver.cpp
// Newton–Kantorovich solver for the critical-collapse boundary value problem.
// Responsibilities:
//   • Hold simulation config and working buffers (spectral/state/grids).
//   • Evaluate equations of motion for given field config.
//   • Assemble finite-difference Jacobian via repeated “shoot” evaluations.
//   • Solve J·dx = −res with LAPACK (row-major), apply line damping, iterate.
//==============================================================================

#include "NewtonSolver.hpp"

//------------------------------------------------------------------------------
// Ctor: cache config scalars & sizes, allocate work arrays, construct helpers.
//------------------------------------------------------------------------------
NewtonSolver::NewtonSolver(SimulationConfig configIn, SimulationConfig& configOut, bool benchmarkIn)
    : config(configIn), result(configOut), Nt(configIn.Nt), Nx(configIn.Nx), Nnewton(configIn.Nt * configIn.Nx/2), maxIts(configIn.MaxIterNewton), 
    Dim(configIn.Dim), Delta(configIn.Delta), EpsNewton(configIn.EpsNewton), TolNewton(configIn.PrecisionNewton),
    Debug(configIn.Debug), Verbose(configIn.Verbose), Converged(configIn.Converged), benchmark(benchmarkIn), 
    F(configIn.F), Om(configIn.Om), Pi(configIn.Pi), Psi(configIn.Psi),
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
        real_t err = 1.0, errOld = 1.0;

        generateGrid();
        packer.pack(F, Om, Pi, Psi, in0);            //Build the vector in0 
        in0[3*Nt*Nx/8 + 2] = Delta;                  //Store Delta in slot for Re(Fc_2) (gauged to zero)

        for (size_t its=0; its<maxIts; ++its)
        {
            
            std::cout << "Newton iteration: " << its+1 << std::endl;
            auto toc_outer = std::chrono::high_resolution_clock::now();
            
            errOld = err;

            EOM(in0, out0);
           
            err = computeL2Norm(out0);
            std::cout << "Mismatch norm: " << err << std::endl;
            std::cout << "Delta = " << in0[3*Nt*Nx/8 + 2] << std::endl;

            if (err<TolNewton)
            {
                Converged = true;
                Delta = in0[3*Nt*Nx/8 + 2];
                in0[3*Nt*Nx/8 + 2] = 0.0;
                packer.NewtonToFields(in0, F, Om, Pi, Psi);
                std::cout << "The solution has converged!" << std::endl << std::endl;
                writeFinalOutput(its, err);
                break;
            }

            if (err >= errOld && its > 0)
            {
                Converged = true;
                Delta = in0[3*Nt*Nx/8 + 2];
                in0[3*Nt*Nx/8 + 2] = 0.0;
                packer.NewtonToFields(in0, F, Om, Pi, Psi);
                writeFinalOutput(its, errOld);
                std::cerr << "Mismatch increased – terminating Newton.\n";
                break;                    
            }

            // Assemble J via finite differences
            mat_real J(Nnewton, vec_real(Nnewton));
            assembleJacobian(in0, out0, J);
            

            // Solve J·dx = −res and update input
            vec_real dx(Nnewton);
            vec_real rhs = out0;
            std::for_each(rhs.begin(), rhs.end(), [](auto& e){ e *= -1.0; });
            std::cerr << "Inverting Jacobian...\n";
            solveLinearSystem(J, rhs, dx);

            // Print error of linear solver
            vec_real residual(Nnewton, 0.0);
            for (size_t i = 0; i < Nnewton; ++i)
            {
                double sum = 0.0;
                for (size_t j = 0; j < Nnewton; ++j)
                {
                    sum += J[i][j] * dx[j];
                }
                residual[i] = sum + out0[i]; 
            }
            double errlinsol = computeL2Norm(residual);
            std::cout << "Inverted with error: " << errlinsol << std::endl;

            
            // Do line search for updating step
            double lambda = 1.0;
            vec_real in_trial(Nnewton);
            vec_real out_trial(Nnewton);

            while (lambda > 1e-6)
            {
                for (size_t i = 0; i < Nnewton; ++i)
                    in_trial[i] = in0[i] + lambda * dx[i];

                EOM(in_trial, out_trial);
                double normF_trial = computeL2Norm(out_trial);

                //accept value of lambda if residual decreses
                if (normF_trial < err)  
                    break;

                //otherwise halve lambda
                lambda *= 0.5; 
            }

            for (size_t i=0; i<Nnewton; ++i) in0[i] += lambda * dx[i];

            auto tic_outer = std::chrono::high_resolution_clock::now();
            std::cout << "Time for Newton Iteration: " << static_cast<real_t>((tic_outer-toc_outer).count()) / 1e9
              << " s." << std::endl << std::endl;

        }

        if (!Converged)
        {
            std::cout << in0[3*Nt*Nx/8 + 2] << std::endl;
            std::cerr << "Newton method did not converge in " << maxIts << " iterations. For dimension: " << Dim << std::endl;
            std::exit(EXIT_FAILURE);
        }
    }
    else
    {
        // Already converged: recompute residual for logging and store output
        generateGrid();
        packer.pack(F, Om, Pi, Psi, in0);
        in0[3*Nt*Nx/8 + 2] = Delta;
        EOM(in0, out0);
        real_t err = computeL2Norm(out0);
        std::cout << "The solution is already marked as converged!" << std::endl;
        writeFinalOutput(0, err);
    }

}



//------------------------------------------------------------------------------
// EOM
// 1) Extract (F, Om, Pi, Psi, Δ) from packed input vector (spectral).
// 2) Evaluate EOM and reduce to outputvector (throw away dependent modes and dealias)
// If `fieldVals` is provided, the shooter exports intermediate field slices. (not working currently)
//------------------------------------------------------------------------------
void NewtonSolver::EOM(vec_real& inputVec, vec_real& outputVec, json* fieldVals)
{
    // Extract Δ stored in the slot of in0 where Re(Fc_2) is
    Delta = inputVec[3*Nt*Nx/8 + 2];
    inputVec[3*Nt*Nx/8 + 2] = 0.0;

    packer.unpack(inputVec, Yin); //Yin is in spectral space, doubled
    inputVec[3*Nt*Nx/8 + 2] = Delta; //write Delta back to original slot

    //Take Yin (in spectral space), evaluate EOM and condense
    evaluator.ComputeResidual(Yin, Delta, xGrid, z_prime, outputVec);
}



//------------------------------------------------------------------------------
// generateGrid
// Build Chebyshev grid and transform to fundamental domain
//------------------------------------------------------------------------------
void NewtonSolver::generateGrid()
{
    // currently they are mapped to x as x=(1-z) / 2 may be changed according to needs
    for (size_t k = 0; k < Nx; ++k)
    {
        real_t z = std::cos(M_PI * k / (static_cast<real_t>(Nx) - 1.0));
        xGrid[k] = (1.0 - z) / 2.0;
        z_prime[k] = - 2.0;
    }
}


//------------------------------------------------------------------------------
// assembleJacobian
// Note: we fill J by columns (jacobian[j][i]) for LAPACK row-major layout later.
//------------------------------------------------------------------------------
void NewtonSolver::assembleJacobian(const vec_real& baseInput, const vec_real& baseOutput, mat_real& jacobian)
{
    Verbose = false;

    std::cout << "Starting to assemble Jacobian: " << std::endl << std::endl;

    for (size_t i=0; i<Nnewton; ++i)
    {
        auto toc_inner = std::chrono::high_resolution_clock::now();

        vec_real perturbedInput = baseInput;
        real_t scale = std::max(1.0, abs(baseInput[i]));

        //std::cout << scale << std::endl;
        perturbedInput[i] += EpsNewton * scale;

        vec_real perturbedOutput(Nnewton);
        EOM(perturbedInput, perturbedOutput);

        for (size_t j=0; j<Nnewton; ++j)
        {
            jacobian[j][i] = (perturbedOutput[j] - baseOutput[j]) / (EpsNewton * scale);
        }

        if (i%1000 == 0)
        {
            auto tic_inner = std::chrono::high_resolution_clock::now();
            std::cout << "Varied " << i << "/" << Nnewton
                      << " parameters. Time per par: " << static_cast<real_t>((tic_inner-toc_inner).count()) / 1e9
                      << " s." << std::endl;
        }
    }

    Verbose = config.Verbose;
}


//------------------------------------------------------------------------------
// solveLinearSystem (Serial)
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
// Fill `result` with solver metadata, final mismatch, and initial data.
//------------------------------------------------------------------------------
void NewtonSolver::writeFinalOutput(size_t newtonIts, real_t mismatchNorm)
{
    result.F = F;
    result.Om = Om;
    result.Pi = Pi;
    result.Psi = Psi;
    result.Delta = Delta;

    result.Nt = Nt;
    result.Nx = Nx;
    result.Dim = Dim;
    result.Converged = Converged;

    result.EpsNewton = EpsNewton;
    result.PrecisionNewton = TolNewton;
    result.MaxIterNewton = maxIts;
    result.IterNewton = newtonIts;
    result.ErrorNorm = mismatchNorm;

    std::cout << "Final result stored in simulation dictionary for dimension D = " << Dim << "." << std::endl << std::endl;
   
}
