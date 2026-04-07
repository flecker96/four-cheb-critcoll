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
    perturbedOutput.resize(Nnewton);
    perturbedOutputMinus.resize(Nnewton);
    Yin.resize(Nt * Nx);

    // xGrid
    xGrid.resize(Nx);
    z_prime.resize(Nx);
}

//------------------------------------------------------------------------------
// run (Serial)
//------------------------------------------------------------------------------
void NewtonSolver::run(std::string method)
{
    if (!Converged)
    {
        real_t err = 1.0, errOld = 1.0;              

        generateGrid();
        packer.pack(F, Om, Pi, Psi, in0);            //Build the vector in0 
        in0[3*Nt*Nx/8 + 2] = Delta;                  //Store Delta in slot for Re(Fc_2) (gauged to zero)

        vec_real in0old = in0;                       //Save always previous step. in case error increases return in0old

        for (size_t its=0; its<maxIts; ++its)
        { 
            std::cout << "=========   Newton iteration: " << its+1 << "  =========" << std::endl;
            auto toc_outer = std::chrono::high_resolution_clock::now();
            
            errOld = err;

            EOM(in0, out0, true);
           
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
                Delta = in0old[3*Nt*Nx/8 + 2];
                in0old[3*Nt*Nx/8 + 2] = 0.0;
                packer.NewtonToFields(in0old, F, Om, Pi, Psi);
                writeFinalOutput(its, errOld);
                std::cerr << "Mismatch increased – terminating Newton.\n";
                break;                    
            }

            // Solve J·dx = − out0 and update input
            vec_real dx(Nnewton);

            if (method == "Newton") {
                // Assemble J via finite differences
                mat_real J(Nnewton, vec_real(Nnewton));
                //EpsNewton = std::max(5E-8, 5E-5 * err);                 //adapt eps according to current residual; decrease during convergence
                std::cout << "EpsNewton = " << EpsNewton << std::endl;
                assembleJacobian(in0, out0, J);
                
                vec_real rhs = out0;
                real_t rcond = 1.0;                                               // Condition number of J
                std::for_each(rhs.begin(), rhs.end(), [](auto& e){ e *= -1.0; });
                std::cerr << "Inverting Jacobian...\n";
                solveLinearSystem(J, rhs, dx, rcond);

                // Print upper bound for relative error in dx, dx/x < rcond * drhs / rhs
                /*vec_real residual(Nnewton, 0.0);
                for (size_t i = 0; i < Nnewton; ++i)
                {
                    real_t sum = 0.0;
                    for (size_t j = 0; j < Nnewton; ++j)
                    {
                        sum += J[i][j] * dx[j];
                    }
                    residual[i] = sum + out0[i]; 
                }
                real_t errlinsol = computeL2Norm(residual);*/
                real_t dxnorm = computeL2Norm(dx);
                //std::cout << "Inversion error (upper bound): " << rcond * errlinsol / err << std::endl;
                std::cout << "Condition number: " << 1.0 / rcond << "\n";
                std::cout << "Norm of dx:       " << dxnorm << "\n";

            }
            else if (method == "Krylov"){

                // ... to be implemented if necessary
            }
            else 
                throw std::runtime_error("Need to provide method for Solver. ");
            
            // Do line search for updating step
            real_t lambda = 1.0;
            vec_real in_trial(Nnewton);
            vec_real out_trial(Nnewton);

            while (lambda > 1e-6)
            {
                for (size_t i = 0; i < Nnewton; ++i)
                    in_trial[i] = in0[i] + lambda * dx[i];

                EOM(in_trial, out_trial);
                real_t normF_trial = computeL2Norm(out_trial);

                //accept value of lambda if residual decreases
                if (normF_trial < err * (1.0))  
                    break;

                //otherwise halve lambda
                lambda *= 0.5; 
            }

            in0old = in0;

            for (size_t i=0; i<Nnewton; ++i) in0[i] += lambda * dx[i];

            auto tic_outer = std::chrono::high_resolution_clock::now();
            std::cout << "Time for Newton Iteration: " << static_cast<real_t>((tic_outer-toc_outer).count()) / 1e9
              << " s." << std::endl << std::endl;

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
        packer.pack(F, Om, Pi, Psi, in0);
        in0[3*Nt*Nx/8 + 2] = Delta;
        EOM(in0, out0, true);
        real_t err = computeL2Norm(out0);
        std::cout << "The solution is already marked as converged!" << std::endl;
        writeFinalOutput(0, err);
    }

}

 

//------------------------------------------------------------------------------
// EOM
// 1) Extract (F, Om, Pi, Psi, Δ) from packed input vector (spectral).
// 2) Evaluate EOM and reduce to outputvector (throw away dependent modes and dealias)
//------------------------------------------------------------------------------
void NewtonSolver::EOM(vec_real& inputVec, vec_real& outputVec, bool print_constr)
{
    // Extract Δ stored in the slot of in0 where Re(Fc_2) is
    Delta = inputVec[3*Nt*Nx/8 + 2];
    inputVec[3*Nt*Nx/8 + 2] = 0.0;

    packer.unpack(inputVec, Yin); //Yin is in spectral space, doubled
    inputVec[3*Nt*Nx/8 + 2] = Delta; //write Delta back to original slot

    //Take Yin (in spectral space), evaluate EOM and condense
    evaluator.ComputeResidual(Yin, Delta, xGrid, z_prime, print_constr, outputVec);

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
        
        perturbedInput[i] += EpsNewton * scale;

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
// build J[in0]*v (J is Jacobian) -- needed for Krylov
//------------------------------------------------------------------------------
void NewtonSolver::Jv(const vec_real& baseInput, const vec_real& trialInput, const vec_real& baseOutput, vec_real& trialOutput)
{

    vec_real perturbedInput = baseInput;
    vec_real perturbedInputMinus = baseInput;

    for (size_t i = 0; i < Nnewton; ++i){
        perturbedInput[i] += EpsNewton * trialInput[i];
        perturbedInputMinus[i] -= EpsNewton * trialInput[i];
    }

    EOM(perturbedInput, perturbedOutput);
    EOM(perturbedInputMinus, perturbedOutputMinus);

    //Build central finite difference (baseOutput actually not needed here.)
    for (size_t i = 0; i < Nnewton; ++i)
        trialOutput[i] = (perturbedOutput[i] - perturbedOutputMinus[i]) / (2.0*EpsNewton);

}

//------------------------------------------------------------------------------
// solveLinearSystem
// row-major dense copy made by memcpy per row.
//------------------------------------------------------------------------------
void NewtonSolver::solveLinearSystem(const mat_real& A_in, vec_real& rhs, vec_real& dx, real_t& rcond)
{
    vec_real A_flat(Nnewton * Nnewton);
    for (size_t i = 0; i < Nnewton; ++i)
        std::memcpy(&A_flat[i * Nnewton], A_in[i].data(), Nnewton * sizeof(real_t));

    //vec_real A_svd = A_flat;

    // compute 1-norm of A before dgesv overwrites it with LU factors
    double anorm = LAPACKE_dlange(LAPACK_ROW_MAJOR, '1',
                                  static_cast<int>(Nnewton),
                                  static_cast<int>(Nnewton),
                                  A_flat.data(),
                                  static_cast<int>(Nnewton));

    std::vector<lapack_int> ipiv(Nnewton);

    lapack_int info = LAPACKE_dgesv(LAPACK_ROW_MAJOR,
                                    static_cast<int>(Nnewton), 1,
                                    A_flat.data(), static_cast<int>(Nnewton),
                                    ipiv.data(), rhs.data(), 1);
    if (info != 0)
    {
        std::cerr << "ERROR: LAPACKE_dgesv failed with error code " << info << std::endl;
        std::exit(EXIT_FAILURE);
    }

    // A_flat now contains LU factors — use to estimate condition number
    LAPACKE_dgecon(LAPACK_ROW_MAJOR, '1',
                   static_cast<int>(Nnewton),
                   A_flat.data(), static_cast<int>(Nnewton),
                   anorm, &rcond);

    // SVD — only need singular values and VT, skip U
    /*int N = static_cast<int>(Nnewton);
    vec_real S(N);
    vec_real VT(N * N);
    vec_real U_dummy(N * N);

    lapack_int info_svd = LAPACKE_dgesdd(LAPACK_ROW_MAJOR, 'A',
                                         N, N,
                                         A_svd.data(), N,
                                         S.data(),
                                         U_dummy.data(), N,
                                         VT.data(), N);
    if (info_svd != 0)
    {
        std::cerr << "ERROR: LAPACKE_dgesdd failed with error code " << info_svd << std::endl;
        std::exit(EXIT_FAILURE);
    }

    std::ofstream outfile2("singvals.txt");
    outfile2 << std::scientific << std::setprecision(16);
    for (int i = 0; i < N; ++i)
        outfile2 << S[i] << "\n";
    outfile2.close();

    std::ofstream outfile("null_vector.txt");
    outfile << std::scientific << std::setprecision(16);
    for (int i = 0; i < N; ++i)
        outfile << VT[(N-1)*N + i] << " " << VT[(N-2)*N + i] << " " << VT[(N-3)*N + i] << "\n";
    
    outfile.close();
    std::cout << "Singular values written to file." << std::endl;
    exit(0);*/

    dx = rhs;
}

//------------------------------------------------------------------------------
// solveLinearSystem (Serial)
// Same LAPACK call; row-major dense copy made by memcpy per row.
//------------------------------------------------------------------------------
/*void NewtonSolver::solveLinearSystem(const mat_real& A_in, vec_real& rhs, vec_real& dx)
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
}*/


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
