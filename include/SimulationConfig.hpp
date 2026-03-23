#pragma once
/**
 * @file SimulationConfig.hpp
 * @brief Lightweight data structures for loading and organizing simulation parameters.
 *
 * @details
 * - **SimulationConfig**: POD-style container that initializes itself from a JSON object
 *   (or file) and exposes all parameters required by the solver stack (Newton/Shooting/IRK).
 * - **SimulationSuite**: Helper to manage a *multi*-dimension input dictionary and derive
 *   the ordered list of dimensions to run (optionally skipping already converged entries
 *   and/or reversing order with look-back seeding).
 */

#include "common.hpp"

/**
 * @struct SimulationConfig
 * @brief Single-simulation configuration (one dimension D, one Ntau, etc.).
 *
 * @section fields Key Fields
 * - `Ntau`        : Number of τ samples per period.
 * - `Dim`         : Physical dimension D.
 * - `XLeft/XMid/XRight` : Spatial match points (domain decomposition).
 * - `EpsNewton`   : Newton step damping/regularization.
 * - `PrecisionNewton` : Newton tolerance on residual.
 * - `SlowError`   : Slowly-updated error metric (diagnostics).
 * - `MaxIterNewton` : Maximum Newton iterations.
 * - `Verbose/Debug/Converged` : Execution flags and previous status.
 * - `DebugNx/DebugNtau` : Save every nth spatial and time point.
 * - `NLeft/NRight`: Spatial grid resolution left/right of the match.
 * - `SchemeIRK`: Scheme (order) for implicit RK method (e.g. 1(1),2(4),3(6)).
 * - `PrecisionIRK`: Tolerance for implicit RK step solves.
 * - `MaxIterIRK`  : Maximum Newton iterations per IRK step.
 * - `Delta`       : Echoing period.
 * - `fc/psic/Up`  : Periodic input arrays for boundary expansions.
 */
struct SimulationConfig
{
    int Nt, Nx;
    real_t Dim;
    real_t EpsNewton, PrecisionNewton, SlowError;
    int    MaxIterNewton, IterNewton;
    bool   Verbose, Debug, Converged;
    size_t DebugNx, DebugNtau;
    real_t Delta, ErrorNorm;
    vec_real f, Om, Pi, Psi;

    
    SimulationConfig(int Nt_, int Nx_)
        : Nt(Nt_), Nx(Nx_)
    {
        f.resize(Nx * Nt);
        Om.resize(Nx * Nt);
        Pi.resize(Nx * Nt);
        Psi.resize(Nx * Nt);
    }

    /**
     * @brief Construct from a JSON object.
     *
     * Expected layout (keys must exist; some initial condition entries may be null):
     * ```
     * {
     *   "Ntau": ..., "Dim": ...,
     *   "XLeft": ..., "XMid": ..., "XRight": ...,
     *   "EpsNewton": ..., "PrecisionNewton": ..., "SlowError": ...,
     *   "MaxIterNewton": ..., "Verbose": ..., "Debug": ..., "DebugNx": ..., "DebugNtau": ...,
     *   "Converged": ..., "NLeft": ..., "NRight": ...,
     *   "PrecisionIRK": ..., "SchemeIRK": ..., "MaxIterIRK": ...,
     *   "Initial_Conditions": {
     *     "Delta": <float or null>,
     *     "fc":    <array or null>,
     *     "psic":  <array or null>,
     *     "Up":    <array or null>
     *   }
     * }
     * ```
     */
    SimulationConfig(H5::H5File& file)
    {
        readAttribute(file, "Delta", Delta);
        readAttribute(file, "Nx", Nx);
        readAttribute(file, "Nt", Nt);
        readAttribute(file, "Dim", Dim);
        readAttribute(file, "MaxIterNewton", MaxIterNewton);
        readAttribute(file, "EpsNewton", EpsNewton);
        readAttribute(file, "TolNewton", PrecisionNewton);
        readAttribute(file, "slowErr", SlowError);
        readAttribute(file, "Converged", Converged);

        f.resize(Nx * Nt);
        Om.resize(Nx * Nt);
        Pi.resize(Nx * Nt);
        Psi.resize(Nx * Nt);

        readDataset(file, "f", f);
        readDataset(file, "Om", Om);
        readDataset(file, "Pi", Pi);
        readDataset(file, "Psi", Psi);
    }

    

    static SimulationConfig loadFromHDF5(const std::string& filename)
    {
        try {
            H5::H5File file(filename, H5F_ACC_RDONLY);
            return SimulationConfig(file);
        } catch (const H5::Exception& e) {
            throw std::runtime_error("Failed to load HDF5 file: " + filename);
        }
    }

    void readDataset(H5::H5File& file, const std::string& name, std::vector<double>& data)
    {
        H5::DataSet dataset = file.openDataSet(name);
        dataset.read(data.data(), H5::PredType::NATIVE_DOUBLE);
    }

    void readAttribute(H5::H5File& file, const std::string& name, int& value)
    {
        H5::Attribute attr = file.openAttribute(name);
        attr.read(H5::PredType::NATIVE_INT, &value);
    }

    void readAttribute(H5::H5File& file, const std::string& name, real_t& value)
    {
        H5::Attribute attr = file.openAttribute(name);
        attr.read(H5::PredType::NATIVE_DOUBLE, &value);
    }

    void readAttribute(H5::H5File& file, const std::string& name, bool& value)
    {
        H5::Attribute attr = file.openAttribute(name);
        attr.read(H5::PredType::NATIVE_HBOOL, &value);
    }

    /**
     * @brief Save a JSON object (e.g., results) to file.
     * @note Writes the JSON as-is; caller controls formatting/indentation.
     */
    static void saveToJson(const std::string& filename, json& simRes)
    {
        std::ofstream file(filename);
        file << simRes;
    }

    /// Print a human-readable configuration summary to stdout.
    void print_config()
    {
        std::cout << "Simulation configuration:" << std::endl;
        std::cout << "Ntau: " << Nt << std::endl;
        std::cout << "Dim: " << Dim << std::endl;
        std::cout << "XLeft: " << Nx << std::endl;
        std::cout << "EpsNewton: " << EpsNewton << std::endl;
        std::cout << "PrecisionNewton: " << PrecisionNewton << std::endl;
        std::cout << "SlowError: " << SlowError << std::endl;
        std::cout << "MaxIterNewton: " << MaxIterNewton << std::endl;
        std::cout << "Verbose: " << Verbose << std::endl;
        std::cout << "Debug: " << Debug << std::endl;
        std::cout << "DebugNx: " << DebugNx << std::endl;
        std::cout << "DebugNtau: " << DebugNtau << std::endl;
        std::cout << "Converged: " << Converged << std::endl;
        std::cout << "Delta: " << Delta << std::endl;
        std::cout << "f is not empty: " << !f.empty() << std::endl;
        std::cout << "Om is not empty: " << !Om.empty() << std::endl;
        std::cout << "Pi is not empty: " << !Pi.empty() << std::endl;
        std::cout << "Psi is not empty: " << !Psi.empty() << std::endl;
    }
};


