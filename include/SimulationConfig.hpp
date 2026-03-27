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
 * - `Nx`          : Number of x samples between [0,1].
 * - `Dim`         : Physical dimension D.
 * - `EpsNewton`   : Newton step damping/regularization.
 * - `PrecisionNewton` : Newton tolerance on residual.
 * - `SlowError`   : Damping in outer Newton loop.
 * - `MaxIterNewton` : Maximum Newton iterations.
 * - `Verbose/Debug/Converged` : Execution flags and previous status.
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
    vec_real F, Om, Pi, Psi;

    template <typename T>
    static H5::DataType get_hdf5_type();
    
    SimulationConfig(int Nt_, int Nx_)
        : Nt(Nt_), Nx(Nx_)
    {
        F.resize(Nx * Nt);
        Om.resize(Nx * Nt);
        Pi.resize(Nx * Nt);
        Psi.resize(Nx * Nt);
    }

    
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

        F.resize(Nx * Nt);
        Om.resize(Nx * Nt);
        Pi.resize(Nx * Nt);
        Psi.resize(Nx * Nt);

        readDataset(file, "F", F);
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

    void writeToHdf5(const std::string& filename) const {
        // Create a new HDF5 file (overwrite if exists)
        H5::H5File file(filename, H5F_ACC_TRUNC);

        hsize_t dims[2] = {
            static_cast<hsize_t>(Nx),
            static_cast<hsize_t>(Nt)
        };
        H5::DataSpace dataspace(2, dims);
        H5::DataSpace scalar(H5S_SCALAR);

        // datasets
        writeDataset(file, dataspace, "F", F);
        writeDataset(file, dataspace, "Om", Om);
        writeDataset(file, dataspace, "Pi", Pi);
        writeDataset(file, dataspace, "Psi", Psi);

        // attributes
        writeAttribute(file, scalar, "Nx", Nx);
        writeAttribute(file, scalar, "Nt", Nt);
        writeAttribute(file, scalar, "Delta", Delta);
        writeAttribute(file, scalar, "Dim", Dim);
        writeAttribute(file, scalar, "MaxIterNewton", MaxIterNewton);
        writeAttribute(file, scalar, "EpsNewton", EpsNewton);
        writeAttribute(file, scalar, "Converged", Converged);
        writeAttribute(file, scalar, "mismatchNorm", ErrorNorm);
        writeAttribute(file, scalar, "IterNewton", IterNewton);

        std::cout << "Output written to file. " << std::endl;
    }

    void readDataset(H5::H5File& file, const std::string& name, std::vector<double>& data)
    {
        H5::DataSet dataset = file.openDataSet(name);
        dataset.read(data.data(), H5::PredType::NATIVE_DOUBLE);
    }

    template <typename T>
    void readAttribute(H5::H5File& file, const std::string& name, T& value)
    {
        H5::Attribute attr = file.openAttribute(name);
        auto type = get_hdf5_type<T>();
        attr.read(type, &value);
    }

    static void writeDataset(H5::H5File& file,
                              const H5::DataSpace& dataspace,
                              const std::string& name,
                              const std::vector<double>& data)
    {
        H5::DataSet ds = file.createDataSet(
            name,
            H5::PredType::IEEE_F64LE,
            dataspace
        );
        ds.write(data.data(), H5::PredType::IEEE_F64LE);
    }

    template <typename T>
    static void writeAttribute(H5::H5File& file,
                                const H5::DataSpace& scalar,
                                const std::string& name,
                                const T& value)
    {
        auto type = get_hdf5_type<T>();
        H5::Attribute attr = file.createAttribute(name, type, scalar);
        attr.write(type, &value);
    }

    /// Print a human-readable configuration summary to stdout.
    void print_config()
    {
        std::cout << "Simulation configuration:" << std::endl;
        std::cout << "Ntau: " << Nt << std::endl;
        std::cout << "Dim: " << Dim << std::endl;
        std::cout << "Nx: " << Nx << std::endl;
        std::cout << "EpsNewton: " << EpsNewton << std::endl;
        std::cout << "PrecisionNewton: " << PrecisionNewton << std::endl;
        std::cout << "SlowError: " << SlowError << std::endl;
        std::cout << "MaxIterNewton: " << MaxIterNewton << std::endl;
        std::cout << "Verbose: " << Verbose << std::endl;
        std::cout << "Converged: " << Converged << std::endl;
        std::cout << "Delta: " << Delta << std::endl;
        std::cout << "f is not empty: " << !F.empty() << std::endl;
        std::cout << "Om is not empty: " << !Om.empty() << std::endl;
        std::cout << "Pi is not empty: " << !Pi.empty() << std::endl;
        std::cout << "Psi is not empty: " << !Psi.empty() << std::endl;
    }
};

template <>
inline H5::DataType SimulationConfig::get_hdf5_type<int>()
{
    return H5::PredType::NATIVE_INT;
}

// double
template <>
inline H5::DataType SimulationConfig::get_hdf5_type<double>()
{
    return H5::PredType::IEEE_F64LE;
}

// bool
template <>
inline H5::DataType SimulationConfig::get_hdf5_type<bool>()
{
    return H5::PredType::NATIVE_HBOOL;
}

