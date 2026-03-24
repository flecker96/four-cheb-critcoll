#include "common.hpp"
#include "SimulationConfig.hpp"
#include "NewtonSolver.hpp"

void readDataset(H5::H5File& file, const std::string& name, std::vector<double>& data);
void readAttribute(H5::H5File& file, const std::string& name, int& value);
void readAttribute(H5::H5File& file, const std::string& name, real_t& value);

int main(){

    std::string inputPath = "../to_hdf5/data.h5";
    std::string outputPath = "../to_hdf5/data_out.h5";

    SimulationConfig config = SimulationConfig::loadFromHDF5(inputPath);
    SimulationConfig result(config.Nt, config.Nx);

    NewtonSolver solver(config, result);
    solver.run();
    
    try {
        // Create a new HDF5 file (overwrite if exists)
        H5::H5File SimData(outputPath, H5F_ACC_TRUNC);

        // Define the dimensions of the dataset
        hsize_t dims[2] = { static_cast<hsize_t>(result.Nx), static_cast<hsize_t>(result.Nt) };
        H5::DataSpace dataspace(2, dims);

        // Create the dataset
        H5::DataSet datasetF = SimData.createDataSet(
            "F",                   // dataset name
            H5::PredType::NATIVE_DOUBLE, // data type
            dataspace
        );

        datasetF.write(result.F.data(), H5::PredType::NATIVE_DOUBLE);

        // Create the dataset
        H5::DataSet datasetOm = SimData.createDataSet(
            "Om",                   // dataset name
            H5::PredType::NATIVE_DOUBLE, // data type
            dataspace
        );

        datasetOm.write(result.Om.data(), H5::PredType::NATIVE_DOUBLE);

        // Create the dataset
        H5::DataSet datasetPi = SimData.createDataSet(
            "Pi",                   // dataset name
            H5::PredType::NATIVE_DOUBLE, // data type
            dataspace
        );

        datasetPi.write(result.Pi.data(), H5::PredType::NATIVE_DOUBLE);

        // Create the dataset
        H5::DataSet datasetPsi = SimData.createDataSet(
            "Psi",                   // dataset name
            H5::PredType::NATIVE_DOUBLE, // data type
            dataspace
        );

        datasetPsi.write(result.Psi.data(), H5::PredType::NATIVE_DOUBLE);

        // Write some attributes
        H5::DataSpace scalar_space(H5S_SCALAR);

        // Nx
        H5::Attribute attrNx = SimData.createAttribute(
            "Nx",
            H5::PredType::NATIVE_INT,
            scalar_space
        );
        attrNx.write(H5::PredType::NATIVE_INT, &result.Nx);

        // Nt
        H5::Attribute attrNt = SimData.createAttribute(
            "Nt",
            H5::PredType::NATIVE_INT,
            scalar_space
        );
        attrNt.write(H5::PredType::NATIVE_INT, &result.Nt);

        // Delta
        H5::Attribute attrDelta = SimData.createAttribute(
            "Delta",
            H5::PredType::IEEE_F64LE,
            scalar_space
        );
        attrDelta.write(H5::PredType::IEEE_F64LE, &result.Delta);

        // Dim
        H5::Attribute attrDim = SimData.createAttribute(
            "Dim",
            H5::PredType::IEEE_F64LE,
            scalar_space
        );
        attrDim.write(H5::PredType::IEEE_F64LE, &result.Dim);

        // MaxitsNewton
        H5::Attribute attrMaxIterNewton = SimData.createAttribute(
            "MaxIterNewton",
            H5::PredType::NATIVE_INT,
            scalar_space
        );
        attrMaxIterNewton.write(H5::PredType::NATIVE_INT, &result.MaxIterNewton);

        // epsNewton
        H5::Attribute attrepsNewton = SimData.createAttribute(
            "EpsNewton",
            H5::PredType::IEEE_F64LE,
            scalar_space
        );
        attrepsNewton.write(H5::PredType::IEEE_F64LE, &result.EpsNewton);

        // converged
        H5::Attribute attrconverged = SimData.createAttribute(
            "Converged",
            H5::PredType::NATIVE_HBOOL,
            scalar_space
        );
        attrconverged.write(H5::PredType::NATIVE_HBOOL, &result.Converged);

        std::cout << "Field successfully written to field.h5\n";
    }
    catch (H5::Exception& e) {
        std::cerr << "Error writing HDF5 file: " << e.getCDetailMsg() << "\n";
        return 1;
    }



    //*************************************************** */
    /*try {
        // Open the file and read in fields + parameters
        H5::H5File file("../to_hdf5/data.h5", H5F_ACC_RDONLY);

        readAttribute(file, "Delta", Delta);
        readAttribute(file, "Nx", Nx);
        readAttribute(file, "Nt", Nt);
        readAttribute(file, "Dim", Dim);

        f.resize(Nx * Nt);
        Om.resize(Nx * Nt);
        Pi.resize(Nx * Nt);
        Psi.resize(Nx * Nt);

        readDataset(file, "f", f);
        readDataset(file, "Om", Om);
        readDataset(file, "Pi", Pi);
        readDataset(file, "Psi", Psi);

    }
    catch(H5::Exception& e) {
        std::cerr << "HDF5 error: " << e.getCDetailMsg() << "\n";
        return 1;
    }*/

    

    /*Z.resize(Nt * Nx);

    SpectralTransformer st(Nt, Nx);


    for(int i=0;i<(Nt * Nx);i++){
        Z[i] = complex_t(f[i], Pi[i]);
    }


    st.forwardDCT_space(Z, Zout);
    st.backwardDCT_space2(Zout);


    for(int i=0;i<(Nt);i++){
        std::cout << std::setprecision(15)<< (Z[Nx-1+i*Nx] - Zout[Nx-1+i*Nx]).real() << ", " << (Z[Nx-1+i*Nx] - Zout[Nx-1+i*Nx]).imag() << std::endl;
    }
    */
    /*st.differentiate_x(Zout, Zout2);*/


    

    

    /*std::string inputPath{"data/simulation_4D_512.json"};
    //--------------------------------------------------------------------------
    // Derive data directory from input file (absolute path, parent folder).
    // Used to store outputs (e.g., benchmark_*.json).
    //--------------------------------------------------------------------------
    

    json result;
    SimulationConfig config = SimulationConfig::loadFromJson(inputPath);
    if (ignoreConverged) config.Converged = false;

    NewtonSolver solver(config, dataPath);
    result = solver.run();

    std::cout << "Result stored in file: " << inputPath << "\n\n";
    OutputWriter::writeJsonToFile(inputPath, result);

    std::cout << "Simulation finished successfully.\n\n";*/

}

void write_field(const std::string& filename,
                 const std::vector<double>& field,
                 int Nx, int Ny)
{
    H5::H5File file(filename, H5F_ACC_TRUNC);

    hsize_t dims[2];
    dims[0] = Nx;
    dims[1] = Ny;

    H5::DataSpace dataspace(2, dims);

    H5::DataSet dataset = file.createDataSet(
        "field",
        H5::PredType::NATIVE_DOUBLE,
        dataspace
    );

    dataset.write(field.data(), H5::PredType::NATIVE_DOUBLE);
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