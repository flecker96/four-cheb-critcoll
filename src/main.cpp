#include "common.hpp"
#include "SimulationConfig.hpp"
#include "NewtonSolver.hpp"

int main(){

    std::string inputPath = "../to_hdf5/data.h5";
    std::string outputPath = "../to_hdf5/data_out.h5";

    SimulationConfig config = SimulationConfig::loadFromHDF5(inputPath);

    config.print_config();
    /*std::ifstream inFile("../Ftest.txt");
    for (double& val : config.F) {
        inFile >> val;
        if (inFile.fail()) {
            std::cerr << "Error reading value from file!\n";
            break;
        }
    }
    inFile.close();*/

    SimulationConfig result(config.Nt, config.Nx);
    NewtonSolver solver(config, result);
    solver.run();
    
    result.writeToHdf5(outputPath);

}
