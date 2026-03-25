#include "common.hpp"
#include "SimulationConfig.hpp"
#include "NewtonSolver.hpp"

int main(){

    std::string inputPath = "../to_hdf5/data_128_4D.h5";
    std::string outputPath = "../to_hdf5/data_out.h5";

    SimulationConfig config = SimulationConfig::loadFromHDF5(inputPath);
    SimulationConfig result(config.Nt, config.Nx);
    //for (int j=0; j<config.Nx; ++j) std::cout << std::setprecision(15) << config.F[j] << "," << std::endl;
    //exit(0);
    NewtonSolver solver(config, result);
    solver.run();
    
    result.writeToHdf5(outputPath);

}
