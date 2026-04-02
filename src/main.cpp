#include "common.hpp"
#include "Sampler.hpp"
#include "SimulationConfig.hpp"
#include "NewtonSolver.hpp"

int main(){

    std::string inputPath = "../to_hdf5/Nt64/data_128_14D.h5";
    std::string outputPath = "../to_hdf5/data_out.h5";

    SimulationConfig configIn = SimulationConfig::loadFromHDF5(inputPath);

    //Upsample functions (padding with zeros)
    real_t fact = 9.0/8.0, facx = 1.0;
    SimulationConfig config = changeModes(configIn, fact, facx);
    
    //config.Dim = 5.5;
    config.EpsNewton = 1E-6;
    config.Converged = false;
    config.PrecisionNewton = 1E-13;
    config.print_config();

    SimulationConfig result(config.Nt, config.Nx);
    NewtonSolver solver(config, result);
    solver.run();
    
    result.writeToHdf5(outputPath);

}
