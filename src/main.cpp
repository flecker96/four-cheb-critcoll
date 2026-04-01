#include "common.hpp"
#include "Sampler.hpp"
#include "SimulationConfig.hpp"
#include "NewtonSolver.hpp"

int main(){

    std::string inputPath = "../to_hdf5/data_128_14D.h5";
    std::string outputPath = "../to_hdf5/data_out.h5";

    SimulationConfig config = SimulationConfig::loadFromHDF5(inputPath);

    //Double Nt
    SimulationConfig confignew(2*config.Nt, config.Nx);
    {
        Sampler samp(config.Nt, config.Nx);
        samp.double_Nt(config, confignew);
    }

    //config.Dim = 14.0;
    confignew.EpsNewton = 1E-5;
    confignew.Converged = false;
    config.PrecisionNewton = 1E-2;
    confignew.print_config();

    SimulationConfig result(confignew.Nt, confignew.Nx);
    NewtonSolver solver(confignew, result);
    solver.run();
    
    result.writeToHdf5(outputPath);

}
