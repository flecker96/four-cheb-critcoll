#include "common.hpp"
#include "Sampler.hpp"
#include "SimulationConfig.hpp"
#include "NewtonSolver.hpp"

int main(){

    std::string inputPath = "../sample_data/data_128_8D.h5";
    std::string outputPath = "../sample_data/data_out.h5";

    SimulationConfig config = SimulationConfig::loadFromHDF5(inputPath);

    /*//Double Nt
    SimulationConfig confignew(2*config.Nt, config.Nx);
    {
        Sampler samp(config.Nt, config.Nx);
        samp.double_Nt(config, confignew);
    }*/

    //config.Dim = 14.0;
    config.EpsNewton = 1E-8;
    config.Converged = false;
    config.PrecisionNewton = 1E-12;
    config.print_config();

    SimulationConfig result(config.Nt, config.Nx);
    NewtonSolver solver(config, result);
    solver.run();
    
    result.writeToHdf5(outputPath);

}
