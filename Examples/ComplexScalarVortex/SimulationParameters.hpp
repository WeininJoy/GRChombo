/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef SIMULATIONPARAMETERS_HPP_
#define SIMULATIONPARAMETERS_HPP_

// General includes
#include "GRParmParse.hpp"
#include "SimulationParametersBase.hpp"

// Problem specific includes:
#include "ComplexStaticVortex.hpp"
#include "ScalarPotential.hpp"
#include "InitialConditions.hpp"

class SimulationParameters : public SimulationParametersBase
{
  public:
    SimulationParameters(GRParmParse &pp) : SimulationParametersBase(pp)
    {
        // read the problem specific params
        read_params(pp);
    }

    void read_params(GRParmParse &pp)
    {
        pp.load("regrid_length", regrid_length, L);
        pp.load("G_Newton", G_Newton, 1.0); // for now the example neglects backreaction
        // vortex data
        pp.load("vortex_amplitude", vortex_params.Amp, 1.0);
        vortex_params.center = center;
        pp.load("winding_n", vortex_params.n, 1);
        // initial data
        pp.load("scalar_mass", initial_params.scalar_mass, 1.0);
        initial_params.field_amplitude = vortex_params.Amp;
        initial_params.center = center;
        // potential data
        potential_params.scalar_mass = initial_params.scalar_mass;
    }

    // Problem specific parameters
    double regrid_length, G_Newton;
    ComplexStaticVortex::params_t vortex_params;
    ScalarPotential::params_t potential_params;
    InitialConditions::params_t initial_params;
};

#endif /* SIMULATIONPARAMETERS_HPP_ */
