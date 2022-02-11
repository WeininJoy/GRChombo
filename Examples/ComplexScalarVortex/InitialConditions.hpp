/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef INITIALCONDITIONS_HPP_
#define INITIALCONDITIONS_HPP_

#include <cmath> 
#include "Cell.hpp"
#include "Coordinates.hpp"
#include "MatterCCZ4RHS.hpp"
#include "ComplexScalarField.hpp"
#include "Tensor.hpp"
#include "UserVariables.hpp" //This files needs NUM_VARS - total no. components
#include "VarsTools.hpp"
#include "simd.hpp"

//! Class which creates the initial conditions

class InitialConditions
{

  public:
    //! A structure for the input params for scalar field properties and initial
    //! conditions
    struct params_t
    {
        double field_amplitude; //!< Amplitude of bump in initial SF bubble
        double scalar_mass;
        std::array<double, CH_SPACEDIM>
            center;   //!< Centre of perturbation in initial SF bubble
    };

    //! The constructor for the class
    InitialConditions(params_t a_params, double a_dx)
        : m_dx(a_dx), m_params(a_params)
    {
    }
    

    //! Function to compute the value of all the initial vars on the grid
    template <class data_t> void compute(Cell<data_t> current_cell) const
    {
        // where am i?
        Coordinates<data_t> coords(current_cell, m_dx, m_params.center);
        const data_t x = coords.x;
        const double y = coords.y;
        const double z = coords.z;
        data_t rho2 = simd_max(x * x + y * y, 1e-12);
	
        data_t phi_Re =   m_params.field_amplitude  * x * exp(-rho2 / (0.5*0.5));
        data_t phi_Im =   m_params.field_amplitude  * y * exp(-rho2 / (0.5*0.5));
        data_t Pi_Re  =   m_params.field_amplitude  * m_params.scalar_mass * y ;
        data_t Pi_Im  = - m_params.field_amplitude  * m_params.scalar_mass * x ;

        // Store the initial values of the variables
        current_cell.store_vars(phi_Re, c_phi_Re);
        current_cell.store_vars(Pi_Re, c_Pi_Re);
        current_cell.store_vars(phi_Im, c_phi_Im);
        current_cell.store_vars(Pi_Im, c_Pi_Im);
    }

  protected:
    double m_dx;
    const params_t m_params; //!< The matter initial condition params
};

#endif /* INITIALCONDITIONS_HPP_ */
