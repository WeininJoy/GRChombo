/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef REALSCALARCIRCULATION_HPP_
#define REALSCALARCIRCULATION_HPP_

#include "ADMFixedBGVars.hpp"
#include "CCZ4Geometry.hpp"
#include "Cell.hpp"
#include "Coordinates.hpp"
#include "FourthOrderDerivatives.hpp"
#include "GRInterval.hpp"
#include "Tensor.hpp"
#include "UserVariables.hpp" //This files needs NUM_VARS - total number of components
#include "VarsTools.hpp"
#include "simd.hpp"
#include "ScalarPotential.hpp"


//! Calculates the circulation: (-y*v_x + x*v_y) with type matter_t and writes it to the grid
template <class matter_t, class background_t> class ComplexScalarCirculation
{
    // Use the variable definition in the matter class
    template <class data_t>
    using MatterVars = typename matter_t::template Vars<data_t>;

    // Now the non grid ADM vars
    template <class data_t> using MetricVars = ADMFixedBGVars::Vars<data_t>;

  protected:
    const FourthOrderDerivatives
        m_deriv; //!< An object for calculating derivatives of the variables
    const matter_t m_matter;         //!< The matter object
    const double m_dx;               //!< The matter object
    const background_t m_background; //!< The matter object
    const ScalarPotential m_potential;        //!< The params object
    const std::array<double, CH_SPACEDIM> m_center;

  public:
    ComplexScalarCirculation(matter_t a_matter, background_t a_background, ScalarPotential a_potential, const double a_dx,
                   const std::array<double, CH_SPACEDIM> a_center)
        : m_matter(a_matter), m_potential(a_potential), m_deriv(a_dx), m_dx(a_dx),
          m_background(a_background), m_center(a_center)
    {
    }

    template <class data_t> void compute(Cell<data_t> current_cell) const
    {
        // copy data from chombo gridpoint into local variables, and calc 1st
        // derivs
        const auto vars = current_cell.template load_vars<MatterVars>();
        const auto d1 = m_deriv.template diff1<MatterVars>(current_cell);
	      Coordinates<data_t> coords(current_cell, m_dx, m_center);

        // get the metric vars
        MetricVars<data_t> metric_vars;
        m_background.compute_metric_background(metric_vars, current_cell);

        using namespace TensorAlgebra;
        const auto gamma_UU = compute_inverse_sym(metric_vars.gamma);
        const auto chris_phys =
            compute_christoffel(metric_vars.d1_gamma, gamma_UU);
        const emtensor_t<data_t> emtensor = m_matter.compute_emtensor(
            vars, metric_vars, d1, gamma_UU, chris_phys.ULL);

        // assign values of density in output box
        data_t v_x = 1.0/( pow(vars.phi_Re,2.0)+ pow(vars.phi_Im,2.0)) * (vars.phi_Im*d1.phi_Re[0] - vars.phi_Re*d1.phi_Im[0]);
        data_t v_y = 1.0/( pow(vars.phi_Re,2.0)+ pow(vars.phi_Im,2.0)) * (vars.phi_Im*d1.phi_Re[1] - vars.phi_Re*d1.phi_Im[1]);
        current_cell.store_vars(v_x, c_vx);
        current_cell.store_vars(v_y, c_vy);
        current_cell.store_vars(-coords.y*v_x + coords.x*v_y , c_cir);
    }
};

#endif /* REALSCALARCIRCULATION_HPP_ */
