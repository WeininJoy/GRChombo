/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef COMPLEXSTATICVORTEX_HPP_
#define COMPLEXSTATICVORTEX_HPP_

#include <cmath> 
#include "ADMConformalVars.hpp"
#include "Cell.hpp"
#include "Coordinates.hpp"
#include "CoordinateTransformations.hpp"
#include "DimensionDefinitions.hpp"
#include "Tensor.hpp"
#include "TensorAlgebra.hpp"
#include "UserVariables.hpp" //This files needs NUM_VARS - total number of components
#include "VarsTools.hpp"
#include "simd.hpp"

//! Class which computes the initial conditions for a static vortex
// starting with cylindrical coords
// ds^2 = -dt^2 + exp(-4*pi*G*C^2*r^2n)*(dr^2 + r^2 d\psi^2) + dz^2
// \phi = C*r^n * exp(i*n*\psi)
// Use coordinate transform to transfer g, K, shift in cylindrical coords
// to h, A, shift in Cartesian coords 

class ComplexStaticVortex
{
    // Use the variable definition in CCZ4
    template <class data_t> using Vars = ADMConformalVars::VarsWithGauge<data_t>;    
  
  public:
    //! Struct for the params of the static vortex
    struct params_t
    {
        double Amp;                      //!<< The amplitude (C) of the static vortex
        int n;                             //!<< The winding number (n) of the static vortex
        std::array<double, CH_SPACEDIM> center; //!< The center of the vortex
        double G_Newton;
    };

  protected:
    const params_t m_params;
    const double m_dx;

  public:
    ComplexStaticVortex(params_t a_params, double a_dx)
        : m_params(a_params), m_dx(a_dx) {}

    template <class data_t> void compute(Cell<data_t> current_cell) const
    {
        // set up vars for the metric and extrinsic curvature, shift and lapse in
        // cylindrical coords
        Tensor<2, data_t> cylindrical_g;
        Tensor<2, data_t> cylindrical_K;
        Tensor<1, data_t> cylindrical_shift;
        data_t vortex_lapse;
        
        // The cartesian variables and coords
        Vars<data_t> vars;
        Coordinates<data_t> coords(current_cell, m_dx, m_params.center);
        
        // Compute the components in cylindrical coords
        compute_vortex(cylindrical_g, cylindrical_K, cylindrical_shift, vortex_lapse, coords);
        
        // work out where we are on the grid
        data_t x = coords.x;
        double y = coords.y;
        double z = coords.z;
        // the radius in xy plane, subject to a floor
        data_t rho2 = simd_max(x * x + y * y, 1e-12);
        data_t rho = sqrt(rho2);

        using namespace CoordinateTransformations;
        // Convert cylindrical components to cartesian components using coordinate
        // transforms
        vars.h = cylindrical_to_cartesian_LL(cylindrical_g, x, y, z);
        vars.A = cylindrical_to_cartesian_LL(cylindrical_K, x, y, z);
        vars.shift = cylindrical_to_cartesian_U(cylindrical_shift, x, y, z);

        using namespace TensorAlgebra;
        // Convert to BSSN vars
        data_t deth = compute_determinant(vars.h);
        auto h_UU = compute_inverse_sym(vars.h);
        vars.chi = pow(deth, -1. / 3.);

        // transform extrinsic curvature into A and TrK - note h is still non
        // conformal version which is what we need here
        vars.K = compute_trace(vars.A, h_UU);
        make_trace_free(vars.A, vars.h, h_UU);
        
        // Make conformal
        FOR(i, j)
        {
            vars.h[i][j] *= vars.chi;
            vars.A[i][j] *= vars.chi;
        }

        // use a pre collapsed lapse, could also use analytic one
        // vars.lapse = vortex_lapse;
        vars.lapse = pow(vars.chi, 0.5);
        
        // Populate the variables on the grid
        // NB We stil need to set Gamma^i which is NON ZERO
        // but we do this via a separate class/compute function
        // as we need the gradients of the metric which are not yet available
        current_cell.store_vars(vars);
    }

  protected:

    template <class data_t>
    void compute_vortex(Tensor<2, data_t> &cylindrical_g,
                        Tensor<2, data_t> &cylindrical_K,
                        Tensor<1, data_t> &cylindrical_shift,
                        data_t &vortex_lapse,
                        const Coordinates<data_t> coords) const
    {
        // Static vortex params - amplitude Amp and winding number n
        double Amp = m_params.Amp;
        double n = m_params.n;
        double G_Newton = m_params.G_Newton;

        // work out where we are on the grid
        data_t x = coords.x;
        double y = coords.y;
        double z = coords.z;

        // the radius in xy plane, subject to a floor
        data_t rho2 = simd_max(x * x + y * y, 1e-12);
        data_t rho = sqrt(rho2);

        // Metric in static vortex coordinates, rho, phi and z
        FOR2(i, j) { cylindrical_g[i][j] = 0.0; }
        cylindrical_g[0][0] = exp(-4.0*M_PI*G_Newton*Amp*Amp* pow(rho2,n));        // gamma_rr
        cylindrical_g[1][1] = rho2 * exp(-4.0*M_PI*G_Newton*Amp*Amp* pow(rho2,n)); // gamma_psipsi
        // cylindrical_g[0][0] = 1.0;        // gamma_rr
        // cylindrical_g[1][1] = rho2;
        cylindrical_g[2][2] = 1.0;                        // gamma_zz

        // Extrinsic curvature. K_ij are 0 since the system is static, and 
        // K_ij = \partial_t \gamma_ij =0 
        FOR2(i, j) { cylindrical_K[i][j] = 0.0; }

        // set the analytic lapse
        vortex_lapse = 1.0;

        // set the shift
        FOR1(i) { cylindrical_shift[i] = 0.0; }
    }    
};

#endif /* COMPLEXSTATICVORTEX_HPP_ */
