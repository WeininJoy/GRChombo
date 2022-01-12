/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef COMPLEXSTATICVORTEXBG_HPP_
#define COMPLEXSTATICVORTEXBG_HPP_

#include <cmath> 
#include "ADMFixedBGVars.hpp"
#include "Cell.hpp"
#include "Coordinates.hpp"
#include "CoordinateTransformations.hpp"
#include "DimensionDefinitions.hpp"
#include "Tensor.hpp"
#include "TensorAlgebra.hpp"
#include "UserVariables.hpp" //This files needs NUM_VARS - total number of components
#include "simd.hpp"

//! Class which computes the initial conditions for a static vortex
// starting with cylindrical coords
// ds^2 = -dt^2 + exp(-4*pi*G*C^2*r^2n)*(dr^2 + r^2 d\psi^2) + dz^2
// \phi = C*r^n * exp(i*n*\psi)
// Use coordinate transform to transfer g, K, shift in cylindrical coords
// to h, A, shift in Cartesian coords 

class ComplexStaticVortexBG
{
  public:
    //! Struct for the params of the static vortex
    struct params_t
    {
        double Amp;                      //!<< The amplitude (C) of the static vortex
        int n;                             //!<< The winding number (n) of the static vortex
        std::array<double, CH_SPACEDIM> center; //!< The center of the vortex
    };

    template <class data_t> using Vars = ADMFixedBGVars::Vars<data_t>;

    const params_t m_params;
    const double m_dx;

    ComplexStaticVortexBG(params_t a_params, double a_dx)
        : m_params(a_params), m_dx(a_dx) {}

    /// This just calculates chi which helps with regridding, debug etc
    /// it is only done once on setup as the BG is fixed

    template <class data_t> void compute(Cell<data_t> current_cell) const
    {
        // get position and set vars

        Vars<data_t> metric_vars;
        compute_metric_background(metric_vars, current_cell);

        // calculate and save chi
        // data_t chi = TensorAlgebra::compute_determinant_sym(metric_vars.gamma);
        // chi = pow(chi, -1.0 / 3.0);
        data_t chi = metric_vars.gamma[0][0];
        current_cell.store_vars(chi, c_chi);

        // save gamma[0][0], K, shift[0], d1_gamma[0][0][0]
        current_cell.store_vars(metric_vars.gamma[0][0], c_gamma_00);
        current_cell.store_vars(metric_vars.K_tensor[0][0], c_K_00);
        current_cell.store_vars(metric_vars.K, c_K);
        current_cell.store_vars(metric_vars.shift[0], c_shift_0);
        current_cell.store_vars(metric_vars.d1_gamma[0][0][0], c_d1_gamma_000);
    }

    template <class data_t, template <typename> class vars_t>
    void compute_metric_background(vars_t<data_t> &vars,
                                   const Cell<data_t> &current_cell) const
    {
        // set up vars for the metric and extrinsic curvature, shift and lapse in
        // cylindrical coords
        Tensor<2, data_t> cylindrical_g;
        Tensor<2, data_t> cylindrical_K;
        Tensor<1, data_t> cylindrical_shift;
        data_t lapse;
        Tensor<2, Tensor<1, data_t>> cylindrical_d1_g;
        Tensor<1, data_t> cylindrical_d1_lapse;
        Tensor<2, data_t> cylindrical_d1_shift;
        
        // get position
        Coordinates<data_t> coords(current_cell, m_dx, m_params.center);
        
        // Compute the components in cylindrical coords
        compute_metric_cylindrical(cylindrical_g, cylindrical_K, cylindrical_shift, lapse, 
                            cylindrical_d1_g, cylindrical_d1_lapse, cylindrical_d1_shift, coords);
        
        // work out where we are on the grid
        data_t x = coords.x;
        double y = coords.y;
        double z = coords.z;
        // the radius in xy plane, subject to a floor
        data_t rho2 = simd_max(x * x + y * y, 1e-1);
        data_t rho = sqrt(rho2);

        using namespace CoordinateTransformations;
        // Convert cylindrical components to cartesian components using coordinate
        // transforms
        vars.gamma = cylindrical_to_cartesian_LL(cylindrical_g, x, y, z);
        vars.shift = cylindrical_to_cartesian_U(cylindrical_shift, x, y, z);
        vars.K_tensor = cylindrical_to_cartesian_LL(cylindrical_K, x, y, z);
        vars.d1_gamma = cylindrical_to_cartesian_LLL(cylindrical_d1_g, x, y, z);
        vars.d1_shift = cylindrical_to_cartesian_UL(cylindrical_d1_shift, x, y, z);
        vars.d1_lapse = cylindrical_to_cartesian_L(cylindrical_d1_lapse, x, y, z);
        vars.lapse = lapse;


        using namespace TensorAlgebra;
        // Convert to BSSN vars
        data_t deth = compute_determinant(vars.gamma);
        const auto gamma_UU = compute_inverse_sym(vars.gamma);

        // transform extrinsic curvature into A and TrK - note h is still non
        // conformal version which is what we need here
        vars.K = compute_trace(vars.K_tensor, gamma_UU);
        make_trace_free(vars.K_tensor, vars.gamma, gamma_UU);
        

        // Populate the variables on the grid
        // NB We stil need to set Gamma^i which is NON ZERO
        // but we do this via a separate class/compute function
        // as we need the gradients of the metric which are not yet available
        
    }

    template <class data_t>
    void compute_metric_cylindrical(Tensor<2, data_t> &cylindrical_g,
                               Tensor<2, data_t> &cylindrical_K,
                               Tensor<1, data_t> &cylindrical_shift,
                               data_t &lapse,
                               Tensor<2, Tensor<1, data_t>> &cylindrical_d1_g,
                               Tensor<1, data_t> &cylindrical_d1_lapse,
                               Tensor<2, data_t> &cylindrical_d1_shift,
                               const Coordinates<data_t> coords) const
    {
        // Static vortex params - amplitude Amp and winding number n
        double Amp = m_params.Amp;
        double n = m_params.n;

        // work out where we are on the grid
        data_t x = coords.x;
        double y = coords.y;
        double z = coords.z;

        // the radius in xy plane, subject to a floor
        data_t rho2 = simd_max(x * x + y * y, 1e-4);
        data_t rho = sqrt(rho2);

        // Metric in static vortex coordinates, rho, phi and z
        FOR2(i, j) { cylindrical_g[i][j] = 0.0; }
        cylindrical_g[0][0] = exp(-4.0*M_PI*rho2);        // gamma_rr
        cylindrical_g[1][1] = rho2 * exp(-4.0*M_PI*rho2); // gamma_psipsi
        cylindrical_g[2][2] = 1.0;                        // gamma_zz

        // Extrinsic curvature. K_ij are 0 since the system is static, and 
        // K_ij = \partial_t \gamma_ij =0 
        FOR2(i, j) { cylindrical_K[i][j] = 0.0; }

        // set the analytic lapse
        lapse = 1.0;

        // set the shift
        FOR1(i) { cylindrical_shift[i] = 0.0; }
        
        // Calculate partial derivative of spatial metric
        FOR3(i, j, k) { cylindrical_d1_g[i][j][k] = 0.0; }
        cylindrical_d1_g[0][0][0] = -8.0*M_PI*rho* cylindrical_g[0][0];
        cylindrical_d1_g[1][1][0] = (2.0/rho -8.0*M_PI*rho)* cylindrical_g[1][1];

        // calculate derivs of lapse and shift
        FOR1(i) {cylindrical_d1_lapse[i] = 0.0; }
        FOR2(i, j) {cylindrical_d1_shift[i][j] = 0.0; }
    }    
};

#endif /* COMPLEXSTATICVORTEXBG_HPP_ */